import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union


def _find_structure_file(work_dir: Path) -> Optional[Path]:
    """Busca o arquivo de estrutura (md_clean.gro prioritário, complex.pdb, complex.gro ou md.gro) no work_dir ou subdiretórios (suporta prefixos)."""
    candidates = ["md_clean.gro", "complex.pdb", "complex.gro", "md.gro"]
    for cand in candidates:
        direct = work_dir / cand
        if direct.exists():
            return direct
        prefixed = list(work_dir.glob(f"*_{cand}"))
        if prefixed:
            return prefixed[0]
    # Busca recursiva rasa (1 nível)
    for cand in candidates:
        matches = list(work_dir.glob(f"*/{cand}")) or list(work_dir.glob(f"*/*_{cand}"))
        if matches:
            return matches[0]
    # Busca em diretório md_files adjacente/pai
    for parent in [work_dir] + list(work_dir.parents):
        candidate_md = parent / "md_files"
        if candidate_md.exists():
            for cand in candidates:
                cand_file = candidate_md / cand
                if cand_file.exists():
                    return cand_file
                prefixed = list(candidate_md.glob(f"*_{cand}"))
                if prefixed:
                    return prefixed[0]
    return None


def _normalize_residues(residues: Optional[Sequence[Any]]) -> List[str]:
    """
    Normaliza sequências de identificadores de resíduos (int, str, dict do PLIP ou objetos)
    em uma lista ordenada e única de strings numéricas/identificadores para seleção no PyMOL.
    Exemplos aceitos:
      - [57, 102, 195] -> ['57', '102', '195']
      - ['57', '102'] -> ['57', '102']
      - ['His57', 'Asp102'] -> ['57', '102']
      - [{'resnr': 57, 'resname': 'HIS'}, ...] -> ['57']
    """
    if not residues:
        return []
    result: List[str] = []
    for item in residues:
        if item is None:
            continue
        if isinstance(item, dict):
            res = item.get("resnr") or item.get("resi") or item.get("residue")
            if res is not None:
                result.append(str(res).strip())
        elif isinstance(item, (int, float)):
            result.append(str(int(item)))
        elif isinstance(item, str):
            parts = [p.strip() for p in re.split(r"[,+\s]+", item.strip()) if p.strip()]
            for part in parts:
                m = re.search(r"(\d+[A-Za-z]?)", part)
                if m:
                    result.append(m.group(1))
                else:
                    result.append(part)
        else:
            res = getattr(item, "resnr", None) or getattr(item, "resi", None)
            if res is not None:
                result.append(str(res).strip())

    seen: Set[str] = set()
    unique: List[str] = []
    for r in result:
        if r not in seen:
            seen.add(r)
            unique.append(r)

    try:
        def _norm_sort_key(x: str):
            m = re.match(r"^(\d+)([A-Za-z]?)$", str(x).strip())
            if m:
                num, ic = m.groups()
                return (int(num), ic)
            return (999999, str(x))

        unique.sort(key=_norm_sort_key)
    except Exception:
        pass

    return unique


def _extract_interactions_and_residues(
    work_dir: Path, interactions_data: Optional[Dict[str, Any]] = None
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Any], Dict[str, List[Dict[str, Any]]]]:
    """
    Recupera todas as interações (pontes de hidrogênio, contatos hidrofóbicos, pontes salinas,
    pi-stacking, pi-cátion, ligações de halogênio, complexos metálicos e pontes de água)
    e a lista consolidada de resíduos-chave a partir de interactions_data em memória
    ou localizando o arquivo interactions.json.
    """
    all_interactions: Dict[str, List[Dict[str, Any]]] = {
        "hydrogen_bonds": [],
        "hydrophobic_contacts": [],
        "salt_bridges": [],
        "pi_stacks": [],
        "pi_cation_interactions": [],
        "halogen_bonds": [],
        "metal_complexes": [],
        "water_bridges": [],
    }
    key_residue_numbers: Set[Any] = set()

    inter_dict: Optional[Dict[str, Any]] = None
    if interactions_data is not None:
        inter_dict = interactions_data
    else:
        interactions_file = work_dir / "interactions.json"
        if not interactions_file.exists():
            matches = (
                list(work_dir.glob("*_interactions.json"))
                or list(work_dir.glob("*/interactions.json"))
                or list(work_dir.glob("*/*_interactions.json"))
            )
            if matches:
                interactions_file = matches[0]

        if interactions_file and interactions_file.exists():
            try:
                with open(interactions_file, "r", encoding="utf-8") as f:
                    inter_dict = json.load(f)
            except Exception:
                inter_dict = None

    if inter_dict:
        for cat in all_interactions:
            items = inter_dict.get(cat, [])
            if isinstance(items, list):
                all_interactions[cat] = items
                for item in items:
                    resnr = item.get("resnr")
                    if resnr is not None and str(resnr).strip() not in ("", "0"):
                        key_residue_numbers.add(resnr)

    hbonds = all_interactions["hydrogen_bonds"]
    hcontacts = all_interactions["hydrophobic_contacts"]

    def _res_sort_key(val: Any):
        s = str(val).strip()
        m = re.match(r"^(\d+)([A-Za-z]?)$", s)
        if m:
            num, ic = m.groups()
            return (int(num), ic)
        return (999999, s)

    return (
        hbonds,
        hcontacts,
        sorted(list(key_residue_numbers), key=_res_sort_key),
        all_interactions,
    )


def generate_pymol_script(
    work_dir: Union[str, Path],
    key_residues: Optional[Sequence[Any]] = None,
    catalytic_residues: Optional[Sequence[Any]] = None,
    interactions_data: Optional[Dict[str, Any]] = None,
) -> Path:
    """
    Gera um script automatizado do PyMOL (show_complex.pml) com padrão estético de publicação científica:
    - Fundo branco (bg_color white), ray_shadows 0, antialias 2, depth_cue 0, specular 0.1
    - Proteína em cartoon branco, sobreposta com superfície também branca e transparência 0.85
    - Ligante em bastões (sticks) com stick_radius 0.25 (carbonos em magenta, oxigênios em vermelho, nitrogênios em azul)
    - Resíduos-chave dinâmicos (PLIP) em bastões com carbonos em gray80 (sem valores fixos/chumbados no código)
    - Resíduos catalíticos opcionais destacados em verde nos carbonos
    - Interações de pontes de hidrogênio tracejadas em deepblue (dash_width 4.0, dash_gap 0.3) sem rótulos numéricos
    - Legendas de resíduos em preto, fonte limpa 7, tamanho 26, posição [0, 0, 1.5] e formato contínuo (ex: 'His57')
    - Enquadramento final centralizado e focado no ligante com zoom 4.5

    :param work_dir: Diretório contendo a estrutura (complex.pdb / md_clean.pdb) e opcionalmente interactions.json.
    :param key_residues: Lista opcional de resíduos-chave (ex: ints, strings 'His57', ou dicts do PLIP). Se None, extrai do interactions.json.
    :param catalytic_residues: Lista opcional de resíduos catalíticos para destacar os carbonos em verde.
    :param interactions_data: Dicionário opcional contendo 'hydrogen_bonds' e 'hydrophobic_contacts'.
    :return: Path para o arquivo show_complex.pml gerado.
    """
    work_dir = Path(work_dir).resolve()
    if not work_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {work_dir}")

    struct_file = _find_structure_file(work_dir)
    if not struct_file:
        raise FileNotFoundError(
            f"Nenhum arquivo de estrutura ('md_clean.gro' ou 'complex.pdb') encontrado no diretório: {work_dir}"
        )

    # 0. Tenta gerar versões PDB sem água para carregamento e renderização instantâneos no PyMOL
    gmx_bin = shutil.which("gmx")
    if not gmx_bin:
        try:
            from docking.md_prep import find_executable

            gmx_bin = find_executable("gmx")
        except Exception:
            gmx_bin = None

    tpr_file = work_dir / "md.tpr"
    if not tpr_file.exists():
        matches_tpr = list(work_dir.glob("*_md.tpr"))
        if matches_tpr:
            tpr_file = matches_tpr[0]

    index_file = work_dir / "index.ndx"
    if not index_file.exists():
        matches_ndx = list(work_dir.glob("*_index.ndx"))
        if matches_ndx:
            index_file = matches_ndx[0]

    nowat_pdb = work_dir / "md_clean_nowat.pdb"
    if (
        not nowat_pdb.exists()
        and gmx_bin
        and struct_file
        and struct_file.exists()
        and tpr_file.exists()
        and index_file.exists()
    ):
        try:
            env = os.environ.copy()
            subprocess.run(
                [
                    gmx_bin,
                    "trjconv",
                    "-s",
                    str(tpr_file.name),
                    "-f",
                    str(struct_file.name),
                    "-n",
                    str(index_file.name),
                    "-o",
                    "md_clean_nowat.pdb",
                ],
                cwd=str(work_dir),
                input=b"Protein_LIG\n",
                capture_output=True,
                env=env,
                check=False,
            )
        except Exception:
            pass

    clean_pdb = work_dir / "md_clean.pdb"
    if (
        not clean_pdb.exists()
        and not nowat_pdb.exists()
        and gmx_bin
        and struct_file
        and struct_file.exists()
    ):
        try:
            env = os.environ.copy()
            subprocess.run(
                [
                    gmx_bin,
                    "editconf",
                    "-f",
                    str(struct_file.name),
                    "-o",
                    "md_clean.pdb",
                ],
                cwd=str(work_dir),
                capture_output=True,
                env=env,
                check=False,
            )
        except Exception:
            pass

    # Verifica se a estrutura medóide mais representativa do cluster está presente
    medoid_file = work_dir / "cluster_medoid.gro"
    if not medoid_file.exists():
        matches_medoid = list(work_dir.glob("*cluster_medoid.gro")) or list(
            work_dir.glob("*/cluster_medoid.gro")
        )
        if matches_medoid:
            medoid_file = matches_medoid[0]

    medoid_pdb = work_dir / "cluster_medoid.pdb"
    if not medoid_pdb.exists() and gmx_bin and medoid_file and medoid_file.exists():
        try:
            env = os.environ.copy()
            subprocess.run(
                [
                    gmx_bin,
                    "editconf",
                    "-f",
                    str(medoid_file.name),
                    "-o",
                    "cluster_medoid.pdb",
                ],
                cwd=str(work_dir),
                capture_output=True,
                env=env,
                check=False,
            )
        except Exception:
            pass

    # Trajetória ajustada (PBC Corrigido & Fit rot+trans)
    fit_xtc_file = work_dir / "md_fit.xtc"
    if not fit_xtc_file.exists():
        matches_xtc = list(work_dir.glob("*md_fit.xtc")) or list(
            work_dir.glob("*/md_fit.xtc")
        )
        if matches_xtc:
            fit_xtc_file = matches_xtc[0]

    nowat_xtc = work_dir / "md_fit_nowat.xtc"
    if (
        not nowat_xtc.exists()
        and gmx_bin
        and fit_xtc_file
        and fit_xtc_file.exists()
        and tpr_file.exists()
        and index_file.exists()
    ):
        try:
            env = os.environ.copy()
            subprocess.run(
                [
                    gmx_bin,
                    "trjconv",
                    "-s",
                    str(tpr_file.name),
                    "-f",
                    str(fit_xtc_file.name),
                    "-n",
                    str(index_file.name),
                    "-o",
                    "md_fit_nowat.xtc",
                    "-dt",
                    "100",
                ],
                cwd=str(work_dir),
                input=b"Protein_LIG\n",
                capture_output=True,
                env=env,
                check=False,
            )
        except Exception:
            pass

    # Escolhe o melhor arquivo de estrutura primária e trajetória compatível
    primary_struct = (
        "md_clean_nowat.pdb"
        if nowat_pdb.exists()
        else ("md_clean.pdb" if clean_pdb.exists() else struct_file.name)
    )
    primary_xtc = (
        "md_fit_nowat.xtc"
        if (nowat_pdb.exists() and nowat_xtc.exists())
        else (
            fit_xtc_file.name
            if fit_xtc_file and fit_xtc_file.exists() and not nowat_pdb.exists()
            else None
        )
    )
    primary_medoid = (
        "cluster_medoid.pdb"
        if medoid_pdb.exists()
        else (medoid_file.name if medoid_file and medoid_file.exists() else None)
    )

    # Extração e normalização de dados de interação
    hbonds, hcontacts, extracted_resnrs, all_interactions = (
        _extract_interactions_and_residues(work_dir, interactions_data)
    )

    # Resíduos-chave dinâmicos: argumento explícito > extração PLIP
    if key_residues is not None:
        key_res_list = _normalize_residues(key_residues)
    else:
        key_res_list = _normalize_residues(extracted_resnrs)

    # Resíduos catalíticos opcionais
    cat_res_list = _normalize_residues(catalytic_residues)

    work_dir_posix = str(work_dir).replace("\\", "/")

    # Monta comandos do script PyMOL
    pml_lines = [
        "# ==============================================================================",
        "# PyMOL Automated Visualization Script",
        "# Generated automatically by Molecular Docking Pipeline",
        "# Publication Quality Preset",
        "# ==============================================================================",
        "",
        "# 0. Garantia de Diretório de Trabalho Autônomo (Python API do PyMOL)",
        "python",
        "import os",
        "try:",
        f"    os.chdir(r'{work_dir_posix}')",
        "except Exception:",
        "    pass",
        "python end",
        "",
        "# 1. Iluminação, Fundo e Configurações de Renderização",
        "reinitialize",
        "bg_color white",
        "set ray_shadows, 0",
        "set antialias, 2",
        "set depth_cue, 0",
        "set specular, 0.1",
        "set cartoon_fancy_helices, 1",
        "set cartoon_smooth_loops, 1",
        "",
        "# 2. Carregamento do Complexo Receptor-Ligante",
        f"load {primary_struct}, complex",
    ]

    if primary_xtc:
        pml_lines.extend(
            [
                "",
                "# 2.1 Carregamento da Trajetória Ajustada (PBC Corrigido & Fit rot+trans)",
                f"load_traj {primary_xtc}, complex",
            ]
        )

    if primary_medoid:
        pml_lines.extend(
            [
                "",
                "# 2.2 Estrutura Representativa do Cluster de Equilíbrio (GROMOS Medoid)",
                f"load {primary_medoid}, rep_cluster",
                "remove (rep_cluster and (resn SOL or resn HOH or resn NA or resn CL or resn TIP3 or resn ION))",
                "hide everything, rep_cluster",
                "show cartoon, rep_cluster and polymer",
                "color warmpink, rep_cluster and polymer",
                "set cartoon_transparency, 0.45, rep_cluster and polymer",
                "select rep_ligand, rep_cluster and (resn LIG or resn UNK or resn UNL or resn MOL or resn ligand_md or (not polymer and not solvent))",
                "show sticks, rep_ligand",
                "color orange, rep_ligand",
                "util.cnc rep_ligand",
                "set stick_radius, 0.22, rep_ligand",
                "# Disponível para alternar visibilidade (overlay)",
                "disable rep_cluster",
            ]
        )

    pml_lines.extend(
        [
            "",
            "# 3. Limpeza de Solvente e Representação da Proteína (Cartoon Branco + Superfície Translúcida)",
            "remove resn SOL or resn HOH or resn TIP3 or resn NA or resn CL or resn ION",
            "hide everything, complex",
            "show cartoon, complex and polymer",
            "color white, complex and polymer",
            "show surface, complex and polymer",
            "color white, complex and polymer",
            "set transparency, 0.85, complex and polymer",
            "",
            "# 4. Representação do Ligante (Sticks)",
            "select ligand, complex and (resn LIG or resn UNK or resn UNL or resn MOL or resn ligand_md or (not polymer and not solvent))",
            "show sticks, ligand",
            "set stick_radius, 0.25, ligand",
            "color magenta, ligand",
            "color magenta, ligand and elem C",
            "color red, ligand and elem O",
            "color blue, ligand and elem N",
            "",
        ]
    )

    # 5. Resíduos Dinâmicos (PLIP) e Resíduos Catalíticos Opcionais
    if key_res_list:
        resi_selection = "+".join(key_res_list)
        pml_lines.extend(
            [
                "# 5. Resíduos Chave de Interação Mapeados pelo PLIP (Carbonos em gray80)",
                f"select key_residues, polymer and resi {resi_selection}",
                "show sticks, key_residues",
                "set stick_radius, 0.20, key_residues",
                "color gray80, key_residues and elem C",
                "color red, key_residues and elem O",
                "color blue, key_residues and elem N",
                "",
            ]
        )
    else:
        pml_lines.extend(
            [
                "# 5. Resíduos do Sítio de Ligação (Raio de Proximidade 5Å)",
                "select key_residues, polymer within 5.0 of ligand",
                "show sticks, key_residues",
                "set stick_radius, 0.20, key_residues",
                "color gray80, key_residues and elem C",
                "color red, key_residues and elem O",
                "color blue, key_residues and elem N",
                "",
            ]
        )

    if cat_res_list:
        cat_selection = "+".join(cat_res_list)
        pml_lines.extend(
            [
                "# 5.1 Resíduos Catalíticos (Carbonos em Verde)",
                f"select catalytic_residues, polymer and resi {cat_selection}",
                "show sticks, catalytic_residues",
                "set stick_radius, 0.20, catalytic_residues",
                "color green, catalytic_residues and elem C",
                "color red, catalytic_residues and elem O",
                "color blue, catalytic_residues and elem N",
                "",
            ]
        )

    # 5.2 Legendas e Rótulos Contínuos (ex: 'His57')
    if key_res_list and cat_res_list:
        label_target = "(key_residues or catalytic_residues)"
    elif key_res_list:
        label_target = "key_residues"
    elif cat_res_list:
        label_target = "catalytic_residues"
    else:
        label_target = "key_residues"

    pml_lines.extend(
        [
            "# 5.2 Rótulos dos Resíduos em Preto (Tamanho 26, Fonte 7, Posição [0, 0, 1.5], Formato 'His57')",
            f'label {label_target} and name CA, "%s%s" % (resn.capitalize(), resi)',
            "set label_color, black",
            "set label_font_id, 7",
            "set label_size, 26",
            "set label_position, [0, 0, 1.5]",
            "",
        ]
    )

    # 6. Interações e Linhas Tracejadas
    if hbonds:
        pml_lines.append(
            "# 6.1 Pontes de Hidrogênio (Linhas Tracejadas Deepblue sem Rótulos Numéricos)"
        )
        unique_hb_pairs: Set[Any] = set()
        for hb in hbonds:
            resnr = hb.get("resnr")
            resname = str(hb.get("resname", "RES")).capitalize()
            if resnr and resnr not in unique_hb_pairs:
                unique_hb_pairs.add(resnr)
                dist_name = f"hb_{resname}_{resnr}"
                pml_lines.append(f"# H-Bond {resname}{resnr}")
                pml_lines.append(
                    f"distance {dist_name}, (polymer and resi {resnr}), (ligand), 4.2, mode=2"
                )
                pml_lines.append(f"hide labels, {dist_name}")
                pml_lines.append(f"set dash_color, deepblue, {dist_name}")
                pml_lines.append(f"set dash_gap, 0.3, {dist_name}")
                pml_lines.append(f"set dash_width, 4.0, {dist_name}")
                pml_lines.append(f"set dash_radius, 0.05, {dist_name}")
        pml_lines.append("")

    salt_bridges = all_interactions.get("salt_bridges", [])
    if salt_bridges:
        pml_lines.append(
            "# 6.2 Pontes Salinas (Linhas Tracejadas Warmpink)"
        )
        unique_sb_pairs: Set[Any] = set()
        for sb in salt_bridges:
            resnr = sb.get("resnr")
            resname = str(sb.get("resname", "RES")).capitalize()
            if resnr and resnr not in unique_sb_pairs:
                unique_sb_pairs.add(resnr)
                dist_name = f"sb_{resname}_{resnr}"
                pml_lines.append(f"# Salt Bridge {resname}{resnr}")
                pml_lines.append(
                    f"distance {dist_name}, (polymer and resi {resnr}), (ligand), 5.0"
                )
                pml_lines.append(f"hide labels, {dist_name}")
                pml_lines.append(f"set dash_color, warmpink, {dist_name}")
                pml_lines.append(f"set dash_gap, 0.3, {dist_name}")
                pml_lines.append(f"set dash_width, 4.0, {dist_name}")
                pml_lines.append(f"set dash_radius, 0.05, {dist_name}")
        pml_lines.append("")

    pi_interactions = (
        all_interactions.get("pi_cation_interactions", [])
        + all_interactions.get("pi_stacks", [])
    )
    if pi_interactions:
        pml_lines.append(
            "# 6.3 Interações Pi (Cátion-Pi e Pi-Stacking em Forest Green)"
        )
        unique_pi_pairs: Set[Any] = set()
        for pi_item in pi_interactions:
            resnr = pi_item.get("resnr")
            resname = str(pi_item.get("resname", "RES")).capitalize()
            if resnr and resnr not in unique_pi_pairs:
                unique_pi_pairs.add(resnr)
                dist_name = f"pi_{resname}_{resnr}"
                pml_lines.append(f"# Pi Interaction {resname}{resnr}")
                pml_lines.append(
                    f"distance {dist_name}, (polymer and resi {resnr}), (ligand), 5.0"
                )
                pml_lines.append(f"hide labels, {dist_name}")
                pml_lines.append(f"set dash_color, forest, {dist_name}")
                pml_lines.append(f"set dash_gap, 0.3, {dist_name}")
                pml_lines.append(f"set dash_width, 4.0, {dist_name}")
                pml_lines.append(f"set dash_radius, 0.05, {dist_name}")
        pml_lines.append("")

    halogens = all_interactions.get("halogen_bonds", [])
    if halogens:
        pml_lines.append(
            "# 6.4 Ligações de Halogênio (Linhas Tracejadas em Ciano)"
        )
        unique_hal_pairs: Set[Any] = set()
        for hg in halogens:
            resnr = hg.get("resnr")
            resname = str(hg.get("resname", "RES")).capitalize()
            if resnr and resnr not in unique_hal_pairs:
                unique_hal_pairs.add(resnr)
                dist_name = f"hal_{resname}_{resnr}"
                pml_lines.append(f"# Halogen Bond {resname}{resnr}")
                pml_lines.append(
                    f"distance {dist_name}, (polymer and resi {resnr}), (ligand), 4.5"
                )
                pml_lines.append(f"hide labels, {dist_name}")
                pml_lines.append(f"set dash_color, cyan, {dist_name}")
                pml_lines.append(f"set dash_gap, 0.3, {dist_name}")
                pml_lines.append(f"set dash_width, 4.0, {dist_name}")
                pml_lines.append(f"set dash_radius, 0.05, {dist_name}")
        pml_lines.append("")

    # 7. Enquadramento e Foco Final no Ligante
    pml_lines.extend(
        [
            "# 7. Enquadramento e Foco Final no Sítio de Ligação",
            "deselect",
            "center ligand",
            "zoom ligand, 4.5",
            "",
        ]
    )

    output_pml = work_dir / "show_complex.pml"
    with open(output_pml, "w", encoding="utf-8") as f:
        f.write("\n".join(pml_lines) + "\n")

    return output_pml
