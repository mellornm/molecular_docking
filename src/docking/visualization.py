import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Set


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


def generate_pymol_script(work_dir: Path) -> Path:
    """
    Gera um script automatizado do PyMOL (show_complex.pml) no diretório de trabalho,
    configurando representação visual científica de alta fidelidade:
    - Suporte a execução independente de caminho (os.chdir embutido via Python API do PyMOL)
    - Conversão e carregamento de estrutura limpa sem solvente (md_clean_nowat.pdb / md_clean.pdb)
    - Carregamento da trajetória ajustada md_fit.xtc
    - Fundo branco de publicação (set bg_rgb, [1, 1, 1])
    - Proteína em cartoon ciano suave (color cyan, polymer)
    - Ligante (LIG) em bastões coloridos por elemento (color magenta)
    - Seleção e destaque em bastões dos resíduos do sítio ativo mapeados no interactions.json
    - Linhas tracejadas amarelas para pontes de hidrogênio com rótulos de distância
    - Foco e enquadramento automático do ligante (center, zoom, orient)

    :param work_dir: Diretório de trabalho contendo md_clean.gro / complex.pdb e interactions.json.
    :return: Caminho do arquivo show_complex.pml gerado.
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
        from docking.md_prep import find_executable

        gmx_bin = find_executable("gmx")

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

    # Leitura do interactions.json (se existir)
    interactions_file = work_dir / "interactions.json"
    if not interactions_file.exists():
        matches = (
            list(work_dir.glob("*_interactions.json"))
            or list(work_dir.glob("*/interactions.json"))
            or list(work_dir.glob("*/*_interactions.json"))
        )
        if matches:
            interactions_file = matches[0]

    hbonds: List[Dict[str, Any]] = []
    hcontacts: List[Dict[str, Any]] = []

    if interactions_file and interactions_file.exists():
        try:
            with open(interactions_file, "r", encoding="utf-8") as f:
                inter_data = json.load(f)
                hbonds = inter_data.get("hydrogen_bonds", [])
                hcontacts = inter_data.get("hydrophobic_contacts", [])
        except Exception:
            hbonds = []
            hcontacts = []

    # Extrai conjunto de resíduos únicos do sítio ativo
    key_residue_numbers: Set[int] = set()
    for hb in hbonds:
        resnr = hb.get("resnr")
        if resnr:
            key_residue_numbers.add(int(resnr))

    for hc in hcontacts:
        resnr = hc.get("resnr")
        if resnr:
            key_residue_numbers.add(int(resnr))

    sorted_resnrs = sorted(list(key_residue_numbers))

    work_dir_posix = str(work_dir).replace("\\", "/")

    # Monta comandos do script PyMOL
    pml_lines = [
        "# ==============================================================================",
        "# PyMOL Automated Visualization Script",
        "# Generated automatically by Molecular Docking Pipeline",
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
        "# 1. Inicialização e Configurações de Fundo e Renderização",
        "reinitialize",
        "set bg_rgb, [1, 1, 1]",
        "set ray_shadows, 0",
        "set antialias, 2",
        "set depth_cue, 1",
        "set specular, 0.25",
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
            "# 3. Limpeza de Solvente e Representação da Proteína (Cartoon)",
            "remove resn SOL or resn HOH or resn TIP3 or resn NA or resn CL or resn ION",
            "hide everything, complex",
            "show cartoon, complex and polymer",
            "color cyan, complex and polymer",
            "set cartoon_transparency, 0.15, complex and polymer",
            "",
            "# 4. Representação do Ligante (Sticks)",
            "select ligand, complex and (resn LIG or resn UNK or resn UNL or resn MOL or resn ligand_md or (not polymer and not solvent))",
            "show sticks, ligand",
            "color magenta, ligand",
            "util.cnc ligand",
            "set stick_radius, 0.25, ligand",
            "",
        ]
    )

    # 5. Seleção e exibição dos resíduos do sítio ativo
    if sorted_resnrs:
        resi_selection = "+".join(str(r) for r in sorted_resnrs)
        pml_lines.extend(
            [
                "# 5. Resíduos Chave de Interação Mapeados pelo PLIP",
                f"select key_residues, polymer and resi {resi_selection}",
                "show sticks, key_residues",
                "color gray80, key_residues and elem C",
                "util.cnc key_residues",
                "set stick_radius, 0.18, key_residues",
                'label key_residues and name CA, "%s %s" % (resn, resi)',
                "set label_size, 14",
                "set label_color, black",
                "set label_font_id, 7",
                "set label_position, [0, 0, 1.5]",
                "",
            ]
        )
    else:
        pml_lines.extend(
            [
                "# 5. Resíduos do Sítio de Ligação (Raio de Proximidade 5Å)",
                "select key_residues, polymer within 5.0 of ligand",
                "show sticks, key_residues",
                "color gray80, key_residues and elem C",
                "util.cnc key_residues",
                "set stick_radius, 0.18, key_residues",
                'label key_residues and name CA, "%s %s" % (resn, resi)',
                "set label_size, 14",
                "set label_color, black",
                "",
            ]
        )

    # 6. Pontes de Hidrogênio com distâncias e linhas tracejadas amarelas
    if hbonds:
        pml_lines.append("# 6. Pontes de Hidrogênio (Linhas Tracejadas Amarelas)")
        unique_hb_pairs: Set[int] = set()
        for idx, hb in enumerate(hbonds, 1):
            resnr = hb.get("resnr")
            resname = hb.get("resname", "RES")
            if resnr and resnr not in unique_hb_pairs:
                unique_hb_pairs.add(resnr)
                pml_lines.append(f"# H-Bond {resname} {resnr}")
                pml_lines.append(
                    f"distance hb_{resname}_{resnr}, (polymer and resi {resnr}), (ligand), 4.2, mode=2"
                )

        pml_lines.extend(
            [
                "set dash_color, yellow",
                "set dash_gap, 0.25",
                "set dash_width, 3.0",
                "set dash_radius, 0.05",
                "set label_color, black",
                "set label_size, 12",
                "",
            ]
        )

    # 7. Centralização e Foco
    pml_lines.extend(
        [
            "# 7. Centralização e Foco no Sítio Ativo",
            "center ligand",
            "zoom ligand, 8",
            "orient ligand",
            "",
            "# Deselecionar tudo para limpar a visualização",
            "deselect",
            "",
        ]
    )

    output_pml = work_dir / "show_complex.pml"
    with open(output_pml, "w", encoding="utf-8") as f:
        f.write("\n".join(pml_lines) + "\n")

    return output_pml
