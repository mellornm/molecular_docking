import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Set


def _find_structure_file(work_dir: Path) -> Optional[Path]:
    """Busca o arquivo de estrutura (md_clean.gro, complex.pdb, complex.gro ou md.gro) no work_dir ou subdiretórios."""
    candidates = ["md_clean.gro", "complex.pdb", "complex.gro", "md.gro"]
    for cand in candidates:
        direct = work_dir / cand
        if direct.exists():
            return direct
    # Busca recursiva rasa (1 nível)
    for cand in candidates:
        matches = list(work_dir.glob(f"*/{cand}"))
        if matches:
            return matches[0]
    return None


def generate_pymol_script(work_dir: Path) -> Path:
    """
    Gera um script automatizado do PyMOL (show_complex.pml) no diretório de trabalho,
    configurando representação visual científica de alta fidelidade:
    - Carregamento de md_clean.gro (ou complex.pdb) e da trajetória ajustada md_fit.xtc
    - Fundo branco de publicação (set bg_rgb, [1, 1, 1])
    - Proteína em cartoon ciano suave (color cyan, polymer)
    - Ligante (LIG) em bastões coloridos por elemento (color magenta)
    - Seleção e destaque em bastões dos resíduos do sítio ativo mapeados no interactions.json
    - Linhas tracejadas amarelas para pontes de hidrogênio com rótulos de distância
    - Foco e enquadramento automático do sítio ativo (zoom)

    :param work_dir: Diretório de trabalho contendo md_clean.gro / complex.pdb e interactions.json.
    :return: Caminho do arquivo show_complex.pml gerado.
    """
    work_dir = Path(work_dir)
    if not work_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {work_dir}")

    struct_file = _find_structure_file(work_dir)
    if not struct_file:
        raise FileNotFoundError(
            f"Nenhum arquivo de estrutura ('md_clean.gro' ou 'complex.pdb') encontrado no diretório: {work_dir}"
        )

    # Leitura do interactions.json (se existir)
    interactions_file = work_dir / "interactions.json"
    if not interactions_file.exists():
        matches = list(work_dir.glob("*/interactions.json"))
        if matches:
            interactions_file = matches[0]

    hbonds: List[Dict[str, Any]] = []
    hcontacts: List[Dict[str, Any]] = []

    if interactions_file.exists():
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

    # Verifica se a trajetória tratada md_fit.xtc está presente
    fit_xtc_file = work_dir / "md_fit.xtc"
    if not fit_xtc_file.exists():
        matches_xtc = list(work_dir.glob("*/md_fit.xtc"))
        if matches_xtc:
            fit_xtc_file = matches_xtc[0]

    # Monta comandos do script PyMOL
    pml_lines = [
        "# ==============================================================================",
        "# PyMOL Automated Visualization Script",
        "# Generated automatically by Molecular Docking Pipeline",
        "# ==============================================================================",
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
        f"load {struct_file.name}, complex",
    ]

    if fit_xtc_file.exists():
        pml_lines.extend(
            [
                "",
                "# 2.1 Carregamento da Trajetória Ajustada (PBC Corrigido & Fit rot+trans)",
                f"load_traj {fit_xtc_file.name}, complex",
            ]
        )

    pml_lines.extend(
        [
            "",
            "# 3. Limpeza de Solvente e Representação da Proteína (Cartoon)",
            "remove resn SOL or resn HOH or resn TIP3 or resn NA or resn CL or resn ION",
            "hide everything, all",
            "show cartoon, polymer",
            "color cyan, polymer",
            "set cartoon_transparency, 0.15, polymer",
            "",
            "# 4. Representação do Ligante (Sticks)",
            "select ligand, (resn LIG or resn UNK or resn UNL or resn MOL or resn ligand_md or (not polymer and not solvent))",
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
                "select key_residues, polymer within 5.0 of resn LIG",
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
                    f"distance hb_{resname}_{resnr}, (polymer and resi {resnr}), (resn LIG), 4.2, mode=2"
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
            "center resn LIG",
            "zoom resn LIG, 8",
            "orient resn LIG",
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
