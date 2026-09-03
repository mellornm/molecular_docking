#!/usr/bin/env python3
"""
Script de arquivamento otimizado de experimentos de Docking e Dinâmica Molecular.
Preserva todos os dados essenciais para reanálises completas e reprodutibilidade,
suportando arquitetura de Isolamento de Alvos (Target Isolation) e prefixos explícitos.
"""

import shutil
import sys
from pathlib import Path
from typing import List, Optional


def find_target_md_dirs(data_dir: Path) -> List[Path]:
    """Descobre todos os diretórios válidos de simulação MD dentro de data/md_files/ ou data/."""
    md_base = data_dir / "md_files"
    found = []

    if md_base.exists():
        # Subpastas de alvos (ex: data/md_files/7CFN)
        subdirs = [
            d for d in md_base.iterdir() if d.is_dir() and not d.name.startswith(".")
        ]
        if subdirs:
            found.extend(subdirs)
        else:
            # Pasta raiz direta data/md_files
            found.append(md_base)
    return found


def archive_experiment(
    exp_name: Optional[str] = None, target_md_dir: Optional[Path] = None
) -> Path:
    print("=" * 68)
    print("📦 ARQUIVAMENTO OTIMIZADO DE DOCKING & DINÂMICA MOLECULAR")
    print("=" * 68)

    root_dir = Path(__file__).resolve().parent
    data_dir = root_dir / "data"

    # 1. Descoberta e resolução dos diretórios de MD
    candidate_md_dirs = find_target_md_dirs(data_dir)

    selected_md_dir = None
    if target_md_dir:
        selected_md_dir = Path(target_md_dir).resolve()
    elif exp_name:
        for d in candidate_md_dirs:
            if d.name.lower() == exp_name.lower() or exp_name.lower() in d.name.lower():
                selected_md_dir = d
                break
        if not selected_md_dir:
            custom_dir = data_dir / "md_files" / exp_name
            if custom_dir.exists():
                selected_md_dir = custom_dir

    if not selected_md_dir:
        if len(candidate_md_dirs) == 1:
            selected_md_dir = candidate_md_dirs[0]
        elif len(candidate_md_dirs) > 1:
            print("\nDiretórios de Dinâmica Molecular encontrados:")
            for idx, d in enumerate(candidate_md_dirs, 1):
                print(f"  [{idx}] {d.name} ({d})")
            if not exp_name:
                try:
                    choice = input(
                        f"\nEscolha o alvo a arquivar [1-{len(candidate_md_dirs)}] (padrão: 1): "
                    ).strip()
                    choice_idx = (
                        int(choice) - 1
                        if choice.isdigit()
                        and 1 <= int(choice) <= len(candidate_md_dirs)
                        else 0
                    )
                    selected_md_dir = candidate_md_dirs[choice_idx]
                except (EOFError, KeyboardInterrupt):
                    selected_md_dir = candidate_md_dirs[0]
            else:
                selected_md_dir = candidate_md_dirs[0]
        else:
            selected_md_dir = data_dir / "md_files"

    target_id = (
        selected_md_dir.name
        if selected_md_dir.name != "md_files"
        else (exp_name or "DEFAULT")
    )

    if not exp_name:
        default_name = (
            f"DS-{target_id}" if target_id != "DEFAULT" else "EXPERIMENT_NAME"
        )
        try:
            user_input = input(
                f"\nNome da pasta de arquivamento [{default_name}]: "
            ).strip()
            exp_name = user_input if user_input else default_name
        except (EOFError, KeyboardInterrupt):
            exp_name = default_name

    archive_base = root_dir / "archive" / exp_name

    dest_dirs = {
        "docking": archive_base / "docking",
        "topology": archive_base / "topology",
        "trajectory": archive_base / "trajectory",
        "report": archive_base / "report",
    }

    for path in dest_dirs.values():
        path.mkdir(parents=True, exist_ok=True)

    print(f"\n📁 Alvo Selecionado:       {target_id} ({selected_md_dir})")
    print(f"📁 Diretório de Destino:  {archive_base}\n")
    copied_count = 0
    total_bytes = 0
    seen_dest_files = set()
    seen_mirror_keys = set()

    def copy_file(src: Path, dest: Path, label: str) -> bool:
        nonlocal copied_count, total_bytes
        if src.exists() and src.is_file():
            dest_file = dest / src.name
            if dest_file.name in seen_dest_files:
                return False

            # Identifica nome base normalizado sem o prefixo do alvo
            clean_base = src.name
            if target_id and clean_base.startswith(f"{target_id}_"):
                clean_base = clean_base[len(target_id) + 1 :]

            mirror_key = (dest.name, clean_base, src.stat().st_size)
            if mirror_key in seen_mirror_keys:
                return False

            shutil.copy2(src, dest_file)
            seen_dest_files.add(dest_file.name)
            seen_mirror_keys.add(mirror_key)
            size_b = src.stat().st_size
            total_bytes += size_b
            size_mb = size_b / (1024 * 1024)
            if size_mb >= 1.0:
                print(f"  [OK] ({label:14}) {src.name:<34} ({size_mb:.2f} MB)")
            else:
                size_kb = size_b / 1024
                print(f"  [OK] ({label:14}) {src.name:<34} ({size_kb:.1f} KB)")
            copied_count += 1
            return True
        return False

    def copy_dir_tree(src_dir: Path, dest_parent: Path, label: str) -> bool:
        nonlocal copied_count, total_bytes
        if src_dir.exists() and src_dir.is_dir():
            dest = dest_parent / src_dir.name
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(src_dir, dest)
            dir_size = sum(f.stat().st_size for f in dest.rglob("*") if f.is_file())
            total_bytes += dir_size
            print(
                f"  [OK] ({label:14}) {src_dir.name}/ ({dir_size / (1024 * 1024):.2f} MB)"
            )
            copied_count += 1
            return True
        return False

    # 1. Docking, Estruturas Iniciais, Receptores e Ligantes
    print("━" * 68)
    print("1. Arquivando Docking, Estruturas e Interações (PLIP/ADMET)...")
    print("━" * 68)
    for lig in (
        list(data_dir.glob("*.sdf"))
        + list(data_dir.glob("*.pdbqt"))
        + list(data_dir.glob("*.mol2"))
    ):
        copy_file(lig, dest_dirs["docking"], "Ligand-Input")

    screening_base = data_dir / "screening"
    if screening_base.exists():
        for s_item in screening_base.rglob("*"):
            if s_item.is_file() and (
                target_id.lower() in s_item.name.lower()
                or target_id.lower() in str(s_item).lower()
            ):
                rel_p = s_item.relative_to(screening_base).parent
                target_sub = dest_dirs["docking"] / "screening" / rel_p
                target_sub.mkdir(parents=True, exist_ok=True)
                copy_file(s_item, target_sub, "Screening-Doc")

    for rec_dir in data_dir.glob("*"):
        if rec_dir.is_dir() and rec_dir.name not in [
            "md_files",
            "screening",
            "results",
            "archive",
        ]:
            if rec_dir.name.lower() == target_id.lower() or target_id == "DEFAULT":
                for f in rec_dir.rglob("*"):
                    if f.is_file() and f.suffix in [
                        ".pdb",
                        ".pdbqt",
                        ".sdf",
                        ".xml",
                        ".json",
                        ".txt",
                    ]:
                        rel_sub = f.relative_to(rec_dir).parent
                        target_sub = dest_dirs["docking"] / rec_dir.name / rel_sub
                        target_sub.mkdir(parents=True, exist_ok=True)
                        copy_file(f, target_sub, f"Rec-{rec_dir.name}")

    # 2. Topologias e Parâmetros de Campo de Força
    print("\n" + "━" * 68)
    print("2. Arquivando Topologias (GROMACS, ACPYPE e AMBER)...")
    print("━" * 68)
    top_candidates = [
        "topol.top",
        f"{target_id}_topol.top",
        "index.ndx",
        f"{target_id}_index.ndx",
    ]
    for tc in top_candidates:
        copy_file(selected_md_dir / tc, dest_dirs["topology"], "Topology/Index")

    for itp in sorted(selected_md_dir.glob("*.itp")):
        copy_file(itp, dest_dirs["topology"], "ITP-Forcefield")

    for prm in (
        list(selected_md_dir.glob("*.prmtop"))
        + list(selected_md_dir.glob("*.inpcrd"))
        + list(selected_md_dir.glob("*.frcmod"))
    ):
        copy_file(prm, dest_dirs["topology"], "AMBER-Parm")

    for acpype_dir in sorted(selected_md_dir.glob("*.acpype")):
        copy_dir_tree(acpype_dir, dest_dirs["topology"], "ACPYPE-Params")

    # 3. Trajetória e Simulação de Produção (GROMACS)
    print("\n" + "━" * 68)
    print("3. Arquivando Simulação, Trajetórias e Coordenadas...")
    print("━" * 68)
    traj_candidates = [
        # 1. Arquivo Binário de Parâmetros e Execução (TPR)
        "md.tpr",
        f"{target_id}_md.tpr",
        # 2. Coordenadas Finais e Estruturas Limpas de Equilíbrio
        "md.gro",
        f"{target_id}_md.gro",
        "md_clean.gro",
        f"{target_id}_md_clean.gro",
        "md_clean_nowat.pdb",
        f"{target_id}_md_clean_nowat.pdb",
        "cluster_medoid.gro",
        f"{target_id}_cluster_medoid.gro",
        "cluster_medoid.pdb",
        f"{target_id}_cluster_medoid.pdb",
        "complex.gro",
        f"{target_id}_complex.gro",
        # 3. Trajetória Bruta Integral (Raw Source of Truth - 100% dos átomos)
        "md.xtc",
        f"{target_id}_md.xtc",
        # 4. Trajetória Compacta para Visualização 3D no PyMOL (Proteína + Ligante sem Solvente)
        "md_fit_nowat.xtc",
        f"{target_id}_md_fit_nowat.xtc",
        # 5. Energias, Logs e Checkpoints para Continuação
        "md.edr",
        f"{target_id}_md.edr",
        "md.log",
        f"{target_id}_md.log",
        "md.cpt",
        f"{target_id}_md.cpt",
    ]
    for tc in traj_candidates:
        copy_file(selected_md_dir / tc, dest_dirs["trajectory"], "Traj/Coords")

    # 4. Relatórios Executivos, Gráficos, Matrizes CSV e Resultados MM-PBSA
    print("\n" + "━" * 68)
    print("4. Arquivando Relatórios, Gráficos, Matrizes CSV e Termodinâmica...")
    print("━" * 68)
    report_candidates = [
        "report.html",
        f"{target_id}_report.html",
        "show_complex.pml",
        f"{target_id}_show_complex.pml",
        "FINAL_RESULTS_MMPBSA.dat",
        f"{target_id}_FINAL_RESULTS_MMPBSA.dat",
        "FINAL_DECOMP_MMPBSA.dat",
        f"{target_id}_FINAL_DECOMP_MMPBSA.dat",
        "FINAL_RESULTS_MMPBSA.csv",
        f"{target_id}_FINAL_RESULTS_MMPBSA.csv",
        "RESULTS_gmx_MMPBSA.h5",
        f"{target_id}_RESULTS_gmx_MMPBSA.h5",
        "mmpbsa_summary.json",
        f"{target_id}_mmpbsa_summary.json",
        "interactions.json",
        f"{target_id}_interactions.json",
        "hbond_occupancy.json",
        f"{target_id}_hbond_occupancy.json",
        "pharmacokinetics.json",
        f"{target_id}_pharmacokinetics.json",
        "mmpbsa.in",
        f"{target_id}_mmpbsa.in",
        "gmx_MMPBSA.log",
        f"{target_id}_gmx_MMPBSA.log",
        "cluster.log",
        f"{target_id}_cluster.log",
        "rmsd-clust.xpm",
        f"{target_id}_rmsd-clust.xpm",
    ]
    for rc in report_candidates:
        copy_file(selected_md_dir / rc, dest_dirs["report"], "Report/Data")

    for ext_pattern in ["*.png", "*.csv", "*.xvg", "*.json"]:
        for f in sorted(selected_md_dir.glob(ext_pattern)):
            copy_file(f, dest_dirs["report"], "Analysis-Plot")

    # 5. Pacote de Exportação de Cluster (se existir)
    cluster_export_dir = root_dir / "cluster_export" / target_id
    if cluster_export_dir.exists():
        copy_dir_tree(cluster_export_dir, archive_base, "Cluster-Export")

    total_mb = total_bytes / (1024 * 1024)
    total_gb = total_mb / 1024

    print("\n" + "=" * 68)
    print("✨ Arquivamento Concluído com Sucesso!")
    print(f"Total de itens preservados: {copied_count}")
    if total_gb >= 1.0:
        print(f"Tamanho total do arquivo:   {total_gb:.2f} GB")
    else:
        print(f"Tamanho total do arquivo:   {total_mb:.2f} MB")
    print(f"Pasta de destino:           {archive_base.resolve()}")
    print("=" * 68)
    return archive_base


if __name__ == "__main__":
    exp_arg = sys.argv[1] if len(sys.argv) > 1 else None
    archive_experiment(exp_arg)
