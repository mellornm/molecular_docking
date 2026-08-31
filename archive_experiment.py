#!/usr/bin/env python3
"""
Script de arquivamento otimizado de experimentos de Docking e Dinâmica Molecular.
Preserva todos os dados essenciais para reanálises completas e reprodutibilidade,
eliminando arquivos temporários e intermediários redundantes (md_nojump.xtc, md.xtc, step*.pdb).
"""

import shutil
from pathlib import Path


def archive_experiment(exp_name: str = None) -> Path:
    print("=" * 68)
    print("📦 ARQUIVAMENTO OTIMIZADO DE DOCKING & DINÂMICA MOLECULAR")
    print("=" * 68)

    root_dir = Path(__file__).resolve().parent
    data_dir = root_dir / "data"

    default_name = "DS-7CFN"
    if not exp_name:
        try:
            user_input = input(
                f"Nome do experimento para arquivamento [{default_name}]: "
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

    print(f"\n📁 Diretório de Destino: {archive_base}\n")
    copied_count = 0
    total_bytes = 0

    def copy_file(src: Path, dest: Path, label: str) -> bool:
        nonlocal copied_count, total_bytes
        if src.exists() and src.is_file():
            dest_file = dest / src.name
            shutil.copy2(src, dest_file)
            size_b = src.stat().st_size
            total_bytes += size_b
            size_mb = size_b / (1024 * 1024)
            if size_mb >= 1.0:
                print(f"  [OK] ({label:14}) {src.name:<32} ({size_mb:.2f} MB)")
            else:
                size_kb = size_b / 1024
                print(f"  [OK] ({label:14}) {src.name:<32} ({size_kb:.1f} KB)")
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

    md_dir = data_dir / "md_files"

    # 1. Docking, Estruturas Iniciais, Receptores e Ligantes
    print("━" * 68)
    print("1. Arquivando Docking, Estruturas e Interações (PLIP/ADMET)...")
    print("━" * 68)
    # Ligantes na raiz de data/
    for lig in list(data_dir.glob("*.sdf")) + list(data_dir.glob("*.pdbqt")):
        copy_file(lig, dest_dirs["docking"], "Ligand-Input")

    # Pastas de Receptores (ex: data/7CFN, data/1OSV)
    for rec_dir in data_dir.glob("*"):
        if rec_dir.is_dir() and rec_dir.name not in [
            "md_files",
            "screening",
            "results",
        ]:
            for f in rec_dir.rglob("*"):
                if f.is_file() and f.suffix in [
                    ".pdb",
                    ".pdbqt",
                    ".sdf",
                    ".xml",
                    ".json",
                    ".txt",
                ]:
                    # Preserva identificador da subpasta relativa
                    rel_sub = f.relative_to(rec_dir).parent
                    target_sub = dest_dirs["docking"] / rec_dir.name / rel_sub
                    target_sub.mkdir(parents=True, exist_ok=True)
                    copy_file(f, target_sub, f"Rec-{rec_dir.name}")

    # 2. Topologias e Parâmetros de Campo de Força
    print("\n" + "━" * 68)
    print("2. Arquivando Topologias (GROMACS, ACPYPE e AMBER)...")
    print("━" * 68)
    copy_file(md_dir / "topol.top", dest_dirs["topology"], "Topology-TOP")
    copy_file(md_dir / "index.ndx", dest_dirs["topology"], "Index-NDX")

    # ITPs do sistema e posre
    for itp in sorted(md_dir.glob("*.itp")):
        copy_file(itp, dest_dirs["topology"], "ITP-Forcefield")

    # Topologias AMBER geradas para o MM-PBSA
    for prm in (
        list(md_dir.glob("*.prmtop"))
        + list(md_dir.glob("*.inpcrd"))
        + list(md_dir.glob("*.frcmod"))
    ):
        copy_file(prm, dest_dirs["topology"], "AMBER-Parm")

    # Diretório ACPYPE completo (parâmetros GAFF2 do ligante)
    for acpype_dir in sorted(md_dir.glob("*.acpype")):
        copy_dir_tree(acpype_dir, dest_dirs["topology"], "ACPYPE-Params")

    # 3. Trajetória e Simulação de Produção (GROMACS)
    print("\n" + "━" * 68)
    print("3. Arquivando Simulação e Trajetória Bruta Original (md.xtc)...")
    print("━" * 68)
    copy_file(md_dir / "md.tpr", dest_dirs["trajectory"], "TPR-Binary")
    copy_file(md_dir / "md.gro", dest_dirs["trajectory"], "Coords-Final")
    copy_file(md_dir / "md.edr", dest_dirs["trajectory"], "Energies-EDR")
    copy_file(md_dir / "md.log", dest_dirs["trajectory"], "GMX-Log")
    copy_file(md_dir / "md.cpt", dest_dirs["trajectory"], "Checkpoint")

    # Arquiva a trajetória bruta original (md.xtc) e coordenadas de referência/equilíbrio
    copy_file(md_dir / "md.xtc", dest_dirs["trajectory"], "Traj-Raw")
    copy_file(md_dir / "md_clean.gro", dest_dirs["trajectory"], "Structure-Clean")
    copy_file(md_dir / "cluster_medoid.gro", dest_dirs["trajectory"], "Cluster-Medoid")

    # 4. Relatórios Executivos, Gráficos, Matrizes CSV e Resultados MM-PBSA
    print("\n" + "━" * 68)
    print("4. Arquivando Relatórios, Gráficos, Matrizes CSV e Termodinâmica...")
    print("━" * 68)
    analysis_files = (
        [
            md_dir / "report.html",
            md_dir / "show_complex.pml",
            md_dir / "FINAL_RESULTS_MMPBSA.dat",
            md_dir / "FINAL_DECOMP_MMPBSA.dat",
            md_dir / "FINAL_RESULTS_MMPBSA.csv",
            md_dir / "RESULTS_gmx_MMPBSA.h5",
            md_dir / "mmpbsa_summary.json",
            md_dir / "hbond_occupancy.json",
            md_dir / "mmpbsa.in",
            md_dir / "gmx_MMPBSA.log",
            md_dir / "cluster.log",
        ]
        + sorted(md_dir.glob("*.png"))
        + sorted(md_dir.glob("*.xvg"))
        + sorted(md_dir.glob("*.csv"))
        + sorted(md_dir.glob("*.json"))
    )

    # Elimina duplicatas mantendo a ordem
    seen_paths = set()
    unique_analysis_files = []
    for f in analysis_files:
        if f not in seen_paths:
            seen_paths.add(f)
            unique_analysis_files.append(f)

    for f in unique_analysis_files:
        copy_file(f, dest_dirs["report"], "Report/Data")

    # Relatórios adicionais em screening/ se houver
    for s_dir in data_dir.glob("screening/*"):
        if s_dir.is_dir():
            for p in (
                list(s_dir.glob("*.html"))
                + list(s_dir.glob("*.pml"))
                + list(s_dir.glob("*.json"))
            ):
                copy_file(p, dest_dirs["report"], "Screening")

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
    import sys

    exp_arg = sys.argv[1] if len(sys.argv) > 1 else None
    archive_experiment(exp_arg)
