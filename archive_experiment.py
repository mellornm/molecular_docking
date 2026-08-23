#!/usr/bin/env python3
import glob
import shutil
from pathlib import Path


def archive_experiment():
    print("=" * 60)
    print("ARQUIVAMENTO DE DOCKING & DINÂMICA MOLECULAR")
    print("=" * 60)

    # Solicita o nome do experimento
    default_name = "desoxicolato"
    exp_name = input(f"Nome do experimento/ligante [{default_name}]: ").strip()
    if not exp_name:
        exp_name = default_name

    root_dir = Path.cwd()
    data_dir = root_dir / "data"
    archive_base = root_dir / "archive" / exp_name

    # Estrutura de destino
    dest_dirs = {
        "report": archive_base / "report",
        "trajectory": archive_base / "trajectory",
        "topology": archive_base / "topology",
        "docking": archive_base / "docking",
    }

    for path in dest_dirs.values():
        path.mkdir(parents=True, exist_ok=True)

    print(f"\nDestino: {archive_base}\n")

    copied_count = 0

    def copy_file(src_path: Path, dest_dir: Path, label: str):
        nonlocal copied_count
        if src_path.exists() and src_path.is_file():
            shutil.copy2(src_path, dest_dir / src_path.name)
            size_mb = src_path.stat().st_size / (1024 * 1024)
            print(f" [OK] ({label}) {src_path.name} -> {size_mb:.2f} MB")
            copied_count += 1
        return

    def copy_pattern(pattern: str, dest_dir: Path, label: str):
        for file_str in glob.glob(pattern):
            copy_file(Path(file_str), dest_dir, label)

    # 1. Relatórios, Scripts e Dados Termodinâmicos
    screening_dir = data_dir / "screening" / exp_name
    md_dir = data_dir / "md_files"

    copy_file(screening_dir / "report.html", dest_dirs["report"], "Report")
    copy_file(screening_dir / "show_complex.pml", dest_dirs["report"], "PyMOL")
    copy_pattern(str(screening_dir / "*.png"), dest_dirs["report"], "Plot")
    copy_pattern(str(screening_dir / "*.json"), dest_dirs["report"], "JSON")
    copy_pattern(str(screening_dir / "*.dat"), dest_dirs["report"], "Data")

    copy_pattern(str(md_dir / "*.xvg"), dest_dirs["report"], "Plot-XVG")
    copy_pattern(str(md_dir / "*.png"), dest_dirs["report"], "Plot-PNG")
    copy_file(md_dir / "mmpbsa_summary.json", dest_dirs["report"], "MMPBSA-JSON")
    copy_file(md_dir / "FINAL_RESULTS_MMPBSA.dat", dest_dirs["report"], "MMPBSA-DAT")

    # 2. Trajetórias e Simulação de Produção (GROMACS)
    copy_file(md_dir / "md_fit.xtc", dest_dirs["trajectory"], "Traj-Fit")
    copy_file(md_dir / "md.xtc", dest_dirs["trajectory"], "Traj-Raw")
    copy_file(md_dir / "md.gro", dest_dirs["trajectory"], "Structure")
    copy_file(md_dir / "md.tpr", dest_dirs["trajectory"], "TPR-Bin")
    copy_file(md_dir / "md.log", dest_dirs["trajectory"], "Log")
    copy_file(md_dir / "md.edr", dest_dirs["trajectory"], "Energy-EDR")

    # 3. Topologia e Parâmetros
    copy_file(md_dir / "topol.top", dest_dirs["topology"], "Topology")
    copy_file(md_dir / "index.ndx", dest_dirs["topology"], "Index")
    copy_pattern(str(md_dir / "*.itp"), dest_dirs["topology"], "ITP-Forcefield")

    # 4. Docking, Ligante e Receptor
    # Pega ligantes na raiz de data/
    copy_pattern(str(data_dir / f"{exp_name}.*"), dest_dirs["docking"], "Ligand")
    copy_pattern(
        str(data_dir / "results" / f"*{exp_name}*"),
        dest_dirs["docking"],
        "Docking-Result",
    )
    copy_pattern(
        str(data_dir / "results" / "*.pdb*"), dest_dirs["docking"], "Docking-Pose"
    )

    # Pega arquivos de receptor (ex: data/1OSV/* ou data/*.pdb)
    for receptor_folder in data_dir.glob("*"):
        if receptor_folder.is_dir() and receptor_folder.name not in [
            "md_files",
            "results",
            "screening",
        ]:
            copy_pattern(
                str(receptor_folder / "*.pdb*"), dest_dirs["docking"], "Receptor"
            )
            copy_pattern(
                str(receptor_folder / "*.sdf"), dest_dirs["docking"], "Receptor-SDF"
            )

    print("-" * 60)
    print("✨ Arquivamento concluído com sucesso!")
    print(f"Total de arquivos preservados: {copied_count}")
    print(f"Pasta de destino: {archive_base.resolve()}")
    print("-" * 60)


if __name__ == "__main__":
    archive_experiment()
