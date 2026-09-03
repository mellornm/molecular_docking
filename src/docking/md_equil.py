import os
import shutil
import subprocess
from pathlib import Path
from typing import Generator, Optional, Tuple

from docking.md_prep import (
    DependencyError,
    SimulationPrepError,
    find_executable,
    sanitize_target_id,
    verify_tpr_consistency,
)


def run_md_equilibration(
    md_dir: Path, target_id: Optional[str] = None
) -> Generator[Tuple[str, str], None, None]:
    """
    Executa o Equilíbrio Termodinâmico (NVT e NPT) no GROMACS com isolamento estrito por alvo (Target Isolation),
    controle estrito de erro no grompp e checagem de consistência pós-geração do TPR.

    Retorna um gerador (etapa, status).
    """
    md_dir = Path(md_dir)
    if not md_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {md_dir}")

    # Identificação do alvo
    if not target_id:
        target_id = sanitize_target_id(md_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            # Tenta encontrar arquivos com prefixo no diretório
            gro_files = list(md_dir.glob("*_em.gro"))
            if gro_files:
                target_id = gro_files[0].stem.replace("_em", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError(
            "O executável 'gmx' (GROMACS) não foi encontrado no PATH."
        )

    def run_command(cmd, cwd, input_val=None, step_name="", expect_zero_only=True):
        try:
            exec_name = cmd[0]
            exec_path = find_executable(exec_name)
            if not exec_path:
                raise DependencyError(
                    f"O executável '{exec_name}' não foi encontrado no PATH."
                )
            cmd[0] = exec_path

            env = os.environ.copy()
            env.pop("PYTHONPATH", None)
            env.pop("PYTHONHOME", None)
            exec_dir = str(Path(exec_path).parent)
            env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

            if input_val is not None and isinstance(input_val, bytes):
                input_val = input_val.decode("utf-8")

            result = subprocess.run(
                cmd,
                cwd=str(cwd),
                env=env,
                capture_output=True,
                text=True,
                input=input_val,
            )
            if expect_zero_only and result.returncode != 0:
                error_msg = result.stderr.strip() or result.stdout.strip()
                raise SimulationPrepError(
                    f"Erro estrito na {step_name}:\n"
                    f"Comando: {' '.join(cmd)}\n"
                    f"Código de retorno: {result.returncode}\n"
                    f"Erro do GROMACS: {error_msg}"
                )
            return result
        except Exception as e:
            if isinstance(e, (SimulationPrepError, DependencyError)):
                raise e
            raise SimulationPrepError(
                f"Falha ao executar o comando da {step_name}: {e}"
            )

    # Identificação dos arquivos mdp
    project_root = Path(__file__).resolve().parent.parent.parent
    nvt_mdp = project_root / "src" / "templates" / "mdp" / "nvt.mdp"
    npt_mdp = project_root / "src" / "templates" / "mdp" / "npt.mdp"

    if not nvt_mdp.exists():
        nvt_mdp = Path("src/templates/mdp/nvt.mdp").resolve()
        if not nvt_mdp.exists():
            raise FileNotFoundError("Arquivo template nvt.mdp não encontrado.")

    if not npt_mdp.exists():
        npt_mdp = Path("src/templates/mdp/npt.mdp").resolve()
        if not npt_mdp.exists():
            raise FileNotFoundError("Arquivo template npt.mdp não encontrado.")

    # Localiza o arquivo de coordenadas de entrada (em.gro)
    em_gro = md_dir / f"{prefix}em.gro"
    if not em_gro.exists():
        em_gro = md_dir / "em.gro"
    if not em_gro.exists():
        raise FileNotFoundError(
            f"Arquivo de coordenadas minimizadas '{em_gro.name}' não encontrado em {md_dir}. "
            "Execute a etapa de preparação/minimização primeiro."
        )

    topol_top = md_dir / f"{prefix}topol.top"
    if not topol_top.exists():
        topol_top = md_dir / "topol.top"
    if not topol_top.exists():
        raise FileNotFoundError(f"Arquivo 'topol.top' não encontrado em {md_dir}")

    # Etapa A: Geração do Índice (make_ndx)
    yield "A", "start"
    index_name = f"{prefix}index.ndx"
    cmd_make_ndx = [gmx_bin, "make_ndx", "-f", str(em_gro.name), "-o", index_name]
    run_command(
        cmd_make_ndx,
        md_dir,
        input_val="q\n",
        step_name="Etapa A (Geração do Índice - make_ndx)",
    )

    # Processa e anexa os grupos Protein_LIG e Water_and_ions programaticamente no arquivo index.ndx
    try:
        index_path = md_dir / index_name
        if not index_path.exists():
            raise FileNotFoundError(f"Arquivo {index_name} não foi gerado em {md_dir}")

        groups = {}
        current_group = None
        with open(index_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("[") and line.endswith("]"):
                    current_group = line[1:-1].strip()
                    groups[current_group] = []
                elif current_group is not None:
                    groups[current_group].extend(line.split())

        # 1. Criação do grupo Protein_LIG (Protein + ligand_md / LIG)
        protein_atoms = groups.get("Protein", [])
        ligand_atoms = groups.get("ligand_md", [])
        if not ligand_atoms:
            ligand_atoms = groups.get("LIG", [])

        if not protein_atoms:
            raise ValueError("Grupo 'Protein' não encontrado no index.ndx padrão.")
        if not ligand_atoms:
            raise ValueError(
                "Grupo do ligante ('ligand_md' ou 'LIG') não encontrado no index.ndx padrão."
            )

        protein_lig_atoms = protein_atoms + ligand_atoms

        # 2. Criação do grupo Water_and_ions (SOL/Water + Ions/NA/CL)
        sol_atoms = groups.get("SOL", [])
        if not sol_atoms:
            sol_atoms = groups.get("Water", [])

        ions_atoms = groups.get("Ions", [])
        if not ions_atoms:
            ions_atoms = groups.get("NA", []) + groups.get("CL", [])

        water_ions_atoms = sol_atoms + ions_atoms

        def format_group(name, atoms):
            lines = [f"[ {name} ]\n"]
            for i in range(0, len(atoms), 15):
                lines.append(" ".join(atoms[i : i + 15]) + "\n")
            return "".join(lines)

        with open(index_path, "r", encoding="utf-8") as f:
            content = f.read()

        if not content.endswith("\n"):
            content += "\n"

        content += "\n" + format_group("Protein_LIG", protein_lig_atoms)
        content += "\n" + format_group("Water_and_ions", water_ions_atoms)

        with open(index_path, "w", encoding="utf-8") as f:
            f.write(content)

        # Espelho index.ndx
        shutil.copy2(index_path, md_dir / "index.ndx")

    except Exception as e:
        if isinstance(e, SimulationPrepError):
            raise e
        raise SimulationPrepError(f"Erro ao atualizar o arquivo de índice: {e}")

    yield "A", "success"

    # Etapa B: Compilação NVT (grompp)
    yield "B", "start"
    nvt_tpr = md_dir / f"{prefix}nvt.tpr"
    nvt_tpr.unlink(missing_ok=True)
    (md_dir / "nvt.tpr").unlink(missing_ok=True)

    cmd_grompp_nvt = [
        gmx_bin,
        "grompp",
        "-f",
        str(nvt_mdp),
        "-c",
        str(em_gro.name),
        "-r",
        str(em_gro.name),
        "-p",
        str(topol_top.name),
        "-n",
        index_name,
        "-o",
        f"{prefix}nvt.tpr",
    ]
    run_command(cmd_grompp_nvt, md_dir, step_name="Etapa B (Compilação NVT)")
    shutil.copy2(nvt_tpr, md_dir / "nvt.tpr")
    verify_tpr_consistency(nvt_tpr)
    yield "B", "success"

    # Etapa C: Execução NVT (mdrun)
    yield "C", "start"
    cmd_mdrun_nvt = [gmx_bin, "mdrun", "-v", "-deffnm", f"{prefix}nvt"]
    run_command(cmd_mdrun_nvt, md_dir, step_name="Etapa C (Execução NVT)")
    if (md_dir / f"{prefix}nvt.gro").exists():
        shutil.copy2(md_dir / f"{prefix}nvt.gro", md_dir / "nvt.gro")
    yield "C", "success"

    # Etapa D: Compilação NPT (grompp)
    yield "D", "start"
    npt_tpr = md_dir / f"{prefix}npt.tpr"
    npt_tpr.unlink(missing_ok=True)
    (md_dir / "npt.tpr").unlink(missing_ok=True)

    nvt_gro = md_dir / f"{prefix}nvt.gro"
    if not nvt_gro.exists():
        nvt_gro = md_dir / "nvt.gro"

    cmd_grompp_npt = [
        gmx_bin,
        "grompp",
        "-f",
        str(npt_mdp),
        "-c",
        str(nvt_gro.name),
        "-r",
        str(nvt_gro.name),
        "-p",
        str(topol_top.name),
        "-n",
        index_name,
        "-o",
        f"{prefix}npt.tpr",
    ]
    run_command(cmd_grompp_npt, md_dir, step_name="Etapa D (Compilação NPT)")
    shutil.copy2(npt_tpr, md_dir / "npt.tpr")
    verify_tpr_consistency(npt_tpr)
    yield "D", "success"

    # Etapa E: Execução NPT (mdrun)
    yield "E", "start"
    cmd_mdrun_npt = [gmx_bin, "mdrun", "-v", "-deffnm", f"{prefix}npt"]
    run_command(cmd_mdrun_npt, md_dir, step_name="Etapa E (Execução NPT)")
    if (md_dir / f"{prefix}npt.gro").exists():
        shutil.copy2(md_dir / f"{prefix}npt.gro", md_dir / "npt.gro")
    yield "E", "success"
