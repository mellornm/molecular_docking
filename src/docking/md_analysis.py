import subprocess
import sys
from pathlib import Path
from docking.md_prep import find_executable, DependencyError, SimulationPrepError

def run_production_md(working_dir: Path):
    """
    Compila e executa a etapa de Produção de Dinâmica Molecular no GROMACS.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    npt_gro = working_dir / "npt.gro"
    topol_top = working_dir / "topol.top"

    if not npt_gro.exists():
        raise FileNotFoundError(f"Arquivo de equilíbrio 'npt.gro' não encontrado no diretório: {working_dir}")
    if not topol_top.exists():
        raise FileNotFoundError(f"Arquivo de topologia 'topol.top' não encontrado no diretório: {working_dir}")

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError("O executável 'gmx' (GROMACS) não foi encontrado no PATH.")

    # Resolução do template md.mdp
    project_root = Path(__file__).resolve().parent.parent.parent
    md_mdp = project_root / "src" / "templates" / "mdp" / "md.mdp"
    if not md_mdp.exists():
        md_mdp = Path("src/templates/mdp/md.mdp").resolve()
        if not md_mdp.exists():
            raise FileNotFoundError("Arquivo template md.mdp não encontrado.")

    # 1. Compilação do arquivo de produção (grompp)
    cmd_grompp = [
        gmx_bin, "grompp",
        "-f", str(md_mdp),
        "-c", "npt.gro",
        "-t", "npt.cpt",
        "-p", "topol.top",
        "-n", "index.ndx",
        "-o", "md.tpr"
    ]

    def run_cmd(cmd, cwd, input_val=None, step_name=""):
        try:
            import os
            env = os.environ.copy()
            exec_dir = str(Path(gmx_bin).parent)
            env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

            result = subprocess.run(
                cmd,
                cwd=str(cwd),
                env=env,
                capture_output=True,
                text=True,
                input=input_val
            )
            if result.returncode != 0:
                error_msg = result.stderr.strip() or result.stdout.strip()
                raise SimulationPrepError(
                    f"Erro na {step_name}:\n"
                    f"Comando: {' '.join(cmd)}\n"
                    f"Código de retorno: {result.returncode}\n"
                    f"Erro real: {error_msg}"
                )
            return result
        except Exception as e:
            if isinstance(e, (SimulationPrepError, DependencyError)):
                raise e
            raise SimulationPrepError(f"Falha ao executar o comando da {step_name}: {e}")

    run_cmd(cmd_grompp, working_dir, step_name="Compilação de Produção (grompp)")

    # 2. Execução da simulação de produção (mdrun)
    cmd_mdrun = [gmx_bin, "mdrun", "-v", "-deffnm", "md"]
    
    try:
        import os
        env = os.environ.copy()
        exec_dir = str(Path(gmx_bin).parent)
        env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

        process = subprocess.Popen(
            cmd_mdrun,
            cwd=str(working_dir),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        # Captura e imprime o output em tempo real
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            
        process.wait()
        if process.returncode != 0:
            raise SimulationPrepError(
                f"Erro na Execução de Produção (mdrun):\n"
                f"Comando: {' '.join(cmd_mdrun)}\n"
                f"Código de retorno: {process.returncode}"
            )
    except Exception as e:
        if isinstance(e, (SimulationPrepError, DependencyError)):
            raise e
        raise SimulationPrepError(f"Falha ao executar a Produção de Dinâmica Molecular (mdrun) em tempo real: {e}")


def analyze_trajectory(working_dir: Path):
    """
    Executa a análise da trajetória da Dinâmica Molecular no GROMACS.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError("O executável 'gmx' (GROMACS) não foi encontrado no PATH.")

    def run_analysis_cmd(cmd, cwd, input_val, step_name=""):
        try:
            import os
            env = os.environ.copy()
            exec_dir = str(Path(gmx_bin).parent)
            env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

            result = subprocess.run(
                cmd,
                cwd=str(cwd),
                env=env,
                capture_output=True,
                text=True,
                input=input_val
            )
            if result.returncode != 0:
                error_msg = result.stderr.strip() or result.stdout.strip()
                raise SimulationPrepError(
                    f"Erro na análise ({step_name}):\n"
                    f"Comando: {' '.join(cmd)}\n"
                    f"Código de retorno: {result.returncode}\n"
                    f"Erro real: {error_msg}"
                )
            return result
        except Exception as e:
            if isinstance(e, (SimulationPrepError, DependencyError)):
                raise e
            raise SimulationPrepError(f"Falha ao executar comando de análise ({step_name}): {e}")

    # 1. RMSD do esqueleto da proteína
    cmd_rmsd = [gmx_bin, "rms", "-s", "md.tpr", "-f", "md.xtc", "-o", "rmsd.xvg", "-tu", "ns"]
    run_analysis_cmd(cmd_rmsd, working_dir, "4 4\n", "RMSD do esqueleto da proteína")

    # 2. RMSF por resíduo
    cmd_rmsf = [gmx_bin, "rmsf", "-s", "md.tpr", "-f", "md.xtc", "-o", "rmsf.xvg", "-res"]
    run_analysis_cmd(cmd_rmsf, working_dir, "3\n", "RMSF por resíduo")

    # 3. Pontes de hidrogênio entre Proteína e Ligante
    cmd_hbond = [gmx_bin, "hbond", "-s", "md.tpr", "-f", "md.xtc", "-num", "hbond.xvg"]
    run_analysis_cmd(cmd_hbond, working_dir, "Protein ligand_md\n", "Pontes de hidrogênio (Proteína-Ligante)")
