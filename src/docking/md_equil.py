import subprocess
from pathlib import Path
from docking.md_prep import find_executable, DependencyError, SimulationPrepError

def run_md_equilibration(md_dir: Path):
    """
    Executa a etapa de Equilíbrio Termodinâmico (NVT e NPT) do sistema no GROMACS.
    Retorna um gerador yielding (step_code, status).
    """
    md_dir = Path(md_dir)
    if not md_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {md_dir}")

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError("O executável 'gmx' (GROMACS) não foi encontrado no PATH.")

    def run_command(cmd, cwd, input_val=None, step_name=""):
        try:
            exec_name = cmd[0]
            exec_path = find_executable(exec_name)
            if not exec_path:
                raise DependencyError(f"O executável '{exec_name}' não foi encontrado no PATH.")
            cmd[0] = exec_path
            
            import os
            env = os.environ.copy()
            exec_dir = str(Path(exec_path).parent)
            env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

            if input_val is not None and isinstance(input_val, bytes):
                input_val = input_val.decode('utf-8')

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

    # Encontrar caminhos para NVT e NPT mdp dinamicamente
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

    # Etapa A: Geração do Índice (make_ndx)
    yield "A", "start"
    cmd_make_ndx = [
        gmx_bin, "make_ndx",
        "-f", "em.gro",
        "-o", "index.ndx"
    ]
    # Executa apenas para salvar os grupos padrão em index.ndx
    run_command(cmd_make_ndx, md_dir, input_val="q\n", step_name="Etapa A (Geração do Índice - make_ndx)")
    
    # Processa e anexa os grupos Protein_LIG e Water_and_ions programaticamente no arquivo index.ndx
    try:
        index_path = md_dir / "index.ndx"
        if not index_path.exists():
            raise FileNotFoundError(f"Arquivo index.ndx não foi gerado em {md_dir}")
            
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
            raise ValueError("Grupo do ligante ('ligand_md' ou 'LIG') não encontrado no index.ndx padrão.")
            
        protein_lig_atoms = protein_atoms + ligand_atoms
        
        # 2. Criação do grupo Water_and_ions (SOL/Water + Ions/NA/CL)
        sol_atoms = groups.get("SOL", [])
        if not sol_atoms:
            sol_atoms = groups.get("Water", [])
            
        ions_atoms = groups.get("Ions", [])
        if not ions_atoms:
            # Tenta combinar os grupos individuais de íons caso Ions não esteja presente
            ions_atoms = groups.get("NA", []) + groups.get("CL", [])
            
        water_ions_atoms = sol_atoms + ions_atoms
        
        def format_group(name, atoms):
            lines = [f"[ {name} ]\n"]
            for i in range(0, len(atoms), 15):
                lines.append(" ".join(atoms[i:i+15]) + "\n")
            return "".join(lines)
            
        with open(index_path, "r", encoding="utf-8") as f:
            content = f.read()
            
        if not content.endswith("\n"):
            content += "\n"
            
        content += "\n" + format_group("Protein_LIG", protein_lig_atoms)
        content += "\n" + format_group("Water_and_ions", water_ions_atoms)
        
        with open(index_path, "w", encoding="utf-8") as f:
            f.write(content)
            
    except Exception as e:
        raise SimulationPrepError(f"Erro ao atualizar o arquivo index.ndx: {e}")
        
    yield "A", "success"

    # Etapa B: Compilação NVT (grompp)
    yield "B", "start"
    cmd_grompp_nvt = [
        gmx_bin, "grompp",
        "-f", str(nvt_mdp),
        "-c", "em.gro",
        "-r", "em.gro",
        "-p", "topol.top",
        "-n", "index.ndx",
        "-o", "nvt.tpr"
    ]
    run_command(cmd_grompp_nvt, md_dir, step_name="Etapa B (Compilação NVT)")
    yield "B", "success"

    # Etapa C: Execução NVT (mdrun)
    yield "C", "start"
    cmd_mdrun_nvt = [
        gmx_bin, "run" if "gmx" not in gmx_bin else "mdrun",
        "-v",
        "-deffnm", "nvt"
    ]
    # We should ensure we call the correct mdrun
    cmd_mdrun_nvt[0] = "mdrun" # find_executable handles cmd[0] resolution inside run_command anyway.
    cmd_mdrun_nvt = [gmx_bin, "mdrun", "-v", "-deffnm", "nvt"]
    run_command(cmd_mdrun_nvt, md_dir, step_name="Etapa C (Execução NVT)")
    yield "C", "success"

    # Etapa D: Compilação NPT (grompp)
    yield "D", "start"
    cmd_grompp_npt = [
        gmx_bin, "grompp",
        "-f", str(npt_mdp),
        "-c", "nvt.gro",
        "-r", "nvt.gro",
        "-p", "topol.top",
        "-n", "index.ndx",
        "-o", "npt.tpr"
    ]
    run_command(cmd_grompp_npt, md_dir, step_name="Etapa D (Compilação NPT)")
    yield "D", "success"

    # Etapa E: Execução NPT (mdrun)
    yield "E", "start"
    cmd_mdrun_npt = [gmx_bin, "mdrun", "-v", "-deffnm", "npt"]
    run_command(cmd_mdrun_npt, md_dir, step_name="Etapa E (Execução NPT)")
    yield "E", "success"
