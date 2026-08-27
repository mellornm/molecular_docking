import collections
import shutil
import subprocess
from pathlib import Path
from rdkit import Chem
from pdbfixer import PDBFixer
from openmm.app import PDBFile


class DependencyError(Exception):
    """Exceção levantada quando um binário externo (GROMACS ou ACPYPE) não é encontrado no PATH."""

    pass


class SimulationPrepError(Exception):
    """Exceção levantada quando há falha no processamento ou execução das etapas de preparação."""

    pass


def find_executable(name: str) -> str:
    """Procura pelo executável no PATH atual ou em locais comuns do conda 'bioinfo'."""
    path = (
        shutil.which(name) or shutil.which(f"{name}.py") or shutil.which(f"{name}.exe")
    )
    if path:
        return path

    home = Path.home()
    # Caminhos para verificar no Linux/WSL e macOS
    possible_paths = [
        home / "miniforge3" / "envs" / "bioinfo" / "bin" / name,
        home / "miniforge3" / "envs" / "bioinfo" / "bin" / f"{name}.py",
        home / "miniconda3" / "envs" / "bioinfo" / "bin" / name,
        home / "miniconda3" / "envs" / "bioinfo" / "bin" / f"{name}.py",
        Path("/home/rmello/miniforge3/envs/bioinfo/bin") / name,
        Path("/home/rmello/miniforge3/envs/bioinfo/bin") / f"{name}.py",
    ]

    # Caminhos para verificar no Windows
    possible_paths += [
        home / "Miniconda3" / "envs" / "bioinfo" / "Scripts" / name,
        home / "Miniconda3" / "envs" / "bioinfo" / "Scripts" / f"{name}.exe",
        home / "Miniforge3" / "envs" / "bioinfo" / "Scripts" / name,
        home / "Miniforge3" / "envs" / "bioinfo" / "Scripts" / f"{name}.exe",
    ]

    for p in possible_paths:
        if p.exists():
            return str(p)

    return None


def extract_ligand(ligand_sdf: Path, output_dir: Path) -> Path:
    """
    Etapa A: Usa o RDKit para ler a primeira pose de 'ligand_sdf'.
    Salva em 'ligand_md.pdb' dentro de 'output_dir', garantindo átomos com
    nomes únicos (C1, C2, O1...) e resíduo 'LIG' na cadeia 'X'.
    """
    try:
        # Importante: removeHs=False garante que os hidrogênios do SDF sejam lidos
        supplier = Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)
        if not supplier or len(supplier) == 0 or supplier[0] is None:
            raise ValueError(
                f"Não foi possível ler o arquivo SDF ou ele está vazio: {ligand_sdf}"
            )
        mol = supplier[0]
        # Adiciona hidrogênios para evitar problemas de radical ímpar (sqm/antechamber)
        mol = Chem.AddHs(mol, addCoords=True)
    except Exception as e:
        raise SimulationPrepError(f"Falha ao ler o arquivo SDF com RDKit: {e}")

    # Garante nomes de átomos únicos, resíduo 'LIG' e cadeia 'X'
    element_counts = collections.defaultdict(int)
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        element_counts[symbol] += 1
        atom_name = f"{symbol}{element_counts[symbol]}"

        # Formatação para o formato PDB (4 caracteres)
        if len(symbol) == 1:
            padded_name = f" {atom_name:<3}"
        else:
            padded_name = f"{atom_name:<4}"

        info = Chem.AtomPDBResidueInfo()
        info.SetName(padded_name)
        info.SetResidueName("LIG")
        info.SetChainId("X")
        info.SetResidueNumber(1)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)

    ligand_pdb_path = output_dir / "ligand_md.pdb"
    try:
        Chem.MolToPDBFile(mol, str(ligand_pdb_path))
    except Exception as e:
        raise SimulationPrepError(f"Falha ao exportar ligante para PDB: {e}")

    return ligand_pdb_path


def run_acpype(ligand_pdb: Path, output_dir: Path):
    """
    Etapa B: Executa o ACPYPE via subprocesso: 'acpype -i ligand_md.pdb -c bcc'.
    O diretório de trabalho é definido para output_dir para manter os outputs organizados.
    """
    acpype_bin = find_executable("acpype")
    if not acpype_bin:
        raise DependencyError(
            "O executável 'acpype' não foi encontrado no PATH ou no ambiente 'bioinfo'. "
            "Certifique-se de que ele está instalado e configurado."
        )

    # Garante que o caminho para o arquivo PDB do ligante seja relativo a output_dir se estiver nele,
    # ou use o nome do arquivo diretamente se rodando no output_dir.
    ligand_pdb_name = ligand_pdb.name

    # Copia o PATH e adiciona o diretório do binário do acpype (onde o obabel está localizado no Conda)
    import os

    env = os.environ.copy()
    acpype_dir = str(Path(acpype_bin).parent)
    env["PATH"] = f"{acpype_dir}{os.pathsep}{env.get('PATH', '')}"

    cmd_acpype = [acpype_bin, "-i", ligand_pdb_name, "-c", "bcc", "-f"]
    try:
        subprocess.run(
            cmd_acpype,
            cwd=str(output_dir),
            env=env,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        # Mostra stdout ou stderr dependendo de onde o ACPYPE printou o erro (geralmente stdout)
        detailed_error = e.stderr.strip() or e.stdout.strip()
        raise SimulationPrepError(
            f"Falha ao rodar ACPYPE (código {e.returncode}):\n"
            f"Comando: {' '.join(cmd_acpype)}\n"
            f"Erro: {detailed_error}"
        )
    except Exception as e:
        raise SimulationPrepError(f"Erro ao iniciar processo do ACPYPE: {e}")


def run_pdb2gmx(receptor_pdb: Path, output_dir: Path):
    """
    Etapa C: Executa o comando do GROMACS para preparar a proteína via subprocesso:
    'gmx pdb2gmx -f receptor_pdb -o protein_processed.gro -p topol.top -ff amber99sb-ildn -water tip3p'.
    Se falhar devido a átomos ausentes na estrutura (comum em arquivos PDB), tenta novamente com o flag '-missing'.
    """
    # Se o receptor recebido não for o arquivo já curado, realiza a cura estrutural primeiro
    if receptor_pdb.name != "receptor_fixed.pdb":
        try:
            fixer = PDBFixer(filename=str(receptor_pdb))
            fixer.removeHeterogens(keepWater=False)
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()

            receptor_fixed = output_dir / "receptor_fixed.pdb"
            with open(receptor_fixed, "w") as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)
            receptor_pdb = receptor_fixed
        except Exception as e:
            raise SimulationPrepError(f"Erro ao curar o receptor com PDBFixer: {e}")

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError(
            "O executável 'gmx' (GROMACS) não foi encontrado no PATH ou no ambiente 'bioinfo'. "
            "Certifique-se de que o GROMACS está instalado."
        )

    # Resolve o caminho do receptor para absoluto para garantir correto funcionamento em cwd diferente
    abs_receptor_pdb = receptor_pdb.resolve()
    if not abs_receptor_pdb.exists():
        raise FileNotFoundError(
            f"Arquivo do receptor não encontrado em: {abs_receptor_pdb}"
        )

    # Copia o PATH e adiciona o diretório do binário do gmx para o subprocesso
    import os

    env = os.environ.copy()
    gmx_dir = str(Path(gmx_bin).parent)
    env["PATH"] = f"{gmx_dir}{os.pathsep}{env.get('PATH', '')}"

    cmd_gmx = [
        gmx_bin,
        "pdb2gmx",
        "-f",
        str(abs_receptor_pdb),
        "-o",
        "protein_processed.gro",
        "-p",
        "topol.top",
        "-ff",
        "amber99sb-ildn",
        "-water",
        "tip3p",
    ]
    try:
        subprocess.run(
            cmd_gmx,
            cwd=str(output_dir),
            env=env,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        detailed_error = e.stderr.strip() or e.stdout.strip()
        # Se falhar por causa de átomos ausentes (ex: CG de um ASP ausente), tenta novamente com -missing
        if (
            "not found in the input file" in detailed_error
            or "missing" in detailed_error.lower()
        ):
            cmd_gmx_missing = cmd_gmx + ["-missing"]
            try:
                subprocess.run(
                    cmd_gmx_missing,
                    cwd=str(output_dir),
                    env=env,
                    check=True,
                    capture_output=True,
                    text=True,
                )
                return
            except subprocess.CalledProcessError as e_inner:
                detailed_error = e_inner.stderr.strip() or e_inner.stdout.strip()
                raise SimulationPrepError(
                    f"Falha ao rodar GROMACS pdb2gmx (mesmo com -missing, código {e_inner.returncode}):\n"
                    f"Comando: {' '.join(cmd_gmx_missing)}\n"
                    f"Erro: {detailed_error}"
                )
        raise SimulationPrepError(
            f"Falha ao rodar GROMACS pdb2gmx (código {e.returncode}):\n"
            f"Comando: {' '.join(cmd_gmx)}\n"
            f"Erro: {detailed_error}"
        )
    except Exception as e:
        raise SimulationPrepError(f"Erro ao iniciar processo do GROMACS: {e}")


def stitch_topology(output_dir: Path):
    """
    Etapa D: Lógica em Python para ler o arquivo 'topol.top' gerado e injetar o include do ligante
    e adicioná-lo na seção [ molecules ] abaixo da proteína.
    """
    topol_path = output_dir / "topol.top"
    if not topol_path.exists():
        raise SimulationPrepError(
            f"O arquivo 'topol.top' esperado não foi gerado em {output_dir}"
        )

    try:
        with open(topol_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
    except Exception as e:
        raise SimulationPrepError(f"Falha ao ler 'topol.top': {e}")

    # Insere o include do ligante após a inclusão de forcefield.itp
    new_lines = []
    ff_found = False
    for line in lines:
        new_lines.append(line)
        if "forcefield.itp" in line and not ff_found:
            new_lines.append(
                "\n; Inclui a topologia do ligante gerada pelo ACPYPE\n"
                '#include "ligand_md.acpype/ligand_md_GMX.itp"\n'
                "\n; Restrições de posição do ligante para NVT/NPT\n"
                "#ifdef POSRES\n"
                '#include "ligand_md.acpype/posre_ligand_md.itp"\n'
                "#endif\n"
            )
            ff_found = True

    if not ff_found:
        raise SimulationPrepError(
            "Inclusão de 'forcefield.itp' não encontrada em topol.top."
        )

    # Localiza a seção [ molecules ]
    molecules_idx = -1
    for i, line in enumerate(new_lines):
        if line.strip().startswith("[ molecules ]"):
            molecules_idx = i
            break

    if molecules_idx == -1:
        raise SimulationPrepError("Seção '[ molecules ]' não encontrada em topol.top.")

    # Localiza a primeira linha de molécula abaixo de [ molecules ]
    protein_idx = -1
    for i in range(molecules_idx + 1, len(new_lines)):
        line_strip = new_lines[i].strip()
        if (
            line_strip
            and not line_strip.startswith(";")
            and not line_strip.startswith("#")
        ):
            protein_idx = i
            break

    if protein_idx == -1:
        raise SimulationPrepError(
            "Nenhuma molécula (proteína) ativa foi encontrada sob a seção '[ molecules ]'."
        )

    # Garante que a última linha termina com quebra de linha antes de adicionar o ligante
    if new_lines and not new_lines[-1].endswith("\n"):
        new_lines[-1] = new_lines[-1] + "\n"

    # Insere a definição do ligante ligand_md estritamente na última linha
    new_lines.append("ligand_md                 1\n")

    try:
        with open(topol_path, "w", encoding="utf-8") as f:
            f.writelines(new_lines)
    except Exception as e:
        raise SimulationPrepError(
            f"Falha ao salvar as modificações em 'topol.top': {e}"
        )


def prepare_md_system(receptor_pdb: Path, ligand_sdf: Path, output_dir: Path):
    """
    Prepara o sistema completo de Dinâmica Molecular executando sequencialmente as etapas de A a L.
    Retorna um gerador que indica o início (start) e fim (success) de cada subetapa.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    gmx_bin = find_executable("gmx")
    acpype_bin = find_executable("acpype")

    def run_command(cmd, cwd, input_val=None, step_name=""):
        try:
            exec_name = cmd[0]
            exec_path = find_executable(exec_name)
            if not exec_path:
                raise DependencyError(
                    f"O executável '{exec_name}' não foi encontrado no PATH."
                )
            cmd[0] = exec_path

            import os

            env = os.environ.copy()
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
            raise SimulationPrepError(
                f"Falha ao executar o comando da {step_name}: {e}"
            )

    # Etapa A: Cura com PDBFixer
    yield "A", "start"
    try:
        fixer = PDBFixer(filename=str(receptor_pdb))
        fixer.removeHeterogens(keepWater=False)

        # Validação: verifica se há resíduos proteicos válidos
        standard_amino_acids = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
            "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN", "ARN"
        }
        protein_res = [
            r for r in fixer.topology.residues() if r.name.upper() in standard_amino_acids
        ]
        if not protein_res:
            raise SimulationPrepError(
                f"O arquivo fornecido como receptor ('{receptor_pdb}') não contém resíduos de proteína (aminoácidos padrão). "
                f"Certifique-se de selecionar o arquivo 'receptor.pdb' (ex: data/7CFN/processed/receptor.pdb) e NÃO o arquivo do ligante ('ligand.pdb')."
            )

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        receptor_fixed = output_dir / "receptor_fixed.pdb"
        with open(receptor_fixed, "w", encoding="utf-8") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
    except Exception as e:
        if isinstance(e, SimulationPrepError):
            raise e
        raise SimulationPrepError(f"Erro na Etapa A (Cura com PDBFixer): {e}")
    yield "A", "success"

    # Etapa B: Extração do Ligante
    yield "B", "start"
    try:
        supplier = Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)
        if not supplier or len(supplier) == 0 or supplier[0] is None:
            raise ValueError(
                f"Não foi possível ler o arquivo SDF ou ele está vazio: {ligand_sdf}"
            )
        mol = supplier[0]
        # Remove hidrogênios residuais cujas coordenadas possam ter sido corrompidas/deslocadas pelo docking
        mol = Chem.RemoveHs(mol)
        # Recalcula e adiciona todos os hidrogênios com posições geométricas 3D precisas
        mol = Chem.AddHs(mol, addCoords=True)

        element_counts = collections.defaultdict(int)
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            element_counts[symbol] += 1
            atom_name = f"{symbol}{element_counts[symbol]}"
            if len(symbol) == 1:
                padded_name = f" {atom_name:<3}"
            else:
                padded_name = f"{atom_name:<4}"

            info = Chem.AtomPDBResidueInfo()
            info.SetName(padded_name)
            info.SetResidueName("LIG")
            info.SetChainId("X")
            info.SetResidueNumber(1)
            info.SetIsHeteroAtom(True)
            atom.SetMonomerInfo(info)

        ligand_pdb = output_dir / "ligand_md.pdb"
        Chem.MolToPDBFile(mol, str(ligand_pdb))
    except Exception as e:
        raise SimulationPrepError(f"Erro na Etapa B (Extração do Ligante): {e}")
    yield "B", "success"

    # Etapa C: Parametrização
    yield "C", "start"
    if not acpype_bin:
        raise DependencyError("O executável 'acpype' não foi encontrado no PATH.")
    cmd_acpype = [acpype_bin, "-i", "ligand_md.pdb", "-c", "bcc", "-f"]
    run_command(cmd_acpype, output_dir, step_name="Etapa C (Parametrização)")
    yield "C", "success"

    # Etapa D: Topologia da Proteína
    yield "D", "start"
    if not gmx_bin:
        raise DependencyError(
            "O executável 'gmx' (GROMACS) não foi encontrado no PATH."
        )
    cmd_pdb2gmx = [
        gmx_bin,
        "pdb2gmx",
        "-f",
        "receptor_fixed.pdb",
        "-o",
        "protein_processed.gro",
        "-p",
        "topol.top",
        "-ff",
        "amber99sb-ildn",
        "-water",
        "tip3p",
    ]
    try:
        run_command(
            cmd_pdb2gmx, output_dir, step_name="Etapa D (Topologia da Proteína)"
        )
    except SimulationPrepError as e:
        if "not found in the input file" in str(e) or "missing" in str(e).lower():
            cmd_pdb2gmx_missing = cmd_pdb2gmx + ["-missing"]
            run_command(
                cmd_pdb2gmx_missing,
                output_dir,
                step_name="Etapa D (Topologia da Proteína com -missing)",
            )
        else:
            raise e
    yield "D", "success"

    # Etapa E: Fusão de Coordenadas
    yield "E", "start"
    try:
        prot_gro_path = output_dir / "protein_processed.gro"
        lig_gro_path = output_dir / "ligand_md.acpype" / "ligand_md_GMX.gro"

        if not prot_gro_path.exists():
            raise FileNotFoundError(f"Arquivo {prot_gro_path} não encontrado.")
        if not lig_gro_path.exists():
            raise FileNotFoundError(f"Arquivo {lig_gro_path} não encontrado.")

        with open(prot_gro_path, "r", encoding="utf-8") as f:
            prot_lines = f.readlines()
        with open(lig_gro_path, "r", encoding="utf-8") as f:
            lig_lines = f.readlines()

        if len(prot_lines) < 3:
            raise ValueError("protein_processed.gro possui menos de 3 linhas.")
        if len(lig_lines) < 3:
            raise ValueError("ligand_md_GMX.gro possui menos de 3 linhas.")

        prot_title = prot_lines[0]
        prot_atoms = prot_lines[2:-1]
        lig_atoms = lig_lines[2:-1]

        total_atoms = len(prot_atoms) + len(lig_atoms)
        box_vector = prot_lines[-1]

        complex_lines = ["Complex of Protein and Ligand\n", f" {total_atoms}\n"]
        complex_lines.extend(prot_atoms)
        complex_lines.extend(lig_atoms)
        complex_lines.append(box_vector)

        complex_gro_path = output_dir / "complex.gro"
        with open(complex_gro_path, "w", encoding="utf-8") as f:
            f.writelines(complex_lines)

    except Exception as e:
        raise SimulationPrepError(f"Erro na Etapa E (Fusão de Coordenadas): {e}")
    yield "E", "success"

    # Etapa F: Fusão de Topologia
    yield "F", "start"
    try:
        topol_path = output_dir / "topol.top"
        if not topol_path.exists():
            raise FileNotFoundError(f"Arquivo {topol_path} não encontrado.")

        with open(topol_path, "r", encoding="utf-8") as f:
            topol_lines = f.readlines()

        new_topol_lines = []
        ff_found = False
        for line in topol_lines:
            new_topol_lines.append(line)
            if "forcefield.itp" in line and not ff_found:
                new_topol_lines.append(
                    '\n; Inclui a topologia do ligante gerada pelo ACPYPE\n#include "ligand_md.acpype/ligand_md_GMX.itp"\n'
                )
                ff_found = True

        if not ff_found:
            raise ValueError(
                "Inclusão de 'forcefield.itp' não encontrada em topol.top."
            )

        # Garante que a última linha termina com quebra de linha antes de adicionar o ligante
        if new_topol_lines and not new_topol_lines[-1].endswith("\n"):
            new_topol_lines[-1] = new_topol_lines[-1] + "\n"

        # Append ligand molecule name and count strictly as the last line
        new_topol_lines.append("ligand_md             1\n")

        with open(topol_path, "w", encoding="utf-8") as f:
            f.writelines(new_topol_lines)

    except Exception as e:
        raise SimulationPrepError(f"Erro na Etapa F (Fusão de Topologia): {e}")
    yield "F", "success"

    # Etapa G: Definição da Caixa
    yield "G", "start"
    cmd_editconf = [
        gmx_bin,
        "editconf",
        "-f",
        "complex.gro",
        "-o",
        "complex_box.gro",
        "-c",
        "-d",
        "1.0",
        "-bt",
        "dodecahedron",
    ]
    run_command(cmd_editconf, output_dir, step_name="Etapa G (Definição da Caixa)")
    yield "G", "success"

    # Etapa H: Solvatação
    yield "H", "start"
    cmd_solvate = [
        gmx_bin,
        "solvate",
        "-cp",
        "complex_box.gro",
        "-cs",
        "spc216.gro",
        "-o",
        "complex_solv.gro",
        "-p",
        "topol.top",
    ]
    run_command(cmd_solvate, output_dir, step_name="Etapa H (Solvatação)")
    yield "H", "success"

    # Etapa I: Compilação de Íons
    yield "I", "start"
    try:
        project_root = Path(__file__).resolve().parent.parent.parent
        minim_mdp = project_root / "src" / "templates" / "mdp" / "minim.mdp"
        if not minim_mdp.exists():
            minim_mdp = Path("src/templates/mdp/minim.mdp").resolve()
            if not minim_mdp.exists():
                raise FileNotFoundError(
                    "Arquivo minim.mdp não encontrado no caminho esperado."
                )

        cmd_grompp_ions = [
            gmx_bin,
            "grompp",
            "-f",
            str(minim_mdp),
            "-c",
            "complex_solv.gro",
            "-p",
            "topol.top",
            "-o",
            "ions.tpr",
            "-maxwarn",
            "3",
        ]
        run_command(
            cmd_grompp_ions, output_dir, step_name="Etapa I (Compilação de Íons)"
        )
    except Exception as e:
        if isinstance(e, SimulationPrepError):
            raise e
        raise SimulationPrepError(f"Erro na Etapa I (Compilação de Íons): {e}")
    yield "I", "success"

    # Etapa J: Neutralização Automatizada
    yield "J", "start"
    cmd_genion = [
        gmx_bin,
        "genion",
        "-s",
        "ions.tpr",
        "-o",
        "complex_ions.gro",
        "-p",
        "topol.top",
        "-pname",
        "NA",
        "-nname",
        "CL",
        "-neutral",
        "-conc",
        "0.15",
    ]
    run_command(
        cmd_genion,
        output_dir,
        input_val=b"15\n",
        step_name="Etapa J (Neutralização Automatizada)",
    )
    yield "J", "success"

    # Etapa K: Grompp Definitivo
    yield "K", "start"
    try:
        project_root = Path(__file__).resolve().parent.parent.parent
        minim_mdp = project_root / "src" / "templates" / "mdp" / "minim.mdp"
        if not minim_mdp.exists():
            minim_mdp = Path("src/templates/mdp/minim.mdp").resolve()
            if not minim_mdp.exists():
                raise FileNotFoundError(
                    "Arquivo minim.mdp não encontrado no caminho esperado."
                )

        cmd_grompp_em = [
            gmx_bin,
            "grompp",
            "-f",
            str(minim_mdp),
            "-c",
            "complex_ions.gro",
            "-p",
            "topol.top",
            "-o",
            "em.tpr",
            "-maxwarn",
            "2",
        ]
        run_command(cmd_grompp_em, output_dir, step_name="Etapa K (Grompp Definitivo)")
    except Exception as e:
        if isinstance(e, SimulationPrepError):
            raise e
        raise SimulationPrepError(f"Erro na Etapa K (Grompp Definitivo): {e}")
    yield "K", "success"

    # Etapa L: Minimização de Energia
    yield "L", "start"
    cmd_mdrun = [gmx_bin, "mdrun", "-v", "-deffnm", "em"]
    run_command(cmd_mdrun, output_dir, step_name="Etapa L (Minimização de Energia)")
    yield "L", "success"
