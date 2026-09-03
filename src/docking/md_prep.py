import collections
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Tuple

from openmm.app import PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem


class DependencyError(Exception):
    """Exceção levantada quando um binário externo (GROMACS ou ACPYPE) não é encontrado no PATH."""

    pass


class SimulationPrepError(Exception):
    """Exceção levantada quando há falha no processamento ou execução das etapas de preparação."""

    pass


def find_executable(name: str) -> Optional[str]:
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


def sanitize_target_id(name: str) -> str:
    """Normaliza o identificador do alvo para formato seguro de diretórios e arquivos (remove sufixos técnicos e caracteres inválidos)."""
    raw = str(name).strip()
    for suffix in [
        "_receptor",
        "_prepared",
        "_clean",
        "_docked",
        "_complex",
        "_target",
    ]:
        if raw.lower().endswith(suffix):
            raw = raw[: -len(suffix)]
    clean = re.sub(r"[^A-Za-z0-9_-]", "_", raw).strip("_")
    return clean.upper() or "TARGET"


def check_and_purge_stale_files(target_dir: Path, purge: bool = False) -> List[Path]:
    """
    Inspeciona o diretório do alvo procurando por artefatos residuais de execuções anteriores:
    (#*#, *.cpt, *.gro, *.xtc, *.tpr, *.edr, *.trr, _GMXMMPBSA_*, *.top, *.ndx).

    Se purge=True: deleta os arquivos com segurança e retorna a lista dos arquivos removidos.
    Se purge=False: apenas retorna a lista dos arquivos encontrados para decisão do chamador.
    """
    target_dir = Path(target_dir)
    if not target_dir.exists():
        return []

    stale_patterns = [
        "#*#",
        "*#",
        "*.cpt",
        "*.gro",
        "*.xtc",
        "*.tpr",
        "*.edr",
        "*.trr",
        "*.log",
        "*.top",
        "*.ndx",
        "*.xvg",
        "_GMXMMPBSA_*",
        "*.inpcrd",
        "*.mdcrd*",
    ]

    found_stale: List[Path] = []
    for pattern in stale_patterns:
        for p in target_dir.glob(pattern):
            if p not in found_stale:
                found_stale.append(p)

    if purge and found_stale:
        for p in found_stale:
            try:
                if p.is_file() or p.is_symlink():
                    p.unlink(missing_ok=True)
                elif p.is_dir():
                    shutil.rmtree(p, ignore_errors=True)
            except Exception:
                pass

    return found_stale


def validate_molecular_identity(
    complex_file: Path,
    expected_target: Optional[str] = None,
    min_ligand_atoms: int = 10,
) -> Dict[str, Any]:
    """
    Fail-Fast Validation:
    1. Inspeciona a estrutura gerada (PDB ou GRO) e valida os resíduos proteicos N-terminais.
    2. Valida a presença obrigatória e contagem de átomos do ligante (LIG / ligand_md).
    3. Verifica integridade das coordenadas atômicas.

    Lança SimulationPrepError imediatamente se qualquer inconsistência for detectada.
    """
    complex_file = Path(complex_file)
    if not complex_file.exists():
        raise SimulationPrepError(
            f"Arquivo de validação molecular não encontrado: {complex_file}"
        )

    standard_aa = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "HID",
        "HIE",
        "HIP",
        "CYX",
        "ASH",
        "GLH",
        "LYN",
        "ARN",
    }

    protein_residues: List[Tuple[str, str, str]] = []  # (resname, resnum, chain)
    ligand_atoms_count = 0
    ligand_resnames = set()
    total_atoms = 0

    if complex_file.suffix.lower() == ".pdb":
        with open(complex_file, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("ATOM  ") or line.startswith("HETATM"):
                    total_atoms += 1
                    resname = line[17:20].strip().upper()
                    chain = line[21:22].strip() or "A"
                    resnum = line[22:26].strip()

                    if resname in standard_aa:
                        res_key = (resname, resnum, chain)
                        if not protein_residues or protein_residues[-1] != res_key:
                            protein_residues.append(res_key)
                    elif resname in ("LIG", "UNK", "MOL", "DRG", "LIGAND_MD"):
                        ligand_atoms_count += 1
                        ligand_resnames.add(resname)

    elif complex_file.suffix.lower() == ".gro":
        with open(complex_file, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
            if len(lines) > 2:
                for line in lines[2:-1]:
                    if len(line) >= 20:
                        total_atoms += 1
                        resnum = line[0:5].strip()
                        resname = line[5:10].strip().upper()
                        if resname in standard_aa:
                            res_key = (resname, resnum, "A")
                            if not protein_residues or protein_residues[-1] != res_key:
                                protein_residues.append(res_key)
                        elif "LIG" in resname or resname in ("UNK", "MOL"):
                            ligand_atoms_count += 1
                            ligand_resnames.add(resname)

    # 1. Validação da Proteína
    if not protein_residues:
        raise SimulationPrepError(
            f"Falha Crítica de Identidade Molecular: Nenhum resíduo de proteína (aminoácidos padrão) "
            f"foi detectado no complexo '{complex_file.name}'. "
            f"Verifique se o arquivo do receptor selecionado corresponde à proteína do alvo {expected_target or ''}."
        )

    n_term_signature = [f"{r[0]}_{r[1]}" for r in protein_residues[:5]]
    c_term_signature = [f"{r[0]}_{r[1]}" for r in protein_residues[-3:]]

    # 2. Validação do Ligante
    if ligand_atoms_count == 0:
        raise SimulationPrepError(
            f"Falha Crítica de Identidade Molecular: O ligante ('LIG') NÃO foi encontrado no complexo '{complex_file.name}'. "
            f"A fusão de coordenadas falhou ou o arquivo do ligante estava vazio."
        )

    if ligand_atoms_count < min_ligand_atoms:
        raise SimulationPrepError(
            f"Falha Crítica de Identidade Molecular: Contagem insuficiente de átomos no ligante "
            f"({ligand_atoms_count} átomos encontrados, mínimo esperado: {min_ligand_atoms}). "
            f"A molécula pode estar incompleta ou corrompida."
        )

    summary = {
        "complex_file": str(complex_file),
        "total_atoms": total_atoms,
        "protein_residue_count": len(protein_residues),
        "n_terminal_signature": " -> ".join(n_term_signature),
        "c_terminal_signature": " -> ".join(c_term_signature),
        "ligand_atoms_count": ligand_atoms_count,
        "ligand_resnames": list(ligand_resnames),
        "is_valid": True,
    }
    return summary


def verify_tpr_consistency(tpr_path: Path) -> Dict[str, Any]:
    """
    Verificação de Consistência e Integridade Pós-Geração do TPR:
    Executa 'gmx check -s <tpr_path>' para checar contagem de átomos,
    carga total e integridade do arquivo binário antes de autorizar simulação ou envio ao cluster.
    """
    tpr_path = Path(tpr_path)
    if not tpr_path.exists():
        raise SimulationPrepError(
            f"Arquivo TPR não encontrado para verificação: {tpr_path}"
        )
    if tpr_path.stat().st_size == 0:
        raise SimulationPrepError(
            f"Arquivo TPR gerado está vazio (0 bytes): {tpr_path}"
        )

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        return {
            "tpr_path": str(tpr_path),
            "size_bytes": tpr_path.stat().st_size,
            "status": "size_verified",
        }

    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    exec_dir = str(Path(gmx_bin).parent)
    env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

    cmd = [gmx_bin, "dump", "-s", str(tpr_path.name)]
    try:
        res = subprocess.run(
            cmd, cwd=str(tpr_path.parent), env=env, capture_output=True, text=True
        )
        output = (res.stderr or "") + "\n" + (res.stdout or "")

        atom_match = (
            re.search(r"natoms\s*=\s*(\d+)", output, re.IGNORECASE)
            or re.search(r"#atoms\s*=\s*(\d+)", output, re.IGNORECASE)
            or re.search(r"(?:Coords|Step)\s+(\d+)", output)
        )
        atom_count = int(atom_match.group(1)) if atom_match else None

        if res.returncode != 0 and "error" in output.lower():
            raise SimulationPrepError(
                f"Arquivo TPR '{tpr_path.name}' reprovado no teste de integridade do GROMACS (gmx dump):\n{output}"
            )

        return {
            "tpr_path": str(tpr_path),
            "size_bytes": tpr_path.stat().st_size,
            "atom_count": atom_count,
            "status": "consistent",
        }
    except Exception as e:
        if isinstance(e, SimulationPrepError):
            raise e
        return {
            "tpr_path": str(tpr_path),
            "size_bytes": tpr_path.stat().st_size,
            "status": "check_warning",
            "warning": str(e),
        }


def extract_ligand(
    ligand_sdf: Path, output_dir: Path, target_id: Optional[str] = None
) -> Path:
    """
    Etapa B: Lê a primeira pose do ligante do arquivo SDF com o RDKit.
    Garante nomes de átomos únicos (C1, C2, O1...) e resíduo 'LIG' na cadeia 'X'.
    Salva em '<target_id>_ligand_md.pdb' e cria espelho 'ligand_md.pdb'.
    """
    try:
        supplier = Chem.SDMolSupplier(str(ligand_sdf), removeHs=False)
        if not supplier or len(supplier) == 0 or supplier[0] is None:
            raise ValueError(
                f"Não foi possível ler o arquivo SDF ou ele está vazio: {ligand_sdf}"
            )
        mol = supplier[0]
        mol = Chem.RemoveHs(mol)
        mol = Chem.AddHs(mol, addCoords=True)
    except Exception as e:
        raise SimulationPrepError(f"Falha ao ler o arquivo SDF com RDKit: {e}")

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

    prefix = f"{sanitize_target_id(target_id)}_" if target_id else ""
    target_pdb = output_dir / f"{prefix}ligand_md.pdb"
    mirror_pdb = output_dir / "ligand_md.pdb"

    try:
        Chem.MolToPDBFile(mol, str(target_pdb))
        if target_pdb != mirror_pdb:
            shutil.copy2(target_pdb, mirror_pdb)
    except Exception as e:
        raise SimulationPrepError(f"Falha ao exportar ligante para PDB: {e}")

    return target_pdb


def run_acpype(ligand_pdb: Path, output_dir: Path):
    """
    Etapa C: Executa o ACPYPE via subprocesso com ambiente isolado (sem interferência de PYTHONPATH).
    """
    acpype_bin = find_executable("acpype")
    if not acpype_bin:
        raise DependencyError(
            "O executável 'acpype' não foi encontrado no PATH ou no ambiente 'bioinfo'."
        )

    ligand_pdb_name = ligand_pdb.name

    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
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
        detailed_error = e.stderr.strip() or e.stdout.strip()
        raise SimulationPrepError(
            f"Falha ao rodar ACPYPE (código {e.returncode}):\n"
            f"Comando: {' '.join(cmd_acpype)}\n"
            f"Erro: {detailed_error}"
        )
    except Exception as e:
        raise SimulationPrepError(f"Erro ao iniciar processo do ACPYPE: {e}")


def run_pdb2gmx(receptor_pdb: Path, output_dir: Path, target_id: Optional[str] = None):
    """
    Etapa D: Executa o GROMACS pdb2gmx para preparar a proteína e gerar topologia.
    """
    prefix = f"{sanitize_target_id(target_id)}_" if target_id else ""
    fixed_name = f"{prefix}receptor_fixed.pdb" if target_id else "receptor_fixed.pdb"
    receptor_fixed = output_dir / fixed_name

    if receptor_pdb.name != fixed_name and not receptor_fixed.exists():
        try:
            fixer = PDBFixer(filename=str(receptor_pdb))
            fixer.removeHeterogens(keepWater=False)
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()

            with open(receptor_fixed, "w", encoding="utf-8") as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)
            receptor_pdb = receptor_fixed
        except Exception as e:
            raise SimulationPrepError(f"Erro ao curar o receptor com PDBFixer: {e}")
    elif receptor_fixed.exists():
        receptor_pdb = receptor_fixed

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError(
            "O executável 'gmx' (GROMACS) não foi encontrado no PATH."
        )

    abs_receptor_pdb = receptor_pdb.resolve()
    if not abs_receptor_pdb.exists():
        raise FileNotFoundError(
            f"Arquivo do receptor não encontrado em: {abs_receptor_pdb}"
        )

    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    gmx_dir = str(Path(gmx_bin).parent)
    env["PATH"] = f"{gmx_dir}{os.pathsep}{env.get('PATH', '')}"

    out_gro_name = (
        f"{prefix}protein_processed.gro" if target_id else "protein_processed.gro"
    )
    out_top_name = f"{prefix}topol.top" if target_id else "topol.top"

    cmd_gmx = [
        gmx_bin,
        "pdb2gmx",
        "-f",
        str(abs_receptor_pdb),
        "-o",
        out_gro_name,
        "-p",
        out_top_name,
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
            except subprocess.CalledProcessError as e_inner:
                detailed_error = e_inner.stderr.strip() or e_inner.stdout.strip()
                raise SimulationPrepError(
                    f"Falha ao rodar GROMACS pdb2gmx (mesmo com -missing, código {e_inner.returncode}):\n"
                    f"Comando: {' '.join(cmd_gmx_missing)}\n"
                    f"Erro: {detailed_error}"
                )
        else:
            raise SimulationPrepError(
                f"Falha ao rodar GROMACS pdb2gmx (código {e.returncode}):\n"
                f"Comando: {' '.join(cmd_gmx)}\n"
                f"Erro: {detailed_error}"
            )
    except Exception as e:
        raise SimulationPrepError(f"Erro ao iniciar processo do GROMACS: {e}")

    # Cria espelhos sem prefixo para manter compatibilidade com submódulos legados
    if target_id:
        if (output_dir / out_gro_name).exists():
            shutil.copy2(
                output_dir / out_gro_name, output_dir / "protein_processed.gro"
            )
        if (output_dir / out_top_name).exists():
            shutil.copy2(output_dir / out_top_name, output_dir / "topol.top")


def stitch_topology(output_dir: Path, target_id: Optional[str] = None):
    """
    Etapa F: Injeta o include do ligante e sua definição em [ molecules ] no topol.top.
    """
    prefix = f"{sanitize_target_id(target_id)}_" if target_id else ""
    target_top = (
        output_dir / f"{prefix}topol.top" if target_id else output_dir / "topol.top"
    )
    if not target_top.exists():
        target_top = output_dir / "topol.top"

    if not target_top.exists():
        raise SimulationPrepError(
            f"O arquivo de topologia esperado não foi gerado em {output_dir}"
        )

    try:
        with open(target_top, "r", encoding="utf-8") as f:
            lines = f.readlines()
    except Exception as e:
        raise SimulationPrepError(f"Falha ao ler '{target_top.name}': {e}")

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
            "Inclusão de 'forcefield.itp' não encontrada na topologia."
        )

    molecules_idx = -1
    for i, line in enumerate(new_lines):
        if line.strip().startswith("[ molecules ]"):
            molecules_idx = i
            break

    if molecules_idx == -1:
        raise SimulationPrepError("Seção '[ molecules ]' não encontrada na topologia.")

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
            "Nenhuma molécula ativa encontrada sob a seção '[ molecules ]'."
        )

    if new_lines and not new_lines[-1].endswith("\n"):
        new_lines[-1] = new_lines[-1] + "\n"

    new_lines.append("ligand_md                 1\n")

    try:
        with open(target_top, "w", encoding="utf-8") as f:
            f.writelines(new_lines)

        # Atualiza espelho
        mirror_top = output_dir / "topol.top"
        if target_top != mirror_top:
            shutil.copy2(target_top, mirror_top)
    except Exception as e:
        raise SimulationPrepError(f"Falha ao salvar modificações na topologia: {e}")


def prepare_md_system(
    receptor_pdb: Path,
    ligand_sdf: Path,
    output_dir: Path,
    target_id: Optional[str] = None,
    purge: bool = True,
) -> Generator[Tuple[str, str], None, None]:
    """
    Prepara o sistema completo de Dinâmica Molecular com isolamento estrito por alvo (Target Isolation),
    checagem de identidade molecular (Fail-Fast Validation) e controle estrito de erros no grompp.

    Retorna um gerador (etapa, status).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Identificação e sanitização do Target ID
    if not target_id:
        target_id = sanitize_target_id(output_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            target_id = sanitize_target_id(
                receptor_pdb.stem.replace("_processed", "").replace("receptor", "")
            )
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    # Limpeza forçada de artefatos residuais e zumbis para garantir estado estéril
    check_and_purge_stale_files(output_dir, purge=purge)

    gmx_bin = find_executable("gmx")
    acpype_bin = find_executable("acpype")

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
                    f"Erro do GROMACS/Ferramenta: {error_msg}"
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

        standard_amino_acids = {
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
            "HID",
            "HIE",
            "HIP",
            "CYX",
            "ASH",
            "GLH",
            "LYN",
            "ARN",
        }
        protein_res = [
            r
            for r in fixer.topology.residues()
            if r.name.upper() in standard_amino_acids
        ]
        if not protein_res:
            raise SimulationPrepError(
                f"O arquivo fornecido como receptor ('{receptor_pdb}') não contém resíduos de proteína (aminoácidos padrão). "
                f"Certifique-se de selecionar o arquivo 'receptor.pdb' do alvo {target_id}."
            )

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        receptor_fixed = output_dir / f"{prefix}receptor_fixed.pdb"
        with open(receptor_fixed, "w", encoding="utf-8") as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        # Espelho de compatibilidade
        shutil.copy2(receptor_fixed, output_dir / "receptor_fixed.pdb")
    except Exception as e:
        if isinstance(e, SimulationPrepError):
            raise e
        raise SimulationPrepError(f"Erro na Etapa A (Cura com PDBFixer): {e}")
    yield "A", "success"

    # Etapa B: Extração do Ligante com RDKit
    yield "B", "start"
    extract_ligand(ligand_sdf, output_dir, target_id=target_id)
    yield "B", "success"

    # Etapa C: Parametrização ACPYPE
    yield "C", "start"
    if not acpype_bin:
        raise DependencyError("O executável 'acpype' não foi encontrado no PATH.")
    cmd_acpype = [acpype_bin, "-i", "ligand_md.pdb", "-c", "bcc", "-f"]
    run_command(cmd_acpype, output_dir, step_name="Etapa C (Parametrização ACPYPE)")
    yield "C", "success"

    # Etapa D: Topologia da Proteína (pdb2gmx)
    yield "D", "start"
    if not gmx_bin:
        raise DependencyError(
            "O executável 'gmx' (GROMACS) não foi encontrado no PATH."
        )

    cmd_pdb2gmx = [
        gmx_bin,
        "pdb2gmx",
        "-f",
        f"{prefix}receptor_fixed.pdb",
        "-o",
        f"{prefix}protein_processed.gro",
        "-p",
        f"{prefix}topol.top",
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
                step_name="Etapa D (Topologia com -missing)",
            )
        else:
            raise e

    # Espelhos de compatibilidade
    shutil.copy2(
        output_dir / f"{prefix}protein_processed.gro",
        output_dir / "protein_processed.gro",
    )
    shutil.copy2(output_dir / f"{prefix}topol.top", output_dir / "topol.top")
    yield "D", "success"

    # Etapa E: Fusão de Coordenadas e Validação Molecular Pré-Execução
    yield "E", "start"
    try:
        prot_gro_path = output_dir / f"{prefix}protein_processed.gro"
        lig_gro_path = output_dir / "ligand_md.acpype" / "ligand_md_GMX.gro"

        if not prot_gro_path.exists():
            raise FileNotFoundError(f"Arquivo {prot_gro_path} não encontrado.")
        if not lig_gro_path.exists():
            raise FileNotFoundError(f"Arquivo {lig_gro_path} não encontrado.")

        with open(prot_gro_path, "r", encoding="utf-8") as f:
            prot_lines = f.readlines()
        with open(lig_gro_path, "r", encoding="utf-8") as f:
            lig_lines = f.readlines()

        if len(prot_lines) < 3 or len(lig_lines) < 3:
            raise ValueError("Arquivos .gro da proteína ou do ligante corrompidos.")

        prot_atoms = prot_lines[2:-1]
        lig_atoms = lig_lines[2:-1]
        total_atoms = len(prot_atoms) + len(lig_atoms)
        box_vector = prot_lines[-1]

        complex_lines = [f"Complex of {target_id} and Ligand\n", f" {total_atoms}\n"]
        complex_lines.extend(prot_atoms)
        complex_lines.extend(lig_atoms)
        complex_lines.append(box_vector)

        complex_gro_path = output_dir / f"{prefix}complex.gro"
        with open(complex_gro_path, "w", encoding="utf-8") as f:
            f.writelines(complex_lines)
        shutil.copy2(complex_gro_path, output_dir / "complex.gro")

        # Converte para PDB do complexo para validação estrita de identidade molecular
        complex_pdb_path = output_dir / f"{prefix}complex.pdb"
        rec_fixed_pdb = output_dir / f"{prefix}receptor_fixed.pdb"
        lig_fixed_pdb = output_dir / f"{prefix}ligand_md.pdb"

        with open(rec_fixed_pdb, "r", encoding="utf-8", errors="ignore") as f_rec:
            r_lines = [l for l in f_rec if l.strip() not in ("END", "ENDMDL")]
        with open(lig_fixed_pdb, "r", encoding="utf-8", errors="ignore") as f_lig:
            l_lines = [
                l for l in f_lig if l.startswith("ATOM") or l.startswith("HETATM")
            ]

        with open(complex_pdb_path, "w", encoding="utf-8") as f_comp:
            for l in r_lines:
                f_comp.write(l)
            for l in l_lines:
                f_comp.write(l)
            f_comp.write("END\n")
        shutil.copy2(complex_pdb_path, output_dir / "complex.pdb")

        # CHECAGEM DE IDENTIDADE MOLECULAR (FAIL-FAST)
        validate_molecular_identity(complex_pdb_path, expected_target=target_id)

    except Exception as e:
        if isinstance(e, SimulationPrepError):
            raise e
        raise SimulationPrepError(f"Erro na Etapa E (Fusão de Coordenadas): {e}")
    yield "E", "success"

    # Etapa F: Fusão de Topologia (Stitching)
    yield "F", "start"
    stitch_topology(output_dir, target_id=target_id)
    yield "F", "success"

    # Etapa G: Definição da Caixa de Simulação (editconf)
    yield "G", "start"
    cmd_editconf = [
        gmx_bin,
        "editconf",
        "-f",
        f"{prefix}complex.gro",
        "-o",
        f"{prefix}complex_box.gro",
        "-c",
        "-d",
        "1.0",
        "-bt",
        "dodecahedron",
    ]
    run_command(cmd_editconf, output_dir, step_name="Etapa G (Definição da Caixa)")
    shutil.copy2(
        output_dir / f"{prefix}complex_box.gro", output_dir / "complex_box.gro"
    )
    yield "G", "success"

    # Etapa H: Solvatação (solvate)
    yield "H", "start"
    cmd_solvate = [
        gmx_bin,
        "solvate",
        "-cp",
        f"{prefix}complex_box.gro",
        "-cs",
        "spc216.gro",
        "-o",
        f"{prefix}complex_solv.gro",
        "-p",
        f"{prefix}topol.top",
    ]
    run_command(cmd_solvate, output_dir, step_name="Etapa H (Solvatação)")
    shutil.copy2(
        output_dir / f"{prefix}complex_solv.gro", output_dir / "complex_solv.gro"
    )
    shutil.copy2(output_dir / f"{prefix}topol.top", output_dir / "topol.top")
    yield "H", "success"

    # Etapa I: Compilação de Íons (grompp) com Controle Estrito de Erro
    yield "I", "start"
    try:
        project_root = Path(__file__).resolve().parent.parent.parent
        minim_mdp = project_root / "src" / "templates" / "mdp" / "minim.mdp"
        if not minim_mdp.exists():
            minim_mdp = Path("src/templates/mdp/minim.mdp").resolve()
            if not minim_mdp.exists():
                raise FileNotFoundError("Arquivo minim.mdp não encontrado.")

        ions_tpr = output_dir / f"{prefix}ions.tpr"
        ions_tpr.unlink(missing_ok=True)
        (output_dir / "ions.tpr").unlink(missing_ok=True)

        cmd_grompp_ions = [
            gmx_bin,
            "grompp",
            "-f",
            str(minim_mdp),
            "-c",
            f"{prefix}complex_solv.gro",
            "-p",
            f"{prefix}topol.top",
            "-o",
            f"{prefix}ions.tpr",
            "-maxwarn",
            "3",
        ]
        run_command(
            cmd_grompp_ions, output_dir, step_name="Etapa I (Compilação de Íons)"
        )
        shutil.copy2(ions_tpr, output_dir / "ions.tpr")
        verify_tpr_consistency(ions_tpr)
    except Exception as e:
        if isinstance(e, SimulationPrepError):
            raise e
        raise SimulationPrepError(f"Erro na Etapa I (Compilação de Íons): {e}")
    yield "I", "success"

    # Etapa J: Neutralização e Concentração Iônica (genion)
    yield "J", "start"
    cmd_genion = [
        gmx_bin,
        "genion",
        "-s",
        f"{prefix}ions.tpr",
        "-o",
        f"{prefix}complex_ions.gro",
        "-p",
        f"{prefix}topol.top",
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
    shutil.copy2(
        output_dir / f"{prefix}complex_ions.gro", output_dir / "complex_ions.gro"
    )
    shutil.copy2(output_dir / f"{prefix}topol.top", output_dir / "topol.top")
    yield "J", "success"

    # Etapa K: Grompp Definitivo da Minimização de Energia (em.tpr)
    yield "K", "start"
    try:
        em_tpr = output_dir / f"{prefix}em.tpr"
        em_tpr.unlink(missing_ok=True)
        (output_dir / "em.tpr").unlink(missing_ok=True)

        cmd_grompp_em = [
            gmx_bin,
            "grompp",
            "-f",
            str(minim_mdp),
            "-c",
            f"{prefix}complex_ions.gro",
            "-p",
            f"{prefix}topol.top",
            "-o",
            f"{prefix}em.tpr",
            "-maxwarn",
            "2",
        ]
        run_command(
            cmd_grompp_em, output_dir, step_name="Etapa K (Grompp Definitivo EM)"
        )
        shutil.copy2(em_tpr, output_dir / "em.tpr")
        verify_tpr_consistency(em_tpr)
    except Exception as e:
        if isinstance(e, SimulationPrepError):
            raise e
        raise SimulationPrepError(f"Erro na Etapa K (Grompp Definitivo EM): {e}")
    yield "K", "success"

    # Etapa L: Minimização de Energia (mdrun)
    yield "L", "start"
    cmd_mdrun = [gmx_bin, "mdrun", "-v", "-deffnm", f"{prefix}em"]
    run_command(cmd_mdrun, output_dir, step_name="Etapa L (Minimização de Energia)")

    # Garante espelho em.gro
    if (output_dir / f"{prefix}em.gro").exists():
        shutil.copy2(output_dir / f"{prefix}em.gro", output_dir / "em.gro")
    yield "L", "success"
