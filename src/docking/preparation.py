import urllib.request
from pathlib import Path
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
import shutil


def get_executable(name: str) -> str:
    """Retorna o nome do executável disponível no PATH (com ou sem extensão .py)."""
    if shutil.which(name):
        return name
    if shutil.which(f"{name}.py"):
        return f"{name}.py"
    return name


def download_pdb(pdb_id: str, raw_dir: Path) -> Path:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    pdb_file = raw_dir / f"{pdb_id}.pdb"
    if not pdb_file.exists():
        urllib.request.urlretrieve(url, pdb_file)
    return pdb_file


def download_pubchem_sdf(cid: str, out_path: Path) -> Path:
    """Baixa um composto do PubChem em formato SDF 3D."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/record/SDF/?record_type=3d"
    if not out_path.exists():
        urllib.request.urlretrieve(url, out_path)
    return out_path


def split_receptor_ligand(pdb_file: Path, processed_dir: Path):
    receptor_path = processed_dir / "receptor.pdb"
    ligand_path = processed_dir / "ligand.pdb"

    with (
        open(pdb_file, "r") as fin,
        open(receptor_path, "w") as frec,
        open(ligand_path, "w") as flig,
    ):
        for line in fin:
            if line.startswith("ATOM"):
                frec.write(line)
            elif line.startswith("HETATM") and "HOH" not in line and "WAT" not in line:
                flig.write(line)
            elif line.startswith("TER") or line.startswith("END"):
                frec.write(line)
    return receptor_path, ligand_path


def prepare_receptor(pdb_path: Path, out_path: Path):
    """
    Prepara o receptor (PDBQT) a partir de um PDB.
    Utiliza o mk_prepare_receptor do Meeko diretamente.
    """
    # Usamos o arquivo PDB original gerado pelo split
    # Adicionamos flags para lidar com altlocs e resíduos problemáticos
    import subprocess

    exec_name = get_executable("mk_prepare_receptor")
    cmd = [
        exec_name,
        "--read_pdb",
        str(pdb_path),
        "-o",
        str(out_path.with_suffix("")),
        "--default_altloc",
        "A",
        "--allow_bad_res",
        "-p",
    ]

    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Falha na conversão do receptor via {exec_name}.\n"
            f"Comando: {' '.join(cmd)}\n"
            f"Erro (stderr): {e.stderr}\n"
            f"Saída (stdout): {e.stdout}"
        )
    except FileNotFoundError:
        raise RuntimeError(
            f"O comando '{exec_name}' não foi encontrado no PATH.\n"
            f"Certifique-se de que o ambiente virtual está ativo e que o pacote 'meeko' está instalado."
        )

    if not out_path.exists():
        raise RuntimeError(
            f"Falha na conversão do receptor. O arquivo de saída {out_path} não foi gerado."
        )


def prepare_ligand_sdf(sdf_path: Path, out_path: Path):
    """Prepara o ligante (PDBQT) a partir de um SDF usando mk_prepare_ligand."""
    import subprocess

    exec_name = get_executable("mk_prepare_ligand")
    cmd = [exec_name, "-i", str(sdf_path), "-o", str(out_path)]
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Falha na conversão do ligante via {exec_name}.\n"
            f"Comando: {' '.join(cmd)}\n"
            f"Erro (stderr): {e.stderr}\n"
            f"Saída (stdout): {e.stdout}"
        )
    except FileNotFoundError:
        raise RuntimeError(
            f"O comando '{exec_name}' não foi encontrado no PATH.\n"
            f"Certifique-se de que o ambiente virtual está ativo e que o pacote 'meeko' está instalado."
        )


def prepare_ligand(pdb_path: Path, out_path: Path):
    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
    if not mol:
        raise ValueError("RDKit não conseguiu processar o PDB do ligante.")

    frags = Chem.GetMolFrags(mol, asMols=True)
    # Isola o maior fragmento (assume ser o ligante de interesse)
    ligand = max(frags, key=lambda m: m.GetNumAtoms())
    ligand = Chem.AddHs(ligand, addCoords=True)

    # 1. Tentativa com charge_model padrão (Gasteiger)
    prep = MoleculePreparation(charge_model="gasteiger")
    setups = prep.prepare(ligand)
    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setups[0])

    # 2. Fallback caso Gasteiger gere NaN/Inf (comum para fosfatos/grupos inorgânicos sem ordens de ligação explícitas no PDB)
    if not is_ok:
        prep_fallback = MoleculePreparation(charge_model="zero")
        setups_fallback = prep_fallback.prepare(ligand)
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(
            setups_fallback[0]
        )

    if not is_ok or not pdbqt_string.strip():
        raise RuntimeError(f"Falha ao gerar PDBQT do ligante: {error_msg}")

    with open(out_path, "w") as f:
        f.write(pdbqt_string)
