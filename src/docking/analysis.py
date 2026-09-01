import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolAlign


def extract_vina_score(log_path: Path):
    try:
        with open(log_path, "r") as f:
            for line in f:
                if line.startswith("   1 "):
                    return float(line.split()[1])
    except Exception:
        return None


def analyze_results(docked_pdbqt: Path, reference_pdb: Path, results_dir: Path):
    sdf_out = results_dir / "docked_poses.sdf"

    # Exporta para SDF usando o meeko CLI
    from docking.preparation import get_executable

    exec_name = get_executable("mk_export")
    cmd = [exec_name, str(docked_pdbqt), "-s", str(sdf_out)]
    try:
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except Exception as e:
        return None, f"Falha ao gerar o arquivo SDF via {exec_name}: {str(e)}"

    if not sdf_out.exists():
        return None, "Falha ao gerar o arquivo SDF. O arquivo de saída não foi gerado."

    try:
        ref_mol = Chem.MolFromPDBFile(str(reference_pdb), removeHs=True)
        suppl = Chem.SDMolSupplier(str(sdf_out), removeHs=True)
        best_pose = suppl[0]

        if not best_pose:
            return None, "RDKit falhou em ler a pose do SDF."

        rmsd = rdMolAlign.GetBestRMS(best_pose, ref_mol)
        return rmsd, None
    except Exception as e:
        return None, str(e)


def generate_complex_pdb(receptor_pdb: Path, ligand_sdf: Path, output_pdb: Path):
    """
    Lê a primeira pose do ligante do arquivo SDF, altera seu resíduo para 'LIG',
    garante nomes de átomos ÚNICOS (ex: C1, C2, O1) para não quebrar a Dinâmica Molecular,
    e concatena como HETATM ao final do receptor.pdb.
    """
    suppl = Chem.SDMolSupplier(str(ligand_sdf))
    if not suppl:
        raise ValueError(
            f"Não foi possível ler o arquivo SDF do ligante em: {ligand_sdf}"
        )

    mol = next(iter(suppl))
    if mol is None:
        raise ValueError(f"Arquivo SDF inválido ou vazio: {ligand_sdf}")

    # Dicionário contador para gerar nomes únicos por elemento (C1, C2, O1...)
    element_counters = {}

    # Força propriedades PDB para o ligante com nomenclatura estrita
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        element_counters[symbol] = element_counters.get(symbol, 0) + 1

        # Cria um nome único de até 4 caracteres (ex: " C1  ", " O12 ")
        atom_name = f"{symbol}{element_counters[symbol]}"
        formatted_name = f" {atom_name:<3}"[:4]

        info = atom.GetPDBResidueInfo()
        if info is None:
            info = Chem.AtomPDBResidueInfo()

        info.SetName(formatted_name)
        info.SetResidueName("LIG")
        info.SetChainId("X")
        info.SetResidueNumber(1)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)

    # Converte o ligante para formato PDB em memória
    pdb_block = Chem.MolToPDBBlock(mol)

    # Filtra e isola as linhas de coordenadas do ligante
    ligand_lines = []
    for line in pdb_block.splitlines():
        if (line.startswith("HETATM") or line.startswith("ATOM")) and "LIG" in line:
            if line.startswith("ATOM  "):
                line = "HETATM" + line[6:]
            ligand_lines.append(line)

    # Lê o arquivo do receptor limpando travas de final de arquivo
    with open(receptor_pdb, "r") as f:
        receptor_content = f.read()

    receptor_lines = []
    for line in receptor_content.splitlines():
        if line.strip() not in ("END", "ENDMDL"):
            receptor_lines.append(line)

    # Consolida o arquivo do complexo
    with open(output_pdb, "w") as f:
        for line in receptor_lines:
            f.write(line + "\n")
        for line in ligand_lines:
            f.write(line + "\n")
        f.write("END\n")


def run_plip_docker(complex_pdb: Path, output_dir: Path):
    """
    Executa o container Docker 'pharmai/plip' via subprocess.
    Monta 'output_dir' como um volume no container e gera o relatório XML 'report.xml'.
    """
    complex_pdb = complex_pdb.resolve()
    output_dir = output_dir.resolve()

    # Remove relatórios pré-existentes para evitar leitura de dados desatualizados
    for old_xml in list(output_dir.glob("*report.xml")) + [output_dir / "report.xml"]:
        if old_xml.exists():
            try:
                old_xml.unlink()
            except Exception:
                pass

    # Garante que o complexo esteja dentro do volume que será montado
    if complex_pdb.parent != output_dir:
        shutil.copy(complex_pdb, output_dir / complex_pdb.name)
        complex_pdb = output_dir / complex_pdb.name

    # Prepara comando Docker
    cmd = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{output_dir}:/results",
        "-w",
        "/results",
        "pharmai/plip",
        "-f",
        complex_pdb.name,
        "-x",
    ]

    try:
        # Executa de forma síncrona capturando logs de erro se houver falhas
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Padroniza saída: se o PLIP gerou <stem>_report.xml, cria espelho report.xml
        xml_out = output_dir / "report.xml"
        if not xml_out.exists():
            stem_xml = output_dir / f"{complex_pdb.stem}_report.xml"
            if stem_xml.exists():
                shutil.copy2(stem_xml, xml_out)
            else:
                any_xmls = list(output_dir.glob("*_report.xml"))
                if any_xmls:
                    shutil.copy2(any_xmls[0], xml_out)

        return True, result.stdout
    except subprocess.CalledProcessError as e:
        error_msg = f"Falha na execução do PLIP via Docker. Código: {e.returncode}.\nStdout: {e.stdout}\nStderr: {e.stderr}"
        return False, error_msg
    except Exception as e:
        return False, f"Erro ao iniciar o container Docker do PLIP: {str(e)}"


def parse_plip_xml(xml_path: Path):
    """
    Realiza o parsing do relatório XML gerado pelo PLIP.
    Extrai as pontes de hidrogênio e contatos hidrofóbicos detectados.
    """
    interactions = {"hydrogen_bonds": [], "hydrophobic_contacts": []}
    xml_path = Path(xml_path)

    if not xml_path.exists():
        # Tenta buscar qualquer *_report.xml na mesma pasta
        candidates = list(xml_path.parent.glob("*_report.xml")) + list(
            xml_path.parent.glob("*report*.xml")
        )
        if candidates:
            xml_path = candidates[0]
        else:
            print(f"[DEBUG] Arquivo XML não localizado em: {xml_path}")
            return interactions

    print(f"[DEBUG] Arquivo XML localizado com sucesso: {xml_path}")

    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()

        bindingsites = root.findall(".//bindingsite")
        print(f"[DEBUG] {len(bindingsites)} bindingsite(s) detectado(s) no XML.")

        for bindingsite in bindingsites:
            # Extração de Pontes de Hidrogênio
            hbonds_node = bindingsite.find(".//hydrogen_bonds")
            if hbonds_node is not None:
                for hb in hbonds_node.findall("hydrogen_bond"):
                    resnr_el = hb.find("resnr")
                    restype_el = hb.find("restype")
                    dist_d_a_el = hb.find("dist_d-a")  # Distância Doador-Aceitador
                    dist_h_a_el = hb.find(
                        "dist_h-a"
                    )  # Backup: Distância Hidrogênio-Aceitador

                    if resnr_el is not None and restype_el is not None:
                        resname = (
                            restype_el.text.strip()
                            if restype_el.text and restype_el.text.strip()
                            else "UNK"
                        )

                        try:
                            resnr = (
                                int(resnr_el.text.strip())
                                if resnr_el.text and resnr_el.text.strip()
                                else 0
                            )
                        except ValueError:
                            resnr = 0

                        dist = 0.0
                        if (
                            dist_d_a_el is not None
                            and dist_d_a_el.text
                            and dist_d_a_el.text.strip()
                        ):
                            try:
                                dist = float(dist_d_a_el.text.strip())
                            except ValueError:
                                dist = 0.0
                        elif (
                            dist_h_a_el is not None
                            and dist_h_a_el.text
                            and dist_h_a_el.text.strip()
                        ):
                            try:
                                dist = float(dist_h_a_el.text.strip())
                            except ValueError:
                                dist = 0.0

                        interactions["hydrogen_bonds"].append(
                            {
                                "resname": resname,
                                "resnr": resnr,
                                "distance": dist,
                            }
                        )

            # Extração de Contatos Hidrofóbicos
            hydrophobic_node = bindingsite.find(".//hydrophobic_interactions")
            if hydrophobic_node is not None:
                for hc in hydrophobic_node.findall("hydrophobic_interaction"):
                    resnr_el = hc.find("resnr")
                    restype_el = hc.find("restype")
                    dist_el = hc.find("dist")

                    if resnr_el is not None and restype_el is not None:
                        resname = (
                            restype_el.text.strip()
                            if restype_el.text and restype_el.text.strip()
                            else "UNK"
                        )

                        try:
                            resnr = (
                                int(resnr_el.text.strip())
                                if resnr_el.text and resnr_el.text.strip()
                                else 0
                            )
                        except ValueError:
                            resnr = 0

                        dist = 0.0
                        if (
                            dist_el is not None
                            and dist_el.text
                            and dist_el.text.strip()
                        ):
                            try:
                                dist = float(dist_el.text.strip())
                            except ValueError:
                                dist = 0.0

                        interactions["hydrophobic_contacts"].append(
                            {
                                "resname": resname,
                                "resnr": resnr,
                                "distance": dist,
                            }
                        )

        total_hb = len(interactions["hydrogen_bonds"])
        total_hc = len(interactions["hydrophobic_contacts"])
        print(
            f"[DEBUG] Total de interações extraídas com sucesso: {total_hb + total_hc} "
            f"({total_hb} pontes de hidrogênio, {total_hc} contatos hidrofóbicos)"
        )

    except Exception as e:
        import traceback

        print(f"[DEBUG] Ocorreu uma exceção ao fazer o parsing do XML: {e}")
        traceback.print_exc()

    return interactions
