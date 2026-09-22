import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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


def _find_pdb_for_xml(
    xml_path: Path, pdb_path: Optional[Path] = None
) -> Optional[Path]:
    """
    Localiza o arquivo PDB correspondente (complex.pdb ou receptor.pdb)
    para validação e resgate de rótulos de resíduos com códigos de inserção.
    """
    if pdb_path is not None:
        p = Path(pdb_path).resolve()
        if p.exists():
            return p

    xml_path = Path(xml_path).resolve()
    xml_dir = xml_path.parent

    # 1. complex.pdb direto no diretório do XML
    direct_complex = xml_dir / "complex.pdb"
    if direct_complex.exists():
        return direct_complex

    # 2. PDB com stem similar (ex.: 1H1B_complex_report.xml -> 1H1B_complex.pdb)
    stem_clean = (
        xml_path.stem.replace("_report", "").replace("report", "").strip("_")
    )
    if stem_clean:
        stem_pdb = xml_dir / f"{stem_clean}.pdb"
        if stem_pdb.exists():
            return stem_pdb

    # 3. Qualquer *complex*.pdb (evitando plipfixed)
    complex_cands = [
        p
        for p in xml_dir.glob("*complex*.pdb")
        if not p.name.startswith("plipfixed")
    ]
    if complex_cands:
        return complex_cands[0]

    # 4. receptor.pdb no diretório
    rec_pdb = xml_dir / "receptor.pdb"
    if rec_pdb.exists():
        return rec_pdb

    # 5. Qualquer *.pdb válido no diretório
    any_pdbs = [
        p
        for p in xml_dir.glob("*.pdb")
        if not p.name.startswith("plipfixed") and "clean" not in p.name
    ]
    if any_pdbs:
        return any_pdbs[0]

    # 6. Diretório pai ou irmãos estruturais (processed/receptor.pdb, results/complex.pdb)
    parent_cands = [
        xml_dir.parent / "processed" / "receptor.pdb",
        xml_dir.parent / "results" / "complex.pdb",
        xml_dir.parent / "complex.pdb",
    ]
    for pc in parent_cands:
        if pc.exists():
            return pc

    return None


def _load_pdb_residue_mapping(pdb_path: Optional[Path]):
    """
    Lê o arquivo PDB e constrói tabelas de dispersão por coordenadas 3D, número serial de átomo
    e índice sequencial para resgatar o identificador canônico do resíduo (resname, resnr com
    código de inserção como '62B' e chain).
    """
    coords_map: Dict[Tuple[float, float, float], Dict[str, Any]] = {}
    serials_map: Dict[int, Dict[str, Any]] = {}
    index_map: Dict[int, Dict[str, Any]] = {}
    atom_list: List[Tuple[Tuple[float, float, float], Dict[str, Any]]] = []

    if not pdb_path or not Path(pdb_path).exists():
        return coords_map, serials_map, index_map, atom_list

    atom_count = 0
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue

                resname = line[17:20].strip()
                # Ignora ligantes e solventes adicionados ao complexo
                if resname in ("LIG", "HOH", "WAT", "DOD"):
                    continue

                chain = line[21].strip() if len(line) > 21 else ""
                resnum_str = line[22:26].strip() if len(line) > 25 else ""
                icode = line[26].strip() if len(line) > 26 else ""

                if not resnum_str:
                    continue

                try:
                    resnum_int = int(resnum_str)
                    # Preserva a letra de inserção se existir (ex.: '62B'), senão int puro (ex.: 62)
                    resnr_val: Any = f"{resnum_int}{icode}" if icode else resnum_int
                except ValueError:
                    resnr_val = f"{resnum_str}{icode}".strip()

                serial = None
                try:
                    serial = int(line[6:11].strip())
                except ValueError:
                    pass

                coord = None
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coord = (round(x, 3), round(y, 3), round(z, 3))
                except ValueError:
                    pass

                atom_count += 1
                atom_info = {
                    "resname": resname,
                    "resnr": resnr_val,
                    "chain": chain,
                    "serial": serial,
                    "coord": coord,
                }

                if coord is not None:
                    coords_map[coord] = atom_info
                    atom_list.append((coord, atom_info))
                if serial is not None:
                    serials_map[serial] = atom_info
                index_map[atom_count] = atom_info
    except Exception as e:
        print(f"[DEBUG] Erro ao construir mapeamento de resíduos do PDB ({pdb_path}): {e}")

    return coords_map, serials_map, index_map, atom_list


def _resolve_protein_residue(
    inter_node: ET.Element,
    coords_map: Dict[Tuple[float, float, float], Dict[str, Any]],
    serials_map: Dict[int, Dict[str, Any]],
    index_map: Dict[int, Dict[str, Any]],
    atom_list: List[Tuple[Tuple[float, float, float], Dict[str, Any]]],
    fallback_resname: str,
    fallback_resnr_raw: str,
    fallback_chain: str = "",
) -> Tuple[str, Any, str]:
    """
    Resolve o rótulo do resíduo (resname, resnr, chain), resgatando o código de inserção real
    (ex.: 'VAL62B' ou resnr '62B') quando o PLIP gera <resnr>0</resnr> ou descarta a letra.
    """
    matched = None

    # 1. Busca por coordenadas atômicas exatas em <protcoo>
    protcoo = inter_node.find("protcoo")
    if protcoo is not None:
        try:
            x_el = protcoo.find("x")
            y_el = protcoo.find("y")
            z_el = protcoo.find("z")
            if (
                x_el is not None
                and y_el is not None
                and z_el is not None
                and x_el.text
                and y_el.text
                and z_el.text
            ):
                px = float(x_el.text.strip())
                py = float(y_el.text.strip())
                pz = float(z_el.text.strip())
                key = (round(px, 3), round(py, 3), round(pz, 3))
                if key in coords_map:
                    matched = coords_map[key]
                else:
                    # Tolerância euclidiana para diferenças mínimas de arredondamento (< 0.05 Å)
                    min_d2 = 0.05 * 0.05
                    best_info = None
                    for (cx, cy, cz), info in atom_list:
                        d2 = (cx - px) ** 2 + (cy - py) ** 2 + (cz - pz) ** 2
                        if d2 < min_d2:
                            min_d2 = d2
                            best_info = info
                    if best_info is not None:
                        matched = best_info
        except Exception:
            pass

    # 2. Busca por índice/serial de átomo da proteína
    if not matched:
        protisdon = inter_node.find("protisdon")
        candidate_tags = []
        if protisdon is not None and protisdon.text and protisdon.text.strip().lower() == "true":
            candidate_tags = ["donoridx"]
        elif protisdon is not None and protisdon.text and protisdon.text.strip().lower() == "false":
            candidate_tags = ["acceptoridx"]
        candidate_tags.extend(["protcarbonidx", "protidx", "donoridx", "acceptoridx"])

        for tag in candidate_tags:
            el = inter_node.find(tag)
            if el is not None and el.text and el.text.strip():
                try:
                    idx_val = int(el.text.strip())
                    if idx_val in serials_map:
                        matched = serials_map[idx_val]
                        break
                    elif idx_val in index_map:
                        matched = index_map[idx_val]
                        break
                except ValueError:
                    pass

    # 3. Se casou com um átomo do PDB de entrada, resgata a anotação canônica verdadeira
    if matched:
        return matched["resname"], matched["resnr"], matched["chain"]

    # 4. Fallback com o dado bruto do XML
    cleaned_resnr: Any = 0
    if fallback_resnr_raw and fallback_resnr_raw.strip():
        raw_val = fallback_resnr_raw.strip()
        m = re.match(r"^(\d+)([A-Za-z]?)$", raw_val)
        if m:
            num_part, icode_part = m.groups()
            cleaned_resnr = f"{int(num_part)}{icode_part}" if icode_part else int(num_part)
        else:
            try:
                cleaned_resnr = int(raw_val)
            except ValueError:
                cleaned_resnr = raw_val

    return fallback_resname, cleaned_resnr, fallback_chain


def parse_plip_xml(xml_path: Path, pdb_path: Optional[Path] = None):
    """
    Realiza o parsing do relatório XML gerado pelo PLIP.
    Extrai as pontes de hidrogênio e contatos hidrofóbicos detectados,
    cruzando coordenadas e índices com o PDB original para preservar
    códigos de inserção canônicos (ex.: 'VAL62B', 'ASN62A') em vez de resíduos '0'.
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

    # Localiza e mapeia o PDB de referência para recuperar anotações verdadeiras
    ref_pdb = _find_pdb_for_xml(xml_path, pdb_path)
    if ref_pdb:
        print(f"[DEBUG] PDB de referência para resgate de resíduos: {ref_pdb}")
    coords_map, serials_map, index_map, atom_list = _load_pdb_residue_mapping(ref_pdb)

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
                    reschain_el = hb.find("reschain")
                    dist_d_a_el = hb.find("dist_d-a")  # Distância Doador-Aceitador
                    dist_h_a_el = hb.find(
                        "dist_h-a"
                    )  # Backup: Distância Hidrogênio-Aceitador

                    raw_resname = (
                        restype_el.text.strip()
                        if restype_el is not None and restype_el.text and restype_el.text.strip()
                        else "UNK"
                    )
                    raw_resnr = (
                        resnr_el.text.strip()
                        if resnr_el is not None and resnr_el.text and resnr_el.text.strip()
                        else "0"
                    )
                    raw_chain = (
                        reschain_el.text.strip()
                        if reschain_el is not None and reschain_el.text and reschain_el.text.strip()
                        else ""
                    )

                    resname, resnr, chain = _resolve_protein_residue(
                        hb,
                        coords_map,
                        serials_map,
                        index_map,
                        atom_list,
                        raw_resname,
                        raw_resnr,
                        raw_chain,
                    )

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
                            "reschain": chain,
                            "distance": dist,
                        }
                    )

            # Extração de Contatos Hidrofóbicos
            hydrophobic_node = bindingsite.find(".//hydrophobic_interactions")
            if hydrophobic_node is not None:
                for hc in hydrophobic_node.findall("hydrophobic_interaction"):
                    resnr_el = hc.find("resnr")
                    restype_el = hc.find("restype")
                    reschain_el = hc.find("reschain")
                    dist_el = hc.find("dist")

                    raw_resname = (
                        restype_el.text.strip()
                        if restype_el is not None and restype_el.text and restype_el.text.strip()
                        else "UNK"
                    )
                    raw_resnr = (
                        resnr_el.text.strip()
                        if resnr_el is not None and resnr_el.text and resnr_el.text.strip()
                        else "0"
                    )
                    raw_chain = (
                        reschain_el.text.strip()
                        if reschain_el is not None and reschain_el.text and reschain_el.text.strip()
                        else ""
                    )

                    resname, resnr, chain = _resolve_protein_residue(
                        hc,
                        coords_map,
                        serials_map,
                        index_map,
                        atom_list,
                        raw_resname,
                        raw_resnr,
                        raw_chain,
                    )

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
                            "reschain": chain,
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
