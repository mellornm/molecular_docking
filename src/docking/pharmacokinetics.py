from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen


def calculate_admet_descriptors(ligand_sdf: Path) -> dict:
    """
    Calcula descritores físico-químicos e predições farmacocinéticas/toxicológicas avançadas (ADMET)
    usando regras moleculares e QSAR nativos baseados em RDKit.
    Analisa a primeira pose válida contida no arquivo SDF fornecido.
    """
    ligand_sdf = Path(ligand_sdf)
    if not ligand_sdf.exists():
        raise FileNotFoundError(f"Arquivo SDF do ligante não encontrado em: {ligand_sdf}")

    try:
        suppl = Chem.SDMolSupplier(str(ligand_sdf))
        mol = None
        for m in suppl:
            if m is not None:
                mol = m
                break
    except Exception as e:
        raise RuntimeError(f"Erro ao ler o arquivo SDF com o RDKit: {e}")

    if mol is None:
        raise ValueError(f"RDKit não conseguiu identificar uma molécula válida no arquivo: {ligand_sdf}")

    try:
        # Cálculo dos descritores físico-químicos com RDKit
        mw = float(Descriptors.ExactMolWt(mol))
        logp = float(Crippen.MolLogP(mol))
        hbd = int(Descriptors.NumHDonors(mol))
        hba = int(Descriptors.NumHAcceptors(mol))
        tpsa = float(Descriptors.TPSA(mol))
        rotb = int(Descriptors.NumRotatableBonds(mol))
        charge = int(Chem.GetFormalCharge(mol))
    except Exception as e:
        raise RuntimeError(f"Erro ao calcular os descritores moleculares com RDKit: {e}")

    # Validação da Regra de Cinco de Lipinski
    lipinski_violations = []
    if mw > 500.0:
        lipinski_violations.append(f"Peso Molecular elevado ({mw:.2f} > 500)")
    if logp > 5.0:
        lipinski_violations.append(f"LogP elevado ({logp:.2f} > 5)")
    if hbd > 5:
        lipinski_violations.append(f"Doadores de H em excesso ({hbd} > 5)")
    if hba > 10:
        lipinski_violations.append(f"Aceitadores de H em excesso ({hba} > 10)")

    lipinski_pass = len(lipinski_violations) <= 1

    # Validação das Regras de Veber
    veber_violations = []
    if rotb > 10:
        veber_violations.append(f"Ligações rotacionáveis em excesso ({rotb} > 10)")
    if tpsa > 140.0:
        veber_violations.append(f"TPSA elevado ({tpsa:.2f} > 140)")

    veber_pass = len(veber_violations) == 0

    # 1. Absorção Intestinal Humana (HIA) - Filtro de Egan (Egan Egg)
    # Alta absorção se TPSA <= 132 e -1.0 <= LogP <= 5.8
    hia_ok = (tpsa <= 132.0) and (-1.0 <= logp <= 5.8)
    hia_status = "Alta Absorção" if hia_ok else "Baixa Absorção"

    # 2. Permeabilidade da Barreira Hematoencefálica (BBB) - Regra de Clark
    # Se neutra, TPSA < 90 e LogP entre 1 e 5 -> Permeável
    bbb_ok = (charge == 0) and (tpsa < 90.0) and (1.0 <= logp <= 5.0)
    bbb_status = "Permeável" if bbb_ok else "Incompatível/Baixa"

    # 3. Substrato de P-glicoproteína (P-gp) - Modelo baseado em carga/tamanho
    # Moléculas com MW > 400 e TPSA > 80 tendem a ser substratos (efluxo provável)
    pgp_substrate = (mw > 400.0) and (tpsa > 80.0)
    pgp_status = "Substrato (Efluxo provável)" if pgp_substrate else "Não Substrato (Baixo Efluxo)"

    # 4. Alerta de Toxicidade (PAINS e Subestruturas Tóxicas/Reativas)
    toxic_alerts = []
    TOX_ALERTS = {
        "Quinona": "O=C1C=CC(=O)C=C1",
        "Catecol": "Oc1c(O)cccc1",
        "Epóxido (Anel reativo)": "C1OC1",
        "Haleto de Ácido": "C(=O)[Cl,Br,I]",
        "Aldeído Alifático": "[CH1](=O)",
        "Nitrogrupo": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
        "Hidrazina": "[NX3][NX3]",
        "Tiocarbonila": "C=S"
    }
    for name, smarts in TOX_ALERTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            toxic_alerts.append(name)

    # Veredito Final combinando físico-química e biologia
    has_risk = (hia_status == "Baixa Absorção") or (len(toxic_alerts) > 0)
    pass_filters = lipinski_pass and veber_pass and not has_risk

    return {
        "molecular_weight": round(mw, 2),
        "logp": round(logp, 2),
        "hydrogen_bond_donors": hbd,
        "hydrogen_bond_acceptors": hba,
        "tpsa": round(tpsa, 2),
        "rotatable_bonds": rotb,
        "formal_charge": charge,
        "lipinski_violations": lipinski_violations,
        "lipinski_pass": lipinski_pass,
        "veber_violations": veber_violations,
        "veber_pass": veber_pass,
        "hia_status": hia_status,
        "bbb_status": bbb_status,
        "pgp_status": pgp_status,
        "toxic_alerts": toxic_alerts,
        "pass_filters": pass_filters,
    }
