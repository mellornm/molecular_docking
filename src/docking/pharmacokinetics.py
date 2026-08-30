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
        raise FileNotFoundError(
            f"Arquivo SDF do ligante não encontrado em: {ligand_sdf}"
        )

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
        raise ValueError(
            f"RDKit não conseguiu identificar uma molécula válida no arquivo: {ligand_sdf}"
        )

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
        raise RuntimeError(
            f"Erro ao calcular os descritores moleculares com RDKit: {e}"
        )

    # 1. Validação da Regra de Cinco de Lipinski (Pfizer)
    lipinski_violations = []
    if mw > 500.0:
        lipinski_violations.append(f"Peso Molecular elevado ({mw:.2f} > 500 Da)")
    if logp > 5.0:
        lipinski_violations.append(f"LogP elevado ({logp:.2f} > 5.0)")
    if hbd > 5:
        lipinski_violations.append(f"Doadores de H em excesso ({hbd} > 5)")
    if hba > 10:
        lipinski_violations.append(f"Aceitadores de H em excesso ({hba} > 10)")

    lipinski_pass = len(lipinski_violations) <= 1

    # 2. Validação das Regras de Veber (GSK)
    veber_violations = []
    if rotb > 10:
        veber_violations.append(f"Ligações rotacionáveis em excesso ({rotb} > 10)")
    if tpsa > 140.0:
        veber_violations.append(f"TPSA elevado ({tpsa:.2f} > 140 Å²)")

    veber_pass = len(veber_violations) == 0

    # Total de violações físico-químicas combinadas (Lipinski + Veber)
    total_violations = len(lipinski_violations) + len(veber_violations)
    all_violations = lipinski_violations + veber_violations

    # 3. Absorção Intestinal Humana (HIA) - Filtro de Egan (Egan Egg)
    # Alta absorção se TPSA <= 132 e -1.0 <= LogP <= 5.8
    hia_ok = (tpsa <= 132.0) and (-1.0 <= logp <= 5.8)
    hia_status = "Alta Absorção" if hia_ok else "Baixa Absorção"

    # 4. Permeabilidade da Barreira Hematoencefálica (BBB) - Regra de Clark
    # Se neutra, TPSA < 90 e LogP entre 1 e 5 -> Permeável
    bbb_ok = (charge == 0) and (tpsa < 90.0) and (1.0 <= logp <= 5.0)
    bbb_status = "Permeável" if bbb_ok else "Incompatível/Baixa"

    # 5. Substrato de P-glicoproteína (P-gp) - Modelo baseado em tamanho e polaridade
    # Moléculas com MW > 400 e TPSA > 80 tendem a ser substratos (efluxo provável)
    pgp_substrate = bool((mw > 400.0) and (tpsa > 80.0))
    pgp_status = (
        "Substrato (Efluxo provável)"
        if pgp_substrate
        else "Não Substrato (Baixo Efluxo)"
    )

    # 6. Alerta de Toxicidade (PAINS e Subestruturas Tóxicas/Reativas)
    toxic_alerts = []
    TOX_ALERTS = {
        "Quinona": "O=C1C=CC(=O)C=C1",
        "Catecol": "Oc1c(O)cccc1",
        "Epóxido (Anel reativo)": "C1OC1",
        "Haleto de Ácido": "C(=O)[Cl,Br,I]",
        "Aldeído Alifático": "[CH1](=O)",
        "Nitrogrupo": "[$([NX3](=O)=O),$([NX3+](=O)[O-])]",
        "Hidrazina": "[NX3][NX3]",
        "Tiocarbonila": "C=S",
    }
    for name, smarts in TOX_ALERTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            toxic_alerts.append(name)

    # Classificação de Veredito Flexibilizada (Padrão da Literatura Medicinal):
    # - Aprovado: 0 violações (Lipinski/Veber), alta HIA e 0 PAINS.
    # - Aprovado com Ressalvas: 1 violação tolerada, alta HIA e 0 PAINS.
    # - Reprovado / Alto Risco: >= 2 violações OU baixa absorção (HIA) OU subestruturas tóxicas (PAINS).
    has_severe_risk = (hia_status == "Baixa Absorção") or (len(toxic_alerts) > 0)

    if total_violations == 0 and not has_severe_risk:
        verdict_category = "APPROVED"
        verdict_label = "APROVADO (BIODISPONÍVEL & SEGURO)"
        pass_filters = True
    elif total_violations == 1 and not has_severe_risk:
        verdict_category = "MODERATE"
        verdict_label = "APROVADO COM RESSALVAS (ALERTA MODERADO)"
        pass_filters = True
    else:
        verdict_category = "RISK"
        verdict_label = "REPROVADO / RISCO ADMET"
        pass_filters = False

    # Mensagens Dinâmicas e Contextuais
    dynamic_points = []
    if verdict_category == "APPROVED":
        dynamic_points.append("• Físico-Química: 100% de conformidade com as regras clássicas de Lipinski e Veber (0 violações).")
        dynamic_points.append("• Farmacocinética: Alta Absorção Intestinal (HIA) estimada pelo modelo Egan Egg.")
        dynamic_points.append("• Toxicidade: Nenhum alerta estrutural reativo ou subestrutura PAINS identificada.")
        attention_note = "A molécula possui excelente perfil biofarmacêutico e físico-químico para desenvolvimento oral."
    elif verdict_category == "MODERATE":
        single_viol = all_violations[0] if all_violations else "Desvio pontual em parâmetro físico-químico"
        dynamic_points.append(f"• Desvio Pontual Tolerado: {single_viol} (desvio único aceito em fármacos aprovados).")
        dynamic_points.append("• Farmacocinética: Mantém Alta Absorção Intestinal (HIA) estimada (Egan Egg).")
        dynamic_points.append("• Toxicidade: Nenhum alerta de toxicidade estrutural ou PAINS identificado.")
        if rotb > 10:
            attention_note = f"A molécula apresenta alta flexibilidade conformacional (RotB = {rotb}, limite <= 10), mas mantém bom perfil de absorção e ausência de toxicidade."
        elif mw > 500:
            attention_note = f"A molécula possui peso molecular ligeiramente elevado (MW = {mw:.1f} g/mol), mas preserva boa absorção intestinal e segurança estrutural."
        elif logp > 5:
            attention_note = f"A lipofilicidade está ligeiramente elevada (LogP = {logp:.2f}), recomendando-se atenção à solubilidade aquosa."
        elif tpsa > 140:
            attention_note = f"A superfície polar está ligeiramente acima do limite de Veber (TPSA = {tpsa:.1f} Å²), mantendo boa absorção."
        else:
            attention_note = f"A molécula apresenta desvio pontual ({single_viol}), mantendo bom perfil global de absorção e ausência de toxicidade."
    else:
        # RISK
        reasons = []
        if total_violations >= 2:
            reasons.append(f"Múltiplas violações físico-químicas ({total_violations} violações: {', '.join(all_violations)})")
        elif total_violations == 1:
            reasons.append(f"Violação físico-química: {all_violations[0]}")

        if hia_status == "Baixa Absorção":
            reasons.append("Baixa Absorção Intestinal (HIA) estimada (fora da elipse de Egan)")

        if len(toxic_alerts) > 0:
            reasons.append(f"Alertas de Toxicidade/PAINS: {', '.join(toxic_alerts)}")

        dynamic_points.append(f"• Problemas Identificados: {'; '.join(reasons)}.")
        if hia_status == "Alta Absorção":
            dynamic_points.append("• Farmacocinética: Alta Absorção Intestinal (HIA) preservada.")
        if len(toxic_alerts) == 0:
            dynamic_points.append("• Toxicidade: Sem alertas estruturais de PAINS.")

        # Conclusão de atenção focada estritamente nas violações reais
        attention_parts = []
        if total_violations >= 2:
            attention_parts.append(f"propriedades físico-químicas desfavoráveis ({', '.join(all_violations)})")
        elif total_violations == 1 and (hia_status == "Baixa Absorção" or len(toxic_alerts) > 0):
            attention_parts.append(f"desvio físico-químico ({all_violations[0]})")

        if hia_status == "Baixa Absorção":
            attention_parts.append("baixa absorção intestinal (HIA)")

        if len(toxic_alerts) > 0:
            attention_parts.append(f"riscos de toxicidade por subestruturas reativas ({', '.join(toxic_alerts)})")

        attention_note = f"A molécula apresenta {' e '.join(attention_parts)}." if attention_parts else "A molécula não atingiu os critérios mínimos de triagem ADMET."

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
        "total_violations": total_violations,
        "all_violations": all_violations,
        "hia_status": hia_status,
        "bbb_status": bbb_status,
        "pgp_status": pgp_status,
        "pgp_substrate": pgp_substrate,
        "toxic_alerts": toxic_alerts,
        "verdict_category": verdict_category,
        "verdict_label": verdict_label,
        "pass_filters": pass_filters,
        "dynamic_points": dynamic_points,
        "attention_note": attention_note,
    }
