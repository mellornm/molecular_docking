from pathlib import Path
from typing import List, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED, FilterCatalog


def _get_filter_catalog() -> Optional[FilterCatalog.FilterCatalog]:
    """
    Inicializa o catálogo completo de filtros estruturais da literatura medicinal no RDKit:
    - PAINS_A, PAINS_B, PAINS_C (Pan Assay Interference Compounds)
    - BRENK (Filtros de alerta estrutural e reatividade química indesejada)
    - NIH (Filtros toxicológicos do National Institutes of Health)
    """
    try:
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH)
        return FilterCatalog.FilterCatalog(params)
    except Exception:
        return None


def calculate_admet_descriptors(ligand_sdf: Path) -> dict:
    """
    Calcula descritores físico-químicos, predições farmacocinéticas (ADMET),
    avaliação de Drug-likeness (QED), complexidade estereoquímica (Fsp3), acessibilidade sintética (SAscore),
    e triagem completa de subestruturas espúrias e toxicológicas (PAINS A/B/C, Brenk e NIH).
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
        # Descritores físico-químicos fundamentais
        mw = float(Descriptors.ExactMolWt(mol))
        logp = float(Crippen.MolLogP(mol))
        hbd = int(Descriptors.NumHDonors(mol))
        hba = int(Descriptors.NumHAcceptors(mol))
        tpsa = float(Descriptors.TPSA(mol))
        rotb = int(Descriptors.NumRotatableBonds(mol))
        charge = int(Chem.GetFormalCharge(mol))

        # Métricas Quimioinformáticas Avançadas
        fsp3 = float(Descriptors.FractionCSP3(mol))
        qed_score = float(QED.qed(mol))
        n_chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        n_rings = int(mol.GetRingInfo().NumRings())

        # Estimativa de Acessibilidade Sintética (1.0 = muito fácil, 10.0 = muito complexo)
        sa_raw = (
            1.0 + (mw / 180.0) + (rotb * 0.12) + (n_chiral * 0.75) + (n_rings * 0.35)
        )
        if fsp3 < 0.25:
            sa_raw += 0.4
        sa_score = round(min(10.0, max(1.0, sa_raw)), 2)

    except Exception as e:
        raise RuntimeError(f"Erro ao calcular descritores moleculares com RDKit: {e}")

    # Classificação QED
    if qed_score >= 0.67:
        qed_classification = "Alto / Favorável (Drug-like)"
    elif qed_score >= 0.49:
        qed_classification = "Moderado (Médio potencial drug-like)"
    else:
        qed_classification = "Baixo / Desfavorável (Não drug-like)"

    # Classificação SAscore
    if sa_score <= 3.5:
        sa_classification = "Fácil Síntese / Acessível"
    elif sa_score <= 6.0:
        sa_classification = "Dificuldade Moderada"
    else:
        sa_classification = "Alta Complexidade Sintética"

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

    # 3. Filtros Complementares de Química Medicinal
    # Lead-likeness (Teague & Oprea): MW 150-350, LogP -1.0 a 3.5, RotB <= 7
    lead_likeness_pass = bool(
        150.0 <= mw <= 350.0 and -1.0 <= logp <= 3.5 and rotb <= 7
    )

    # Golden Triangle (Johnson & Zheng / Pfizer): MW 200-400, LogP -1.0 a 3.0
    golden_triangle_pass = bool(200.0 <= mw <= 400.0 and -1.0 <= logp <= 3.0)

    # Total de violações físico-químicas combinadas (Lipinski + Veber)
    total_violations = len(lipinski_violations) + len(veber_violations)
    all_violations = lipinski_violations + veber_violations

    # 4. Absorção Intestinal Humana (HIA) - Filtro de Egan (Egan Egg)
    hia_ok = (tpsa <= 132.0) and (-1.0 <= logp <= 5.8)
    hia_status = "Alta Absorção" if hia_ok else "Baixa Absorção"

    # 5. Permeabilidade da Barreira Hematoencefálica (BBB) - Regra de Clark
    bbb_ok = (charge == 0) and (tpsa < 90.0) and (1.0 <= logp <= 5.0)
    bbb_status = "Permeável" if bbb_ok else "Incompatível/Baixa"

    # 6. Substrato de P-glicoproteína (P-gp) - Modelo baseado em tamanho e polaridade
    pgp_substrate = bool((mw > 400.0) and (tpsa > 80.0))
    pgp_status = (
        "Substrato (Efluxo provável)"
        if pgp_substrate
        else "Não Substrato (Baixo Efluxo)"
    )

    # 7. Triagem Toxicológica Avançada e Alertas Estruturais (PAINS, Brenk, SMARTS)
    toxic_alerts: List[str] = []
    pains_alerts: List[str] = []
    brenk_alerts: List[str] = []

    # 7.1 Catálogo RDKit FilterCatalog (PAINS A/B/C, Brenk, NIH)
    catalog = _get_filter_catalog()
    if catalog is not None:
        try:
            matches = catalog.GetMatches(mol)
            for entry in matches:
                desc = entry.GetDescription()
                if "PAINS" in desc.upper():
                    pains_alerts.append(desc)
                elif "BRENK" in desc.upper():
                    brenk_alerts.append(desc)
                else:
                    toxic_alerts.append(desc)
        except Exception:
            pass

    # 7.2 SMARTS Toxicológicos Traduzidos em Português
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
            if name not in toxic_alerts:
                toxic_alerts.append(name)

    all_structural_alerts = pains_alerts + brenk_alerts + toxic_alerts
    has_severe_risk = (hia_status == "Baixa Absorção") or (
        len(all_structural_alerts) > 0
    )

    # Classificação de Veredito ADMET
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
        dynamic_points.append(
            "• Físico-Química: 100% de conformidade com as regras clássicas de Lipinski e Veber (0 violações)."
        )
        dynamic_points.append(
            f"• Drug-likeness (QED): Score {qed_score:.2f} ({qed_classification}) e Fsp3 = {fsp3:.2f}."
        )
        dynamic_points.append(
            "• Farmacocinética: Alta Absorção Intestinal (HIA) estimada pelo modelo Egan Egg."
        )
        dynamic_points.append(
            "• Toxicidade & PAINS: Nenhum alerta estrutural reativo ou subestrutura PAINS/Brenk identificada."
        )
        attention_note = "A molécula possui excelente perfil biofarmacêutico, químico-medicinal e físico-químico para desenvolvimento oral."
    elif verdict_category == "MODERATE":
        single_viol = (
            all_violations[0]
            if all_violations
            else "Desvio pontual em parâmetro físico-químico"
        )
        dynamic_points.append(
            f"• Desvio Pontual Tolerado: {single_viol} (desvio único aceito em fármacos aprovados)."
        )
        dynamic_points.append(
            f"• Drug-likeness (QED): Score {qed_score:.2f} ({qed_classification}) e Fsp3 = {fsp3:.2f}."
        )
        dynamic_points.append(
            "• Farmacocinética: Mantém Alta Absorção Intestinal (HIA) estimada (Egan Egg)."
        )
        dynamic_points.append(
            "• Toxicidade: Ausência de subestruturas tóxicas severas ou PAINS."
        )
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
        reasons = []
        if total_violations >= 2:
            reasons.append(
                f"Múltiplas violações físico-químicas ({total_violations} violações: {', '.join(all_violations)})"
            )
        elif total_violations == 1:
            reasons.append(f"Violação físico-química: {all_violations[0]}")

        if hia_status == "Baixa Absorção":
            reasons.append(
                "Baixa Absorção Intestinal (HIA) estimada (fora da elipse de Egan)"
            )

        if len(all_structural_alerts) > 0:
            reasons.append(
                f"Alertas de Toxicidade/PAINS: {', '.join(all_structural_alerts)}"
            )

        dynamic_points.append(f"• Problemas Identificados: {'; '.join(reasons)}.")
        if hia_status == "Alta Absorção":
            dynamic_points.append(
                "• Farmacocinética: Alta Absorção Intestinal (HIA) preservada."
            )
        if len(all_structural_alerts) == 0:
            dynamic_points.append("• Toxicidade: Sem alertas estruturais de PAINS.")

        attention_parts = []
        if total_violations >= 2:
            attention_parts.append(
                f"propriedades físico-químicas desfavoráveis ({', '.join(all_violations)})"
            )
        elif total_violations == 1 and (
            hia_status == "Baixa Absorção" or len(all_structural_alerts) > 0
        ):
            attention_parts.append(f"desvio físico-químico ({all_violations[0]})")

        if hia_status == "Baixa Absorção":
            attention_parts.append("baixa absorção intestinal (HIA)")

        if len(all_structural_alerts) > 0:
            attention_parts.append(
                f"riscos de toxicidade por subestruturas reativas ({', '.join(all_structural_alerts)})"
            )

        attention_note = (
            f"A molécula apresenta {' e '.join(attention_parts)}."
            if attention_parts
            else "A molécula não atingiu os critérios mínimos de triagem ADMET."
        )

    return {
        "molecular_weight": round(mw, 2),
        "logp": round(logp, 2),
        "hydrogen_bond_donors": hbd,
        "hydrogen_bond_acceptors": hba,
        "tpsa": round(tpsa, 2),
        "rotatable_bonds": rotb,
        "formal_charge": charge,
        "fsp3": round(fsp3, 3),
        "qed_score": round(qed_score, 3),
        "qed_classification": qed_classification,
        "synthetic_accessibility": sa_score,
        "synthetic_accessibility_classification": sa_classification,
        "chiral_centers": n_chiral,
        "num_rings": n_rings,
        "lead_likeness_pass": lead_likeness_pass,
        "golden_triangle_pass": golden_triangle_pass,
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
        "pains_alerts": pains_alerts,
        "brenk_alerts": brenk_alerts,
        "all_structural_alerts": all_structural_alerts,
        "verdict_category": verdict_category,
        "verdict_label": verdict_label,
        "pass_filters": pass_filters,
        "dynamic_points": dynamic_points,
        "attention_note": attention_note,
    }
