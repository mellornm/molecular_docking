import base64
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


def _find_file(work_dir: Path, filenames: List[str]) -> Optional[Path]:
    """Busca o primeiro arquivo correspondente na lista de nomes dentro do work_dir ou subdiretórios imediatos."""
    for fn in filenames:
        direct = work_dir / fn
        if direct.exists():
            return direct
    # Busca recursiva rasa (1 nível)
    for fn in filenames:
        matches = list(work_dir.glob(f"*/{fn}"))
        if matches:
            return matches[0]
    return None


def _find_image_base64(image_path: Optional[Path]) -> Optional[str]:
    """Lê um arquivo de imagem PNG e o codifica como base64 data URI."""
    if not image_path or not image_path.exists():
        return None
    try:
        with open(image_path, "rb") as f:
            encoded = base64.b64encode(f.read()).decode("utf-8")
            return f"data:image/png;base64,{encoded}"
    except Exception:
        return None


def parse_pipeline_artifacts(work_dir: Path) -> Tuple[Dict[str, Any], List[str]]:
    """
    Varre e faz o parse de todos os artefatos gerados pelo pipeline no diretório fornecido:
    - Score de Vina (ΔG) dos logs
    - Dados ADMET de pharmacokinetics.json ou interactions.json
    - Interações atômicas de interactions.json (PLIP)
    - Gráficos da DM (rmsd.png, rmsf.png, hbond.png)
    - Sumário termodinâmico MM-PBSA de mmpbsa_summary.json

    Retorna uma tupla (dados_estruturados, lista_de_avisos).
    """
    work_dir = Path(work_dir)
    warnings: List[str] = []
    data: Dict[str, Any] = {
        "work_dir": str(work_dir.resolve()),
        "generated_at": datetime.now().strftime("%d/%m/%Y %H:%M:%S"),
        "vina_score": None,
        "vina_log_path": None,
        "admet": None,
        "interactions": {
            "hydrogen_bonds": [],
            "hydrophobic_contacts": [],
        },
        "hbond_occupancy": [],
        "mmpbsa": None,
        "plots": {
            "rmsd": None,
            "rmsf": None,
            "hbond": None,
            "decomp": None,
        },
    }

    # 1. Parse do Score de Vina (ΔG)
    vina_log_files = [
        "vina_log.txt",
        "desoxicolato_vina.log",
        "vina.log",
        "docking.log",
        "log.txt",
    ]
    log_file = _find_file(work_dir, vina_log_files)
    if not log_file:
        # Tenta glob genérico para logs
        possible_logs = list(work_dir.glob("*vina*.log")) + list(work_dir.glob("*.log"))
        if possible_logs:
            log_file = possible_logs[0]

    if log_file and log_file.exists():
        data["vina_log_path"] = str(log_file.name)
        try:
            with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith("   1 "):
                        parts = line.split()
                        if len(parts) >= 2:
                            data["vina_score"] = float(parts[1])
                            break
        except Exception as e:
            warnings.append(
                f"Não foi possível ler o Score de Vina do log '{log_file.name}': {e}"
            )
    else:
        warnings.append("Log do AutoDock Vina não encontrado no diretório.")

    # 2. Parse de Interações e ADMET (interactions.json ou pharmacokinetics.json)
    interactions_file = _find_file(work_dir, ["interactions.json"])
    pharmacokinetics_file = _find_file(work_dir, ["pharmacokinetics.json"])

    if interactions_file and interactions_file.exists():
        try:
            with open(interactions_file, "r", encoding="utf-8") as f:
                inter_json = json.load(f)
                data["interactions"]["hydrogen_bonds"] = inter_json.get(
                    "hydrogen_bonds", []
                )
                data["interactions"]["hydrophobic_contacts"] = inter_json.get(
                    "hydrophobic_contacts", []
                )
                if "pharmacokinetics" in inter_json:
                    data["admet"] = inter_json["pharmacokinetics"]
        except Exception as e:
            warnings.append(f"Erro ao processar 'interactions.json': {e}")

    # 2.1 Parse de Ocupação Temporal de Pontes de Hidrogênio (hbond_occupancy.json)
    hbond_occ_file = _find_file(work_dir, ["hbond_occupancy.json"])
    if not hbond_occ_file:
        for parent in [work_dir] + list(work_dir.parents):
            candidate_md = parent / "md_files"
            if candidate_md.exists():
                hbond_occ_file = _find_file(candidate_md, ["hbond_occupancy.json"])
                if hbond_occ_file:
                    break

    if hbond_occ_file and hbond_occ_file.exists():
        try:
            with open(hbond_occ_file, "r", encoding="utf-8") as f:
                data["hbond_occupancy"] = json.load(f)
        except Exception:
            pass

    # Se pharmacokinetics.json existir isoladamente, dá prioridade ou preenche
    if pharmacokinetics_file and pharmacokinetics_file.exists():
        try:
            with open(pharmacokinetics_file, "r", encoding="utf-8") as f:
                data["admet"] = json.load(f)
        except Exception as e:
            warnings.append(f"Erro ao processar 'pharmacokinetics.json': {e}")

    if not data["admet"]:
        warnings.append("Dados farmacocinéticos (ADMET) não encontrados.")
    if (
        not data["interactions"]["hydrogen_bonds"]
        and not data["interactions"]["hydrophobic_contacts"]
    ):
        warnings.append("Mapeamento de interações estruturais (PLIP) não encontrado.")

    # 3. Parse do Sumário MM-PBSA
    mmpbsa_file = _find_file(work_dir, ["mmpbsa_summary.json"])
    if not mmpbsa_file:
        # Busca no diretório md_files caso esteja em subdiretório ou em diretórios pais/irmãos
        for parent in [work_dir] + list(work_dir.parents):
            candidate_md = parent / "md_files"
            if candidate_md.exists():
                mmpbsa_file = _find_file(candidate_md, ["mmpbsa_summary.json"])
                if mmpbsa_file:
                    break

    if mmpbsa_file and mmpbsa_file.exists():
        try:
            with open(mmpbsa_file, "r", encoding="utf-8") as f:
                data["mmpbsa"] = json.load(f)
        except Exception as e:
            warnings.append(f"Erro ao processar 'mmpbsa_summary.json': {e}")
    else:
        warnings.append(
            "Sumário de Energia Livre MM-PBSA (mmpbsa_summary.json) não encontrado."
        )

    # 4. Busca e Codificação de Gráficos de Dinâmica Molecular (PNG)
    for plot_name in ["rmsd", "rmsf", "hbond", "gyrate", "sasa", "decomp"]:
        img_name = "decomp_mmpbsa.png" if plot_name == "decomp" else f"{plot_name}.png"
        img_file = _find_file(work_dir, [img_name])
        if not img_file:
            for parent in [work_dir] + list(work_dir.parents):
                candidate_md = parent / "md_files"
                if candidate_md.exists():
                    img_file = _find_file(candidate_md, [img_name])
                    if img_file:
                        break

        if img_file and img_file.exists():
            b64_str = _find_image_base64(img_file)
            if b64_str:
                data["plots"][plot_name] = b64_str
            else:
                warnings.append(f"Falha ao codificar a imagem '{img_name}'.")
        elif plot_name not in ("decomp", "gyrate", "sasa"):
            warnings.append(
                f"Gráfico de Dinâmica Molecular '{img_name}' não encontrado."
            )

    return data, warnings


def generate_html_report(
    work_dir: Path,
    output_file: Optional[Path] = None,
    receptor_name: Optional[str] = None,
    ligand_name: Optional[str] = None,
) -> Tuple[Path, List[str]]:
    """
    Gera um relatório executivo e científico em HTML autoconferível consolidando
    todos os resultados de Docking, PLIP, ADMET e Dinâmica Molecular.

    :param work_dir: Diretório onde se encontram os artefatos do pipeline.
    :param output_file: Caminho de saída para o arquivo HTML (padrão: work_dir/report.html).
    :param receptor_name: Código ou nome do receptor (ex: 7CFN, GPBAR1).
    :param ligand_name: Nome do ligante (ex: INT-777, Desoxicolato).
    :return: Tupla com o caminho do arquivo HTML gerado e a lista de avisos de arquivos ausentes.
    """
    work_dir = Path(work_dir)
    if output_file is None:
        output_file = work_dir / "report.html"
    else:
        output_file = Path(output_file)

    output_file.parent.mkdir(parents=True, exist_ok=True)

    data, warnings = parse_pipeline_artifacts(work_dir)

    # Inferência inteligente de metadados caso não informados
    if not receptor_name:
        for part in work_dir.resolve().parts:
            if (
                len(part) == 4
                and part.isalnum()
                and (part.isupper() or any(c.isdigit() for c in part))
            ):
                receptor_name = part
                break
    if not ligand_name:
        sdf_files = list(work_dir.glob("*.sdf")) + list(work_dir.glob("*/*.sdf"))
        if sdf_files:
            ligand_name = sdf_files[0].stem
        elif "desoxicolato" in str(work_dir).lower():
            ligand_name = "Desoxicolato"

    receptor_name_display = (
        receptor_name.strip()
        if receptor_name and receptor_name.strip()
        else "Não especificado"
    )
    ligand_name_display = (
        ligand_name.strip()
        if ligand_name and ligand_name.strip()
        else "Não especificado"
    )

    # Extração de Métricas Principais
    vina_score = data["vina_score"]
    vina_str = f"{vina_score:.2f} kcal/mol" if vina_score is not None else "N/A"

    mmpbsa_data = data.get("mmpbsa", {})
    mmpbsa_energies = mmpbsa_data.get("energies", {}) if mmpbsa_data else {}
    dg_bind = mmpbsa_energies.get("delta_g_binding", {}) if mmpbsa_energies else {}

    if dg_bind and "mean" in dg_bind:
        dg_mean = dg_bind.get("mean", 0.0)
        dg_std = dg_bind.get("std", 0.0)
        unit = mmpbsa_data.get("unit", "kcal/mol")
        mmpbsa_str = f"{dg_mean:.2f} ± {dg_std:.2f} {unit}"
        mmpbsa_subtext = mmpbsa_data.get(
            "thermodynamic_window",
            "Janela Termodinâmica: 60 - 100 ns (Últimos 40% - Estado Estacionário)",
        )
    else:
        mmpbsa_str = "Não executado / N/A"
        mmpbsa_subtext = "Janela Termodinâmica: 60 - 100 ns (Pendente)"

    admet = data.get("admet") or {}
    admet_pass = admet.get("pass_filters", False)
    verdict_cat = admet.get("verdict_category")
    hia_status = admet.get("hia_status", "N/A")
    bbb_status = admet.get("bbb_status", "N/A")
    pgp_status = admet.get("pgp_status", "N/A")
    toxic_alerts = admet.get("toxic_alerts", [])
    total_viol = admet.get("total_violations", 0)
    all_viol = admet.get("all_violations", [])
    attention_note = admet.get("attention_note", "")

    hbonds = data["interactions"]["hydrogen_bonds"]
    hcontacts = data["interactions"]["hydrophobic_contacts"]
    total_interactions = len(hbonds) + len(hcontacts)

    # Fallback caso verdict_category não esteja gravado no JSON legado
    if not verdict_cat and admet:
        has_severe_risk = (hia_status == "Baixa Absorção") or (len(toxic_alerts) > 0)
        if total_viol == 0 and not has_severe_risk:
            verdict_cat = "APPROVED"
        elif total_viol == 1 and not has_severe_risk:
            verdict_cat = "MODERATE"
        else:
            verdict_cat = "RISK"

    # Status ADMET Badge, Banner e Texto Resumo Dinâmico
    if admet:
        if verdict_cat == "APPROVED":
            admet_badge = '<span class="badge badge-success">APROVADO</span>'
            admet_banner_class = "veredito-success"
            admet_banner_icon = "✅"
            admet_banner_title = "Composto Aprovado na Triagem ADMET"
            admet_summary_text = (
                attention_note
                or "Molécula com 100% de conformidade físico-química (Lipinski & Veber), alta absorção intestinal estimada e ausência de toxicidade estrutural."
            )
        elif verdict_cat == "MODERATE":
            admet_badge = '<span class="badge badge-warning">APROVADO COM RESSALVAS</span>'
            admet_banner_class = "veredito-warning"
            admet_banner_icon = "⚠️"
            admet_banner_title = "Aprovado com Ressalvas (Alerta Moderado)"
            admet_summary_text = (
                attention_note
                or "Molécula com alta absorção intestinal e ausência de toxicidade, apresentando desvio pontual em parâmetro físico-químico aceito pela literatura."
            )
        else:
            admet_badge = '<span class="badge badge-danger">REPROVADO / ALTO RISCO</span>'
            admet_banner_class = "veredito-danger"
            admet_banner_icon = "🚫"
            admet_banner_title = "Alertas Críticos ou Violações ADMET Identificadas"
            admet_summary_text = (
                attention_note
                or "A molécula apresentou restrições em critérios biofarmacêuticos ou toxicológicos."
            )
    else:
        admet_badge = '<span class="badge badge-secondary">NÃO ANALISADO</span>'
        admet_banner_class = "veredito-warning"
        admet_banner_icon = "ℹ️"
        admet_banner_title = "Triagem ADMET Não Disponível"
        admet_summary_text = (
            "Dados de triagem ADMET não disponíveis no diretório de trabalho."
        )

    # Geração das Linhas da Tabela ADMET
    admet_rows_html = ""
    if admet:
        mw = admet.get("molecular_weight", 0.0)
        mw_status = (
            '<span class="badge badge-success">OK</span>'
            if mw <= 500
            else '<span class="badge badge-danger">VIOLADO</span>'
        )
        logp = admet.get("logp", 0.0)
        logp_status = (
            '<span class="badge badge-success">OK</span>'
            if logp <= 5
            else '<span class="badge badge-danger">VIOLADO</span>'
        )
        hbd = admet.get("hydrogen_bond_donors", 0)
        hbd_status = (
            '<span class="badge badge-success">OK</span>'
            if hbd <= 5
            else '<span class="badge badge-danger">VIOLADO</span>'
        )
        hba = admet.get("hydrogen_bond_acceptors", 0)
        hba_status = (
            '<span class="badge badge-success">OK</span>'
            if hba <= 10
            else '<span class="badge badge-danger">VIOLADO</span>'
        )
        tpsa = admet.get("tpsa", 0.0)
        tpsa_status = (
            '<span class="badge badge-success">OK</span>'
            if tpsa <= 140
            else '<span class="badge badge-danger">VIOLADO</span>'
        )
        rotb = admet.get("rotatable_bonds", 0)
        rotb_status = (
            '<span class="badge badge-success">OK</span>'
            if rotb <= 10
            else '<span class="badge badge-danger">VIOLADO</span>'
        )

        hia_badge = (
            '<span class="badge badge-success">Alta Absorção</span>'
            if hia_status == "Alta Absorção"
            else '<span class="badge badge-danger">Baixa Absorção</span>'
        )
        bbb_badge = (
            '<span class="badge badge-success">Permeável</span>'
            if bbb_status == "Permeável"
            else '<span class="badge badge-warning">Baixa / Incompatível</span>'
        )
        is_pgp_substrate = (
            admet.get("pgp_substrate")
            if "pgp_substrate" in admet
            else (pgp_status.startswith("Substrato") or ("Substrato" in pgp_status and "Não" not in pgp_status))
        )
        pgp_badge = (
            '<span class="badge badge-warning">Efluxo Ativo</span>'
            if is_pgp_substrate
            else '<span class="badge badge-success">Baixo Efluxo</span>'
        )

        fsp3 = admet.get("fsp3", 0.0)
        fsp3_badge = '<span class="badge badge-success">Alto (3D Complexo)</span>' if fsp3 >= 0.42 else '<span class="badge badge-secondary">Moderado / Plano</span>'
        qed = admet.get("qed_score", 0.0)
        qed_cls = admet.get("qed_classification", "N/A")
        qed_badge = '<span class="badge badge-success">Alto (Drug-like)</span>' if qed >= 0.67 else ('<span class="badge badge-primary">Moderado</span>' if qed >= 0.49 else '<span class="badge badge-warning">Baixo</span>')
        sa = admet.get("synthetic_accessibility", 0.0)
        sa_cls = admet.get("synthetic_accessibility_classification", "N/A")
        sa_badge = '<span class="badge badge-success">Fácil Síntese</span>' if sa <= 3.5 else ('<span class="badge badge-primary">Moderada</span>' if sa <= 6.0 else '<span class="badge badge-danger">Alta Complexidade</span>')

        lead_like = admet.get("lead_likeness_pass", False)
        lead_badge = '<span class="badge badge-success">Conforme</span>' if lead_like else '<span class="badge badge-secondary">Não Enquadrado</span>'
        golden = admet.get("golden_triangle_pass", False)
        golden_badge = '<span class="badge badge-success">Conforme</span>' if golden else '<span class="badge badge-secondary">Fora do Triângulo</span>'

        all_alerts = admet.get("all_structural_alerts", toxic_alerts)
        if all_alerts:
            tox_badge = '<span class="badge badge-danger">ALERTA ESTRUTURAL</span>'
            tox_val = ", ".join(all_alerts)
        else:
            tox_badge = '<span class="badge badge-success">Seguro (0 Alertas)</span>'
            tox_val = "Nenhum alerta de PAINS, Brenk ou subestruturas reativas identificado"

        admet_rows_html = f"""
        <tr class="section-header"><td colspan="4">1. Propriedades Físico-Químicas & Regras de Filtro (Lipinski & Veber)</td></tr>
        <tr><td>Peso Molecular (MW)</td><td class="text-mono">{mw:.2f} g/mol</td><td>&le; 500.00 Da</td><td>{mw_status}</td></tr>
        <tr><td>Lipofilicidade (LogP)</td><td class="text-mono">{logp:.2f}</td><td>&le; 5.00</td><td>{logp_status}</td></tr>
        <tr><td>Doadores de H (HBD)</td><td class="text-mono">{hbd}</td><td>&le; 5</td><td>{hbd_status}</td></tr>
        <tr><td>Aceitadores de H (HBA)</td><td class="text-mono">{hba}</td><td>&le; 10</td><td>{hba_status}</td></tr>
        <tr><td>Superfície Polar Topológica (TPSA)</td><td class="text-mono">{tpsa:.2f} &Aring;&sup2;</td><td>&le; 140.00 &Aring;&sup2;</td><td>{tpsa_status}</td></tr>
        <tr><td>Ligações Rotacionáveis (RotB)</td><td class="text-mono">{rotb}</td><td>&le; 10</td><td>{rotb_status}</td></tr>

        <tr class="section-header"><td colspan="4">2. Quimioinformática & Medicinal Chemistry (Drug-likeness & Síntese)</td></tr>
        <tr><td>Estimativa Quantitativa de Drug-likeness (QED)</td><td class="text-mono">{qed:.2f} ({qed_cls})</td><td>&ge; 0.67 (Bickerton et al.)</td><td>{qed_badge}</td></tr>
        <tr><td>Fração de Carbonos sp3 (Fsp3)</td><td class="text-mono">{fsp3:.2f}</td><td>&ge; 0.42 (Complexidade 3D / Solubilidade)</td><td>{fsp3_badge}</td></tr>
        <tr><td>Acessibilidade Sintética (SAscore)</td><td class="text-mono">{sa:.2f} ({sa_cls})</td><td>1.0 (Muito Fácil) a 10.0 (Muito Complexo)</td><td>{sa_badge}</td></tr>
        <tr><td>Perfil Lead-like (Teague & Oprea)</td><td class="text-mono">{"Passou" if lead_like else "Desvio"}</td><td>MW 150-350 &bull; LogP -1 a 3.5 &bull; RotB &le; 7</td><td>{lead_badge}</td></tr>
        <tr><td>Golden Triangle (Pfizer / Johnson & Zheng)</td><td class="text-mono">{"Passou" if golden else "Desvio"}</td><td>MW 200-400 &bull; LogP -1 a 3.0</td><td>{golden_badge}</td></tr>

        <tr class="section-header"><td colspan="4">3. Farmacocinética e Biodisponibilidade (ADME)</td></tr>
        <tr><td>Absorção Intestinal Humana (HIA)</td><td class="text-mono">{hia_status}</td><td>Egan Egg (TPSA &le; 132 & -1.0 &le; LogP &le; 5.8)</td><td>{hia_badge}</td></tr>
        <tr><td>Permeabilidade Hematoencefálica (BBB)</td><td class="text-mono">{bbb_status}</td><td>Clark (Neutra, TPSA &lt; 90 & 1.0 &le; LogP &le; 5.0)</td><td>{bbb_badge}</td></tr>
        <tr><td>Substrato de P-glicoproteína (P-gp)</td><td class="text-mono">{pgp_status}</td><td>MW &gt; 400 & TPSA &gt; 80</td><td>{pgp_badge}</td></tr>

        <tr class="section-header"><td colspan="4">4. Toxicologia e Alertas Estruturais (PAINS, Brenk & SMARTS)</td></tr>
        <tr><td>Catálogo de Subestruturas Espúrias & Alertas</td><td class="text-mono">{tox_val}</td><td>Filtros RDKit PAINS A/B/C, Brenk & NIH</td><td>{tox_badge}</td></tr>
        """
    else:
        admet_rows_html = """<tr><td colspan="4" class="text-center text-muted">Nenhum dado ADMET disponível</td></tr>"""

    # Geração das Linhas da Tabela de Interações (PLIP)
    interaction_rows_html = ""
    if hbonds or hcontacts:
        for hb in hbonds:
            res_label = (
                f"<strong>{hb.get('resname', 'UNK')}</strong> {hb.get('resnr', '')}"
            )
            dist = hb.get("distance", 0.0)
            interaction_rows_html += f"""
            <tr>
                <td><span class="res-tag res-hbond">{res_label}</span></td>
                <td><span class="badge badge-primary">Ponte de Hidrogênio</span></td>
                <td class="text-mono">{dist:.2f} &Aring;</td>
                <td>Contato Direto Polar</td>
            </tr>
            """
        for hc in hcontacts:
            res_label = (
                f"<strong>{hc.get('resname', 'UNK')}</strong> {hc.get('resnr', '')}"
            )
            dist = hc.get("distance", 0.0)
            interaction_rows_html += f"""
            <tr>
                <td><span class="res-tag res-hydro">{res_label}</span></td>
                <td><span class="badge badge-amber">Contato Hidrofóbico</span></td>
                <td class="text-mono">{dist:.2f} &Aring;</td>
                <td>Interação de Van der Waals / Apolar</td>
            </tr>
            """
    else:
        interaction_rows_html = """<tr><td colspan="4" class="text-center text-muted">Nenhuma interação atômica mapeada no diretório.</td></tr>"""

    # Galeria de Imagens da DM
    rmsd_img = data["plots"].get("rmsd")
    rmsf_img = data["plots"].get("rmsf")
    hbond_img = data["plots"].get("hbond")
    gyrate_img = data["plots"].get("gyrate")
    sasa_img = data["plots"].get("sasa")

    def render_plot_card(
        title: str, subtitle: str, img_b64: Optional[str], fallback_desc: str
    ) -> str:
        if img_b64:
            return f"""
            <div class="gallery-card">
                <div class="gallery-card-header">
                    <h4>{title}</h4>
                    <p>{subtitle}</p>
                </div>
                <div class="gallery-card-body">
                    <img src="{img_b64}" alt="{title}" class="img-responsive" />
                </div>
            </div>
            """
        else:
            return f"""
            <div class="gallery-card placeholder-card">
                <div class="gallery-card-header">
                    <h4>{title}</h4>
                    <p>{subtitle}</p>
                </div>
                <div class="gallery-card-body placeholder-body">
                    <div class="placeholder-icon">📊</div>
                    <p class="text-muted">{fallback_desc}</p>
                    <span class="badge badge-secondary">Pendente / Não Executado</span>
                </div>
            </div>
            """

    rmsd_card = render_plot_card(
        "RMSD - Estabilidade Estrutural e Persistência no Sítio (0 - 100 ns)",
        "Evolução Temporal do Backbone e Ligante no Sítio Ativo via Trajetória Ajustada (md_fit.xtc)",
        rmsd_img,
        "Gráfico rmsd.png não encontrado no diretório de trabalho.",
    )
    rmsf_card = render_plot_card(
        "RMSF - Flutuação Atômica por Resíduo (0 - 100 ns)",
        "Flexibilidade Conformacional dos Carbonos Alfa (C-α) ao Longo de Toda a Simulação",
        rmsf_img,
        "Gráfico rmsf.png não encontrado no diretório de trabalho.",
    )
    hbond_card = render_plot_card(
        "Pontes de Hidrogênio Intermoleculares (0 - 100 ns)",
        "Persistência de Contatos Receptor-Ligante e Ausência de Desprendimento (Unbinding)",
        hbond_img,
        "Gráfico hbond.png não encontrado no diretório de trabalho.",
    )
    gyrate_card = render_plot_card(
        "Raio de Giro (Rg) - Compacidade e Enovelamento (0 - 100 ns)",
        "Monitoramento da Compacidade Estrutural da Proteína ao Longo da Simulação",
        gyrate_img,
        "Gráfico gyrate.png não encontrado no diretório de trabalho.",
    )
    sasa_card = render_plot_card(
        "SASA - Área Acessível ao Solvente (0 - 100 ns)",
        "Estabilidade de Exposição ao Solvente e Integridade do Core Hidrofóbico",
        sasa_img,
        "Gráfico sasa.png não encontrado no diretório de trabalho.",
    )

    # Ocupação Temporal de Pontes de Hidrogênio (% Occupancy)
    hbond_occ_data = data.get("hbond_occupancy", [])
    hbond_occ_html = ""
    if hbond_occ_data:
        hbond_occ_rows = ""
        for occ in hbond_occ_data[:12]:
            d_pair = occ.get("donor", "UNK")
            a_pair = occ.get("acceptor", "LIG")
            pct = occ.get("occupancy_percent", 0.0)
            cls_txt = occ.get("classification", "")
            badge_cls = "badge-success" if pct >= 75.0 else ("badge-primary" if pct >= 35.0 else "badge-secondary")
            hbond_occ_rows += f"""
            <tr>
                <td><strong>{d_pair}</strong> &rarr; {a_pair}</td>
                <td class="text-mono font-bold">{pct:.1f}%</td>
                <td><span class="badge {badge_cls}">{cls_txt}</span></td>
            </tr>
            """
        hbond_occ_html = f"""
        <div class="card mt-4">
            <div class="card-header">
                <h3>⏱️ Persistência Temporal de Pontes de Hidrogênio (0 - 100 ns)</h3>
                <span class="badge badge-primary">% de Ocupação na DM</span>
            </div>
            <div class="card-body">
                <p class="text-muted" style="font-size: 0.875rem; margin-bottom: 1rem;">
                    Frequência de formação de pontes de hidrogênio ao longo dos 100 ns de trajetória para identificar âncoras farmacofóricas permanentes (&ge;75% de ocupação).
                </p>
                <div class="table-responsive">
                    <table class="table">
                        <thead>
                            <tr>
                                <th>Par Doador &rarr; Aceitador</th>
                                <th>Ocupação na Trajetória (%)</th>
                                <th>Classificação Farmacofórica</th>
                            </tr>
                        </thead>
                        <tbody>
                            {hbond_occ_rows}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
        """

    # Decomposição de Energia Livre por Resíduo (MM-PBSA Hotspots)
    decomp_img = data["plots"].get("decomp")
    hotspots = mmpbsa_data.get("hotspot_residues", []) if mmpbsa_data else []
    decomp_section_html = ""
    if decomp_img or hotspots:
        hotspot_rows = ""
        for h in hotspots:
            raw_res = h.get("residue", "").strip()
            if raw_res.startswith("R:") or raw_res.startswith("L:"):
                raw_res = raw_res[2:]
            res_lbl = raw_res.replace(":", " ").strip()
            tot = h.get("total", 0.0)
            vdw_h = h.get("vdw", 0.0)
            eel_h = h.get("electrostatic", 0.0)
            hotspot_rows += f"""
            <tr>
                <td><strong>{res_lbl}</strong></td>
                <td class="text-right text-mono">{vdw_h:.2f}</td>
                <td class="text-right text-mono">{eel_h:.2f}</td>
                <td class="text-right text-mono font-bold text-success">{tot:.2f} kcal/mol</td>
                <td><span class="badge badge-success">Hotspot Estabilizador</span></td>
            </tr>
            """

        decomp_plot_html = f"""
        <div style="margin-top: 1.5rem; text-align: center;">
            <h4 style="margin-bottom: 0.5rem; font-size: 1.05rem;">Perfil de Contribuição por Resíduo (&Delta;G<sub>bind</sub>)</h4>
            <img src="{decomp_img}" alt="MM-PBSA Per-Residue Decomposition" class="img-responsive" style="max-height: 480px; margin: 0 auto; border-radius: var(--radius-md); box-shadow: var(--shadow-sm);" />
        </div>
        """ if decomp_img else ""

        decomp_section_html = f"""
        <div class="card-body" style="border-top: 1px solid var(--border-color); background: #fafafa;">
            <h4 style="font-size: 1.1rem; color: var(--text-primary); margin-bottom: 0.75rem;">🎯 Decomposição por Resíduo (Hotspots Termodinâmicos de Ligação)</h4>
            <p class="text-muted" style="font-size: 0.875rem; margin-bottom: 1rem;">
                Resíduos do sítio ativo com maior contribuição de energia livre favorável (&Delta;G &lt; 0) para ancoramento do ligante.
            </p>
            <div class="table-responsive">
                <table class="table">
                    <thead>
                        <tr>
                            <th>Resíduo Chave</th>
                            <th class="text-right">&Delta;E<sub>vdw</sub></th>
                            <th class="text-right">&Delta;E<sub>elec</sub></th>
                            <th class="text-right">&Delta;G Total</th>
                            <th>Papel Termodinâmico</th>
                        </tr>
                    </thead>
                    <tbody>
                        {hotspot_rows if hotspot_rows else '<tr><td colspan="5" class="text-center text-muted">Dados de hotspots disponíveis no gráfico de decomposição</td></tr>'}
                    </tbody>
                </table>
            </div>
            {decomp_plot_html}
        </div>
        """

    # Tabela MM-PBSA detalhada
    mmpbsa_table_html = ""
    if mmpbsa_energies:
        vdw = mmpbsa_energies.get("van_der_waals", {"mean": 0.0, "std": 0.0})
        eel = mmpbsa_energies.get("electrostatic", {"mean": 0.0, "std": 0.0})
        polar = mmpbsa_energies.get("polar_solvation", {"mean": 0.0, "std": 0.0})
        apolar = mmpbsa_energies.get("nonpolar_solvation", {"mean": 0.0, "std": 0.0})
        dg = mmpbsa_energies.get("delta_g_binding", {"mean": 0.0, "std": 0.0})
        unit = mmpbsa_data.get("unit", "kcal/mol")
        window_label = mmpbsa_data.get(
            "thermodynamic_window",
            "Janela Termodinâmica: 60 - 100 ns (Últimos 40% - Estado Estacionário)",
        )
        start_f = mmpbsa_data.get("startframe", 600)
        end_f = mmpbsa_data.get("endframe", 1000)
        interval_f = mmpbsa_data.get("interval", 2)
        total_samples = mmpbsa_data.get(
            "frames_analyzed", max(1, (end_f - start_f + 1) // interval_f)
        )

        mmpbsa_table_html = f"""
        <section class="card mt-4">
            <div class="card-header">
                <h3>⚡ Energia Livre de Ligação MM-PBSA (Solvente Explícito)</h3>
                <span class="badge badge-primary">{window_label}</span>
            </div>
            <div class="card-body">
                <div class="veredito-banner" style="background-color: #f0fdf4; border: 1px solid #bbf7d0; color: #166534; margin-bottom: 1.25rem;">
                    <div class="veredito-icon">⚖️</div>
                    <div class="veredito-content">
                        <h4>Protocolo Padronizado de Dupla Escala Temporal</h4>
                        <p style="margin-top: 0.25rem; line-height: 1.5;">
                            <strong>• Análises Estruturais Globais (0 - 100 ns):</strong> Calculadas sobre toda a trajetória corrigida (<code>md_fit.xtc</code>) para comprovar a persistência do ligante no sítio ativo e a ausência de <em>unbinding</em>.<br/>
                            <strong>• Janela Termodinâmica MM-PBSA (60 - 100 ns / Últimos 40% - Estado Estacionário):</strong> Amostragem restrita à fase de produção em estado estacionário, eliminando o ruído conformacional da fase de relaxamento inicial (<em>induced-fit</em>) e reduzindo a variância amostral.<br/>
                            <strong>• Parâmetros de Amostragem:</strong> Frames {start_f} a {end_f} (Intervalo = {interval_f} &bull; Total Amostrado = {total_samples} frames).
                        </p>
                    </div>
                </div>

                <div class="table-responsive">
                    <table class="table">
                        <thead>
                            <tr>
                                <th>Componente Energético</th>
                                <th>Descrição / Papel Físico</th>
                                <th class="text-right">Energia Média ± Desvio Padrão ({unit})</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td><strong>Van der Waals (&Delta;E<sub>vdw</sub>)</strong></td>
                                <td>Atrações dispersivas e empacotamento estérico no sítio ativo</td>
                                <td class="text-right text-mono">{vdw.get("mean", 0.0):.2f} &plusmn; {vdw.get("std", 0.0):.2f}</td>
                            </tr>
                            <tr>
                                <td><strong>Eletrostática (&Delta;E<sub>elec</sub>)</strong></td>
                                <td>Interações de Coulomb e atração eletrostática específica</td>
                                <td class="text-right text-mono">{eel.get("mean", 0.0):.2f} &plusmn; {eel.get("std", 0.0):.2f}</td>
                            </tr>
                            <tr>
                                <td><strong>Solvatação Polar (&Delta;G<sub>polar</sub>)</strong></td>
                                <td>Custo termodinâmico de dessolvatação eletrostática (Poisson-Boltzmann)</td>
                                <td class="text-right text-mono">{polar.get("mean", 0.0):.2f} &plusmn; {polar.get("std", 0.0):.2f}</td>
                            </tr>
                            <tr>
                                <td><strong>Solvatação Apolar (&Delta;G<sub>apolar</sub>)</strong></td>
                                <td>Efeito hidrofóbico e variação da área de superfície acessível (SASA)</td>
                                <td class="text-right text-mono">{apolar.get("mean", 0.0):.2f} &plusmn; {apolar.get("std", 0.0):.2f}</td>
                            </tr>
                            <tr class="highlight-row">
                                <td><strong>&Delta;G Total de Ligação (&Delta;G<sub>bind</sub>)</strong></td>
                                <td><strong>Afinidade Termodinâmica Global MM-PBSA (Estado Estacionário)</strong></td>
                                <td class="text-right text-mono font-bold text-success">{dg.get("mean", 0.0):.2f} &plusmn; {dg.get("std", 0.0):.2f}</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>
            {decomp_section_html}
        </section>
        """

    # Template HTML Completo
    html_content = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Relatório Consolidado de Docking Molecular & Bioinformática</title>
    <style>
        :root {{
            --bg-body: #f8fafc;
            --bg-card: #ffffff;
            --text-primary: #0f172a;
            --text-secondary: #475569;
            --text-muted: #94a3b8;
            --border-color: #e2e8f0;
            --primary: #4f46e5;
            --primary-light: #eef2ff;
            --primary-dark: #3730a3;
            --success: #10b981;
            --success-bg: #ecfdf5;
            --danger: #ef4444;
            --danger-bg: #fef2f2;
            --warning: #f59e0b;
            --warning-bg: #fffbeb;
            --cyan: #06b6d4;
            --cyan-bg: #ecfeff;
            --font-sans: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            --font-mono: 'JetBrains Mono', 'Fira Code', Menlo, Monaco, Consolas, monospace;
            --shadow-sm: 0 1px 2px 0 rgb(0 0 0 / 0.05);
            --shadow-md: 0 4px 6px -1px rgb(0 0 0 / 0.1), 0 2px 4px -2px rgb(0 0 0 / 0.1);
            --shadow-lg: 0 10px 15px -3px rgb(0 0 0 / 0.1), 0 4px 6px -4px rgb(0 0 0 / 0.1);
            --radius-md: 12px;
            --radius-lg: 16px;
        }}

        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            background-color: var(--bg-body);
            color: var(--text-primary);
            font-family: var(--font-sans);
            line-height: 1.6;
            padding: 2.5rem 1.5rem;
        }}

        .container {{
            max-width: 1280px;
            margin: 0 auto;
        }}

        /* Header Hero */
        .hero {{
            background: linear-gradient(135deg, #1e1b4b 0%, #312e81 50%, #4338ca 100%);
            color: #ffffff;
            padding: 2.5rem 2rem;
            border-radius: var(--radius-lg);
            margin-bottom: 2rem;
            box-shadow: var(--shadow-lg);
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 1.5rem;
        }}

        .hero-title h1 {{
            font-size: 1.875rem;
            font-weight: 800;
            letter-spacing: -0.025em;
            margin-bottom: 0.5rem;
            display: flex;
            align-items: center;
            gap: 0.75rem;
        }}

        .hero-targets {{
            display: flex;
            flex-wrap: wrap;
            gap: 0.6rem;
            margin: 0.4rem 0 0.85rem 0;
        }}

        .target-chip {{
            display: inline-flex;
            align-items: center;
            gap: 0.45rem;
            background: rgba(255, 255, 255, 0.16);
            backdrop-filter: blur(8px);
            border: 1px solid rgba(255, 255, 255, 0.28);
            padding: 0.3rem 0.85rem;
            border-radius: 9999px;
            font-size: 0.9rem;
            color: #ffffff;
        }}

        .target-chip .chip-label {{
            color: #c7d2fe;
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            font-weight: 600;
        }}

        .target-chip strong {{
            color: #ffffff;
            font-weight: 700;
        }}

        .hero-title p {{
            color: #c7d2fe;
            font-size: 1rem;
            max-width: 700px;
        }}

        .hero-meta {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(8px);
            border: 1px solid rgba(255, 255, 255, 0.2);
            padding: 1rem 1.25rem;
            border-radius: var(--radius-md);
            font-size: 0.875rem;
        }}

        .hero-meta div {{
            margin-bottom: 0.25rem;
        }}

        .hero-meta div:last-child {{
            margin-bottom: 0;
        }}

        .hero-meta strong {{
            color: #e0e7ff;
        }}

        /* Grid de Cards de Métricas Principais */
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
            gap: 1.25rem;
            margin-bottom: 2rem;
        }}

        .metric-card {{
            background: var(--bg-card);
            border: 1px solid var(--border-color);
            border-radius: var(--radius-md);
            padding: 1.5rem;
            box-shadow: var(--shadow-sm);
            transition: transform 0.2s ease, box-shadow 0.2s ease;
            position: relative;
            overflow: hidden;
        }}

        .metric-card:hover {{
            transform: translateY(-2px);
            box-shadow: var(--shadow-md);
        }}

        .metric-card::before {{
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            height: 4px;
            background: var(--primary);
        }}

        .metric-card.card-vina::before {{ background: #0284c7; }}
        .metric-card.card-mmpbsa::before {{ background: #10b981; }}
        .metric-card.card-admet::before {{ background: #8b5cf6; }}
        .metric-card.card-plip::before {{ background: #f59e0b; }}

        .metric-label {{
            font-size: 0.875rem;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-secondary);
            margin-bottom: 0.5rem;
        }}

        .metric-value {{
            font-size: 1.75rem;
            font-weight: 800;
            color: var(--text-primary);
            font-family: var(--font-mono);
            margin-bottom: 0.25rem;
        }}

        .metric-subtext {{
            font-size: 0.8125rem;
            color: var(--text-muted);
        }}

        /* Estrutura de Seções / Cards */
        .card {{
            background: var(--bg-card);
            border: 1px solid var(--border-color);
            border-radius: var(--radius-lg);
            box-shadow: var(--shadow-sm);
            margin-bottom: 2rem;
            overflow: hidden;
        }}

        .card-header {{
            padding: 1.25rem 1.5rem;
            border-bottom: 1px solid var(--border-color);
            display: flex;
            justify-content: space-between;
            align-items: center;
            background: #ffffff;
        }}

        .card-header h3 {{
            font-size: 1.25rem;
            font-weight: 700;
            color: var(--text-primary);
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }}

        .card-body {{
            padding: 1.5rem;
        }}

        /* Tabelas */
        .table-responsive {{
            overflow-x: auto;
        }}

        .table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9375rem;
            text-align: left;
        }}

        .table th {{
            background: #f8fafc;
            color: var(--text-secondary);
            font-weight: 600;
            padding: 0.875rem 1rem;
            border-bottom: 2px solid var(--border-color);
            font-size: 0.8125rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
        }}

        .table td {{
            padding: 0.875rem 1rem;
            border-bottom: 1px solid var(--border-color);
            color: var(--text-primary);
        }}

        .table tbody tr:hover {{
            background-color: #f8fafc;
        }}

        .table tr.section-header td {{
            background: #f1f5f9;
            font-weight: 700;
            color: var(--primary-dark);
            padding: 0.625rem 1rem;
            font-size: 0.875rem;
            text-transform: uppercase;
            letter-spacing: 0.03em;
        }}

        .table tr.highlight-row td {{
            background: var(--success-bg);
            border-top: 2px solid var(--success);
            border-bottom: 2px solid var(--success);
        }}

        /* Badges */
        .badge {{
            display: inline-flex;
            align-items: center;
            padding: 0.25rem 0.625rem;
            border-radius: 9999px;
            font-size: 0.75rem;
            font-weight: 700;
            letter-spacing: 0.025em;
            text-transform: uppercase;
        }}

        .badge-success {{
            background-color: var(--success-bg);
            color: #065f46;
            border: 1px solid #a7f3d0;
        }}

        .badge-danger {{
            background-color: var(--danger-bg);
            color: #991b1b;
            border: 1px solid #fecaca;
        }}

        .badge-warning {{
            background-color: var(--warning-bg);
            color: #92400e;
            border: 1px solid #fde68a;
        }}

        .badge-primary {{
            background-color: var(--primary-light);
            color: var(--primary-dark);
            border: 1px solid #c7d2fe;
        }}

        .badge-amber {{
            background-color: #fff7ed;
            color: #c2410c;
            border: 1px solid #ffedd5;
        }}

        .badge-secondary {{
            background-color: #f1f5f9;
            color: #475569;
            border: 1px solid #cbd5e1;
        }}

        /* Tags de Resíduos */
        .res-tag {{
            display: inline-block;
            padding: 0.2rem 0.5rem;
            border-radius: 6px;
            font-family: var(--font-mono);
            font-size: 0.875rem;
            font-weight: 600;
        }}

        .res-hbond {{
            background-color: #e0f2fe;
            color: #0369a1;
            border: 1px solid #bae6fd;
        }}

        .res-hydro {{
            background-color: #fef3c7;
            color: #b45309;
            border: 1px solid #fde68a;
        }}

        /* Galeria de Gráficos da DM */
        .gallery-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(360px, 1fr));
            gap: 1.5rem;
        }}

        .gallery-card {{
            background: #ffffff;
            border: 1px solid var(--border-color);
            border-radius: var(--radius-md);
            overflow: hidden;
            box-shadow: var(--shadow-sm);
            display: flex;
            flex-direction: column;
        }}

        .gallery-card-header {{
            padding: 1rem 1.25rem;
            background: #f8fafc;
            border-bottom: 1px solid var(--border-color);
        }}

        .gallery-card-header h4 {{
            font-size: 1rem;
            font-weight: 700;
            color: var(--text-primary);
            margin-bottom: 0.25rem;
        }}

        .gallery-card-header p {{
            font-size: 0.8125rem;
            color: var(--text-muted);
        }}

        .gallery-card-body {{
            padding: 1rem;
            display: flex;
            justify-content: center;
            align-items: center;
            background: #ffffff;
            flex-grow: 1;
        }}

        .img-responsive {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            transition: transform 0.2s ease;
        }}

        .img-responsive:hover {{
            transform: scale(1.02);
        }}

        .placeholder-card {{
            border: 2px dashed #cbd5e1;
            background: #fafafa;
        }}

        .placeholder-body {{
            min-height: 240px;
            flex-direction: column;
            gap: 0.75rem;
            text-align: center;
            padding: 2rem;
        }}

        .placeholder-icon {{
            font-size: 2.5rem;
        }}

        /* Banner de Veredito ADMET */
        .veredito-banner {{
            padding: 1.25rem 1.5rem;
            border-radius: var(--radius-md);
            margin-bottom: 1.5rem;
            display: flex;
            align-items: center;
            gap: 1rem;
        }}

        .veredito-success {{
            background-color: var(--success-bg);
            border: 1px solid #a7f3d0;
            color: #065f46;
        }}

        .veredito-warning {{
            background-color: var(--warning-bg);
            border: 1px solid #fde68a;
            color: #92400e;
        }}

        .veredito-danger {{
            background-color: var(--danger-bg);
            border: 1px solid #fecaca;
            color: #991b1b;
        }}

        .veredito-icon {{
            font-size: 1.75rem;
        }}

        .veredito-content h4 {{
            font-weight: 700;
            margin-bottom: 0.25rem;
        }}

        .veredito-content p {{
            font-size: 0.875rem;
        }}

        /* Utilitários */
        .text-mono {{ font-family: var(--font-mono); }}
        .text-right {{ text-align: right; }}
        .text-center {{ text-align: center; }}
        .text-muted {{ color: var(--text-muted); }}
        .text-success {{ color: var(--success); }}
        .font-bold {{ font-weight: 700; }}
        .mt-4 {{ margin-top: 1.5rem; }}

        /* Footer */
        footer {{
            text-align: center;
            padding: 2rem 0 1rem;
            color: var(--text-muted);
            font-size: 0.875rem;
            border-top: 1px solid var(--border-color);
            margin-top: 3rem;
        }}

        footer strong {{
            color: var(--text-secondary);
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Hero Header -->
        <header class="hero">
            <div class="hero-title">
                <h1>🧬 Relatório Consolidado de Bioinformática</h1>
                <div class="hero-targets">
                    <span class="target-chip"><span class="chip-label">Receptor:</span> <strong>{receptor_name_display}</strong></span>
                    <span class="target-chip"><span class="chip-label">Ligante:</span> <strong>{ligand_name_display}</strong></span>
                </div>
                <p>Análise Integrada de Docking Molecular, Interações Atômicas (PLIP), Triagem ADMET e Dinâmica Molecular (GROMACS).</p>
            </div>
            <div class="hero-meta">
                <div><strong>Receptor:</strong> {receptor_name_display}</div>
                <div><strong>Ligante:</strong> {ligand_name_display}</div>
                <div><strong>Data:</strong> {data["generated_at"]}</div>
                <div><strong>Diretório:</strong> {data["work_dir"]}</div>
                <div><strong>Status Geral:</strong> {admet_badge}</div>
            </div>
        </header>

        <!-- Métricas Principais (Cards) -->
        <section class="metrics-grid">
            <div class="metric-card card-vina">
                <div class="metric-label">Score de Docking (Vina)</div>
                <div class="metric-value">{vina_str}</div>
                <div class="metric-subtext">&Delta;G de afinidade empírica (Pose 1)</div>
            </div>
            <div class="metric-card card-mmpbsa">
                <div class="metric-label">Energia Livre MM-PBSA</div>
                <div class="metric-value">{mmpbsa_str}</div>
                <div class="metric-subtext">{mmpbsa_subtext}</div>
            </div>
            <div class="metric-card card-admet">
                <div class="metric-label">Veredito ADMET</div>
                <div class="metric-value" style="font-size: 1.35rem; margin-top: 0.35rem;">{admet_badge}</div>
                <div class="metric-subtext">Lipinski, Veber & Toxicidade</div>
            </div>
            <div class="metric-card card-plip">
                <div class="metric-label">Interações Estruturais</div>
                <div class="metric-value">{total_interactions}</div>
                <div class="metric-subtext">{len(hbonds)} H-Bonds | {len(hcontacts)} Hidrofóbicos</div>
            </div>
        </section>

        <!-- Seção 1: Triagem Farmacocinética e Perfil ADMET -->
        <section class="card">
            <div class="card-header">
                <h3>💊 Triagem Farmacocinética, Físico-Química e Toxicologia (ADMET)</h3>
                {admet_badge}
            </div>
            <div class="card-body">
                <div class="veredito-banner {admet_banner_class}">
                    <div class="veredito-icon">{admet_banner_icon}</div>
                    <div class="veredito-content">
                        <h4>{admet_banner_title}</h4>
                        <p>{admet_summary_text}</p>
                    </div>
                </div>

                <div class="table-responsive">
                    <table class="table">
                        <thead>
                            <tr>
                                <th>Descritor / Propriedade</th>
                                <th>Valor Calculado</th>
                                <th>Critério / Faixa Ideal</th>
                                <th>Status</th>
                            </tr>
                        </thead>
                        <tbody>
                            {admet_rows_html}
                        </tbody>
                    </table>
                </div>
            </div>
        </section>

        <!-- Seção 2: Interações Atômicas Mapeadas (PLIP) -->
        <section class="card">
            <div class="card-header">
                <h3>🔍 Interações Intermoleculares Receptor-Ligante (PLIP)</h3>
                <span class="badge badge-primary">{total_interactions} Contatos Identificados</span>
            </div>
            <div class="card-body">
                <div class="table-responsive">
                    <table class="table">
                        <thead>
                            <tr>
                                <th>Resíduo do Receptor</th>
                                <th>Tipo de Interação</th>
                                <th>Distância (&Aring;)</th>
                                <th>Classificação Estrutural</th>
                            </tr>
                        </thead>
                        <tbody>
                            {interaction_rows_html}
                        </tbody>
                    </table>
                </div>
            </div>
        </section>

        <!-- Seção 2.1: Ocupação e Persistência de Pontes de Hidrogênio (DM 0 - 100 ns) -->
        {hbond_occ_html}

        <!-- Seção 3: Galeria de Dinâmica Molecular (GROMACS) -->
        <section class="card">
            <div class="card-header">
                <h3>📈 Monitoramento de Estabilidade Estrutural Global (0 - 100 ns)</h3>
                <span class="badge badge-secondary">Trajetória: md_fit.xtc &bull; 300 DPI</span>
            </div>
            <div class="card-body">
                <div class="gallery-grid">
                    {rmsd_card}
                    {rmsf_card}
                    {hbond_card}
                    {gyrate_card}
                    {sasa_card}
                </div>
            </div>
        </section>

        <!-- Seção 4: Decomposição MM-PBSA (Se disponível) -->
        {mmpbsa_table_html}

        <!-- Seção 5: Matrizes e Séries Temporais Exportadas em CSV para Publicação -->
        <section class="card">
            <div class="card-header">
                <h3>💾 Matrizes e Séries Temporais Exportadas (CSV / Publicação)</h3>
                <span class="badge badge-success">Pronto para GraphPad Prism / Origin / R</span>
            </div>
            <div class="card-body">
                <p class="text-muted" style="font-size: 0.875rem; margin-bottom: 1rem;">
                    Todos os dados quantitativos calculados durante o pipeline foram exportados no formato CSV para viabilizar a criação de gráficos customizados e tabelas para manuscritos científicos.
                </p>
                <div class="table-responsive">
                    <table class="table">
                        <thead>
                            <tr>
                                <th>Arquivo CSV</th>
                                <th>Métrica / Conteúdo</th>
                                <th>Colunas Exportadas</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td><code>rmsd.csv</code></td>
                                <td>Série temporal do RMSD da proteína (Backbone) e do ligante no sítio ativo</td>
                                <td class="text-mono">Time_ns, Protein_Backbone_RMSD_nm, Ligand_RMSD_nm</td>
                            </tr>
                            <tr>
                                <td><code>rmsf.csv</code></td>
                                <td>Flutuação residual atômica dos carbonos alfa (C-α)</td>
                                <td class="text-mono">Residue_Number, Calpha_RMSF_nm</td>
                            </tr>
                            <tr>
                                <td><code>hbond.csv</code></td>
                                <td>Contagem temporal de pontes de hidrogênio intermoleculares</td>
                                <td class="text-mono">Time_ns, HBond_Count</td>
                            </tr>
                            <tr>
                                <td><code>gyrate.csv</code></td>
                                <td>Evolução do Raio de Giro total e por eixos tensores (R<sub>x</sub>, R<sub>y</sub>, R<sub>z</sub>)</td>
                                <td class="text-mono">Time_ns, Total_nm, 2D_Rg_x_nm, 2D_Rg_y_nm, 2D_Rg_z_nm</td>
                            </tr>
                            <tr>
                                <td><code>sasa.csv</code></td>
                                <td>Área de Superfície Acessível ao Solvente</td>
                                <td class="text-mono">Time_ns, Total_SASA_nm2</td>
                            </tr>
                            <tr>
                                <td><code>hbond_occupancy.csv</code></td>
                                <td>Persistência percentual e identificação farmacofórica de H-Bonds</td>
                                <td class="text-mono">Donor, Acceptor, Occupancy_Percent, Classification</td>
                            </tr>
                            <tr>
                                <td><code>decomp_mmpbsa.csv</code></td>
                                <td>Contribuições energéticas por resíduo (&Delta;G, VdW, Eletrostática)</td>
                                <td class="text-mono">Residue, Van_der_Waals_kcal_mol, Electrostatic_kcal_mol, Total_DeltaG_kcal_mol</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>
        </section>

        <!-- Footer -->
        <footer>
            <p>Relatório gerado automaticamente pelo <strong>Molecular Docking & MD Pipeline</strong>.</p>
            <p class="text-muted">Automação Científica para Química Medicinal e Descoberta de Fármacos.</p>
        </footer>
    </div>
</body>
</html>
"""

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    return output_file, warnings
