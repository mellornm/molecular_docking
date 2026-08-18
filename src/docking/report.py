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
        "mmpbsa": None,
        "plots": {
            "rmsd": None,
            "rmsf": None,
            "hbond": None,
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
            warnings.append(f"Não foi possível ler o Score de Vina do log '{log_file.name}': {e}")
    else:
        warnings.append("Log do AutoDock Vina não encontrado no diretório.")

    # 2. Parse de Interações e ADMET (interactions.json ou pharmacokinetics.json)
    interactions_file = _find_file(work_dir, ["interactions.json"])
    pharmacokinetics_file = _find_file(work_dir, ["pharmacokinetics.json"])

    if interactions_file and interactions_file.exists():
        try:
            with open(interactions_file, "r", encoding="utf-8") as f:
                inter_json = json.load(f)
                data["interactions"]["hydrogen_bonds"] = inter_json.get("hydrogen_bonds", [])
                data["interactions"]["hydrophobic_contacts"] = inter_json.get("hydrophobic_contacts", [])
                if "pharmacokinetics" in inter_json:
                    data["admet"] = inter_json["pharmacokinetics"]
        except Exception as e:
            warnings.append(f"Erro ao processar 'interactions.json': {e}")

    # Se pharmacokinetics.json existir isoladamente, dá prioridade ou preenche
    if pharmacokinetics_file and pharmacokinetics_file.exists():
        try:
            with open(pharmacokinetics_file, "r", encoding="utf-8") as f:
                data["admet"] = json.load(f)
        except Exception as e:
            warnings.append(f"Erro ao processar 'pharmacokinetics.json': {e}")

    if not data["admet"]:
        warnings.append("Dados farmacocinéticos (ADMET) não encontrados.")
    if not data["interactions"]["hydrogen_bonds"] and not data["interactions"]["hydrophobic_contacts"]:
        warnings.append("Mapeamento de interações estruturais (PLIP) não encontrado.")

    # 3. Parse do Sumário MM-PBSA
    mmpbsa_file = _find_file(work_dir, ["mmpbsa_summary.json"])
    if not mmpbsa_file:
        # Busca no diretório md_files caso esteja aninhado
        md_dir = work_dir / "md_files"
        if md_dir.exists():
            mmpbsa_file = _find_file(md_dir, ["mmpbsa_summary.json"])

    if mmpbsa_file and mmpbsa_file.exists():
        try:
            with open(mmpbsa_file, "r", encoding="utf-8") as f:
                data["mmpbsa"] = json.load(f)
        except Exception as e:
            warnings.append(f"Erro ao processar 'mmpbsa_summary.json': {e}")
    else:
        warnings.append("Sumário de Energia Livre MM-PBSA (mmpbsa_summary.json) não encontrado.")

    # 4. Busca e Codificação de Gráficos de Dinâmica Molecular (PNG)
    for plot_name in ["rmsd", "rmsf", "hbond"]:
        img_file = _find_file(work_dir, [f"{plot_name}.png"])
        if not img_file:
            md_dir = work_dir / "md_files"
            if md_dir.exists():
                img_file = _find_file(md_dir, [f"{plot_name}.png"])

        if img_file and img_file.exists():
            b64_str = _find_image_base64(img_file)
            if b64_str:
                data["plots"][plot_name] = b64_str
            else:
                warnings.append(f"Falha ao codificar a imagem '{plot_name}.png'.")
        else:
            warnings.append(f"Gráfico de Dinâmica Molecular '{plot_name}.png' não encontrado.")

    return data, warnings


def generate_html_report(work_dir: Path, output_file: Optional[Path] = None) -> Tuple[Path, List[str]]:
    """
    Gera um relatório executivo e científico em HTML autoconferível consolidando
    todos os resultados de Docking, PLIP, ADMET e Dinâmica Molecular.

    :param work_dir: Diretório onde se encontram os artefatos do pipeline.
    :param output_file: Caminho de saída para o arquivo HTML (padrão: work_dir/report.html).
    :return: Tupla com o caminho do arquivo HTML gerado e a lista de avisos de arquivos ausentes.
    """
    work_dir = Path(work_dir)
    if output_file is None:
        output_file = work_dir / "report.html"
    else:
        output_file = Path(output_file)

    output_file.parent.mkdir(parents=True, exist_ok=True)

    data, warnings = parse_pipeline_artifacts(work_dir)

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
    else:
        mmpbsa_str = "Não executado / N/A"

    admet = data.get("admet") or {}
    admet_pass = admet.get("pass_filters", False)
    hia_status = admet.get("hia_status", "N/A")
    bbb_status = admet.get("bbb_status", "N/A")
    pgp_status = admet.get("pgp_status", "N/A")
    toxic_alerts = admet.get("toxic_alerts", [])

    hbonds = data["interactions"]["hydrogen_bonds"]
    hcontacts = data["interactions"]["hydrophobic_contacts"]
    total_interactions = len(hbonds) + len(hcontacts)

    # Status ADMET Badge
    if admet:
        if admet_pass:
            admet_badge = '<span class="badge badge-success">APROVADO</span>'
            admet_summary_text = "Molécula com perfil físico-químico favorável, alta absorção e sem alertas de toxicidade."
        else:
            admet_badge = '<span class="badge badge-danger">REPROVADO / RISCO</span>'
            admet_summary_text = "A molécula apresentou violações de regras físico-químicas, baixa absorção ou alertas de toxicidade."
    else:
        admet_badge = '<span class="badge badge-secondary">NÃO ANALISADO</span>'
        admet_summary_text = "Dados de triagem ADMET não disponíveis no diretório de trabalho."

    # Geração das Linhas da Tabela ADMET
    admet_rows_html = ""
    if admet:
        mw = admet.get("molecular_weight", 0.0)
        mw_status = '<span class="badge badge-success">OK</span>' if mw <= 500 else '<span class="badge badge-danger">VIOLADO</span>'
        logp = admet.get("logp", 0.0)
        logp_status = '<span class="badge badge-success">OK</span>' if logp <= 5 else '<span class="badge badge-danger">VIOLADO</span>'
        hbd = admet.get("hydrogen_bond_donors", 0)
        hbd_status = '<span class="badge badge-success">OK</span>' if hbd <= 5 else '<span class="badge badge-danger">VIOLADO</span>'
        hba = admet.get("hydrogen_bond_acceptors", 0)
        hba_status = '<span class="badge badge-success">OK</span>' if hba <= 10 else '<span class="badge badge-danger">VIOLADO</span>'
        tpsa = admet.get("tpsa", 0.0)
        tpsa_status = '<span class="badge badge-success">OK</span>' if tpsa <= 140 else '<span class="badge badge-danger">VIOLADO</span>'
        rotb = admet.get("rotatable_bonds", 0)
        rotb_status = '<span class="badge badge-success">OK</span>' if rotb <= 10 else '<span class="badge badge-danger">VIOLADO</span>'

        hia_badge = '<span class="badge badge-success">Alta Absorção</span>' if hia_status == "Alta Absorção" else '<span class="badge badge-danger">Baixa Absorção</span>'
        bbb_badge = '<span class="badge badge-success">Permeável</span>' if bbb_status == "Permeável" else '<span class="badge badge-warning">Baixa / Incompatível</span>'
        pgp_badge = '<span class="badge badge-warning">Efluxo Ativo</span>' if "Substrato" in pgp_status else '<span class="badge badge-success">Baixo Efluxo</span>'

        if toxic_alerts:
            tox_badge = '<span class="badge badge-danger">ALERTA PAINS</span>'
            tox_val = ", ".join(toxic_alerts)
        else:
            tox_badge = '<span class="badge badge-success">Seguro</span>'
            tox_val = "Nenhum alerta estrutural identificado"

        admet_rows_html = f"""
        <tr class="section-header"><td colspan="4">1. Propriedades Físico-Químicas (Lipinski & Veber)</td></tr>
        <tr><td>Peso Molecular (MW)</td><td class="text-mono">{mw:.2f} g/mol</td><td>&le; 500.00 Da</td><td>{mw_status}</td></tr>
        <tr><td>Lipofilicidade (LogP)</td><td class="text-mono">{logp:.2f}</td><td>&le; 5.00</td><td>{logp_status}</td></tr>
        <tr><td>Doadores de H (HBD)</td><td class="text-mono">{hbd}</td><td>&le; 5</td><td>{hbd_status}</td></tr>
        <tr><td>Aceitadores de H (HBA)</td><td class="text-mono">{hba}</td><td>&le; 10</td><td>{hba_status}</td></tr>
        <tr><td>Superfície Polar Topológica (TPSA)</td><td class="text-mono">{tpsa:.2f} &Aring;&sup2;</td><td>&le; 140.00 &Aring;&sup2;</td><td>{tpsa_status}</td></tr>
        <tr><td>Ligações Rotacionáveis (RotB)</td><td class="text-mono">{rotb}</td><td>&le; 10</td><td>{rotb_status}</td></tr>

        <tr class="section-header"><td colspan="4">2. Farmacocinética e Biodisponibilidade (ADME)</td></tr>
        <tr><td>Absorção Intestinal Humana (HIA)</td><td class="text-mono">{hia_status}</td><td>Egan Egg (TPSA &le; 132 & -1.0 &le; LogP &le; 5.8)</td><td>{hia_badge}</td></tr>
        <tr><td>Permeabilidade Hematoencefálica (BBB)</td><td class="text-mono">{bbb_status}</td><td>Clark (Neutra, TPSA &lt; 90 & 1.0 &le; LogP &le; 5.0)</td><td>{bbb_badge}</td></tr>
        <tr><td>Substrato de P-glicoproteína (P-gp)</td><td class="text-mono">{pgp_status}</td><td>MW &gt; 400 & TPSA &gt; 80</td><td>{pgp_badge}</td></tr>

        <tr class="section-header"><td colspan="4">3. Toxicologia e Alertas Estruturais (T)</td></tr>
        <tr><td>Subestruturas Reativas / PAINS</td><td class="text-mono">{tox_val}</td><td>Filtro de Subestruturas Tóxicas</td><td>{tox_badge}</td></tr>
        """
    else:
        admet_rows_html = """<tr><td colspan="4" class="text-center text-muted">Nenhum dado ADMET disponível</td></tr>"""

    # Geração das Linhas da Tabela de Interações (PLIP)
    interaction_rows_html = ""
    if hbonds or hcontacts:
        for hb in hbonds:
            res_label = f"<strong>{hb.get('resname', 'UNK')}</strong> {hb.get('resnr', '')}"
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
            res_label = f"<strong>{hc.get('resname', 'UNK')}</strong> {hc.get('resnr', '')}"
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
    rmsd_img = data["plots"]["rmsd"]
    rmsf_img = data["plots"]["rmsf"]
    hbond_img = data["plots"]["hbond"]

    def render_plot_card(title: str, subtitle: str, img_b64: Optional[str], fallback_desc: str) -> str:
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
        "RMSD - Desvio Quadrático Médio",
        "Estabilidade Estrutural do Backbone da Proteína ao Longo do Tempo (ns)",
        rmsd_img,
        "Gráfico rmsd.png não encontrado no diretório de trabalho."
    )
    rmsf_card = render_plot_card(
        "RMSF - Flutuação Atômica por Resíduo",
        "Flexibilidade Conformacional dos Carbonos Alfa (C-α)",
        rmsf_img,
        "Gráfico rmsf.png não encontrado no diretório de trabalho."
    )
    hbond_card = render_plot_card(
        "Pontes de Hidrogênio Intermoleculares",
        "Monitoramento Temporal de Contatos Específicos Receptor-Ligante",
        hbond_img,
        "Gráfico hbond.png não encontrado no diretório de trabalho."
    )

    # Tabela MM-PBSA detalhada
    mmpbsa_table_html = ""
    if mmpbsa_energies:
        vdw = mmpbsa_energies.get("van_der_waals", {"mean": 0.0, "std": 0.0})
        eel = mmpbsa_energies.get("electrostatic", {"mean": 0.0, "std": 0.0})
        polar = mmpbsa_energies.get("polar_solvation", {"mean": 0.0, "std": 0.0})
        apolar = mmpbsa_energies.get("nonpolar_solvation", {"mean": 0.0, "std": 0.0})
        dg = mmpbsa_energies.get("delta_g_binding", {"mean": 0.0, "std": 0.0})
        unit = mmpbsa_data.get("unit", "kcal/mol")

        mmpbsa_table_html = f"""
        <div class="card mt-4">
            <div class="card-header">
                <h3>Decomposição de Energia Livre MM-PBSA (Solvente Explícito)</h3>
                <span class="badge badge-primary">{unit}</span>
            </div>
            <div class="card-body">
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
                                <td>Atrações dispersivas e empacotamento estérico</td>
                                <td class="text-right text-mono">{vdw.get('mean', 0.0):.2f} &plusmn; {vdw.get('std', 0.0):.2f}</td>
                            </tr>
                            <tr>
                                <td><strong>Eletrostática (&Delta;E<sub>elec</sub>)</strong></td>
                                <td>Interações de Coulomb e pares iônicos</td>
                                <td class="text-right text-mono">{eel.get('mean', 0.0):.2f} &plusmn; {eel.get('std', 0.0):.2f}</td>
                            </tr>
                            <tr>
                                <td><strong>Solvatação Polar (&Delta;G<sub>polar</sub>)</strong></td>
                                <td>Custo de dessolvatação eletrostática (Poisson-Boltzmann)</td>
                                <td class="text-right text-mono">{polar.get('mean', 0.0):.2f} &plusmn; {polar.get('std', 0.0):.2f}</td>
                            </tr>
                            <tr>
                                <td><strong>Solvatação Apolar (&Delta;G<sub>apolar</sub>)</strong></td>
                                <td>Efeito hidrofóbico e área de superfície acessível ao solvente (SASA)</td>
                                <td class="text-right text-mono">{apolar.get('mean', 0.0):.2f} &plusmn; {apolar.get('std', 0.0):.2f}</td>
                            </tr>
                            <tr class="highlight-row">
                                <td><strong>&Delta;G Total de Ligação (&Delta;G<sub>bind</sub>)</strong></td>
                                <td><strong>Afinidade Termodinâmica Global MM-PBSA</strong></td>
                                <td class="text-right text-mono font-bold text-success">{dg.get('mean', 0.0):.2f} &plusmn; {dg.get('std', 0.0):.2f}</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
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
                <p>Análise Integrada de Docking Molecular, Interações Atômicas (PLIP), Triagem ADMET e Dinâmica Molecular (GROMACS).</p>
            </div>
            <div class="hero-meta">
                <div><strong>Data:</strong> {data['generated_at']}</div>
                <div><strong>Diretório:</strong> {data['work_dir']}</div>
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
                <div class="metric-subtext">&Delta;G<sub>bind</sub> solvatação explícita</div>
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
                <div class="veredito-banner {'veredito-success' if admet_pass else 'veredito-danger'}">
                    <div class="veredito-icon">{'✅' if admet_pass else '⚠️'}</div>
                    <div class="veredito-content">
                        <h4>{'Composto Aprovado na Triagem ADMET' if admet_pass else 'Atenção: Alertas ou Violações ADMET Identificadas'}</h4>
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

        <!-- Seção 3: Galeria de Dinâmica Molecular (GROMACS) -->
        <section class="card">
            <div class="card-header">
                <h3>📈 Monitoramento de Estabilidade na Dinâmica Molecular (GROMACS)</h3>
                <span class="badge badge-secondary">Publicação 300 DPI</span>
            </div>
            <div class="card-body">
                <div class="gallery-grid">
                    {rmsd_card}
                    {rmsf_card}
                    {hbond_card}
                </div>
            </div>
        </section>

        <!-- Seção 4: Decomposição MM-PBSA (Se disponível) -->
        {mmpbsa_table_html}

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
