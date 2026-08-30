#!/usr/bin/env python

import json
import time
from pathlib import Path

import questionary
import typer
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table

from docking import (
    analysis,
    box_utils,
    md_analysis,
    md_equil,
    md_prep,
    notifier,
    pharmacokinetics,
    preparation,
    report,
    vina_runner,
    visualization,
)
from docking.preparation import get_executable

app = typer.Typer(help="Pipeline de Docking Molecular Automatizado")
console = Console()

# Definição de caminhos base
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
VINA_BIN = None


def get_vina_bin() -> Path:
    import os
    import platform
    import shutil
    import stat
    import urllib.request

    bin_dir = BASE_DIR / "bin"
    bin_dir.mkdir(exist_ok=True)

    if platform.system() == "Windows":
        vina_bin = bin_dir / "vina.exe"
        if not vina_bin.exists():
            console.print("[yellow]Baixando AutoDock Vina para Windows...[/yellow]")
            url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_win64.exe"
            urllib.request.urlretrieve(url, vina_bin)
        return vina_bin
    else:
        vina_bin = bin_dir / "vina"
        if not vina_bin.exists():
            # Tenta verificar se já existe vina globalmente no PATH
            system_vina = shutil.which("vina")
            if system_vina:
                return Path(system_vina)

            console.print("[yellow]Baixando AutoDock Vina para Linux...[/yellow]")
            url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
            try:
                urllib.request.urlretrieve(url, vina_bin)
                st = os.stat(vina_bin)
                os.chmod(vina_bin, st.st_mode | stat.S_IEXEC)
            except Exception as e:
                raise RuntimeError(
                    f"AutoDock Vina para Linux não pôde ser baixado e não está instalado no PATH global. "
                    f"Erro original: {e}"
                )
        return vina_bin


def render_interactions_table(interactions: dict):
    """
    Exibe uma tabela no terminal com os aminoácidos do receptor que fizeram
    contatos estáticos (pontes de hidrogênio e contatos hidrofóbicos) na Pose 1.
    """
    table = Table(
        title="[bold magenta]Interações Estáticas Receptor-Ligante (Pose 1)[/bold magenta]",
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Resíduo", style="yellow")
    table.add_column("Tipo de Interação", style="green")
    table.add_column("Distância (Å)", justify="right", style="white")

    # Adiciona as pontes de hidrogênio
    hbonds = interactions.get("hydrogen_bonds", [])
    for hb in hbonds:
        res = f"{hb['resname']} {hb['resnr']}"
        table.add_row(res, "Ponte de Hidrogênio", f"{hb['distance']:.2f}")

    # Adiciona os contatos hidrofóbicos
    hcontacts = interactions.get("hydrophobic_contacts", [])
    for hc in hcontacts:
        res = f"{hc['resname']} {hc['resnr']}"
        table.add_row(res, "Contato Hidrofóbico", f"{hc['distance']:.2f}")

    # Mensagem caso não existam interações mapeadas
    if not hbonds and not hcontacts:
        table.add_row("Nenhuma interação mapeada", "-", "-")

    console.print(table)


def render_admet_table(admet: dict):
    """
    Exibe uma tabela detalhada e organizada em seções no terminal com os
    descritores físico-químicos, predições farmacocinéticas e alertas de toxicidade.
    """
    if "error" in admet:
        console.print(f"[bold red]Erro ao calcular ADMET:[/bold red] {admet['error']}")
        return

    table = Table(
        title="[bold magenta]Triagem ADMET e Perfil Farmacocinético[/bold magenta]",
        show_header=True,
        header_style="bold cyan",
    )
    table.add_column("Propriedade/Parâmetro", style="yellow")
    table.add_column("Valor/Resultado", justify="right", style="white")
    table.add_column("Critério/Limite", style="blue")
    table.add_column("Status/Veredito", justify="center")

    # --- SEÇÃO 1: Físico-Química (Lipinski & Veber) ---
    table.add_row(
        "[bold cyan]1. Parâmetros Físico-Químicos (Lipinski/Veber)[/bold cyan]",
        "",
        "",
        "",
    )

    mw = admet.get("molecular_weight", 0.0)
    mw_status = (
        "[bold green]OK[/bold green]" if mw <= 500 else "[bold red]VIOLADO[/bold red]"
    )
    table.add_row("  Peso Molecular (MW)", f"{mw:.2f} g/mol", "<= 500.00", mw_status)

    logp = admet.get("logp", 0.0)
    logp_status = (
        "[bold green]OK[/bold green]" if logp <= 5 else "[bold red]VIOLADO[/bold red]"
    )
    table.add_row("  Lipofilicidade (LogP)", f"{logp:.2f}", "<= 5.00", logp_status)

    hbd = admet.get("hydrogen_bond_donors", 0)
    hbd_status = (
        "[bold green]OK[/bold green]" if hbd <= 5 else "[bold red]VIOLADO[/bold red]"
    )
    table.add_row("  Doadores de H (HBD)", str(hbd), "<= 5", hbd_status)

    hba = admet.get("hydrogen_bond_acceptors", 0)
    hba_status = (
        "[bold green]OK[/bold green]" if hba <= 10 else "[bold red]VIOLADO[/bold red]"
    )
    table.add_row("  Aceitadores de H (HBA)", str(hba), "<= 10", hba_status)

    tpsa = admet.get("tpsa", 0.0)
    tpsa_status = (
        "[bold green]OK[/bold green]" if tpsa <= 140 else "[bold red]VIOLADO[/bold red]"
    )
    table.add_row(
        "  Superfície Polar (TPSA)", f"{tpsa:.2f} Å²", "<= 140.00", tpsa_status
    )

    rotb = admet.get("rotatable_bonds", 0)
    rotb_status = (
        "[bold green]OK[/bold green]" if rotb <= 10 else "[bold red]VIOLADO[/bold red]"
    )
    table.add_row("  Ligações Rotacionáveis", str(rotb), "<= 10", rotb_status)

    # --- SEÇÃO 2: Predições Farmacocinéticas (ADME) ---
    table.add_section()
    table.add_row(
        "[bold cyan]2. Predições Farmacocinéticas (ADME)[/bold cyan]", "", "", ""
    )

    hia = admet.get("hia_status", "N/A")
    hia_status = (
        "[bold green]Alta[/bold green]"
        if hia == "Alta Absorção"
        else "[bold red]Baixa[/bold red]"
    )
    table.add_row(
        "  Absorção Intestinal (HIA)",
        hia,
        "Egan Egg (TPSA<=132 & -1.0<=LogP<=5.8)",
        hia_status,
    )

    bbb = admet.get("bbb_status", "N/A")
    bbb_status = (
        "[bold green]Permeável[/bold green]"
        if bbb == "Permeável"
        else "[bold yellow]Incompatível/Baixa[/bold yellow]"
    )
    table.add_row(
        "  Permeabilidade SNC (BBB)",
        bbb,
        "Clark (Neutra, TPSA<90 & 1.0<=LogP<=5.0)",
        bbb_status,
    )

    pgp = admet.get("pgp_status", "N/A")
    is_pgp_substrate = (
        admet.get("pgp_substrate")
        if "pgp_substrate" in admet
        else (pgp.startswith("Substrato") or ("Substrato" in pgp and "Não" not in pgp))
    )
    pgp_status = (
        "[bold yellow]Efluxo Ativo[/bold yellow]"
        if is_pgp_substrate
        else "[bold green]Baixo Efluxo[/bold green]"
    )
    table.add_row("  Perfil de Efluxo (P-gp)", pgp, "MW > 400 & TPSA > 80", pgp_status)

    # --- SEÇÃO 3: Triagem de Toxicidade (T) ---
    table.add_section()
    table.add_row("[bold cyan]3. Triagem de Toxicidade (T)[/bold cyan]", "", "", "")

    toxic_alerts = admet.get("toxic_alerts", [])
    if toxic_alerts:
        tox_status = "[bold red]ALERTA[/bold red]"
        tox_val = ", ".join(toxic_alerts)
    else:
        tox_status = "[bold green]Seguro[/bold green]"
        tox_val = "Nenhum alerta estrutural encontrado"
    table.add_row(
        "  Alertas Estruturais (PAINS)",
        tox_val,
        "Subestruturas Reativas / PAINS",
        tox_status,
    )

    console.print(table)

    # Veredito Geral Integrado (3 Tiers: Aprovado, Aprovado com Ressalvas, Reprovado/Risco)
    verdict_cat = admet.get("verdict_category")
    total_viol = admet.get("total_violations", 0)
    all_viol = admet.get("all_violations", [])
    hia_status_val = admet.get("hia_status", "")
    dynamic_points = admet.get("dynamic_points", [])
    attention_note = admet.get("attention_note", "")

    # Fallback caso os campos estendidos não estejam presentes
    if not verdict_cat:
        has_severe_risk = (hia_status_val == "Baixa Absorção") or (len(toxic_alerts) > 0)
        if total_viol == 0 and not has_severe_risk:
            verdict_cat = "APPROVED"
        elif total_viol == 1 and not has_severe_risk:
            verdict_cat = "MODERATE"
        else:
            verdict_cat = "RISK"

    if verdict_cat == "APPROVED":
        veredito = "[bold white on green]  APROVADO (BIODISPONÍVEL & SEGURO)  [/bold white on green]"
        points_str = "\n".join(dynamic_points) if dynamic_points else (
            "• Físico-Química: 100% de conformidade com as regras clássicas de Lipinski e Veber (0 violações).\n"
            "• Farmacocinética: Alta Absorção Intestinal (HIA) estimada (Egan Egg).\n"
            "• Toxicidade: Nenhum alerta estrutural reativo ou PAINS identificado."
        )
        message = (
            f"Veredito de Triagem ADMET: {veredito}\n"
            f"{points_str}\n"
            f"[green]Nota:[/green] {attention_note or 'A molécula possui excelente perfil biofarmacêutico para desenvolvimento oral.'}"
        )
        border_style = "green"

    elif verdict_cat == "MODERATE":
        veredito = "[bold black on yellow]  APROVADO COM RESSALVAS (ALERTA MODERADO)  [/bold black on yellow]"
        points_str = "\n".join(dynamic_points) if dynamic_points else (
            f"• Desvio Pontual Tolerado: {all_viol[0] if all_viol else '1 desvio físico-químico'} (desvio único aceito em fármacos comerciais).\n"
            "• Farmacocinética: Mantém Alta Absorção Intestinal (HIA) estimada (Egan Egg).\n"
            "• Toxicidade: Nenhum alerta estrutural reativo ou PAINS identificado."
        )
        message = (
            f"Veredito de Triagem ADMET: {veredito}\n"
            f"{points_str}\n"
            f"[yellow]Atenção:[/yellow] {attention_note or 'A molécula mantém bom perfil de absorção e ausência de toxicidade, recomendando-se atenção ao desvio pontual.'}"
        )
        border_style = "yellow"

    else:
        veredito = "[bold white on red]  REPROVADO / RISCO ADMET  [/bold white on red]"
        points_str = "\n".join(dynamic_points) if dynamic_points else (
            f"• Problemas Identificados: {', '.join(all_viol) if all_viol else 'Critérios ADMET não atendidos.'}"
        )
        message = (
            f"Veredito de Triagem ADMET: {veredito}\n"
            f"{points_str}\n"
            f"[red]Atenção:[/red] {attention_note or 'A molécula possui restrições que podem comprometer sua biodisponibilidade ou segurança.'}"
        )
        border_style = "red"

    console.print(Panel(message, border_style=border_style))


@app.command(name="validate")
def validate(
    pdb_id: str = typer.Option(
        "4HG7", "--pdb", help="ID do PDB para baixar e analisar"
    ),
    exhaustiveness: int = typer.Option(
        16, "--ex", help="Exaustividade do Vina (padrão: 16)"
    ),
):
    """
    POSITIVE CONTROL / REDOCKING:
    Baixa um PDB, separa o ligante nativo, prepara e executa o docking para validar o RMSD.
    """
    global VINA_BIN
    if VINA_BIN is None:
        VINA_BIN = get_vina_bin()

    run_dir = DATA_DIR / pdb_id
    raw_dir = run_dir / "raw"
    processed_dir = run_dir / "processed"
    results_dir = run_dir / "results"

    for folder in [raw_dir, processed_dir, results_dir]:
        folder.mkdir(parents=True, exist_ok=True)

    console.print(
        Panel.fit(
            f"[bold blue]Pipeline de Validação (Redocking)[/bold blue]\n"
            f"PDB ID: {pdb_id} | Exhaustiveness: {exhaustiveness}\n"
            f"Output: {run_dir}",
            border_style="blue",
        )
    )

    start_time = time.time()
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task1 = progress.add_task(
                description="Baixando e separando PDB...", total=1
            )
            pdb_file = preparation.download_pdb(pdb_id, raw_dir)
            rec_pdb, lig_pdb = preparation.split_receptor_ligand(
                pdb_file, processed_dir
            )
            progress.update(task1, completed=1)

            task2 = progress.add_task(
                description="Preparando receptor e ligante (PDBQT)...", total=1
            )
            rec_pdbqt = processed_dir / "receptor.pdbqt"
            lig_pdbqt = processed_dir / "ligand.pdbqt"
            preparation.prepare_receptor(rec_pdb, rec_pdbqt)
            preparation.prepare_ligand(lig_pdb, lig_pdbqt)
            progress.update(task2, completed=1)

            task3 = progress.add_task(description="Calculando grid box...", total=1)
            box_params = box_utils.calculate_centroid(lig_pdb)
            progress.update(task3, completed=1)

            task4 = progress.add_task(description="Rodando AutoDock Vina...", total=1)
            docked_out = results_dir / "docked.pdbqt"
            vina_log = results_dir / "vina_log.txt"
            vina_runner.run_vina(
                VINA_BIN,
                rec_pdbqt,
                lig_pdbqt,
                box_params,
                docked_out,
                vina_log,
                exhaustiveness,
            )
            progress.update(task4, completed=1)

            task5 = progress.add_task(description="Analisando resultados...", total=1)
            score = analysis.extract_vina_score(vina_log)
            rmsd, error = analysis.analyze_results(docked_out, lig_pdb, results_dir)
            progress.update(task5, completed=1)

            # Executa fluxo do PLIP imediatamente após a análise inicial
            task_plip = progress.add_task(
                description="Executando PLIP (Docker)...", total=1
            )
            complex_pdb = results_dir / "complex.pdb"
            analysis.generate_complex_pdb(
                rec_pdb, results_dir / "docked_poses.sdf", complex_pdb
            )

            plip_ok, plip_msg = analysis.run_plip_docker(complex_pdb, results_dir)
            if not plip_ok:
                raise RuntimeError(plip_msg)

            interactions = analysis.parse_plip_xml(results_dir / "complex_report.xml")

            # Triagem ADMET
            try:
                admet = pharmacokinetics.calculate_admet_descriptors(
                    results_dir / "docked_poses.sdf"
                )
            except Exception as admet_err:
                admet = {"error": str(admet_err), "pass_filters": False}

            interactions["pharmacokinetics"] = admet  # type: ignore

            # Salva o arquivo JSON consolidado
            with open(results_dir / "interactions.json", "w") as f:
                json.dump(interactions, f, indent=4)

            progress.update(task_plip, completed=1)

        duration = time.time() - start_time
        console.print("\n[bold green]✓ Validação concluída![/bold green]")
        console.print(f"[bold]Energia de Afinidade (Score):[/bold] {score} kcal/mol")

        veredito_str = "Indefinido"
        if rmsd is not None:
            color = "green" if rmsd <= 2.0 else "red"
            veredito_str = "SUCESSO" if rmsd <= 2.0 else "FRACASSO"
            console.print(
                f"[bold]RMSD (vs Cristal):[/bold] [{color}]{rmsd:.3f} Å[/{color}]"
            )
            console.print(f"[bold]Veredito:[/bold] [{color}]{veredito_str}[/{color}]")
        else:
            console.print(f"[bold red]Erro na análise:[/bold red] {error}")

        # Renderiza a tabela de contatos estáticos no terminal
        render_interactions_table(interactions)
        render_admet_table(interactions.get("pharmacokinetics", {}))  # type: ignore

        # Notificação por E-mail
        details = {
            "PDB ID": pdb_id,
            "Afinidade Vina (Score)": f"{score} kcal/mol",
            "Exaustividade": str(exhaustiveness),
            "RMSD vs Cristal": f"{rmsd:.3f} Å" if rmsd is not None else "N/A",
            "Veredito": veredito_str,
            "Diretório de Resultados": str(results_dir),
        }
        admet_info = interactions.get("pharmacokinetics", {})
        if isinstance(admet_info, dict) and "pass_filters" in admet_info:
            v_cat = admet_info.get("verdict_category", "")
            if v_cat == "APPROVED":
                v_label = "Aprovado (Biodisponível & Seguro)"
            elif v_cat == "MODERATE":
                v_label = "Aprovado com Ressalvas (Alerta Moderado)"
            elif v_cat == "RISK":
                v_label = "Reprovado / Risco ADMET"
            else:
                v_label = "Aprovado" if admet_info.get("pass_filters") else "Reprovado / Risco"
            details["Triagem ADMET"] = v_label
        notifier.send_email_alert(
            step_name=f"Validação / Redocking ({pdb_id})",
            status="success"
            if (rmsd is not None and rmsd <= 2.0)
            else ("warning" if rmsd else "error"),
            duration_seconds=duration,
            details=details,
            console_logger=console,
        )

    except Exception as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name=f"Validação / Redocking ({pdb_id})",
            status="error",
            duration_seconds=duration,
            error_message=str(e),
            console_logger=console,
        )
        console.print(f"\n[bold red]FATAL ERROR:[/bold red] {e}")
        raise typer.Exit(code=1)


@app.command(name="screen")
def screen(
    receptor: Path = typer.Option(
        ..., "--receptor", help="Caminho para o receptor preparado (.pdbqt)"
    ),
    ligand: Path = typer.Option(
        ..., "--ligand", help="Caminho para o novo composto preparado (.pdbqt)"
    ),
    cx: float = typer.Option(..., "--cx", help="Coordenada X do centro do sítio ativo"),
    cy: float = typer.Option(..., "--cy", help="Coordenada Y do centro do sítio ativo"),
    cz: float = typer.Option(..., "--cz", help="Coordenada Z do centro do sítio ativo"),
    size: float = typer.Option(22.0, "--size", help="Tamanho da caixa (A)"),
    exhaustiveness: int = typer.Option(16, "--ex", help="Exaustividade do Vina"),
):
    """
    TRIAGEM VIRTUAL (VIRTUAL SCREENING):
    Executa o docking de um novo xenobiótico em um receptor já preparado em coordenadas específicas.
    """
    global VINA_BIN
    if VINA_BIN is None:
        VINA_BIN = get_vina_bin()

    # Isolamento de output pelo nome do ligante
    ligand_name = ligand.stem
    results_dir = DATA_DIR / "screening" / ligand_name
    results_dir.mkdir(parents=True, exist_ok=True)

    console.print(
        Panel.fit(
            f"[bold cyan]Triagem Virtual[/bold cyan]\n"
            f"Receptor: {receptor.name} | Ligante: {ligand_name}\n"
            f"Box: Center({cx}, {cy}, {cz}) | Size({size})",
            border_style="cyan",
        )
    )

    # Montagem dos parâmetros da caixa
    box_params = {
        "center_x": cx,
        "center_y": cy,
        "center_z": cz,
        "size_x": size,
        "size_y": size,
        "size_z": size,
    }

    start_time = time.time()
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            # Execução do Vina
            task1 = progress.add_task(
                description=f"Rodando Vina para {ligand_name}...", total=1
            )
            docked_out = results_dir / f"{ligand_name}_docked.pdbqt"
            vina_log = results_dir / f"{ligand_name}_vina.log"

            vina_runner.run_vina(
                VINA_BIN,
                receptor,
                ligand,
                box_params,
                docked_out,
                vina_log,
                exhaustiveness,
            )
            progress.update(task1, completed=1)

            # Extração de score
            task2 = progress.add_task(description="Extraindo score...", total=1)
            score = analysis.extract_vina_score(vina_log)
            progress.update(task2, completed=1)

            # Exporta para SDF e executa o fluxo do PLIP
            task_plip = progress.add_task(
                description="Executando PLIP (Docker)...", total=1
            )
            sdf_out = results_dir / "docked_poses.sdf"
            import subprocess

            exec_name = get_executable("mk_export")
            subprocess.run([exec_name, str(docked_out), "-s", str(sdf_out)], check=True)

            # Resolve o receptor PDB correspondente
            receptor_pdb = receptor.with_suffix(".pdb")
            if not receptor_pdb.exists():
                pdbs = list(receptor.parent.glob("*.pdb"))
                if pdbs:
                    receptor_pdb = pdbs[0]
                else:
                    raise FileNotFoundError(
                        f"Não foi possível encontrar o arquivo receptor PDB em {receptor.parent}"
                    )

            complex_pdb = results_dir / "complex.pdb"
            analysis.generate_complex_pdb(receptor_pdb, sdf_out, complex_pdb)

            plip_ok, plip_msg = analysis.run_plip_docker(complex_pdb, results_dir)
            if not plip_ok:
                raise RuntimeError(plip_msg)

            interactions = analysis.parse_plip_xml(results_dir / "complex_report.xml")

            # Triagem ADMET
            try:
                admet = pharmacokinetics.calculate_admet_descriptors(sdf_out)
            except Exception as admet_err:
                admet = {"error": str(admet_err), "pass_filters": False}

            interactions["pharmacokinetics"] = admet  # type: ignore

            # Salva o arquivo JSON consolidado
            with open(results_dir / "interactions.json", "w") as f:
                json.dump(interactions, f, indent=4)

            progress.update(task_plip, completed=1)

        duration = time.time() - start_time
        console.print(
            f"\n[bold green]✓ Triagem concluída para {ligand_name}![/bold green]"
        )
        console.print(
            f"[bold]Energia de Afinidade (Score):[/bold] [yellow]{score}[/yellow] kcal/mol"
        )
        console.print(f"[bold]Resultado salvo em:[/bold] {docked_out}")

        # Renderiza a tabela de contatos estáticos no terminal
        render_interactions_table(interactions)
        render_admet_table(interactions.get("pharmacokinetics", {}))  # type: ignore

        # Notificação por E-mail
        details = {
            "Receptor": receptor.name,
            "Ligante": ligand_name,
            "Afinidade Vina (Score)": f"{score} kcal/mol",
            "Exaustividade": str(exhaustiveness),
            "Diretório de Resultados": str(results_dir),
        }
        admet_info = interactions.get("pharmacokinetics", {})
        if isinstance(admet_info, dict) and "pass_filters" in admet_info:
            v_cat = admet_info.get("verdict_category", "")
            if v_cat == "APPROVED":
                v_label = "Aprovado (Biodisponível & Seguro)"
            elif v_cat == "MODERATE":
                v_label = "Aprovado com Ressalvas (Alerta Moderado)"
            elif v_cat == "RISK":
                v_label = "Reprovado / Risco ADMET"
            else:
                v_label = "Aprovado" if admet_info.get("pass_filters") else "Reprovado / Risco"
            details["Veredito ADMET"] = v_label

        notifier.send_email_alert(
            step_name=f"Triagem Virtual ({ligand_name})",
            status="success",
            duration_seconds=duration,
            details=details,
            console_logger=console,
        )

    except Exception as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name=f"Triagem Virtual ({ligand_name})",
            status="error",
            duration_seconds=duration,
            error_message=str(e),
            console_logger=console,
        )
        console.print(f"\n[bold red]FATAL ERROR during screening:[/bold red] {e}")
        raise typer.Exit(code=1)


@app.command(name="interactive")
def interactive():
    """Interface interativa (TUI) para facilitar o uso do pipeline."""
    while True:
        choice = questionary.select(
            "O que você deseja fazer?",
            choices=[
                "1. Validação (Redocking)",
                "2. Download de Ligante (PubChem)",
                "3. Preparação de Ligante (SDF -> PDBQT)",
                "4. Triagem Virtual (Screening)",
                "5. Preparar Dinâmica Molecular (GROMACS)",
                "6. Rodar Equilíbrio da Dinâmica (NVT/NPT)",
                "7. Compilar TPR de Produção (grompp -> md.tpr)",
                "8. Executar Produção da Dinâmica (100 ns)",
                "9. Pós-processamento, Gráficos e MM-PBSA da DM",
                "10. Gerar Relatório Executivo (HTML) e Script PyMOL (3D)",
                "11. Testar Configuração de E-mail de Alerta",
                "12. Sair",
            ],
        ).ask()

        if choice == "1. Validação (Redocking)":
            pdb_id = questionary.text(
                "Digite o ID do PDB (ex: 4HG7):", default="4HG7"
            ).ask()
            ex = questionary.text("Exaustividade (ex: 16):", default="16").ask()
            validate(pdb_id=pdb_id, exhaustiveness=int(ex))

        elif choice == "2. Download de Ligante (PubChem)":
            cid = questionary.text("Digite o CID do composto no PubChem:").ask()
            name = questionary.text(
                "Digite o nome do arquivo (ex: desoxicolato.sdf):"
            ).ask()
            if not name.endswith(".sdf"):
                name += ".sdf"
            out_path = DATA_DIR / name
            preparation.download_pubchem_sdf(cid, out_path)
            console.print(f"[bold green]✓ Download concluído:[/bold green] {out_path}")

        elif choice == "3. Preparação de Ligante (SDF -> PDBQT)":
            sdf_file = questionary.path("Caminho para o arquivo SDF:").ask()
            name = Path(sdf_file).stem
            out_pdbqt = questionary.text(
                "Caminho de saída (PDBQT):", default=f"data/{name}.pdbqt"
            ).ask()
            preparation.prepare_ligand_sdf(Path(sdf_file), Path(out_pdbqt))
            console.print(
                f"[bold green]✓ Preparação concluída:[/bold green] {out_pdbqt}"
            )

        elif choice == "4. Triagem Virtual (Screening)":
            receptor = questionary.path("Caminho para o receptor (.pdbqt):").ask()
            ligand = questionary.path("Caminho para o ligante (.pdbqt):").ask()
            cx = questionary.text("Coordenada X:").ask()
            cy = questionary.text("Coordenada Y:").ask()
            cz = questionary.text("Coordenada Z:").ask()
            size = questionary.text("Tamanho da caixa (A):", default="22.0").ask()
            ex = questionary.text("Exaustividade (ex: 16):", default="16").ask()

            screen(
                receptor=Path(receptor),
                ligand=Path(ligand),
                cx=float(cx),
                cy=float(cy),
                cz=float(cz),
                size=float(size),
                exhaustiveness=int(ex),
            )

        elif choice == "5. Preparar Dinâmica Molecular (GROMACS)":
            rec_default = ""
            for candidate in [
                Path("data/7CFN/processed/receptor.pdb"),
                Path("data/1OSV/processed/receptor.pdb"),
            ]:
                if candidate.exists():
                    rec_default = str(candidate)
                    break
            if not rec_default:
                found_rec = list(Path("data").glob("*/processed/receptor.pdb"))
                if found_rec:
                    rec_default = str(found_rec[0])

            sdf_default = ""
            for candidate_sdf in [
                Path("data/7CFN/results/docked_poses.sdf"),
                Path("data/screening/desoxicolato/docked_poses.sdf"),
            ]:
                if candidate_sdf.exists():
                    sdf_default = str(candidate_sdf)
                    break
            if not sdf_default:
                found_sdf = list(Path("data").glob("**/docked_poses.sdf"))
                if found_sdf:
                    sdf_default = str(found_sdf[0])

            out_default = "data/md_files"

            receptor_path = questionary.path(
                "Caminho para o PDB original da proteína (Receptor):",
                default=rec_default,
            ).ask()

            sdf_path = questionary.path(
                "Caminho para o arquivo docked_poses.sdf (Ligante):",
                default=sdf_default,
            ).ask()

            output_dir = questionary.text(
                "Diretório de saída para a Dinâmica Molecular:", default=out_default
            ).ask()

            if not receptor_path or not sdf_path or not output_dir:
                console.print(
                    "[bold red]Operação cancelada: todos os caminhos devem ser preenchidos.[/bold red]"
                )
                continue

            start_time = time.time()
            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Preparação e Minimização de Energia de Dinâmica Molecular (GROMACS)[/bold blue]\n"
                        f"Receptor: {receptor_path}\n"
                        f"Ligante (SDF): {sdf_path}\n"
                        f"Diretório de Saída: {output_dir}",
                        border_style="blue",
                    )
                )

                with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    console=console,
                ) as progress:
                    tasks = {
                        "A": progress.add_task(
                            description="[A] Cura do receptor com PDBFixer",
                            total=1,
                            start=False,
                        ),
                        "B": progress.add_task(
                            description="[B] Extração do Ligante com RDKit",
                            total=1,
                            start=False,
                        ),
                        "C": progress.add_task(
                            description="[C] Parametrização do Ligante (ACPYPE)",
                            total=1,
                            start=False,
                        ),
                        "D": progress.add_task(
                            description="[D] Topologia da Proteína (pdb2gmx)",
                            total=1,
                            start=False,
                        ),
                        "E": progress.add_task(
                            description="[E] Fusão de Coordenadas (complex.gro)",
                            total=1,
                            start=False,
                        ),
                        "F": progress.add_task(
                            description="[F] Fusão de Topologia (Stitching)",
                            total=1,
                            start=False,
                        ),
                        "G": progress.add_task(
                            description="[G] Definição da Caixa de Simulação (editconf)",
                            total=1,
                            start=False,
                        ),
                        "H": progress.add_task(
                            description="[H] Solvatação do Sistema (solvate)",
                            total=1,
                            start=False,
                        ),
                        "I": progress.add_task(
                            description="[I] Compilação de Íons (grompp)",
                            total=1,
                            start=False,
                        ),
                        "J": progress.add_task(
                            description="[J] Neutralização e Concentração (genion)",
                            total=1,
                            start=False,
                        ),
                        "K": progress.add_task(
                            description="[K] Grompp Definitivo (em.tpr)",
                            total=1,
                            start=False,
                        ),
                        "L": progress.add_task(
                            description="[L] Minimização de Energia (mdrun)",
                            total=1,
                            start=False,
                        ),
                    }

                    for step, status in md_prep.prepare_md_system(
                        Path(receptor_path), Path(sdf_path), Path(output_dir)
                    ):
                        task_id = tasks[step]
                        if status == "start":
                            progress.start_task(task_id)
                        elif status == "success":
                            progress.update(task_id, completed=1)

                duration = time.time() - start_time
                console.print(
                    "\n[bold green]✓ Preparação e Minimização de Energia concluídas com sucesso![/bold green]"
                )
                console.print(
                    f"Arquivos e outputs gerados em: [cyan]{output_dir}[/cyan]"
                )

                notifier.send_email_alert(
                    step_name="Opção 5: Preparação de DM e Minimização de Energia",
                    status="success",
                    duration_seconds=duration,
                    details={
                        "Receptor": str(receptor_path),
                        "Ligante": str(sdf_path),
                        "Diretório de Saída": str(output_dir),
                        "Minimização (EM)": "em.gro e topol.top gerados com sucesso",
                        "Próxima Ação": "Rodar Equilíbrio NVT/NPT (Opção 6)",
                    },
                    console_logger=console,
                )

            except md_prep.DependencyError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 5: Preparação de DM e Minimização de Energia",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro de Dependência: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro de Dependência GROMACS/ACPYPE:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha de Dependência",
                    )
                )
            except md_prep.SimulationPrepError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 5: Preparação de DM e Minimização de Energia",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro de Preparação: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro de Preparação na Dinâmica Molecular:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha na Preparação",
                    )
                )
            except FileNotFoundError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 5: Preparação de DM e Minimização de Energia",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Arquivo Não Encontrado: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Arquivo Não Encontrado:[/bold red]\n{e}",
                        border_style="red",
                        title="Erro de Caminho",
                    )
                )
            except Exception as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 5: Preparação de DM e Minimização de Energia",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro Inesperado: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro Inesperado no Fluxo de Dinâmica:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha Crítica",
                    )
                )

        elif choice == "6. Rodar Equilíbrio da Dinâmica (NVT/NPT)":
            md_dir_default = "data/md_files"
            md_dir = questionary.text(
                "Diretório de trabalho da Dinâmica Molecular (onde contêm em.gro e topol.top):",
                default=md_dir_default,
            ).ask()

            if not md_dir:
                console.print(
                    "[bold red]Operação cancelada: o diretório de trabalho deve ser preenchido.[/bold red]"
                )
                continue

            start_time = time.time()
            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Equilíbrio Termodinâmico da Dinâmica Molecular (GROMACS)[/bold blue]\n"
                        f"Diretório de Trabalho: {md_dir}",
                        border_style="blue",
                    )
                )

                with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    console=console,
                ) as progress:
                    tasks = {
                        "A": progress.add_task(
                            description="[A] Geração de Grupos e Índices (make_ndx)",
                            total=1,
                            start=False,
                        ),
                        "B": progress.add_task(
                            description="[B] Compilação da Caixa NVT (grompp)",
                            total=1,
                            start=False,
                        ),
                        "C": progress.add_task(
                            description="[C] Execução do Equilíbrio NVT (mdrun)",
                            total=1,
                            start=False,
                        ),
                        "D": progress.add_task(
                            description="[D] Compilação da Caixa NPT (grompp)",
                            total=1,
                            start=False,
                        ),
                        "E": progress.add_task(
                            description="[E] Execução do Equilíbrio NPT (mdrun)",
                            total=1,
                            start=False,
                        ),
                    }

                    for step, status in md_equil.run_md_equilibration(Path(md_dir)):
                        task_id = tasks[step]
                        if status == "start":
                            progress.start_task(task_id)
                        elif status == "success":
                            progress.update(task_id, completed=1)

                duration = time.time() - start_time
                console.print(
                    "\n[bold green]✓ Equilíbrio NVT/NPT concluído com sucesso![/bold green]"
                )
                console.print(
                    f"Estruturas equilibradas geradas em: [cyan]{md_dir}[/cyan]"
                )

                notifier.send_email_alert(
                    step_name="Opção 6: Equilíbrio da Dinâmica (NVT/NPT)",
                    status="success",
                    duration_seconds=duration,
                    details={
                        "Diretório de Trabalho": str(md_dir),
                        "Equilíbrio NVT": "Concluído (nvt.gro gerado)",
                        "Equilíbrio NPT": "Concluído (npt.gro gerado)",
                        "Próxima Ação": "Compilar TPR ou Executar Produção (Opção 7/8)",
                    },
                    console_logger=console,
                )

            except md_prep.DependencyError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 6: Equilíbrio da Dinâmica (NVT/NPT)",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro de Dependência: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro de Dependência GROMACS:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha de Dependência",
                    )
                )
            except md_prep.SimulationPrepError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 6: Equilíbrio da Dinâmica (NVT/NPT)",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro de Execução no GROMACS: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro de Execução no Equilíbrio GROMACS:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha no Equilíbrio",
                    )
                )
            except FileNotFoundError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 6: Equilíbrio da Dinâmica (NVT/NPT)",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Arquivo Não Encontrado: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Arquivo/Diretório Não Encontrado:[/bold red]\n{e}",
                        border_style="red",
                        title="Erro de Caminho",
                    )
                )
            except Exception as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 6: Equilíbrio da Dinâmica (NVT/NPT)",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro Inesperado: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro Inesperado no Fluxo de Equilíbrio:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha Crítica",
                    )
                )

        elif choice == "7. Compilar TPR de Produção (grompp -> md.tpr)":
            md_dir_default = "data/md_files"
            md_dir = questionary.text(
                "Diretório de trabalho da Dinâmica Molecular (onde contêm npt.gro e topol.top):",
                default=md_dir_default,
            ).ask()

            if not md_dir:
                console.print(
                    "[bold red]Operação cancelada: o diretório de trabalho deve ser preenchido.[/bold red]"
                )
                continue

            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Compilação do Arquivo de Produção (md.tpr)[/bold blue]\n"
                        f"Diretório de Trabalho: {md_dir}",
                        border_style="blue",
                    )
                )

                console.print(
                    "[yellow]Compilando md.tpr via GROMACS (grompp)...[/yellow]"
                )
                tpr_path = md_analysis.compile_production_tpr(Path(md_dir))
                console.print(
                    f"[bold green]✓ Arquivo 'md.tpr' gerado com sucesso em:[/bold green] [cyan]{tpr_path}[/cyan]"
                )
                console.print(
                    "[dim]Agora você pode rodar a produção (opção 8) ou transferir o md.tpr para um servidor/cluster com GPU.[/dim]"
                )

            except md_prep.DependencyError as e:
                console.print(
                    Panel(
                        f"[bold red]Erro de Dependência GROMACS:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha de Dependência",
                    )
                )
            except md_prep.SimulationPrepError as e:
                console.print(
                    Panel(
                        f"[bold red]Erro na Compilação do md.tpr:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha no grompp",
                    )
                )
            except FileNotFoundError as e:
                console.print(
                    Panel(
                        f"[bold red]Arquivo/Diretório Não Encontrado:[/bold red]\n{e}",
                        border_style="red",
                        title="Erro de Caminho",
                    )
                )
            except Exception as e:
                console.print(
                    Panel(
                        f"[bold red]Erro Inesperado na Compilação:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha Crítica",
                    )
                )

        elif choice == "8. Executar Produção da Dinâmica (100 ns)":
            md_dir_default = "data/md_files"
            md_dir = questionary.text(
                "Diretório de trabalho da Dinâmica Molecular (onde contêm npt.gro e topol.top):",
                default=md_dir_default,
            ).ask()

            if not md_dir:
                console.print(
                    "[bold red]Operação cancelada: o diretório de trabalho deve ser preenchido.[/bold red]"
                )
                continue

            start_time = time.time()
            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Produção de Dinâmica Molecular (GROMACS)[/bold blue]\n"
                        f"Diretório de Trabalho: {md_dir}",
                        border_style="blue",
                    )
                )

                console.print(
                    "[yellow]Iniciando a produção da Dinâmica Molecular (grompp & mdrun)...[/yellow]"
                )
                md_analysis.run_production_md(Path(md_dir))
                duration = time.time() - start_time
                console.print(
                    "[bold green]✓ Produção concluída com sucesso![/bold green]"
                )
                console.print(
                    "[cyan]Execute a opção 9 para tratamento de PBC, gráficos e MM-PBSA.[/cyan]"
                )

                notifier.send_email_alert(
                    step_name="Opção 8: Produção de Dinâmica Molecular (100 ns)",
                    status="success",
                    duration_seconds=duration,
                    details={
                        "Diretório de Trabalho": str(md_dir),
                        "Trajetória Gerada": "md.xtc e md.gro gerados",
                        "Próxima Ação Recomendada": "Executar Pós-processamento e MM-PBSA (Opção 9)",
                    },
                    console_logger=console,
                )

            except md_prep.DependencyError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 8: Produção de Dinâmica Molecular (100 ns)",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro de Dependência: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro de Dependência GROMACS:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha de Dependência",
                    )
                )
            except md_prep.SimulationPrepError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 8: Produção de Dinâmica Molecular (100 ns)",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro na Simulação GROMACS: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro de Execução no GROMACS:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha na Simulação",
                    )
                )
            except FileNotFoundError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 8: Produção de Dinâmica Molecular (100 ns)",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Arquivo Não Encontrado: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Arquivo/Diretório Não Encontrado:[/bold red]\n{e}",
                        border_style="red",
                        title="Erro de Caminho",
                    )
                )
            except Exception as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 8: Produção de Dinâmica Molecular (100 ns)",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro Inesperado: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro Inesperado no Fluxo de Produção:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha Crítica",
                    )
                )

        elif choice == "9. Pós-processamento, Gráficos e MM-PBSA da DM":
            md_dir_default = "data/md_files"
            md_dir = questionary.text(
                "Diretório de trabalho da Dinâmica Molecular (onde contêm md.tpr e md.xtc):",
                default=md_dir_default,
            ).ask()

            if not md_dir:
                console.print(
                    "[bold red]Operação cancelada: o diretório de trabalho deve ser preenchido.[/bold red]"
                )
                continue

            run_mmpbsa = questionary.confirm(
                "Deseja executar o cálculo de Energia Livre de Ligação MM-PBSA?",
                default=True,
            ).ask()

            start_time = time.time()
            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Pós-processamento, Gráficos e MM-PBSA da Dinâmica Molecular[/bold blue]\n"
                        f"Diretório de Trabalho: {md_dir}",
                        border_style="blue",
                    )
                )

                with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    console=console,
                ) as progress:
                    task_pbc = progress.add_task(
                        description="[A] Tratamento de PBC e Fit rot+trans (md_center.xtc, md_fit.xtc, md_clean.gro)",
                        total=1,
                    )
                    md_analysis.fix_pbc(Path(md_dir))
                    progress.update(task_pbc, completed=1)

                    task_traj = progress.add_task(
                        description="[B] Análise da Trajetória Ajustada (RMSD Backbone, RMSF C-α, HBond)",
                        total=1,
                    )
                    md_analysis.analyze_trajectory(Path(md_dir))
                    progress.update(task_traj, completed=1)

                    task_plot = progress.add_task(
                        description="[C] Geração de Gráficos Científicos de Publicação (300 DPI)",
                        total=1,
                    )
                    plots = md_analysis.plot_md_results(Path(md_dir))
                    progress.update(task_plot, completed=1)

                    mmpbsa_res = None
                    if run_mmpbsa:
                        task_mmpbsa = progress.add_task(
                            description="[D] Cálculo de Energia Livre MM-PBSA (gmx_MMPBSA)",
                            total=1,
                        )
                        mmpbsa_res = md_analysis.calculate_mmpbsa(Path(md_dir))
                        progress.update(task_mmpbsa, completed=1)

                duration = time.time() - start_time
                console.print(
                    "\n[bold green]✓ Pós-processamento concluído com sucesso![/bold green]"
                )

                if plots:
                    console.print(
                        "\n[bold cyan]Gráficos de Publicação Gerados:[/bold cyan]"
                    )
                    for p_name, p_path in plots.items():
                        console.print(
                            f"  • [bold]{p_name.upper()}:[/bold] [green]{p_path}[/green]"
                        )

                details = {
                    "Diretório": str(md_dir),
                    "Tratamento PBC": "Concluído (md_fit.xtc e md_clean.gro gerados)",
                    "Análise de Trajetória": "RMSD, RMSF e HBond gerados",
                    "Gráficos Gerados": f"{len(plots)} gráficos salvos"
                    if plots
                    else "Nenhum",
                }

                if mmpbsa_res and "energies" in mmpbsa_res:
                    table = Table(
                        title="Resumo Termodinâmico de Energia Livre MM-PBSA",
                        show_header=True,
                        header_style="bold magenta",
                    )
                    table.add_column("Componente Energético", style="dim", width=32)
                    table.add_column(
                        f"Energia ({mmpbsa_res.get('unit', 'kcal/mol')})",
                        justify="right",
                    )

                    energies = mmpbsa_res["energies"]
                    table.add_row(
                        "Van der Waals (ΔE_vdw)",
                        f"{energies['van_der_waals']['mean']:.2f} ± {energies['van_der_waals']['std']:.2f}",
                    )
                    table.add_row(
                        "Eletrostática (ΔE_elec)",
                        f"{energies['electrostatic']['mean']:.2f} ± {energies['electrostatic']['std']:.2f}",
                    )
                    table.add_row(
                        "Solvatação Polar (ΔG_polar)",
                        f"{energies['polar_solvation']['mean']:.2f} ± {energies['polar_solvation']['std']:.2f}",
                    )
                    table.add_row(
                        "Solvatação Apolar (ΔG_apolar)",
                        f"{energies['nonpolar_solvation']['mean']:.2f} ± {energies['nonpolar_solvation']['std']:.2f}",
                    )
                    table.add_section()
                    table.add_row(
                        "[bold]ΔG Total de Ligação (ΔG_bind)[/bold]",
                        f"[bold green]{energies['delta_g_binding']['mean']:.2f} ± {energies['delta_g_binding']['std']:.2f}[/bold green]",
                    )

                    console.print(table)
                    console.print(
                        f"Sumário JSON salvo em: [cyan]{Path(md_dir) / 'mmpbsa_summary.json'}[/cyan]"
                    )

                    dg_bind = f"{energies['delta_g_binding']['mean']:.2f} ± {energies['delta_g_binding']['std']:.2f} {mmpbsa_res.get('unit', 'kcal/mol')}"
                    details["MM-PBSA ΔG Total (Ligação)"] = dg_bind
                    details["Sumário JSON"] = str(Path(md_dir) / "mmpbsa_summary.json")

                # Envio do E-mail de Alerta
                notifier.send_email_alert(
                    step_name="Opção 9: Pós-processamento e MM-PBSA da DM",
                    status="success",
                    duration_seconds=duration,
                    details=details,
                    console_logger=console,
                )

            except md_prep.DependencyError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 9: Pós-processamento e MM-PBSA da DM",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro de Dependência: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro de Dependência GROMACS/gmx_MMPBSA:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha de Dependência",
                    )
                )
            except md_prep.SimulationPrepError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 9: Pós-processamento e MM-PBSA da DM",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro no Pós-processamento: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro no Pós-processamento / MM-PBSA:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha na Análise",
                    )
                )
            except FileNotFoundError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 9: Pós-processamento e MM-PBSA da DM",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Arquivo Não Encontrado: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Arquivo/Diretório Não Encontrado:[/bold red]\n{e}",
                        border_style="red",
                        title="Erro de Caminho",
                    )
                )
            except Exception as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name="Opção 9: Pós-processamento e MM-PBSA da DM",
                    status="error",
                    duration_seconds=duration,
                    error_message=f"Erro Inesperado: {e}",
                    console_logger=console,
                )
                console.print(
                    Panel(
                        f"[bold red]Erro Inesperado no Pós-processamento:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha Crítica",
                    )
                )

        elif choice == "10. Gerar Relatório Executivo (HTML) e Script PyMOL (3D)":
            # Sugere um diretório inteligente padrão se disponível
            default_dir = "data/screening/desoxicolato"
            if not Path(default_dir).exists():
                default_dir = "data/1OSV/results"
                if not Path(default_dir).exists():
                    default_dir = "data/md_files"
                    if not Path(default_dir).exists():
                        default_dir = "data"

            work_dir = questionary.text(
                "Diretório contendo os artefatos do pipeline (docking, PLIP, ADMET, DM):",
                default=default_dir,
            ).ask()

            if not work_dir:
                console.print(
                    "[bold red]Operação cancelada: o diretório de trabalho deve ser informado.[/bold red]"
                )
                continue

            receptor_code = questionary.text(
                "Código/Nome do Receptor (ex: 7CFN, GPBAR1):",
                default="7CFN",
            ).ask()

            ligand_name = questionary.text(
                "Nome do Ligante (ex: Desoxicolato, INT-777):",
                default="Desoxicolato",
            ).ask()

            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Geração de Relatório Executivo e Visualização 3D[/bold blue]\n"
                        f"Diretório de Trabalho: {work_dir}\n"
                        f"Receptor: {receptor_code or 'Não informado'} | Ligante: {ligand_name or 'Não informado'}",
                        border_style="blue",
                    )
                )

                with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    console=console,
                ) as progress:
                    task_html = progress.add_task(
                        description="[1/2] Compilando e gerando Relatório Executivo HTML...",
                        total=1,
                    )
                    html_path, missing = report.generate_html_report(
                        Path(work_dir),
                        receptor_name=receptor_code,
                        ligand_name=ligand_name,
                    )
                    progress.update(task_html, completed=1)

                    task_pymol = progress.add_task(
                        description="[2/2] Gerando Script de Visualização PyMOL (.pml)...",
                        total=1,
                    )
                    try:
                        pml_path = visualization.generate_pymol_script(Path(work_dir))
                    except Exception as e:
                        pml_path = None
                        missing.append(f"PyMOL Script: {e}")
                    progress.update(task_pymol, completed=1)

                console.print(
                    "\n[bold green]✓ Entrega de Resultados Concluída com Sucesso![/bold green]"
                )
                console.print(
                    f"  • [bold]Relatório HTML Consolidado:[/bold] [cyan]{html_path}[/cyan]"
                )
                if pml_path:
                    console.print(
                        f"  • [bold]Script PyMOL 3D:[/bold] [cyan]{pml_path}[/cyan]"
                    )
                    console.print(
                        "    [dim]Para abrir a cena 3D no PyMOL, execute: pymol show_complex.pml[/dim]"
                    )

                if missing:
                    console.print(
                        "\n[bold yellow]Avisos / Artefatos Ausentes:[/bold yellow]"
                    )
                    for item in missing:
                        console.print(f"  [yellow]⚠[/yellow] {item}")

            except FileNotFoundError as e:
                console.print(
                    Panel(
                        f"[bold red]Arquivo/Diretório Não Encontrado:[/bold red]\n{e}",
                        border_style="red",
                        title="Erro de Caminho",
                    )
                )
            except Exception as e:
                console.print(
                    Panel(
                        f"[bold red]Erro Inesperado na Geração de Relatórios:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha Crítica",
                    )
                )

        elif choice == "11. Testar Configuração de E-mail de Alerta":
            console.print(
                Panel.fit(
                    "[bold blue]Teste de Configuração de E-mail de Alerta (SMTP_SSL)[/bold blue]",
                    border_style="blue",
                )
            )
            success, msg = notifier.test_email_connection(console_logger=console)
            if success:
                console.print(
                    "\n[bold green]✓ Teste concluído com sucesso! Verifique sua caixa de entrada.[/bold green]"
                )
            else:
                console.print(
                    f"\n[bold red]✗ Falha no envio do e-mail de teste:[/bold red] {msg}"
                )

        elif choice == "12. Sair":
            break


@app.command(name="md-compile")
def md_compile_command(
    working_dir: Path = typer.Option(
        ...,
        "--dir",
        help="Diretório de trabalho contendo os arquivos do equilíbrio (npt.gro, topol.top, etc.)",
    ),
):
    """
    COMPILAÇÃO DO ARQUIVO DE PRODUÇÃO (md.tpr):
    Executa o 'gmx grompp' para compilar o arquivo de produção md.tpr sem iniciar a simulação.
    """
    console.print(
        Panel.fit(
            f"[bold blue]Compilação do Arquivo de Produção (md.tpr)[/bold blue]\n"
            f"Diretório de Trabalho: {working_dir}",
            border_style="blue",
        )
    )

    try:
        console.print("[yellow]Compilando md.tpr via GROMACS (grompp)...[/yellow]")
        tpr_path = md_analysis.compile_production_tpr(working_dir)
        console.print(
            f"[bold green]✓ Arquivo 'md.tpr' gerado com sucesso em:[/bold green] [cyan]{tpr_path}[/cyan]"
        )
        console.print(
            "[dim]Você pode transferir este arquivo para execução em GPU/cluster ou rodar 'md-run'.[/dim]"
        )
    except md_prep.DependencyError as e:
        console.print(f"\n[bold red]Erro de Dependência:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except md_prep.SimulationPrepError as e:
        console.print(f"\n[bold red]Erro na Compilação:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        console.print(f"\n[bold red]Arquivo Não Encontrado:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"\n[bold red]Erro Inesperado:[/bold red]\n{e}")
        raise typer.Exit(code=1)


@app.command(name="md-run")
def md_run_command(
    working_dir: Path = typer.Option(
        ...,
        "--dir",
        help="Diretório de trabalho onde estão os arquivos do equilíbrio (nvt/npt)",
    ),
):
    """
    PRODUÇÃO DE DINÂMICA MOLECULAR:
    Compila e executa a produção da Dinâmica Molecular no GROMACS (grompp e mdrun).
    """
    console.print(
        Panel.fit(
            f"[bold blue]Produção de Dinâmica Molecular (GROMACS)[/bold blue]\n"
            f"Diretório de Trabalho: {working_dir}",
            border_style="blue",
        )
    )

    start_time = time.time()
    try:
        console.print(
            "[yellow]Iniciando a produção da Dinâmica Molecular (grompp & mdrun)...[/yellow]"
        )
        md_analysis.run_production_md(working_dir)
        duration = time.time() - start_time
        console.print("[bold green]✓ Produção concluída com sucesso![/bold green]")
        console.print(
            f"[cyan]Execute 'md-postprocess --dir {working_dir}' para realizar o tratamento de PBC, gráficos e MM-PBSA.[/cyan]"
        )

        notifier.send_email_alert(
            step_name="Produção de Dinâmica Molecular (GROMACS)",
            status="success",
            duration_seconds=duration,
            details={
                "Diretório": str(working_dir),
                "Trajetória Gerada": "md.xtc e md.gro",
                "Próxima Ação Recomendada": f"Executar md-postprocess --dir {working_dir}",
            },
            console_logger=console,
        )

    except md_prep.DependencyError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Produção de Dinâmica Molecular (GROMACS)",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro de Dependência: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro de Dependência:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except md_prep.SimulationPrepError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Produção de Dinâmica Molecular (GROMACS)",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro na Dinâmica: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro na Dinâmica:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Produção de Dinâmica Molecular (GROMACS)",
            status="error",
            duration_seconds=duration,
            error_message=f"Arquivo Não Encontrado: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Arquivo Não Encontrado:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except Exception as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Produção de Dinâmica Molecular (GROMACS)",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro Inesperado: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro Inesperado:[/bold red]\n{e}")
        raise typer.Exit(code=1)


@app.command(name="md-postprocess")
def md_postprocess_command(
    working_dir: Path = typer.Option(
        ...,
        "--dir",
        help="Diretório de trabalho contendo os arquivos de simulação (md.tpr, md.xtc, etc.)",
    ),
    skip_mmpbsa: bool = typer.Option(
        False, "--skip-mmpbsa", help="Pular o cálculo de energia livre MM-PBSA"
    ),
):
    """
    PÓS-PROCESSAMENTO, GRÁFICOS E MM-PBSA DA DINÂMICA MOLECULAR:
    1. Tratamento de Condições Periódicas de Contorno (PBC: remoção de saltos e centralização).
    2. Análise de trajetória (RMSD, RMSF, HBond) utilizando a trajetória corrigida (md_fit.xtc).
    3. Geração automatizada de gráficos científicos (.png a 300 DPI) para publicação.
    4. Cálculo de energia livre de ligação MM-PBSA (gmx_MMPBSA) e exportação de mmpbsa_summary.json.
    """
    working_dir = Path(working_dir)
    console.print(
        Panel.fit(
            f"[bold blue]Pós-processamento, Gráficos e MM-PBSA da Dinâmica Molecular[/bold blue]\n"
            f"Diretório de Trabalho: {working_dir}",
            border_style="blue",
        )
    )

    start_time = time.time()
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task_pbc = progress.add_task(
                description="[1/4] Tratamento de PBC e Fit rot+trans (md_center.xtc, md_fit.xtc, md_clean.gro)...",
                total=1,
            )
            md_analysis.fix_pbc(working_dir)
            progress.update(task_pbc, completed=1)

            task_traj = progress.add_task(
                description="[2/4] Análise da Trajetória Ajustada (RMSD Backbone, RMSF C-α, HBond)...",
                total=1,
            )
            md_analysis.analyze_trajectory(working_dir)
            progress.update(task_traj, completed=1)

            task_plot = progress.add_task(
                description="[3/4] Geração de Gráficos Científicos de Publicação (300 DPI)...",
                total=1,
            )
            plots = md_analysis.plot_md_results(working_dir)
            progress.update(task_plot, completed=1)

            mmpbsa_res = None
            if not skip_mmpbsa:
                task_mmpbsa = progress.add_task(
                    description="[4/4] Cálculo de Energia Livre de Ligação MM-PBSA (gmx_MMPBSA)...",
                    total=1,
                )
                mmpbsa_res = md_analysis.calculate_mmpbsa(working_dir)
                progress.update(task_mmpbsa, completed=1)

        duration = time.time() - start_time
        console.print(
            "\n[bold green]✓ Pós-processamento concluído com sucesso![/bold green]"
        )

        if plots:
            console.print("\n[bold cyan]Gráficos Gerados:[/bold cyan]")
            for p_name, p_path in plots.items():
                console.print(
                    f"  • [bold]{p_name.upper()}:[/bold] [green]{p_path}[/green]"
                )

        details = {
            "Diretório": str(working_dir),
            "Tratamento PBC": "Concluído (md_fit.xtc e md_clean.gro gerados)",
            "Análises de Trajetória": "RMSD, RMSF e HBond gerados",
            "Gráficos Gerados": f"{len(plots)} gráficos salvos" if plots else "Nenhum",
        }

        if mmpbsa_res and "energies" in mmpbsa_res:
            table = Table(
                title="Resultados Termodinâmicos MM-PBSA",
                show_header=True,
                header_style="bold magenta",
            )
            table.add_column("Componente Energético", style="dim", width=32)
            table.add_column(
                f"Energia ({mmpbsa_res.get('unit', 'kcal/mol')})", justify="right"
            )

            energies = mmpbsa_res["energies"]
            table.add_row(
                "Van der Waals (ΔE_vdw)",
                f"{energies['van_der_waals']['mean']:.2f} ± {energies['van_der_waals']['std']:.2f}",
            )
            table.add_row(
                "Eletrostática (ΔE_elec)",
                f"{energies['electrostatic']['mean']:.2f} ± {energies['electrostatic']['std']:.2f}",
            )
            table.add_row(
                "Solvatação Polar (ΔG_polar)",
                f"{energies['polar_solvation']['mean']:.2f} ± {energies['polar_solvation']['std']:.2f}",
            )
            table.add_row(
                "Solvatação Apolar (ΔG_apolar)",
                f"{energies['nonpolar_solvation']['mean']:.2f} ± {energies['nonpolar_solvation']['std']:.2f}",
            )
            table.add_section()
            table.add_row(
                "[bold]ΔG Total de Ligação (ΔG_bind)[/bold]",
                f"[bold green]{energies['delta_g_binding']['mean']:.2f} ± {energies['delta_g_binding']['std']:.2f}[/bold green]",
            )

            console.print(table)
            console.print(
                f"Sumário JSON salvo em: [cyan]{working_dir / 'mmpbsa_summary.json'}[/cyan]"
            )

            dg_bind = f"{energies['delta_g_binding']['mean']:.2f} ± {energies['delta_g_binding']['std']:.2f} {mmpbsa_res.get('unit', 'kcal/mol')}"
            details["MM-PBSA ΔG Total (Ligação)"] = dg_bind
            details["Sumário JSON"] = str(working_dir / "mmpbsa_summary.json")

        # Notificação por E-mail
        notifier.send_email_alert(
            step_name="Pós-processamento e MM-PBSA da DM",
            status="success",
            duration_seconds=duration,
            details=details,
            console_logger=console,
        )

    except md_prep.DependencyError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Pós-processamento e MM-PBSA da DM",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro de Dependência: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro de Dependência:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except md_prep.SimulationPrepError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Pós-processamento e MM-PBSA da DM",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro no Pós-processamento: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro no Pós-processamento:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Pós-processamento e MM-PBSA da DM",
            status="error",
            duration_seconds=duration,
            error_message=f"Arquivo Não Encontrado: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Arquivo Não Encontrado:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except Exception as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Pós-processamento e MM-PBSA da DM",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro Inesperado: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro Inesperado:[/bold red]\n{e}")
        raise typer.Exit(code=1)


@app.command(name="md-prep")
def md_prep_command(
    receptor: Path = typer.Option(
        ..., "--receptor", help="Caminho para o PDB original da proteína"
    ),
    sdf: Path = typer.Option(
        ..., "--sdf", help="Caminho para o arquivo docked_poses.sdf gerado no docking"
    ),
    out: Path = typer.Option(
        ..., "--out", help="Diretório de saída para a Dinâmica Molecular"
    ),
):
    """
    PREPARAÇÃO E MINIMIZAÇÃO DE ENERGIA DE DINÂMICA MOLECULAR:
    Prepara o receptor e o ligante, combina suas topologias e executa a minimização de energia no GROMACS.
    """
    console.print(
        Panel.fit(
            f"[bold blue]Preparação e Minimização de Energia de Dinâmica Molecular (GROMACS)[/bold blue]\n"
            f"Receptor: {receptor}\n"
            f"Ligante (SDF): {sdf}\n"
            f"Diretório de Saída: {out}",
            border_style="blue",
        )
    )

    start_time = time.time()
    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            tasks = {
                "A": progress.add_task(
                    description="[A] Cura do receptor com PDBFixer",
                    total=1,
                    start=False,
                ),
                "B": progress.add_task(
                    description="[B] Extração do Ligante com RDKit",
                    total=1,
                    start=False,
                ),
                "C": progress.add_task(
                    description="[C] Parametrização do Ligante (ACPYPE)",
                    total=1,
                    start=False,
                ),
                "D": progress.add_task(
                    description="[D] Topologia da Proteína (pdb2gmx)",
                    total=1,
                    start=False,
                ),
                "E": progress.add_task(
                    description="[E] Fusão de Coordenadas (complex.gro)",
                    total=1,
                    start=False,
                ),
                "F": progress.add_task(
                    description="[F] Fusão de Topologia (Stitching)",
                    total=1,
                    start=False,
                ),
                "G": progress.add_task(
                    description="[G] Definição da Caixa de Simulação (editconf)",
                    total=1,
                    start=False,
                ),
                "H": progress.add_task(
                    description="[H] Solvatação do Sistema (solvate)",
                    total=1,
                    start=False,
                ),
                "I": progress.add_task(
                    description="[I] Compilação de Íons (grompp)", total=1, start=False
                ),
                "J": progress.add_task(
                    description="[J] Neutralização e Concentração (genion)",
                    total=1,
                    start=False,
                ),
                "K": progress.add_task(
                    description="[K] Grompp Definitivo (em.tpr)", total=1, start=False
                ),
                "L": progress.add_task(
                    description="[L] Minimização de Energia (mdrun)",
                    total=1,
                    start=False,
                ),
            }

            for step, status in md_prep.prepare_md_system(receptor, sdf, out):
                task_id = tasks[step]
                if status == "start":
                    progress.start_task(task_id)
                elif status == "success":
                    progress.update(task_id, completed=1)

        duration = time.time() - start_time
        console.print(
            "\n[bold green]✓ Preparação e Minimização de Energia concluídas com sucesso![/bold green]"
        )
        console.print(f"Arquivos e outputs gerados em: [cyan]{out}[/cyan]")

        notifier.send_email_alert(
            step_name="Preparação de DM e Minimização de Energia",
            status="success",
            duration_seconds=duration,
            details={
                "Receptor": str(receptor),
                "Ligante": str(sdf),
                "Diretório de Saída": str(out),
                "Minimização (EM)": "em.gro e topol.top gerados com sucesso",
            },
            console_logger=console,
        )

    except md_prep.DependencyError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Preparação de DM e Minimização de Energia",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro de Dependência: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro de Dependência:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except md_prep.SimulationPrepError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Preparação de DM e Minimização de Energia",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro de Preparação: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro de Preparação:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Preparação de DM e Minimização de Energia",
            status="error",
            duration_seconds=duration,
            error_message=f"Arquivo Não Encontrado: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Arquivo Não Encontrado:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except Exception as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name="Preparação de DM e Minimização de Energia",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro Inesperado: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro Inesperado:[/bold red]\n{e}")
        raise typer.Exit(code=1)


@app.command(name="test-email")
def test_email_command():
    """
    TESTAR CONFIGURAÇÃO DE E-MAIL:
    Envia um e-mail de teste para validar o servidor SMTP (porta 465 SSL) e credenciais do .env.
    """
    console.print(
        Panel.fit(
            "[bold blue]Teste de Configuração de E-mail de Alerta (SMTP_SSL)[/bold blue]",
            border_style="blue",
        )
    )
    success, msg = notifier.test_email_connection(console_logger=console)
    if success:
        console.print(
            "\n[bold green]✓ Teste concluído com sucesso! Verifique sua caixa de entrada.[/bold green]"
        )
    else:
        console.print(f"\n[bold red]✗ Falha no teste:[/bold red] {msg}")
        raise typer.Exit(code=1)


@app.command(name="report")
def report_command(
    work_dir: Path = typer.Option(
        ...,
        "--dir",
        help="Diretório de trabalho contendo os artefatos do pipeline (docking, PLIP, ADMET, DM)",
    ),
    output_html: Path = typer.Option(
        None,
        "--out",
        help="Caminho personalizado para o relatório HTML (padrão: report.html dentro do diretório)",
    ),
    receptor: str = typer.Option(
        None,
        "--receptor",
        help="Código ou nome do receptor (ex: 7CFN, GPBAR1)",
    ),
    ligand: str = typer.Option(
        None,
        "--ligand",
        help="Nome ou identificador do ligante (ex: INT-777, Desoxicolato)",
    ),
):
    """
    RELATÓRIO EXECUTIVO (HTML) E VISUALIZAÇÃO 3D (PyMOL):
    1. Consolida resultados de docking (Vina), PLIP, ADMET, Gráficos de DM e MM-PBSA em um relatório HTML moderno e autoconferível.
    2. Gera o script PyMOL (.pml) automatizado para visualização 3D do complexo ligante-receptor e suas interações.
    """
    work_dir = Path(work_dir)
    console.print(
        Panel.fit(
            f"[bold blue]Geração de Relatório Executivo e Visualização 3D[/bold blue]\n"
            f"Diretório de Trabalho: {work_dir}\n"
            f"Receptor: {receptor or 'Não informado'} | Ligante: {ligand or 'Não informado'}",
            border_style="blue",
        )
    )

    if not work_dir.exists():
        console.print(
            f"[bold red]Erro:[/bold red] Diretório de trabalho não encontrado: {work_dir}"
        )
        raise typer.Exit(code=1)

    try:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=console,
        ) as progress:
            task_html = progress.add_task(
                description="[1/2] Compilando e gerando Relatório Executivo HTML...",
                total=1,
            )
            html_path, missing = report.generate_html_report(
                work_dir, output_html, receptor_name=receptor, ligand_name=ligand
            )
            progress.update(task_html, completed=1)

            task_pymol = progress.add_task(
                description="[2/2] Gerando Script de Visualização PyMOL (.pml)...",
                total=1,
            )
            try:
                pml_path = visualization.generate_pymol_script(work_dir)
            except Exception as e:
                pml_path = None
                missing.append(f"PyMOL Script: {e}")
            progress.update(task_pymol, completed=1)

        console.print(
            "\n[bold green]✓ Relatório e Script gerados com sucesso![/bold green]"
        )
        console.print(f"  • [bold]Relatório HTML:[/bold] [cyan]{html_path}[/cyan]")
        if pml_path:
            console.print(
                f"  • [bold]Script PyMOL (3D):[/bold] [cyan]{pml_path}[/cyan]"
            )
            console.print(
                "    [dim]Para abrir a cena 3D, execute no terminal: pymol show_complex.pml[/dim]"
            )

        if missing:
            console.print("\n[bold yellow]Avisos / Artefatos Ausentes:[/bold yellow]")
            for item in missing:
                console.print(f"  [yellow]⚠[/yellow] {item}")

    except Exception as e:
        console.print(
            f"\n[bold red]FATAL ERROR ao gerar relatório/script:[/bold red] {e}"
        )
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
