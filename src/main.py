#!/usr/bin/env python

import json
import time
from pathlib import Path

try:
    import questionary
except ImportError:
    questionary = None  # type: ignore[assignment]
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
    target: str = typer.Option(
        None, "--target", help="Identificador único do alvo (ex: 2AGV, 7CFN, 4HG7)"
    ),
):
    """
    TRIAGEM VIRTUAL (VIRTUAL SCREENING):
    Executa o docking de um novo xenobiótico em um receptor já preparado em coordenadas específicas.
    Aplica isolamento estrito de diretório por alvo: data/screening/<PDB_ID>/<LIGAND_NAME>/
    """
    global VINA_BIN
    if VINA_BIN is None:
        VINA_BIN = get_vina_bin()

    # Identificação do alvo e do ligante com isolamento estrito
    receptor = Path(receptor)
    ligand = Path(ligand)

    if target:
        target_id = md_prep.sanitize_target_id(target)
    else:
        if receptor.parent.name in ("processed", "raw", "results") and receptor.parent.parent.name not in ("data", "screening", ""):
            target_id = receptor.parent.parent.name
        elif receptor.parent.name not in ("data", "screening", "processed", "results", "temp", "tmp", ""):
            target_id = receptor.parent.name
        else:
            stem_clean = receptor.stem
            for sfx in ["_receptor", "_prepared", "_clean", "_docked", "_complex", "_target"]:
                if stem_clean.lower().endswith(sfx):
                    stem_clean = stem_clean[:-len(sfx)]
            target_id = stem_clean

        target_id = md_prep.sanitize_target_id(target_id)
        if target_id in ("RECEPTOR", "PROTEIN", "TARGET", "COMPLEX") and receptor.parent.parent.name not in ("data", "screening", "", ".", "/"):
            target_id = md_prep.sanitize_target_id(receptor.parent.parent.name)

    ligand_name = ligand.stem

    results_dir = DATA_DIR / "screening" / target_id / ligand_name
    results_dir.mkdir(parents=True, exist_ok=True)

    console.print(
        Panel.fit(
            f"[bold cyan]Triagem Virtual (Target Isolation)[/bold cyan]\n"
            f"Alvo / Receptor: {target_id} ({receptor.name}) | Ligante: {ligand_name}\n"
            f"Diretório Exclusivo: {results_dir}\n"
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
            # Execução do Vina com prefixo do alvo
            task1 = progress.add_task(
                description=f"Rodando Vina para {target_id} + {ligand_name}...", total=1
            )
            docked_out = results_dir / f"{target_id}_{ligand_name}_docked.pdbqt"
            vina_log = results_dir / f"{target_id}_{ligand_name}_vina.log"

            vina_runner.run_vina(
                VINA_BIN,
                receptor,
                ligand,
                box_params,
                docked_out,
                vina_log,
                exhaustiveness,
            )
            # Cria espelhos
            import shutil
            shutil.copy2(docked_out, results_dir / "docked.pdbqt")
            shutil.copy2(vina_log, results_dir / "vina.log")
            progress.update(task1, completed=1)

            # Extração de score
            task2 = progress.add_task(description="Extraindo score...", total=1)
            score = analysis.extract_vina_score(vina_log)
            progress.update(task2, completed=1)

            # Exporta para SDF e executa o fluxo do PLIP
            task_plip = progress.add_task(
                description="Executando PLIP (Docker)...", total=1
            )
            sdf_out = results_dir / f"{target_id}_{ligand_name}_docked_poses.sdf"
            import subprocess

            exec_name = get_executable("mk_export")
            subprocess.run([exec_name, str(docked_out), "-s", str(sdf_out)], check=True)
            shutil.copy2(sdf_out, results_dir / "docked_poses.sdf")

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

            complex_pdb = results_dir / f"{target_id}_{ligand_name}_complex.pdb"
            analysis.generate_complex_pdb(receptor_pdb, sdf_out, complex_pdb)
            shutil.copy2(complex_pdb, results_dir / "complex.pdb")

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

            # Salva o arquivo JSON consolidado com prefixo e espelho
            interactions_file = results_dir / f"{target_id}_{ligand_name}_interactions.json"
            with open(interactions_file, "w", encoding="utf-8") as f:
                json.dump(interactions, f, indent=4, ensure_ascii=False)
            shutil.copy2(interactions_file, results_dir / "interactions.json")

            progress.update(task_plip, completed=1)

        duration = time.time() - start_time
        console.print(
            f"\n[bold green]✓ Triagem concluída para {target_id} + {ligand_name}![/bold green]"
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
            "Alvo": target_id,
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
            step_name=f"Triagem Virtual ({target_id} + {ligand_name})",
            status="success",
            duration_seconds=duration,
            details=details,
            console_logger=console,
        )

    except Exception as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name=f"Triagem Virtual ({target_id} + {ligand_name})",
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
    if questionary is None:
        console.print(
            "[bold red]Erro ao carregar a interface interativa:[/bold red] O pacote 'questionary' / 'prompt_toolkit' não está disponível.\n"
            "[yellow]Execute o pipeline através do gerenciador de pacotes do projeto:[/yellow]\n"
            "  [cyan]uv run src/main.py interactive[/cyan]\n"
            "Ou instale as dependências com: [cyan]uv sync[/cyan]"
        )
        raise typer.Exit(code=1)

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
                "7. Compilar TPR de Produção & Pacote para Cluster (SSH/tmux)",
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
            rec_default = ""
            for candidate in [
                Path("data/7CFN/processed/receptor.pdbqt"),
                Path("data/2AGV/processed/receptor.pdbqt"),
                Path("data/1OSV/processed/receptor.pdbqt"),
            ]:
                if candidate.exists():
                    rec_default = str(candidate)
                    break
            if not rec_default:
                found_rec = list(Path("data").glob("*/processed/receptor.pdbqt"))
                if found_rec:
                    rec_default = str(found_rec[0])

            receptor = questionary.path(
                "Caminho para o receptor (.pdbqt):",
                default=rec_default,
            ).ask()
            if not receptor:
                console.print("[bold red]Operação cancelada: caminho do receptor não fornecido.[/bold red]")
                continue

            rec_p = Path(receptor)
            target_default = "TARGET"
            if rec_p.parent.name in ("processed", "raw", "results") and rec_p.parent.parent.name not in ("data", "screening", ""):
                target_default = rec_p.parent.parent.name
            elif rec_p.parent.name not in ("data", "screening", "processed", "results", "temp", "tmp", ""):
                target_default = rec_p.parent.name
            else:
                stem_clean = rec_p.stem
                for sfx in ["_receptor", "_prepared", "_clean", "_docked", "_complex", "_target"]:
                    if stem_clean.lower().endswith(sfx):
                        stem_clean = stem_clean[:-len(sfx)]
                if stem_clean.lower() not in ("receptor", "protein", "target", "complex"):
                    target_default = stem_clean
            target_default = md_prep.sanitize_target_id(target_default)

            target = questionary.text(
                "Identificador único do Alvo (Target ID, ex: 2AGV, 7CFN):",
                default=target_default,
            ).ask()
            target = md_prep.sanitize_target_id(target) if target else target_default

            lig_default = ""
            found_ligs = list(Path("data").glob("*.pdbqt"))
            if found_ligs:
                lig_default = str(found_ligs[0])

            ligand = questionary.path(
                "Caminho para o ligante (.pdbqt):",
                default=lig_default,
            ).ask()
            if not ligand:
                console.print("[bold red]Operação cancelada: caminho do ligante não fornecido.[/bold red]")
                continue

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
                target=target,
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
                Path("data/screening/7CFN/desoxicolato/7CFN_desoxicolato_docked_poses.sdf"),
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

            receptor_path = questionary.path(
                "Caminho para o PDB original da proteína (Receptor):",
                default=rec_default,
            ).ask()

            sdf_path = questionary.path(
                "Caminho para o arquivo docked_poses.sdf (Ligante):",
                default=sdf_default,
            ).ask()

            # Inferência inteligente de Target ID
            target_default = "7CFN"
            if receptor_path:
                rec_p = Path(receptor_path)
                if rec_p.parent.name == "processed" and rec_p.parent.parent.name not in ("data", ""):
                    target_default = rec_p.parent.parent.name
                else:
                    target_default = rec_p.stem.replace("_receptor", "").replace("_clean", "")
            target_default = md_prep.sanitize_target_id(target_default)

            target_id = questionary.text(
                "Identificador único do Alvo (ex: 7CFN, 1OSV, 4HG7):",
                default=target_default,
            ).ask()
            target_id = md_prep.sanitize_target_id(target_id)

            out_default = f"data/md_files/{target_id}"
            output_dir = questionary.text(
                "Diretório exclusivo de saída para a Dinâmica Molecular (Target Isolation):",
                default=out_default,
            ).ask()

            if not receptor_path or not sdf_path or not output_dir or not target_id:
                console.print(
                    "[bold red]Operação cancelada: todos os campos devem ser preenchidos.[/bold red]"
                )
                continue

            out_path = Path(output_dir)
            # Verificação e sanitização de resíduos pré-execução
            if out_path.exists():
                stale_files = md_prep.check_and_purge_stale_files(out_path, purge=False)
                if stale_files:
                    console.print(
                        f"[bold yellow]Aviso de Segurança:[/bold yellow] Foram detectados {len(stale_files)} arquivos residuais em [cyan]{out_path}[/cyan]."
                    )
                    purge_confirm = questionary.confirm(
                        "Deseja purgar com segurança os arquivos residuais antigos para evitar contaminação cruzada?",
                        default=True,
                    ).ask()
                    if purge_confirm:
                        md_prep.check_and_purge_stale_files(out_path, purge=True)
                        console.print("[bold green]✓ Limpeza de resíduos concluída com sucesso.[/bold green]")

            start_time = time.time()
            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Preparação e Minimização de Energia de Dinâmica Molecular (GROMACS)[/bold blue]\n"
                        f"Alvo / Target ID: [green]{target_id}[/green]\n"
                        f"Receptor: {receptor_path}\n"
                        f"Ligante (SDF): {sdf_path}\n"
                        f"Diretório Exclusivo: {output_dir}",
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
                            description=f"[E] Fusão de Coordenadas ({target_id}_complex.gro)",
                            total=1,
                            start=False,
                        ),
                        "F": progress.add_task(
                            description=f"[F] Fusão de Topologia ({target_id}_topol.top)",
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
                            description=f"[K] Grompp Definitivo ({target_id}_em.tpr)",
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
                        Path(receptor_path), Path(sdf_path), out_path, target_id=target_id, purge=False
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
                    f"Arquivos e outputs isolados em: [cyan]{output_dir}[/cyan]"
                )

                notifier.send_email_alert(
                    step_name=f"Opção 5: Preparação de DM ({target_id})",
                    status="success",
                    duration_seconds=duration,
                    details={
                        "Alvo": target_id,
                        "Receptor": str(receptor_path),
                        "Ligante": str(sdf_path),
                        "Diretório de Saída": str(output_dir),
                        "Minimização (EM)": f"{target_id}_em.gro e {target_id}_topol.top gerados com sucesso",
                        "Próxima Ação": "Rodar Equilíbrio NVT/NPT (Opção 6)",
                    },
                    console_logger=console,
                )

            except md_prep.DependencyError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name=f"Opção 5: Preparação de DM ({target_id})",
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
                    step_name=f"Opção 5: Preparação de DM ({target_id})",
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
                    step_name=f"Opção 5: Preparação de DM ({target_id})",
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
                    step_name=f"Opção 5: Preparação de DM ({target_id})",
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
            candidate_dirs = [str(d) for d in (DATA_DIR / "md_files").glob("*") if d.is_dir()]
            md_dir_default = candidate_dirs[0] if candidate_dirs else "data/md_files/7CFN"

            md_dir = questionary.text(
                "Diretório de trabalho da Dinâmica Molecular (onde contêm em.gro e topol.top):",
                default=md_dir_default,
            ).ask()

            if not md_dir:
                console.print(
                    "[bold red]Operação cancelada: o diretório de trabalho deve ser preenchido.[/bold red]"
                )
                continue

            md_path = Path(md_dir)
            target_id = md_prep.sanitize_target_id(md_path.name)
            if target_id.lower() in ("md_files", "screening", "data"):
                candidates = list(md_path.glob("*_em.gro")) or list(md_path.glob("*_topol.top"))
                if candidates:
                    target_id = candidates[0].stem.replace("_em", "").replace("_topol", "")
            target_id = md_prep.sanitize_target_id(target_id)

            start_time = time.time()
            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Equilíbrio Termodinâmico da Dinâmica Molecular (GROMACS)[/bold blue]\n"
                        f"Alvo / Target ID: [green]{target_id}[/green]\n"
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
                            description=f"[A] Geração de Grupos e Índices ({target_id}_index.ndx)",
                            total=1,
                            start=False,
                        ),
                        "B": progress.add_task(
                            description=f"[B] Compilação da Caixa NVT ({target_id}_nvt.tpr)",
                            total=1,
                            start=False,
                        ),
                        "C": progress.add_task(
                            description=f"[C] Execução do Equilíbrio NVT ({target_id}_nvt)",
                            total=1,
                            start=False,
                        ),
                        "D": progress.add_task(
                            description=f"[D] Compilação da Caixa NPT ({target_id}_npt.tpr)",
                            total=1,
                            start=False,
                        ),
                        "E": progress.add_task(
                            description=f"[E] Execução do Equilíbrio NPT ({target_id}_npt)",
                            total=1,
                            start=False,
                        ),
                    }

                    for step, status in md_equil.run_md_equilibration(md_path, target_id=target_id):
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
                    step_name=f"Opção 6: Equilíbrio da Dinâmica ({target_id})",
                    status="success",
                    duration_seconds=duration,
                    details={
                        "Alvo": target_id,
                        "Diretório de Trabalho": str(md_dir),
                        "Equilíbrio NVT": f"Concluído ({target_id}_nvt.gro gerado)",
                        "Equilíbrio NPT": f"Concluído ({target_id}_npt.gro gerado)",
                        "Próxima Ação": "Compilar TPR & Exportar Pacote (Opção 7) ou Rodar Produção (Opção 8)",
                    },
                    console_logger=console,
                )

            except md_prep.DependencyError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name=f"Opção 6: Equilíbrio da Dinâmica ({target_id})",
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
                    step_name=f"Opção 6: Equilíbrio da Dinâmica ({target_id})",
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
                    step_name=f"Opção 6: Equilíbrio da Dinâmica ({target_id})",
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
                    step_name=f"Opção 6: Equilíbrio da Dinâmica ({target_id})",
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

        elif choice == "7. Compilar TPR de Produção & Pacote para Cluster (SSH/tmux)":
            candidate_dirs = [str(d) for d in (DATA_DIR / "md_files").glob("*") if d.is_dir()]
            md_dir_default = candidate_dirs[0] if candidate_dirs else "data/md_files/7CFN"

            md_dir = questionary.text(
                "Diretório de trabalho da Dinâmica Molecular (onde contêm npt.gro e topol.top):",
                default=md_dir_default,
            ).ask()

            if not md_dir:
                console.print(
                    "[bold red]Operação cancelada: o diretório de trabalho deve ser preenchido.[/bold red]"
                )
                continue

            md_path = Path(md_dir)
            target_id = md_prep.sanitize_target_id(md_path.name)
            if target_id.lower() in ("md_files", "screening", "data"):
                candidates = list(md_path.glob("*_npt.gro")) or list(md_path.glob("*_topol.top"))
                if candidates:
                    target_id = candidates[0].stem.replace("_npt", "").replace("_topol", "")
            target_id = md_prep.sanitize_target_id(target_id)

            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Compilação de Produção & Exportação de Pacote para Cluster[/bold blue]\n"
                        f"Alvo / Target ID: [green]{target_id}[/green]\n"
                        f"Diretório de Trabalho: {md_dir}",
                        border_style="blue",
                    )
                )

                console.print(
                    f"[yellow]Compilando {target_id}_md.tpr via GROMACS (grompp) e validando integridade...[/yellow]"
                )
                tpr_path = md_analysis.compile_production_tpr(md_path, target_id=target_id)
                console.print(
                    f"[bold green]✓ Arquivo '{tpr_path.name}' gerado e validado com sucesso em:[/bold green] [cyan]{tpr_path}[/cyan]"
                )

                # Exporta pacote standalone para cluster/SSH
                export_dir = md_analysis.export_cluster_package(md_path, target_id=target_id)
                console.print(
                    f"\n[bold green]✓ Pacote Modular para Cluster exportado com sucesso em:[/bold green] [cyan]{export_dir}[/cyan]"
                )
                console.print(
                    Panel(
                        f"[bold cyan]Instruções para Execução em Servidor/Cluster (SSH / tmux):[/bold cyan]\n\n"
                        f"1. Envie a pasta do pacote para seu servidor remoto:\n"
                        f"   [yellow]rsync -avP cluster_export/{target_id}/ user@cluster:/path/to/simulations/{target_id}/[/yellow]\n\n"
                        f"2. Conecte-se ao servidor e abra uma sessão tmux persistente:\n"
                        f"   [yellow]ssh user@cluster[/yellow]\n"
                        f"   [yellow]tmux new -s md_{target_id}[/yellow]\n"
                        f"   [yellow]cd /path/to/simulations/{target_id}[/yellow]\n\n"
                        f"3. Inicie a produção (com detecção automática de GPU e auto-retomada):\n"
                        f"   [yellow]./run_local.sh[/yellow]\n\n"
                        f"4. Desanexe da sessão com [bold]Ctrl+B[/bold], depois [bold]D[/bold].\n"
                        f"   Para acompanhar os logs a qualquer momento: [yellow]tail -f {target_id}_md.log[/yellow]",
                        title=f"Manual de Execução Remota - {target_id}",
                        border_style="green",
                    )
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
                        f"[bold red]Erro na Compilação do TPR:[/bold red]\n{e}",
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
            candidate_dirs = [str(d) for d in (DATA_DIR / "md_files").glob("*") if d.is_dir()]
            md_dir_default = candidate_dirs[0] if candidate_dirs else "data/md_files/7CFN"

            md_dir = questionary.text(
                "Diretório de trabalho da Dinâmica Molecular (onde contêm npt.gro e topol.top):",
                default=md_dir_default,
            ).ask()

            if not md_dir:
                console.print(
                    "[bold red]Operação cancelada: o diretório de trabalho deve ser preenchido.[/bold red]"
                )
                continue

            md_path = Path(md_dir)
            target_id = md_prep.sanitize_target_id(md_path.name)
            if target_id.lower() in ("md_files", "screening", "data"):
                candidates = list(md_path.glob("*_npt.gro")) or list(md_path.glob("*_topol.top"))
                if candidates:
                    target_id = candidates[0].stem.replace("_npt", "").replace("_topol", "")
            target_id = md_prep.sanitize_target_id(target_id)

            start_time = time.time()
            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Produção de Dinâmica Molecular (GROMACS)[/bold blue]\n"
                        f"Alvo / Target ID: [green]{target_id}[/green]\n"
                        f"Diretório de Trabalho: {md_dir}",
                        border_style="blue",
                    )
                )

                console.print(
                    f"[yellow]Iniciando a produção da Dinâmica Molecular para {target_id} (grompp & mdrun)...[/yellow]"
                )
                md_analysis.run_production_md(md_path, target_id=target_id)
                duration = time.time() - start_time
                console.print(
                    f"[bold green]✓ Produção concluída com sucesso para {target_id}![/bold green]"
                )
                console.print(
                    "[cyan]Execute a opção 9 para tratamento de PBC, gráficos e MM-PBSA.[/cyan]"
                )

                notifier.send_email_alert(
                    step_name=f"Opção 8: Produção de Dinâmica Molecular ({target_id})",
                    status="success",
                    duration_seconds=duration,
                    details={
                        "Alvo": target_id,
                        "Diretório de Trabalho": str(md_dir),
                        "Trajetória Gerada": f"{target_id}_md.xtc e {target_id}_md.gro gerados",
                        "Próxima Ação Recomendada": "Executar Pós-processamento e MM-PBSA (Opção 9)",
                    },
                    console_logger=console,
                )

            except md_prep.DependencyError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name=f"Opção 8: Produção de Dinâmica Molecular ({target_id})",
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
                    step_name=f"Opção 8: Produção de Dinâmica Molecular ({target_id})",
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
                    step_name=f"Opção 8: Produção de Dinâmica Molecular ({target_id})",
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
                    step_name=f"Opção 8: Produção de Dinâmica Molecular ({target_id})",
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
            candidate_dirs = [str(d) for d in (DATA_DIR / "md_files").glob("*") if d.is_dir()]
            md_dir_default = candidate_dirs[0] if candidate_dirs else "data/md_files/7CFN"

            md_dir = questionary.text(
                "Diretório de trabalho da Dinâmica Molecular (onde contêm md.tpr e md.xtc):",
                default=md_dir_default,
            ).ask()

            if not md_dir:
                console.print(
                    "[bold red]Operação cancelada: o diretório de trabalho deve ser preenchido.[/bold red]"
                )
                continue

            md_path = Path(md_dir)
            target_id = md_prep.sanitize_target_id(md_path.name)
            if target_id.lower() in ("md_files", "screening", "data"):
                candidates = list(md_path.glob("*_md.tpr")) or list(md_path.glob("*_md.xtc"))
                if candidates:
                    target_id = candidates[0].stem.replace("_md", "")
            target_id = md_prep.sanitize_target_id(target_id)

            run_mmpbsa = questionary.confirm(
                "Deseja executar o cálculo de Energia Livre de Ligação MM-PBSA (Janela: 60 - 100 ns / Últimos 40%)?",
                default=True,
            ).ask()

            start_time = time.time()
            try:
                console.print(
                    Panel.fit(
                        f"[bold blue]Pós-processamento, Gráficos e MM-PBSA da Dinâmica Molecular[/bold blue]\n"
                        f"Alvo / Target ID: [green]{target_id}[/green]\n"
                        f"Diretório de Trabalho: {md_dir}\n"
                        f"[dim]Protocolo: Dupla Escala Temporal (Estrutural: 0 - 100 ns | Termodinâmica: 60 - 100 ns)[/dim]",
                        border_style="blue",
                    )
                )

                with Progress(
                    SpinnerColumn(),
                    TextColumn("[progress.description]{task.description}"),
                    console=console,
                ) as progress:
                    task_pbc = progress.add_task(
                        description=f"[1/4] Tratamento Automatizado de Trajetória ({target_id}_md_fit.xtc, {target_id}_md_clean.gro)...",
                        total=1,
                    )
                    md_analysis.fix_pbc(md_path, target_id=target_id)
                    progress.update(task_pbc, completed=1)

                    task_traj = progress.add_task(
                        description="[2/4] Análises Estruturais Globais (0 - 100 ns: RMSD, RMSF, HBond, Rg, SASA, Clustering & CSVs)...",
                        total=1,
                    )
                    md_analysis.analyze_trajectory(md_path, target_id=target_id)
                    progress.update(task_traj, completed=1)

                    task_plot = progress.add_task(
                        description="[3/4] Geração de Gráficos Científicos de Publicação (Janela Completa: 0 - 100 ns, 300 DPI)...",
                        total=1,
                    )
                    plots = md_analysis.plot_md_results(md_path, target_id=target_id)
                    progress.update(task_plot, completed=1)

                    mmpbsa_res = None
                    if run_mmpbsa:
                        task_mmpbsa = progress.add_task(
                            description="[4/4] Cálculo de Energia Livre MM-PBSA (Janela: 60 - 100 ns / Últimos 40% - Estado Estacionário)...",
                            total=1,
                        )
                        mmpbsa_res = md_analysis.calculate_mmpbsa(md_path, target_id=target_id)
                        progress.update(task_mmpbsa, completed=1)

                duration = time.time() - start_time
                console.print(
                    f"\n[bold green]✓ Pós-processamento concluído com sucesso para {target_id}![/bold green]"
                )

                if plots:
                    console.print(
                        "\n[bold cyan]Gráficos Científicos Gerados (Janela: 0 - 100 ns):[/bold cyan]"
                    )
                    for p_name, p_path in plots.items():
                        console.print(
                            f"  • [bold]{p_name.upper()}:[/bold] [green]{p_path}[/green]"
                        )

                details = {
                    "Alvo": target_id,
                    "Diretório": str(md_dir),
                    "Tratamento PBC": f"Concluído ({target_id}_md_fit.xtc e {target_id}_md_clean.gro gerados)",
                    "Análise Estrutural": "Janela Completa 0 - 100 ns (RMSD Backbone + Ligante, RMSF, HBond, Rg, SASA)",
                    "Gráficos Gerados": f"{len(plots)} gráficos salvos (300 DPI)"
                    if plots
                    else "Nenhum",
                }

                if mmpbsa_res and "energies" in mmpbsa_res:
                    table = Table(
                        title=f"Resumo Termodinâmico MM-PBSA ({target_id}) [Janela: 60 - 100 ns (Últimos 40%)]",
                        show_header=True,
                        header_style="bold magenta",
                    )
                    table.add_column("Componente Energético", style="dim", width=34)
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
                        f"Sumário JSON salvo em: [cyan]{md_path / f'{target_id}_mmpbsa_summary.json'}[/cyan]"
                    )

                    dg_bind = f"{energies['delta_g_binding']['mean']:.2f} ± {energies['delta_g_binding']['std']:.2f} {mmpbsa_res.get('unit', 'kcal/mol')}"
                    details["MM-PBSA Janela"] = "60 - 100 ns (Últimos 40% - Estado Estacionário)"
                    details["MM-PBSA ΔG Total (Ligação)"] = dg_bind
                    details["Sumário JSON"] = str(md_path / f"{target_id}_mmpbsa_summary.json")

                # Envio do E-mail de Alerta
                notifier.send_email_alert(
                    step_name=f"Opção 9: Pós-processamento e MM-PBSA da DM ({target_id})",
                    status="success",
                    duration_seconds=duration,
                    details=details,
                    console_logger=console,
                )

            except md_prep.DependencyError as e:
                duration = time.time() - start_time
                notifier.send_email_alert(
                    step_name=f"Opção 9: Pós-processamento e MM-PBSA da DM ({target_id})",
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
                    step_name=f"Opção 9: Pós-processamento e MM-PBSA da DM ({target_id})",
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
                    step_name=f"Opção 9: Pós-processamento e MM-PBSA da DM ({target_id})",
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
                    step_name=f"Opção 9: Pós-processamento e MM-PBSA da DM ({target_id})",
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
            default_dirs = [
                "data/md_files/7CFN",
                "data/screening/7CFN/desoxicolato",
                "data/screening/desoxicolato",
                "data/1OSV/results",
                "data/md_files",
                "data",
            ]
            default_dir = next((d for d in default_dirs if Path(d).exists()), "data")

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


@app.command(name="md-prep")
def md_prep_command(
    receptor: Path = typer.Option(
        ..., "--receptor", help="Caminho para o PDB original da proteína"
    ),
    sdf: Path = typer.Option(
        ..., "--sdf", help="Caminho para o arquivo docked_poses.sdf gerado no docking"
    ),
    out: Path = typer.Option(
        None, "--out", help="Diretório de saída para a Dinâmica Molecular (padrão: data/md_files/<target_id>)"
    ),
    target: str = typer.Option(
        None, "--target", help="Identificador único do alvo (ex: 7CFN, 1OSV, 4HG7)"
    ),
    purge: bool = typer.Option(
        True, "--purge/--no-purge", help="Purgar arquivos residuais antigos do diretório de saída antes de iniciar"
    ),
):
    """
    PREPARAÇÃO E MINIMIZAÇÃO DE ENERGIA DE DINÂMICA MOLECULAR:
    Prepara o receptor e o ligante, combina suas topologias e executa a minimização de energia no GROMACS.
    Gera arquivos isolados com prefixo do alvo (<target_id>_*.gro, <target_id>_*.top).
    """
    receptor = Path(receptor)
    sdf = Path(sdf)

    if not target:
        if receptor.parent.name == "processed" and receptor.parent.parent.name not in ("data", ""):
            target = receptor.parent.parent.name
        else:
            target = receptor.stem.replace("_receptor", "").replace("_clean", "")
    target = md_prep.sanitize_target_id(target)

    if out is None:
        out = DATA_DIR / "md_files" / target
    out = Path(out)

    console.print(
        Panel.fit(
            f"[bold blue]Preparação e Minimização de Energia de Dinâmica Molecular (GROMACS)[/bold blue]\n"
            f"Alvo / Target ID: [green]{target}[/green]\n"
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
                    description=f"[E] Fusão de Coordenadas ({target}_complex.gro)",
                    total=1,
                    start=False,
                ),
                "F": progress.add_task(
                    description=f"[F] Fusão de Topologia ({target}_topol.top)",
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
                    description=f"[K] Grompp Definitivo ({target}_em.tpr)", total=1, start=False
                ),
                "L": progress.add_task(
                    description="[L] Minimização de Energia (mdrun)",
                    total=1,
                    start=False,
                ),
            }

            for step, status in md_prep.prepare_md_system(receptor, sdf, out, target_id=target, purge=purge):
                task_id = tasks[step]
                if status == "start":
                    progress.start_task(task_id)
                elif status == "success":
                    progress.update(task_id, completed=1)

        duration = time.time() - start_time
        console.print(
            f"\n[bold green]✓ Preparação e Minimização de Energia concluídas com sucesso para {target}![/bold green]"
        )
        console.print(f"Arquivos e outputs isolados em: [cyan]{out}[/cyan]")

        notifier.send_email_alert(
            step_name=f"Preparação de DM ({target})",
            status="success",
            duration_seconds=duration,
            details={
                "Alvo": target,
                "Receptor": str(receptor),
                "Ligante": str(sdf),
                "Diretório de Saída": str(out),
                "Minimização (EM)": f"{target}_em.gro e {target}_topol.top gerados com sucesso",
            },
            console_logger=console,
        )

    except md_prep.DependencyError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name=f"Preparação de DM ({target})",
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
            step_name=f"Preparação de DM ({target})",
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
            step_name=f"Preparação de DM ({target})",
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
            step_name=f"Preparação de DM ({target})",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro Inesperado: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro Inesperado:[/bold red]\n{e}")
        raise typer.Exit(code=1)


@app.command(name="md-equil")
def md_equil_command(
    working_dir: Path = typer.Option(
        ...,
        "--dir",
        help="Diretório de trabalho contendo os arquivos da minimização (em.gro, topol.top, etc.)",
    ),
    target: str = typer.Option(
        None, "--target", help="Identificador único do alvo (ex: 7CFN, 1OSV, 4HG7)"
    ),
):
    """
    EQUILÍBRIO TERMODINÂMICO (NVT/NPT):
    Gera restrições de posição, compila e executa o equilíbrio termodinâmico NVT e NPT no GROMACS.
    """
    working_dir = Path(working_dir)
    if not target:
        target = md_prep.sanitize_target_id(working_dir.name)
        if target.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_em.gro")) or list(working_dir.glob("*_topol.top"))
            if candidates:
                target = candidates[0].stem.replace("_em", "").replace("_topol", "")
    target = md_prep.sanitize_target_id(target)

    console.print(
        Panel.fit(
            f"[bold blue]Equilíbrio Termodinâmico da Dinâmica Molecular (GROMACS)[/bold blue]\n"
            f"Alvo / Target ID: [green]{target}[/green]\n"
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
            tasks = {
                "A": progress.add_task(
                    description=f"[A] Geração de Grupos e Índices ({target}_index.ndx)",
                    total=1,
                    start=False,
                ),
                "B": progress.add_task(
                    description=f"[B] Compilação da Caixa NVT ({target}_nvt.tpr)",
                    total=1,
                    start=False,
                ),
                "C": progress.add_task(
                    description=f"[C] Execução do Equilíbrio NVT ({target}_nvt)",
                    total=1,
                    start=False,
                ),
                "D": progress.add_task(
                    description=f"[D] Compilação da Caixa NPT ({target}_npt.tpr)",
                    total=1,
                    start=False,
                ),
                "E": progress.add_task(
                    description=f"[E] Execução do Equilíbrio NPT ({target}_npt)",
                    total=1,
                    start=False,
                ),
            }

            for step, status in md_equil.run_md_equilibration(working_dir, target_id=target):
                task_id = tasks[step]
                if status == "start":
                    progress.start_task(task_id)
                elif status == "success":
                    progress.update(task_id, completed=1)

        duration = time.time() - start_time
        console.print(
            f"\n[bold green]✓ Equilíbrio NVT/NPT concluído com sucesso para {target}![/bold green]"
        )
        console.print(f"Estruturas equilibradas geradas em: [cyan]{working_dir}[/cyan]")

        notifier.send_email_alert(
            step_name=f"Equilíbrio de DM ({target})",
            status="success",
            duration_seconds=duration,
            details={
                "Alvo": target,
                "Diretório": str(working_dir),
                "Equilíbrio NVT": f"Concluído ({target}_nvt.gro gerado)",
                "Equilíbrio NPT": f"Concluído ({target}_npt.gro gerado)",
            },
            console_logger=console,
        )

    except md_prep.DependencyError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name=f"Equilíbrio de DM ({target})",
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
            step_name=f"Equilíbrio de DM ({target})",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro no GROMACS: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro no Equilíbrio:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name=f"Equilíbrio de DM ({target})",
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
            step_name=f"Equilíbrio de DM ({target})",
            status="error",
            duration_seconds=duration,
            error_message=f"Erro Inesperado: {e}",
            console_logger=console,
        )
        console.print(f"\n[bold red]Erro Inesperado:[/bold red]\n{e}")
        raise typer.Exit(code=1)


@app.command(name="md-compile")
def md_compile_command(
    working_dir: Path = typer.Option(
        ...,
        "--dir",
        help="Diretório de trabalho contendo os arquivos do equilíbrio (npt.gro, topol.top, etc.)",
    ),
    target: str = typer.Option(
        None, "--target", help="Identificador único do alvo (ex: 7CFN, 1OSV, 4HG7)"
    ),
):
    """
    COMPILAÇÃO DO ARQUIVO DE PRODUÇÃO (md.tpr) & EXPORTAÇÃO PARA CLUSTER:
    Executa o 'gmx grompp' com validação estrita de integridade e exporta o pacote modular cluster_export/<PDB_ID>/.
    """
    working_dir = Path(working_dir)
    if not target:
        target = md_prep.sanitize_target_id(working_dir.name)
        if target.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_npt.gro")) or list(working_dir.glob("*_topol.top"))
            if candidates:
                target = candidates[0].stem.replace("_npt", "").replace("_topol", "")
    target = md_prep.sanitize_target_id(target)

    console.print(
        Panel.fit(
            f"[bold blue]Compilação do Arquivo de Produção & Pacote para Cluster[/bold blue]\n"
            f"Alvo / Target ID: [green]{target}[/green]\n"
            f"Diretório de Trabalho: {working_dir}",
            border_style="blue",
        )
    )

    try:
        console.print(f"[yellow]Compilando {target}_md.tpr via GROMACS (grompp) e validando integridade...[/yellow]")
        tpr_path = md_analysis.compile_production_tpr(working_dir, target_id=target)
        console.print(
            f"[bold green]✓ Arquivo '{tpr_path.name}' gerado e validado com sucesso em:[/bold green] [cyan]{tpr_path}[/cyan]"
        )

        export_dir = md_analysis.export_cluster_package(working_dir, target_id=target)
        console.print(
            f"\n[bold green]✓ Pacote Modular para Cluster exportado com sucesso em:[/bold green] [cyan]{export_dir}[/cyan]"
        )
        console.print(
            Panel(
                f"[bold cyan]Instruções para Execução em Servidor/Cluster (SSH / tmux):[/bold cyan]\n\n"
                f"1. Envie a pasta do pacote para seu servidor remoto:\n"
                f"   [yellow]rsync -avP cluster_export/{target}/ user@cluster:/path/to/simulations/{target}/[/yellow]\n\n"
                f"2. Conecte-se ao servidor e abra uma sessão tmux persistente:\n"
                f"   [yellow]ssh user@cluster[/yellow]\n"
                f"   [yellow]tmux new -s md_{target}[/yellow]\n"
                f"   [yellow]cd /path/to/simulations/{target}[/yellow]\n\n"
                f"3. Inicie a produção (com detecção automática de GPU e auto-retomada):\n"
                f"   [yellow]./run_local.sh[/yellow]\n\n"
                f"4. Desanexe da sessão com [bold]Ctrl+B[/bold], depois [bold]D[/bold].\n"
                f"   Para acompanhar os logs a qualquer momento: [yellow]tail -f {target}_md.log[/yellow]",
                title=f"Manual de Execução Remota - {target}",
                border_style="green",
            )
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


@app.command(name="md-export")
def md_export_command(
    working_dir: Path = typer.Option(
        ...,
        "--dir",
        help="Diretório de trabalho contendo o arquivo md.tpr ou os dados da dinâmica",
    ),
    target: str = typer.Option(
        None, "--target", help="Identificador único do alvo (ex: 7CFN, 1OSV, 4HG7)"
    ),
    output_export_dir: Path = typer.Option(
        None, "--out", help="Diretório de exportação de destino (padrão: cluster_export/<target_id>)"
    ),
):
    """
    EMPACOTAMENTO MODULAR PARA CLUSTER (SSH / tmux / sem Slurm):
    Empacota o arquivo md.tpr validado, o script de execução autônomo run_local.sh e o README.md explicativo.
    """
    working_dir = Path(working_dir)
    if not target:
        target = md_prep.sanitize_target_id(working_dir.name)
        if target.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_md.tpr")) or list(working_dir.glob("*_npt.gro"))
            if candidates:
                target = candidates[0].stem.replace("_md", "").replace("_npt", "")
    target = md_prep.sanitize_target_id(target)

    console.print(
        Panel.fit(
            f"[bold blue]Exportação de Pacote Modular para Cluster (SSH / tmux)[/bold blue]\n"
            f"Alvo / Target ID: [green]{target}[/green]\n"
            f"Diretório de Trabalho: {working_dir}",
            border_style="blue",
        )
    )

    try:
        export_dir = md_analysis.export_cluster_package(
            working_dir, target_id=target, output_export_dir=output_export_dir
        )
        console.print(
            f"\n[bold green]✓ Pacote Modular exportado com sucesso em:[/bold green] [cyan]{export_dir}[/cyan]"
        )
        console.print(
            Panel(
                f"[bold cyan]Instruções para Execução em Servidor/Cluster (SSH / tmux):[/bold cyan]\n\n"
                f"1. Envie a pasta do pacote para seu servidor remoto:\n"
                f"   [yellow]rsync -avP {export_dir}/ user@cluster:/path/to/simulations/{target}/[/yellow]\n\n"
                f"2. Conecte-se ao servidor e abra uma sessão tmux persistente:\n"
                f"   [yellow]ssh user@cluster[/yellow]\n"
                f"   [yellow]tmux new -s md_{target}[/yellow]\n"
                f"   [yellow]cd /path/to/simulations/{target}[/yellow]\n\n"
                f"3. Inicie a produção (com detecção automática de GPU e auto-retomada):\n"
                f"   [yellow]./run_local.sh[/yellow]\n\n"
                f"4. Desanexe da sessão com [bold]Ctrl+B[/bold], depois [bold]D[/bold].\n"
                f"   Para acompanhar os logs a qualquer momento: [yellow]tail -f {target}_md.log[/yellow]",
                title=f"Manual de Execução Remota - {target}",
                border_style="green",
            )
        )
    except Exception as e:
        console.print(f"\n[bold red]Erro na Exportação do Pacote:[/bold red] {e}")
        raise typer.Exit(code=1)


@app.command(name="md-run")
def md_run_command(
    working_dir: Path = typer.Option(
        ...,
        "--dir",
        help="Diretório de trabalho onde estão os arquivos do equilíbrio (nvt/npt)",
    ),
    target: str = typer.Option(
        None, "--target", help="Identificador único do alvo (ex: 7CFN, 1OSV, 4HG7)"
    ),
):
    """
    PRODUÇÃO DE DINÂMICA MOLECULAR:
    Compila e executa a produção da Dinâmica Molecular no GROMACS (grompp e mdrun).
    """
    working_dir = Path(working_dir)
    if not target:
        target = md_prep.sanitize_target_id(working_dir.name)
        if target.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_npt.gro")) or list(working_dir.glob("*_topol.top"))
            if candidates:
                target = candidates[0].stem.replace("_npt", "").replace("_topol", "")
    target = md_prep.sanitize_target_id(target)

    console.print(
        Panel.fit(
            f"[bold blue]Produção de Dinâmica Molecular (GROMACS)[/bold blue]\n"
            f"Alvo / Target ID: [green]{target}[/green]\n"
            f"Diretório de Trabalho: {working_dir}",
            border_style="blue",
        )
    )

    start_time = time.time()
    try:
        console.print(
            f"[yellow]Iniciando a produção da Dinâmica Molecular para {target} (grompp & mdrun)...[/yellow]"
        )
        md_analysis.run_production_md(working_dir, target_id=target)
        duration = time.time() - start_time
        console.print(f"[bold green]✓ Produção concluída com sucesso para {target}![/bold green]")
        console.print(
            f"[cyan]Execute 'md-postprocess --dir {working_dir}' para realizar o tratamento de PBC, gráficos e MM-PBSA.[/cyan]"
        )

        notifier.send_email_alert(
            step_name=f"Produção de Dinâmica Molecular ({target})",
            status="success",
            duration_seconds=duration,
            details={
                "Alvo": target,
                "Diretório": str(working_dir),
                "Trajetória Gerada": f"{target}_md.xtc e {target}_md.gro",
                "Próxima Ação Recomendada": f"Executar md-postprocess --dir {working_dir}",
            },
            console_logger=console,
        )

    except md_prep.DependencyError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name=f"Produção de Dinâmica Molecular ({target})",
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
            step_name=f"Produção de Dinâmica Molecular ({target})",
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
            step_name=f"Produção de Dinâmica Molecular ({target})",
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
            step_name=f"Produção de Dinâmica Molecular ({target})",
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
    target: str = typer.Option(
        None, "--target", help="Identificador único do alvo (ex: 7CFN, 1OSV, 4HG7)"
    ),
    skip_mmpbsa: bool = typer.Option(
        False, "--skip-mmpbsa", help="Pular o cálculo de energia livre MM-PBSA"
    ),
):
    """
    PÓS-PROCESSAMENTO, GRÁFICOS E MM-PBSA DA DINÂMICA MOLECULAR:
    1. Tratamento de Condições Periódicas de Contorno (PBC: remoção de saltos e centralização via <target_id>_md_fit.xtc e <target_id>_md_clean.gro).
    2. Análises Estruturais Globais (0 - 100 ns: RMSD Backbone + Ligante, RMSF C-α, HBond, Rg, SASA, Clustering & CSVs).
    3. Geração automatizada de gráficos científicos (.png a 300 DPI) para publicação cobrindo a janela completa de 0 - 100 ns.
    4. Cálculo de energia livre de ligação MM-PBSA (gmx_MMPBSA) na Janela Termodinâmica de estado estacionário (60 - 100 ns / Últimos 40%).
    """
    working_dir = Path(working_dir)
    if not target:
        target = md_prep.sanitize_target_id(working_dir.name)
        if target.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_md.tpr")) or list(working_dir.glob("*_md.xtc"))
            if candidates:
                target = candidates[0].stem.replace("_md", "")
    target = md_prep.sanitize_target_id(target)

    console.print(
        Panel.fit(
            f"[bold blue]Pós-processamento, Gráficos e MM-PBSA da Dinâmica Molecular[/bold blue]\n"
            f"Alvo / Target ID: [green]{target}[/green]\n"
            f"Diretório de Trabalho: {working_dir}\n"
            f"[dim]Protocolo: Dupla Escala Temporal (Estrutural: 0 - 100 ns | Termodinâmica: 60 - 100 ns)[/dim]",
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
                description=f"[1/4] Tratamento Automatizado de Trajetória ({target}_md_fit.xtc, {target}_md_clean.gro)...",
                total=1,
            )
            md_analysis.fix_pbc(working_dir, target_id=target)
            progress.update(task_pbc, completed=1)

            task_traj = progress.add_task(
                description=f"[2/4] Análises Estruturais Globais (0 - 100 ns: RMSD, RMSF, HBond, Rg, SASA, Clustering & CSVs)...",
                total=1,
            )
            md_analysis.analyze_trajectory(working_dir, target_id=target)
            progress.update(task_traj, completed=1)

            task_plot = progress.add_task(
                description="[3/4] Geração de Gráficos Científicos de Publicação (Janela Completa: 0 - 100 ns, 300 DPI)...",
                total=1,
            )
            plots = md_analysis.plot_md_results(working_dir, target_id=target)
            progress.update(task_plot, completed=1)

            mmpbsa_res = None
            if not skip_mmpbsa:
                task_mmpbsa = progress.add_task(
                    description="[4/4] Cálculo de Energia Livre MM-PBSA (Janela: 60 - 100 ns / Últimos 40% - Estado Estacionário)...",
                    total=1,
                )
                mmpbsa_res = md_analysis.calculate_mmpbsa(working_dir, target_id=target)
                progress.update(task_mmpbsa, completed=1)

        duration = time.time() - start_time
        console.print(
            f"\n[bold green]✓ Pós-processamento concluído com sucesso para {target}![/bold green]"
        )

        if plots:
            console.print(
                "\n[bold cyan]Gráficos Científicos Gerados (Janela: 0 - 100 ns):[/bold cyan]"
            )
            for p_name, p_path in plots.items():
                console.print(
                    f"  • [bold]{p_name.upper()}:[/bold] [green]{p_path}[/green]"
                )

        details = {
            "Alvo": target,
            "Diretório": str(working_dir),
            "Tratamento PBC": f"Concluído ({target}_md_fit.xtc e {target}_md_clean.gro gerados)",
            "Análise Estrutural": "Janela Completa 0 - 100 ns (RMSD Backbone + Ligante, RMSF, HBond, Rg, SASA)",
            "Gráficos Gerados": f"{len(plots)} gráficos salvos (300 DPI)"
            if plots
            else "Nenhum",
        }

        if mmpbsa_res and "energies" in mmpbsa_res:
            table = Table(
                title=f"Resumo Termodinâmico MM-PBSA ({target}) [Janela: 60 - 100 ns (Últimos 40%)]",
                show_header=True,
                header_style="bold magenta",
            )
            table.add_column("Componente Energético", style="dim", width=34)
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
                f"Sumário JSON salvo em: [cyan]{working_dir / f'{target}_mmpbsa_summary.json'}[/cyan]"
            )

            dg_bind = f"{energies['delta_g_binding']['mean']:.2f} ± {energies['delta_g_binding']['std']:.2f} {mmpbsa_res.get('unit', 'kcal/mol')}"
            details["MM-PBSA Janela"] = "60 - 100 ns (Últimos 40% - Estado Estacionário)"
            details["MM-PBSA ΔG Total (Ligação)"] = dg_bind
            details["Sumário JSON"] = str(working_dir / f"{target}_mmpbsa_summary.json")

        # Notificação por E-mail
        notifier.send_email_alert(
            step_name=f"Pós-processamento e MM-PBSA da DM ({target})",
            status="success",
            duration_seconds=duration,
            details=details,
            console_logger=console,
        )

    except md_prep.DependencyError as e:
        duration = time.time() - start_time
        notifier.send_email_alert(
            step_name=f"Pós-processamento e MM-PBSA da DM ({target})",
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
            step_name=f"Pós-processamento e MM-PBSA da DM ({target})",
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
            step_name=f"Pós-processamento e MM-PBSA da DM ({target})",
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
            step_name=f"Pós-processamento e MM-PBSA da DM ({target})",
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
