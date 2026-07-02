#!/usr/bin/env python

import json
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
    preparation,
    vina_runner,
    pharmacokinetics,
    md_prep,
    md_equil,
    md_analysis,
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
    pgp_status = (
        "[bold yellow]Efluxo Ativo[/bold yellow]"
        if "Substrato" in pgp
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

    # Veredito Geral Integrado
    pass_filters = admet.get("pass_filters", False)
    lipinski_pass = admet.get("lipinski_pass", False)
    veber_pass = admet.get("veber_pass", False)
    hia_status_val = admet.get("hia_status", "")

    if pass_filters:
        veredito = "[bold white on green]  APROVADO (BIODISPONÍVEL & SEGURO)  [/bold white on green]"
        message = (
            f"Veredito de Triagem ADMET: {veredito}\n"
            f"• Físico-Química: A molécula atende às regras clássicas de Lipinski e Veber.\n"
            f"• Farmacocinética: Alta Absorção Intestinal (HIA) estimada.\n"
            f"• Toxicidade: Nenhum alerta estrutural reativo ou PAINS foi identificado."
        )
        border_style = "green"
    else:
        reasons = []
        if not lipinski_pass:
            reasons.append("Violou regras de Lipinski")
        if not veber_pass:
            reasons.append("Violou regras de Veber")
        if hia_status_val == "Baixa Absorção":
            reasons.append("Baixa Absorção Intestinal (HIA)")
        if len(toxic_alerts) > 0:
            reasons.append(f"Alertas de Toxicidade/PAINS: {', '.join(toxic_alerts)}")

        veredito = "[bold white on red]  REPROVADO / RISCO ADMET  [/bold white on red]"
        reasons_str = "; ".join(reasons)
        message = (
            f"Veredito de Triagem ADMET: {veredito}\n"
            f"Problemas identificados: [red]{reasons_str}[/red]\n"
            f"Atenção: A molécula possui propriedades físico-químicas desfavoráveis, baixa absorção intestinal ou riscos de toxicidade estrutural."
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

        console.print("\n[bold green]✓ Validação concluída![/bold green]")
        console.print(f"[bold]Energia de Afinidade (Score):[/bold] {score} kcal/mol")

        if rmsd is not None:
            color = "green" if rmsd <= 2.0 else "red"
            veredito = "SUCESSO" if rmsd <= 2.0 else "FRACASSO"
            console.print(
                f"[bold]RMSD (vs Cristal):[/bold] [{color}]{rmsd:.3f} Å[/{color}]"
            )
            console.print(f"[bold]Veredito:[/bold] [{color}]{veredito}[/{color}]")
        else:
            console.print(f"[bold red]Erro na análise:[/bold red] {error}")

        # Renderiza a tabela de contatos estáticos no terminal
        render_interactions_table(interactions)
        render_admet_table(interactions.get("pharmacokinetics", {}))  # type: ignore

    except Exception as e:
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

    except Exception as e:
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
                "7. Executar Produção e Análise da Dinâmica (100 ns)",
                "8. Sair",
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
            rec_default = "data/1OSV/processed/receptor.pdb"
            if not Path(rec_default).exists():
                rec_default = ""

            sdf_default = "data/screening/desoxicolato/docked_poses.sdf"
            if not Path(sdf_default).exists():
                sdf_default = ""

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

                console.print(
                    "\n[bold green]✓ Preparação e Minimização de Energia concluídas com sucesso![/bold green]"
                )
                console.print(
                    f"Arquivos e outputs gerados em: [cyan]{output_dir}[/cyan]"
                )

            except md_prep.DependencyError as e:
                console.print(
                    Panel(
                        f"[bold red]Erro de Dependência GROMACS/ACPYPE:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha de Dependência",
                    )
                )
            except md_prep.SimulationPrepError as e:
                console.print(
                    Panel(
                        f"[bold red]Erro de Preparação na Dinâmica Molecular:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha na Preparação",
                    )
                )
            except FileNotFoundError as e:
                console.print(
                    Panel(
                        f"[bold red]Arquivo Não Encontrado:[/bold red]\n{e}",
                        border_style="red",
                        title="Erro de Caminho",
                    )
                )
            except Exception as e:
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

                console.print(
                    "\n[bold green]✓ Equilíbrio NVT/NPT concluído com sucesso![/bold green]"
                )
                console.print(
                    f"Estruturas equilibradas geradas em: [cyan]{md_dir}[/cyan]"
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
                        f"[bold red]Erro de Execução no Equilíbrio GROMACS:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha no Equilíbrio",
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
                        f"[bold red]Erro Inesperado no Fluxo de Equilíbrio:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha Crítica",
                    )
                )

        elif choice == "7. Executar Produção e Análise da Dinâmica (100 ns)":
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
                        f"[bold blue]Produção e Análise de Dinâmica Molecular (GROMACS)[/bold blue]\n"
                        f"Diretório de Trabalho: {md_dir}",
                        border_style="blue",
                    )
                )

                console.print("[yellow]Iniciando a produção da Dinâmica Molecular (grompp & mdrun)...[/yellow]")
                md_analysis.run_production_md(Path(md_dir))
                console.print("[bold green]✓ Produção concluída com sucesso![/bold green]")
                
                console.print("[yellow]Iniciando a análise da trajetória (RMSD, RMSF, Pontes de Hidrogênio)...[/yellow]")
                md_analysis.analyze_trajectory(Path(md_dir))
                console.print("[bold green]✓ Análise da trajetória concluída com sucesso![/bold green]")
                console.print(f"Arquivos de análise gerados em: [cyan]{md_dir}[/cyan]")

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
                        f"[bold red]Erro de Execução/Análise no GROMACS:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha na Simulação/Análise",
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
                        f"[bold red]Erro Inesperado no Fluxo de Produção:[/bold red]\n{e}",
                        border_style="red",
                        title="Falha Crítica",
                    )
                )

        elif choice == "8. Sair":
            break


@app.command(name="md-run")
def md_run_command(
    working_dir: Path = typer.Option(
        ..., "--dir", help="Diretório de trabalho onde estão os arquivos do equilíbrio (nvt/npt)"
    ),
):
    """
    PRODUÇÃO E ANÁLISE DE DINÂMICA MOLECULAR:
    Compila e executa a produção da Dinâmica Molecular no GROMACS e gera as análises de trajetória (RMSD, RMSF, HBond).
    """
    console.print(
        Panel.fit(
            f"[bold blue]Produção e Análise de Dinâmica Molecular (GROMACS)[/bold blue]\n"
            f"Diretório de Trabalho: {working_dir}",
            border_style="blue",
        )
    )

    try:
        console.print("[yellow]Iniciando a produção da Dinâmica Molecular (grompp & mdrun)...[/yellow]")
        md_analysis.run_production_md(working_dir)
        console.print("[bold green]✓ Produção concluída com sucesso![/bold green]")
        
        console.print("[yellow]Iniciando a análise da trajetória (RMSD, RMSF, Pontes de Hidrogênio)...[/yellow]")
        md_analysis.analyze_trajectory(working_dir)
        console.print("[bold green]✓ Análise da trajetória concluída com sucesso![/bold green]")
        console.print(f"Arquivos de análise gerados em: [cyan]{working_dir}[/cyan]")
    except md_prep.DependencyError as e:
        console.print(f"\n[bold red]Erro de Dependência:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except md_prep.SimulationPrepError as e:
        console.print(f"\n[bold red]Erro na Dinâmica/Análise:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        console.print(f"\n[bold red]Arquivo Não Encontrado:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except Exception as e:
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

        console.print(
            "\n[bold green]✓ Preparação e Minimização de Energia concluídas com sucesso![/bold green]"
        )
        console.print(f"Arquivos e outputs gerados em: [cyan]{out}[/cyan]")
    except md_prep.DependencyError as e:
        console.print(f"\n[bold red]Erro de Dependência:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except md_prep.SimulationPrepError as e:
        console.print(f"\n[bold red]Erro de Preparação:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        console.print(f"\n[bold red]Arquivo Não Encontrado:[/bold red]\n{e}")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"\n[bold red]Erro Inesperado:[/bold red]\n{e}")
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
