import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib  # type: ignore[reportMissingImports]

matplotlib.use("Agg")  # Backend não-interativo para renderização de imagens
import matplotlib.pyplot as plt  # type: ignore[reportMissingImports]

try:
    import seaborn as sns  # type: ignore[reportMissingImports]

    sns.set_theme(style="whitegrid")
except ImportError:
    sns = None
    try:
        plt.style.use("seaborn-v0_8-whitegrid")
    except Exception:
        pass

from docking.md_prep import DependencyError, SimulationPrepError, find_executable


def get_index_groups(index_file: Path) -> List[str]:
    """
    Lê um arquivo index.ndx e retorna a lista ordenada de nomes de grupos presentes.
    """
    groups_list: List[str] = []
    if not index_file.exists():
        return groups_list
    with open(index_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                group_name = line[1:-1].strip()
                groups_list.append(group_name)
    return groups_list


def identify_complex_groups(groups_list: List[str]) -> Tuple[str, str, int, int]:
    """
    Identifica o nome e os índices numéricos dos grupos de receptor (proteína) e ligante
    a partir da lista de grupos do index.ndx.

    No GROMACS e gmx_MMPBSA, a numeração é 0-based correspondente à ordem em index.ndx
    (0 = System, 1 = Protein, ..., 13 = LIG).

    Retorna: (prot_name, lig_name, prot_idx, lig_idx)
    """
    # 1. Identificação do receptor
    if "Protein" in groups_list:
        prot_name = "Protein"
    elif "Protein-H" in groups_list:
        prot_name = "Protein-H"
    else:
        prot_name = "Protein"

    prot_idx = groups_list.index(prot_name) if prot_name in groups_list else 1

    # 2. Identificação do ligante
    lig_name = None
    for candidate in ["ligand_md", "LIG", "UNL", "UNK", "MOL", "Other"]:
        if candidate in groups_list:
            lig_name = candidate
            break

    if not lig_name:
        standard_groups = {
            "System",
            "Protein",
            "Protein-H",
            "C-alpha",
            "Backbone",
            "MainChain",
            "MainChain+Cb",
            "MainChain+H",
            "SideChain",
            "SideChain-H",
            "Prot-Masses",
            "non-Protein",
            "Other",
            "NA",
            "CL",
            "Ion",
            "Water",
            "SOL",
            "non-Water",
            "Water_and_ions",
            "Protein_LIG",
        }
        for g in groups_list:
            if g not in standard_groups:
                lig_name = g
                break

    if not lig_name:
        lig_name = "LIG" if "LIG" in groups_list else "ligand_md"

    lig_idx = groups_list.index(lig_name) if lig_name in groups_list else 13

    return prot_name, lig_name, prot_idx, lig_idx


def compile_production_tpr(working_dir: Path) -> Path:
    """
    Compila o arquivo de entrada binário para a produção da Dinâmica Molecular (md.tpr)
    usando o comando 'gmx grompp' a partir de npt.gro, npt.cpt, topol.top, index.ndx e md.mdp.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    npt_gro = working_dir / "npt.gro"
    npt_cpt = working_dir / "npt.cpt"
    topol_top = working_dir / "topol.top"
    index_ndx = working_dir / "index.ndx"

    if not npt_gro.exists():
        raise FileNotFoundError(
            f"Arquivo de equilíbrio 'npt.gro' não encontrado no diretório: {working_dir}"
        )
    if not topol_top.exists():
        raise FileNotFoundError(
            f"Arquivo de topologia 'topol.top' não encontrado no diretório: {working_dir}"
        )

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError(
            "O executável 'gmx' (GROMACS) não foi encontrado no PATH."
        )

    # Resolução do template md.mdp
    project_root = Path(__file__).resolve().parent.parent.parent
    md_mdp = project_root / "src" / "templates" / "mdp" / "md.mdp"
    if not md_mdp.exists():
        md_mdp = Path("src/templates/mdp/md.mdp").resolve()
        if not md_mdp.exists():
            raise FileNotFoundError("Arquivo template md.mdp não encontrado.")

    # 1. Compilação do arquivo de produção (grompp)
    cmd_grompp = [
        gmx_bin,
        "grompp",
        "-f",
        str(md_mdp),
        "-c",
        "npt.gro",
        "-p",
        "topol.top",
        "-o",
        "md.tpr",
    ]
    if npt_cpt.exists():
        cmd_grompp.extend(["-t", "npt.cpt"])
    if index_ndx.exists():
        cmd_grompp.extend(["-n", "index.ndx"])

    try:
        env = os.environ.copy()
        exec_dir = str(Path(gmx_bin).parent)
        env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

        result = subprocess.run(
            cmd_grompp, cwd=str(working_dir), env=env, capture_output=True, text=True
        )
        if result.returncode != 0:
            error_msg = result.stderr.strip() or result.stdout.strip()
            raise SimulationPrepError(
                f"Erro na Compilação de Produção (grompp):\n"
                f"Comando: {' '.join(cmd_grompp)}\n"
                f"Código de retorno: {result.returncode}\n"
                f"Erro real: {error_msg}"
            )
    except Exception as e:
        if isinstance(e, (SimulationPrepError, DependencyError)):
            raise e
        raise SimulationPrepError(
            f"Falha ao executar a Compilação de Produção (grompp): {e}"
        )

    md_tpr = working_dir / "md.tpr"
    if not md_tpr.exists():
        raise FileNotFoundError(f"Arquivo 'md.tpr' não foi gerado em: {working_dir}")

    return md_tpr


def run_production_md(working_dir: Path):
    """
    Compila e executa a etapa de Produção de Dinâmica Molecular no GROMACS.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError(
            "O executável 'gmx' (GROMACS) não foi encontrado no PATH."
        )

    # 1. Compilação do arquivo de produção (grompp -> md.tpr)
    compile_production_tpr(working_dir)

    # 2. Execução da simulação de produção (mdrun)
    cmd_mdrun = [gmx_bin, "mdrun", "-v", "-deffnm", "md"]

    try:
        env = os.environ.copy()
        exec_dir = str(Path(gmx_bin).parent)
        env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

        process = subprocess.Popen(
            cmd_mdrun,
            cwd=str(working_dir),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )

        # Captura e imprime o output em tempo real
        for line in process.stdout:  # type: ignore
            sys.stdout.write(line)
            sys.stdout.flush()

        process.wait()
        if process.returncode != 0:
            raise SimulationPrepError(
                f"Erro na Execução de Produção (mdrun):\n"
                f"Comando: {' '.join(cmd_mdrun)}\n"
                f"Código de retorno: {process.returncode}"
            )
    except Exception as e:
        if isinstance(e, (SimulationPrepError, DependencyError)):
            raise e
        raise SimulationPrepError(
            f"Falha ao executar a Produção de Dinâmica Molecular (mdrun) em tempo real: {e}"
        )


def fix_pbc(working_dir: Path, force: bool = False) -> Path:
    """
    Executa o tratamento completo e automatizado de Condições Periódicas de Contorno (PBC)
    e remoção de movimentos de corpo rígido (fit rot+trans) via GROMACS (gmx trjconv):

    1. Etapa 1 - Centralização e Correção de PBC:
       gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -pbc mol -center
       (Input stdin: '1\\n0\\n' -> Protein para centralizar, System para salvar)

    2. Etapa 2 - Ajuste de Rotação e Translação (Fit de Mínimos Quadrados):
       gmx trjconv -s md.tpr -f md_center.xtc -o md_fit.xtc -fit rot+trans
       (Input stdin: '4\\n0\\n' -> Backbone para fit de RMSD, System para salvar)

    3. Etapa 3 - Estrutura Estática Limpa para PyMOL:
       gmx trjconv -s md.tpr -f md.gro -o md_clean.gro -pbc mol -center
       (Input stdin: '1\\n0\\n' -> Protein para centralizar, System para salvar)

    Se 'md_fit.xtc' e 'md_clean.gro' já existirem no diretório e force=False, reutiliza os arquivos existentes.
    Retorna o caminho do arquivo de trajetória corrigida md_fit.xtc.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    md_fit_path = working_dir / "md_fit.xtc"
    md_clean_path = working_dir / "md_clean.gro"

    if not force and md_fit_path.exists() and md_clean_path.exists():
        return md_fit_path

    tpr_file = working_dir / "md.tpr"
    xtc_file = working_dir / "md.xtc"

    if not tpr_file.exists():
        raise FileNotFoundError(
            f"Arquivo de topologia compilada 'md.tpr' não encontrado em: {working_dir}"
        )
    if not xtc_file.exists():
        raise FileNotFoundError(
            f"Arquivo de trajetória de produção 'md.xtc' não encontrado em: {working_dir}"
        )

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError(
            "O executável 'gmx' (GROMACS) não foi encontrado no PATH."
        )

    def run_trjconv_cmd(cmd: List[str], input_val: str, step_name: str):
        try:
            env = os.environ.copy()
            exec_dir = str(Path(gmx_bin).parent)
            env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

            result = subprocess.run(
                cmd,
                cwd=str(working_dir),
                env=env,
                capture_output=True,
                text=True,
                input=input_val,
            )
            if result.returncode != 0:
                error_msg = result.stderr.strip() or result.stdout.strip()
                raise SimulationPrepError(
                    f"Erro no tratamento de PBC ({step_name}):\n"
                    f"Comando: {' '.join(cmd)}\n"
                    f"Código de retorno: {result.returncode}\n"
                    f"Erro real: {error_msg}"
                )
            return result
        except Exception as e:
            if isinstance(e, (SimulationPrepError, DependencyError, FileNotFoundError)):
                raise e
            raise SimulationPrepError(
                f"Falha ao executar o tratamento de PBC ({step_name}): {e}"
            )

    # Etapa 1: Centralização da Proteína e Correção de Moléculas Quebradas em PBC
    # Seleção de grupos: 1 (Protein para centralizar) e 0 (System para salvar)
    cmd_center = [
        gmx_bin,
        "trjconv",
        "-s",
        "md.tpr",
        "-f",
        "md.xtc",
        "-o",
        "md_center.xtc",
        "-pbc",
        "mol",
        "-center",
    ]
    run_trjconv_cmd(
        cmd_center,
        input_val="1\n0\n",
        step_name="Etapa 1 (Centralização e Correção de PBC - pbc mol -center)",
    )

    # Etapa 2: Ajuste de Rotação e Translação (Fit de Mínimos Quadrados via Backbone)
    # Seleção de grupos: 4 (Backbone para fit de RMSD) e 0 (System para salvar)
    cmd_fit = [
        gmx_bin,
        "trjconv",
        "-s",
        "md.tpr",
        "-f",
        "md_center.xtc",
        "-o",
        "md_fit.xtc",
        "-fit",
        "rot+trans",
    ]
    run_trjconv_cmd(
        cmd_fit,
        input_val="4\n0\n",
        step_name="Etapa 2 (Ajuste de Rotação e Translação - fit rot+trans)",
    )

    # Etapa 3: Estrutura Estática Limpa para Visualização no PyMOL
    # Seleção de grupos: 1 (Protein para centralizar) e 0 (System para salvar)
    gro_candidates = ["md.gro", "npt.gro", "em.gro", "complex.gro"]
    gro_source = None
    for cand in gro_candidates:
        if (working_dir / cand).exists():
            gro_source = cand
            break

    if gro_source:
        cmd_clean_gro = [
            gmx_bin,
            "trjconv",
            "-s",
            "md.tpr",
            "-f",
            gro_source,
            "-o",
            "md_clean.gro",
            "-pbc",
            "mol",
            "-center",
        ]
        run_trjconv_cmd(
            cmd_clean_gro,
            input_val="1\n0\n",
            step_name="Etapa 3 (Estrutura Estática Limpa para PyMOL - md_clean.gro)",
        )

    if not md_fit_path.exists():
        raise FileNotFoundError(
            f"Arquivo de trajetória corrigida 'md_fit.xtc' não foi gerado em: {working_dir}"
        )
    if not md_clean_path.exists():
        raise FileNotFoundError(
            f"Arquivo de estrutura estática limpa 'md_clean.gro' não foi gerado em: {working_dir}"
        )

    return md_fit_path


def analyze_trajectory(working_dir: Path):
    """
    Executa a análise quantitativa da trajetória de Dinâmica Molecular no GROMACS (Janela Completa: 0 - 100 ns):
    1. RMSD do Backbone da Proteína e do Ligante (gmx rms -s md.tpr -f md_fit.xtc -o rmsd.xvg -tu ns)
    2. RMSF por resíduo dos Carbonos Alfa (gmx rmsf -s md.tpr -f md_fit.xtc -o rmsf.xvg -res)
    3. Monitoramento de Pontes de Hidrogênio Proteína-Ligante (gmx hbond -s md.tpr -f md_fit.xtc -num hbond.xvg -tu ns)

    Utiliza estritamente a trajetória ajustada e sem artefatos de PBC (md_fit.xtc).
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    tpr_file = working_dir / "md.tpr"
    fit_xtc = working_dir / "md_fit.xtc"

    if not tpr_file.exists():
        raise FileNotFoundError(
            f"Arquivo 'md.tpr' não encontrado no diretório: {working_dir}"
        )
    if not fit_xtc.exists():
        raise FileNotFoundError(
            f"Trajetória corrigida 'md_fit.xtc' não encontrada em: {working_dir}. "
            "Certifique-se de executar o tratamento de PBC (fix_pbc) antes da análise."
        )

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError(
            "O executável 'gmx' (GROMACS) não foi encontrado no PATH."
        )

    def run_analysis_cmd(cmd: List[str], cwd: Path, input_val: str, step_name: str = ""):
        try:
            env = os.environ.copy()
            exec_dir = str(Path(gmx_bin).parent)
            env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

            result = subprocess.run(
                cmd,
                cwd=str(cwd),
                env=env,
                capture_output=True,
                text=True,
                input=input_val,
            )
            if result.returncode != 0:
                error_msg = result.stderr.strip() or result.stdout.strip()
                raise SimulationPrepError(
                    f"Erro na análise ({step_name}):\n"
                    f"Comando: {' '.join(cmd)}\n"
                    f"Código de retorno: {result.returncode}\n"
                    f"Erro real: {error_msg}"
                )
            return result
        except Exception as e:
            if isinstance(e, (SimulationPrepError, DependencyError)):
                raise e
            raise SimulationPrepError(
                f"Falha ao executar comando de análise ({step_name}): {e}"
            )

    index_file = working_dir / "index.ndx"
    prot_name = "Protein"
    lig_name = "LIG"
    prot_idx = 1
    lig_idx = 13
    if index_file.exists():
        groups_list = get_index_groups(index_file)
        prot_name, lig_name, prot_idx, lig_idx = identify_complex_groups(groups_list)

    # 1.1 RMSD do esqueleto da proteína (Backbone/Backbone -> 4 4)
    cmd_rmsd = [
        gmx_bin,
        "rms",
        "-s",
        "md.tpr",
        "-f",
        "md_fit.xtc",
        "-o",
        "rmsd.xvg",
        "-tu",
        "ns",
    ]
    run_analysis_cmd(cmd_rmsd, working_dir, "4\n4\n", "RMSD do esqueleto da proteína (Backbone)")

    # 1.2 RMSD do Ligante no sítio ativo (Fit: Backbone 4, Calc: Ligand) para monitorar persistência e ausência de unbinding
    if index_file.exists() and lig_idx is not None:
        try:
            cmd_rmsd_lig = [
                gmx_bin,
                "rms",
                "-s",
                "md.tpr",
                "-f",
                "md_fit.xtc",
                "-o",
                "rmsd_lig.xvg",
                "-tu",
                "ns",
                "-n",
                "index.ndx",
            ]
            run_analysis_cmd(
                cmd_rmsd_lig,
                working_dir,
                f"4\n{lig_idx}\n",
                f"RMSD do ligante ({lig_name}) no sítio ativo",
            )
        except Exception:
            pass

    # 2. RMSF por resíduo dos Carbonos Alfa (C-alpha -> 3)
    cmd_rmsf = [
        gmx_bin,
        "rmsf",
        "-s",
        "md.tpr",
        "-f",
        "md_fit.xtc",
        "-o",
        "rmsf.xvg",
        "-res",
    ]
    run_analysis_cmd(cmd_rmsf, working_dir, "3\n", "RMSF por resíduo (C-alpha)")

    # 3. Pontes de hidrogênio entre Proteína e Ligante ao longo dos 100 ns
    cmd_hbond = [
        gmx_bin,
        "hbond",
        "-s",
        "md.tpr",
        "-f",
        "md_fit.xtc",
        "-num",
        "hbond.xvg",
        "-tu",
        "ns",
    ]
    if index_file.exists():
        cmd_hbond.extend(["-n", "index.ndx"])
        hbond_input = f"{prot_name}\n{lig_name}\n"
    else:
        hbond_input = "Protein\nLIG\n"

    run_analysis_cmd(
        cmd_hbond,
        working_dir,
        hbond_input,
        "Pontes de hidrogênio (Proteína-Ligante)",
    )


def parse_xvg_with_meta(
    file_path: Path,
) -> Tuple[List[float], List[float], Dict[str, str]]:
    """
    Faz o parse de arquivos .xvg gerados pelo GROMACS, extraindo dados e metadados (@ xaxis label, etc.).
    Retorna (x_vals, y_vals, metadata_dict).
    """
    x_vals: List[float] = []
    y_vals: List[float] = []
    meta: Dict[str, str] = {}

    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo XVG não encontrado: {file_path}")

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("@"):
                line_clean = line[1:].strip()
                if "xaxis" in line_clean and "label" in line_clean:
                    parts = line_clean.split("label", 1)
                    if len(parts) > 1:
                        meta["xaxis_label"] = parts[1].strip().strip('"')
                elif "yaxis" in line_clean and "label" in line_clean:
                    parts = line_clean.split("label", 1)
                    if len(parts) > 1:
                        meta["yaxis_label"] = parts[1].strip().strip('"')
                elif "title" in line_clean:
                    parts = line_clean.split("title", 1)
                    if len(parts) > 1:
                        meta["title"] = parts[1].strip().strip('"')
                continue
            if line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) >= 2:
                try:
                    x_vals.append(float(parts[0]))
                    y_vals.append(float(parts[1]))
                except ValueError:
                    continue

    return x_vals, y_vals, meta


def parse_xvg(file_path: Path) -> Tuple[List[float], List[float]]:
    """
    Faz o parse de arquivos .xvg gerados pelo GROMACS, ignorando metadados iniciados com '@' e '#'.
    Retorna uma tupla contendo duas listas de floats (valores do eixo X e eixo Y).
    """
    x_vals, y_vals, _ = parse_xvg_with_meta(file_path)
    return x_vals, y_vals


def plot_md_results(working_dir: Path) -> Dict[str, Path]:
    """
    Lê os arquivos .xvg gerados na análise (rmsd.xvg, rmsd_lig.xvg, rmsf.xvg, hbond.xvg) e gera gráficos
    com padrão estético científico de publicação (300 DPI) cobrindo a extensão completa da simulação (0 - 100 ns).

    Salva diretamente no working_dir:
    - rmsd.png: Tempo (ns) vs RMSD (nm) para Backbone da Proteína e Ligante no sítio ativo (0 - 100 ns)
    - rmsf.png: Número do Resíduo vs Flutuação RMSF (nm) (0 - 100 ns)
    - hbond.png: Tempo (ns) vs Número de Pontes de Hidrogênio (0 - 100 ns)

    Retorna um dicionário mapeando o nome da análise ao caminho do arquivo .png gerado.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    # Configuração de estilo científico de alta qualidade
    if sns is not None:
        sns.set_theme(style="whitegrid")  # type: ignore
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "figure.titlesize": 14,
            "figure.autolayout": True,
        }
    )

    generated_plots: Dict[str, Path] = {}

    # 1. Gráfico de RMSD (Janela Completa: 0 - 100 ns - Backbone & Ligante)
    rmsd_file = working_dir / "rmsd.xvg"
    rmsd_lig_file = working_dir / "rmsd_lig.xvg"

    if rmsd_file.exists():
        x_time, y_rmsd, meta_rmsd = parse_xvg_with_meta(rmsd_file)
        if x_time and y_rmsd:
            sorted_rmsd = sorted(zip(x_time, y_rmsd), key=lambda p: p[0])
            x_time = [p[0] for p in sorted_rmsd]
            y_rmsd = [p[1] for p in sorted_rmsd]

            # Conversão automática de unidades (ps para ns) se necessário
            xaxis_lbl = meta_rmsd.get("xaxis_label", "").lower()
            if (
                "ps" in xaxis_lbl
                or "(ps)" in xaxis_lbl
                or (max(x_time) > 1000 and "ns" not in xaxis_lbl)
            ):
                x_time = [t / 1000.0 for t in x_time]

            fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=300)
            ax.plot(
                x_time, y_rmsd, color="#1f77b4", linewidth=1.6, label="Protein Backbone"
            )

            # Adiciona linha de média da proteína
            mean_rmsd = sum(y_rmsd) / len(y_rmsd)
            ax.axhline(
                mean_rmsd,
                color="#1f77b4",
                linestyle="--",
                alpha=0.6,
                label=f"Protein Mean: {mean_rmsd:.3f} nm",
            )

            max_x = max(x_time) if x_time else 100.0

            # Plota RMSD do ligante se disponível para demonstrar persistência no sítio
            if rmsd_lig_file.exists():
                x_lig, y_lig, meta_lig = parse_xvg_with_meta(rmsd_lig_file)
                if x_lig and y_lig:
                    sorted_lig = sorted(zip(x_lig, y_lig), key=lambda p: p[0])
                    x_lig = [p[0] for p in sorted_lig]
                    y_lig = [p[1] for p in sorted_lig]
                    xaxis_lbl_lig = meta_lig.get("xaxis_label", "").lower()
                    if (
                        "ps" in xaxis_lbl_lig
                        or "(ps)" in xaxis_lbl_lig
                        or (max(x_lig) > 1000 and "ns" not in xaxis_lbl_lig)
                    ):
                        x_lig = [t / 1000.0 for t in x_lig]
                    ax.plot(
                        x_lig,
                        y_lig,
                        color="#e76f51",
                        linewidth=1.4,
                        alpha=0.85,
                        label="Ligand (Site Persistence)",
                    )
                    mean_lig = sum(y_lig) / len(y_lig)
                    ax.axhline(
                        mean_lig,
                        color="#e76f51",
                        linestyle=":",
                        alpha=0.7,
                        label=f"Ligand Mean: {mean_lig:.3f} nm",
                    )
                    max_x = max(max_x, max(x_lig))

            ax.set_xlabel("Time (ns)", fontweight="bold")
            ax.set_ylabel("RMSD (nm)", fontweight="bold")
            ax.set_title(
                "Structural Stability & Ligand Residence (0 - 100 ns)",
                fontweight="bold",
                pad=12,
            )
            ax.set_xlim(left=0, right=max_x)
            ax.set_ylim(bottom=0)
            ax.legend(loc="upper left", frameon=True, framealpha=0.9)

            out_rmsd_png = working_dir / "rmsd.png"
            fig.savefig(out_rmsd_png, dpi=300, bbox_inches="tight")  # type: ignore
            plt.close(fig)
            generated_plots["rmsd"] = out_rmsd_png

    # 2. Gráfico de RMSF (Janela Completa: 0 - 100 ns)
    rmsf_file = working_dir / "rmsf.xvg"
    if rmsf_file.exists():
        x_res, y_rmsf, _ = parse_xvg_with_meta(rmsf_file)
        if x_res and y_rmsf:
            fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=300)

            # Tratamento para ordenação e suporte a múltiplos segmentos/cadeias
            if len(set(x_res)) == len(x_res):
                sorted_data = sorted(zip(x_res, y_rmsf), key=lambda p: p[0])
                x_plot = [p[0] for p in sorted_data]
                y_plot = [p[1] for p in sorted_data]
                chain_boundaries = []
                x_label = "Residue Number"
            else:
                segments: List[List[Tuple[float, float]]] = []
                current_segment: List[Tuple[float, float]] = []

                for res, val in zip(x_res, y_rmsf):
                    if current_segment and res <= current_segment[-1][0]:
                        segments.append(current_segment)
                        current_segment = []
                    current_segment.append((res, val))
                if current_segment:
                    segments.append(current_segment)

                x_plot = []
                y_plot = []
                chain_boundaries = []
                accum_offset = 0

                for seg_idx, seg in enumerate(segments):
                    sorted_seg = sorted(seg, key=lambda p: p[0])
                    max_res_in_seg = int(sorted_seg[-1][0])
                    for res, val in sorted_seg:
                        x_plot.append(res + accum_offset)
                        y_plot.append(val)
                    accum_offset += max_res_in_seg
                    if seg_idx < len(segments) - 1:
                        chain_boundaries.append(accum_offset)

                x_label = "Residue Index (Continuous)"

            ax.plot(x_plot, y_plot, color="#2a9d8f", linewidth=1.4, label="C-α RMSF")
            ax.fill_between(x_plot, y_plot, color="#2a9d8f", alpha=0.25)  # type: ignore

            for cb in chain_boundaries:
                ax.axvline(
                    x=cb + 0.5,
                    color="#6c757d",
                    linestyle=":",
                    linewidth=1.0,
                    alpha=0.6,
                )

            ax.set_xlabel(x_label, fontweight="bold")
            ax.set_ylabel("RMSF (nm)", fontweight="bold")
            ax.set_title(
                "Root Mean Square Fluctuation per Residue (0 - 100 ns)", fontweight="bold", pad=12
            )
            ax.set_xlim(
                left=min(x_plot) if x_plot else 0, right=max(x_plot) if x_plot else 1
            )
            ax.set_ylim(bottom=0)
            ax.legend(loc="upper right", frameon=True, framealpha=0.9)

            out_rmsf_png = working_dir / "rmsf.png"
            fig.savefig(out_rmsf_png, dpi=300, bbox_inches="tight")  # type: ignore
            plt.close(fig)
            generated_plots["rmsf"] = out_rmsf_png

    # 3. Gráfico de Pontes de Hidrogênio (HBond) (Janela Completa: 0 - 100 ns)
    hbond_file = working_dir / "hbond.xvg"
    if hbond_file.exists():
        x_time, y_hb, meta_hb = parse_xvg_with_meta(hbond_file)
        if x_time and y_hb:
            sorted_hb = sorted(zip(x_time, y_hb), key=lambda p: p[0])
            x_time = [p[0] for p in sorted_hb]
            y_hb = [p[1] for p in sorted_hb]

            xaxis_lbl = meta_hb.get("xaxis_label", "").lower()
            if (
                "ps" in xaxis_lbl
                or "(ps)" in xaxis_lbl
                or (max(x_time) > 1000 and "ns" not in xaxis_lbl)
            ):
                x_time = [t / 1000.0 for t in x_time]

            fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=300)
            ax.plot(
                x_time,
                y_hb,
                color="#e76f51",
                linewidth=1.0,
                alpha=0.75,
                label="H-Bonds count",
            )

            if len(y_hb) >= 10:
                window_size = max(5, len(y_hb) // 25)
                smoothed = []
                for i in range(len(y_hb)):
                    start_i = max(0, i - window_size // 2)
                    end_i = min(len(y_hb), i + window_size // 2 + 1)
                    smoothed.append(sum(y_hb[start_i:end_i]) / (end_i - start_i))
                ax.plot(
                    x_time,
                    smoothed,
                    color="#9d0208",
                    linewidth=1.8,
                    label=f"Running Average (window {window_size})",
                )

            ax.set_xlabel("Time (ns)", fontweight="bold")
            ax.set_ylabel("H-Bonds count", fontweight="bold")
            ax.set_title("Protein–Ligand Hydrogen Bonds (0 - 100 ns)", fontweight="bold", pad=12)
            ax.set_xlim(left=0, right=max(x_time) if x_time else 100.0)
            ax.set_ylim(bottom=0)
            ax.legend(loc="upper right", frameon=True, framealpha=0.9)

            out_hbond_png = working_dir / "hbond.png"
            fig.savefig(out_hbond_png, dpi=300, bbox_inches="tight")  # type: ignore
            plt.close(fig)
            generated_plots["hbond"] = out_hbond_png

    if not generated_plots:
        raise FileNotFoundError(
            f"Nenhum arquivo de análise (.xvg) foi encontrado em {working_dir} para geração dos gráficos."
        )

    return generated_plots


def parse_mmpbsa_dat(dat_path: Path) -> Dict[str, Any]:
    """
    Parse estruturado do arquivo de resultados FINAL_RESULTS_MMPBSA.dat gerado pelo gmx_MMPBSA.
    Extrai contribuições de Van der Waals, Eletrostática, Solvatação Polar e Apolar, e o Delta G total
    tanto para o modelo Generalized Born (MM-GBSA) quanto Poisson Boltzmann (MM-PBSA).
    """
    summary: Dict[str, Any] = {
        "unit": "kcal/mol",
        "energies": {
            "van_der_waals": {"mean": 0.0, "std": 0.0},
            "electrostatic": {"mean": 0.0, "std": 0.0},
            "polar_solvation": {"mean": 0.0, "std": 0.0},
            "nonpolar_solvation": {"mean": 0.0, "std": 0.0},
            "delta_g_binding": {"mean": 0.0, "std": 0.0},
        },
        "models": {
            "gb": {
                "van_der_waals": {"mean": 0.0, "std": 0.0},
                "electrostatic": {"mean": 0.0, "std": 0.0},
                "polar_solvation": {"mean": 0.0, "std": 0.0},
                "nonpolar_solvation": {"mean": 0.0, "std": 0.0},
                "gas_energy": {"mean": 0.0, "std": 0.0},
                "solvation_energy": {"mean": 0.0, "std": 0.0},
                "delta_g_binding": {"mean": 0.0, "std": 0.0},
            },
            "pb": {
                "van_der_waals": {"mean": 0.0, "std": 0.0},
                "electrostatic": {"mean": 0.0, "std": 0.0},
                "polar_solvation": {"mean": 0.0, "std": 0.0},
                "nonpolar_solvation": {"mean": 0.0, "std": 0.0},
                "gas_energy": {"mean": 0.0, "std": 0.0},
                "solvation_energy": {"mean": 0.0, "std": 0.0},
                "delta_g_binding": {"mean": 0.0, "std": 0.0},
            },
        },
    }

    if not dat_path.exists():
        return summary

    with open(dat_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    current_model = None
    in_delta_section = False

    for line in lines:
        line_clean = line.strip()
        if not line_clean:
            continue

        if "GENERALIZED BORN:" in line_clean:
            current_model = "gb"
            in_delta_section = False
            continue
        elif "POISSON BOLTZMANN:" in line_clean:
            current_model = "pb"
            in_delta_section = False
            continue

        if "Delta (Complex - Receptor - Ligand):" in line_clean:
            in_delta_section = True
            continue

        if in_delta_section and current_model:
            cleaned_line = line_clean.replace("Δ", "").replace("Delta", "").strip()
            parts = cleaned_line.split()
            if len(parts) >= 3:
                comp_name = parts[0].upper()
                try:
                    mean_val = float(parts[1])
                    std_val = float(parts[3]) if len(parts) >= 4 else float(parts[2])
                except ValueError:
                    continue

                target_dict = summary["models"][current_model]
                if comp_name == "VDWAALS":
                    target_dict["van_der_waals"] = {"mean": mean_val, "std": std_val}
                elif comp_name == "EEL":
                    target_dict["electrostatic"] = {"mean": mean_val, "std": std_val}
                elif comp_name in ("EGB", "EPB"):
                    target_dict["polar_solvation"] = {"mean": mean_val, "std": std_val}
                elif comp_name in ("ESURF", "ENPOLAR"):
                    target_dict["nonpolar_solvation"] = {
                        "mean": mean_val,
                        "std": std_val,
                    }
                elif comp_name == "GGAS":
                    target_dict["gas_energy"] = {"mean": mean_val, "std": std_val}
                elif comp_name == "GSOLV":
                    target_dict["solvation_energy"] = {"mean": mean_val, "std": std_val}
                elif comp_name == "TOTAL":
                    target_dict["delta_g_binding"] = {"mean": mean_val, "std": std_val}

    gb = summary["models"]["gb"]
    pb = summary["models"]["pb"]
    primary = gb if gb["delta_g_binding"]["mean"] != 0.0 else pb
    summary["energies"] = {
        "van_der_waals": primary["van_der_waals"],
        "electrostatic": primary["electrostatic"],
        "polar_solvation": primary["polar_solvation"],
        "nonpolar_solvation": primary["nonpolar_solvation"],
        "delta_g_binding": primary["delta_g_binding"],
    }

    return summary


def _ensure_gmx_mmpbsa_cys_patched() -> None:
    """
    Verifica e aplica auto-patch no arquivo make_top.py do pacote gmx_MMPBSA caso detecte
    o bug de indexação de pontes dissulfeto em sistemas proteicos multicadeia.
    """
    try:
        import GMXMMPBSA.make_top as mt

        make_top_path = Path(mt.__file__)
        if make_top_path.exists():
            with open(make_top_path, "r", encoding="utf-8") as f:
                code = f.read()

            target = (
                "if str_name == 'COM':\n"
                "                                        cys1 = c\n"
                "                                        cys2 = structure.residues.index(bondedatm.residue) + 1\n"
                "                                    else:\n"
                "                                        cys1 = residue.number\n"
                "                                        cys2 = bondedatm.residue.number"
            )
            replacement = (
                "if str_name:\n"
                "                                    cys1 = c\n"
                "                                    cys2 = structure.residues.index(bondedatm.residue) + 1"
            )

            if target in code:
                code = code.replace(target, replacement)
                with open(make_top_path, "w", encoding="utf-8") as f:
                    f.write(code)
    except Exception:
        pass


def calculate_mmpbsa(working_dir: Path) -> Dict[str, Any]:
    """
    Executa o cálculo de Energia Livre de Ligação MM-PBSA via gmx_MMPBSA (Janela Termodinâmica: 60 - 100 ns / Últimos 40%):
    1. Extrai o número total de frames e configura dinamicamente o arquivo 'mmpbsa.in' para os últimos 40% (estado estacionário).
       Ex: para 1000 frames totais, define startframe=600, endframe=1000, interval=2.
    2. Identifica os grupos do receptor (Protein) e ligante (ligand_md / LIG) em index.ndx.
    3. Executa o gmx_MMPBSA via subprocesso não-interativo com captura de stderr.
    4. Extrai contribuições energéticas (Van der Waals, Eletrostática, Solvatação Polar e Apolar) e Delta G total.
    5. Salva os resultados estruturados no arquivo mmpbsa_summary.json.

    Retorna o dicionário com o sumário dos resultados termodinâmicos.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    tpr_file = working_dir / "md.tpr"
    fit_xtc = working_dir / "md_fit.xtc"
    index_file = working_dir / "index.ndx"

    if not tpr_file.exists():
        raise FileNotFoundError(f"Arquivo 'md.tpr' não encontrado em: {working_dir}")
    if not fit_xtc.exists():
        raise FileNotFoundError(
            f"Trajetória corrigida 'md_fit.xtc' não encontrada em: {working_dir}. "
            "Execute fix_pbc antes de calcular o MM-PBSA."
        )
    if not index_file.exists():
        raise FileNotFoundError(f"Arquivo 'index.ndx' não encontrado em: {working_dir}")

    mmpbsa_bin = find_executable("gmx_MMPBSA")
    if not mmpbsa_bin:
        raise DependencyError("O executável 'gmx_MMPBSA' não foi encontrado no PATH.")

    # Garante que o patch de suporte a sistemas multicadeia esteja ativo
    _ensure_gmx_mmpbsa_cys_patched()

    # 1. Identificação do número de frames e cálculo da janela de equilíbrio (Últimos 40%)
    gmx_bin = find_executable("gmx")
    total_frames = 1000
    if gmx_bin:
        try:
            env = os.environ.copy()
            exec_dir = str(Path(gmx_bin).parent)
            env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

            chk_res = subprocess.run(
                [gmx_bin, "check", "-f", "md_fit.xtc"],
                cwd=str(working_dir),
                env=env,
                capture_output=True,
                text=True,
            )
            chk_out = (chk_res.stderr or "") + "\n" + (chk_res.stdout or "")
            frame_matches = re.findall(
                r"(?:Last\s+frame|Step|Coords|Time|Found)\s+(\d+)", chk_out, re.IGNORECASE
            )
            if frame_matches:
                f_count = max(int(m) for m in frame_matches)
                if f_count > 0:
                    total_frames = f_count
        except Exception:
            pass

    # Protocolo de Janela Termodinâmica: Últimos 40% (ex: 600 a 1000 para 1000 frames totais)
    startframe = max(1, int(round(total_frames * 0.60)))
    endframe = total_frames
    frames_in_window = max(1, endframe - startframe + 1)
    interval = max(1, frames_in_window // 200)
    if frames_in_window <= 400 and interval > 2:
        interval = 2
    elif interval < 1:
        interval = 1

    # 2. Criação do arquivo de entrada mmpbsa.in
    mmpbsa_in_path = working_dir / "mmpbsa.in"
    mmpbsa_in_content = f"""&general
sys_name="Protein_Ligand_Complex",
startframe={startframe},
endframe={endframe},
interval={interval},
verbose=2,
/
&gb
igb=5,
saltcon=0.150,
/
&pb
istrng=0.150,
fillratio=4.0,
radiopt=1,
/
"""
    with open(mmpbsa_in_path, "w", encoding="utf-8") as f:
        f.write(mmpbsa_in_content)

    # 3. Identificação dos índices dos grupos no index.ndx
    groups_list = get_index_groups(index_file)
    _, _, rec_idx, lig_idx = identify_complex_groups(groups_list)

    # 4. Execução do gmx_MMPBSA
    cmd_mmpbsa = [
        mmpbsa_bin,
        "-O",
        "-nogui",
        "-i",
        "mmpbsa.in",
        "-cs",
        "md.tpr",
        "-ct",
        "md_fit.xtc",
        "-ci",
        "index.ndx",
        "-cg",
        str(rec_idx),
        str(lig_idx),
        "-o",
        "FINAL_RESULTS_MMPBSA.dat",
        "-eo",
        "FINAL_RESULTS_MMPBSA.csv",
    ]

    # Identifica o arquivo mol2 parametrizado do ligante (gerado pelo ACPYPE)
    mol2_candidates = (
        list(working_dir.glob("*/*_bcc_gaff2.mol2"))
        + list(working_dir.glob("*_bcc_gaff2.mol2"))
        + list(working_dir.glob("*/*_AC.mol2"))
        + list(working_dir.glob("*/*.mol2"))
        + list(working_dir.glob("*.mol2"))
    )
    if mol2_candidates:
        mol2_rel = mol2_candidates[0].resolve().relative_to(working_dir.resolve())
        cmd_mmpbsa.extend(["-lm", str(mol2_rel).replace("\\", "/")])

    try:
        env = os.environ.copy()
        exec_dir = str(Path(mmpbsa_bin).parent)
        system_paths = ["/usr/bin", "/bin", "/usr/local/bin"]
        env["PATH"] = (
            f"{exec_dir}{os.pathsep}{os.pathsep.join(system_paths)}{os.pathsep}{env.get('PATH', '')}"
        )

        result = subprocess.run(
            cmd_mmpbsa, cwd=str(working_dir), env=env, capture_output=True, text=True
        )
        if result.returncode != 0:
            error_msg = result.stderr.strip() or result.stdout.strip()
            raise SimulationPrepError(
                f"Erro na execução do gmx_MMPBSA:\n"
                f"Comando: {' '.join(cmd_mmpbsa)}\n"
                f"Código de retorno: {result.returncode}\n"
                f"Erro real: {error_msg}"
            )
    except Exception as e:
        if isinstance(e, (SimulationPrepError, DependencyError, FileNotFoundError)):
            raise e
        raise SimulationPrepError(
            f"Falha ao executar o cálculo MM-PBSA via gmx_MMPBSA: {e}"
        )

    # 5. Parse dos resultados e geração de mmpbsa_summary.json
    dat_output = working_dir / "FINAL_RESULTS_MMPBSA.dat"
    summary_data = parse_mmpbsa_dat(dat_output)
    summary_data["thermodynamic_window"] = "60 - 100 ns (Últimos 40% - Estado Estacionário)"
    summary_data["startframe"] = startframe
    summary_data["endframe"] = endframe
    summary_data["interval"] = interval
    summary_data["total_frames_estimated"] = total_frames
    summary_data["frames_analyzed"] = max(1, (endframe - startframe + 1) // interval)
    summary_data["protocol"] = "Dupla Escala Temporal (MM-PBSA 60-100 ns / Trajetória 0-100 ns)"
    summary_data["raw_output_file"] = "FINAL_RESULTS_MMPBSA.dat"

    summary_json_path = working_dir / "mmpbsa_summary.json"
    with open(summary_json_path, "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=2, ensure_ascii=False)

    return summary_data
