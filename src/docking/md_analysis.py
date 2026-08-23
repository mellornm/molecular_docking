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


def fix_pbc(working_dir: Path) -> Path:
    """
    Executa o tratamento de Condições Periódicas de Contorno (PBC) via GROMACS (gmx trjconv):
    - Passo 1: Remove saltos e quebras através da caixa de simulação (-pbc nojump -> md_nojump.xtc)
    - Passo 2: Centraliza o complexo na proteína e compacta a caixa (-center -pbc mol -ur compact -> md_fit.xtc)

    Retorna o caminho do arquivo de trajetória corrigida md_fit.xtc.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

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

    # Passo 1: Remover saltos periódicos (nojump)
    # Seleção de grupo: 0 (System)
    cmd_nojump = [
        gmx_bin,
        "trjconv",
        "-s",
        "md.tpr",
        "-f",
        "md.xtc",
        "-o",
        "md_nojump.xtc",
        "-pbc",
        "nojump",
    ]
    run_trjconv_cmd(
        cmd_nojump, input_val="0\n", step_name="Passo 1 (Remoção de Saltos - nojump)"
    )

    # Passo 2: Centralizar na proteína e compactar a caixa
    # Seleção de grupos: 1 (Protein para centralizar) e 0 (System para salvar na trajetória final)
    cmd_fit = [
        gmx_bin,
        "trjconv",
        "-s",
        "md.tpr",
        "-f",
        "md_nojump.xtc",
        "-o",
        "md_fit.xtc",
        "-center",
        "-pbc",
        "mol",
        "-ur",
        "compact",
    ]
    run_trjconv_cmd(
        cmd_fit,
        input_val="1\n0\n",
        step_name="Passo 2 (Centralização e Compactação - mol compact)",
    )

    md_fit_path = working_dir / "md_fit.xtc"
    if not md_fit_path.exists():
        raise FileNotFoundError(
            f"Arquivo de trajetória corrigida 'md_fit.xtc' não foi gerado em: {working_dir}"
        )

    return md_fit_path


def analyze_trajectory(working_dir: Path):
    """
    Executa a análise da trajetória da Dinâmica Molecular no GROMACS utilizando
    estritamente a trajetória corrigida com PBC (md_fit.xtc).
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

    def run_analysis_cmd(cmd, cwd, input_val, step_name=""):
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

    # 1. RMSD do esqueleto da proteína (Backbone/Backbone -> 4 4)
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
    run_analysis_cmd(cmd_rmsd, working_dir, "4 4\n", "RMSD do esqueleto da proteína")

    # 2. RMSF por resíduo (C-alpha -> 3)
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
    run_analysis_cmd(cmd_rmsf, working_dir, "3\n", "RMSF por resíduo")

    # 3. Pontes de hidrogênio entre Proteína e Ligante
    cmd_hbond = [
        gmx_bin,
        "hbond",
        "-s",
        "md.tpr",
        "-f",
        "md_fit.xtc",
        "-num",
        "hbond.xvg",
    ]
    index_file = working_dir / "index.ndx"
    if index_file.exists():
        cmd_hbond.extend(["-n", "index.ndx"])
        groups_list = get_index_groups(index_file)
        prot_name, lig_name, _, _ = identify_complex_groups(groups_list)
        hbond_input = f"{prot_name}\n{lig_name}\n"
    else:
        hbond_input = "Protein\nLIG\n"

    run_analysis_cmd(
        cmd_hbond,
        working_dir,
        hbond_input,
        "Pontes de hidrogênio (Proteína-Ligante)",
    )


def parse_xvg(file_path: Path) -> Tuple[List[float], List[float]]:
    """
    Faz o parse de arquivos .xvg gerados pelo GROMACS, ignorando metadados iniciados com '@' e '#'.
    Retorna uma tupla contendo duas listas de floats (valores do eixo X e eixo Y).
    """
    x_vals: List[float] = []
    y_vals: List[float] = []

    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo XVG não encontrado: {file_path}")

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("@") or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    x_vals.append(float(parts[0]))
                    y_vals.append(float(parts[1]))
                except ValueError:
                    continue

    return x_vals, y_vals


def plot_md_results(working_dir: Path) -> Dict[str, Path]:
    """
    Lê os arquivos .xvg gerados na análise (rmsd.xvg, rmsf.xvg, hbond.xvg) e gera gráficos
    com padrão estético científico de publicação (300 DPI) utilizando matplotlib e seaborn.

    Salva diretamente no working_dir:
    - rmsd.png: Tempo (ns) vs RMSD (nm)
    - rmsf.png: Número do Resíduo vs Flutuação RMSF (nm)
    - hbond.png: Tempo (ns) vs Número de Pontes de Hidrogênio

    Retorna um dicionário mapeando o nome da análise ao caminho do arquivo .png gerado.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    # Configuração de estilo científico de alta qualidade
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

    # 1. Gráfico de RMSD
    rmsd_file = working_dir / "rmsd.xvg"
    if rmsd_file.exists():
        x_time, y_rmsd = parse_xvg(rmsd_file)
        if x_time and y_rmsd:
            fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=300)
            ax.plot(
                x_time, y_rmsd, color="#1f77b4", linewidth=1.5, label="Backbone RMSD"
            )

            # Adiciona linha de média
            mean_rmsd = sum(y_rmsd) / len(y_rmsd)
            ax.axhline(
                mean_rmsd,
                color="#d62728",
                linestyle="--",
                alpha=0.7,
                label=f"Média: {mean_rmsd:.3f} nm",
            )

            ax.set_xlabel("Tempo (ns)", fontweight="bold")
            ax.set_ylabel("RMSD (nm)", fontweight="bold")
            ax.set_title(
                "Evolução Temporal do RMSD do Esqueleto Protéico",
                fontweight="bold",
                pad=12,
            )
            ax.set_xlim(left=0, right=max(x_time) if x_time else 1)
            ax.set_ylim(bottom=0)
            ax.legend(loc="lower right", frameon=True, framealpha=0.9)

            out_rmsd_png = working_dir / "rmsd.png"
            fig.savefig(out_rmsd_png, dpi=300, bbox_inches="tight")  # type: ignore
            plt.close(fig)
            generated_plots["rmsd"] = out_rmsd_png

    # 2. Gráfico de RMSF
    rmsf_file = working_dir / "rmsf.xvg"
    if rmsf_file.exists():
        x_res, y_rmsf = parse_xvg(rmsf_file)
        if x_res and y_rmsf:
            fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=300)
            ax.plot(x_res, y_rmsf, color="#2a9d8f", linewidth=1.4, label="RMSF C-α")
            ax.fill_between(x_res, y_rmsf, color="#2a9d8f", alpha=0.25)  # type: ignore

            ax.set_xlabel("Número do Resíduo", fontweight="bold")
            ax.set_ylabel("Flutuação RMSF (nm)", fontweight="bold")
            ax.set_title(
                "Flutuação Atômica por Resíduo (RMSF)", fontweight="bold", pad=12
            )
            ax.set_xlim(
                left=min(x_res) if x_res else 0, right=max(x_res) if x_res else 1
            )
            ax.set_ylim(bottom=0)
            ax.legend(loc="upper right", frameon=True, framealpha=0.9)

            out_rmsf_png = working_dir / "rmsf.png"
            fig.savefig(out_rmsf_png, dpi=300, bbox_inches="tight")  # type: ignore
            plt.close(fig)
            generated_plots["rmsf"] = out_rmsf_png

    # 3. Gráfico de Pontes de Hidrogênio (HBond)
    hbond_file = working_dir / "hbond.xvg"
    if hbond_file.exists():
        x_time, y_hb = parse_xvg(hbond_file)
        if x_time and y_hb:
            fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=300)
            ax.plot(
                x_time,
                y_hb,
                color="#e76f51",
                linewidth=1.0,
                alpha=0.75,
                label="Pontes de H",
            )

            # Cálculo de média simples ou móvel se houver pontos suficientes
            if len(y_hb) >= 10:
                window_size = max(5, len(y_hb) // 25)
                # Média móvel manual pura em python
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
                    label=f"Média Suavizada (janela {window_size})",
                )

            ax.set_xlabel("Tempo (ns)", fontweight="bold")
            ax.set_ylabel("Número de H-Bonds", fontweight="bold")
            ax.set_title(
                "Pontes de Hidrogênio (Proteína - Ligante)", fontweight="bold", pad=12
            )
            ax.set_xlim(left=0, right=max(x_time) if x_time else 1)
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


def calculate_mmpbsa(working_dir: Path) -> Dict[str, Any]:
    """
    Executa o cálculo de Energia Livre de Ligação MM-PBSA via gmx_MMPBSA:
    1. Gera o arquivo de configuração mmpbsa.in com 100 frames distribuídos uniformemente ao longo de md_fit.xtc.
    2. Identifica os grupos do receptor (Protein) e ligante (ligand_md / LIG) em index.ndx.
    3. Executa o gmx_MMPBSA via subprocesso.
    4. Extrai contribuições energéticas (Van der Waals, Eletrostática, Solvatação Polar e Apolar) e Delta G.
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

    # 1. Identificação do número de frames e cálculo de intervalo para 100 frames uniformes
    gmx_bin = find_executable("gmx")
    total_frames = 100
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
                r"(?:Last\s+frame|Step|Coords|Time)\s+(\d+)", chk_out, re.IGNORECASE
            )
            if frame_matches:
                f_count = max(int(m) for m in frame_matches)
                if f_count > 0:
                    total_frames = f_count
        except Exception:
            pass

    interval = max(1, total_frames // 100)

    # 2. Criação do arquivo de entrada mmpbsa.in
    mmpbsa_in_path = working_dir / "mmpbsa.in"
    mmpbsa_in_content = f"""&general
sys_name="Protein_Ligand_Complex",
startframe=1,
endframe=999999,
interval={interval},
/
&gb
igb=5,
saltcon=0.150,
/
&pb
istrng=0.150,
fillratio=4.0,
radiopt=0,
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
        env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

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
    summary_data["frames_extracted_interval"] = interval
    summary_data["total_frames_estimated"] = total_frames
    summary_data["raw_output_file"] = "FINAL_RESULTS_MMPBSA.dat"

    summary_json_path = working_dir / "mmpbsa_summary.json"
    with open(summary_json_path, "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=2, ensure_ascii=False)

    return summary_data
