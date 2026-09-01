import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

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

from docking.md_prep import (
    DependencyError,
    SimulationPrepError,
    find_executable,
    sanitize_target_id,
    verify_tpr_consistency,
)


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


def export_cluster_package(
    working_dir: Path,
    target_id: Optional[str] = None,
    output_export_dir: Optional[Path] = None,
) -> Path:
    """
    Gera um diretório autocontido e isolado para execução remota/manual de Dinâmica Molecular:
    cluster_export/<target_id>/
    ├── <target_id>_md.tpr
    ├── run_local.sh (script executável pré-configurado para tmux/nohup)
    └── README.md (guia rápido de execução e monitoramento)
    """
    working_dir = Path(working_dir)
    if not target_id:
        target_id = sanitize_target_id(working_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            tpr_candidates = list(working_dir.glob("*_md.tpr"))
            if tpr_candidates:
                target_id = tpr_candidates[0].stem.replace("_md", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    # Localiza o TPR de produção
    tpr_source = working_dir / f"{prefix}md.tpr"
    if not tpr_source.exists():
        tpr_source = working_dir / "md.tpr"
    if not tpr_source.exists():
        raise FileNotFoundError(
            f"Arquivo TPR de produção não encontrado em {working_dir}. "
            "Compile o TPR primeiro via compile_production_tpr."
        )

    # Valida integridade do TPR
    verify_tpr_consistency(tpr_source)

    if output_export_dir is None:
        project_root = Path(__file__).resolve().parent.parent.parent
        export_root = project_root / "cluster_export"
    else:
        export_root = Path(output_export_dir)

    target_export_dir = export_root / target_id
    target_export_dir.mkdir(parents=True, exist_ok=True)

    # 1. Copia o TPR
    dest_tpr = target_export_dir / f"{prefix}md.tpr"
    shutil.copy2(tpr_source, dest_tpr)

    # 2. Gera o script executável run_local.sh
    run_script_content = f"""#!/usr/bin/env bash
# ==============================================================================
# Script de Execução Autônoma de Dinâmica Molecular (GROMACS)
# Receptor / Alvo: {target_id}
# Arquivo TPR: {dest_tpr.name}
# ==============================================================================

set -euo pipefail

TARGET="{target_id}"
TPR_FILE="${{TARGET}}_md.tpr"
DEFFNM="${{TARGET}}_md"

echo "======================================================================"
echo "  Iniciando Produção de Dinâmica Molecular: ${{TARGET}}"
echo "  Data / Hora de Início: $(date)"
echo "  Arquivo TPR: ${{TPR_FILE}}"
echo "  Diretório Atual: $(pwd)"
echo "======================================================================"

# 1. Identificação do binário do GROMACS
if command -v gmx &> /dev/null; then
    GMX_BIN="gmx"
elif [ -f "$HOME/miniforge3/envs/bioinfo/bin/gmx" ]; then
    GMX_BIN="$HOME/miniforge3/envs/bioinfo/bin/gmx"
elif [ -f "$HOME/miniconda3/envs/bioinfo/bin/gmx" ]; then
    GMX_BIN="$HOME/miniconda3/envs/bioinfo/bin/gmx"
else
    echo "ERRO: Executável 'gmx' (GROMACS) não encontrado no PATH nem no Conda bioinfo." >&2
    exit 1
fi

echo "-> Executável GROMACS: ${{GMX_BIN}}"

# 2. Verificação do arquivo TPR
if [ ! -f "${{TPR_FILE}}" ]; then
    echo "ERRO: Arquivo '${{TPR_FILE}}' não encontrado no diretório atual!" >&2
    exit 1
fi

# 3. Verificação de retomada automática via Checkpoint (cpi)
CPT_FLAG=""
if [ -f "${{DEFFNM}}.cpt" ]; then
    echo "-> Checkpoint prévio detectado (${{DEFFNM}}.cpt). Retomando simulação..."
    CPT_FLAG="-cpi ${{DEFFNM}}.cpt -append"
fi

# 4. Execução do mdrun
echo "-> Disparando mdrun com prefixo exclusivo '${{DEFFNM}}'..."
${{GMX_BIN}} mdrun -v \\
    -s "${{TPR_FILE}}" \\
    -deffnm "${{DEFFNM}}" \\
    ${{CPT_FLAG}}

echo "======================================================================"
echo "  Dinâmica Molecular Finalizada com Sucesso: ${{TARGET}}"
echo "  Data / Hora de Término: $(date)"
echo "  Arquivos Gerados: ${{DEFFNM}}.xtc, ${{DEFFNM}}.edr, ${{DEFFNM}}.gro, ${{DEFFNM}}.log"
echo "======================================================================"
"""

    script_path = target_export_dir / "run_local.sh"
    with open(script_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(run_script_content)

    try:
        os.chmod(script_path, 0o755)
    except Exception:
        pass

    # 3. Gera README.md com instruções de transferência e execução
    readme_content = f"""# Pacote de Execução Autônoma - Alvo: {target_id}

Este diretório contém os arquivos necessários para transferir e executar a simulação de Dinâmica Molecular (100 ns) em um servidor remoto, estação de trabalho ou GPU cluster sem depender de gerenciadores de fila (Slurm/PBS).

---

## 📦 Arquivos do Pacote
* `{dest_tpr.name}`: Arquivo binário de entrada da simulação (compilado e verificado).
* `run_local.sh`: Script bash otimizado para execução com tolerância a falhas e retomada via checkpoint.

---

## 🚀 Como Executar no Servidor Remoto

### 1. Transferir para o Servidor (SSH / SCP):
```bash
scp -r cluster_export/{target_id} usuario@servidor.remoto:/caminho/de/destino/
```

### 2. Acessar o Diretório:
```bash
cd /caminho/de/destino/{target_id}
```

### 3. Disparar a Simulação:

#### Opção A (Recomendado - via tmux):
```bash
tmux new -s md_{target_id}
bash run_local.sh
# Pressione Ctrl+B e depois D para desanexar do terminal
```

#### Opção B (via nohup em background):
```bash
nohup bash run_local.sh > simulation_{target_id}.log 2>&1 &
```

### 4. Acompanhar o Progresso:
```bash
tail -f {prefix}md.log
```

---

## 📥 Pós-Simulação
Após o término, copie de volta os arquivos `{prefix}md.xtc` e `{prefix}md.tpr` para `data/md_files/{target_id}/` e execute a Opção 9 do pipeline para análise, gráficos e MM-PBSA.
"""

    readme_path = target_export_dir / "README.md"
    with open(readme_path, "w", encoding="utf-8") as f:
        f.write(readme_content)

    return target_export_dir


def compile_production_tpr(
    working_dir: Path, target_id: Optional[str] = None
) -> Path:
    """
    Compila o arquivo de entrada binário para a produção da Dinâmica Molecular (<target_id>_md.tpr)
    usando o comando 'gmx grompp' a partir de npt.gro, npt.cpt, topol.top, index.ndx e md.mdp.

    Aplica:
    1. Isolamento estrito com prefixo do alvo.
    2. Deleção prévia de TPR antigo para evitar reutilização indevida.
    3. Controle estrito de saída do grompp (código 0).
    4. Checagem de consistência pós-geração via gmx check.
    5. Empacotamento automático para execução remota em cluster_export/<target_id>/.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    # Identificação do alvo
    if not target_id:
        target_id = sanitize_target_id(working_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            gro_files = list(working_dir.glob("*_npt.gro")) or list(working_dir.glob("*_em.gro"))
            if gro_files:
                target_id = gro_files[0].stem.replace("_npt", "").replace("_em", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    npt_gro = working_dir / f"{prefix}npt.gro"
    if not npt_gro.exists():
        npt_gro = working_dir / "npt.gro"
    if not npt_gro.exists():
        raise FileNotFoundError(
            f"Arquivo de equilíbrio 'npt.gro' não encontrado no diretório: {working_dir}"
        )

    npt_cpt = working_dir / f"{prefix}npt.cpt"
    if not npt_cpt.exists():
        npt_cpt = working_dir / "npt.cpt"

    topol_top = working_dir / f"{prefix}topol.top"
    if not topol_top.exists():
        topol_top = working_dir / "topol.top"
    if not topol_top.exists():
        raise FileNotFoundError(
            f"Arquivo de topologia 'topol.top' não encontrado no diretório: {working_dir}"
        )

    index_ndx = working_dir / f"{prefix}index.ndx"
    if not index_ndx.exists():
        index_ndx = working_dir / "index.ndx"

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError("O executável 'gmx' (GROMACS) não foi encontrado no PATH.")

    # Resolução do template md.mdp
    project_root = Path(__file__).resolve().parent.parent.parent
    md_mdp = project_root / "src" / "templates" / "mdp" / "md.mdp"
    if not md_mdp.exists():
        md_mdp = Path("src/templates/mdp/md.mdp").resolve()
        if not md_mdp.exists():
            raise FileNotFoundError("Arquivo template md.mdp não encontrado.")

    # Deleção prévia de TPRs antigos para impedir reaproveitamento indevido
    target_tpr = working_dir / f"{prefix}md.tpr"
    target_tpr.unlink(missing_ok=True)
    (working_dir / "md.tpr").unlink(missing_ok=True)

    # 1. Compilação do arquivo de produção (grompp)
    cmd_grompp = [
        gmx_bin,
        "grompp",
        "-f",
        str(md_mdp),
        "-c",
        str(npt_gro.name),
        "-p",
        str(topol_top.name),
        "-o",
        f"{prefix}md.tpr",
    ]
    if npt_cpt.exists():
        cmd_grompp.extend(["-t", str(npt_cpt.name)])
    if index_ndx.exists():
        cmd_grompp.extend(["-n", str(index_ndx.name)])

    try:
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env.pop("PYTHONHOME", None)
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

    if not target_tpr.exists():
        raise FileNotFoundError(f"Arquivo '{target_tpr.name}' não foi gerado em: {working_dir}")

    # Espelho md.tpr
    shutil.copy2(target_tpr, working_dir / "md.tpr")

    # Checagem de consistência pós-geração
    verify_tpr_consistency(target_tpr)

    # Empacotamento automático para exportação de cluster
    try:
        export_cluster_package(working_dir, target_id=target_id)
    except Exception:
        pass

    return target_tpr


def run_production_md(
    working_dir: Path, target_id: Optional[str] = None
):
    """
    Compila e executa a etapa de Produção de Dinâmica Molecular no GROMACS.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    # Identificação do alvo
    if not target_id:
        target_id = sanitize_target_id(working_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            gro_files = list(working_dir.glob("*_npt.gro")) or list(working_dir.glob("*_em.gro"))
            if gro_files:
                target_id = gro_files[0].stem.replace("_npt", "").replace("_em", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError("O executável 'gmx' (GROMACS) não foi encontrado no PATH.")

    # 1. Compilação do arquivo de produção (grompp -> <target_id>_md.tpr)
    compile_production_tpr(working_dir, target_id=target_id)

    # 2. Execução da simulação de produção (mdrun)
    cmd_mdrun = [gmx_bin, "mdrun", "-v", "-deffnm", f"{prefix}md"]

    try:
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env.pop("PYTHONHOME", None)
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

        # Espelhos para compatibilidade
        if (working_dir / f"{prefix}md.xtc").exists():
            shutil.copy2(working_dir / f"{prefix}md.xtc", working_dir / "md.xtc")
        if (working_dir / f"{prefix}md.gro").exists():
            shutil.copy2(working_dir / f"{prefix}md.gro", working_dir / "md.gro")

    except Exception as e:
        if isinstance(e, (SimulationPrepError, DependencyError)):
            raise e
        raise SimulationPrepError(
            f"Falha ao executar a Produção de Dinâmica Molecular (mdrun) em tempo real: {e}"
        )


def fix_pbc(
    working_dir: Path, target_id: Optional[str] = None, force: bool = False
) -> Path:
    """
    Executa o tratamento completo e automatizado de Condições Periódicas de Contorno (PBC)
    e remoção de movimentos de corpo rígido (fit rot+trans) via GROMACS (gmx trjconv):

    1. Etapa 1 - Centralização e Correção de PBC:
       gmx trjconv -s <target_id>_md.tpr -f <target_id>_md.xtc -o <target_id>_md_center.xtc -pbc mol -center
       (Input stdin: '1\\n0\\n' -> Protein para centralizar, System para salvar)

    2. Etapa 2 - Ajuste de Rotação e Translação (Fit de Mínimos Quadrados):
       gmx trjconv -s <target_id>_md.tpr -f <target_id>_md_center.xtc -o <target_id>_md_fit.xtc -fit rot+trans
       (Input stdin: '4\\n0\\n' -> Backbone para fit de RMSD, System para salvar)

    3. Etapa 3 - Estrutura Estática Limpa para PyMOL:
       gmx trjconv -s <target_id>_md.tpr -f <target_id>_md.gro -o <target_id>_md_clean.gro -pbc mol -center
       (Input stdin: '1\\n0\\n' -> Protein para centralizar, System para salvar)

    Retorna o caminho do arquivo de trajetória corrigida md_fit.xtc.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    # Identificação do alvo
    if not target_id:
        target_id = sanitize_target_id(working_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_md.tpr")) or list(working_dir.glob("*_md.xtc"))
            if candidates:
                target_id = candidates[0].stem.replace("_md", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    md_fit_path = working_dir / f"{prefix}md_fit.xtc"
    md_clean_path = working_dir / f"{prefix}md_clean.gro"

    if not force and md_fit_path.exists() and md_clean_path.exists():
        if not (working_dir / "md_fit.xtc").exists():
            shutil.copy2(md_fit_path, working_dir / "md_fit.xtc")
        if not (working_dir / "md_clean.gro").exists():
            shutil.copy2(md_clean_path, working_dir / "md_clean.gro")
        return md_fit_path

    if not force and (working_dir / "md_fit.xtc").exists() and (working_dir / "md_clean.gro").exists():
        return working_dir / "md_fit.xtc"

    tpr_file = working_dir / f"{prefix}md.tpr"
    if not tpr_file.exists():
        tpr_file = working_dir / "md.tpr"

    xtc_file = working_dir / f"{prefix}md.xtc"
    if not xtc_file.exists():
        xtc_file = working_dir / "md.xtc"

    if not tpr_file.exists():
        raise FileNotFoundError(
            f"Arquivo de topologia compilada '{tpr_file.name}' não encontrado em: {working_dir}"
        )
    if not xtc_file.exists():
        raise FileNotFoundError(
            f"Arquivo de trajetória de produção '{xtc_file.name}' não encontrado em: {working_dir}"
        )

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError("O executável 'gmx' (GROMACS) não foi encontrado no PATH.")

    def run_trjconv_cmd(cmd: List[str], input_val: str, step_name: str):
        try:
            env = os.environ.copy()
            env.pop("PYTHONPATH", None)
            env.pop("PYTHONHOME", None)
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
    center_xtc_name = f"{prefix}md_center.xtc"
    cmd_center = [
        gmx_bin,
        "trjconv",
        "-s",
        str(tpr_file.name),
        "-f",
        str(xtc_file.name),
        "-o",
        center_xtc_name,
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
    fit_xtc_name = f"{prefix}md_fit.xtc"
    cmd_fit = [
        gmx_bin,
        "trjconv",
        "-s",
        str(tpr_file.name),
        "-f",
        center_xtc_name,
        "-o",
        fit_xtc_name,
        "-fit",
        "rot+trans",
    ]
    run_trjconv_cmd(
        cmd_fit,
        input_val="4\n0\n",
        step_name="Etapa 2 (Ajuste de Rotação e Translação - fit rot+trans)",
    )

    # Etapa 3: Estrutura Estática Limpa para Visualização no PyMOL
    gro_candidates = [
        f"{prefix}md.gro", "md.gro",
        f"{prefix}npt.gro", "npt.gro",
        f"{prefix}em.gro", "em.gro",
        f"{prefix}complex.gro", "complex.gro"
    ]
    gro_source = None
    for cand in gro_candidates:
        if (working_dir / cand).exists():
            gro_source = cand
            break

    clean_gro_name = f"{prefix}md_clean.gro"
    if gro_source:
        cmd_clean_gro = [
            gmx_bin,
            "trjconv",
            "-s",
            str(tpr_file.name),
            "-f",
            gro_source,
            "-o",
            clean_gro_name,
            "-pbc",
            "mol",
            "-center",
        ]
        run_trjconv_cmd(
            cmd_clean_gro,
            input_val="1\n0\n",
            step_name="Etapa 3 (Estrutura Estática Limpa para PyMOL - md_clean.gro)",
        )

    # Cria espelhos
    if (working_dir / fit_xtc_name).exists():
        shutil.copy2(working_dir / fit_xtc_name, working_dir / "md_fit.xtc")
    if (working_dir / clean_gro_name).exists():
        shutil.copy2(working_dir / clean_gro_name, working_dir / "md_clean.gro")

    if not md_fit_path.exists() and not (working_dir / "md_fit.xtc").exists():
        raise FileNotFoundError(
            f"Arquivo de trajetória corrigida 'md_fit.xtc' não foi gerado em: {working_dir}"
        )

    return md_fit_path if md_fit_path.exists() else working_dir / "md_fit.xtc"


def analyze_trajectory(
    working_dir: Path, target_id: Optional[str] = None
):
    """
    Executa a análise quantitativa da trajetória de Dinâmica Molecular no GROMACS (Janela Completa: 0 - 100 ns):
    1. RMSD do Backbone da Proteína e do Ligante
    2. RMSF por resíduo dos Carbonos Alfa
    3. Monitoramento de Pontes de Hidrogênio Proteína-Ligante
    4. Raio de Giro (Rg)
    5. Área de Superfície Acessível ao Solvente (SASA)
    6. Ocupação Temporal de Pontes de Hidrogênio
    7. Agrupamento Conformacional (Clustering)
    8. Exportação automatizada de todas as matrizes brutas em formato CSV para publicação.

    Utiliza estritamente a trajetória ajustada e sem artefatos de PBC.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    # Identificação do alvo
    if not target_id:
        target_id = sanitize_target_id(working_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_md.tpr")) or list(working_dir.glob("*_md_fit.xtc"))
            if candidates:
                target_id = candidates[0].stem.replace("_md_fit", "").replace("_md", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    tpr_file = working_dir / f"{prefix}md.tpr"
    if not tpr_file.exists():
        tpr_file = working_dir / "md.tpr"

    fit_xtc = working_dir / f"{prefix}md_fit.xtc"
    if not fit_xtc.exists():
        fit_xtc = working_dir / "md_fit.xtc"

    if not tpr_file.exists():
        raise FileNotFoundError(
            f"Arquivo '{tpr_file.name}' não encontrado no diretório: {working_dir}"
        )
    if not fit_xtc.exists():
        raise FileNotFoundError(
            f"Trajetória corrigida '{fit_xtc.name}' não encontrada em: {working_dir}. "
            "Certifique-se de executar o tratamento de PBC (fix_pbc) antes da análise."
        )

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        raise DependencyError("O executável 'gmx' (GROMACS) não foi encontrado no PATH.")

    def run_analysis_cmd(cmd: List[str], cwd: Path, input_val: str, step_name: str = ""):
        try:
            env = os.environ.copy()
            env.pop("PYTHONPATH", None)
            env.pop("PYTHONHOME", None)
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

    index_file = working_dir / f"{prefix}index.ndx"
    if not index_file.exists():
        index_file = working_dir / "index.ndx"

    prot_name = "Protein"
    lig_name = "LIG"
    prot_idx = 1
    lig_idx = 13
    if index_file.exists():
        groups_list = get_index_groups(index_file)
        prot_name, lig_name, prot_idx, lig_idx = identify_complex_groups(groups_list)

    # 1.1 RMSD do esqueleto da proteína (Backbone/Backbone -> 4 4)
    rmsd_xvg = f"{prefix}rmsd.xvg"
    cmd_rmsd = [
        gmx_bin,
        "rms",
        "-s",
        str(tpr_file.name),
        "-f",
        str(fit_xtc.name),
        "-o",
        rmsd_xvg,
        "-tu",
        "ns",
    ]
    run_analysis_cmd(cmd_rmsd, working_dir, "4\n4\n", "RMSD do esqueleto da proteína (Backbone)")
    shutil.copy2(working_dir / rmsd_xvg, working_dir / "rmsd.xvg")

    # 1.2 RMSD do Ligante no sítio ativo (Fit: Backbone 4, Calc: Ligand)
    if index_file.exists() and lig_idx is not None:
        try:
            rmsd_lig_xvg = f"{prefix}rmsd_lig.xvg"
            cmd_rmsd_lig = [
                gmx_bin,
                "rms",
                "-s",
                str(tpr_file.name),
                "-f",
                str(fit_xtc.name),
                "-o",
                rmsd_lig_xvg,
                "-tu",
                "ns",
                "-n",
                str(index_file.name),
            ]
            run_analysis_cmd(
                cmd_rmsd_lig,
                working_dir,
                f"4\n{lig_idx}\n",
                f"RMSD do ligante ({lig_name}) no sítio ativo",
            )
            shutil.copy2(working_dir / rmsd_lig_xvg, working_dir / "rmsd_lig.xvg")
        except Exception:
            pass

    # 2. RMSF por resíduo dos Carbonos Alfa (C-alpha -> 3)
    rmsf_xvg = f"{prefix}rmsf.xvg"
    cmd_rmsf = [
        gmx_bin,
        "rmsf",
        "-s",
        str(tpr_file.name),
        "-f",
        str(fit_xtc.name),
        "-o",
        rmsf_xvg,
        "-res",
    ]
    run_analysis_cmd(cmd_rmsf, working_dir, "3\n", "RMSF por resíduo (C-alpha)")
    shutil.copy2(working_dir / rmsf_xvg, working_dir / "rmsf.xvg")

    # 3. Pontes de hidrogênio entre Proteína e Ligante ao longo dos 100 ns
    hbond_xvg = f"{prefix}hbond.xvg"
    cmd_hbond = [
        gmx_bin,
        "hbond",
        "-s",
        str(tpr_file.name),
        "-f",
        str(fit_xtc.name),
        "-num",
        hbond_xvg,
        "-hbn",
        f"{prefix}hbond.ndx",
        "-logfile",
        f"{prefix}hbond.log",
        "-tu",
        "ns",
    ]
    if index_file.exists():
        cmd_hbond.extend(["-n", str(index_file.name)])
        hbond_input = f"{prot_name}\n{lig_name}\n"
    else:
        hbond_input = "Protein\nLIG\n"

    try:
        run_analysis_cmd(
            cmd_hbond,
            working_dir,
            hbond_input,
            "Pontes de hidrogênio (Proteína-Ligante)",
        )
    except Exception:
        cmd_hbond_fallback = [
            gmx_bin,
            "hbond",
            "-s",
            str(tpr_file.name),
            "-f",
            str(fit_xtc.name),
            "-num",
            hbond_xvg,
            "-tu",
            "ns",
        ]
        if index_file.exists():
            cmd_hbond_fallback.extend(["-n", str(index_file.name)])
        run_analysis_cmd(
            cmd_hbond_fallback,
            working_dir,
            hbond_input,
            "Pontes de hidrogênio (Proteína-Ligante)",
        )
    shutil.copy2(working_dir / hbond_xvg, working_dir / "hbond.xvg")

    # 4. Raio de Giro (Rg - Compacidade e Estabilidade Global de Enovelamento)
    gyrate_xvg = f"{prefix}gyrate.xvg"
    cmd_gyrate = [
        gmx_bin,
        "gyrate",
        "-s",
        str(tpr_file.name),
        "-f",
        str(fit_xtc.name),
        "-o",
        gyrate_xvg,
    ]
    try:
        run_analysis_cmd(cmd_gyrate, working_dir, "1\n", "Raio de Giro (Rg - Proteína)")
        shutil.copy2(working_dir / gyrate_xvg, working_dir / "gyrate.xvg")
    except Exception:
        pass

    # 5. Área de Superfície Acessível ao Solvente (SASA)
    sasa_xvg = f"{prefix}sasa.xvg"
    cmd_sasa = [
        gmx_bin,
        "sasa",
        "-s",
        str(tpr_file.name),
        "-f",
        str(fit_xtc.name),
        "-o",
        sasa_xvg,
        "-tu",
        "ns",
    ]
    try:
        run_analysis_cmd(cmd_sasa, working_dir, "1\n", "Área de Superfície Acessível ao Solvente (SASA)")
        shutil.copy2(working_dir / sasa_xvg, working_dir / "sasa.xvg")
    except Exception:
        pass

    # 6. Quantificação de persistência temporal (% H-Bond Occupancy)
    parse_hbond_occupancy(working_dir)

    # 7. Agrupamento Conformacional GROMOS (gmx cluster) para extração da pose representativa
    calculate_clusters(working_dir, target_id=target_id)

    # 8. Exportação automatizada de todas as matrizes brutas em formato CSV para publicação
    export_analysis_csv(working_dir, target_id=target_id)


def parse_xvg_multicolumn(
    file_path: Path,
) -> Tuple[List[float], Dict[str, List[float]], Dict[str, str]]:
    """
    Faz o parse de arquivos .xvg multicolunas gerados pelo GROMACS (ex: gyrate.xvg, sasa.xvg).
    Retorna (x_vals, y_series_dict, metadata_dict).
    """
    x_vals: List[float] = []
    series_data: List[List[float]] = []
    meta: Dict[str, str] = {}
    legends: Dict[int, str] = {}

    if not file_path.exists():
        raise FileNotFoundError(f"Arquivo XVG não encontrado: {file_path}")

    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith("@"):
                line_clean = line_str[1:].strip()
                if "xaxis" in line_clean and "label" in line_clean:
                    parts = line_clean.split("label", 1)
                    if len(parts) > 1:
                        meta["xaxis_label"] = parts[1].strip().strip('"')
                elif "yaxis" in line_clean and "label" in line_clean:
                    parts = line_clean.split("label", 1)
                    if len(parts) > 1:
                        meta["yaxis_label"] = parts[1].strip().strip('"')
                elif "legend" in line_clean:
                    m = re.search(r's(\d+)\s+legend\s+"([^"]+)"', line_clean)
                    if m:
                        legends[int(m.group(1))] = m.group(2)
                continue
            if line_str.startswith("#"):
                continue

            parts = line_str.split()
            if len(parts) >= 2:
                try:
                    x_val = float(parts[0])
                    y_row = [float(p) for p in parts[1:]]
                    x_vals.append(x_val)
                    if not series_data:
                        series_data = [[] for _ in range(len(y_row))]
                    for s_idx, val in enumerate(y_row):
                        if s_idx < len(series_data):
                            series_data[s_idx].append(val)
                except ValueError:
                    continue

    y_series_dict: Dict[str, List[float]] = {}
    for s_idx, col_vals in enumerate(series_data):
        label = legends.get(s_idx, f"Series_{s_idx + 1}")
        y_series_dict[label] = col_vals

    return x_vals, y_series_dict, meta


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


def plot_md_results(
    working_dir: Path, target_id: Optional[str] = None
) -> Dict[str, Path]:
    """
    Lê os arquivos .xvg gerados na análise e gera gráficos
    com padrão estético científico de publicação (300 DPI) cobrindo a extensão completa da simulação (0 - 100 ns).

    Salva no working_dir os gráficos com prefixo do alvo (<target_id>_*.png) e espelhos (*.png).
    Retorna um dicionário mapeando o nome da análise ao caminho do arquivo .png gerado.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    # Identificação do alvo
    if not target_id:
        target_id = sanitize_target_id(working_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_rmsd.xvg")) or list(working_dir.glob("*_md.tpr"))
            if candidates:
                target_id = candidates[0].stem.replace("_rmsd", "").replace("_md", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

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
    rmsd_file = working_dir / f"{prefix}rmsd.xvg"
    if not rmsd_file.exists():
        rmsd_file = working_dir / "rmsd.xvg"

    rmsd_lig_file = working_dir / f"{prefix}rmsd_lig.xvg"
    if not rmsd_lig_file.exists():
        rmsd_lig_file = working_dir / "rmsd_lig.xvg"

    if rmsd_file.exists():
        x_time, y_rmsd, meta_rmsd = parse_xvg_with_meta(rmsd_file)
        if x_time and y_rmsd:
            sorted_rmsd = sorted(zip(x_time, y_rmsd), key=lambda p: p[0])
            x_time = [p[0] for p in sorted_rmsd]
            y_rmsd = [p[1] for p in sorted_rmsd]

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

            mean_rmsd = sum(y_rmsd) / len(y_rmsd)
            ax.axhline(
                mean_rmsd,
                color="#1f77b4",
                linestyle="--",
                alpha=0.6,
                label=f"Protein Mean: {mean_rmsd:.3f} nm",
            )

            max_x = max(x_time) if x_time else 100.0

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
                f"Structural Stability & Ligand Residence ({target_id} | 0 - 100 ns)",
                fontweight="bold",
                pad=12,
            )
            ax.set_xlim(left=0, right=max_x)
            ax.set_ylim(bottom=0)
            ax.legend(loc="upper left", frameon=True, framealpha=0.9)

            out_rmsd_png = working_dir / f"{prefix}rmsd.png"
            fig.savefig(out_rmsd_png, dpi=300, bbox_inches="tight")  # type: ignore
            plt.close(fig)
            shutil.copy2(out_rmsd_png, working_dir / "rmsd.png")
            generated_plots["rmsd"] = out_rmsd_png

    # 2. Gráfico de RMSF (Janela Completa: 0 - 100 ns)
    rmsf_file = working_dir / f"{prefix}rmsf.xvg"
    if not rmsf_file.exists():
        rmsf_file = working_dir / "rmsf.xvg"

    if rmsf_file.exists():
        x_res, y_rmsf, _ = parse_xvg_with_meta(rmsf_file)
        if x_res and y_rmsf:
            fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=300)

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
                f"Root Mean Square Fluctuation per Residue ({target_id} | 0 - 100 ns)",
                fontweight="bold",
                pad=12,
            )
            ax.set_xlim(
                left=min(x_plot) if x_plot else 0, right=max(x_plot) if x_plot else 1
            )
            ax.set_ylim(bottom=0)
            ax.legend(loc="upper right", frameon=True, framealpha=0.9)

            out_rmsf_png = working_dir / f"{prefix}rmsf.png"
            fig.savefig(out_rmsf_png, dpi=300, bbox_inches="tight")  # type: ignore
            plt.close(fig)
            shutil.copy2(out_rmsf_png, working_dir / "rmsf.png")
            generated_plots["rmsf"] = out_rmsf_png

    # 3. Gráfico de Pontes de Hidrogênio (HBond) (Janela Completa: 0 - 100 ns)
    hbond_file = working_dir / f"{prefix}hbond.xvg"
    if not hbond_file.exists():
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
            ax.set_title(f"Protein–Ligand Hydrogen Bonds ({target_id} | 0 - 100 ns)", fontweight="bold", pad=12)
            ax.set_xlim(left=0, right=max(x_time) if x_time else 100.0)
            ax.set_ylim(bottom=0)
            ax.legend(loc="upper right", frameon=True, framealpha=0.9)

            out_hbond_png = working_dir / f"{prefix}hbond.png"
            fig.savefig(out_hbond_png, dpi=300, bbox_inches="tight")  # type: ignore
            plt.close(fig)
            shutil.copy2(out_hbond_png, working_dir / "hbond.png")
            generated_plots["hbond"] = out_hbond_png

    # 4. Gráfico de Raio de Giro (Rg - Compacidade e Enovelamento)
    gyrate_file = working_dir / f"{prefix}gyrate.xvg"
    if not gyrate_file.exists():
        gyrate_file = working_dir / "gyrate.xvg"

    if gyrate_file.exists():
        try:
            x_time, y_dict, meta_gyrate = parse_xvg_multicolumn(gyrate_file)
            if x_time and y_dict:
                xaxis_lbl = meta_gyrate.get("xaxis_label", "").lower()
                if "ps" in xaxis_lbl or (max(x_time) > 1000 and "ns" not in xaxis_lbl):
                    x_time = [t / 1000.0 for t in x_time]

                fig, ax = plt.subplots(figsize=(8.0, 4.5), dpi=300)
                total_key = next((k for k in y_dict.keys() if "total" in k.lower() or "rg" in k.lower()), list(y_dict.keys())[0])
                y_total = y_dict[total_key]
                ax.plot(x_time, y_total, color="#457b9d", linewidth=1.4, alpha=0.85, label=f"Total Rg ({total_key})")

                mean_rg = sum(y_total) / len(y_total)
                ax.axhline(mean_rg, color="#1d3557", linestyle="--", alpha=0.7, label=f"Mean Rg: {mean_rg:.3f} nm")

                ax.set_xlabel("Time (ns)", fontweight="bold")
                ax.set_ylabel("Radius of Gyration (nm)", fontweight="bold")
                ax.set_title(f"Protein Compactness & Folding Stability ({target_id} | 0 - 100 ns)", fontweight="bold", pad=12)
                ax.set_xlim(left=0, right=max(x_time) if x_time else 100.0)
                ax.legend(loc="upper right", frameon=True, framealpha=0.9)

                out_gyrate_png = working_dir / f"{prefix}gyrate.png"
                fig.savefig(out_gyrate_png, dpi=300, bbox_inches="tight")
                plt.close(fig)
                shutil.copy2(out_gyrate_png, working_dir / "gyrate.png")
                generated_plots["gyrate"] = out_gyrate_png
        except Exception:
            pass

    # 5. Gráfico de Área de Superfície Acessível ao Solvente (SASA)
    sasa_file = working_dir / f"{prefix}sasa.xvg"
    if not sasa_file.exists():
        sasa_file = working_dir / "sasa.xvg"

    if sasa_file.exists():
        try:
            x_time, y_dict, meta_sasa = parse_xvg_multicolumn(sasa_file)
            if x_time and y_dict:
                xaxis_lbl = meta_sasa.get("xaxis_label", "").lower()
                if "ps" in xaxis_lbl or (max(x_time) > 1000 and "ns" not in xaxis_lbl):
                    x_time = [t / 1000.0 for t in x_time]

                fig, ax = plt.subplots(figsize=(8.0, 4.5), dpi=300)
                total_key = list(y_dict.keys())[0]
                y_sasa = y_dict[total_key]
                ax.plot(x_time, y_sasa, color="#2b9348", linewidth=1.4, alpha=0.85, label="Total SASA")

                mean_sasa = sum(y_sasa) / len(y_sasa)
                ax.axhline(mean_sasa, color="#007f5f", linestyle="--", alpha=0.7, label=f"Mean SASA: {mean_sasa:.2f} nm²")

                ax.set_xlabel("Time (ns)", fontweight="bold")
                ax.set_ylabel(r"SASA ($\mathrm{nm}^2$)", fontweight="bold")
                ax.set_title(f"Solvent Accessible Surface Area ({target_id} | 0 - 100 ns)", fontweight="bold", pad=12)
                ax.set_xlim(left=0, right=max(x_time) if x_time else 100.0)
                ax.legend(loc="upper right", frameon=True, framealpha=0.9)

                out_sasa_png = working_dir / f"{prefix}sasa.png"
                fig.savefig(out_sasa_png, dpi=300, bbox_inches="tight")
                plt.close(fig)
                shutil.copy2(out_sasa_png, working_dir / "sasa.png")
                generated_plots["sasa"] = out_sasa_png
        except Exception:
            pass

    # 6. Gráfico de Decomposição MM-PBSA por Resíduo (Hotspots Energéticos) se disponível
    decomp_dat_file = working_dir / f"{prefix}FINAL_DECOMP_MMPBSA.dat"
    if not decomp_dat_file.exists():
        decomp_dat_file = working_dir / "FINAL_DECOMP_MMPBSA.dat"

    if decomp_dat_file.exists():
        decomp_data = parse_mmpbsa_decomp(decomp_dat_file)
        if decomp_data:
            out_decomp_png = plot_mmpbsa_decomp(decomp_data, working_dir, target_id=target_id)
            if out_decomp_png:
                generated_plots["decomp"] = out_decomp_png

    if not generated_plots:
        raise FileNotFoundError(
            f"Nenhum arquivo de análise (.xvg) foi encontrado em {working_dir} para geração dos gráficos."
        )

    return generated_plots


def calculate_clusters(
    working_dir: Path, target_id: Optional[str] = None, cutoff: float = 0.15
) -> Optional[Path]:
    """
    Executa o agrupamento conformacional da trajetória via algoritmo GROMOS no GROMACS (gmx cluster).
    Gera a estrutura medóide mais representativa do estado estacionário (<target_id>_cluster_medoid.gro).
    """
    working_dir = Path(working_dir)
    if not target_id:
        target_id = sanitize_target_id(working_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_md.tpr")) or list(working_dir.glob("*_md_fit.xtc"))
            if candidates:
                target_id = candidates[0].stem.replace("_md_fit", "").replace("_md", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    tpr_file = working_dir / f"{prefix}md.tpr"
    if not tpr_file.exists():
        tpr_file = working_dir / "md.tpr"

    fit_xtc = working_dir / f"{prefix}md_fit.xtc"
    if not fit_xtc.exists():
        fit_xtc = working_dir / "md_fit.xtc"

    if not tpr_file.exists() or not fit_xtc.exists():
        return None

    gmx_bin = find_executable("gmx")
    if not gmx_bin:
        return None

    medoid_name = f"{prefix}cluster_medoid.gro"
    log_name = f"{prefix}cluster.log"
    dist_name = f"{prefix}clust-dist.xvg"

    cmd_cluster = [
        gmx_bin,
        "cluster",
        "-s",
        str(tpr_file.name),
        "-f",
        str(fit_xtc.name),
        "-method",
        "gromos",
        "-cutoff",
        str(cutoff),
        "-cl",
        medoid_name,
        "-g",
        log_name,
        "-dist",
        dist_name,
    ]

    try:
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env.pop("PYTHONHOME", None)
        exec_dir = str(Path(gmx_bin).parent)
        env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

        subprocess.run(
            cmd_cluster,
            cwd=str(working_dir),
            env=env,
            capture_output=True,
            text=True,
            input="1\n1\n",
        )
        medoid_path = working_dir / medoid_name
        if medoid_path.exists():
            shutil.copy2(medoid_path, working_dir / "cluster_medoid.gro")
            return medoid_path
        if (working_dir / "cluster_medoid.gro").exists():
            return working_dir / "cluster_medoid.gro"
    except Exception:
        pass

    return None


def parse_hbond_occupancy(working_dir: Path, total_frames: int = 1000) -> List[Dict[str, Any]]:
    """
    Faz o parsing das pontes de hidrogênio geradas pelo GROMACS (hbond.log / hbond.ndx)
    para quantificar a persistência temporal e ocupação percentual por par de resíduos.
    Salva os resultados estruturados em 'hbond_occupancy.json'.
    """
    working_dir = Path(working_dir)
    log_file = working_dir / "hbond.log"
    occupancy_list: List[Dict[str, Any]] = []

    if log_file.exists():
        try:
            with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()

            for line in lines:
                line_str = line.strip()
                if "%" in line_str and ("-" in line_str or "atom" in line_str.lower() or "res" in line_str.lower() or "side" in line_str.lower() or "main" in line_str.lower()):
                    parts = line_str.split()
                    pct_str = [p for p in parts if "%" in p]
                    if pct_str:
                        try:
                            pct_val = float(pct_str[0].replace("%", "").replace(",", "."))
                            donor = parts[0] if len(parts) > 0 else "UNK"
                            acceptor = parts[1] if len(parts) > 1 else "LIG"
                            occupancy_list.append({
                                "donor": donor,
                                "acceptor": acceptor,
                                "occupancy_percent": round(pct_val, 2),
                                "classification": "Permanente / Âncora Farmacofórica" if pct_val >= 75.0 else ("Moderada" if pct_val >= 35.0 else "Transitória")
                            })
                        except ValueError:
                            continue
        except Exception:
            pass

    json_path = working_dir / "hbond_occupancy.json"
    if occupancy_list:
        occupancy_list.sort(key=lambda x: x["occupancy_percent"], reverse=True)
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(occupancy_list, f, indent=2, ensure_ascii=False)

    return occupancy_list


def parse_mmpbsa_decomp(decomp_path: Path) -> List[Dict[str, Any]]:
    """
    Parse estruturado do arquivo de decomposição por resíduo (FINAL_DECOMP_MMPBSA.dat).
    Extrai as contribuições energéticas individuais por resíduo (Van der Waals, Eletrostática, Polar, Apolar e Total).
    Retorna uma lista ordenada pela energia total (mais favoráveis / estabilizadoras primeiro).
    """
    if not decomp_path.exists():
        return []

    residues_data: List[Dict[str, Any]] = []

    with open(decomp_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    in_total_decomp = False
    for line in lines:
        line_clean = line.strip()
        if not line_clean or line_clean.startswith("#") or line_clean.startswith("---") or line_clean.startswith("|"):
            continue

        if "Total Energy Decomposition:" in line_clean or "DELTAS:" in line_clean or "Energy Decomposition Analysis" in line_clean or "TDC" in line_clean:
            in_total_decomp = True
            continue

        # Sai da seção Total Energy Decomposition se encontrar outra seção como Sidechain ou Backbone
        if in_total_decomp and ("Sidechain Energy Decomposition:" in line_clean or "Backbone Energy Decomposition:" in line_clean):
            break

        if in_total_decomp:
            if line_clean.startswith("Residue") or line_clean.startswith(",Avg") or line_clean.startswith("==="):
                continue

            # 1. Formato CSV do gmx_MMPBSA (Residue, Internal, vdW, Elec, Polar, Non-Polar, Total)
            if "," in line_clean:
                parts = [p.strip() for p in line_clean.split(",")]
                if len(parts) >= 17:
                    res_raw = parts[0]
                    if not res_raw or res_raw.lower().startswith("residue") or res_raw.lower().startswith("avg"):
                        continue
                    try:
                        vdw_mean = float(parts[4])
                        eel_mean = float(parts[7])
                        polar_mean = float(parts[10])
                        apolar_mean = float(parts[13])
                        total_mean = float(parts[16])
                        total_std = float(parts[17]) if len(parts) > 17 and parts[17] else 0.0
                        total_sem = float(parts[18]) if len(parts) > 18 and parts[18] else 0.0

                        residues_data.append({
                            "residue": res_raw,
                            "vdw": round(vdw_mean, 3),
                            "electrostatic": round(eel_mean, 3),
                            "polar": round(polar_mean, 3),
                            "nonpolar": round(apolar_mean, 3),
                            "total": round(total_mean, 3),
                            "std": round(total_std, 3),
                            "sem": round(total_sem, 3),
                        })
                        continue
                    except ValueError:
                        pass

            # 2. Formato clássico com pipe (|)
            if "|" in line_clean:
                cols = [c.strip() for c in line_clean.split("|")]
                if len(cols) >= 4:
                    res_raw = cols[0]
                    try:
                        total_col = cols[-1]
                        total_parts = total_col.replace("±", "+/-").split("+/-")
                        total_mean = float(total_parts[0].strip())
                        total_std = float(total_parts[1].strip()) if len(total_parts) > 1 else 0.0

                        vdw_mean = 0.0
                        if len(cols) > 2:
                            vdw_parts = cols[2].replace("±", "+/-").split("+/-")
                            try:
                                vdw_mean = float(vdw_parts[0].strip())
                            except ValueError:
                                vdw_mean = 0.0

                        eel_mean = 0.0
                        if len(cols) > 3:
                            eel_parts = cols[3].replace("±", "+/-").split("+/-")
                            try:
                                eel_mean = float(eel_parts[0].strip())
                            except ValueError:
                                eel_mean = 0.0

                        residues_data.append({
                            "residue": res_raw,
                            "vdw": round(vdw_mean, 3),
                            "electrostatic": round(eel_mean, 3),
                            "total": round(total_mean, 3),
                            "std": round(total_std, 3),
                            "sem": 0.0,
                        })
                    except ValueError:
                        continue
            else:
                parts = line_clean.split()
                if len(parts) >= 2:
                    res_raw = parts[0]
                    try:
                        total_mean = float(parts[-1])
                        residues_data.append({
                            "residue": res_raw,
                            "vdw": 0.0,
                            "electrostatic": 0.0,
                            "total": round(total_mean, 3),
                            "std": 0.0,
                            "sem": 0.0,
                        })
                    except ValueError:
                        continue

    # Ordena os resíduos: mais estabilizadores primeiro (delta G mais negativo)
    residues_data.sort(key=lambda x: x["total"])
    return residues_data


def plot_mmpbsa_decomp(
    decomp_data: List[Dict[str, Any]],
    working_dir: Path,
    target_id: Optional[str] = None,
) -> Optional[Path]:
    """
    Gera gráfico de barras de publicação (decomp_mmpbsa.png a 300 DPI) destacando os resíduos chave (hotspots)
    que mais contribuem para a energia livre de ligação MM-PBSA.
    """
    if not decomp_data:
        return None

    working_dir = Path(working_dir)
    if not target_id:
        target_id = sanitize_target_id(working_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_FINAL_DECOMP_MMPBSA.dat")) or list(working_dir.glob("*_md.tpr"))
            if candidates:
                target_id = candidates[0].stem.replace("_FINAL_DECOMP_MMPBSA", "").replace("_md", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    # Filtra os resíduos mais estabilizadores (< 0) e principais desestabilizadores (> 0.5 kcal/mol)
    stabilizing = [r for r in decomp_data if r["total"] < 0][:15]
    destabilizing = [r for r in decomp_data if r["total"] > 0.3][:5]
    plot_items = stabilizing + destabilizing
    if not plot_items:
        plot_items = decomp_data[:15]

    plot_items.sort(key=lambda x: x["total"])

    labels = []
    for r in plot_items:
        raw_res = r.get("residue", "").strip()
        parts = raw_res.split(":")
        if len(parts) == 4 and parts[0] == "R":
            chain, resname, resnum = parts[1], parts[2], parts[3]
            labels.append(f"{resname} {resnum} (Chain {chain})")
        elif len(parts) == 4 and parts[0] == "L":
            resname, resnum = parts[2], parts[3]
            labels.append(f"Ligand ({resname} {resnum})")
        else:
            clean = raw_res
            if clean.startswith("R:") or clean.startswith("L:"):
                clean = clean[2:]
            labels.append(clean.replace(":", " ").strip())

    totals = [r["total"] for r in plot_items]
    sems = [r.get("sem", 0.0) for r in plot_items]
    colors = ["#2a9d8f" if v < 0 else "#e76f51" for v in totals]

    fig, ax = plt.subplots(figsize=(8.5, max(4.5, len(labels) * 0.35)), dpi=300)
    y_pos = list(range(len(labels)))
    has_errors = any(s > 0 for s in sems)
    bars = ax.barh(
        y_pos,
        totals,
        xerr=sems if has_errors else None,
        color=colors,
        alpha=0.88,
        edgecolor="black",
        linewidth=0.6,
        capsize=3,
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontweight="bold")
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.9, alpha=0.7)
    ax.set_xlabel(r"Per-Residue $\Delta G_{\mathrm{bind}}$ Contribution (kcal/mol)", fontweight="bold")
    ax.set_title(f"MM-PBSA Per-Residue Free Energy Decomposition ({target_id} | Hotspots)", fontweight="bold", pad=12)

    for idx, bar in enumerate(bars):
        val = totals[idx]
        offset = 0.12 if val >= 0 else -0.12
        ha = "left" if val >= 0 else "right"
        ax.text(val + offset, bar.get_y() + bar.get_height() / 2, f"{val:.2f}", va="center", ha=ha, fontsize=8.5, fontweight="bold")

    out_png = working_dir / f"{prefix}decomp_mmpbsa.png"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    shutil.copy2(out_png, working_dir / "decomp_mmpbsa.png")
    return out_png


def export_analysis_csv(
    working_dir: Path, target_id: Optional[str] = None
) -> Dict[str, Path]:
    """
    Exporta todas as séries temporais e dados calculados de DM para matrizes CSV limpas
    com prefixo do alvo (<target_id>_*.csv) e espelhos (*.csv).
    """
    working_dir = Path(working_dir)
    if not target_id:
        target_id = sanitize_target_id(working_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_rmsd.xvg")) or list(working_dir.glob("*_md.tpr"))
            if candidates:
                target_id = candidates[0].stem.replace("_rmsd", "").replace("_md", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    exported_csvs: Dict[str, Path] = {}

    # 1. RMSD (Backbone e Ligante)
    rmsd_file = working_dir / f"{prefix}rmsd.xvg"
    if not rmsd_file.exists():
        rmsd_file = working_dir / "rmsd.xvg"

    rmsd_lig_file = working_dir / f"{prefix}rmsd_lig.xvg"
    if not rmsd_lig_file.exists():
        rmsd_lig_file = working_dir / "rmsd_lig.xvg"

    if rmsd_file.exists():
        x_time, y_prot, meta = parse_xvg_with_meta(rmsd_file)
        if "ps" in meta.get("xaxis_label", "").lower() or (x_time and max(x_time) > 1000 and "ns" not in meta.get("xaxis_label", "").lower()):
            x_time = [t / 1000.0 for t in x_time]

        y_lig_map = {}
        if rmsd_lig_file.exists():
            x_l, y_l, meta_l = parse_xvg_with_meta(rmsd_lig_file)
            if "ps" in meta_l.get("xaxis_label", "").lower() or (x_l and max(x_l) > 1000 and "ns" not in meta_l.get("xaxis_label", "").lower()):
                x_l = [t / 1000.0 for t in x_l]
            for xl, yl in zip(x_l, y_l):
                y_lig_map[round(xl, 3)] = yl

        rmsd_csv_path = working_dir / f"{prefix}rmsd.csv"
        with open(rmsd_csv_path, "w", encoding="utf-8") as f:
            f.write("Time_ns,Protein_Backbone_RMSD_nm,Ligand_RMSD_nm\n")
            for xt, yp in zip(x_time, y_prot):
                yl = y_lig_map.get(round(xt, 3), "")
                f.write(f"{xt:.3f},{yp:.4f},{yl if yl != '' else ''}\n")
        shutil.copy2(rmsd_csv_path, working_dir / "rmsd.csv")
        exported_csvs["rmsd"] = rmsd_csv_path

    # 2. RMSF
    rmsf_file = working_dir / f"{prefix}rmsf.xvg"
    if not rmsf_file.exists():
        rmsf_file = working_dir / "rmsf.xvg"

    if rmsf_file.exists():
        x_res, y_rmsf, _ = parse_xvg_with_meta(rmsf_file)
        rmsf_csv_path = working_dir / f"{prefix}rmsf.csv"
        with open(rmsf_csv_path, "w", encoding="utf-8") as f:
            f.write("Residue_Number,Calpha_RMSF_nm\n")
            for xr, yr in zip(x_res, y_rmsf):
                f.write(f"{int(xr)},{yr:.4f}\n")
        shutil.copy2(rmsf_csv_path, working_dir / "rmsf.csv")
        exported_csvs["rmsf"] = rmsf_csv_path

    # 3. HBond
    hbond_file = working_dir / f"{prefix}hbond.xvg"
    if not hbond_file.exists():
        hbond_file = working_dir / "hbond.xvg"

    if hbond_file.exists():
        x_time, y_hb, meta = parse_xvg_with_meta(hbond_file)
        if "ps" in meta.get("xaxis_label", "").lower() or (x_time and max(x_time) > 1000 and "ns" not in meta.get("xaxis_label", "").lower()):
            x_time = [t / 1000.0 for t in x_time]
        hbond_csv_path = working_dir / f"{prefix}hbond.csv"
        with open(hbond_csv_path, "w", encoding="utf-8") as f:
            f.write("Time_ns,HBond_Count\n")
            for xt, yb in zip(x_time, y_hb):
                f.write(f"{xt:.3f},{int(yb)}\n")
        shutil.copy2(hbond_csv_path, working_dir / "hbond.csv")
        exported_csvs["hbond"] = hbond_csv_path

    # 4. Raio de Giro (gyrate.xvg)
    gyrate_file = working_dir / f"{prefix}gyrate.xvg"
    if not gyrate_file.exists():
        gyrate_file = working_dir / "gyrate.xvg"

    if gyrate_file.exists():
        try:
            x_time, y_dict, meta = parse_xvg_multicolumn(gyrate_file)
            if "ps" in meta.get("xaxis_label", "").lower() or (x_time and max(x_time) > 1000 and "ns" not in meta.get("xaxis_label", "").lower()):
                x_time = [t / 1000.0 for t in x_time]
            headers = ["Time_ns"] + [k.replace(" ", "_") + "_nm" for k in y_dict.keys()]
            gyrate_csv_path = working_dir / f"{prefix}gyrate.csv"
            with open(gyrate_csv_path, "w", encoding="utf-8") as f:
                f.write(",".join(headers) + "\n")
                keys = list(y_dict.keys())
                for idx, xt in enumerate(x_time):
                    row_vals = [f"{xt:.3f}"] + [f"{y_dict[k][idx]:.4f}" if idx < len(y_dict[k]) else "" for k in keys]
                    f.write(",".join(row_vals) + "\n")
            shutil.copy2(gyrate_csv_path, working_dir / "gyrate.csv")
            exported_csvs["gyrate"] = gyrate_csv_path
        except Exception:
            pass

    # 5. SASA (sasa.xvg)
    sasa_file = working_dir / f"{prefix}sasa.xvg"
    if not sasa_file.exists():
        sasa_file = working_dir / "sasa.xvg"

    if sasa_file.exists():
        try:
            x_time, y_dict, meta = parse_xvg_multicolumn(sasa_file)
            if "ps" in meta.get("xaxis_label", "").lower() or (x_time and max(x_time) > 1000 and "ns" not in meta.get("xaxis_label", "").lower()):
                x_time = [t / 1000.0 for t in x_time]
            headers = ["Time_ns"] + [k.replace(" ", "_") + "_nm2" for k in y_dict.keys()]
            sasa_csv_path = working_dir / f"{prefix}sasa.csv"
            with open(sasa_csv_path, "w", encoding="utf-8") as f:
                f.write(",".join(headers) + "\n")
                keys = list(y_dict.keys())
                for idx, xt in enumerate(x_time):
                    row_vals = [f"{xt:.3f}"] + [f"{y_dict[k][idx]:.4f}" if idx < len(y_dict[k]) else "" for k in keys]
                    f.write(",".join(row_vals) + "\n")
            shutil.copy2(sasa_csv_path, working_dir / "sasa.csv")
            exported_csvs["sasa"] = sasa_csv_path
        except Exception:
            pass

    # 6. Decomposição MM-PBSA
    decomp_file = working_dir / f"{prefix}FINAL_DECOMP_MMPBSA.dat"
    if not decomp_file.exists():
        decomp_file = working_dir / "FINAL_DECOMP_MMPBSA.dat"

    if decomp_file.exists():
        decomp_data = parse_mmpbsa_decomp(decomp_file)
        if decomp_data:
            decomp_csv_path = working_dir / f"{prefix}decomp_mmpbsa.csv"
            with open(decomp_csv_path, "w", encoding="utf-8") as f:
                f.write(
                    "Residue,Van_der_Waals_kcal_mol,Electrostatic_kcal_mol,Polar_Solvation_kcal_mol,Nonpolar_Solvation_kcal_mol,Total_DeltaG_kcal_mol,Std_kcal_mol,SEM_kcal_mol\n"
                )
                for d in decomp_data:
                    raw_res = d.get("residue", "").strip()
                    f.write(
                        f'"{raw_res}",{d.get("vdw", 0.0):.3f},{d.get("electrostatic", 0.0):.3f},{d.get("polar", 0.0):.3f},{d.get("nonpolar", 0.0):.3f},{d.get("total", 0.0):.3f},{d.get("std", 0.0):.3f},{d.get("sem", 0.0):.3f}\n'
                    )
            shutil.copy2(decomp_csv_path, working_dir / "decomp_mmpbsa.csv")
            exported_csvs["decomp"] = decomp_csv_path

    # 7. Ocupação de Pontes de Hidrogênio
    hbond_occ_file = working_dir / f"{prefix}hbond_occupancy.json"
    if not hbond_occ_file.exists():
        hbond_occ_file = working_dir / "hbond_occupancy.json"

    if hbond_occ_file.exists():
        try:
            with open(hbond_occ_file, "r", encoding="utf-8") as f:
                occ_data = json.load(f)
            if occ_data:
                occ_csv_path = working_dir / f"{prefix}hbond_occupancy.csv"
                with open(occ_csv_path, "w", encoding="utf-8") as f:
                    f.write("Donor,Acceptor,Occupancy_Percent,Classification\n")
                    for row in occ_data:
                        f.write(f'"{row.get("donor", "")}","{row.get("acceptor", "")}",{row.get("occupancy_percent", 0.0):.2f},"{row.get("classification", "")}"\n')
                shutil.copy2(occ_csv_path, working_dir / "hbond_occupancy.csv")
                exported_csvs["hbond_occupancy"] = occ_csv_path
        except Exception:
            pass

    return exported_csvs

    # 7. Ocupação de Pontes de Hidrogênio
    hbond_occ_file = working_dir / "hbond_occupancy.json"
    if hbond_occ_file.exists():
        try:
            with open(hbond_occ_file, "r", encoding="utf-8") as f:
                occ_data = json.load(f)
            if occ_data:
                occ_csv_path = working_dir / "hbond_occupancy.csv"
                with open(occ_csv_path, "w", encoding="utf-8") as f:
                    f.write("Donor,Acceptor,Occupancy_Percent,Classification\n")
                    for row in occ_data:
                        f.write(f'"{row.get("donor", "")}","{row.get("acceptor", "")}",{row.get("occupancy_percent", 0.0):.2f},"{row.get("classification", "")}"\n')
                exported_csvs["hbond_occupancy"] = occ_csv_path
        except Exception:
            pass

    return exported_csvs


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


def _ensure_gmx_mmpbsa_patched() -> None:
    """
    Verifica e aplica auto-patch no pacote gmx_MMPBSA:
    1. make_top.py: Corrige o bug de indexação de pontes dissulfeto em sistemas proteicos multicadeia.
    2. utils.py: Envolve os.remove em try/except para evitar erros de FileNotFoundError (race condition em MPI).
    """
    # Patch 1: Suporte multicadeia CYS em make_top.py
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

    # Patch 2: Remoção segura de arquivos em utils.py para evitar FileNotFoundError em MPI concorrente
    try:
        import GMXMMPBSA.utils as mu

        utils_path = Path(mu.__file__)
        if utils_path.exists():
            with open(utils_path, "r", encoding="utf-8") as f:
                code_u = f.read()

            if "os.remove(fil)" in code_u:
                target_str = "            os.remove(fil)"
                safe_str = (
                    "            try:\n"
                    "                os.remove(fil)\n"
                    "            except OSError:\n"
                    "                pass"
                )
                if target_str in code_u and "try:" not in code_u.split("def remove(")[1].split("def ")[0]:
                    code_u = code_u.replace(target_str, safe_str)
                    with open(utils_path, "w", encoding="utf-8") as f:
                        f.write(code_u)
    except Exception:
        pass


def _ensure_ligand_mol2(working_dir: Path, lig_idx: int) -> Optional[Path]:
    """
    Garante a disponibilidade de um arquivo .mol2 parametrizado com GAFF/GAFF2
    para o ligante, necessário para a geração de topologias AMBER via tleap no gmx_MMPBSA.

    1. Busca por arquivos .mol2 existentes no diretório de trabalho e subdiretórios (ex: ligand_md.acpype/).
    2. Busca no diretório pai ou na árvore de dados.
    3. Se não encontrar, tenta gerar automaticamente via ACPYPE a partir de ligand_md.pdb
       ou extraindo as coordenadas do ligante via gmx trjconv.
    """
    working_dir = Path(working_dir)

    # 1. Procura em working_dir
    local_candidates = [
        p
        for p in working_dir.glob("**/*.mol2")
        if not p.name.startswith("_GMXMMPBSA_")
        and not p.parent.name.startswith("_GMXMMPBSA_")
    ]
    for pattern in ["*_bcc_gaff2.mol2", "*_AC.mol2", "*.mol2"]:
        for p in local_candidates:
            if p.match(pattern):
                return p

    # 2. Procura no diretório pai (ex: data/ ou workspace)
    try:
        parent_candidates = [
            p
            for p in working_dir.parent.glob("**/*.mol2")
            if not p.name.startswith("_GMXMMPBSA_")
            and not p.parent.name.startswith("_GMXMMPBSA_")
        ]
        for pattern in ["*_bcc_gaff2.mol2", "*_AC.mol2", "*.mol2"]:
            for p in parent_candidates:
                if p.match(pattern):
                    target_dir = working_dir / "ligand_md.acpype"
                    target_dir.mkdir(parents=True, exist_ok=True)
                    dest = target_dir / p.name
                    if not dest.exists():
                        shutil.copy2(p, dest)
                    return dest
    except Exception:
        pass

    # 3. Geração automatizada via ACPYPE
    acpype_bin = find_executable("acpype")
    if not acpype_bin:
        return None

    ligand_pdb_path = working_dir / "ligand_md.pdb"
    if not ligand_pdb_path.exists():
        # 3.1 Tenta extrair a partir de md.tpr + index.ndx usando gmx trjconv
        gmx_bin = find_executable("gmx")
        tpr_file = working_dir / "md.tpr"
        index_file = working_dir / "index.ndx"
        if gmx_bin and tpr_file.exists() and index_file.exists():
            try:
                env = os.environ.copy()
                exec_dir = str(Path(gmx_bin).parent)
                env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"
                subprocess.run(
                    [
                        gmx_bin,
                        "trjconv",
                        "-f",
                        "md.tpr",
                        "-s",
                        "md.tpr",
                        "-o",
                        "ligand_md.pdb",
                        "-n",
                        "index.ndx",
                        "-dump",
                        "0",
                    ],
                    cwd=str(working_dir),
                    input=f"{lig_idx}\n".encode(),
                    env=env,
                    capture_output=True,
                    check=True,
                )
            except Exception:
                pass

        # 3.2 Se ainda não existe, tenta converter a partir de algum arquivo de ligante (.sdf/.pdbqt/.pdb)
        if not ligand_pdb_path.exists():
            for ext in ["*.sdf", "*.mol2", "*.pdbqt"]:
                for parent_lig in working_dir.parent.glob(ext):
                    try:
                        from docking.md_prep import export_ligand_pdb

                        export_ligand_pdb(parent_lig, working_dir)
                        break
                    except Exception:
                        pass
                if ligand_pdb_path.exists():
                    break

    if ligand_pdb_path.exists():
        try:
            env = os.environ.copy()
            env.pop("PYTHONPATH", None)
            env.pop("PYTHONHOME", None)
            acpype_dir = str(Path(acpype_bin).parent)
            env["PATH"] = f"{acpype_dir}{os.pathsep}{env.get('PATH', '')}"
            subprocess.run(
                [acpype_bin, "-i", ligand_pdb_path.name, "-c", "bcc", "-f"],
                cwd=str(working_dir),
                env=env,
                capture_output=True,
                text=True,
                check=True,
            )
            acpype_folder = working_dir / f"{ligand_pdb_path.stem}.acpype"
            if acpype_folder.exists():
                for f in acpype_folder.glob("*_bcc_gaff2.mol2"):
                    return f
                for f in acpype_folder.glob("*.mol2"):
                    return f
        except Exception:
            pass

    return None


def calculate_mmpbsa(
    working_dir: Path, target_id: Optional[str] = None, force: bool = False
) -> Dict[str, Any]:
    """
    Executa o cálculo de Energia Livre de Ligação MM-PBSA via gmx_MMPBSA (Janela Termodinâmica: 60 - 100 ns / Últimos 40%):
    1. Se os resultados já foram previamente calculados e force=False, reutiliza os dados existentes sem recalcular.
    2. Realiza a limpeza em processo único de quaisquer artefatos e arquivos residuais antes da invocação do MPI.
    3. Extrai o número total de frames e configura dinamicamente o arquivo 'mmpbsa.in' para os últimos 40% (estado estacionário).
    4. Identifica os grupos do receptor (Protein) e ligante (ligand_md / LIG) em index.ndx.
    5. Executa o gmx_MMPBSA de forma paralela via MPI com mitigação de race condition e tolerância a falhas.
    6. Extrai contribuições energéticas (Van der Waals, Eletrostática, Solvatação Polar e Apolar) e Delta G total.
    7. Gera gráficos científicos (.png a 300 DPI), matrizes CSV e salva mmpbsa_summary.json.

    Retorna o dicionário com o sumário dos resultados termodinâmicos.
    """
    working_dir = Path(working_dir)
    if not working_dir.exists():
        raise FileNotFoundError(f"Diretório de trabalho não encontrado: {working_dir}")

    # Identificação do alvo
    if not target_id:
        target_id = sanitize_target_id(working_dir.name)
        if target_id.lower() in ("md_files", "screening", "data"):
            candidates = list(working_dir.glob("*_md.tpr")) or list(working_dir.glob("*_md_fit.xtc"))
            if candidates:
                target_id = candidates[0].stem.replace("_md_fit", "").replace("_md", "")
    target_id = sanitize_target_id(target_id)
    prefix = f"{target_id}_"

    tpr_file = working_dir / f"{prefix}md.tpr"
    if not tpr_file.exists():
        tpr_file = working_dir / "md.tpr"

    fit_xtc = working_dir / f"{prefix}md_fit.xtc"
    if not fit_xtc.exists():
        fit_xtc = working_dir / "md_fit.xtc"

    index_file = working_dir / f"{prefix}index.ndx"
    if not index_file.exists():
        index_file = working_dir / "index.ndx"

    dat_output = working_dir / f"{prefix}FINAL_RESULTS_MMPBSA.dat"
    if not dat_output.exists() and (working_dir / "FINAL_RESULTS_MMPBSA.dat").exists():
        dat_output = working_dir / "FINAL_RESULTS_MMPBSA.dat"

    decomp_dat_output = working_dir / f"{prefix}FINAL_DECOMP_MMPBSA.dat"
    if not decomp_dat_output.exists() and (working_dir / "FINAL_DECOMP_MMPBSA.dat").exists():
        decomp_dat_output = working_dir / "FINAL_DECOMP_MMPBSA.dat"

    if not tpr_file.exists():
        raise FileNotFoundError(f"Arquivo '{tpr_file.name}' não encontrado em: {working_dir}")
    if not fit_xtc.exists():
        raise FileNotFoundError(
            f"Trajetória corrigida '{fit_xtc.name}' não encontrada em: {working_dir}. "
            "Execute fix_pbc antes de calcular o MM-PBSA."
        )
    if not index_file.exists():
        raise FileNotFoundError(f"Arquivo '{index_file.name}' não encontrado em: {working_dir}")

    # Reutilização imediata de cálculos já concluídos com sucesso (evita reexecução desnecessária de ~1h)
    if not force and dat_output.exists() and dat_output.stat().st_size > 0:
        existing_summary = parse_mmpbsa_dat(dat_output)
        if existing_summary.get("energies", {}).get("delta_g_binding", {}).get("mean", 0.0) != 0.0:
            decomp_data = []
            if decomp_dat_output.exists() and decomp_dat_output.stat().st_size > 0:
                decomp_data = parse_mmpbsa_decomp(decomp_dat_output)
            if not decomp_data and dat_output.exists():
                decomp_data = parse_mmpbsa_decomp(dat_output)

            existing_summary["thermodynamic_window"] = "60 - 100 ns (Últimos 40% - Estado Estacionário)"
            existing_summary["protocol"] = "Dupla Escala Temporal (MM-PBSA 60-100 ns / Trajetória 0-100 ns)"
            existing_summary["raw_output_file"] = dat_output.name

            if decomp_data:
                existing_summary["per_residue_decomposition"] = decomp_data
                existing_summary["hotspot_residues"] = [r for r in decomp_data if r["total"] < 0][:10]
                plot_mmpbsa_decomp(decomp_data, working_dir, target_id=target_id)

            export_analysis_csv(working_dir, target_id=target_id)

            summary_json_path = working_dir / f"{prefix}mmpbsa_summary.json"
            with open(summary_json_path, "w", encoding="utf-8") as f:
                json.dump(existing_summary, f, indent=2, ensure_ascii=False)
            shutil.copy2(summary_json_path, working_dir / "mmpbsa_summary.json")

            return existing_summary

    mmpbsa_bin = find_executable("gmx_MMPBSA")
    if not mmpbsa_bin:
        raise DependencyError("O executável 'gmx_MMPBSA' não foi encontrado no PATH.")

    # Garante que os patches no gmx_MMPBSA (multicadeia CYS e safe unlink MPI) estejam aplicados
    _ensure_gmx_mmpbsa_patched()

    # 1. Identificação do número de frames e cálculo da janela de equilíbrio (Últimos 40%)
    total_frames = None
    for xvg_name in [f"{prefix}rmsd.xvg", "rmsd.xvg", f"{prefix}gyrate.xvg", "gyrate.xvg"]:
        xvg_file = working_dir / xvg_name
        if xvg_file.exists():
            try:
                x_vals, _ = parse_xvg(xvg_file)
                if x_vals and len(x_vals) > 0:
                    total_frames = len(x_vals)
                    break
            except Exception:
                pass

    if total_frames is None:
        gmx_bin = find_executable("gmx")
        if gmx_bin:
            try:
                env = os.environ.copy()
                env.pop("PYTHONPATH", None)
                env.pop("PYTHONHOME", None)
                exec_dir = str(Path(gmx_bin).parent)
                env["PATH"] = f"{exec_dir}{os.pathsep}{env.get('PATH', '')}"

                chk_res = subprocess.run(
                    [gmx_bin, "check", "-f", str(fit_xtc.name)],
                    cwd=str(working_dir),
                    env=env,
                    capture_output=True,
                    text=True,
                )
                chk_out = (chk_res.stderr or "") + "\n" + (chk_res.stdout or "")
                m_table = re.search(r"(?:Coords|Step)\s+(\d+)", chk_out)
                if m_table:
                    total_frames = int(m_table.group(1))
                else:
                    m_last = re.search(r"Last\s+frame\s+(\d+)", chk_out)
                    if m_last:
                        total_frames = int(m_last.group(1)) + 1
            except Exception:
                pass

    if total_frames is None or total_frames <= 0:
        total_frames = 1000

    # 2. Limpeza atômica em processo único de arquivos de execuções anteriores
    stale_patterns = [
        "_GMXMMPBSA_*",
        "*FINAL_RESULTS_MMPBSA.*",
        "*FINAL_DECOMP_MMPBSA.*",
        "*RESULTS_gmx_MMPBSA.h5",
        "COM.prmtop",
        "REC.prmtop",
        "LIG.prmtop",
        "COM_traj_*",
        "reference.pdb",
        "*.inpcrd",
        "*.mdcrd*",
    ]
    for pattern in stale_patterns:
        for stale_item in working_dir.glob(pattern):
            try:
                if stale_item.is_file():
                    stale_item.unlink(missing_ok=True)
                elif stale_item.is_dir():
                    shutil.rmtree(stale_item, ignore_errors=True)
            except Exception:
                pass

    # Protocolo de Janela Termodinâmica: Últimos 40%
    startframe = max(1, int(round(total_frames * 0.60)))
    endframe = total_frames
    frames_in_window = max(1, endframe - startframe + 1)
    interval = max(1, frames_in_window // 100)
    if frames_in_window <= 200 and interval > 2:
        interval = 2
    elif interval < 1:
        interval = 1

    # 3. Criação do arquivo de entrada mmpbsa.in
    mmpbsa_in_path = working_dir / "mmpbsa.in"
    mmpbsa_in_content = f"""&general
sys_name="{target_id}_Complex",
startframe={startframe},
endframe={endframe},
interval={interval},
verbose=2,
/
&gb
igb=5,
saltcon=0.150,
/
&decomp
idecomp=2,
dec_verbose=1,
print_res="within 6.0",
/
"""
    with open(mmpbsa_in_path, "w", encoding="utf-8") as f:
        f.write(mmpbsa_in_content)

    # 4. Identificação dos índices dos grupos no index.ndx
    groups_list = get_index_groups(index_file)
    _, _, rec_idx, lig_idx = identify_complex_groups(groups_list)

    # Identifica ou gera o arquivo mol2 parametrizado do ligante
    ligand_mol2 = _ensure_ligand_mol2(working_dir, lig_idx)

    # 5. Execução do gmx_MMPBSA
    cmd_mmpbsa = []
    mpirun_bin = find_executable("mpirun")
    if mpirun_bin:
        num_cores = max(1, min(8, (os.cpu_count() or 4) - 1))
        if num_cores > 1:
            cmd_mmpbsa.extend([mpirun_bin, "-np", str(num_cores)])

    out_dat_name = f"{prefix}FINAL_RESULTS_MMPBSA.dat"
    out_decomp_name = f"{prefix}FINAL_DECOMP_MMPBSA.dat"
    out_csv_name = f"{prefix}FINAL_RESULTS_MMPBSA.csv"

    cmd_mmpbsa.extend([
        mmpbsa_bin,
        "-O",
        "-nogui",
        "-i",
        "mmpbsa.in",
        "-cs",
        str(tpr_file.name),
        "-ct",
        str(fit_xtc.name),
        "-ci",
        str(index_file.name),
        "-cg",
        str(rec_idx),
        str(lig_idx),
        "-o",
        out_dat_name,
        "-do",
        out_decomp_name,
        "-eo",
        out_csv_name,
    ])

    if ligand_mol2 and ligand_mol2.exists():
        try:
            mol2_rel = ligand_mol2.resolve().relative_to(working_dir.resolve())
            mol2_arg = str(mol2_rel).replace("\\", "/")
        except ValueError:
            dest_dir = working_dir / "ligand_md.acpype"
            dest_dir.mkdir(parents=True, exist_ok=True)
            dest_mol2 = dest_dir / ligand_mol2.name
            if not dest_mol2.exists():
                shutil.copy2(ligand_mol2, dest_mol2)
            mol2_arg = f"ligand_md.acpype/{dest_mol2.name}"
        cmd_mmpbsa.extend(["-lm", mol2_arg])

    try:
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)
        env.pop("PYTHONHOME", None)
        exec_dir = str(Path(mmpbsa_bin).parent)
        system_paths = ["/usr/bin", "/bin", "/usr/local/bin"]
        env["PATH"] = (
            f"{exec_dir}{os.pathsep}{os.pathsep.join(system_paths)}{os.pathsep}{env.get('PATH', '')}"
        )

        result = subprocess.run(
            cmd_mmpbsa, cwd=str(working_dir), env=env, capture_output=True, text=True
        )

        final_dat = working_dir / out_dat_name
        is_successful = False
        if final_dat.exists() and final_dat.stat().st_size > 0:
            parsed_test = parse_mmpbsa_dat(final_dat)
            if parsed_test.get("energies", {}).get("delta_g_binding", {}).get("mean", 0.0) != 0.0:
                is_successful = True

        if not is_successful and result.returncode != 0:
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

    # Cria espelhos
    final_dat = working_dir / out_dat_name
    final_decomp = working_dir / out_decomp_name
    if final_dat.exists():
        shutil.copy2(final_dat, working_dir / "FINAL_RESULTS_MMPBSA.dat")
    if final_decomp.exists():
        shutil.copy2(final_decomp, working_dir / "FINAL_DECOMP_MMPBSA.dat")

    # 6. Parse dos resultados globais e decomposição por resíduo
    summary_data = parse_mmpbsa_dat(final_dat)
    summary_data["thermodynamic_window"] = "60 - 100 ns (Últimos 40% - Estado Estacionário)"
    summary_data["startframe"] = startframe
    summary_data["endframe"] = endframe
    summary_data["interval"] = interval
    summary_data["total_frames_estimated"] = total_frames
    summary_data["frames_analyzed"] = max(1, (endframe - startframe + 1) // interval)
    summary_data["protocol"] = "Dupla Escala Temporal (MM-PBSA 60-100 ns / Trajetória 0-100 ns)"
    summary_data["raw_output_file"] = out_dat_name

    decomp_data = []
    if final_decomp.exists() and final_decomp.stat().st_size > 0:
        decomp_data = parse_mmpbsa_decomp(final_decomp)
    if not decomp_data and final_dat.exists():
        decomp_data = parse_mmpbsa_decomp(final_dat)

    if decomp_data:
        summary_data["per_residue_decomposition"] = decomp_data
        summary_data["hotspot_residues"] = [r for r in decomp_data if r["total"] < 0][:10]
        plot_mmpbsa_decomp(decomp_data, working_dir, target_id=target_id)

    export_analysis_csv(working_dir, target_id=target_id)

    summary_json_path = working_dir / f"{prefix}mmpbsa_summary.json"
    with open(summary_json_path, "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=2, ensure_ascii=False)
    shutil.copy2(summary_json_path, working_dir / "mmpbsa_summary.json")

    return summary_data
