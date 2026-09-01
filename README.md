# Manual do Pipeline Integrado de Docking Molecular e Dinâmica Molecular (GROMACS)

Este projeto automatiza de ponta a ponta o fluxo de trabalho de Docking Molecular, Triagem Virtual, Dinâmica Molecular (GROMACS) e Termodinâmica MM-PBSA.

---

## 🔒 Arquitetura de Segurança e Isolamento Estrito (Target Isolation)

Para eliminar o risco de contaminação cruzada entre diferentes alvos biológicos e garantir reprodutibilidade estrita:

1. **Hierarquia Obrigatória por Alvo:**
   - **Triagem Virtual:** `data/screening/<PDB_ID>/<LIGAND_NAME>/`
   - **Dinâmica Molecular:** `data/md_files/<PDB_ID>/`
   - **Exportação para Cluster Remoto:** `cluster_export/<PDB_ID>/`
2. **Nomenclatura Explícita:** Todos os arquivos de simulação e análise possuem o prefixo do alvo (ex: `<PDB_ID>_complex.gro`, `<PDB_ID>_em.tpr`, `<PDB_ID>_nvt.tpr`, `<PDB_ID>_md.tpr`, `<PDB_ID>_md_fit.xtc`, `<PDB_ID>_mmpbsa_summary.json`).
3. **Validação Fail-Fast:** Antes de executar etapas pesadas, o pipeline valida a presença da proteína e do ligante ($\ge 10$ átomos pesados) no complexo, e executa a verificação de integridade binária do TPR (`gmx dump -s`) contra arquivos corrompidos.
4. **Purga Segura de Resíduos:** Detecta automaticamente arquivos residuais (`#*#`, `*.cpt`, `*.tpr`) antes de inicializar novas etapas para evitar leitura de checkpoints obsoletos.

---

## 🚀 Comandos Disponíveis (CLI)

O pipeline pode ser executado via interface interativa TUI ou por comandos individuais no terminal:

### 0. Modo Interativo (Recomendado)
```bash
uv run src/main.py interactive
```

### 1. Triagem Virtual (Virtual Screening)
```bash
uv run src/main.py screen --receptor data/7CFN/processed/7CFN_receptor.pdbqt --ligand ligante.pdbqt --target 7CFN --cx 12.5 --cy 8.2 --cz -15.4
```

### 2. Preparação e Minimização de DM (`md-prep`)
Prepara o complexo receptor-ligante, parametriza com ACPYPE/GAFF2, solvata e minimiza a energia:
```bash
uv run src/main.py md-prep --receptor data/7CFN/processed/7CFN_receptor.pdb --sdf data/screening/7CFN/Desoxicolato/docked_poses.sdf --target 7CFN
```

### 3. Equilíbrio Termodinâmico NVT / NPT (`md-equil`)
Executa o equilíbrio com restrição de posição nos átomos pesados:
```bash
uv run src/main.py md-equil --dir data/md_files/7CFN --target 7CFN
```

### 4. Compilação de Produção & Pacote para Cluster (`md-compile`)
Compila `<target_id>_md.tpr`, valida a consistência e gera o pacote modular em `cluster_export/<target_id>/`:
```bash
uv run src/main.py md-compile --dir data/md_files/7CFN --target 7CFN
```

### 5. Exportação Modular para Cluster (`md-export`)
Empacota a simulação para execução via SSH em servidores remotos / tmux sem depender de Slurm:
```bash
uv run src/main.py md-export --dir data/md_files/7CFN --target 7CFN
```

### 6. Produção Local de DM (`md-run`)
Executa a dinâmica de produção localmente (auto-detecta GPU/CPU):
```bash
uv run src/main.py md-run --dir data/md_files/7CFN --target 7CFN
```

### 7. Pós-Processamento, Gráficos e MM-PBSA (`md-postprocess`)
Ajusta condições periódicas de contorno (PBC), calcula RMSD/RMSF/HBond/Rg/SASA e executa o cálculo de energia livre MM-PBSA (60 - 100 ns):
```bash
uv run src/main.py md-postprocess --dir data/md_files/7CFN --target 7CFN
```

### 8. Relatório Executivo HTML e Visualização 3D (`report`)
Gera o relatório consolidado de publicação e script PyMOL 3D:
```bash
uv run src/main.py report --dir data/md_files/7CFN --receptor 7CFN --ligand Desoxicolato
```

---

## 🖥️ Execução em Servidor/Cluster Remoto (SSH / tmux)

O pacote exportado em `cluster_export/<PDB_ID>/` é 100% autossuficiente:

```bash
# 1. Enviar para o servidor
rsync -avP cluster_export/7CFN/ user@cluster:/path/to/simulations/7CFN/

# 2. Conectar e abrir sessão tmux
ssh user@cluster
tmux new -s md_7CFN
cd /path/to/simulations/7CFN

# 3. Executar o script inteligente (auto-detecta GPU e suporta retomada -cpi)
./run_local.sh

# 4. Desanexar do tmux: Ctrl+B, depois D
# 5. Acompanhar em tempo real:
tail -f 7CFN_md.log
```

---

## 📦 Estrutura de Diretórios Gerada

```text
molecular_docking/
├── cluster_export/
│   └── 7CFN/
│       ├── 7CFN_md.tpr
│       ├── run_local.sh
│       └── README.md
├── data/
│   ├── 7CFN/
│   │   ├── raw/
│   │   └── processed/
│   ├── screening/
│   │   └── 7CFN/
│   │       └── Desoxicolato/
│   └── md_files/
│       └── 7CFN/
│           ├── 7CFN_complex.gro
│           ├── 7CFN_topol.top
│           ├── 7CFN_em.tpr / 7CFN_em.gro
│           ├── 7CFN_nvt.tpr / 7CFN_nvt.gro
│           ├── 7CFN_npt.tpr / 7CFN_npt.gro
│           ├── 7CFN_md.tpr / 7CFN_md.xtc
│           ├── 7CFN_md_fit.xtc / 7CFN_md_clean.gro
│           ├── 7CFN_rmsd.png / 7CFN_rmsf.png / 7CFN_hbond.png
│           ├── 7CFN_FINAL_RESULTS_MMPBSA.dat
│           └── 7CFN_mmpbsa_summary.json
```
