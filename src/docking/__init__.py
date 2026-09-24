import logging
from pathlib import Path

# Garante que o logger do RDKit possua ao menos um handler para evitar IndexError no Meeko
logging.getLogger("rdkit").addHandler(logging.NullHandler())

# Garante que subprocessos Python (como os scripts CLI do Meeko: mk_export, mk_prepare_receptor, etc.)
# também possuam o NullHandler através de sitecustomize.py no site-packages ativo
try:
    import site

    for sp in site.getsitepackages():
        sp_path = Path(sp)
        if sp_path.exists() and sp_path.is_dir():
            sc = sp_path / "sitecustomize.py"
            if not sc.exists():
                sc.write_text(
                    "# Auto-generated to fix Meeko IndexError: list index out of range\n"
                    "import logging\n"
                    "logging.getLogger('rdkit').addHandler(logging.NullHandler())\n",
                    encoding="utf-8",
                )
except Exception:
    pass

from . import (  # noqa: E402
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

__all__ = [
    "analysis",
    "box_utils",
    "md_analysis",
    "md_equil",
    "md_prep",
    "notifier",
    "pharmacokinetics",
    "preparation",
    "report",
    "vina_runner",
    "visualization",
]
