import logging

# Garante que o logger do RDKit possua ao menos um handler para evitar IndexError no Meeko
logging.getLogger("rdkit").addHandler(logging.NullHandler())

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
