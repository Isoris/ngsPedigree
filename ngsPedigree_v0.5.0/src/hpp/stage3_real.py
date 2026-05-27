"""
Stage 3 real adapter — swap target for ``stage3_placeholder``.

Activated when ngsPedigree Stage 3 (``STEP_PED_03_inheritance_map.py``)
ships. The function signatures and return types must remain identical
to ``stage3_placeholder`` so downstream HPP code does not change. The
adapter's job is to read whatever schema Stage 3 actually emits and
return ``DyadSegment`` / ``TriadSegment`` / ``parent_phase`` records.

HANDOFF.md open question #1 is the main risk: where does
parent-heterozygous-site phase live in Stage 3's output? Until that is
locked, this module raises ``Stage3RealNotReadyError``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from .stage3_placeholder import (
    DyadSegment,
    PhaseKey,
    Stage3Inputs,
    TriadSegment,
)


class Stage3RealNotReadyError(NotImplementedError):
    """Stage 3 real adapter is not yet implemented. See HANDOFF.md."""


def load_dyad_map(path: str | Path) -> List[DyadSegment]:
    raise Stage3RealNotReadyError(
        "stage3_real.load_dyad_map: ngsPedigree Stage 3 has not shipped yet. "
        "Use stage3_placeholder.load_dyad_map until it does."
    )


def load_triad_map(path: str | Path) -> List[TriadSegment]:
    raise Stage3RealNotReadyError(
        "stage3_real.load_triad_map: ngsPedigree Stage 3 has not shipped yet."
    )


def load_parent_phase(path: str | Path) -> Dict[PhaseKey, str]:
    raise Stage3RealNotReadyError(
        "stage3_real.load_parent_phase: open question — where does "
        "parent-het phase live in Stage 3's output? See HANDOFF.md #1."
    )


def load(
    dyad_map_path: Optional[str | Path] = None,
    triad_map_path: Optional[str | Path] = None,
    parent_phase_path: Optional[str | Path] = None,
) -> Stage3Inputs:
    raise Stage3RealNotReadyError(
        "stage3_real.load: ngsPedigree Stage 3 not ready."
    )
