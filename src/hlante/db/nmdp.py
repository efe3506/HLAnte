"""
hlante.db.nmdp
=============

NMDP (National Marrow Donor Program / Be The Match) allele-frequency
client. Supplements AFND with NMDP registry-derived frequencies when
AFND data is unavailable.

A built-in frequency table for ~280 common alleles across four NMDP
ethnic groups (European, African American, Asian/Pacific Islander,
Hispanic) is bundled with the package. For full registry coverage,
download the official NMDP tables from https://frequency.nmdp.org/
and save them as ``~/.hlante/nmdp/nmdp_frequencies.tsv``.

TSV format expected (same as AFND):
- ``Allele``          — e.g. ``B*57:01``
- ``Population``      — source study population label
- ``Population Group`` — geographic group (``European``, ``African``, ...)
- ``Frequency``       — 0.0–1.0
- ``Sample Size``     — integer
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from hlante.db.afnd import (
    AFNDClient,
    AFNDDatabaseError,
    DEFAULT_MIN_SAMPLE_SIZE,
    DEFAULT_POPULATION_GROUP,
)

logger: logging.Logger = logging.getLogger(__name__)


NMDP_DEFAULT_DIR: Path = Path.home() / ".hlante" / "nmdp"
NMDP_TSV_FILENAME: str = "nmdp_frequencies.tsv"
BUILTIN_NMDP_TSV: Path = Path(__file__).parent / "nmdp_builtin.tsv"


class NMDPClient(AFNDClient):
    """
    NMDP allele-frequency client.

    Drop-in supplement for :class:`~hlante.db.afnd.AFNDClient`. Uses
    NMDP registry-derived frequencies and tags all results with
    ``source="NMDP"``.

    Parameters
    ----------
    local_dir : Path, optional
        Directory holding the NMDP TSV. Defaults to
        ``~/.hlante/nmdp/``.
    population_group : str, optional
        Default population group (same codes as AFND).
    min_sample_size : int, optional
        Minimum per-study sample size (default 50).
    """

    def __init__(
        self,
        local_dir: Optional[Path] = None,
        *,
        population_group: str = DEFAULT_POPULATION_GROUP,
        min_sample_size: int = DEFAULT_MIN_SAMPLE_SIZE,
    ) -> None:
        super().__init__(
            local_dir or NMDP_DEFAULT_DIR,
            population_group=population_group,
            min_sample_size=min_sample_size,
            source_label="NMDP",
        )

    def load(self) -> None:
        """
        Load the NMDP frequency TSV (idempotent).

        Falls back to the bundled built-in table when no user-supplied
        TSV is found. Never raises — NMDP is a secondary source.
        """
        if self._loaded:
            return
        path = self._locate_tsv()
        if path is None:
            if BUILTIN_NMDP_TSV.is_file():
                path = BUILTIN_NMDP_TSV
                logger.info(
                    "NMDP: no local TSV found; using built-in data (%s). "
                    "For full registry coverage download from "
                    "https://frequency.nmdp.org/ and save as "
                    "%s/%s",
                    path,
                    NMDP_DEFAULT_DIR,
                    NMDP_TSV_FILENAME,
                )
            else:
                logger.debug("NMDP: built-in data not found; skipping.")
                self._rows = []
                self._loaded = True
                return
        try:
            self._rows = self._parse_tsv(path)
        except AFNDDatabaseError as exc:
            logger.warning("NMDP: failed to load TSV (%s): %s", path, exc)
            self._rows = []
        self._reindex()
        self._loaded = True
        logger.info("NMDP loaded: %d frequency row(s) (%s)", len(self._rows), path)

    def _locate_tsv(self) -> Optional[Path]:
        candidates = [
            self.local_dir / NMDP_TSV_FILENAME,
            self.local_dir / "nmdp.tsv",
        ]
        for p in candidates:
            if p.is_file():
                return p
        return None


__all__ = [
    "NMDPClient",
    "NMDP_DEFAULT_DIR",
    "NMDP_TSV_FILENAME",
    "BUILTIN_NMDP_TSV",
]
