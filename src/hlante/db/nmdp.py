"""
hlante.db.nmdp
=============

NMDP (National Marrow Donor Program / Be The Match) allele-frequency
client. Supplements AFND with NMDP registry-derived frequencies when
AFND data is unavailable.

NMDP frequency data are **not** redistributed with HLAnte: the
resource is licensed by NMDP/Be The Match and its terms do not permit
redistribution. This client is therefore inert by default. To enable
it, obtain an extract from https://frequency.nmdp.org/ under the terms
of that resource and save it as ``~/.hlante/nmdp/nmdp_frequencies.tsv``;
HLAnte will then use it as a secondary frequency source.

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

        NMDP frequency data are **not** redistributed with HLAnte. The
        client stays inert unless the user supplies their own extract
        (see :data:`NMDP_DEFAULT_DIR`). Never raises — NMDP is an
        optional secondary source.
        """
        if self._loaded:
            return
        path = self._locate_tsv()
        if path is None:
            logger.info(
                "NMDP: no local frequency table found; skipping this source. "
                "NMDP data are not redistributed with HLAnte. To enable it, "
                "obtain an extract from https://frequency.nmdp.org/ under the "
                "terms of that resource and save it as %s/%s",
                NMDP_DEFAULT_DIR,
                NMDP_TSV_FILENAME,
            )
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
]
