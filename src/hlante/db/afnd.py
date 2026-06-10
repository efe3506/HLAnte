"""
hlante.db.afnd
=============

Allele Frequency Net Database (AFND) client.

Update mechanism
----------------
``AFNDClient.update()`` downloads the machine-readable AFND mirror
published by Slowikowski (2024) at :data:`SLOWIKOWSKI_URL`
(``https://raw.githubusercontent.com/slowkow/allelefrequencies/main/afnd.tsv``).
The downloaded file uses Slowikowski's 7-column schema; this module
automatically transforms it into HLAnte's internal 5-column TSV and
saves it to ``~/.hlante/afnd/afnd_frequencies.tsv``.

When ``update()`` has not been run, HLAnte falls back to a compact
built-in frequency table (7 loci × 5 populations) bundled with the
package, derived from Gonzalez-Galarza et al. (2020) and
Pappas et al. (2016). Running ``hlante db-update --db afnd`` replaces
this with the full Slowikowski dataset (~6 MB; >3,000 population
studies, 8 HLA loci).

Expected TSV schema (internal, post-transform)
-----------------------------------------------
Column names are case-insensitive; flexible lookup is applied:

- ``Allele`` — e.g. ``B*57:01`` (``HLA-`` prefix optional)
- ``Population`` — source study population label
- ``Population Group`` — geographic group assigned during transform
  (``"European"``, ``"African"``, ``"East Asian"``, ``"American"``,
  ``"Oceanian"``, ``"South Asian"``, ``"Middle Eastern"``, ``"Unknown"``)
- ``Frequency`` — allele frequency in the 0.0–1.0 range
- ``Sample Size`` — integer

Population selection
--------------------
The user passes a population-group code (``EUR``, ``AFR``, ``EAS``,
``SAS``, ``MID``, ``AMR``, ``OCE``, or ``global``; ``ASN`` is accepted
as a backward-compatible alias for ``EAS``). The code maps to a
list of geographic keywords in :data:`POPULATION_GROUPS`. A row is
considered a match when any keyword is a case-insensitive substring of
the row's population label. No country is special-cased.
"""

from __future__ import annotations

import csv
import logging
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from hlante.db.gwas import (
    _colon_group_count,
    _field_label,
    _get_resolution_levels,
    _truncate_to_fields,
)

logger: logging.Logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_LOCAL_DIR: Path = Path.home() / ".hlante" / "afnd"
DEFAULT_TSV_FILENAME: str = "afnd_frequencies.tsv"
DEFAULT_POPULATION_GROUP: str = "global"
DEFAULT_MIN_SAMPLE_SIZE: int = 50

#: Primary AFND data source — machine-readable mirror maintained by
#: Slowikowski (2024), scraped from allelefrequencies.net.
#: Schema: group / gene / allele / population / indivs_over_n /
#:          alleles_over_2n / n
SLOWIKOWSKI_URL: str = (
    "https://raw.githubusercontent.com/slowkow/allelefrequencies/main/afnd.tsv"
)

# Built-in fallback: literature-derived 7-locus × 5-population frequency table
# bundled with the package so that AFND works offline without a manual download.
# Source: AFND (Gonzalez-Galarza et al. 2020, NAR) and Pappas et al. 2016.
# Replaced by the full Slowikowski dataset after `hlante db-update --db afnd`.
BUILTIN_AFND_TSV: Path = Path(__file__).parent / "afnd_builtin.tsv"

# ---------------------------------------------------------------------------
# Geographic classification
# ---------------------------------------------------------------------------

#: Priority-ordered rules for mapping a free-text population label to a
#: :data:`POPULATION_GROUPS` key.  Earlier entries take precedence; ``EUR``
#: is last because its keyword list is broad.
_GEO_CLASSIFICATION_RULES: List[Tuple[str, List[str]]] = [
    ("AFR", [
        "africa", "african", "nigeria", "nigerian", "ghana", "ghanaian",
        "senegal", "kenyan", "kenya", "ethiopia", "cameroon", "mali",
        "burkina", "ivory coast", "zimbabwe", "mozambique", "tanzania",
        "cape verde", "sub-saharan", "black south",
    ]),
    # Keyed "EAS" to match the 1000 Genomes super-population convention used
    # throughout the manuscript and benchmark. "ASN" is accepted as a
    # backward-compatible alias (see :data:`_GROUP_ALIASES`).
    ("EAS", [
        "japan", "japanese", "china", "chinese", "korea", "korean",
        "taiwan", "taiwanese", "vietnam", "thai", "thailand", "myanmar",
        "cambodia", "laos", "singapore", "malaysia", "indonesian",
        "indonesia", "philippines", "filipino", "east asian",
    ]),
    ("SAS", [
        "india", "indian", "pakistan", "pakistani", "bangladesh",
        "sri lanka", "nepal", "south asian", "gujarati", "punjabi",
        "dravidian", "sindhi",
    ]),
    ("MID", [
        "iran", "iranian", "iraq", "iraqi", "arab", "arabic", "saudi",
        "jordan", "jordanian", "lebanon", "lebanese", "syria", "syrian",
        "turkey", "turkish", "azerbaija", "middle east", "north africa",
        "maghreb", "morocco", "moroccan", "algeria", "algerian",
        "tunisia", "tunisian", "libya", "libyan", "egypt", "egyptian",
        "persian", "kuwaiti", "omani", "yemeni", "bahrain", "qatar",
    ]),
    ("AMR", [
        "american", "hispanic", "latino", "mestizo", "mexico", "mexican",
        "brazil", "brazilian", "argentina", "argentinean", "colombia",
        "colombian", "chile", "chilean", "peru", "peruvian", "venezuela",
        "ecuadorian", "ecuador", "bolivia", "paraguay", "uruguay",
        "costa rica", "cuba", "cuban", "puerto rico", "dominican",
        "caribbean", "native american", "amerindian", "indigenous",
        "guatemal", "hondur", "nicaragu",
    ]),
    ("OCE", [
        "oceania", "oceanian", "pacific", "australia", "australian",
        "new zealand", "papua", "fiji", "samoa", "polynesia", "melanesia",
        "micronesia", "maori", "aboriginal",
    ]),
    ("EUR", [
        "european", "caucasian", "white", "spain", "spanish", "france",
        "french", "germany", "german", "italy", "italian", "greece",
        "greek", "sweden", "swedish", "norway", "norwegian", "denmark",
        "danish", "finland", "finnish", "netherlands", "dutch", "belgium",
        "austrian", "austria", "switzerland", "swiss", "poland", "polish",
        "czech", "hungary", "hungarian", "romania", "romanian", "bulgaria",
        "bulgarian", "croatia", "serbian", "serbia", "ukraine", "ukrainian",
        "russia", "russian", "portugal", "portuguese", "ireland", "irish",
        "scotland", "scottish", "england", "english", "british", "nordic",
        "scandinavian", "basque", "catalan", "iberian", "galician",
        "sardinian", "sicilian", "maltese", "georgian", "armenian",
        "slovak", "latvian", "lithuanian", "estonian",
    ]),
]


#: Accepted aliases for population-group codes. The canonical key for the
#: East/Southeast Asian group is "EAS" (matching the 1000 Genomes
#: super-population convention used in the manuscript and benchmark). "ASN" is
#: accepted as a backward-compatible alias so earlier ``--population`` scripts
#: keep working.
_GROUP_ALIASES: Dict[str, str] = {"ASN": "EAS"}


def _canonical_group(target_group: str) -> str:
    """Resolve a population-group code through :data:`_GROUP_ALIASES`."""
    return _GROUP_ALIASES.get(target_group.strip().upper(), target_group)


def _classify_population_group(label: str) -> str:
    """
    Map a free-text AFND population label to a :data:`POPULATION_GROUPS`
    key using :data:`_GEO_CLASSIFICATION_RULES`.

    Returns ``"Unknown"`` when no rule matches.
    """
    lower = label.lower()
    for group, keywords in _GEO_CLASSIFICATION_RULES:
        if any(kw in lower for kw in keywords):
            return group
    return "Unknown"


def _transform_slowikowski_tsv(src: Path, dest: Path) -> int:
    """
    Convert Slowikowski's ``afnd.tsv`` schema to HLAnte's internal format.

    Input columns (Slowikowski):
        ``group / gene / allele / population / indivs_over_n /
        alleles_over_2n / n``

    Output columns (HLAnte):
        ``Allele / Population / Population Group / Frequency / Sample Size``

    Only rows with ``group == "hla"`` are kept.  Returns the number of
    rows written.
    """
    rows_written = 0
    with (
        src.open("r", encoding="utf-8", newline="") as sf,
        dest.open("w", encoding="utf-8", newline="") as df,
    ):
        reader = csv.DictReader(sf, delimiter="\t")
        writer = csv.writer(df, delimiter="\t")
        writer.writerow(
            ["Allele", "Population", "Population Group", "Frequency", "Sample Size"]
        )
        for row in reader:
            if (row.get("group") or "").strip().lower() != "hla":
                continue
            allele = (row.get("allele") or "").strip()
            if not allele or "*" not in allele:
                continue
            freq_str = (row.get("alleles_over_2n") or "").strip()
            if not freq_str:
                continue
            try:
                freq = float(freq_str)
            except ValueError:
                continue
            if not 0.0 <= freq <= 1.0:
                continue
            try:
                n = int(float(row.get("n") or 0))
            except ValueError:
                n = 0
            population = (row.get("population") or "").strip()
            writer.writerow(
                [allele, population, _classify_population_group(population), freq, n]
            )
            rows_written += 1
    return rows_written

#: Universal population-group taxonomy used for AFND filtering.
#:
#: No country-level entries appear in this mapping — AFND population
#: rows that carry a country name in their ``Population`` column are
#: matched via substring against the geographic keywords below (for
#: example a row with ``Population Group = "East Asian"`` matches
#: ``EAS`` regardless of the country label). Keeping the taxonomy
#: purely geographic avoids privileging any particular country or
#: ethnic group.
POPULATION_GROUPS: Dict[str, List[str]] = {
    "EUR": [
        "European",
        "Caucasian",
        "White",
        "Western European",
        "Northern European",
        "Southern European",
        "Eastern European",
    ],
    "AFR": [
        "African",
        "Sub-Saharan African",
        "West African",
        "East African",
        "Black",
    ],
    "EAS": [
        "Asian",
        "East Asian",
        "Southeast Asian",
    ],
    "SAS": [
        "South Asian",
    ],
    "MID": [
        "Middle Eastern",
        "Arab",
        "North African",
        "West Asian",
    ],
    "AMR": [
        "American",
        "Hispanic",
        "Latino",
        "Mestizo",
    ],
    "OCE": [
        "Oceanian",
        "Pacific Islander",
    ],
    "global": [],  # No filter — aggregate everything (weighted by sample size).
}


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class AFNDDatabaseError(Exception):
    """
    Raised for AFND-related errors.
    """

    pass


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass
class AllelFrequency:
    """
    A single aggregated allele × population-group frequency record.

    Attributes
    ----------
    allele : str
        Allele in IPD-IMGT/HLA canonical form (no ``HLA-`` prefix).
    frequency : float
        Allele frequency in the 0.0–1.0 range.
    population_group : str
        Population-group code selected at query time (e.g., ``"EUR"``
        or ``"global"``).
    sample_size : int
        Aggregate sample size behind the frequency estimate.
    source_resolution : int
        Colon-group count at which the match was found (1–4).
    is_estimated : bool
        ``True`` when fallback resolution had to be used.
    populations_aggregated : int
        Number of source rows combined into this record.
    """

    allele: str
    frequency: float
    population_group: str
    sample_size: int
    source_resolution: int
    is_estimated: bool = False
    populations_aggregated: int = 1
    source: str = "AFND"

    @property
    def resolution_label(self) -> str:
        """
        The reporter-compatible ``"N-field"`` label.
        """
        return _field_label(self.source_resolution)

    # Backward-compatibility alias: legacy callers read ``.population``.
    @property
    def population(self) -> str:
        return self.population_group


# ---------------------------------------------------------------------------
# AFND client
# ---------------------------------------------------------------------------


class AFNDClient:
    """
    AFND local-TSV client.

    Parameters
    ----------
    local_dir : Path, optional
        Directory holding the TSV. Defaults to :data:`DEFAULT_LOCAL_DIR`.
    population_group : str, optional
        Default population group used in lookups. One of
        :data:`POPULATION_GROUPS` keys (``EUR``, ``AFR``, ``EAS``,
        ``SAS``, ``MID``, ``AMR``, ``OCE``, ``global``; ``ASN`` is an
        accepted alias for ``EAS``). Custom strings
        are passed through and interpreted as a case-insensitive
        substring filter against the row population label.
    min_sample_size : int, optional
        Drop rows below this sample size (default ``50``).
    """

    def __init__(
        self,
        local_dir: Optional[Path] = None,
        *,
        population_group: str = DEFAULT_POPULATION_GROUP,
        min_sample_size: int = DEFAULT_MIN_SAMPLE_SIZE,
        source_label: str = "AFND",
        # Backward-compatibility alias — previous releases used
        # ``population``; silently accept it.
        population: Optional[str] = None,
    ) -> None:
        if population is not None and population_group == DEFAULT_POPULATION_GROUP:
            population_group = population
        self.local_dir: Path = Path(local_dir) if local_dir else DEFAULT_LOCAL_DIR
        self.population_group: str = population_group
        self.min_sample_size: int = min_sample_size
        self._source_label: str = source_label
        self._rows: List[Dict[str, str]] = []
        self._by_allele: Dict[str, List[Dict[str, str]]] = {}
        self._loaded: bool = False

    # ---- Public API ----

    def update(self, source_url: Optional[str] = None) -> Path:
        """
        Download and install the AFND frequency TSV.

        By default, fetches the machine-readable AFND mirror maintained
        by Slowikowski (2024) from :data:`SLOWIKOWSKI_URL`.  If the
        downloaded file uses Slowikowski's 7-column schema (detected by
        the presence of ``alleles_over_2n`` in the header), it is
        automatically transformed to HLAnte's 5-column internal format
        before saving.  A user-supplied URL pointing to a file already in
        HLAnte's format is accepted and copied as-is.

        Parameters
        ----------
        source_url : str, optional
            Override the default :data:`SLOWIKOWSKI_URL`.

        Returns
        -------
        Path
            Path to the saved ``afnd_frequencies.tsv``.

        Raises
        ------
        AFNDDatabaseError
            When the download or transform fails.
        """
        self.local_dir.mkdir(parents=True, exist_ok=True)
        url = source_url or SLOWIKOWSKI_URL
        dest = self.local_dir / DEFAULT_TSV_FILENAME
        tmp = self.local_dir / "_afnd_download.tmp"

        logger.info("Downloading AFND data: %s", url)
        try:
            with urllib.request.urlopen(url, timeout=180) as resp:  # noqa: S310
                tmp.write_bytes(resp.read())
        except Exception as exc:  # noqa: BLE001
            if tmp.is_file():
                tmp.unlink()
            raise AFNDDatabaseError(f"AFND download failed ({url}): {exc}") from exc

        # Detect Slowikowski format by inspecting the first header line.
        try:
            with tmp.open("r", encoding="utf-8") as tf:
                first_line = tf.readline()
        except OSError as exc:
            tmp.unlink(missing_ok=True)
            raise AFNDDatabaseError(f"AFND temp file unreadable: {exc}") from exc

        if "alleles_over_2n" in first_line:
            logger.info("Detected Slowikowski format — transforming to HLAnte schema")
            try:
                n_rows = _transform_slowikowski_tsv(tmp, dest)
            except Exception as exc:  # noqa: BLE001
                tmp.unlink(missing_ok=True)
                raise AFNDDatabaseError(
                    f"AFND transform failed: {exc}"
                ) from exc
            tmp.unlink(missing_ok=True)
            logger.info(
                "AFND: %d HLA frequency rows saved → %s", n_rows, dest
            )
        else:
            # Already in HLAnte format; replace in place atomically.
            tmp.replace(dest)
            logger.info("AFND: custom TSV saved → %s", dest)

        self._loaded = False
        return dest

    def load(self) -> None:
        """
        Load the local TSV into memory (idempotent).

        Falls back to the built-in literature-derived frequency table
        (:data:`BUILTIN_AFND_TSV`) when no user-supplied TSV is found,
        so that population-frequency scoring works out-of-the-box
        without a manual AFND download.
        """
        if self._loaded:
            return
        path = self._locate_tsv()
        if path is None:
            if BUILTIN_AFND_TSV.is_file():
                path = BUILTIN_AFND_TSV
                logger.info(
                    "AFND: no local TSV found; using built-in fallback (%s). "
                    "For full population coverage download from "
                    "http://www.allelefrequencies.net/datasets.asp",
                    path,
                )
            else:
                raise AFNDDatabaseError(
                    f"AFND TSV not found in {self.local_dir} and built-in "
                    "fallback is missing. Download the file manually or call update()."
                )
        self._rows = self._parse_tsv(path)
        self._reindex()
        self._loaded = True
        logger.info("AFND loaded: %d frequency row(s) (%s)", len(self._rows), path)

    def get_frequency(
        self,
        allele: str,
        population_group: Optional[str] = None,
    ) -> Optional[AllelFrequency]:
        """
        Aggregate a single frequency for an allele × population group.

        Parameters
        ----------
        allele : str
            HLA allele (``HLA-`` prefix optional).
        population_group : str, optional
            Population-group code. When ``None`` the instance default
            is used.

        Returns
        -------
        AllelFrequency or None
            ``None`` when no row matches the filter.
        """
        try:
            self.load()
        except AFNDDatabaseError as exc:
            logger.debug("AFND could not be loaded: %s", exc)
            return None

        group = population_group if population_group is not None else self.population_group
        bare = _strip_hla(allele)
        resolution = _colon_group_count(bare)
        matched = self._match_rows(bare, group)
        if not matched:
            return None
        return self._aggregate(
            matched,
            allele=bare,
            source_resolution=resolution,
            population_group_label=group,
            is_estimated=False,
        )

    def get_frequency_with_fallback(
        self,
        allele: str,
        population_group: Optional[str] = None,
        min_resolution: int = 2,
    ) -> Optional[AllelFrequency]:
        """
        Same as :meth:`get_frequency` but with resolution fallback.

        Fallback order mirrors
        :func:`hlante.db.gwas._get_resolution_levels` (from the current
        colon-group count down to 1). The first matching level is
        returned. When the matching level is coarser than the original
        input, ``is_estimated`` is set to ``True``.

        Parameters
        ----------
        allele : str
            Allele to look up.
        population_group : str, optional
            Population-group filter.
        min_resolution : int, optional
            Minimum digit-based level to descend to. Default ``2``
            (allows descent to ``A*02``).

        Returns
        -------
        AllelFrequency or None
        """
        bare = _strip_hla(allele)
        original_resolution = _colon_group_count(bare)
        group = population_group if population_group is not None else self.population_group

        # Translate digit-based (2/4/6/8) → colon-group floor.
        min_colon = max(1, min_resolution // 2)
        for n_colon in _get_resolution_levels(bare):
            if n_colon < min_colon:
                break
            truncated = _truncate_to_fields(bare, n_colon)
            freq = self.get_frequency(truncated, population_group=group)
            if freq is not None:
                freq.is_estimated = n_colon != original_resolution
                return freq
        return None

    # ---- Matching logic ----

    def _matches_group(self, row_population: str, target_group: str) -> bool:
        """
        Return ``True`` when a population row belongs to the target group.

        The match is case-insensitive substring against the keyword list
        defined in :data:`POPULATION_GROUPS`. Custom group strings that
        are not predefined codes fall through to a direct substring
        filter against ``row_population``.

        Parameters
        ----------
        row_population : str
            Population label from the AFND TSV row.
        target_group : str
            Either a :data:`POPULATION_GROUPS` key (``EUR``, ``AFR``,
            ...), ``"global"``, or a custom substring.

        Returns
        -------
        bool
        """
        if target_group == "global":
            return True
        target_group = _canonical_group(target_group)
        if target_group in POPULATION_GROUPS:
            # Resolve the row to its single canonical group by precedence
            # (AFR before AMR), then compare for equality. A naive substring
            # match would mis-assign "African American" to AMR via the
            # "american" substring, double-counting it into AMR (and AFR) and
            # skewing AMR frequencies.
            return _classify_population_group(row_population) == target_group
        # Custom (non-predefined) string — direct substring filter.
        return target_group.lower() in row_population.lower()

    # ---- Internal helpers ----

    def _locate_tsv(self) -> Optional[Path]:
        candidates = [
            self.local_dir / DEFAULT_TSV_FILENAME,
            self.local_dir / "afnd.tsv",
            self.local_dir / "allele_frequencies.tsv",
        ]
        for p in candidates:
            if p.is_file():
                return p
        for match in self.local_dir.rglob("*.tsv"):
            return match
        return None

    def _parse_tsv(self, path: Path) -> List[Dict[str, str]]:
        rows: List[Dict[str, str]] = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                raise AFNDDatabaseError(f"AFND TSV has no header: {path}")
            cols = {(n or "").strip().lower(): n for n in reader.fieldnames}

            def pick(*candidates: str) -> Optional[str]:
                for c in candidates:
                    real = cols.get(c.lower())
                    if real:
                        return real
                return None

            f_allele = pick("allele", "alleles")
            f_pop = pick("population", "pop")
            f_group = pick("population group", "population_group", "group", "region")
            f_freq = pick("frequency", "allele frequency", "%", "freq")
            f_size = pick(
                "sample size",
                "sample_size",
                "n",
                "individuals",
                "alleles_no",
            )
            if f_allele is None or f_freq is None:
                raise AFNDDatabaseError(
                    f"AFND TSV is missing expected columns: {reader.fieldnames}"
                )
            for row in reader:
                try:
                    freq = float(row.get(f_freq, "") or "nan")
                except ValueError:
                    continue
                if not (0.0 <= freq <= 1.0):
                    # Values may be expressed as a percentage (0–100); normalize.
                    if 0.0 <= freq <= 100.0:
                        freq = freq / 100.0
                    else:
                        continue
                try:
                    size = int(float(row.get(f_size, 0) or 0)) if f_size else 0
                except ValueError:
                    size = 0
                rows.append(
                    {
                        "allele": _strip_hla((row.get(f_allele) or "").strip()),
                        "population": (row.get(f_pop) or "").strip() if f_pop else "",
                        "group": (row.get(f_group) or "").strip() if f_group else "",
                        "frequency": str(freq),
                        "sample_size": str(size),
                    }
                )
        return rows

    def _reindex(self) -> None:
        """
        Rebuild the ``allele -> rows`` index used by :meth:`_match_rows`.

        Must be called whenever ``self._rows`` is (re)assigned — including by
        subclasses such as :class:`~hlante.db.nmdp.NMDPClient` that override
        :meth:`load` — so per-allele lookups stay O(1).
        """
        index: Dict[str, List[Dict[str, str]]] = {}
        for row in self._rows:
            index.setdefault(row["allele"], []).append(row)
        self._by_allele = index

    def _match_rows(self, allele: str, target_group: str) -> List[Dict[str, str]]:
        """
        Apply the population filter and the min_sample_size threshold.
        """
        matches: List[Dict[str, str]] = []
        # O(1) lookup via the allele index built at load time, instead of an
        # O(N) scan of every row on each query.
        for row in self._by_allele.get(allele, []):
            try:
                size = int(row["sample_size"])
            except ValueError:
                size = 0
            if size < self.min_sample_size:
                continue

            # Check both the row's Population Group column and its
            # Population label. Either must match.
            row_text = f"{row['group']} {row['population']}"
            if self._matches_group(row_text, target_group):
                matches.append(row)
        return matches

    def _aggregate(
        self,
        rows: List[Dict[str, str]],
        *,
        allele: str,
        source_resolution: int,
        population_group_label: str,
        is_estimated: bool,
    ) -> AllelFrequency:
        """
        Combine multiple rows using a sample-size-weighted mean.
        """
        total_size = 0
        weighted_sum = 0.0
        for r in rows:
            try:
                f = float(r["frequency"])
                n = int(r["sample_size"])
            except ValueError:
                continue
            if n <= 0:
                continue
            weighted_sum += f * n
            total_size += n
        if total_size <= 0:
            # Fall back to a plain arithmetic mean when sizes are missing.
            total_size = len(rows)
            weighted_sum = sum(float(r["frequency"]) for r in rows)
        aggregated_freq = weighted_sum / max(total_size, 1)
        return AllelFrequency(
            allele=allele,
            frequency=round(aggregated_freq, 6),
            population_group=population_group_label,
            sample_size=total_size,
            source_resolution=source_resolution,
            is_estimated=is_estimated,
            populations_aggregated=len(rows),
            source=self._source_label,
        )


def _strip_hla(allele: str) -> str:
    if allele.upper().startswith("HLA-"):
        return allele[4:]
    return allele


__all__ = [
    "AFNDClient",
    "AllelFrequency",
    "AFNDDatabaseError",
    "POPULATION_GROUPS",
    "SLOWIKOWSKI_URL",
    "DEFAULT_LOCAL_DIR",
    "DEFAULT_TSV_FILENAME",
    "DEFAULT_POPULATION_GROUP",
    "DEFAULT_MIN_SAMPLE_SIZE",
]
