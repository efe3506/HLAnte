"""
HLAnte
=====

HLA genotype annotation toolkit.

HLAnte parses HLA typing tool outputs (ARCAS-HLA, T1K, HLA-HD,
OptiType) and produces a unified clinical annotation report by
querying IPD-IMGT/HLA, GWAS Catalog, PharmGKB, and AFND. It
is distributed as a command-line tool.

Attributes
----------
__version__ : str
    Semantic (PEP 440) package version.
__author__ : str
    Package author / maintainer.
__license__ : str
    License under which the package is released.
"""

from __future__ import annotations

__version__: str = "0.1.0"
__author__: str = "Efe Dallı"
__license__: str = "MIT"

__all__ = [
    "__version__",
    "__author__",
    "__license__",
]
