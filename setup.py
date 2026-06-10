"""
Backward-compatible setup.py for HLAnte.

The primary package configuration lives in ``pyproject.toml``; this
file is kept around only for ``pip install -e .`` and legacy toolchains.
"""

from __future__ import annotations

from setuptools import setup


if __name__ == "__main__":
    setup()
