"""File to setup the package."""

from __future__ import annotations

from setuptools import find_packages, setup

setup(
    name="phokimo",
    packages=find_packages(),
    entry_points={
        "console_scripts": ["phokimo = phokimo.main:main"],
    },
)
