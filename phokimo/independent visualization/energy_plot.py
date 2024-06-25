"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import os
import numpy as np
from phokimo.src.toml_reader import TomlReader

def energyplot():
    # Get the directory of the currently executing Python script
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Assume the TOML file is in the same directory as the Python script
    # Modify "s1_dynamics.toml" for suitable set-up toml file
    toml_file_path = os.path.join(current_dir, "s1_dynamics.toml")

    toml_data = TomlReader(toml_file_path)