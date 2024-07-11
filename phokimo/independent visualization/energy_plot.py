"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import os
import numpy as np
from matplotlib import pyplot as plt
from phokimo.src.toml_reader import TomlReader

def energyplot():
    # Get the directory of the currently executing Python script
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Assume the TOML file is in the same directory as the Python script
    # Modify "s1_dynamics.toml" for suitable set-up toml file
    toml_relative_file_path = os.path.join(current_dir, "..", "s1_dynamics.toml")
    toml_absolute_file_path = os.path.abspath(toml_relative_file_path)

    toml_data = TomlReader(toml_absolute_file_path)
    list = toml_data.visualize_state_list_name()
    name_list = ['TAB*', 'TAB', 'CAB', 'S1 min', 'S1/S0 reactive', 'S1/S0 unreactive', 'S1 bar']
    energy_list = [2.95, 0.0, 0.59, 2.33, 2.27, 2.79, 2.47]

    x = ["S1/S0 unreactive", "TAB", "TAB", "S1 min", "S1 bar", "S1/S0 reactive", "CAB"]
    y = [2.79, 0.0, 2.95, 2.33, 2.47, 2.27, 0.59]

    plt.scatter(x, y, s=900, marker="_", linewidth=2, zorder=3)

    pairs = [(0, 1), (0, 2), (2, 3), (3, 4), (4, 5), (5, 1), (5, 6)]

    for i, j in pairs:
        plt.plot([x[i], x[j]], [y[i], y[j]], linestyle='--', color='blue')

    plt.show()

energyplot()