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
    toml_relative_file_path = os.path.join(current_dir, "..", "azobenzene_s1_dynamics.toml")
    toml_absolute_file_path = os.path.abspath(toml_relative_file_path)

    toml_data = TomlReader(toml_absolute_file_path)
    list = toml_data.visualize_state_list_name()
    print(list)
    name_list = ['TAB*', 'TAB', 'CAB', 'S1 min', 'S1/S0 reactive', 'S1/S0 unreactive', 'S1 bar']
    energy_list = [2.95, 0.0, 0.59, 2.33, 2.27, 2.79, 2.47]

    x = ["S1/S0 unreactive", "TAB", "TAB", "S1 min", "S1 bar", "S1/S0 reactive", "CAB"]
    y = [2.79, 0, 2.95, 2.33, 2.47, 2.27, 0.59]

    pairs = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 5), (5, 1), (5, 6)]

    unique_x = []
    for item in x:
        if item not in unique_x:
            unique_x.append(item)
    print(unique_x)
    x_numeric_map = {label: i for i, label in enumerate(unique_x)}
    print(x_numeric_map)
    x_numeric = np.array([x_numeric_map[label] for label in x])
    print(x_numeric)

    plt.scatter(x_numeric, y, s=900, marker="_", linewidth=2, zorder=3, color = 'maroon')

    for i, j in pairs:
        offset = 0.1
        if x[i] == x[j]:
            start_x = x_numeric[j]
            start_y = y[j]
            end_x = x_numeric[i]
            end_y = y[i]
        else: 
            if i > j:
                start_x = x_numeric[j] + offset
                start_y = y[j]
                end_x = x_numeric[i] - offset
                end_y = y[i]
            else:
                start_x = x_numeric[i] + offset
                start_y = y[i]
                end_x = x_numeric[j] - offset
                end_y = y[j]
        plt.plot([start_x, end_x], [start_y, end_y], linestyle='--', color='rosybrown')

        # Calculate and annotate the difference
        difference = np.round(y[j] - y[i], 2)
        if (i-j) * difference < 0:
            plt.text((start_x + end_x) / 2 - 0.05, (y[i] + y[j]) / 2 + 0.05, f'{difference}', fontsize=10, ha='right', va='bottom', font = 'Arial', weight = 'bold')
        else:
            plt.text((start_x + end_x) / 2 + 0.05, (y[i] + y[j]) / 2 + 0.05, f'{difference}', fontsize=10, ha='left', va='bottom', font = 'Arial', weight = 'bold')
        # plt.plot([x[i], x[j]], [y[i], y[j]], linestyle='--', color='gray')

    plt.xticks(range(len(unique_x)), unique_x, fontweight = 'bold')

    for i, label in enumerate(y):
        print(i, toml_data.reference_state())
        if toml_data.reference_state() == i:
            plt.text(x_numeric[i], y[i] - 0.02, label, ha='center', va='top', fontsize=10, font = 'Arial', weight = 'bold')
        else:
            plt.text(x_numeric[i], y[i] + 0.01, label, ha='center', va='bottom', fontsize=10, font = 'Arial', weight = 'bold')

    plt.ylabel('Relative energy (eV)', font = 'Arial', fontweight = 'bold', fontsize = 12)
    plt.show()

energyplot()