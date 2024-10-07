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

    ['TAB*', 'TAB', 'CAB', 'S1 min', 'S1/S0 reactive', 'S1/S0 unreactive', 'S1 bar']
    [2.95, 0.0, 0.59, 2.33, 2.27, 2.79, 2.47]

    ['TAB*', 'TAB', 'CAB', 'S1 min', 'S1/S0 reactive', 'S1/S0 unreactive', 'S1 bar', 'S2 min', 'S2/S1 planar', 'S2/S1 twisted', 'S2 bar']
    [4.18, 0.0, 0.59, 2.33, 2.27, 2.78, 2.47, 3.68, 3.79, 2.64, 3.83]
    [4.18, 0.0, 0.59, 2.33, 2.27, 2.78, 2.47, 3.68, 3.79, 2.64, 4.19]

    ['Ethylene*', 'S1/S0 tw-py', 'S1/S0 ethy', 'C2H3 + H', 'C2H2 + H2', 'C2H2 + H + H', 'S1 min', 'Ethylene']
    [7.56, 5.0, 5.15, 4.94, -6.49, 0.99, 7.33, 0.0]

    #x = ["S1/S0 unreactive", "TAB", "TAB", "S1 min", "S1 bar", "S1/S0 reactive", "CAB"]
    #y = [2.79, 0, 2.95, 2.33, 2.47, 2.27, 0.59]
    #x = [" CAB ", " S1/S0 reactive ", "S2/S1 twisted", "S2 barrier", "S2 min", "TAB", "TAB", "S2/S1 planar", "S1/S0 unreactive", "S1 min", "S1 bar", "S1/S0 reactive", "CAB"]
    #y = [0.59, 2.27, 2.64, 3.83, 3.68, 0.0, 4.18, 3.79, 2.78, 2.33, 2.47, 2.27, 0.59]
    #y = [0.59, 2.27, 2.64, 4.19, 3.68, 0.0, 4.18, 3.79, 2.78, 2.33, 2.47, 2.27, 0.59]
    #x = ["S1/S0 tw-py", "C2H3 + H",  "Ethylene", "Ethylene", "C2H2 + H + H", "S1/S0 ethy", "C2H2 + H2"]
    #x = ["Ethylene", "Ethylene", "S1/S0 tw-py", "C2H3 + H", "C2H2 + H + H", "C2H2 + H2", "S1/S0 ethy"]
    #y = [0.0, 7.56, 5.0, 4.94, 0.99, -6.49, 5.15]
    x = ["C2H2 + H2", " C2H2 + H + H ", " C2H3 + H ", "S1/S0 ethy", "Ethylene", "Ethylene", "S1/S0 tw-py", "C2H3 + H", "C2H2 + H + H"]
    y = [-6.49, 0.99, 4.94, 5.15, 0.0, 7.56, 5.0, 4.94, 0.99]

    #pairs = [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 5), (5, 1), (5, 6)]
    #pairs = [(1,0), (2,1), (3,2), (4,3), (5,6), (6,4), (6,7), (7,8), (8,5), (7,9), (9,10), (10,11), (11,12)]
    #pairs = [(0,1), (1,2), (1,6), (2,3), (2,4), (6,3), (6,4), (6,5)]
    pairs = [(4,5), (5,3), (3,2), (3,1), (3,0), (5,6), (6,7), (6,8)]

    unique_x = []
    for item in x:
        if item not in unique_x or item == "CAB" or item == "S1/S0 reactive" or item == "C2H2 + H + H" or item == "C2H3 + H":
            unique_x.append(item)
    print(unique_x)

    x_numeric_map = {label: i for i, label in enumerate(unique_x)}
    print(x_numeric_map)
    x_numeric = np.array([x_numeric_map[label] for label in x])
    print(x_numeric)

    plt.scatter(x_numeric, y, s=900, marker="_", linewidth=2, zorder=3, color = 'maroon')

    for i, j in pairs:
        offset = 0.123
        if x[i] == x[j]:
            start_x = x_numeric[j]
            start_y = y[j]
            end_x = x_numeric[i]
            end_y = y[i]
            plt.annotate('', xy=(start_x, start_y), xytext=(end_x, end_y), arrowprops=dict(arrowstyle='->', linestyle='--', color='rosybrown', lw=1.5))

        else: 
            if i > j:
                start_x = x_numeric[j] + offset
                start_y = y[j]
                end_x = x_numeric[i] - offset
                end_y = y[i]
                plt.annotate('', xy=(start_x, start_y), xytext=(end_x, end_y), arrowprops=dict(arrowstyle='->', linestyle='--', color='rosybrown', lw=1.5))
            else:
                start_x = x_numeric[i] + offset
                start_y = y[i]
                end_x = x_numeric[j] - offset
                end_y = y[j]
                plt.annotate('', xy=(end_x, end_y), xytext=(start_x, start_y), arrowprops=dict(arrowstyle='->', linestyle='--', color='rosybrown', lw=1.5))
        #plt.plot([start_x, end_x], [start_y, end_y], linestyle='--', color='rosybrown')
        #plt.arrow(start_x, start_y, end_x - start_x, end_y - start_y, linestyle='--', head_width=0.05, head_length=0.1, length_includes_head=True, color = 'rosybrown')



        # Calculate and annotate the difference
        difference = np.round(y[j] - y[i], 2)
        if x_numeric[i] == x_numeric[j]:
            plt.text((start_x + end_x) / 2 + 0.05, (y[i] + y[j]) / 2, f'{difference}', fontsize=18, ha='left', va='bottom', font = 'Arial')
        else:
            """
            if (i, j) == (7,8):
                plt.text((start_x + end_x) / 2 - 0.08, (y[i] + y[j]) / 2 - 0.1, f'{difference}', fontsize=18, ha='right', va='bottom', font = 'Arial')
                
                elif (i, j) == (4,3):
                plt.text((start_x + end_x) / 2 - 0.05, (y[i] + y[j]) / 2 - 0.2, f'{difference}', fontsize=18, ha='right', va='bottom', font = 'Arial')
                
            else:
            """
            if (i-j) * difference < 0:
                plt.text((start_x + end_x) / 2 - 0.05, (y[i] + y[j]) / 2 + 0.05, f'{difference}', fontsize=18, ha='right', va='bottom', font = 'Arial')
            else:
                plt.text((start_x + end_x) / 2 + 0.05, (y[i] + y[j]) / 2 + 0.05, f'{difference}', fontsize=18, ha='left', va='bottom', font = 'Arial')
        # plt.plot([x[i], x[j]], [y[i], y[j]], linestyle='--', color='gray')

    plt.xticks(range(len(unique_x)), unique_x, rotation = 30, fontsize = 15)

    """
    for i, label in enumerate(y):
        print(i, toml_data.reference_state())
        if toml_data.reference_state() == i:
            plt.text(x_numeric[i], y[i] - 0.02, label, ha='center', va='top', fontsize=10, font = 'Arial')
        else:
            plt.text(x_numeric[i], y[i] + 0.01, label, ha='center', va='bottom', fontsize=10, font = 'Arial')

    """
    """
    ax = plt.gca()
    for i, label in enumerate(ax.get_xticklabels()):
        if i % 2 == 0:
            label.set_y(0)  # Lower position for even-indexed labels
        else:
            label.set_y(-0.01)  # Higher position for odd-indexed labels
    """

    plt.tight_layout()  # Adjust layout to make room for modified labels   

    plt.ylabel('rel. energy (eV)', font = 'Arial', fontsize = 25)
    plt.show()

energyplot()