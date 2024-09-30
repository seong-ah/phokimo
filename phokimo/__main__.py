"""The main entry point for the application."""

from __future__ import annotations

import os
import sys
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from phokimo.src.additional_functions import product_fraction, graph_builder, expfitting, product_ratio, rate_dict, toml_update
from phokimo.src.ode_builder import construct_ode
from phokimo.src.rate_constants import RateCalculator
from phokimo.src.terachem_values import Reactions, State_Values
from phokimo.src.toml_reader import TomlReader
from tcgm_lib.convert.converter import energy_unit

def main() -> None:
    """Run the application."""
    """ Read data from toml file. """

    # INPUT
    if len(sys.argv) < 2:
        print("Usage: python3 __main__.py <name_of_toml_file>")
        sys.exit(1)

    toml_file = sys.argv[1]  # Get the TOML file from the command line
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Assume the TOML file is in the same directory as the __main__.py
    # Modify toml_file_path for suitable set-up toml file
    toml_file_path = os.path.join(current_dir, toml_file)

    toml_data = TomlReader(toml_file_path)

    # Absolute path of calculation folders: assume that all folders have same structure in parallel (calculation_path/sp/tc.out)
    calculation_path = toml_data.calculation_path()

    num_states = toml_data.num_states()

    visualize_state_list_name = toml_data.visualize_state_list_name()
    start_conc = toml_data.initial_conc()

    state_data = State_Values(toml_data)
    state_list_energy = state_data.state_relative_list_energy()


    rate_formula = RateCalculator()
    reactions = Reactions(toml_data, rate_formula, state_list_energy)

    graph_group = toml_data.graph_teq_group()
    print("graph group", graph_group)

    # state_list_energy[10] = 385600.0
    # state_list_energy[10] = 364000.0

    # rates = reactions.rates(state_list_energy)

    spin_list_dict = toml_data.spin_list_dict()

    product_list_name = toml_data.product_list_name()
    product_list_num = toml_data.product_list_num()

    " Simple print setting for debugging "

    def custom_formatter(x):
        if x == 0:
            return "0"
        else:
            return f"{x:.4f}"

    np.set_printoptions(formatter={"float_kind": custom_formatter})
    np.set_printoptions(suppress=False, precision=2)

    """ Plot Energies(eV) """
    
    relative_energy_ev = energy_unit(state_list_energy, "j/mol", "ev")
    visualize_state_list_ev = [np.round(x, 2) for x in relative_energy_ev]
    print(visualize_state_list_name)
    print(visualize_state_list_ev)

    plt.scatter(visualize_state_list_name, visualize_state_list_ev, s=900, marker="_", linewidth=2, zorder=3)
    [plt.text(x, y, str(y), ha="center", va="bottom", fontsize=15) for x, y in zip(range(len(visualize_state_list_name)), visualize_state_list_ev)]
    plt.ylabel("rel.energy (eV)", fontsize = 25)
    plt.xticks(fontsize=18, rotation = 45)
    plt.yticks(fontsize=18)
    plt.tight_layout()

    plt.show()

    rates = reactions.rates(state_list_energy)

    dict_rate = rate_dict(rates, visualize_state_list_name)

    """ Solving ode """
    table = toml_data.reactions_num
    print(table)

    duration = toml_data.duration()
    spacing = 100000
    time = np.linspace(0, duration, spacing)

    conc = np.zeros((spacing, num_states))

    func = partial(construct_ode, table = table, rates=rates)
    conc = odeint(func, start_conc, time)
    
    x_axis = [i * 10 ** (15) for i in time] #fs

    plt.plot(x_axis, conc, linewidth = 5)
    plt.legend(visualize_state_list_name, fontsize = 18, ncol = 2)
    plt.xlabel("time [fs]", fontsize = 25)
    plt.ylabel("Concentration", fontsize = 25)
    plt.show()

    """
    visualize_conc_list = []
    visualize_name = []
    
    for i in range(num_states):
        if i == reactant_num or i in product_list_num:
            visualize_conc_list.append(conc_transpose[i])
            visualize_name.append(visualize_state_list_name[i])
    visualize_conc = np.transpose(np.array(visualize_conc_list))

    plt.plot(x_axis, visualize_conc)
    plt.legend(visualize_name)
    plt.xlabel("time [fs]", fontsize = 20)
    plt.ylabel("Concentration", fontsize = 20)
    plt.show()
    c
    """

    expfitting_result = expfitting(time, x_axis, spin_list_dict, conc)
    # timeconstant(product_list_num, time, conc)

    """ Plotting fractions """
    # product_fraction(spacing, time, x_axis, conc, product_list_name, product_list_num)
    frac_list = product_ratio(spacing, conc, product_list_name, product_list_num)

    """ Add results in TOML file """
    # toml_update(toml_file, frac_list, expfitting_result, dict_rate)

if __name__ == "__main__":
    main()
