"""The main entry point for the application."""

from __future__ import annotations

import os
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from phokimo.src.additioanl_graphing import fraction, graph_builder
from phokimo.src.ode_builder import construct_ode
from phokimo.src.rate_constants import RateCalculator
from phokimo.src.terachem_values import Reactions, State_Values
from phokimo.src.toml_reader import TomlReader
from tcgm_lib.convert.converter import energy_unit


def main() -> None:
    """Run the application."""
    """ Read data from toml file. """

    # Get the directory of the currently executing Python script
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Assume the TOML file is in the same directory as the Python script
    # Modify "s1_dynamics.toml" for suitable set-up toml file
    toml_file_path = os.path.join(current_dir, "s1_dynamics.toml")

    # Absolute path of calculation folders: assume that all folders have same structure in parallel (calculation_path/sp/tc.out)
    calculation_path = "/home/guests/schoi/kinetic/azobenzene/"

    toml_data = TomlReader(toml_file_path)
    num_states = toml_data.num_states()

    state_list_name = toml_data.state_list_name()
    state_list_num = toml_data.state_list_num()
    start_conc = toml_data.start_conc()

    rate_formula = RateCalculator()
    reactions = Reactions(toml_data, rate_formula)

    state_data = State_Values(toml_data)
    state_list_hartree = state_data.state_list_hartree(calculation_path)
    print(state_list_hartree)
    state_list_energy = state_data.state_list_energy(calculation_path)
    print(state_list_energy)

    graph_table_name = reactions.graph_table_name()
    graph_table_num = reactions.graph_table_num()

    rates = reactions.rates(state_list_energy)

    " Simple print setting for debugging "

    def custom_formatter(x):
        if x == 0:
            return "0"
        else:
            return f"{x:.4f}"

    np.set_printoptions(formatter={"float_kind": custom_formatter})
    np.set_printoptions(suppress=False, precision=2)

    print(rates)

    """ Plot Energies(eV) """

    relative_energy = [(x - state_list_hartree[1]) for x in state_list_hartree] # Eh
    relative_energy_numpy = np.asarray(relative_energy)
    relative_energy_ev = energy_unit(relative_energy_numpy, "eh", "ev") # Relative energy from TAB in eV
    visualize_state_list_ev = [np.round(x, 3) for x in relative_energy_ev]
    print(relative_energy_ev)

    plt.scatter(state_list_name, visualize_state_list_ev, marker="o")

    for i, y in enumerate(visualize_state_list_ev):
        plt.text(state_list_name[i], y, str(y), ha="left", va="bottom", fontsize=10)

    plt.show()

    """ Solving ode """

    #table = graph_builder(state_list_name, state_list_num, graph_table_name, graph_table_num)
    table = toml_data.reaction_list()
    print(table)

    spacing = 1000
    time = np.linspace(0, 10 ** (-12), spacing)

    func = partial(construct_ode, table=table, rates=rates)
    conc = odeint(func, start_conc, time)

    labels = []
    for i in range(toml_data.num_states()):
        labels.append(toml_data.state_name(i))

    plt.plot(time, conc)
    plt.legend(labels)
    plt.show()

    """ Plotting fractions """
    number_of_products = 2
    fraction(spacing, number_of_products, time, conc)

if __name__ == "__main__":
    main()
