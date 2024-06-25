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
    visualize_state_list_name = toml_data.visualize_state_list_name()
    start_conc = toml_data.start_conc()

    rate_formula = RateCalculator()
    reactions = Reactions(toml_data, rate_formula)

    state_data = State_Values(toml_data)
    state_list_hartree = state_data.state_list_hartree(calculation_path)
    state_list_energy = state_data.state_list_energy(calculation_path)

    rates = reactions.rates(state_list_energy)

    reactant_name = toml_data.reactant_name()
    reactant_num = toml_data.reactant_num()
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

    print(rates)

    """ Plot Energies(eV) """

    relative_energy = [(x - state_list_hartree[1]) for x in state_list_hartree] # Eh
    relative_energy_numpy = np.asarray(relative_energy)
    relative_energy_ev = energy_unit(relative_energy_numpy, "eh", "ev") # Relative energy from TAB in eV
    visualize_state_list_ev = [np.round(x, 2) for x in relative_energy_ev]
    print(visualize_state_list_name)
    print(visualize_state_list_ev)

    plt.scatter(visualize_state_list_name, visualize_state_list_ev, s=900, marker="_", linewidth=2, zorder=3)
    [plt.text(x, y, str(y), ha="left", va="bottom", fontsize=10) for x, y in zip(range(len(visualize_state_list_name)), visualize_state_list_ev)]

    plt.show()

    """ Solving ode """
    table = toml_data.reaction_list()

    spacing = 10000
    time = np.linspace(0, 10 ** (-11), spacing)

    func = partial(construct_ode, table=table, rates=rates)
    conc = odeint(func, start_conc, time)

    plt.plot(time, conc)
    plt.legend(visualize_state_list_name)
    plt.show()

    """ Plotting fractions """
    fraction(spacing, time, conc, product_list_name, product_list_num)

if __name__ == "__main__":
    main()
