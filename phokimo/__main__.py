"""The main entry point for the application."""

from __future__ import annotations

import os
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from phokimo.src.additional_graphing import product_fraction, graph_builder, expfitting
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

    # Assume the TOML file is in the same directory as the __main__.py
    # Modify toml_file_path for suitable set-up toml file
    toml_file_path = os.path.join(current_dir, "ethylene_s1_dynamics.toml")

    toml_data = TomlReader(toml_file_path)

    # Absolute path of calculation folders: assume that all folders have same structure in parallel (calculation_path/sp/tc.out)
    calculation_path = toml_data.calculation_path()

    num_states = toml_data.num_states()

    visualize_state_list_name = toml_data.visualize_state_list_name()
    start_conc = toml_data.start_conc()

    rate_formula = RateCalculator()
    reactions = Reactions(toml_data, rate_formula)

    state_data = State_Values(toml_data)
    state_list_energy = state_data.state_list_energy(calculation_path)

    graph_group = toml_data.graph_group()
    print("graph group", graph_group)

    # state_list_energy[10] = 385600.0
    # state_list_energy[10] = 367800.0

    # rates = reactions.rates(state_list_energy)

    spin_list_dict = toml_data.spin_list_dict(1)
    reference_state = toml_data.reference_state()
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

    """ Plot Energies(eV) """
    
    relative_energy_numpy = np.asarray(state_list_energy)
    relative_energy_ev = energy_unit(relative_energy_numpy, "j/mol", "ev") # Relative energy in eV
    # relative_energy_ev_zpe = [a + b for a, b in zip(relative_energy_ev, zpe)]
    # relative_relative_energy_ev_zpe = np.asarray([x - relative_energy_ev_zpe[reference_state] for x in relative_energy_ev_zpe])
    # relative_relative_energy_jmol_zpe = energy_unit(relative_relative_energy_ev_zpe, "ev", "j/mol")
    visualize_state_list_ev = [np.round(x, 2) for x in relative_energy_ev]
    print(visualize_state_list_name)
    print(visualize_state_list_ev)

    plt.scatter(visualize_state_list_name, visualize_state_list_ev, s=900, marker="_", linewidth=2, zorder=3)
    [plt.text(x, y, str(y), ha="center", va="bottom", fontsize=10) for x, y in zip(range(len(visualize_state_list_name)), visualize_state_list_ev)]
    plt.ylabel("rel.energy (eV)")

    plt.show()

    rates = reactions.rates(state_list_energy)

    for i in range(num_states):
        for j in range(num_states):
            if rates[i][j] != 0:
                print(visualize_state_list_name[i], visualize_state_list_name[j], f"{rates[i][j]:.2e}")

    """ Solving ode """
    table = toml_data.reaction_list()
    print(table)

    duration = toml_data.duration()
    spacing = 100000
    time = np.linspace(0, duration, spacing)

    func = partial(construct_ode, table=table, rates=rates)
    conc = odeint(func, start_conc, time)

    x_axis = [i * 10 ** (15) for i in time] #fs

    plt.plot(x_axis, conc)
    plt.legend(visualize_state_list_name)
    plt.xlabel("time (fs)")
    plt.ylabel("Concentration")
    plt.show()

    expfitting(time, x_axis, spin_list_dict, conc)
    # timeconstant(product_list_num, time, conc)

    """ Plotting fractions """
    product_fraction(spacing, time, x_axis, conc, product_list_name, product_list_num)

if __name__ == "__main__":
    main()
