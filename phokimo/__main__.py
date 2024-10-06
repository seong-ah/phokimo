"""The main entry point for the application."""

from __future__ import annotations

import os
import sys
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from phokimo.src.additional_functions import expfitting, product_ratio, dict_generator, toml_generator
from phokimo.src.ode_builder import construct_ode
from phokimo.src.rate_constants import RateCalculator
from phokimo.src.terachem_values import Reactions, State_Values
from phokimo.src.toml_reader import TomlReader
from tcgm_lib.convert.converter import energy_unit

def main() -> None:
    """Run the application."""
    """ Read data from toml file. """
    
    if len(sys.argv) < 2:
        print("Usage: python3 -m phokimo filename.toml")
        sys.exit(1)

    toml_file = sys.argv[1]
    current_dir = os.getcwd()
    toml_file_path = os.path.join(current_dir, toml_file)
    toml_data = TomlReader(toml_file_path)

    num_states = toml_data.num_states()
    visualize_state_list_name = toml_data.visualize_state_list_name()
    state_dict = toml_data.name_to_num
    start_conc = toml_data.initial_conc()

    state_data = State_Values(toml_data)
    state_list_energy = state_data.state_relative_list_energy()

    rate_formula = RateCalculator()
    reactions = Reactions(toml_data, rate_formula, state_list_energy)

    spin_list_dict = toml_data.spin_list_dict()

    product_list_num = toml_data.product_nums
    product_list_name_vis = toml_data.product_names_vis

    """ Plot Energies(eV) """
    
    relative_energy_ev = energy_unit(state_list_energy, "j/mol", "ev")
    visualize_state_list_ev = [np.round(x, 2) for x in relative_energy_ev]

    plt.figure(figsize=(12, 8))
    plt.scatter(visualize_state_list_name, visualize_state_list_ev, s=900, marker="_", linewidth=2, zorder=3)
    [plt.text(x, y, str(y), ha="center", va="bottom") for x, y in zip(range(len(visualize_state_list_name)), visualize_state_list_ev)]
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("rel.energy [eV]")
    plt.tight_layout()
    plt.savefig('phokimo_state_energy.png')
    plt.show()
    plt.clf()

    """ Solving ode """
    dEs = reactions.dEs(state_list_energy)
    rates = reactions.rates(state_list_energy)

    table = toml_data.reactions_num
    table_name = toml_data.reactions_vis

    duration = toml_data.duration()
    spacing = 100000
    time = np.linspace(0, duration, spacing)
    conc = np.zeros((spacing, num_states))
    func = partial(construct_ode, table = table, rates=rates)
    conc = odeint(func, start_conc, time)
    x_axis = [i * 10 ** (15) for i in time] #fs

    plt.plot(x_axis, conc)
    plt.legend(visualize_state_list_name)
    plt.xlabel("time [fs]")
    plt.ylabel("Concentration")
    plt.tight_layout()
    plt.savefig('phokimo_kinetics_state.png')
    plt.show()
    plt.clf()

    """ Generate output TOML file """

    value_dict = dict_generator(visualize_state_list_name, state_list_energy, table_name, dEs, rates)
    exp_dict = expfitting(time, x_axis, state_dict, spin_list_dict, conc)
    frac_dict = product_ratio(spacing, conc, product_list_name_vis, product_list_num)
    toml_generator(value_dict, exp_dict, frac_dict)

if __name__ == "__main__":
    main()
