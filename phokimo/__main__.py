"""The main entry point for the application."""

from __future__ import annotations

import os
from functools import partial

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.integrate import odeint

from phokimo.src.ode_builder import construct_ode
from phokimo.src.rate_constants import RateCalculator
from phokimo.src.terachem_values import Reactions, State_Values
from phokimo.src.toml_reader import TomlReader


def main() -> None:
    """ Run the application. """
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
    state_list_energy = state_data.state_list_energy(calculation_path)

    graph_table_name = reactions.graph_table_name()
    graph_table_num = reactions.graph_table_num()

    dEs = reactions.dEs(state_list_energy)
    dEs_ev = [x * 0.0000103643 for x in dEs]
    rates = reactions.rates(state_list_energy)

    " Simple print setting for debugging "
    def custom_formatter(x):
        if x == 0:
            return '0'
        else:
            return f'{x:.4f}'
    
    np.set_printoptions(formatter={'float_kind': custom_formatter})
    print(dEs)
    print(rates)

    """ Plot Energies(eV) """

    state_list_ev = [(x - state_list_hartree[1]) * 27.2114 for x in state_list_hartree] # Relative energy from TAB in eV
    visualize_state_list_ev = [np.round(x, 3) for x in state_list_ev]

    plt.scatter(state_list_name, visualize_state_list_ev, marker="o")

    for i, y in enumerate(visualize_state_list_ev):
        plt.text(state_list_name[i], y, str(y), ha="left", va="bottom", fontsize=10)

    plt.show()

    """ Graph creation """

    # Create a graph
    name_graph = nx.Graph()
    num_graph = nx.Graph()

    # Add nodes
    name_graph.add_nodes_from(state_list_name)
    num_graph.add_nodes_from(state_list_num)

    # Add edges
    name_graph.add_edges_from(graph_table_name)
    num_graph.add_edges_from(graph_table_num)

    # Convert the graph to a dictionary of lists
    table = nx.to_dict_of_lists(num_graph)

    # Visualize the graph
    nx.draw(name_graph, with_labels=True, font_weight="bold")
    plt.show()

    nx.draw(num_graph, with_labels=True, font_weight="bold")
    plt.show()

    """ Solving ode """

    time = np.linspace(0, 10**(-12), 1000)

    func = partial(construct_ode, table=table, rates=rates)
    conc = odeint(func, start_conc, time)

    labels = []
    for i in range(toml_data.num_states()):
        labels.append(toml_data.state_name(i))

    plt.plot(time, conc)
    plt.legend(labels)
    plt.show()

    """ Graphing fractions """

    dim = (1000, 2)

    fractions = np.zeros(dim)  # [TAB, CAB]
    for i in range(1, 1000):
        denominator = conc[i][1] + conc[i][2]
        tab = conc[i][1]
        cab = conc[i][2]
        fractions[i][0] = tab / denominator
        fractions[i][1] = cab / denominator

    plt.plot(time, fractions)
    plt.show()

if __name__ == "__main__":
    main()
