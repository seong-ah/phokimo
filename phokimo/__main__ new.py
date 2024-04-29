"""The main entry point for the application."""

from __future__ import annotations

from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import toml
import os

from scipy.integrate import odeint

from phokimo.src.ode_builder import general_ode, construct_ode

def main() -> None:
    """Run the application."""

    """ Read data from toml file and set-up ode calculation """

    # Get the directory of the currently executing Python script
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Assume the TOML file is in the same directory as the Python script
    toml_file_path = os.path.join(current_dir, "sample.toml")

    # Read the TOML file
    with open(toml_file_path, "r") as f:
        toml_data = toml.load(f)
    print(toml_data)

    # Generate start_conc
    num_states = len(toml_data['states'])
    start_conc = np.zeros(num_states)
    for state, conc in toml_data['concentrations'].items():
            order = int(state[-1])
            start_conc[order] = conc
    
    print(start_conc)

    # Generate state_list from toml
    state_list_name = []
    state_list_num = []

    for state, state_name in toml_data['states'].items():
        order = int(state[-1])
        state_list_name.append(state_name)
        state_list_num.append(order)

    print(state_list_name)
    print(state_list_num)

    time = np.linspace(0, 10, 1000)

    #Generate reaction types table from toml
    graph_table_name = []
    graph_table_num = []

    for reaction, linkage in toml_data['reactions'].items():
        string1 = str(linkage[0])
        print(linkage[0])
        print(string1)
        string2 = str(linkage[1])
        target1 = "state" + string1
        target2 = "state" + string2
        print(target1)
        print(target2)
        species_name1 = toml_data['states'][target1]
        species_name2 = toml_data['states'][target2]
        reaction_species = tuple([species_name1, species_name2])
        graph_table_name.append(reaction_species)
        graph_table_num.append(tuple(linkage))

    print(graph_table_name)
    print(graph_table_num)

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
    print(table)

    # Visualize the graph
    nx.draw(name_graph, with_labels=True, font_weight='bold')
    plt.show()

    nx.draw(num_graph, with_labels=True, font_weight='bold')
    plt.show()

    """ Solving ode """

    #Generate reaction constant matrix from toml
    dim = (num_states, num_states)
    rates = np.zeros(dim)
    
    for init, final in table.items():
        for target_final in final:
            start = int(init)
            end = int(target_final)
            for const, value in toml_data['rate_constants'].items():
                if start == int(const[1]) and end == int(const[2]):
                    rates[start][end] = value
    print(rates)

    func = partial(construct_ode, table=table, rates=rates)
    conc = odeint(func, start_conc, time)

    plt.plot(time, conc)
    plt.show()   

if __name__ == "__main__":
    main()
