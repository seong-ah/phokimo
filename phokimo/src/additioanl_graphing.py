"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def graph_builder(state_list_name: list, state_list_num: list, graph_table_name: list, graph_table_num: list) -> list[tuple, ]:
    """A tool to build graph represents the reaction relations.

    Args:
        state_list_name (list): list of name of each state
        state_list_num (list): list of numbering of each state
        graph_table_name (list): list of graph edge tuples that represent the reaction with state names
        graph_table_num (list): list of graph edge tuples that represent the reaction with state numbering

    Returns:
        list: list of tuples that represent graph edges (initial/transition state, transition/final state)
    """
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

    return table


def fraction(spacing: int, number_of_products: int, time: np.ndarray, conc: np.ndarray) -> None:
    """A tool to plot fraction of products.

    Args:
        spacing (int): spacing of time
        number_of_products (int): total number of products
        time (np.ndarray): time
        conc (np.ndarray): concentration of corresponding time
    """
    
    dim = (spacing, number_of_products)  # [spacing, number of products]

    fractions = np.zeros(dim)  # [TAB, CAB]
    for i in range(1, spacing):
        denominator = conc[i][1] + conc[i][2]
        tab = conc[i][1]
        cab = conc[i][2]
        fractions[i][0] = tab / denominator
        fractions[i][1] = cab / denominator

    labels = ["TAB", "CAB"]

    plt.plot(time, fractions)
    plt.legend(labels)
    plt.show()
    print(fractions[spacing-1][0], fractions[spacing-1][1])
