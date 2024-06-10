"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def graph_builder(state_list_name: list, state_list_num: list, graph_table_name: list, graph_table_num: list, reactant_name: str, reactant_num: int, product_list_name: list[str], product_list_num: list[int]) -> list[tuple]:
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
    name_graph = nx.DiGraph(directed=True)
    num_graph = nx.DiGraph(directed=True)

    # Add nodes
    name_graph.add_nodes_from(state_list_name)
    num_graph.add_nodes_from(state_list_num)

    # Add edges
    name_graph.add_edges_from(graph_table_name)
    num_graph.add_edges_from(graph_table_num)

    pos_name = nx.kamada_kawai_layout(name_graph, scale = 2)
    pos_num = nx.kamada_kawai_layout(num_graph, scale = 5)

    # Define starting & ending nodes
        
    pos_name[reactant_name] = [0.0, 1.0]
    pos_num[reactant_num] = [0.0, 1.0]
    
    for i in range(len(product_list_num)):
        num_states = float(len(product_list_num))
        product_name = product_list_name[i]
        product_num = product_list_num[i]
        spacing = 1.0 / (num_states - 1)
        pos_name[product_name] = [0.0 + (spacing * i), 0.0]
        pos_num[product_num] = [0.0 + (spacing * i), 0.0]
    
    longest_path = []
    for i in range(len(product_list_num)):
        path = nx.dag_longest_path(num_graph, reactant_num, product_list_name[i])
        if len(path) > len(longest_path):
            longest_path[:] = path

    # Visualize the graph
    nx.draw(name_graph, pos_name, with_labels=True)
    plt.show()

    nx.draw(num_graph, pos_num, with_labels=True)
    plt.show()

def fraction(spacing: int, time: np.ndarray, conc: np.ndarray, product_list_name: list[str], product_list_num: list[int]) -> None:
    """A tool to plot fraction of products.

    Args:
        spacing (int): spacing of time
        time (np.ndarray): time
        conc (np.ndarray): concentration of corresponding time
        product_list_name (list[str]): list of the name of the products
        product_list_num (list[int]): list of the numbering of the products
    """
    
    number_of_products = len(product_list_num)
    dim = (spacing, number_of_products)  # [spacing, number of products]

    fractions = np.zeros(dim)  # [TAB, CAB]

    for i in range(1, spacing):
        denominator = 0
        for x in product_list_num:
            denominator += conc[i][x]
        for j in range(number_of_products):
            state = product_list_num[j]
            fractions[i][j] = conc[i][state] / denominator

    plt.plot(time, fractions)
    plt.legend(product_list_name)
    plt.show()
    print("{:.1f}, {:.1f}".format(fractions[spacing-1][0] * 100, fractions[spacing-1][1] * 100))
