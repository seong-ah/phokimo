"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import graphviz as gp
from graphviz import Source
import pydot as pd


def graph_builder(state_list_name: list, state_list_num: list, graph_table_name: list, graph_table_num: list) -> list[tuple]:
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
    name_tree = gp.Digraph()


    # Add nodes
    gp_nodes = []
    for i in range(len(state_list_name)):
        list = [str(state_list_num[i]), state_list_name[i]]
        gp_nodes.append(tuple(list))
    for node, label in gp_nodes:
        name_tree.node(node, label)
    print(gp_nodes)

    # Add edges
    gp_edges = []
    for init, fin in graph_table_num:
        list = [str(init), str(fin)]
        gp_edges.append(tuple(list))
    for edge in gp_edges:
        name_tree.edge(*edge)
    print(gp_edges)

    name_tree.render('example_graph', format='svg', view=True)

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
    plt.xlabel("time (s)")
    plt.ylabel("Concentration fraction")
    plt.show()
    print("{:.0f}, {:.0f}".format(fractions[spacing-1][0] * 100, fractions[spacing-1][1] * 100))
