"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import graphviz as gp
from graphviz import Source
from scipy.optimize import curve_fit


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

def product_fraction(spacing: int, time: np.ndarray, x_axis: np.ndarray, conc: np.ndarray, product_list_name: list[str], product_list_num: list[int]) -> None:
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

    plt.plot(x_axis, fractions)
    plt.legend(product_list_name)
    plt.xlabel("time (s)")
    plt.ylabel("Concentration fraction")
    for i in range(len(product_list_num)):
        print(str(product_list_name[i]), "{:.0f}".format(fractions[spacing-1][i] * 100)+"%")

def expfitting(time: np.ndarray, x_axis: np.bdarray, spin_list_dict: dict, conc: np.ndarray):
    def exponential_func(x, a, k, c):
        return a * np.exp(k * x) + c
    # initial_guess = [-0.3, 10 ** 13, 0]
    x = time
    for i in range(len(spin_list_dict)):
        y = np.zeros(len(time))
        list = spin_list_dict[i]
        for j in range(len(time)):
            for k in list:
                y[j] += conc[j][k]
        label = "S" + str(i)
        plt.scatter(x_axis, y, label = label)
        try:
            popt, pcov = curve_fit(exponential_func, x, y)
            # popt, pcov = curve_fit(exponential_func, x, y, p0 = initial_guess, bounds = (0, np.inf))
            a, k, c = popt
            print(a, k, c)
            print("Time constant ", label, " : ", np.round(1 / abs(k) * 10 ** (12), 2), "ps")
            y_fit = exponential_func(x, a, k, c)
            plt.plot(x_axis, y_fit)
        except RuntimeError as e:
            print("Optimal parameters not found:", e)
        
        print(label, ": ", np.round(y[len(time)-1], 2) * 100, "%")
        y = np.zeros(len(time))

    plt.legend()
    plt.xlabel("time (fs)")
    plt.ylabel("Concentration fraction")
    plt.show()