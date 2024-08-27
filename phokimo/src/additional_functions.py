"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import graphviz as gp
from graphviz import Source
from scipy.optimize import curve_fit
import toml
import os

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

def rate_dict(rates: np.ndarray, visualize_state_list_name: list):
    rate_dict = {}
    for i in range(len(visualize_state_list_name)):
        for j in range(len(visualize_state_list_name)):
            if rates[i][j] != 0:
                temporary_dict = {}
                temporary_dict[visualize_state_list_name[j]] = rates[i][j]
                rate_dict[visualize_state_list_name[i]] = {}
                rate_dict[visualize_state_list_name[i]] = temporary_dict
                print(visualize_state_list_name[i], visualize_state_list_name[j], f"{rates[i][j]:.2e}")
    return rate_dict


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
    plt.legend(product_list_name, fontsize = 18)
    plt.xlabel("time (s)", fontsize = 25)
    plt.ylabel("Concentration fraction", fontsize = 25)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    for i in range(len(product_list_num)):
        print(str(product_list_name[i]), "{:.0f}".format(fractions[spacing-1][i] * 100)+"%")

def product_ratio(spacing: int, conc: np.ndarray, product_list_name: list[str], product_list_num: list[int]) -> None:
    frac_dict = {}
    
    for i in range(len(product_list_name)):
        key = product_list_name[i]
        if key in frac_dict.keys():
            frac_dict[key] += conc[spacing-1][product_list_num[i]]
        else:
            frac_dict[key] = conc[spacing-1][product_list_num[i]]

    sum = 0.0
    for item in frac_dict.values():
        print(item, type(item))
        sum += item

    for key, item in frac_dict.items():
        frac_dict[key] = item / sum * 100

    print(frac_dict)
    return frac_dict

def expfitting(time: np.ndarray, x_axis: np.ndarray, spin_list_dict: dict, conc: np.ndarray):
    def exponential_func(x, a, k, c):
        return a * np.exp(k * x) + c
    # initial_guess = [-0.3, 10 ** 13, 0]
    x = time
    y_list = []
    y_fit_list = []
    label_list = []
    expfitting_result = {}

    for i in range(len(spin_list_dict)):
        y = np.zeros(len(time))
        list = spin_list_dict[i]
        for j in range(len(time)):
            for k in list:
                y[j] += conc[j][k]
        y_list.append(y)
        label = "S" + str(i)
        expfitting_result[label] = {}
        label_list.append(label)
        plt.scatter(x_axis, y, label = label)
        try:
            popt, pcov = curve_fit(exponential_func, x, y)
            # popt, pcov = curve_fit(exponential_func, x, y, p0 = initial_guess, bounds = (0, np.inf))
            a, k, c = popt
            #print(a, k, c)
            time_constant = 1 / abs(k) * 10 ** (15)
            expfitting_result[label]['Time constant'] = time_constant
            print("Time constant ", label, " : ", np.round(1 / abs(k) * 10 ** (15), 4), "fs")
            y_fit = exponential_func(x, a, k, c)
            y_fit_list.append(y_fit)
        except RuntimeError as e:
            print("Optimal parameters not found:", e)
        expfitting_result[label]['Fraction'] = y[len(time)-1]
        print(label, ": ", np.round(y[len(time)-1], 4) * 100, "%")
        y = np.zeros(len(time))

    plt.legend(fontsize = 18)
    plt.xlabel("time [fs]", fontsize = 25)
    plt.ylabel("Concentration fraction", fontsize = 25)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    plt.show()

    y_fit_transpose = np.transpose(np.array(y_fit_list))
    
    for y in y_list:
        plt.scatter(x_axis, y)
    plt.plot(x_axis, y_fit_transpose)
    plt.legend(label_list, fontsize = 18)
    plt.xlabel("time [fs]", fontsize = 25)
    plt.ylabel("Concentration fraction", fontsize = 25)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.show()    

    print(expfitting_result)

def toml_update(toml_file, frac_list: dict, expfitting_result: dict, dict_rate: dict) -> None:

    current_dir = os.path.dirname(os.path.abspath(__file__))

    toml_file_path = os.path.join(current_dir, '..', toml_file)
    
    with open(toml_file_path, 'r') as f:
        config = toml.load(f)
    
    if 'Rates' not in config:
        config['Rates']= {}

    config['Rates'] = dict_rate

    if 'Exponential fitting' not in config:
        config['Exponential fitting'] = {}
    
    config['Exponential fitting'] = expfitting_result
    
    if 'Product ratio' not in config:
        config['Product ratio'] = {}
    config['Product ratio'] = frac_list

    with open(toml_file_path, 'w') as f:
        toml.dump(config, f)
