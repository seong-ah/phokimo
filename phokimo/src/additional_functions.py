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

def product_ratio(spacing: int, conc: np.ndarray, product_list_name_vis: list[str], product_list_num: list[int]) -> None:
    """Calculate the product ratio and store in dictionary to generate TOML file.

    Args:
        spacing (int): number of timesteps
        conc (np.ndarray): concentration values as (spacing) x (number of total states) adjacency matrix
        product_list_name_vis (list[str]): list of visualizing name of products
        product_list_num (list[int]): list of product numbering

    Returns:
        dict: dictionary of product ratio values to generate TOML file
    """
    frac_dict = {}
    
    for i in range(len(product_list_name_vis)):
        key = product_list_name_vis[i]
        if key in frac_dict.keys():
            frac_dict[key] += conc[spacing-1][product_list_num[i]]
        else:
            frac_dict[key] = conc[spacing-1][product_list_num[i]]

    sum = 0.0
    for item in frac_dict.values():
        sum += item

    for key, item in frac_dict.items():
        frac_dict[key] = item / sum * 100
        print("Product ratio: ", key, ": ", frac_dict[key], " %")

    frac_dict_return = {}
    frac_dict_return["Product ratios %"] = frac_dict
    return frac_dict_return

def expfitting(time: np.ndarray, x_axis: np.ndarray, state_dict: dict, spin_list_dict: dict, conc: np.ndarray) -> dict:
    """Execute exponential fitting per spin state to calculate time constant.
    Store values as dictionary to generate TOML file.

    Args:
        time (np.ndarray): time
        x_axis (np.ndarray): time in fs format for visualization
        state_dict (dict): mapping between state name and numbering
        spin_list_dict (dict): mapping between spin state and corresponding states
        conc (np.ndarray): starting concentration

    Returns:
        dict: dictionary of values related to the exponential to generate TOML file
    """
    def exponential_func(x, a, k, c):
        return a * np.exp(k * x) + c
    
    x = time
    spacing = len(time)
    y_list =  np.zeros((spacing, len(spin_list_dict))) # time * label
    y_fit_list = []
    label_list = ["name"] * len(spin_list_dict)
    expfitting_result = {}

    counter = 0

    for spin_pair in spin_list_dict:
        if spin_pair[0] == 1:
            label = "S" + str(spin_pair[1])
        elif spin_pair[0] == 2:
            label = "D" + str(spin_pair[1])
        elif spin_pair[0] == 3:
            label = "T" + str(spin_pair[1])
        label_list[counter] = label
        expfitting_result[label] = {}

        for state in spin_list_dict[spin_pair]:
            state_num = state_dict[state]
            for t in range(len(time)):
                y_list[t][counter] += conc[t][state_num]
        counter += 1
    plt.plot(x_axis, y_list)
    plt.legend(label_list)
    plt.xlabel("time [fs]")
    plt.ylabel("Concentration fraction")
    plt.tight_layout()
    plt.savefig('phokimo_kinetics_spin.png')
    plt.show()
    plt.clf()

    y_transpose = np.transpose(y_list) # label * time

    for i in range(len(spin_list_dict)):
        y = y_transpose[i]
        try:
            popt, pcov = curve_fit(exponential_func, x, y)
            a, k, c = popt
            time_constant = 1 / abs(k) * 10 ** (15) #fs
            label = label_list[i]
            expfitting_result[label]['Time constant [fs]'] = time_constant
            print("Time constant ", label, " : ", np.round(1 / abs(k) * 10 ** (15), 4), "fs")
            y_fit = exponential_func(x, a, k, c)
            y_fit_list.append(y_fit)
        except RuntimeError as e:
            print("Optimal parameters not found:", e)
        expfitting_result[label]['Fraction'] = y[len(time)-1]
        print(label, ": ", np.round(y[len(time)-1], 4) * 100, "%")
        
    y_fit_transpose = np.transpose(np.array(y_fit_list))
    fit_label_list = [x + '_fit' for x in label_list]
    
    plt.plot(x_axis, y_list, label = label_list)
    plt.plot(x_axis, y_fit_transpose, label = fit_label_list)
    plt.legend()
    plt.xlabel("time [fs]")
    plt.ylabel("Concentration fraction")
    plt.tight_layout()
    plt.savefig('phokimo_expfitting.png')
    plt.show()
    plt.clf()    

    expfitting_return = {}
    expfitting_return["Exponential fitting values"] = expfitting_result
    return expfitting_return


def dict_generator(visualize_state_list_name: list[str], state_list_energy: np.ndarray, table_name: list[tuple], dEs: np.ndarray, rates: np.ndarray) -> dict:
    """Store numerical values of states and reactions (e.g. energies and rate constants) in dictionary format to generate TOML file.

    Args:
        visualize_state_list_name (list[str]): visualizing name of each state in a list
        state_list_energy (np.ndarray): energy (J/mol) per each state
        table_name (list[tuple]): reaction connection tuples with visualizing name
        dEs (np.ndarray): energy differences of the reaction as N x N matrix (initial -> final, dEs[init][fin] = energy difference)
        rates (np.ndarray): reaction rates as N x N adjacency matrix. N: number of states. (initial -> final, rates[init][fin] = rate constant)
    Returns:
        dict: dictionary of values to generate TOML file
    """
    value_dict = {}
    value_dict["relative energy (J/mol)"] = {}
    value_dict["Reactions"] = {}
    for i in range(len(visualize_state_list_name)):
        state = visualize_state_list_name[i]
        value_dict["relative energy (J/mol)"][state] = state_list_energy[i]

    for reaction in table_name:
        init, fin = reaction
        init_num, fin_num = visualize_state_list_name.index(init), visualize_state_list_name.index(fin)
        if init not in value_dict["Reactions"]:
            value_dict["Reactions"][init] = {}
        if fin not in value_dict["Reactions"][init]:
            value_dict["Reactions"][init][fin] = {}
        value_dict["Reactions"][init][fin]["dE (J/mol)"] = dEs[init_num][fin_num]
        value_dict["Reactions"][init][fin]["rate constant"] = rates[init_num][fin_num]
    return value_dict

def toml_generator(value_dict: dict, exp_dict: dict, frac_dict: dict) -> None:
    """ Generate TOML file by merging the dictionaries.

    Args:
        value_dict (dict): values of states and reactions as dictionary
        exp_dict (dict): values from exponential fitting as dictionary
        frac_dict (dict): product fraction values as dictionary
    """
    data = {**value_dict, **exp_dict, **frac_dict}

    with open('phokimo.toml', 'w') as toml_file:
        toml.dump(data, toml_file)