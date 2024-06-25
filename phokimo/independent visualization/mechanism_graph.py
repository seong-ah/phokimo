"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import os
import numpy as np
import graphviz as gp
from phokimo.src.toml_reader import TomlReader

def graph_builder():
    """A tool to build graph represents the reaction relations.

    Args:
        state_list_name (list): list of name of each state
        state_list_num (list): list of numbering of each state
        graph_table_num (list): list of graph edge tuples that represent the reaction with state numbering

    Returns:
        list: list of tuples that represent graph edges (initial/transition state, transition/final state)
    """

    # Get the directory of the currently executing Python script
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Assume the TOML file is in the same directory as the Python script
    # Modify "s1_dynamics.toml" for suitable set-up toml file
    toml_file_path = os.path.join(current_dir, "s1_dynamics.toml")

    toml_data = TomlReader(toml_file_path)

    # Create a graph
    name_tree = gp.Digraph()

    name_tree.attr(fontname = 'Arial')

    # Add nodes
    gp_nodes = []
    for i in range(len(state_list_name)):
        list = [str(state_list_num[i]), state_list_name[i]]
        gp_nodes.append(tuple(list))
    for node, label in gp_nodes:
        name_tree.node(node, label, shape = 'plaintext')
    print(gp_nodes)

    # Add edges
    gp_edges = []
    for init, fin in graph_table_num:
        list = [str(init), str(fin)]
        gp_edges.append(tuple(list))
    for edge in gp_edges:
        name_tree.edge(*edge, splines ='line')
    print(gp_edges)

    name_tree.render('example_graph', format='svg', view=True)

state_list_name = toml_data.
state_list_num = result_list = [i for i in range(len(state_list_name))]
graph_table_name = [('TAB*', 's1min'), ('TAB*', 's1s0_unreactive'), ('s1min', 's1s0_unreactive'), ('s1s0_unreactive', 'TAB'), ('s1min', 's1bar'), ('s1bar', 's1s0_reactive'), ('s1s0_reactive', 'TAB'), ('s1s0_reactive', 'CAB'), ('s1s0_unreactive', 'TAB')]
graph_table_num = [(0, 3), (0, 5), (3, 5), (5, 1), (3, 6), (6, 4), (4, 1), (4, 2)]

graph_builder(state_list_name, state_list_num, graph_table_num)