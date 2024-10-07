"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import os
import numpy as np
import graphviz as gp
from phokimo.src.toml_reader import TomlReader
from collections import Counter

def remove_duplicates(lst):
    seen = set()
    new_list = []
    for item in lst:
        if item not in seen:
            new_list.append(item)
            seen.add(item)
    return new_list

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
    toml_relative_file_path = os.path.join(current_dir, "..", "azobenzene_s2_dynamics.toml")
    toml_file_path = os.path.abspath(toml_relative_file_path)

    toml_data = TomlReader(toml_file_path)

    state_list_name = toml_data.visualize_state_list_name()
    state_list_num = toml_data.state_list_num()
    graph_table_num = remove_duplicates(toml_data.graph_table_num())
    reactant_num = toml_data.reactant_num()
    product_list_num = toml_data.product_list_num()

    # Create a graph
    name_tree = gp.Digraph()

    name_tree.node_attr.update(fontname='Arial')

    # Add edges
    gp_edges = []
    for init, fin in graph_table_num:
        list = [str(init), str(fin)]
        gp_edges.append(tuple(list))
    for edge in gp_edges:
        name_tree.edge(*edge, splines ='line')
    print(gp_edges)

    # Add nodes
    gp_nodes = []
    for i in range(len(state_list_name)):
        list = [str(state_list_num[i]), state_list_name[i]]
        gp_nodes.append(tuple(list))
    edges_num = {item for sublist in gp_edges for item in sublist}
    remove_node = [str(item[0]) for item in gp_nodes if str(item[0]) in edges_num]
    new_nodes = []
    for i in remove_node:
        for j in range(len(gp_nodes)):
            if gp_nodes[j][0] == i:
                new_nodes.append(gp_nodes[j])


    for node, label in new_nodes:
        if int(node) == reactant_num:
            name_tree.node(node, label, shape = 'box', style = 'filled', fillcolor = 'slategray1')
        elif int(node) in product_list_num:
            name_tree.node(node, label, shape = 'box', style = 'filled', fillcolor = 'lavender' )
        else:
            name_tree.node(node, label, shape = 'box', style = 'filled', fillcolor = 'lightyellow')
    print(gp_nodes)

    name_tree.render('example_graph', format='png', view=True)

graph_builder()