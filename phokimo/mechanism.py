"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import os
import sys
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
    """A tool to build graph diagram that represents the reaction mechanism.

    Returns:
        list: list of tuples that represent graph edges (initial/transition state, transition/final state)
    """

    # INPUT
    if len(sys.argv) < 2:
        print("Usage: python3 -m __main__.py <name_of_toml_file>")
        sys.exit(1)

    toml_file = sys.argv[1]  # Get the TOML file from the command line
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Assume the TOML file is in the same directory as the __main__.py
    # Modify toml_file_path for suitable set-up toml file
    toml_file_path = os.path.join(current_dir, toml_file)

    toml_data = TomlReader(toml_file_path)

    state_list_name = toml_data.visualize_state_list_name()
    state_list_num = toml_data.state_list_num()
    graph_table_num = remove_duplicates(toml_data.graph_table_num())
    reactant_nums = toml_data.reactant_nums()
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
        if int(node) in reactant_nums:
            name_tree.node(node, label, shape = 'box', style = 'filled', fillcolor = 'slategray1')
        elif int(node) in product_list_num:
            name_tree.node(node, label, shape = 'box', style = 'filled', fillcolor = 'lavender' )
        else:
            name_tree.node(node, label, shape = 'box', style = 'filled', fillcolor = 'lightyellow')
    print(gp_nodes)

    name_tree.render('example_graph', format='png', view=True)

graph_builder()