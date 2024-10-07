"""Tool to generate graph about reaction connection."""

from __future__ import annotations

import os
import sys
import numpy as np
import graphviz as gp
from phokimo.src.toml_reader import TomlReader
from collections import Counter

def main():
    """ A tool to build graph diagram that represents the reaction mechanism """
    
    if len(sys.argv) < 2:
        print("Usage: python3 -m phokimo filename.toml")
        sys.exit(1)

    toml_file = sys.argv[1]
    current_dir = os.getcwd()
    toml_file_path = os.path.join(current_dir, toml_file)
    toml_data = TomlReader(toml_file_path)

    name_to_num = toml_data.name_to_num
    visualize_name = toml_data.vis_name_list
    graph_table_name = toml_data.graph_name
    reactant_name = toml_data.reactant_names
    product_list_name = toml_data.product_names

    name_tree = gp.Digraph()
    name_tree.node_attr.update(fontname='Arial')

    gp_edges = []
    for init, fin in graph_table_name:
        init_visualize = visualize_name[name_to_num[init]]
        fin_visualize = visualize_name[name_to_num[fin]]
        reaction = (init_visualize, fin_visualize)
        gp_edges.append(reaction)
    for edge in gp_edges:
        name_tree.edge(*edge, splines ='line')

    nodes = {n for edge in gp_edges for n in edge}

    for node in name_to_num:
        if visualize_name[name_to_num[node]] in nodes:
            if node in reactant_name:
                name_tree.node(visualize_name[name_to_num[node]], shape = 'box', style = 'filled', fillcolor = 'slategray1')
            elif node in product_list_name:
                name_tree.node(visualize_name[name_to_num[node]], shape = 'box', style = 'filled', fillcolor = 'lavender' )
            else:
                name_tree.node(visualize_name[name_to_num[node]], shape = 'box', style = 'filled', fillcolor = 'lightyellow')
    name_tree.render('mechanism_diagram.png', format='png', view=True)

if __name__ == "__main__":
    main()