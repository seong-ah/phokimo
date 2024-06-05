""" Tool to generate graph about reaction connection. """
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

def graph_builder(state_list_name, state_list_num, graph_table_name, graph_table_num):
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

    return(table)
