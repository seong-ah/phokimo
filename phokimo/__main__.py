"""The main entry point for the application."""

from __future__ import annotations

from functools import partial

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import toml
import os

from scipy.integrate import odeint

from phokimo.src.io.terachem import TeraChemOutputReader

from phokimo.src.rate_constants import RateCalculator, ReactionTheory, EyringEquation, RelaxationTheory, AdhocRelaxation

from phokimo.src.toml_reader import TomlReader

from phokimo.src.ode_builder import general_ode, construct_ode

def main() -> None:
    """Run the application."""

    """ Read data from toml file and set-up ode calculation """

    # Get the directory of the currently executing Python script
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Assume the TOML file is in the same directory as the Python script
    toml_file_path = os.path.join(current_dir, "sample.toml")

    toml_data = TomlReader(toml_file_path)

    # Generate start_conc and state_list from toml
    total_atoms = float(toml_data._total_atoms())
    normal_modes = float(toml_data._normal_modes())

    num_states = len(toml_data.data['state'])
    start_conc = np.zeros(num_states)

    state_list_name = ["name"] * num_states
    state_list_num = np.zeros(num_states, dtype=int)
    state_list_energy = np.zeros(num_states)
    state_list_oscil = np.zeros(num_states)

    dim = (num_states, num_states)
    rates = np.zeros(dim)

    graph_table_name = []
    graph_table_num = []

    for i in range(num_states):
        state = str(i)
        start_conc[i] = toml_data._conc(i)

        init_name = toml_data._state_name(i)
        state_list_name[i] = init_name
        init_num = toml_data._state_num(i)
        state_list_num[i] = init_num

        if init_name.endswith('*'):
            folder_name = init_name[:-1]
            file_path = os.path.join(current_dir, folder_name, 'tc.out')
        else:
            file_path = os.path.join(current_dir, init_name, 'tc.out')

        terachem_file = TeraChemOutputReader(file_path)
        result = terachem_file.ci_energy(max_roots = 2) #Assume that only singlet states
        target_spin_state = toml_data._target_spin_state(i)
        hartree_energy = result[0][target_spin_state]
        energy = (2625.5 * (10 ** 3)) * hartree_energy #J/mol
        state_list_energy[i] = energy
        if target_spin_state != 0:
            oscilstr = result[0][target_spin_state - 1]
            state_list_oscil[i] = oscilstr

    for i in range(num_states):            
        for j in range(num_states):
            init_name = toml_data._state_name(i)
            init_num = toml_data._state_num(i)

            rate_formula = RateCalculator()

            if toml_data._ts_existence(i, j) == True:
                ts_name = toml_data._ts_name(i, j)
                ts_num = toml_data._ts_num(i, j)
                ts_final_name = toml_data._ts_final_name(i, j)
                ts_final_num = toml_data._ts_final_num(i, j)

                graph_table_name.append(toml_data._graph_edge(init_name, ts_name))
                graph_table_num.append(toml_data._graph_edge(init_num, ts_num))
                graph_table_name.append(toml_data._graph_edge(ts_name, ts_final_name))
                graph_table_num.append(toml_data._graph_edge(ts_num, ts_final_num))
                
                rate_formula = RateCalculator()
                dE = state_list_energy[ts_num] - state_list_energy[init_num]
                rate_constant = rate_formula.reaction_theory.compute_rate(dE)
                rates[init_num][ts_final_num] = rate_constant

            if toml_data._final_existence(i, j) == True:
                final_name = toml_data._final_name(i, j)
                final_num = toml_data._final_num(i, j)

                graph_table_name.append(toml_data._graph_edge(init_name, final_name))
                graph_table_num.append(toml_data._graph_edge(init_num, final_num))

                if toml_data._reaction_type(init_num, final_num) == 'relaxation':
                    dE = state_list_energy[final_num] - state_list_energy[init_num]
                    rate_constant = rate_formula.relaxation_theory.compute_rate(dE, normal_modes, total_atoms)
                    rates[init_num][final_num] = rate_constant

    time = np.linspace(0, 10, 1000)

    #Graph creation

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
    nx.draw(name_graph, with_labels=True, font_weight='bold')
    plt.show()

    nx.draw(num_graph, with_labels=True, font_weight='bold')
    plt.show()

    # Solving ode

    func = partial(construct_ode, table=table, rates=rates)
    conc = odeint(func, start_conc, time)

    plt.plot(time, conc)
    plt.show()

if __name__ == "__main__":
    main()
