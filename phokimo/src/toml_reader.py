"""A tool to read toml file for the rate constant calculation"""

from __future__ import annotations

import os
import networkx as nx
from collections import defaultdict

import numpy as np
import toml


class TomlReader:
    """Extract information from the toml file."""

    def __init__(self, fpath: str) -> None:
        """Extract information from toml file.

        Reads generated toml file.
        If not exist, AssertionError will be raised.

        Args:
            fpath (_type_): absolute file path of the output file.
        """
        self.fpath = fpath
        assert fpath.endswith(".toml")

        self.data = None

        with open(self.fpath) as file:
            self.data = toml.load(file)

        self.calculation_dir = self.calculation_path()

        self.len_states = self.num_states()

        self.name_to_num = self.state_name_num_dict()
        self.reverse_name_to_num = self.reverse_state_name_num_dict()
        self.vis_name_list = self.visualize_state_list_name()
        self.ref_name = self.reference_state()

        self.reactant_names = self.reactant_list_name()
        self.reactant_nums = self.reactant_list_num()
        self.product_names = self.product_list_name()
        self.product_nums = self.product_list_num()

        self.spins = self.spin_list()
        self.spin_dict = self.spin_list_dict()

        self.reactions_name = self.reaction_list_name()
        self.reactions_num = self.reaction_list_num()

        self.teq_graph = self.graph_teq_group()

        if 'substate' in self.data:
            self.sub_name_to_num = self.substate_name_num_dict()
            self.reverse_sub_name_to_num = self.reverse_substate_name_num_dict()
        
        self.file_path = self.file_path_dict()


    def total_atoms(self) -> int:
        """Extract the number of total atoms in toml file.

        Returns:
            int: number of total atoms
        """
        return self.data["molecule"]["total_atoms"]

    def total_normal_modes(self) -> int:
        """Calculate the normal modes from the number of total atoms.

        .. math::
            normal mode = \\3 \\cdot \\N \\- \\6

        Returns:
            int: number of normal modes
        """
        normal_modes = self._total_atoms() * 3 - 6
        return normal_modes
    
    def duration(self) -> float:
        """Get the duration of the modling. In toml file, time is written as fs

        Returns:
            float: duration of the modeling (s)
        """
        return float(self.data["molecule"]["duration"]) * 10 ** (-15) #fs

    def num_states(self) -> int:
        """Get the total number of states of the modeling.

        Returns:
            int: total number of states
        """
        return len(self.data["state"])

    def mult(self, state: str, substate: bool = False) -> int:
        """Extract the spin multiplicity

        If the searching state is intersection, extract the lower one.

        Args:
            state (str): name of the searching state
            substate (bool): default setting is False, True when targeting substate

        Returns:
            int: spin multiplicity
        """
        if substate:
            return self.data["substate"][state]["spin_multiplicity"]
        else:
            return self.data["state"][state]["spin_multiplicity"]

    def target_spin_state(self, state: str, substate: bool = False) -> int:
        """Extract the target spin state.

        Args:
            state (str): name of the searching state
            substate (bool): default setting is False, True when targeting substate

        Returns:
            int: target spin state
        """
        if not substate:
            spin_state = self.data["state"][state]["target_spin_state"]
        else:
            spin_state = self.data["substate"][state]["target_spin_state"]
            
        if type(spin_state) == int:
            return [spin_state, spin_state]
        else:
            return spin_state

    def state_name(self, state: str, substate: bool = False) -> str:
        """Extract the state name.

        Args:
            state (str): name of the searching state
            substate (bool): default setting is False, True when targeting substate

        Returns:
            str: state name
        """
        if substate == False:
            return self.data["state"][state]["name"]
        else:
            return self.data["substate"][state]["name"]
    
    def visualize_state_name(self, state: str) -> str:
        """Extract the state name for visualization.

        Args:
            state (str): name of the searching state

        Returns:
            str: state name for visualization
        """
        return self.data["state"][state]["visualize_name"]

    def state_name_num_dict(self) -> dict:
        """Generate number label of each state as dictionary format.

        Returns:
            dict: dictionary that connects state name and number label
        """
        name_num_dict = {}

        counter = 0

        for name in self.data['state']:

            if name not in name_num_dict:
                name_num_dict[name] = counter
                counter += 1

        return name_num_dict

    def reverse_state_name_num_dict(self) -> dict:
        """Reverse the name_to_num dict.

        Returns:
            dict: dictionary that connects number label to name
        """
        num_to_name = {num: name for name, num in self.name_to_num.items()}
        return num_to_name

    def substate_name_num_dict(self) -> dict:
        """Generate number label of each state as dictionary format for substates.
        substate should exist in TOML file, otherwise AssertionError will be raised.

        Returns:
            dict: dictionary that connects state name and number label
        """
        assert 'substate' in self.data

        name_num_dict = {}

        counter = 0

        for name in self.data["substate"]:
            if name not in name_num_dict:
                name_num_dict[name] = counter
                counter += 1

        return name_num_dict

    def reverse_substate_name_num_dict(self) -> dict:
        """Reverse the name_to_num dict of substates.

        Returns:
            dict: dictionary that connects number label to name
        """
        num_to_name = {num: name for name, num in self.sub_name_to_num.items()}
        return num_to_name

    def state_num(self, state: str) -> int:
        """Extract the numbering of the state.

        Args:
            state (str): name of the searching state

        Returns:
            int: numbering of the state
        """
        return self.name_to_num[state]

    def conc(self, state: str) -> float:
        """Extract the starting concentration.

        Args:
            state (str): name of the searching state

        Returns:
            float: starting concentration
        """
        return self.data["state"][state]["conc"]
    
    def normal_mode(self, init: str, final: str) -> float:
        """Extract the normal mode of the reaction.

        Args:
            init (str): name of the initial state
            fin (str): name of the final state

        Returns:
            float: corresponding normal mode of the reaction
        """
        assert "final" in self.data["state"][init]
        return self.data["state"][init]["final"][final]["normal_mode"]
    
    def substate_existence(self, state: str) -> bool:
        """Check substate exists or not.

        Args:
            state (str): name of the searching state

        Returns:
            bool: True if substate exists and False for else
        """
        if 'substate' in self.data["state"][state]:
            return True
        else:
            return False
        
    def substate_list(self, state: str) -> list:
        """Return list of substates.

        substate_existence should be True, otherwise AssertionError will be raised.
        
        Args:
            state (str): name of the searching state

        Returns:
            list: list of substates
        """
        assert self.substate_existence(state) 

        return self.data["state"][state]["substate"]
    
    def theory_level(self, state: str, substate: bool = False, ground: bool = False) -> str:
        """Return theory level of corresponding state.

        Args:
            state (str): name of the searching state
            substate (bool): default setting is False, True when targeting substate
            ground (bool): default setting is excited, True when theory level is ground state calculation

        Returns:
            str: theory level (excited or ground)
        """
        if state == self.reference_state():
            if ground == False:
                return "excited"
            else:
                return "ground"
        else:
            if substate:
                return self.data["substate"][state]["theory_level"]
            else:
                return self.data["state"][state]["theory_level"]
            
    def final_existence(self, init: str, fin: str) -> bool:
        """Check the existence of the reaction without transition state (relaxation).

        Args:
            init (str): name of the initial state
            fin (str): name of the final state

        Returns:
            bool: True if exists otherwise False
        """
        if "final" in self.data["state"][init]:
            if fin in self.data["state"][init]["final"]:
                return True
        else:
            return False

    def final_name(self, init: str, fin: str) -> str:
        """Extract the name of the final state of a reaction without transition state (relaxation).

        final_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (str): name of the initial state
            fin (str): name of the final state

        Returns:
            str: name of the final state
        """
        assert self.final_existence(init, fin)
        return self.state_name(fin)

    def final_num(self, init: str, fin: str) -> int:
        """Extract the numbering of the final state of a reaction without transition state (relaxation).

        final_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (str): name of the initial state
            fin (str): name of the final state

        Returns:
            int: numbering of the final state
        """
        assert self.final_existence(init, fin)
        return self.state_num(fin)

    def ts_existence(self, init: str, ts: str) -> bool:
        """Check the existence of the reaction with transition state.

        Args:
            init (str): name of the initial state
            ts (str): name of the transition state

        Returns:
            bool: True if exists otherwise False
        """
        if "ts" in self.data["state"][init]:
            if ts in self.data["state"][init]["ts"]:
                return True

    def ts_name(self, init: str, ts: str) -> str:
        """Extract the name of the transition state.

        ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (str): name of the initial state
            ts (str): name of the transition state

        Returns:
            str: name of the transition state
        """
        assert self.ts_existence(init, ts) 
        return self.state_name(ts)

    def ts_num(self, init: str, ts: str) -> int:
        """Extract the numbering of the transition state.

        ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (str): name of the initial state
            ts (str): name of the transition state

        Returns:
            int: numbering of the transition state
        """
        assert self.ts_existence(init, ts)
        return self.state_num(ts)

    def ts_final_name(self, init: str, ts: str) -> str:
        """Extract the name of the final state of a reaction with transition state.

        ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (str): name of the initial state
            ts (str): name of the transition state

        Returns:
            int: numbering of the final state of a reaction with transition state
        """
        assert self.ts_existence(init, ts) 
        final = int(self.data["state"][init]["ts"][ts]["final"])
        return self.state_name(final)

    def ts_final_num(self, init: str, ts: str) -> int:
        """Extract the numbering of the final state in a reaction with transition state.

        ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (str): name of the initial state
            ts (str): name of the transition state

        Returns:
            int: numbering of the final state
        """
        assert self.ts_existence(init, ts) 
        final = int(self.data["state"][init]["ts"][ts]["final"])
        return self.state_num(final)

    def graph_edge(self, init, fin) -> tuple:
        """Generate edge of the mechanism graph.

        Args:
            init (int of str): numbering or name of the initial state
            fin (int or str): numbering or name of the transition state

        Returns:
            tuple: edge connection of initial and final state of reaction (initial, final)
        """
        edge = tuple([init, fin])
        return edge

    def reaction_type(self, init: str, fin: str) -> str:
        """Check the reaction type of the reaction without transition state (relaxation).

        final_existence should be true, otherwise Assertionerror will be raised.

        Args:
            init (str): name of the initial state
            fin (str): name of the final state

        Returns:
            str: type of the reaction
        """
        if self.final_existence(init, fin):
            return self.data["state"][init]["final"][fin]["reaction_type"]
        
    def calculation_path(self) -> str:
        """Get the calculation path of root directory.

        Returns:
            str: root directory of calculation directories
        """
        return self.data["molecule"]["calculation_path"]

    def file_path_dict(self) -> dict:
        """Generate a dictionary of file paths
        
        if there are multiple substates for a state, dict[state] = dict of file paths
        if no substates, dict[state] = file path

        Returns:
            dict: dict of file path except reference state
        """
        file_path_dict = {}
        
        for state in self.name_to_num:
            sub_dict = {}
            if "substate" in self.data["state"][state]:
                for substate in self.data["state"][state]["substate"]:
                    for folder in os.listdir(self.calculation_dir):
                        folder_name = self.data["substate"][substate]["folder_name"]
                        if folder.endswith(folder_name):
                            if self.theory_level(state, substate = True) == "ground":
                                sub_dict[substate] = os.path.join(self.calculation_dir, folder, "ground_sp", "tc.out")
                            else:
                                sub_dict[substate] = os.path.join(self.calculation_dir, folder, "sp", "tc.out")
                    file_path_dict[state] = sub_dict
            else:
                folder_name = self.data["state"][state]["folder_name"]
                for folder in os.listdir(self.calculation_dir):
                    if folder.endswith(folder_name):
                        if self.theory_level(state) == "ground":
                            file_path = os.path.join(self.calculation_dir, folder, "ground_sp", "tc.out")
                        else:
                            file_path = os.path.join(self.calculation_dir, folder, "sp", "tc.out")
                        file_path_dict[state] = file_path
        return file_path_dict

    def ref_file_path(self, ground = False) -> str:
        """Get the file path of the reference state.

        Args:
            ground (bool): default setting is excited, True when theory level is ground state calculation

        Returns:
            str: directory of reference state output calculation file
        """
        folder_name = self.data["state"][self.ref_name]["folder_name"]
        for folder in os.listdir(self.calculation_dir):
            if folder.endswith(folder_name):
                if ground:
                    return os.path.join(self.calculation_dir, folder, "ground_sp", "tc.out")
                else:
                    return os.path.join(self.calculation_dir, folder, "sp", "tc.out")
        
    def initial_conc(self) -> np.ndarray:
        """Generate a list with initial concentration of each state.

        Returns:
            np.ndarray: initial concentration of each state
        """
        start_conc = np.zeros(self.len_states)
        for state in self.name_to_num:
            state_num = self.state_num(state)
            start_conc[state_num] = self.conc(state)
        return start_conc
    
    def visualize_state_list_name(self) -> list:
        """Generate a list with a name of each state for visualization.

        Returns:
            list: name of each state for visualization
        """
        visualize_state_list_name = ["name"] * self.len_states
        for state, num in self.name_to_num.items():
            name = self.visualize_state_name(state)
            visualize_state_list_name[num] = name
        return visualize_state_list_name
    
    def reaction_list_name(self) -> list[tuple]:
        """Generate a list of reaction linkage.

        linkage of init -> final, ignoring transition state

        Returns:
            list: list of tuples of initial state(name) and final state(name) (initial, final)
        """
        reaction_list = []
        for init in self.name_to_num:
            for next in self.name_to_num:
                if self.ts_existence(init, next):
                    ts_final_name = self.ts_final_name(init, next)
                    edge = (init, ts_final_name)
                    reaction_list.append(edge)
                if self.final_existence(init, next):
                    edge = (init, next)
                    reaction_list.append(edge)
        return reaction_list

    def reaction_list_num(self) -> list[tuple]:
        """Generate a list of reaction linkage with numbering.

        linkage of init(int) -> final(int), ignoring transition state

        Returns:
            list: list of tuples of initial state(int) and final state(int) (initial, final)
        """
        reaction_list_num = []
        for edge in self.reactions_name:
            init, fin = edge
            num_edge = (self.state_num(init), self.state_num(fin))
            reaction_list_num.append(num_edge)
        return reaction_list_num
        
    
    def reaction_types(self) -> list:
        """Generate a list of reaction types.

        Each number represents the reaction type of reaction_types[initial_state][final_state]:
            1: reaction with transition state
            2: non-radiative relaxation
            3: radiative emission

        Returns:
            list: reaction type of corresponding reactions (reaction_types[initial state][final state])
        """
        dim = (self.len_states, self.len_states)
        reaction_types = np.zeros(dim)
        for i in range(self.len_states):
            init_num = self.state_num(i)
            for j in range(self.len_states):
                if self.ts_existence(i, j):
                    ts_final_num = self.ts_final_num(i, j)
                    reaction_types[init_num][ts_final_num] = 1
                if self.final_existence(i, j):
                    final_num = self.final_num(i, j)
                    if self.reaction_type(init_num, final_num) == "relaxation":
                        reaction_types[init_num][final_num] = 2
                    if self.reaction_type(init_num, final_num) == "emission":
                        reaction_types[init_num][final_num] = 3
        return reaction_types
    
    def graph_table_name(self) -> list[tuple]:
        """Generate a list with reaction connection with state names.

        Returns:
            list: graph edge tuples that represent the reaction connections with state names
        """
        graph_table_name = []
        for init in self.name_to_num:
            for next in self.name_to_num:
                if self.ts_existence(init, next):
                    ts_name = next
                    ts_final_name = self.ts_final_name(init, next)
                    graph_table_name.append(self.graph_edge(init, ts_name))
                    graph_table_name.append(self.graph_edge(ts_name, ts_final_name))
                if self.final_existence(init, next):
                    final_name = next
                    graph_table_name.append(self.graph_edge(init, final_name))
        return graph_table_name

    def graph_table_num(self) -> list:
        """Generate a list with reaction connection with state numberings.

        Returns:
            list: graph edge tuples that represent the reaction connections with state numberings
        """
        graph_table_num = []
        for init in self.name_to_num:
            init_num = self.state_num(init)
            for next in self.name_to_num:
                if self.ts_existence(init, next):
                    ts_num = self.ts_num(init, next)
                    ts_final_num = self.ts_final_num(init, next)
                    graph_table_num.append(self.graph_edge(init_num, ts_num))
                    graph_table_num.append(self.graph_edge(ts_num, ts_final_num))
                if self.final_existence(init, next):
                    final_num = self.final_num(init, next)
                    graph_table_num.append(self.graph_edge(init_num, final_num))
        return graph_table_num
    
    def _condition(self, state = str) -> str:
        """Extract the condition of corresponding state.

        Args:
            state (str): name of the searching state
        
        Returns:
            str: condition of each state (e.g. reactant, product, intermediate, etc)
        """
        return self.data["state"][state]["condition"]
    
    def reactant_list_name(self) -> list[str]: 
        """Extract the reactant names.

        Returns:
            list[str]: list of reactant names
        """
        reactant_list_name = []
        for state in self.name_to_num:
            if self._condition(state) == "reactant":
                reactant_list_name.append(state)
        return reactant_list_name
            
    def reactant_list_num(self) -> list[int]:
        """Extract the numbering of the reactants.

        Returns:
            list[int]: list of the numbering of the reactants
        """
        return [self.state_num(state) for state in self.reactant_names]

    def product_list_name(self) -> list[str]:
        """Generate a list of name of the products.

        Returns:
            list[int]: list of name of the products
        """
        product_list_name = []
        for state in self.name_to_num:
            if self._condition(state) == "product":
                product_list_name.append(state)
        return product_list_name

    def product_list_num(self) -> list[int]:
        """Generate a list of numbering of the products.

        Returns:
            list[int]: list of numbering of the products
        """
        return [self.state_num(state) for state in self.product_names]
    
    def reference_state(self) -> str:
        """Get the name of reference energy state.

        Returns:
            str: name of reference energy state
        """
        for state in self.name_to_num:
            if "reference_state" in self.data["state"][state]:
                return state
            
    def spin_list(self) -> list:
        """Get the list of existing spin states.

        Returns:
            list[tuple]: list of tuples that represents spin states as (spin, target) (e.g. S0 = (1, 0), S1 = (1, 1), and T2 = (3, 2))
        """
        spin_list = []
        for state in self.name_to_num:
            if type(self.data["state"][state]["target_spin_state"]) != list:
                spin = (self.mult(state), self.target_spin_state(state)[0])
                if spin not in spin_list:
                    spin_list.append(spin)
        return spin_list
    
    def spin_list_dict(self) -> dict:
        """Map the spins and corresponding states.

        Returns:
            dict: map spin (tuple) to corresponding states as list items
        """
        spin_list_dict = {}
        for spin in self.spins:
            spin_states = []
            mult, target = spin
            for state in self.name_to_num:
                if self.data["state"][state]["spin_multiplicity"] == mult and self.data["state"][state]["target_spin_state"] == target:
                    spin_states.append(state)
            spin_list_dict[spin] = spin_states
        return spin_list_dict
    
    def graph_teq_group(self) -> dict:
        """Group the reactions that has same equilibriated temperature.

        Returns:
            dict: map reactant to parent state of T_eq and parent state to corresponding states as list items
        """       
        def get_descendants(G, node):
            descendants = []
            for child in G.successors(node):
                descendants.append((node, child))
                descendants.extend(get_descendants(G, child))
            return descendants
        
        roots = self.reactant_names
        whole_groups = {}

        for root in roots:
            G = nx.DiGraph()
            edges = self.reactions_name
            G.add_edges_from(edges)

            parent_groups = defaultdict(list)
            direct_children = list(G.successors(root))

            for parent in direct_children:
                descendants = get_descendants(G, parent)
                if descendants:
                    parent_groups[parent].extend(descendants)
            whole_groups[root] = parent_groups
        return whole_groups
