"""A tool to read toml file for the rate constant calculation"""

from __future__ import annotations

import os

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

    def total_atoms(self) -> int:
        """Extract the number of total atoms in toml file.

        Returns:
            int: number of total atoms
        """
        return self.data["molecule"]["total_atoms"]

    def normal_modes(self) -> int:
        """Calculate the normal modes from the number of total atoms.

        .. math::
            normal mode = \\3 \\cdot \\N \\- \\6

        Returns:
            int: number of normal modes
        """
        normal_modes = self._total_atoms() * 3 - 6
        return normal_modes

    def num_states(self) -> int:
        """Get the total number of states of the modeling.

        Args:
            num (int): numbering of the searching state

        Returns:
            int: total number of states
        """
        return len(self.data["state"])

    def mult(self, num: int) -> int:
        """Extract the spin multiplicity

        If the searching state is intersection, extract the lower one.

        Args:
            num (int): numbering of the searching state

        Returns:
            int: spin multiplicity
        """
        return self.data["state"][str(num)]["spin_multiplicity"]

    def target_spin_state(self, num: int) -> int:
        """Extract the target spin state.

        Args:
            num (int): numbering of the searching state

        Returns:
            int: target spin state
        """
        state = self.data["state"][str(num)]["target_spin_state"]
        if type(state) == int:
            return state
        if type(state) == list:
            return state[1]

    def state_name(self, num: int) -> str:
        """Extract the state name.

        Args:
            num (int): numbering of the searching state

        Returns:
            str: state name
        """
        return self.data["state"][str(num)]["name"]

    def state_num(self, num: int) -> int:
        """Extract the numbering of the state.

        Args:
            num (int): numbering of the searching state

        Returns:
            int: numbering of the state
        """
        return num

    def conc(self, num: int) -> float:
        """Extract the starting concentration.

        Args:
            num (int): numbering of the searching state

        Returns:
            float: starting concentration
        """
        return self.data["state"][str(num)]["conc"]

    def initial_name(self, num: int) -> str:
        """Extract the name of the initial state of the reaction.

        Args:
            num (int): numbering of the searching state

        Returns:
            str: name of initial state
        """
        return self.state_name(num)

    def initial_num(self, num: int) -> int:
        """Extract the numbering of the initial state of the reaction.

        Args:
            num (int): numbering of the searching state

        Returns:
            int: numbering of initial state
        """
        return self.state_num(num)

    def final_existence(self, init: int, fin: int) -> bool:
        """Check the existence of the reaction without transition state.

        Args:
            init (int): numbering of the initial state
            fin (int): numbering of the final state

        Returns:
            bool: True if exists otherwise False
        """
        if "final" in self.data["state"][str(init)]:
            if str(fin) in self.data["state"][str(init)]["final"]:
                return True
        else:
            return False

    def final_name(self, init: int, fin: int) -> str:
        """Extract the name of the final state of a reaction without transition state.

        final_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (int): numbering of the initial state
            fin (int): numbering of the final state

        Returns:
            str: name of the final state
        """
        assert self.final_existence(init, fin) == True
        return self.state_name(fin)

    def final_num(self, init: int, fin: int) -> int:
        """Extract the numbering of the final state of a reaction without transition state.

        final_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (int): numbering of the initial state
            fin (int): numbering of the final state

        Returns:
            int: numbering of the final state
        """
        assert self.final_existence(init, fin) == True
        return self.state_num(fin)

    def ts_existence(self, init: int, ts: int):
        """Check the existence of the reaction with transition state.

        Args:
            init (int): numbering of the initial state
            ts (int): numbering of the transition state

        Returns:
            bool: True if exists otherwise False
        """
        if "ts" in self.data["state"][str(init)]:
            if str(ts) in self.data["state"][str(init)]["ts"]:
                return True

    def ts_name(self, init: int, ts: int) -> str:
        """Extract the name of the transition state.

        ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (int): numbering of the initial state
            ts (int): numbering of the transition state

        Returns:
            str: name of the transition state
        """
        assert self.ts_existence(init, ts) == True
        return self.state_name(ts)

    def ts_num(self, init: int, ts: int) -> int:
        """Extract the numbering of the final state of a reaction with transition state.

        ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (int): numbering of the initial state
            ts (int): numbering of the transition state

        Returns:
            int: numbering of the final state
        """
        assert self.ts_existence(init, ts) == True
        return self.state_num(ts)

    def ts_final_name(self, init: int, ts: int) -> str:
        """Extract the name of the final state of a reaction with transition state.

        ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (int): numbering of the initial state
            ts (int): numbering of the transition state

        Returns:
            int: numbering of the final state
        """
        assert self.ts_existence(init, ts) == True
        state = int(self.data["state"][str(init)]["ts"][str(ts)]["final"])
        return self.state_name(state)

    def ts_final_num(self, init: int, ts: int) -> int:
        """Extract the numbering of the final state in a reaction with transition state.

        ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (int): numbering of the initial state
            ts (int): numbering of the transition state

        Returns:
            int: numbering of the final state
        """
        assert self.ts_existence(init, ts) == True
        state = int(self.data["state"][str(init)]["ts"][str(ts)]["final"])
        return self.state_num(state)

    def graph_edge(self, init: int, fin: int) -> tuple:
        """Extract the numbering of the final state of a reaction with transition state.

        ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init (int): numbering of the initial/transition state
            fin (int): numbering of the final/transition state

        Returns:
            tuple: connection of initial and final state of reaction (initial, final)
        """
        edge = tuple([init, fin])
        return edge

    def reaction_type(self, init: int, fin: int) -> str:
        """Check the reaction type of the reaciontion without transition state.

        final_existence should be true, otherwise Assertionerror will be raised.

        Args:
            init (int): numbering of the initial state
            fin (int): numbering of the final state

        Returns:
            str: type of the reaction
        """
        if self.final_existence(init, fin):
            return self.data["state"][str(init)]["final"][str(fin)]["reaction_type"]

    def file_path(self, num: int, calculation_path: str):
        """Get the file path of the given state.

        For the optimized(minimized) states, using the corresponding folder of the given state name.
        For the Frank-Condon points (vertical excitations), using the result of corresponding optimized state.

        Args:
            num (int): numbering of the searching state
            calculation_path (str): absolute path of calculation directory (calculation folders should have same parallel structure here)

        Returns:
            str: file path
        """
        for folder in os.listdir(calculation_path):
            state_name = self.state_name(num)
            if state_name.endswith("*"):
                target_folder_name = state_name[:-1]
            else:
                target_folder_name = state_name

            if folder.endswith(target_folder_name):
                file_path = os.path.join(calculation_path, folder, "sp", "tc.out")
        return file_path

    def start_conc(self) -> list:
        """Generate a list with initial concentration of each state.

        Returns:
            list: initial concentrations of each state
        """
        start_conc = np.zeros(self.num_states())
        for i in range(self.num_states()):
            start_conc[i] = self.conc(i)
        return start_conc

    def state_list_name(self) -> list:
        """Generate a list with a name of each state.

        Returns:
            list: name of each state
        """
        state_list_name = ["name"] * self.num_states()
        for i in range(self.num_states()):
            init_name = self.state_name(i)
            state_list_name[i] = init_name
        return state_list_name

    def state_list_num(self) -> list:
        """Generate a list with name of each state.

        Returns:
            list: numbering of each state
        """
        state_list_num = np.zeros(self.num_states(), dtype=int)
        for i in range(self.num_states()):
            init_num = self.state_num(i)
            state_list_num[i] = init_num
        return state_list_num
