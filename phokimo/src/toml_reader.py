"""A tool to read toml file for the rate constant calculation"""

import toml
import os
import numpy as np

from phokimo.src.io.terachem import TeraChemOutputReader

from phokimo.src.rate_constants import RateCalculator

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

        with open(self.fpath, 'r') as file:
            self.data = toml.load(file)

    def _total_atoms(self) -> int:
        """Extract the number of total atoms in toml file.

        Returns:
            int: number of total atoms
        """
        return self.data['molecule']['total_atoms']
    
    def _normal_modes(self) -> int:
        """Calculate the normal modes from the number of total atoms.

        .. math::
            normal mode = \\3 \\cdot \\N \\- \\6

        Returns:
            int: number of normal modes
        """
        normal_modes = self._total_atoms() * 3 - 6
        return normal_modes
    
    def _mult(self, num: int) -> int:
        """Extract the spin multiplicity

        If the searching state is intersection, extract the lower one.

        Args:
            num: numbering of the searching state

        Returns:
            int: spin multiplicity
        """
        return self.data['state'][str(num)]['spin_multiplicity']
    
    def _target_spin_state(self, num: int) -> int:
        """Extract the target spin state.

        Args:
            num: numbering of the searching state

        Returns:
            int: target spin state
        """
        state = self.data['state'][str(num)]['target_spin_state']
        if type(state) == int:
            return state
        if type(state) == list:
            return state[1]

    def _state_name(self, num: int) -> str:
        """Extract the state name.

        Args:
            num: numbering of the searching state

        Returns:
            str: state name
        """
        return self.data['state'][str(num)]['name']
    
    def _state_num(self, num: int) -> int:
        """Extract the numbering of the state.

        Args:
            num: numbering of the searching state

        Returns:
            int: numbering of the state
        """
        return num
    
    def _conc(self, num: int) -> float:
        """Extract the starting concentration.

        Args:
            num: numbering of the searching state

        Returns:
            float: starting concentration
        """
        return self.data['state'][str(num)]['conc']
    
    def _initial_name(self, num: int) -> str:
        """Extract the name of the initial state of the reaction.

        Args:
            num: numbering of the searching state

        Returns:
            str: name of initial state
        """
        return self._state_name(num)
    
    def _initial_num(self, num: int) -> int:
        """Extract the numbering of the initial state of the reaction.

        Args:
            num: numbering of the searching state

        Returns:
            int: numbering of initial state
        """
        return self._state_num(num)
    
    def _final_existence(self, init: int, fin: int) -> bool:
        """Check the existence of the reaction without transition state.

        Args:
            init: numbering of the initial state
            fin: numbering of the final state

        Returns:
            bool: True if exists otherwise False
        """
        if 'final' in self.data['state'][str(init)]:
            if str(fin) in self.data['state'][str(init)]['final']:
                return True
        else:
            return False

    def _final_name(self, init: int, fin: int) -> str:
        """Extract the name of the final state of a reaction without transition state.

        _final_existence should be True, otherwise AssertionError will be raised.

        Args:
            init: numbering of the initial state
            fin: numbering of the final state

        Returns:
            str: name of the final state
        """
        assert self._final_existence(init, fin) == True
        return self._state_name(fin)

    def _final_num(self, init: int, fin: int) -> int:
        """Extract the numbering of the final state of a reaction without transition state.

        _final_existence should be True, otherwise AssertionError will be raised.

        Args:
            init: numbering of the initial state
            fin: numbering of the final state

        Returns:
            int: numbering of the final state
        """
        assert self._final_existence(init, fin) == True
        return self._state_num(fin)
        
    def _ts_existence(self, init: int, ts: int):
        """Check the existence of the reaction with transition state.

        Args:
            init: numbering of the initial state
            ts: numbering of the transition state

        Returns:
            bool: True if exists otherwise False
        """
        if 'ts' in self.data['state'][str(init)]:
            if str(ts) in self.data['state'][str(init)]['ts']:
                return True

    def _ts_name(self, init: int, ts: int) -> str:
        """Extract the name of the transition state.

        _ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init: numbering of the initial state
            ts: numbering of the transition state

        Returns:
            str: name of the transition state
        """
        assert self._ts_existence(init, ts) == True
        return self._state_name(ts)
        
    def _ts_num(self, init: int, ts: int) -> int:
        """Extract the numbering of the final state of a reaction with transition state.

        _ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init: numbering of the initial state
            ts: numbering of the transition state

        Returns:
            int: numbering of the final state
        """
        assert self._ts_existence(init, ts) == True
        return self._state_num(ts)
    
    def _ts_final_name(self, init: int, ts: int) -> str:
        """Extract the name of the final state of a reaction with transition state.

        _ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init: numbering of the initial state
            ts: numbering of the transition state

        Returns:
            int: numbering of the final state
        """
        assert self._ts_existence(init, ts) == True
        state = int(self.data['state'][str(init)]['ts'][str(ts)]['final'])
        return self._state_name(state)
    
    def _ts_final_num(self, init: int, ts: int) -> int:
        """Extract the numbering of the final state in a reaction with transition state.

        _ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init: numbering of the initial state
            ts: numbering of the transition state

        Returns:
            int: numbering of the final state
        """
        assert self._ts_existence(init, ts) == True
        state = int(self.data['state'][str(init)]['ts'][str(ts)]['final'])
        return self._state_num(state)
    
    def _graph_edge(self, init, fin) -> list:
        """Extract the numbering of the final state of a reaction with transition state.

        _ts_existence should be True, otherwise AssertionError will be raised.

        Args:
            init: numbering of the initial state
            fin: numbering of the final state

        Returns:
            int: numbering of the final state
        """
        edge = tuple([init, fin])
        return edge
    
    def _reaction_type(self, init, fin) -> str:
        """Check the reaction type of the reaciontion without transition state.

        _final_existence should be true, otherwise Assertionerror will be raised.

        Args:
            init: numbering of the initial state
            fin: numbering of the final state

        Returns:
            str: type of the reaction
        """
        if self._final_existence(init, fin) == True:
            return self.data['state'][str(init)]['final'][str(fin)]['reaction_type']
