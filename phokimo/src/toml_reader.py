"""A tool to read toml file for the rate constant calculation"""

import toml
import os
import numpy as np

from phokimo.src.io.terachem import TeraChemOutputReader

from phokimo.src.rate_constants import RateCalculator

class TomlReader:
    """Extract information from the toml file."""

    def __init__(self, fpath: str) -> None:
        self.fpath = fpath
        assert fpath.endswith(".toml")

        self.data = None

        try:
            with open(self.fpath, 'r') as file:
                self.data = toml.load(file)
        except FileNotFoundError:
            print(f"File '{self.fpath}' not found.")
        except Exception as e:
            print(f"An error occurred while reading '{self.fpath}': {e}")

    def _total_atoms(self) -> int:
        return self.data['molecule']['total_atoms']
    
    def _normal_modes(self) -> int:
        normal_modes = self._total_atoms() * 3 - 6
        return normal_modes
    
    def _mult(self, num: int) -> int:
        return self.data['state'][str(num)]['spin_multiplicity']
    
    def _target_spin_state(self, num: int) -> int:
        state = self.data['state'][str(num)]['target_spin_state']
        if type(state) == int:
            return state
        if type(state) == list:
            return state[0]

    def _conical(self, num: int):
        if self.data['state'][num]['type'] == 'conical':
            return 'conical'
        else:
            return 'Not_conical'

    def _state_name(self, num: int) -> str:
        return self.data['state'][str(num)]['name']
    
    def _state_num(self, num: int) -> int:
        return num
    
    def _conc(self, num: int) -> float:
        return self.data['state'][str(num)]['conc']
    
    def _initial_name(self, num: int) -> str:
        return self._state_name(num)
    
    def _initial_num(self, num: int) -> int:
        return self._state_num(num)
    
    def _final_existence(self, init: int, fin: int):
        if 'final' in self.data['state'][str(init)]:
            if str(fin) in self.data['state'][str(init)]['final']:
                return True

    def _final_name(self, init: int, fin: int) -> str:
        assert self._final_existence(init, fin) == True
        return self._state_name(fin)

    def _final_num(self, init: int, fin: int) -> int:
        assert self._final_existence(init, fin) == True
        return self._state_num(fin)
        
    def _ts_existence(self, init: int, ts: int):
        if 'ts' in self.data['state'][str(init)]:
            if str(ts) in self.data['state'][str(init)]['ts']:
                return True

    def _ts_name(self, init: int, ts: int) -> str:
        assert self._ts_existence(init, ts) == True
        return self._state_name(ts)
        
    def _ts_num(self, init: int, ts: int) -> int:
        assert self._ts_existence(init, ts) == True
        return self._state_num(ts)
    
    def _ts_final_name(self, init: int, ts: int) -> str:
        assert self._ts_existence(init, ts) == True
        state = int(self.data['state'][str(init)]['ts'][str(ts)]['final'])
        return self._state_name(state)
    
    def _ts_final_num(self, init: int, ts: int) -> int:
        assert self._ts_existence(init, ts) == True
        state = int(self.data['state'][str(init)]['ts'][str(ts)]['final'])
        return self._state_num(state)
    
    def _graph_edge(self, init, fin) -> list:
        edge = tuple([init, fin])
        return edge
    
    def _reaction_type(self, init, state) -> str:
        #if self._ts_existence(init, state) == True:
            #return self.data['state'][str(int)]['ts'][str(state)]['reaction_type']
        if self._final_existence(init, state) == True:
            return self.data['state'][str(init)]['final'][str(state)]['reaction_type']
