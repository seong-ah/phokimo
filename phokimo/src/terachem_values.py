from __future__ import annotations

import os

import numpy as np

from phokimo.src.io.terachem import TeraChemOutputReader
from tcgm_lib.convert.converter import energy_unit


class State_Values:
    def __init__(self, toml) -> None:
        """Extract data from terachem output based on the kinetic model.

        Reads given toml file for the kinetic modeling and get data from terachem based on it.

        Args:
            toml: TomlReader(toml_file_path) that reads toml file
        """
        self.toml = toml

    def terachem_output(self, num: int, calculation_path: str, max_roots: int = 5, substate: bool = False, hhtda: bool = True) -> tuple:
        """Get the energy and oscillation strength from terachem output.

        Args:
            num (int): numbering of the searching state
            calculation_path (str): absolute path of calculation directory
            max_roots (int, optional): Number of roots to extract from. Defaults to 5.

        Returns:
            tuple: energies (Eh) and oscillator strengths (-)
        """
        file_path = self.toml.file_path(num, substate, hhtda)
        assert os.path.exists(file_path)
        terachem = TeraChemOutputReader(file_path)
        if self.toml.mult(num) == 2 or hhtda == False:
            with open(file_path, 'r') as file:
                for line in file:
                    if "FINAL ENERGY:" in line:
                        energy_value = float(line.split()[2])
                        break
        else:
            energy_value = terachem.ci_energy(max_roots)
        return energy_value

    def state_list_hartree(self, calculation_path: str, hhtda = True) -> np.ndarray:
        """Generate a list with a hartree energy(Eh) of each state.

        Args:
            calculation_path (str): absolute path of calculation directory

        Returns:
            list: energy(Eh) of each state
        """
        state_list_hartree = np.zeros(self.toml.num_states())
        for i in range(self.toml.num_states()):
            hartree_energy = 0.0
            if self.toml.substate_existence(i):
                for j in self.toml.substate_list(i):
                    if hhtda == True and self.toml.mult(j, substate = True) != 2:
                        hartree_energy += self.terachem_output(j, calculation_path, substate = True)[0][self.toml.target_spin_state(j, substate = True)[0]]
                    else:
                        hartree_energy += self.terachem_output(j, calculation_path, substate = True, hhtda = False)
            else:
                if hhtda == True and self.toml.mult(i) != 2:
                    hartree_energy = (self.terachem_output(i, calculation_path)[0][self.toml.target_spin_state(i)[0]] + self.terachem_output(i, calculation_path)[0][self.toml.target_spin_state(i)[1]]) / 2
                else:
                    hartree_energy = self.terachem_output(i, calculation_path, hhtda = False)
            state_list_hartree[i] = hartree_energy
        return state_list_hartree
    
    def state_relative_list_hartree(self, calculation_path: str) -> np.ndarray:
        num_states = self.toml.num_states()
        reference_state = self.toml.reference_state()
        state_list_hartree = self.state_list_hartree(calculation_path)
        state_relative_list_energy = [(x - state_list_hartree[reference_state]) for x in state_list_hartree]
        for i in range(num_states):
            if self.toml.mult(i) == 2:
                reference_energy = self.terachem_output(reference_state, calculation_path, hhtda = False)
                absolute_energy = state_list_hartree[i]
                state_relative_list_energy[i] = absolute_energy - reference_energy
        return np.asarray(state_relative_list_energy)

    def state_list_energy(self, calculation_path: str) -> np.ndarray:
        """Generate a list with a energy(J/mol) of each state.

        Returns:
            list: energy(J/mol) of each state
        """
        state_list_energy = energy_unit(self.state_relative_list_hartree(calculation_path), "eh", "j/mol")
        return state_list_energy

    def oscilstr(self, calculation_path: str) -> list[tuple]:
        """Generate a list with an oscillation strength of each state.

        Oscillation strength only appears for excited states, so would be zero for the ground state.

        Args:
            calculation_path (str): absolute path of calculation directory

        Returns:
            list: oscillation strength for each state
        """
        state_list_oscil = np.zeros(self.num_states)
        for i in range(self.num_states):
            if self.toml.target_spin_state(i)[0] != 0:
                oscilstr = self.terachem_output(i, calculation_path)[1][self.toml.target_spin_state(i)[0] - 1]
                state_list_oscil[i] = oscilstr
        return state_list_oscil


class Reactions:
    def __init__(self, toml, rate_constant):
        """Figure out reaction connection and calculate rate constant.

        Use toml data for reaction and terachem output to calculate rate constants.

        Args:
            toml: TomlReader(toml_file_path) that reads toml file
            rate_constant: RateCalculator() that calculates rate constants
        """
        self.toml = toml
        self.rate_constant = rate_constant

    
    def dEs(self, state_list_energy: list) -> np.ndarray:
        """Calculate the energy differences of each reaction.

        dEs[initial_state][final_state] = final energy - initial energy (for reactions with transition state, excitation energy)

        Args:
            state_list_energy (list): energy(J/mol) of each state

        Returns:
            np.ndarray: rate constants
        """
        dim = (self.toml.num_states(), self.toml.num_states())
        dEs = np.zeros(dim)
        
        for i in range(self.toml.num_states()):
            init_num = self.toml.state_num(i)
            for j in range(self.toml.num_states()):
                if self.toml.ts_existence(i, j):
                    ts_num = self.toml.ts_num(i, j)
                    ts_final_num = self.toml.ts_final_num(i, j)
                    dE = state_list_energy[ts_num] - state_list_energy[init_num]
                    dEs[init_num][ts_final_num] = dE

                if self.toml.final_existence(i, j):
                    final_num = self.toml.final_num(i, j)
                    if self.toml.reaction_type(init_num, final_num) == "relaxation":
                        dE = state_list_energy[final_num] - state_list_energy[init_num]
                        dEs[init_num][final_num] = dE
        return dEs

    def rates(self, state_list_energy: list) -> np.ndarray:
        """Calculate the rate constants of each reaction.

        Args:
            state_list_energy (list): energy(J/mol) of each state

        Returns:
            np.ndarray: rate constants
        """
        dim = (self.toml.num_states(), self.toml.num_states())
        rates = np.zeros(dim)

        reaction_types = self.toml.reaction_types()
        dEs = self.dEs(state_list_energy)

        total_atoms = float(self.toml.total_atoms())
        normal_modes = 1.0  # Should be extended later for each reaction

        for i in range(self.toml.num_states()):
            for j in range(self.toml.num_states()):
                if reaction_types[i][j] == 1:
                    rate_constant = self.rate_constant.reaction_theory.compute_rate(dEs[i][j])
                    rates[i][j] = rate_constant
                if reaction_types[i][j] == 2:
                    rate_constant = self.rate_constant.relaxation_theory.compute_rate(dEs[i][j], normal_modes, total_atoms)
                    rates[i][j] = rate_constant
        return rates
