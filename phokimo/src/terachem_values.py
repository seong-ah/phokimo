from __future__ import annotations

import os

import numpy as np

from phokimo.src.io.terachem import TeraChemOutputReader
from scipy import constants
from tcgm_lib.convert.converter import energy_unit


class State_Values:
    def __init__(self, toml) -> None:
        """Extract data from terachem output based on the kinetic model.

        Reads given toml file for the kinetic modeling and extract energy values from terachem outputs.

        Args:
            toml: TomlReader(toml_file_path) that reads toml file
        """
        self.toml = toml
        self.hartree_array = self.state_relative_list_hartree()
        self.jmol_array = self.state_relative_list_energy()

    def terachem_output(self, state: str, max_roots: int = 5, substate_name: str = None, ground: bool = False) -> float:
        """Get the energy value of target state.

        Args:
            state (str): name of the searching state
            max_roots (int, optional): Number of roots to extract from. Defaults to 5.
            substate_name (str): default setting is None, name if targeting substate exists
            ground (bool): default setting is excited, True when theory level is ground state calculation

        Returns:
            float: energy (Eh)
        """
        file_paths = self.toml.file_path
        target_spin_state = self.toml.target_spin_state(state)
        if self.toml.reference_state() == state:
            file_path = self.toml.ref_file_path(ground)
            terachem = TeraChemOutputReader(file_path)
            if not ground:
                energy_tuple = terachem.ci_energy(max_roots)
                energy_value = (energy_tuple[0][target_spin_state[0]] + energy_tuple[0][target_spin_state[1]]) / 2
            else:
                with open(file_path, 'r') as file:
                    for line in file:
                        if "FINAL ENERGY:" in line:
                            energy_value = float(line.split()[2])
                            break

        elif substate_name == None:
            file_path = file_paths[state]
            terachem = TeraChemOutputReader(file_path)
            if self.toml.theory_level(state) == "excited":
                energy_tuple = terachem.ci_energy(max_roots)
                energy_value = (energy_tuple[0][target_spin_state[0]] + energy_tuple[0][target_spin_state[1]]) / 2
            else:
                with open(file_path, 'r') as file:
                    for line in file:
                        if "FINAL ENERGY:" in line:
                            energy_value = float(line.split()[2])
                            break
        else:
            file_path = file_paths[state][substate_name]
            terachem = TeraChemOutputReader(file_path)
            if self.toml.theory_level(substate_name, substate = True) == "excited":
                energy_tuple = terachem.ci_energy(max_roots)
                energy_value = (energy_tuple[0][target_spin_state[0]] + energy_tuple[0][target_spin_state[1]]) / 2
            else:
                with open(file_path, 'r') as file:
                    for line in file:
                        if "FINAL ENERGY:" in line:
                            energy_value = float(line.split()[2])
                            break
        return energy_value

    def state_relative_list_hartree(self) -> np.ndarray:
        """Generate an array of relative energies.

        Returns:
            np.ndarray: array of energies relative to reference state (Eh)
        """
        num_states = self.toml.len_states
        reference_state = self.toml.ref_name
        state_relative_list_energy = np.zeros(num_states)
        for state in self.toml.name_to_num:
            state_num = self.toml.state_num(state)
            hartree_energy = 0.0
            if state == reference_state:
                state_relative_list_energy[state_num] = 0.0
            else:
                if self.toml.substate_existence(state):
                    for substate in self.toml.substate_list(state):
                        hartree_energy += self.terachem_output(state, substate_name = substate)
                else:
                    hartree_energy = self.terachem_output(state)
                if self.toml.theory_level(state) == "excited":
                    reference_energy = self.terachem_output(reference_state)
                else:
                    reference_energy = self.terachem_output(reference_state, ground = True)
                state_relative_list_energy[state_num] = hartree_energy - reference_energy
        return state_relative_list_energy

    def state_relative_list_energy(self) -> np.ndarray:
        """Converts unit to J/mol from relative energies in an array.

        Returns:
            np.ndarray: array of energies relative to reference state (J/mol)
        """
        state_list_energy = energy_unit(self.hartree_array, "eh", "j/mol")
        return state_list_energy

class Reactions:
    def __init__(self, toml, rate_constant, state_list_energy: np.ndarray):
        """Figure out reaction connection and calculate rate constant.

        Use toml data for reaction and terachem output to calculate rate constants.

        Args:
            toml: TomlReader(toml_file_path) that reads toml file
            rate_constant: RateCalculator() that calculates rate constants
            state_list_energy (np.ndarray): relative energy(J/mol) of each state
        """
        self.toml = toml
        self.rate_constant = rate_constant

        self.energy_differences = self.dEs(state_list_energy)
    
    def dEs(self, state_list_energy: np.ndarray) -> np.ndarray:
        """Calculate the energy differences of each reaction.

        dEs[initial_state][final_state] = final energy - initial energy
        if reaction has transition state, dEs[initial_state][final_state] = excitation energy = transition state energy - initial energy

        Args:
            state_list_energy (np.ndarray): relative energy(J/mol) of each state

        Returns:
            np.ndarray: rate constant values in matrix form (k = dEs[initial state][final state])
        """
        dim = (self.toml.len_states, self.toml.len_states)
        dEs = np.zeros(dim)
        
        for init in self.toml.name_to_num:
            init_num = self.toml.state_num(init)
            for next in self.toml.name_to_num:
                if self.toml.ts_existence(init, next):
                    ts_num = self.toml.ts_num(init, next)
                    ts_final_num = self.toml.ts_final_num(init, next)
                    dE = state_list_energy[ts_num] - state_list_energy[init_num]
                    dEs[init_num][ts_final_num] = dE

                if self.toml.final_existence(init, next):
                    final_num = self.toml.final_num(init, next)
                    if self.toml.reaction_type(init, next) == "vibrational relaxation":
                        dE = state_list_energy[final_num] - state_list_energy[init_num]
                        dEs[init_num][final_num] = dE
        return dEs
    
    def T_eq(self, state_list_energy: np.ndarray, reactant: str, temp_state: str) -> float:
        """Calculate the equilibriated temperature of the reaction.

        .. math::
            T_{\text{eq}} = T_{\text{neq}} + \frac{n \Delta \epsilon}{(3N - 6) R}
            T_{\text{neq}} = 300K


        Args:
            state_list_energy (np.ndarray): relative energy(J/mol) of each state
            reactant (str): corresponding reactant of the reaction
            temp_state (str) : corresponding subsequent state from the reactant that defines :sub:`T_{eq}`

        Returns:
            np.ndarray: rate constant values in matrix form (k = dEs[initial state][final state])
        """
        R_GAS = constants.R

        total_atoms = self.toml.total_atoms()
        reactant_num = self.toml.state_num(reactant)
        temp_state_num = self.toml.state_num(temp_state)

        return 300 + (state_list_energy[reactant_num] - state_list_energy[temp_state_num]) / ((3 * total_atoms - 6) * R_GAS)

    def rates(self, state_list_energy: np.ndarray) -> np.ndarray:
        """Calculate the rate constants of each reaction.

        Args:
            state_list_energy (np.ndarray): energy(J/mol) of each state

        Returns:
            np.ndarray: rate constants
        """
        
        graph_teq_group = self.toml.teq_graph
        dim = (self.toml.len_states, self.toml.len_states)
        rates = np.zeros(dim)
        
        reactants = self.toml.reactant_list_name()

        dEs = self.dEs(state_list_energy)

        total_atoms = float(self.toml.total_atoms())
        
        for reactant in reactants:
            print(graph_teq_group)
            graph_teq_reactant = graph_teq_group[reactant]
            print(graph_teq_reactant)
            for init in self.toml.name_to_num:
                init_num = self.toml.state_num(init)
                for next in self.toml.name_to_num:
                    next_num = self.toml.state_num(next)
                    if init == reactant:
                        if self.toml.reaction_type(init, next) == "vibrational relaxation":
                            normal_mode = self.toml.normal_mode(init, next)
                            rate_constant = self.rate_constant.relaxation_theory.compute_rate(dEs[init_num][next_num], normal_mode, total_atoms, T = 300)
                            rates[init_num][next_num] = rate_constant
                    else:
                        for temp_state in graph_teq_reactant:
                            T_eq = self.T_eq(state_list_energy, reactant, temp_state)
                            if (init, next) in graph_teq_reactant.get(temp_state, []) or self.toml.ts_existence(init, next):
                                if self.toml.ts_existence(init, next):
                                    print("transition check", init, next)
                                    fin = self.toml.ts_final_name(init, next)
                                    fin_num = self.toml.state_num(fin)
                                    rate_constant = self.rate_constant.reaction_theory.compute_rate(dEs[init_num][fin_num], T = 329.05)
                                    rates[init_num][fin_num] = rate_constant
                                elif self.toml.reaction_type(init, next) == "vibrational relaxation":
                                    normal_mode = self.toml.normal_mode(init, next)
                                    print("relax", init, next, T_eq, normal_mode)
                                    rate_constant = self.rate_constant.relaxation_theory.compute_rate(dEs[init_num][next_num], normal_mode, total_atoms, T = T_eq)
                                    rates[init_num][next_num] = rate_constant
        return rates