from __future__ import annotations

import os

import numpy as np

from phokimo.src.io.terachem import TeraChemOutputReader


class State_Values:
    def __init__(self, toml) -> None:
        """Extract data from terachem output based on the kinetic model.

        Reads given toml file for the kinetic modeling and get data from terachem based on it.

        Args:
            toml: TomlReader(toml_file_path) that reads toml file
        """
        self.toml = toml

    def terachem_output(self, num: int, calculation_path: str, max_roots: int = 3) -> tuple:
        """Get the energy and oscillation strength from terachem output.

        Args:
            num (int): numbering of the searching state
            calculation_path (str): absolute path of calculation directory
            max_roots (int, optional): Number of roots to extract from. Defaults to 3.

        Returns:
            tuple: energies (Eh) and oscillator strengths (-)
        """
        file_path = self.toml.file_path(num, calculation_path)
        assert os.path.exists(file_path)
        terachem = TeraChemOutputReader(file_path)
        return terachem.ci_energy(max_roots)

    def state_list_hartree(self, calculation_path: str) -> list:
        """Generate a list with a hartree energy(Eh) of each state.

        Args:
            calculation_path (str): absolute path of calculation directory

        Returns:
            list: energy(Eh) of each state
        """
        state_list_hartree = np.zeros(self.toml.num_states())
        for i in range(self.toml.num_states()):
            hartree_energy = self.terachem_output(i, calculation_path)[0][self.toml.target_spin_state(i)]
            state_list_hartree[i] = hartree_energy
        return state_list_hartree

    def state_list_energy(self, calculation_path: str) -> list:
        """Generate a list with a energy(J/mol) of each state.

        Returns:
            list: energy(J/mol) of each state
        """
        state_list_energy = [x * 2625.5 * (10**3) for x in self.state_list_hartree(calculation_path)]
        return state_list_energy

    def oscilstr(self, calculation_path: str) -> list:
        """Generate a list with an oscillation strength of each state.

        Oscillation strength only appears for excited states, so would be zero for the ground state.

        Args:
            calculation_path (str): absolute path of calculation directory

        Returns:
            list: oscillation strength for each state
        """
        state_list_oscil = np.zeros(self.num_states)
        for i in range(self.num_states):
            if self.toml.target_spin_state != 0:
                oscilstr = self.terachem_output(i, calculation_path)[1][self.toml.target_spin_state - 1]
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

    def graph_table_name(self) -> list:
        """Generate a list with reaction connection with state names.

        Returns:
            list: graph edge tuples that represent the reaction with state names
        """
        graph_table_name = []
        for i in range(self.toml.num_states()):
            init_name = self.toml.state_name(i)
            for j in range(self.toml.num_states()):
                if self.toml.ts_existence(i, j):
                    ts_name = self.toml.ts_name(i, j)
                    ts_final_name = self.toml.ts_final_name(i, j)
                    graph_table_name.append(self.toml.graph_edge(init_name, ts_name))
                    graph_table_name.append(self.toml.graph_edge(ts_name, ts_final_name))
                if self.toml.final_existence(i, j):
                    final_name = self.toml.final_name(i, j)
                    graph_table_name.append(self.toml.graph_edge(init_name, final_name))
        return graph_table_name

    def graph_table_num(self) -> list:
        """Generate a list with reaction connection with state numberings.

        Returns:
            list: graph edge tuples that represent the reaction with state numberings
        """
        graph_table_num = []
        for i in range(self.toml.num_states()):
            init_num = self.toml.state_num(i)
            for j in range(self.toml.num_states()):
                if self.toml.ts_existence(i, j):
                    ts_num = self.toml.ts_num(i, j)
                    ts_final_num = self.toml.ts_final_num(i, j)
                    graph_table_num.append(self.toml.graph_edge(init_num, ts_num))
                    graph_table_num.append(self.toml.graph_edge(ts_num, ts_final_num))
                if self.toml.final_existence(i, j):
                    final_num = self.toml.final_num(i, j)
                    graph_table_num.append(self.toml.graph_edge(init_num, final_num))
        return graph_table_num

    def rates(self, state_list_energy: list) -> np.ndarray:
        """Calculate the rate constants of each reaction.

        Args:
            state_list_energy (list): energy(J/mol) of each state

        Returns:
            np.ndarray: rate constants
        """
        dim = (self.toml.num_states(), self.toml.num_states())
        rates = np.zeros(dim)

        total_atoms = float(self.toml.total_atoms())
        normal_modes = 1.0  # Should be extended later for each reaction

        for i in range(self.toml.num_states()):
            init_num = self.toml.state_num(i)
            for j in range(self.toml.num_states()):
                if self.toml.ts_existence(i, j):
                    ts_num = self.toml.ts_num(i, j)
                    ts_final_num = self.toml.ts_final_num(i, j)
                    dE = state_list_energy[ts_num] - state_list_energy[init_num]
                    rate_constant = self.rate_constant.reaction_theory.compute_rate(dE)
                    rates[init_num][ts_final_num] = rate_constant
                if self.toml.final_existence(i, j):
                    final_num = self.toml.final_num(i, j)
                    if self.toml.reaction_type(init_num, final_num) == "relaxation":
                        dE = state_list_energy[final_num] - state_list_energy[init_num]
                        rate_constant = self.rate_constant.relaxation_theory.compute_rate(dE, normal_modes, total_atoms)
                        rates[init_num][final_num] = rate_constant
        return rates
