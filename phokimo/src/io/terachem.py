import numpy as np
import pandas as pd
import re as re


class TeraChemOutputReader:
    def __init__(self, fname: str) -> None:
        """Extract information from TeraChem outputs.

        Reads in the file of the TeraChem and checks if it finished.
        Otherwise AssertionError will be raised.

        Args:
            fname (_type_): filename of the output file.
        """
        assert fname.endswith(".out"), "Wrong filetype given. Must be .out"

        with open(fname) as file:
            self.lines: str = file.readlines()
        file.close()

        assert self._check_finish(), f"{fname} did not finish!"

    def _check_finish(self) -> bool:
        """Checks if a calculation is finished.

        Searches in the last 200 lines of the file if 'Job finished:' is part of a line.
        Otherwise it is considered that the computation did not finish.

        Returns:
            bool: True if finished otherwise False
        """

        finished = False
        substring = " Job finished:"

        max_lines = min(200, len(self.lines))

        for i in range(1, max_lines):

            if substring in self.lines[-i]:
                finished = True
                break

        return finished

    def check_convergence(self) -> bool:
        """Checks if a calculation is converged.

        Searches in the last 200 lines of the file if 'Converged!' is part of a line.
        Otherwise it is considered that the computation did not finish.

        Returns:
            bool: True if converged otherwise False
        """

        substring = "Converged!"
        converged = False

        max_lines = min(200, len(self.lines))

        for i in range(1, max_lines):

            if substring in self.lines[-i]:
                converged = True
                break

        return converged

    def _search_latest_str(self, substring: str) -> int:
        """Tool to search a substring in the lines of the output.

        This tool starts at the end of the file and stops at the first entry where the substring is in.

        Args:
            substring (str): _description_

        Returns:
            int: _description_
        """
        j = 0
        for i in reversed(range(len(self.lines))):
            if substring in self.lines[i]:
                j = i
                break
        return j

    def energy(self) -> float:
        """Extract the final energy from the output.

        Searches for the keyword 'FINAL ENERGY' in the file and extracts from this line the energy.

        Returns:
            float: energy
        """
        substring = "FINAL ENERGY"

        j = self._search_latest_str(substring)

        energy = float(self.lines[j].split()[2])

        return energy

    def ci_energy(self, max_roots: int = 3) -> tuple:
        """Reads configurational interaction energies and oscillator strength.

        !CAUTION! There is no distinguishment between singlet and triplet, or other, states. 

        Args:
            max_roots (int, optional): Number of roots to extract from. Defaults to 3.

        Returns:
            tuple: energies (Eh) and oscillator strengths (-)
        """

        energies = []
        foscs = []
        substring = "Root   Mult.   Total Energy (a.u.) "

        j = self._search_latest_str(substring)
        j += 2

        if j != 2:
            for i in range(max_roots):

                if str(i + 1) in self.lines[j + i]:
                    energies.append(float(self.lines[j + i].split()[2]))
                    if i != 0:
                        foscs.append(float(self.lines[j + i].split()[-1]))
                else:
                    break
        else:
            energies, foscs = [np.nan], [np.nan]

        return np.array(energies), np.array(foscs)

    def state_dipole_moment(self, max_roots: int = 996) -> np.ndarray:
        """Reads state diple moments."""

        substring = f" state dipole moments:"
        state_dipole_moments = []

        j = self._search_latest_str(substring)

        for idx in range(4, max_roots + 4):
            if self.lines[j + idx] == "\n":
                break
            else:
                temp = self.lines[j + idx].split()
                temp.pop(0)
                temp.pop(-1)
                temp = [float(val) for val in temp]
                state_dipole_moments.append(temp)

        return np.array(state_dipole_moments)

    def transition_dipole_moment(self, max_roots: int = 996) -> pd.DataFrame:
        """Reads transition diple moments."""

        substring = f"Transition dipole moments "

        state_dipole_moments = []

        j = self._search_latest_str(substring)

        for idx in range(4, max_roots + 4):
            if self.lines[j + idx] == "\n":
                break
            else:
                temp = self.lines[j + idx].split()
                temp.pop(1)
                temp[:2] = [int(val) for val in temp[:2]]
                temp[2:] = [float(val) for val in temp[2:]]
                state_dipole_moments.append(temp)
        state_dipole_moments = pd.DataFrame(
            state_dipole_moments, columns=["init", "final", "x", "y", "z", "norm"]
        )

        return state_dipole_moments

    def scf_iterations(self, max_steps=101) -> pd.DataFrame:
        """Reads the latest scf iteration with a set number of steps."""

        line_idx = self._search_latest_str(
            "                      *** Start SCF Iterations ***"
        )
        scf_information = []
        for idx in range(line_idx + 6, line_idx + max_steps * 2):
            if "Reached max number of SCF iterations" in self.lines[idx]:
                break
            elif "-" * 20 in self.lines[idx]:
                break

            else:

                splitted_line = self.lines[idx].split()

                if splitted_line[1].isdigit():

                    scf_information.append(map(float, splitted_line[1:]))

        scf_information = pd.DataFrame(
            scf_information,
            columns=[
                "Iter",
                "DIIS Error",
                "Energy change",
                "Electrons",
                "XC Energy",
                "Energy",
                "E_PCM",
                "Time(s)",
            ],
        )
        scf_information.set_index("Iter")

        return scf_information

    def transition_electric_dipole_moment(self, max_roots=996):

        energies, fosc = self.ci_energy(max_roots=100)

        energies = energies[1:] - energies[0]

        trans_dipoles = self.transition_dipole_moment(max_roots)

        trans_dipoles = trans_dipoles[trans_dipoles["init"] == 1]

        trans_dipoles = trans_dipoles.rename(
            columns={
                "final": "state",
                "x": "tx",
                "y": "ty",
                "z": "tz",
                "norm": "t_norm",
            }
        )

        trans_dipoles["energy"] = energies
        trans_dipoles["fosc"] = fosc

        return trans_dipoles

    def _gaussian(self, x, y, xmin, xmax, xstep, sigma):
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(y)):
                yi[i] = yi[i] + y[k] * np.e ** (-((xi[i] - x[k]) ** 2) / (2 * sigma**2))
        return xi, yi

    def plot_spectrum(self, ax, X, Y, xmin=None, xmax=None, xstep=None, gamma=10):

        if xmin is None:
            xmin = np.min(X) * 0.8

        if xmax is None:
            xmax = np.max(X) * 1.2

        if xmin is None:
            xstep = (xmax - xmin) / 300

        xi, yi = self._gaussian(
            X, Y, xmin, xmax, xstep, gamma / np.sqrt(4 * 2 * np.log(2))
        )

        ax.plot(xi, yi / np.sum(yi), label="Gaussian")
        ax.set_xlim([xmin, xmax])
