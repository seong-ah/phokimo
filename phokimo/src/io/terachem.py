"""A tool to extract information from terachem outputs."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from matplotlib import axis


class TeraChemOutputReader:
    """Extract information from TeraChem outputs."""

    def __init__(self, fname: str) -> None:
        """Extract information from TeraChem outputs.

        Reads in the file of the TeraChem and checks if it finished.
        Otherwise AssertionError will be raised.

        Args:
            fname (_type_): filename of the output file.
        """
        assert fname.endswith(".out"), "Wrong filetype given. Must be .out"

        with open(fname) as file:
            self.lines: list[str] = file.readlines()
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

        return float(self.lines[j].split()[2])

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

        if j != 2:  # noQA: PLR2004
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
        """Extract the state dipole moments.

        Args:
            max_roots (int, optional): Number of dipole moments to be extracted. Defaults to 996.

        Returns:
            np.ndarray: State dipole moments
        """
        substring = " state dipole moments:"
        state_dipole_moments = []

        j = self._search_latest_str(substring)

        for idx in range(4, max_roots + 4):
            if self.lines[j + idx] != "\n":
                temp: list = self.lines[j + idx].split()
                temp.pop(0)
                temp.pop(-1)
                temp = [float(val) for val in temp]
                state_dipole_moments.append(temp)
            else:
                break

        return np.array(state_dipole_moments)

    def transition_dipole_moment(self, max_roots: int = 996) -> pd.DataFrame:
        """Extract transition dipole moments.

        Args:
            max_roots (int, optional): Number of roots to consider for the transition dipole moments. Defaults to 996.

            Every state pair is put into the dataframe. In each case the 'intial' and 'final' states are given.
            The dipole moment are given in cartesian coordinates (access via 'x', 'y', 'z') and the norm via 'norm'

        Returns:
            pd.DataFrame: returns a data frame with initial and final states as dipole and its norm.
        """
        substring = "Transition dipole moments "

        state_dipole_moments = []

        j = self._search_latest_str(substring)

        for idx in range(4, max_roots + 4):
            if self.lines[j + idx] != "\n":
                temp: list = self.lines[j + idx].split()
                temp.pop(1)
                temp[:2] = [int(val) for val in temp[:2]]
                temp[2:] = [float(val) for val in temp[2:]]
                state_dipole_moments.append(temp)
            else:
                break

        return pd.DataFrame(
            state_dipole_moments,
            columns=["init", "final", "x", "y", "z", "norm"],
        )

    def scf_iterations(self, max_steps: int = 101) -> pd.DataFrame:
        """Returns the last scf iteration.

        Args:
            max_steps (int, optional): Number of steps that should be used. Defaults to 101.

        Returns:
            pd.DataFrame: Table of information for the scf.
        """
        line_idx = self._search_latest_str(
            "                      *** Start SCF Iterations ***",
        )
        scf_information = []
        for idx in range(line_idx + 6, line_idx + max_steps * 2):
            if not ("Reached max number of SCF iterations" in self.lines[idx] or "-" * 20 in self.lines[idx]):
                splitted_line = self.lines[idx].split()

                if splitted_line[1].isdigit():
                    scf_information.append(map(float, splitted_line[1:]))
            else:
                break

        scf_information_df = pd.DataFrame(
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
        scf_information_df.set_index("Iter")

        return scf_information_df

    def transition_electric_dipole_moment(self, max_roots: int = 996) -> pd.DataFrame:
        """Extract electric transition dipole moments, energies and oscillator strengths.

        Args:
            max_roots (int, optional): Number of roots to consider for the transition dipole moments. Defaults to 996.

            Every state pair is put into the dataframe. In each case the 'final' state is given.
            The electric transition dipole moment are given in cartesian
            coordinates (access via 'x', 'y', 'z') and the norm via 'norm'

        Returns:
            pd.DataFrame: table of information of electric transition dipole moment.
        """
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
            },
        )

        trans_dipoles["energy"] = energies
        trans_dipoles["fosc"] = fosc

        return trans_dipoles

    def _gaussian(  # noQA: PLR0913
        self,
        x: np.ndarray,
        y: np.ndarray,
        xmin: float,
        xmax: float,
        xstep: int,
        sigma: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Applies gaussian smearing on x and y values.

        Args:
            x (np.ndarray): x values for smearing
            y (np.ndarray): y values for smearing
            xmin (float): min value to start with
            xmax (float): max value to end with
            xstep (int): stepsize
            sigma (float): standard deviation/ smearing constant

        Returns:
            tuple[np.ndarray,np.ndarray]: smeared x and y values.
        """
        xi = np.arange(xmin, xmax, xstep)
        yi = np.zeros(len(xi))
        for i in range(len(xi)):
            for k in range(len(y)):
                yi[i] = yi[i] + y[k] * np.e ** (-((xi[i] - x[k]) ** 2) / (2 * sigma**2))
        return xi, yi

    def plot_spectrum(  # noQA: PLR0913
        self,
        ax: axis.Axis,
        X: np.ndarray,
        Y: np.ndarray,
        xmin: float = -np.inf,
        xmax: float = np.inf,
        xstep: int = 1,
        gamma: float = 10,
    ) -> None:
        """Plots a spectrum of gaussian smeared variables.

        Args:
            ax (axis.Axis): Axis on which the spectrum is plotted.
            X (np.ndarray): x values
            Y (np.ndarray): y values
            xmin (float | None, optional): minimal value for x values . Defaults to None.
            xmax (float | None, optional): maximal value for x values. Defaults to None.
            xstep (int | None, optional): stepsize for x . Defaults to None.
            gamma (float, optional): smearing constant. Defaults to 10.
        """
        if xmin is -np.inf:
            xmin = np.min(X) * 0.8

        if xmax is np.inf:
            xmax = np.max(X) * 1.2

        xi, yi = self._gaussian(
            X,
            Y,
            xmin,
            xmax,
            xstep,
            gamma / np.sqrt(4 * 2 * np.log(2)),
        )

        ax.plot(xi, yi / np.sum(yi), label="Gaussian")
        ax.set_xlim([xmin, xmax])
