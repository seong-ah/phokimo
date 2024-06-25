"""A library to store a variety of rate constants."""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from scipy import constants

R_GAS = constants.R
H_PLANCK = constants.Planck
K_BOLTZ = constants.Boltzmann
ELEM_CHRG = constants.elementary_charge
ELEC_MASS = constants.electron_mass
VAC_PERM = constants.epsilon_0
SPEED_OF_LIGHT = constants.speed_of_light
PI = np.pi


class RelaxationTheory(ABC):
    """This is an abstract class for the declaration of common operations for the relaxation theories."""

    def __init__(self) -> None:
        """Instantiate the relaxation theory object."""
        super().__init__()

    @abstractmethod
    def compute_rate(
        self,
        dE: float,
        n_modes: float,
        n_atoms: float,
        T: float = 298.15,
    ) -> float:
        """Abstract method for a relaxtion rate.

        Args:
            dE (float): Energy difference
            n_modes (float): number of involved normal modes
            n_atoms (float): number of atoms
            T (float, optional): temperature. Defaults to 298.15.

        Returns:
            float: _description_
        """


class EmissionTheory(ABC):
    """This is an abstract class for the declaration of common operations for the emission theories."""

    def __init__(self) -> None:
        """Instantiate the emission theory object."""
        super().__init__()

    @abstractmethod
    def compute_rate(
        self,
        nu: float,
        f12: float,
        g1: float = 1,
        g2: float = 1,
    ) -> float:
        """Abstract method for the computation of emission rates.

        Args:
            nu (float): wave length (nm)
            f12 (float): oscillator strength (-)
            g1 (float, optional): degeneracy state 1. Defaults to 1.
            g2 (float, optional): degeneracy state 2. Defaults to 1.

        Returns:
            float: rate
        """


class ReactionTheory(ABC):
    """This is an abstract class for the declaration of common operations for the reaction theories."""

    def __init__(self) -> None:
        """Instantiate the reaction theory object."""
        super().__init__()

    @abstractmethod
    def compute_rate(
        self,
        dG: float,
        T: float = 298.15,
        kappa: float = 1,
    ) -> float:
        """Abstract method for the computation of reaction rates.

        Args:
            dG (float): _description_
            T (float, optional): _description_. Defaults to 298.15.
            kappa (float, optional): _description_. Defaults to 1.

        Returns:
            float: reaction rate
        """


class EyringEquation(ReactionTheory):
    """A theory for the computation of reaction rates."""

    def compute_rate(
        self,
        dG: float,
        T: float = 298.15,
        kappa: float = 1,
    ) -> float:
        """Compute the rate constant using the Eyring equation.

        Eyring: https://doi.org/10.1039/TF9353100875 \n
        Eyring and Polyani: https://doi.org/10.1063/1.1749604

        .. math::
            k = \\kappa \\cdot \\dfrac{k_\\mathrm{b}T}{h}\\exp-\\dfrac{\\Delta G}{RT}

        Args:
            dG (float): Gibbs energy of activation in J/mol
            T (float, optional): Temperature. Defaults to 298.15.
            kappa (float, optional): transimission coefficient. Defaults to 1.

        Returns:
            float: Eyring eq. derived rate constant
        """
        return kappa * K_BOLTZ * T / H_PLANCK * np.exp(-dG / (R_GAS * T))


class EinsteinCoeffientA12(EmissionTheory):
    """A theory for the computation of emission rates."""

    def compute_rate(
        self,
        nu: float,
        f12: float,
        g1: float = 1,
        g2: float = 1,
    ) -> float:
        """Rate constant for spontaneous emission derived by Einstein coeffs.

        Equation from Table I: https://doi.org/10.1119/1.12937 \n
        Original work: A. Einstein, Zur Quantentheorie der Strahlung, \n
        Physikalische Zeitschrift, 1917 ,18, 121.

        .. math::
            A_{21} = \\dfrac{g_1}{g_2} \\cdot \\dfrac{\\pi \\nu  e^2 }{\\epsilon_0 m_\\mathrm{el}}

        Args:
            nu (float): Excitation energy in nm
            f12 (float): Oscillator strength
            g1 (float, optional): Degeneracy factor. Defaults to 1.
            g2 (float, optional): Degeneracy factor. Defaults to 1.

        Returns:
            float: Spontaneous emission coefficient
        """
        prefactor = 2 * PI * nu**2 * ELEM_CHRG**2 / (VAC_PERM * ELEC_MASS)
        return prefactor * g1 / g2 * f12


class AdhocRelaxation(RelaxationTheory):
    """A theory for the computation of relaxtion rates."""

    def compute_rate(
        self,
        dE: float,
        n_modes: float,
        n_atoms: float,
        T: float = 298.15,
    ) -> float:
        """Ad hoc ansatz for the approximation of the relaxation of structure to the ground state.

        This is a self developed approach for the description of the relaxation of vibronically excited states.
        It is based on the Eyring equation and considers #
        that a relaxation occurs through a relaxation of the normal mode.
        It depends on the excess energy between the ground state and its initial state.

        .. math::
            k_\\mathrm{relax} = \\dfrac{k_\\mathrm{B}T}{h}
            \\exp\\dfrac{-n\\Delta\\epsilon_\\mathrm{init}^\\mathrm{final}}{RT}

        Args:
            dE (float): Energy difference between excited and ground structure in J/mol.
            n_modes (float): Number of normal modes involved in the relaxation.
            T (float): Temperature. Defaults to 298.15.
            n_atoms (float): Number of atoms of the structure

        Returns:
            float: Rate of relaxation.
        """
        return K_BOLTZ * T / H_PLANCK * np.exp(-(n_modes * dE) / ((3 * n_atoms - 6) * R_GAS * T))

class RateCalculator:
    """Interface to define the varying theories that should be used for the computation of the rates."""

    def __init__(
        self,
        reaction_theory: ReactionTheory = EyringEquation(),  # noQA: B008
        emission_theory: EmissionTheory = EinsteinCoeffientA12(),  # noQA: B008
        relaxation_theory: RelaxationTheory = AdhocRelaxation(),  # noQA: B008
    ) -> None:
        """Sets the default theories for the computation of the rate constants of different systems.

        Args:
            reaction_theory (ReactionTheory, optional): Theory for Reactions. Defaults to EyringEquation().
            emission_theory (EmissionTheory, optional): Theory for Emission. Defaults to EinsteinCoeffientA12().
            relaxation_theory (RelaxationTheory, optional): Theory for Relaxation. Defaults to AdhocRelaxation().
        """
        super().__init__()

        self._reaction_theory = reaction_theory
        self._emission_theory = emission_theory
        self._relaxation_theory = relaxation_theory

    @property
    def reaction_theory(self) -> ReactionTheory:
        """Property attribute of the reaction theory in the rate calculator.

        Returns:
            ReactionTheory: reaction theory
        """
        return self._reaction_theory

    @reaction_theory.setter
    def reaction_theory(self, reaction_theory: ReactionTheory) -> None:
        """Setter for the relaxation theory.

        Args:
            reaction_theory (RelaxationTheory): New relaxation theory
        """
        self._reaction_theory = reaction_theory

    @property
    def emission_theory(self) -> EmissionTheory:
        """Property attribute of the emission theory in the rate calculator.

        Returns:
            EmissionTheory: emission theory
        """
        return self._emission_theory

    @emission_theory.setter
    def emission_theory(self, emission_theory: EmissionTheory) -> None:
        """Setter for the relaxation theory.

        Args:
            emission_theory (RelaxationTheory): New relaxation theory
        """
        self._emission_theory = emission_theory

    @property
    def relaxation_theory(self) -> RelaxationTheory:
        """Property attribute of the relaxation theory in the rate calculator.

        Returns:
            RelaxationTheory: relaxation theory
        """
        return self._relaxation_theory

    @relaxation_theory.setter
    def relaxation_theory(self, relaxation_theory: RelaxationTheory) -> None:
        """Setter for the relaxation theory.

        Args:
            relaxation_theory (RelaxationTheory): New relaxation theory
        """
        self._relaxation_theory = relaxation_theory
