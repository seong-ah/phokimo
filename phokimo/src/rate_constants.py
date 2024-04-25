"""A library to store a variety of rate constants."""

from __future__ import annotations

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


def eyring_equation(
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


def einstein_coefficient_A12(  # noQA: N802
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


def n_modes_relaxation(
    dE: float,
    n_modes: float,
    n_atoms: float,
    T: float = 298.15,
) -> float:
    """Ad hoc ansatz for the approximation of the relaxation of structure to the ground state.

    This is a self developed approach for the description of the relaxation of vibronically excited states.
    It is based on the Eyring equation and considers, that a relaxation occurs through a relaxation of the normal mode.
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
    return K_BOLTZ * T / H_PLANCK * np.exp(-(n_modes * dE) / (3 * n_atoms - 6) * R_GAS * T)
