"""Tools to build a ode equatin function."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np


def hard_coded_ode(c: np.ndarray, t: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, ...]:  # noQA: ARG001
    """An example for an ordinary differential equation.

    This a hardcoded example for a differential equation. Necessary inputs are the concenatraion
    the rates and time. For the equation a simple first order decay is given.
    This functions needs to be given to the odeint solver.

    Args:
        c (np.ndarray): Starting concentration
        t (np.ndarray): Time
        k (np.ndarray): Rates

    Returns:
        tuple[np.ndarray,...]: ODE equations
    """
    c0 = -k[0] * c[0]
    c1 = k[0] * c[0]
    return c0, c1


def construct_ode(concentration: np.ndarray, t: np.ndarray, table: dict, rates: np.ndarray) -> tuple[float, ...]:  # noQA: ARG001
    """A tool to automatically construct ODEs.

    Args:
        concentration (np.ndarray): starting concentration
        t (np.ndarray): time
        table (dict): table with information of possible reaction
        rates (np.ndarray): reaction rates as N x N adjacency matrix. N: number of structures.

    Returns:
        tuple: concentration profile for ODEs.
    """
    odes = [0.0 for _ in concentration]

    for init, final in table.items(): # init = key, final = value (dic) #originally forward and reverse independently but can be together
        for target_final in final:
            odes[init] -= rates[init, target_final] * concentration[init]
            odes[target_final] += rates[init, target_final] * concentration[init]
    return tuple(odes)

def general_ode(concentration: np.ndarray, t: np.ndarray, table: dict, rates: np.ndarray) -> tuple[float, ...]:
    """A tool to automatically construct ODEs with considering reverse reactions.

    Args:
        concentration (np.ndarray): starting concentration
        t (np.ndarray): time
        table (dict): table with information of possible reaction
        rates (np.ndarray): reaction rates as N x N adjacency matrix. N: number of structures.

    Returns:
        tuple: concentration profile for ODEs.
    """
    odes = [0.0 for _ in concentration]
    
    for init, final in table.items():
        odes[init] -= (rates[init, final] * concentration[init] - rates[final, init] * concentration[final])
        odes[final] += (rates[init, final] * concentration[init] - rates[final, init] * concentration[final])

    return tuple(odes)