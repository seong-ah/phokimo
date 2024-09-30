"""Tools to build a ode equatin function."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np
    from scipy.integrate import odeint


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


def construct_ode(concentration: np.ndarray, t: np.ndarray, table: list[tuple], rates: np.ndarray) -> tuple[float, ...]:  # noQA: ARG001
    """A tool to automatically construct ODEs.

    Args:
        concentration (np.ndarray): starting concentration
        t (np.ndarray): time
        table (list[tuple]): table with information of elementary reactions
        rates (np.ndarray): reaction rates as N x N adjacency matrix. N: number of states.

    Returns:
        tuple: concentration profile for ODEs.
    """
    odes = [0.0 for _ in concentration]

    for reaction in table:
        init = reaction[0]
        final = reaction[1]
        odes[init] -= rates[init, final] * concentration[init]
        odes[final] += rates[init, final] * concentration[init]
    return tuple(odes)