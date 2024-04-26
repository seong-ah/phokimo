"""The main entry point for the application."""

from __future__ import annotations

from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from phokimo.src.ode_builder import construct_ode, hard_coded_ode
from phokimo.src.rate_constants import RateCalculator


def main() -> None:
    """Run the application."""
    start_conc = np.array([1.0, 0.0])
    time = np.linspace(0, 1e1, 1000)

    # Example 1
    rates = np.array([2.0])
    func = partial(hard_coded_ode, k=rates)
    conc = odeint(func, start_conc, time)

    plt.plot(time, conc)
    plt.show()

    # Example 2
    table = {0: 1, 1: 0}

    rates = np.array([[0.0, 1.0], [1.0, 0.0]])

    # An example on how to use the Rate setter
    rate_comp = RateCalculator()

    rates[0, 1] = rate_comp.reaction_theory.compute_rate(dG=70e3)
    rates[1, 0] = rate_comp.reaction_theory.compute_rate(dG=70e3)

    func = partial(construct_ode, table=table, rates=rates)
    conc = odeint(func, start_conc, time)

    plt.plot(time, conc)
    plt.show()


if __name__ == "__main__":
    main()
