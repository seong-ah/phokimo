"""The main entry point for the application."""

from __future__ import annotations

from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from phokimo.src.ode_builder import construct_ode, hard_coded_ode, general_ode


def main() -> None:
    """Run the application."""
    start_conc = np.array([0.0, 0.05, 5.0, 0.1, 0.15, 0.2])
    time = np.linspace(0, 10, 1000)
    '''
    # Example 1
    #A > B only forward reaction?
    
    rates = np.array([3.0])
    func = partial(hard_coded_ode, k=rates)
    conc = odeint(func, start_conc, time) 

    plt.plot(time, conc)
    plt.show()
    '''

    # Example 2
    #table = {0: 1} #possible reactions 0 to 1 and 1 to 0

    #rates = np.array([[0.0, 3.0], [3.0, 0.0]]) #Forward & reverse reaction constants as a matrix form

    table = {2: 4, 4: 0, 2: 5, 5: 3, 3: 0, 3: 1}
    rates = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 4.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 6.0, 1.0], [3.0, 3.0, 0.0, 0.0, 0.0, 0.0], [6.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 3.0, 0.0, 0.0]])
    func = partial(construct_ode, table=table, rates=rates)
    conc = odeint(func, start_conc, time)

    plt.plot(time, conc)
    plt.show()
    

    #initial_simulation
    #table = {"ems1": "bs1", "ems1": "cis2", "bs1": "cis1", "cis1": "gs2", "cis1": "gs1", "cis2": "gs1"}
    
    func = partial(general_ode, table=table, rates=rates)
    conc = odeint(func, start_conc, time)

    plt.plot(time, conc)
    plt.show()   


if __name__ == "__main__":
    main()
