"""Implement the Kermack-McKendrick models of disease spread

Kermack, W. O. and McKendrick, A. G. "A Contribution to the
Mathematical Theory of Epidemics." Proc. Roy. Soc. Lond. A 115,
700-721, 1927.

Variable Definitions:
s - Susceptible Population
i - Infected Population
r - Recovered Population
n - Total Population
beta - Infection rate (Viral coefficient per Population)
gamma - Recovery rate
mu - Death rate
alpha - Birth rate
f - Immunity loss rate
dt - Numerical integration timestep in days
"""

import numpy as np


def calc_ds(s, i, r, beta, alpha=0, mu=0, f=0):
    """ Calculate the differential number of susceptible people for a single timestep.

    Args:
        s (float) - number of susceptible people
        i (float) - number of infected people
        r (float) - number of recovered people
        beta (float) - Infection rate
        alpha (float) - Birth rate
        mu (float) - Death rate
        f (float) - Immunity loss rate

    Returns:
        float
    """
    return -1 * beta * s * i + alpha * (s + i + r) - mu * s + f * r


def calc_di(s, i, beta, gamma=0, mu=0):
    """ Calculate the differential number of infected people for a single timestep.

    Args:
        s (float) - number of susceptible people
        i (float) - number of infected people
        beta (float) - Infection rate
        gamma (float) - Recovery rate
        mu (float) - Death rate

    Returns:
        float
    """
    return beta * s * i - gamma * i - mu * i


def calc_dr(i, r, gamma=0, mu=0, f=0):
    """ Calculate the differential number of recovered people for a single timestep.

    Args:
        i (float) - number of infected people
        r (float) - number of recovered people
        gamma (float) - Recovery rate
        mu (float) - Death rate
        f (float) - Immunity loss rate

    Returns:
        float
    """
    return gamma * i - mu * r - f * r


def eulerstep(s, i, r, beta, gamma=0, alpha=0, mu=0, f=0, dt=0.001):
    """ Calculate the susceptible, infected, and recovered numbersfor a single timestep.

    Use the forward Euler method to keep it simple.  This can go unstable
    for certain parameters at long times. Try making the step size smaller
    in that case.

    Args:
        s (float) - number of susceptible people
        i (float) - number of infected people
        r (float) - number of recovered people
        beta (float) - Infection rate
        gamma (float) - Recovery rate, default 0
        alpha (float) - Birth rate, default 0
        mu (float) - Death rate, default 0
        f (float) - Immunity loss rate, default 0
        dt (float) - Timestep in days, default 0.001

    Returns:
        np.array: length 3 including s, i, r
    """
    s1 = s + dt * calc_ds(s, i, r, beta, alpha, mu, f)
    i1 = i + dt * calc_di(s, i, beta, gamma, mu)
    r1 = r + dt * calc_dr(i, r, gamma, mu, f)

    # Round down to 0 to get the dynamics right in the cases of epidemics dying out.
    if i1 < 0.5:
        i1 = 0

    return np.ceil(np.array([s1, i1, r1]))


def run(s0, i0, beta, r0=0, gamma=0, alpha=0, mu=0, f=0, dt=0.001, t=100):
    """ Run the model for t days.

    Args:
        s (float) - number of susceptible people
        i (float) - number of infected people
        beta (float) - Infection rate
        r (float) - number of recovered people, default 0
        gamma (float) - Recovery rate, default 0
        alpha (float) - Birth rate, default 0
        mu (float) - Death rate, default 0
        f (float) - Immunity loss rate, default 0
        dt (float) - Timestep in days, default 0.001

    Returns:
        np.array: length t / dt, s at each timestep
        np.array: length t / dt, i at each timestep
        np.array: length t / dt, r at each timestep
    """
    def yield_steps():
        for _ in np.arange(0, T, dt):
            yield eulerstep(s0, i0, r0, beta, gamma, alpha, mu, f, dt)


    result = np.array(list(yield_steps())).T
    return result[0], result[1], result[2]
