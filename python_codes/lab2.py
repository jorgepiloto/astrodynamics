""" This is the laboratory session number 2. """

import numpy as np
from poliastro.bodies import Earth
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def task1():
    r""" Comparation between all the different maneuvers. """

    rc_ra = np.linspace(1, 30, 50)
    rb_ra = np.linspace(1, 100, 50)

    hohmann_cost = []
    bielliptic_cost = []

    for R1, R2 in zip(rc_ra, rb_ra):

        r_i = 1
        r_f = R1 * r_i
        r_b = R2 * r_i

        dv_a, dv_b, dv_c, t = bielliptic(Earth.k.value, r_i, r_b, r_f)
        deltaV = np.abs(dv_a) + np.abs(dv_b) + np.abs(dv_c)
        bielliptic_cost.append(deltaV)

        dv_a, dv_b, t = hohmann(Earth.k.value, r_i, r_f)
        deltaV = np.abs(dv_a) + np.abs(dv_b)
        hohmann_cost.append(deltaV)

    plt.plot(rc_ra, hohmann_cost, 'r')
    plt.plot(rb_ra, bielliptic_cost, 'b')
    plt.show()

def hohmann(k, r_i, r_f):
    """ This function solves for Hohmann transfer.

    Parameters
    ----------
    k: float
        Gravitational parameter
    r_i: float
        Magnitude of the initial radius
    r_f: float
        Magnitude of the final radius

    Returns
    -------
    dv_a: float
        First increment of the velocity
    dv_b: float
        Second increment of the velocity
    t: float
        Time for the maneuver
    """

    a_trans = (r_f + r_i) / 2
    dv_a = np.sqrt(2 * k / r_i - k / a_trans) - np.sqrt(2 * k / r_i - k / r_i)
    dv_b = np.sqrt(k / r_f) - np.sqrt(2 * k / r_f - k / a_trans)
    t = np.pi * np.sqrt(a_trans ** 3 / k)

    return dv_a, dv_b, t

def bielliptic(k, r_i, r_b, r_f):
    """ This function solves for the Bielliptic transfer.

    Parameters
    ----------
    k: float
        Gravitational parameter
    r_i: float
        Magnitude of the initial radius
    r_b: float
        Magnitude of the second maneuver radius
    r_f: float
        Magnitude of the final desired altitude radius

    Returns
    -------
    dv_a: float
        First increment of the velocity
    dv_b: float
        Second increment of the velocity
    dv_c: float
        Last increment of the velocity
    t: float
        Time for the maneuver
    """

    a_1 = (r_i + r_b) / 2
    a_2 = (r_b + r_f) / 2

    dv_a = np.sqrt(2 * k / r_i - k / a_1) - np.sqrt(2 * k / r_i - k / r_i)
    dv_b = np.sqrt(2 * k / r_b - k / a_2) - np.sqrt(2 * k / r_b - k / a_1)
    dv_c = np.sqrt(k / r_f) - np.sqrt(2 * k / r_f - k / a_2)

    t = np.pi * np.sqrt(a_1 ** 3 / k) + np.pi * np.sqrt(a_2 ** 3 / k)

    return dv_a, dv_b, dv_c, t




if __name__ == '__main__':
    task1()
