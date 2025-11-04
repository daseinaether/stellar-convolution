"""
Functions to handle calculating the star formation and meta;licity distributions

TODO: Force these functions to return a certain astropy unit
TODO: go through all papers and find what unit they use
TODO: add https://arxiv.org/abs/2111.13704 (MW models)
TODO: add https://ui.adsabs.harvard.edu/abs/2020ApJ...898...71B/abstract (MW models)
TODO: https://arxiv.org/abs/2208.05938
TODO: https://ui.adsabs.harvard.edu/abs/2012ARA%26A..50..531K/abstract
TODO: https://arxiv.org/abs/1208.4256
TODO: https://arxiv.org/abs/1901.11321
TODO:  https://www.aanda.org/articles/aa/full_html/2014/11/aa24441-14/aa24441-14.html


DH0001_file
"""

import astropy.units as u
import numpy as np


def madau_dickinson_sfr(redshifts, a, b, c, d):
    """
    Cosmological star formation rate density from Madau & Dickinson (https://ui.adsabs.harvard.edu/abs/2014ARA%26A..52..415M/abstract)
        {'a': 0.015, 'b': 2.7, 'c': 2.9, 'd': 5.6}
    as a function of redshift.

    Used in Neijsel et al 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3740N/abstract)
    """

    rate = a * np.power(1 + redshifts, b) / (1 + np.power((1 + redshifts) / c, d))

    return rate * (u.Msun / (u.yr * (u.Mpc**3)))


def starformation_rate_distribution_vanSon2023(redshifts):
    """
    Cosmological star formation rate density used in van Son 2021. Based on Madau & Dickinson SFR.
    """

    return madau_dickinson_sfr(redshifts=redshifts, a=0.02, b=1.48, c=4.45, d=5.90)


def mor19_sfr(config, lookback_time):
    """
    Star-formation rate as a function of time (years) since the birth of the milky way, based on Gaia DR2 from Mor et al. 2019 (https://ui.adsabs.harvard.edu/abs/2019A%26A...624L...1M/abstract)

    Input: time since birth of the Galaxy / years
    """

    # print(lookback_time)

    # lookback_time_in_Gyr = lookback_time.to(u.Gyr)

    # hence age in Gyr
    age_Gyr = config["cosmology"].age(0) - lookback_time

    # data only goes back 10Gyr
    lower_age_limit = config["cosmology"].age(0) - 10 * u.Gyr

    # if its further back we just set it to 0
    age_Gyr[age_Gyr < lower_age_limit] = 0

    age_Gyr = age_Gyr.value

    # rough fit to Mor et al. (2019) : Msun/Gyr/pc^2
    sfr = (0.7946) * np.exp((0.2641) * age_Gyr) + (7.3643) * np.exp(
        -((age_Gyr - (2.5566)) ** 2) / (3.036)
    )

    # convert to Msun/year assuming (as Mor+ do)
    # 1Msun/year now == 1.58Msun/Gyr/pc^2 (data point for now)
    sfr *= 1.0 / 1.58

    return sfr * u.Msun / u.yr
