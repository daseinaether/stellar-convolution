"""
Some cosmology utility functions
"""

import astropy.cosmology as c
import astropy.units as u


##############
# Functions that handle the conversion
def age_of_universe_to_redshift(age_of_universe, cosmology):
    """
    Function to turn age of universe to redshift
    """

    # Calculate redshift
    redshift_value = c.z_at_value(
        cosmology.age, age_of_universe * u.Gyr, zmin=0, zmax=10000
    )

    return redshift_value


def lookback_time_to_redshift(lookback_time, cosmology):
    """
    Function to calculate the redshift corresponding to a lookback time

    Parameters
    ----------
    lookback_time : float
        Lookback time to be converted to redshift. Assumed to be in Gyr.
    cosmology : astropy_cosmology
        astropy cosmology object used to convert redshift and time.

    Returns
    -------
    TYPE
        Description
    """

    if lookback_time == 0:
        lookback_time = 1e-6

    # Calculate redshift
    redshift_value = c.z_at_value(
        cosmology.age, cosmology.age(0) - lookback_time * u.Gyr, zmin=0, zmax=10000
    )

    return redshift_value


def redshift_to_lookback_time(redshift, cosmology):
    """
    Function to calculate the lookback time corresponding to a certain redshift
    """

    return cosmology.age(0) - cosmology.age(redshift)


def redshift_to_age_of_universe(redshift, cosmology):
    """
    Function to calculate the age of the universe corresponding to a certain redshift
    """

    return cosmology.age(redshift)
