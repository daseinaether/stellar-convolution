"""
Module for calculating origin redshift arrays.

This module provides functions for calculating the redshift at the origin (birth) of stellar systems or events,
based on their current redshift and delay times.
"""

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve.cosmology_utils import redshift_to_lookback_time


####
# Array functions
def calculate_origin_redshift_array(
    config,
    convolution_redshift_value,
    data_dict,
):
    """
    Calculate the birth redshift array for merging/formation events.

    This function computes the origin (birth) redshift of stellar systems or events by:
    - Calculating the lookback time (age of the universe at the current redshift).
    - Adding the delay time (time between formation and the event) to the lookback time.
    - Converting the resulting birth time back to redshift using cosmological parameters.

    Parameters:

    - config (dict): Configuration dictionary containing:
        - "logger": Logger object for debug information.
        - "cosmology": Astropy cosmology object for calculations.
        - "redshift_interpolator_max_redshift": Maximum redshift for the interpolator.
        - "interpolators": Dictionary with "lookback_time_to_redshift_interpolator" for conversions.
    - convolution_redshift_value (float): The current redshift value for convolution.
    - data_dict (dict): Dictionary containing event data, including:
        - "delay_time": Astropy Quantity array of delay times (time between formation and event).

    Returns:
    - origin_redshift_values (np.ndarray): Array of calculated birth redshifts for the events.
      Values outside the valid range are set to -1.

    Example:
    ```python
    config = {
        "logger": logger,
        "cosmology": cosmology,
        "redshift_interpolator_max_redshift": 10.0,
        "interpolators": {
            "lookback_time_to_redshift_interpolator": interpolator
        },
    }
    convolution_redshift_value = 0.5
    data_dict = {"delay_time": np.array([1.0, 2.0]) * u.Gyr}
    birth_redshifts = calculate_origin_redshift_array(config, convolution_redshift_value, data_dict)
    ```

    """

    config["logger"].debug(
        "Calculating origin redshift of systems by converting to lookback time, adding delay time and converting back to redshift."
    )

    # With the current redshift, we calculate the lookback time, subtract the merger time and formation time and turn back into
    current_lookback_value = redshift_to_lookback_time(
        convolution_redshift_value, cosmology=config["cosmology"]
    )

    # Calculate lookback time of first starformation (which is the same as the upper redshift of the interpolator)
    lookback_time_of_first_starformation = redshift_to_lookback_time(
        config["redshift_interpolator_max_redshift"],
        cosmology=config["cosmology"],
    )

    # Calculate the lookback time of the event, given the delay time of the event (i.e. duration between birth and event-type) and the current time.
    origin_lookback_time_values_in_gyr = (
        current_lookback_value.to(u.yr) + data_dict["delay_time"].to(u.yr)
    ).to(u.Gyr)

    # Get the indices where the event falls inside the correct starformation time range
    indices_within_first_starformation = (
        origin_lookback_time_values_in_gyr < lookback_time_of_first_starformation
    )
    indices_outside_first_starformation = (
        origin_lookback_time_values_in_gyr >= lookback_time_of_first_starformation
    )

    # Create redshift values array of the event
    origin_redshift_values = np.ones(data_dict["delay_time"].shape)
    origin_redshift_values[indices_within_first_starformation] = config[
        "interpolators"
    ]["lookback_time_to_redshift_interpolator"](
        origin_lookback_time_values_in_gyr[indices_within_first_starformation]
    )
    origin_redshift_values[indices_outside_first_starformation] = -1

    #
    return origin_redshift_values
