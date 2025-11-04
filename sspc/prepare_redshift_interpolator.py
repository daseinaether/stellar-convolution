"""
Function to prepare the redshift interpolator
"""

import os
import pickle

import numpy as np
from scipy import interpolate


##################
# Functions to create the datasets to interpolate on
def create_interpolation_datasets(config):
    """
    Function that generates two datasets:
    - First a redshift to lookback time dataset
    - Then we turn it around and use the result to set up a lookback time to redshift dataset (this is much faster, we can increase the resolution while still being off faster)
    """

    #
    config["logger"].debug("Preparing redshift interpolator")

    #
    min_redshift = config["redshift_interpolator_min_redshift"]
    max_redshift = config["redshift_interpolator_max_redshift"]
    redshift_stepsize = config["redshift_interpolator_stepsize"]

    #
    if config["redshift_interpolator_use_log"]:
        if min_redshift == 0:
            min_redshift += config["redshift_interpolator_min_redshift_if_log"]
        redshift_range = 10 ** np.arange(
            np.log10(min_redshift),
            np.log10(max_redshift + redshift_stepsize),
            redshift_stepsize,
        )
    else:
        redshift_range = np.arange(
            min_redshift,
            max_redshift + redshift_stepsize,
            redshift_stepsize,
        )

    #
    redshift_to_time_array = config["cosmology"].age(0) - config["cosmology"].age(
        redshift_range
    )

    # Generate dictionary
    interpolation_data_dict = {
        "redshift_data": redshift_range,
        "lookback_time_data": redshift_to_time_array,
        "min_redshift": min_redshift,
        "max_redshift": max_redshift,
        "redshift_stepsize": redshift_stepsize,
        "interpolate_log": config["redshift_interpolator_use_log"],
    }

    # Write to file
    with open(config["redshift_interpolator_data_output_filename"], "wb") as f:
        pickle.dump(interpolation_data_dict, f)


def load_interpolation_data(config):
    """
    Function to load the interpolation dataset and return the loaded interpolators

    if rebuild_when_settings_not_match: we build the dataset based on the settings we passed

    returns:
        - redshift_to_lookback_time interpolator
        - lookback_time_to_redshift interpolator
    """

    if not os.path.isfile(config["redshift_interpolator_data_output_filename"]):
        config["logger"].debug(
            "Creating new interpolation dict, since it was not found in the place it should be"
        )

        # Rebuild data
        create_interpolation_datasets(config)

    # load data
    with open(config["redshift_interpolator_data_output_filename"], "rb") as f:
        interpolation_data_dict = pickle.load(f)

    # Check if we need to rebuild
    rebuild_interpolation_data = False
    if config["redshift_interpolator_rebuild_when_settings_mismatch"]:
        if (
            interpolation_data_dict["min_redshift"]
            != config["redshift_interpolator_min_redshift"]
        ):
            rebuild_interpolation_data = True
        if (
            interpolation_data_dict["max_redshift"]
            != config["redshift_interpolator_max_redshift"]
        ):
            rebuild_interpolation_data = True
        if (
            interpolation_data_dict["redshift_stepsize"]
            != config["redshift_interpolator_stepsize"]
        ):
            rebuild_interpolation_data = True
        if (
            interpolation_data_dict["interpolate_log"]
            != config["redshift_interpolator_use_log"]
        ):
            rebuild_interpolation_data = True

    # Or if we override it
    if config["redshift_interpolator_force_rebuild"]:
        rebuild_interpolation_data = True

    #
    if rebuild_interpolation_data:
        config["logger"].debug(
            "Creating new interpolation dict, since it was not found in the place it should be"
        )

        # Rebuild data
        create_interpolation_datasets(config)

        # Reload dict
        with open(config["redshift_interpolator_data_output_filename"], "rb") as f:
            interpolation_data_dict = pickle.load(f)

    ###########
    # Create the interpolators

    # redshift to lookback time interpolator
    redshift_to_lookback_time_interpolator = interpolate.interp1d(
        interpolation_data_dict["redshift_data"],
        interpolation_data_dict["lookback_time_data"],
        bounds_error=False,
        fill_value=0,
    )

    #
    lookback_time_to_redshift_interpolator = interpolate.interp1d(
        interpolation_data_dict["lookback_time_data"],
        interpolation_data_dict["redshift_data"],
        bounds_error=False,
        fill_value=0,
    )

    return {
        "redshift_to_lookback_time_interpolator": redshift_to_lookback_time_interpolator,
        "lookback_time_to_redshift_interpolator": lookback_time_to_redshift_interpolator,
    }


def prepare_redshift_interpolator(config):
    """
    Function to set up the redshift interpolator
    """

    # Load the dict for dict with the interpolators
    if config["time_type"] == "redshift":
        config["interpolators"] = load_interpolation_data(config)

    return config
