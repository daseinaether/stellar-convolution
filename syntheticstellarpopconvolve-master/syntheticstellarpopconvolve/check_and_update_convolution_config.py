"""
Function to handle checking the convolution configuration

TODO: handle logic of SFR
TODO: in the SFR dict use ONLY the (dP/dz)*dZ array. let the rest be constructed from there
TODO: make sure that redshift is supported in sampling (we need to use the time bin sizes instead of lookback time bin sizes)
"""

import numpy as np
import voluptuous as vol

from syntheticstellarpopconvolve.check_and_update_convolution_instruction import (
    check_and_update_convolution_instructions,
)
from syntheticstellarpopconvolve.check_and_update_sfr_dict import (
    check_and_update_sfr_dicts,
)
from syntheticstellarpopconvolve.default_convolution_config import (
    default_convolution_config_dict,
)
from syntheticstellarpopconvolve.general_functions import is_time_unit


def update_convolution_config(config):
    """
    Function to calculate some extra quantities based on input
    """

    #
    config["logger"].debug("Updating configuration")

    ########
    # Calculate some extra convolution-bin info if needed (the other two dont need any convolution-bin info)
    requires_convolution_bin_info = any(
        [
            # convolution_instruction["convolution_type"] == "integrate" # NOTE: before, 'integrate' implied 'backward' convolution but was the only one
            convolution_instruction["convolution_direction"] == "backward"
            for convolution_instruction in config["convolution_instructions"]
        ]
    )

    # update lookback time info
    if requires_convolution_bin_info and config["time_type"] == "lookback_time":
        config["convolution_lookback_time_bin_centers"] = (
            config["convolution_lookback_time_bin_edges"][1:]
            + config["convolution_lookback_time_bin_edges"][:-1]
        ) / 2
        config["convolution_lookback_time_bin_sizes"] = np.diff(
            config["convolution_lookback_time_bin_edges"]
        )

        config["convolution_time_bin_edges"] = config[
            "convolution_lookback_time_bin_edges"
        ]
        config["convolution_time_bin_centers"] = config[
            "convolution_lookback_time_bin_centers"
        ]
        config["convolution_time_bin_sizes"] = config[
            "convolution_lookback_time_bin_sizes"
        ]

        config["logger"].debug(
            "Updated bin data: convolution_lookback_time_bin_centers: {} convolution_lookback_time_bin_sizes: {}".format(
                config["convolution_lookback_time_bin_centers"],
                config["convolution_lookback_time_bin_sizes"],
            )
        )

    # update redshift info
    if requires_convolution_bin_info and config["time_type"] == "redshift":
        config["convolution_redshift_bin_centers"] = (
            config["convolution_redshift_bin_edges"][1:]
            + config["convolution_redshift_bin_edges"][:-1]
        ) / 2
        config["convolution_redshift_bin_sizes"] = np.diff(
            config["convolution_redshift_bin_edges"]
        )

        config["convolution_time_bin_edges"] = config["convolution_redshift_bin_edges"]
        config["convolution_time_bin_centers"] = config[
            "convolution_redshift_bin_centers"
        ]
        config["convolution_time_bin_sizes"] = config["convolution_redshift_bin_sizes"]

        config["logger"].debug(
            "Updated bin data: convolution_redshift_bin_centers: {} convolution_redshift_bin_sizes: {}".format(
                config["convolution_redshift_bin_centers"],
                config["convolution_redshift_bin_sizes"],
            )
        )

    return config


def check_convolution_config_general(config):  # DH0001
    """
    Function to handle the general checking of convolution config input.

    This will not check the sfr-dicts, the convolution-instructions, or perform complicated logical checks
    """

    ##########
    # Skip the convolution
    if not config["check_convolution_config"]:
        return

    ##########
    # from the main dictionary, create a validation scheme
    validation_dict = {
        key: value["validation"]
        for key, value in default_convolution_config_dict.items()
        if "validation" in value
    }
    validation_schema = vol.Schema(validation_dict, extra=vol.ALLOW_EXTRA)

    ##########
    # do the validation
    for parameter, parameter_dict in config.items():

        ##########
        # Custom rules. we can decide to skip checking the input on some occasions

        # Skip checking the interpolator
        if parameter == "redshift_interpolator_data_output_filename":
            if config["time_type"] != "redshift":
                continue

        # skip these parameters in general for now
        if parameter in [
            "convolution_redshift_bin_edges",
            "convolution_lookback_time_bin_edges",
        ]:
            continue

        #
        validation_schema({parameter: parameter_dict})


def check_convolution_config_other(config):  # DH0001
    """
    Function to perform some extra custom checks on the convolution
    config. This can contain more complicated logic
    """

    # determine whether any of the config requires
    requires_convolution_bin_info = any(
        [
            convolution_instruction["convolution_type"] == "integrate"
            for convolution_instruction in config["convolution_instructions"]
        ]
    )

    if requires_convolution_bin_info:

        # extract time-type from general config
        time_type = config["time_type"]

        ######
        # Check the convolution time or redshift: TODO: this relies on convolution_type as well.
        if time_type == "lookback_time":
            # check if that is present
            if config.get("convolution_lookback_time_bin_edges", None) is None:
                raise ValueError(
                    "Please provide 'convolution_lookback_time_bin_edges' when using 'lookback-time' as 'time-type'"
                )

            if not is_time_unit(config["convolution_lookback_time_bin_edges"]):
                # if not config["convolution_lookback_time_bin_edges"].unit == u.yr:
                raise ValueError(
                    "Please express 'convolution_lookback_time_bin_edges' in units of time"
                )

            config["convolution_time_bin_edges"] = config[
                "convolution_lookback_time_bin_edges"
            ]
        elif time_type == "redshift":
            if config.get("convolution_redshift_bin_edges", None) is None:
                raise ValueError(
                    "Please provide 'convolution_redshift_bin_edges' when using 'redshift' as 'time-type'"
                )
            config["convolution_time_bin_edges"] = config[
                "convolution_redshift_bin_edges"
            ]
        else:
            raise ValueError("unsupported time-type")


def check_convolution_config(config):
    """
    Function to handle checking the convolution config

    TODO: move this to another file
    """

    #
    config["logger"].debug("Checking configuration")

    ############
    # General check of the convolution config. Will not check
    check_convolution_config_general(config=config)

    #######
    # check the convolution instructions
    check_and_update_convolution_instructions(convolution_config=config)

    #######
    # check the sfr dicts
    check_and_update_sfr_dicts(config=config)

    ##########
    # Perform other custom checks
    check_convolution_config_other(config=config)


def check_and_update_convolution_config(config):  # DH0001
    """
    Main function to check and update the convolution config
    """

    # check convolution config
    check_convolution_config(config=config)

    # update
    update_convolution_config(config=config)
