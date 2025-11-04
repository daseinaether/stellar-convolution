"""
Functions to check and update the SFR dict

TODO: add checks for the dimensions of the input arrays
"""

import numpy as np

from syntheticstellarpopconvolve.cosmology_utils import redshift_to_lookback_time
from syntheticstellarpopconvolve.general_functions import (
    has_unit,
    is_time_unit,
    pad_function,
)
from syntheticstellarpopconvolve.store_redshift_shell_info import (
    store_redshift_shell_info,
)


def pad_sfr_dict(config, sfr_dict):
    """
    Function to pad the entries in the sfr dictionary with empty bins.

    These functions update all the sfr properties and adds new entries that are prepended with 'padded_'
    """

    #
    config["logger"].debug("Padding SFR dict")

    max_pad = 1.0e13

    ##########
    # pad lookback time/redshift array
    if config["time_type"] == "lookback_time":
        # pad lookback time bins
        sfr_dict["padded_lookback_time_bin_edges"] = pad_function(
            array=sfr_dict["lookback_time_bin_edges"],
            left_val=-max_pad,
            right_val=max_pad,
            relative_to_edge_val=True,
        )
        sfr_dict["padded_time_bin_edges"] = sfr_dict["padded_lookback_time_bin_edges"]
        sfr_dict["time_bin_edges"] = sfr_dict["lookback_time_bin_edges"]

        #
        config["logger"].debug(
            "Padded lookback time bin edges {} to {}".format(
                sfr_dict["lookback_time_bin_edges"],
                sfr_dict["padded_lookback_time_bin_edges"],
            )
        )

        # add bin sizes
        sfr_dict["lookback_time_bin_sizes"] = np.abs(
            np.diff(sfr_dict["lookback_time_bin_edges"])
        )
        sfr_dict["time_bin_sizes"] = sfr_dict["lookback_time_bin_sizes"]

        # Pad time-bin sizes
        sfr_dict["padded_lookback_time_bin_sizes"] = pad_function(
            array=sfr_dict["lookback_time_bin_sizes"],
            left_val=0,
            right_val=0,
            relative_to_edge_val=False,
        )
        sfr_dict["padded_time_bin_sizes"] = sfr_dict["padded_lookback_time_bin_sizes"]

        # log the binsizes
        config["logger"].debug(
            "Created lookback time bin sizes {}".format(
                sfr_dict["lookback_time_bin_sizes"],
            )
        )
    elif config["time_type"] == "redshift":
        #
        sfr_dict["padded_redshift_bin_edges"] = pad_function(
            array=sfr_dict["redshift_bin_edges"],
            left_val=-max_pad,
            right_val=max_pad,
            relative_to_edge_val=True,
        )

        #
        sfr_dict["padded_time_bin_edges"] = sfr_dict["padded_redshift_bin_edges"]
        sfr_dict["time_bin_edges"] = sfr_dict["redshift_bin_edges"]

        #
        config["logger"].debug(
            "Padded redshift bin edges {} to {}".format(
                sfr_dict["redshift_bin_edges"],
                sfr_dict["padded_redshift_bin_edges"],
            )
        )

        # create redshift time-bin size
        sfr_dict["redshift_bin_sizes"] = np.abs(np.diff(sfr_dict["redshift_bin_edges"]))

        # TODO: put these steps in a dedicated function that makes an array of a list of astropy-united values
        lookback_time_at_redshift_bin_edges = [
            redshift_to_lookback_time(
                redshift=redshift_bin_edge, cosmology=config["cosmology"]
            )
            for redshift_bin_edge in sfr_dict["redshift_bin_edges"]
        ]
        lookback_time_at_redshift_bin_edges_unit = lookback_time_at_redshift_bin_edges[
            0
        ].unit
        lookback_time_at_redshift_bin_edges = [
            el.value for el in lookback_time_at_redshift_bin_edges
        ]
        lookback_time_at_redshift_bin_edges = (
            np.array(lookback_time_at_redshift_bin_edges)
            * lookback_time_at_redshift_bin_edges_unit
        )

        sfr_dict["time_bin_sizes"] = np.abs(
            np.diff(lookback_time_at_redshift_bin_edges)
        )

        # Pad time-bin sizes
        sfr_dict["padded_redshift_bin_sizes"] = pad_function(
            array=sfr_dict["redshift_bin_sizes"],
            left_val=0,
            right_val=0,
            relative_to_edge_val=False,
        )

        # log the binsizes
        config["logger"].debug(
            "Created redshift bin sizes {} and the corresponding lookback time-bin sizes {} ".format(
                sfr_dict["redshift_bin_sizes"], sfr_dict["time_bin_sizes"]
            )
        )
    else:
        raise ValueError("Invalid time-type ({})".format(config["time_type"]))

    #########
    # Pad time-bin sizes
    sfr_dict["padded_time_bin_sizes"] = pad_function(
        array=sfr_dict["time_bin_sizes"],
        left_val=0,
        right_val=0,
        relative_to_edge_val=False,
    )

    # log the binsizes
    config["logger"].debug(
        "Padded the time bin sizes {}".format(
            sfr_dict["padded_time_bin_sizes"],
        )
    )

    ##########
    # pad SFR rate array
    if "starformation_rate_array" in sfr_dict:  # it should be present always
        #
        sfr_dict["padded_starformation_rate_array"] = pad_function(
            array=sfr_dict["starformation_rate_array"],
            left_val=0,
            right_val=0,
            relative_to_edge_val=False,
        )

        #
        config["logger"].debug(
            "Padded starformation array {} to {}".format(
                sfr_dict["starformation_rate_array"],
                sfr_dict["padded_starformation_rate_array"],
            )
        )

    ##########
    # pad metallicity bins
    if "metallicity_bin_edges" in sfr_dict:
        #
        sfr_dict["padded_metallicity_bin_edges"] = pad_function(
            array=sfr_dict["metallicity_bin_edges"],
            left_val=-1e-20,
            right_val=2,
            relative_to_edge_val=False,
        )

        #
        config["logger"].debug(
            "Padded metallicity bin edges {} to {}".format(
                sfr_dict["metallicity_bin_edges"],
                sfr_dict["padded_metallicity_bin_edges"],
            )
        )

    ##########
    # pad metallicity distribution
    if "metallicity_distribution_array" in sfr_dict:
        #
        sfr_dict["padded_metallicity_distribution_array"] = pad_function(
            array=sfr_dict["metallicity_distribution_array"],
            left_val=0,
            right_val=0,
            relative_to_edge_val=False,
        )

        #
        sfr_dict["padded_metallicity_distribution_array"] = pad_function(
            array=sfr_dict["padded_metallicity_distribution_array"],
            left_val=0,
            right_val=0,
            relative_to_edge_val=False,
            axis=1,
        )

    ##########
    # pad metallicity weighted SFR rate bins
    if "metallicity_weighted_starformation_rate_array" in sfr_dict:
        #
        sfr_dict["padded_metallicity_weighted_starformation_rate_array"] = pad_function(
            array=sfr_dict["metallicity_weighted_starformation_rate_array"],
            left_val=0,
            right_val=0,
            relative_to_edge_val=False,
        )

        #
        sfr_dict["padded_metallicity_weighted_starformation_rate_array"] = pad_function(
            array=sfr_dict["padded_metallicity_weighted_starformation_rate_array"],
            left_val=0,
            right_val=0,
            relative_to_edge_val=False,
            axis=1,
        )

    return sfr_dict


def update_sfr_dict(sfr_dict, config):
    """
    Function to update the SFR dict
    - provides padding
    - adds redshift shell info
    - updates metallicity-based info
    """

    #
    config["logger"].debug("Updating SFR dict")

    sfr_dict["include_metallicity_info"] = False
    if "metallicity_distribution_array" in sfr_dict.keys():
        sfr_dict["include_metallicity_info"] = True

        # determine centers
        sfr_dict["metallicity_bin_centers"] = (
            sfr_dict["metallicity_bin_edges"][1:]
            + sfr_dict["metallicity_bin_edges"][:-1]
        ) / 2

        # determine sizes
        sfr_dict["metallicity_bin_sizes"] = (
            sfr_dict["metallicity_bin_edges"][1:]
            - +sfr_dict["metallicity_bin_edges"][:-1]
        )
        sfr_dict["metallicity_bin_sizes"] = np.diff(sfr_dict["metallicity_bin_edges"])

        # construct the combined array: multiplies the SFR, dp/dZ, and delta Z
        sfr_dict["metallicity_weighted_starformation_rate_array"] = (
            sfr_dict["starformation_rate_array"][:, np.newaxis]
            * sfr_dict["metallicity_bin_sizes"][np.newaxis, :]
            * sfr_dict["metallicity_distribution_array"]
        )

    # add bin centers
    if config["time_type"] == "lookback_time":
        sfr_dict["lookback_time_bin_centers"] = (
            sfr_dict["lookback_time_bin_edges"][1:]
            + sfr_dict["lookback_time_bin_edges"][:-1]
        ) / 2
        sfr_dict["time_bin_centers"] = sfr_dict["lookback_time_bin_centers"]
    elif config["time_type"] == "redshift":
        sfr_dict["redshift_bin_centers"] = (
            sfr_dict["redshift_bin_edges"][1:] + sfr_dict["redshift_bin_edges"][:-1]
        ) / 2
        sfr_dict["time_bin_centers"] = sfr_dict["redshift_bin_centers"]

    # Pad the SFR dict with the empty bins around
    sfr_dict = pad_sfr_dict(config=config, sfr_dict=sfr_dict)

    # Add redshift shell info to dict.
    sfr_dict = store_redshift_shell_info(config=config, sfr_dict=sfr_dict)

    return sfr_dict


def check_sfr_dict(
    sfr_dict, config, requires_name, requires_metallicity_info, time_type
):
    """
    Function to check the sfr dictionary
    """

    ##########
    # Check if the name exists if the sfr dict requires it
    if requires_name:
        if "name" not in sfr_dict:
            raise ValueError("Name is required in the sfr dictionary")

    ##########
    # Check if the correct time bins are present
    if "starformation_rate_array" not in sfr_dict:
        raise ValueError("starformation_rate_array is required in the sfr dictionary")

    # check if starformation array has any unit
    try:
        sfr_dict["starformation_rate_array"].unit
    except AttributeError:
        raise AttributeError("starformation_rate_array requires an astropy unit")

    ##########
    # Check if the correct time bins are present
    if time_type == "lookback_time":
        if "lookback_time_bin_edges" not in sfr_dict:
            raise ValueError(
                "lookback_time_bin_edges is required in the sfr dictionary"
            )

        # check if has time-units
        if not is_time_unit(sfr_dict["lookback_time_bin_edges"]):
            raise ValueError(
                "Please express 'lookback_time_bin_edges' in units of time"
            )

        # check if length is correct:
        if (
            not len(sfr_dict["lookback_time_bin_edges"])
            == len(sfr_dict["starformation_rate_array"]) + 1
        ):
            raise ValueError(
                "Please ensure the length of the `starformation_rate_array` ({}) is one element shorter than the length of `lookback_time_bin_edges` ({})".format(
                    len(sfr_dict["starformation_rate_array"]),
                    len(sfr_dict["lookback_time_bin_edges"]),
                )
            )

    elif time_type == "redshift":
        if "redshift_bin_edges" not in sfr_dict:
            raise ValueError("redshift_bin_edges is required in the sfr dictionary")

        # check if length is correct:
        if (
            not len(sfr_dict["redshift_bin_edges"])
            == len(sfr_dict["starformation_rate_array"]) + 1
        ):
            raise ValueError(
                "Please ensure the length of the `starformation_rate_array` ({}) is one element shorter than the length of `redshift_bin_edges` ({})".format(
                    len(sfr_dict["starformation_rate_array"]),
                    len(sfr_dict["redshift_bin_edges"]),
                )
            )

    ##########
    # check if metallicity information is present
    if requires_metallicity_info:
        # Check if the metallicity bins are present
        if "metallicity_bin_edges" not in sfr_dict:
            raise ValueError("metallicity_bin_edges is required in the SFR dictionary.")

        # Check for metallicity distribution
        # NOTE: from 2024-06-08 I have decided to require only the
        # "metallicity_distribution_array". "metallicity_weighted_starformation_rate_array"
        # will be created from this.
        if "metallicity_distribution_array" not in sfr_dict:
            raise ValueError(
                "metallicity_distribution_array is required in the sfr dictionary"
            )

        # check if the array has units. (not allowed)
        if has_unit(sfr_dict["metallicity_distribution_array"]):
            raise ValueError(
                "metallicity_distribution_array should not contain any units"
            )

        # check if length in the metallicity direction is correct:
        if (
            not sfr_dict["metallicity_distribution_array"].shape[1]
            == len(sfr_dict["metallicity_bin_edges"]) - 1
        ):
            raise ValueError(
                "Please ensure the length of the `metallicity_distribution_array.shape[1]` ({}) is one element shorter than the length of `metallicity_bin_edges` ({})".format(
                    sfr_dict["metallicity_distribution_array"].shape[1],
                    len(sfr_dict["metallicity_bin_edges"]),
                )
            )

        # check if length in the time direction is correct:
        if time_type == "lookback_time":
            if (
                not sfr_dict["metallicity_distribution_array"].shape[0]
                == len(sfr_dict["lookback_time_bin_edges"]) - 1
            ):
                raise ValueError(
                    "Please ensure the length of the metallicity_distribution_array.shape[0]` ({}) is one element shorter than the length of `lookback_time_bin_edges` ({})".format(
                        sfr_dict["metallicity_distribution_array"].shape[0],
                        len(sfr_dict["lookback_time_bin_edges"]),
                    )
                )

        elif time_type == "redshift":
            if (
                not sfr_dict["metallicity_distribution_array"].shape[0]
                == len(sfr_dict["redshift_bin_edges"]) - 1
            ):
                raise ValueError(
                    "Please ensure the length of the metallicity_distribution_array.shape[0]` ({}) is one element shorter than the length of `redshift_bin_edges` ({})".format(
                        sfr_dict["metallicity_distribution_array"].shape[0],
                        len(sfr_dict["redshift_bin_edges"]),
                    )
                )


def check_and_update_sfr_dict(  # DH0001
    sfr_dict, config, requires_name, requires_metallicity_info, time_type
):  # DH0001
    """
    Function to check the SFR dict for the appropriate content and update
    """

    # check sfr dict
    check_sfr_dict(
        sfr_dict=sfr_dict,
        config=config,
        requires_name=requires_name,
        requires_metallicity_info=requires_metallicity_info,
        time_type=time_type,
    )

    # update sfr dict
    sfr_dict = update_sfr_dict(sfr_dict=sfr_dict, config=config)

    return sfr_dict


def check_and_update_sfr_dicts(config):  # DH0001
    """
    Function to check the SFR dict for the appropriate content and update
    """

    # determine whether any of the convolution instructions use metallicity
    requires_metallicity_info = any(
        [
            "metallicity" in convolution_instruction["data_column_dict"].keys()
            for convolution_instruction in config["convolution_instructions"]
        ]
    )

    #######
    # check the SFR information
    if "SFR_info" not in config or not config["SFR_info"]:
        raise ValueError(
            'please provide a non-empty list or dictionary to config["SFR_INFO"]'
        )

    if isinstance(config["SFR_info"], dict):
        config["SFR_info"] = check_and_update_sfr_dict(
            sfr_dict=config["SFR_info"],
            requires_name=False,
            requires_metallicity_info=requires_metallicity_info,
            time_type=config["time_type"],
            config=config,
        )
    elif isinstance(config["SFR_info"], list):
        # check all sfr dicts
        for sfr_dict in config["SFR_info"]:
            sfr_dict = check_and_update_sfr_dict(
                sfr_dict=sfr_dict,
                requires_name=True,
                requires_metallicity_info=requires_metallicity_info,
                time_type=config["time_type"],
                config=config,
            )

    # TODO: handle a generator type
