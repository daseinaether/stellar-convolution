"""
File containing a selection of functions to calculate the star formation rates.

TODO: reconsider if this whole padding is really necessary
"""

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve.calculate_birth_redshift_array import (
    calculate_origin_redshift_array,
)
from syntheticstellarpopconvolve.convolve_binned_data import calculate_overlap_fractions


def general_sfr_digitise_function(
    config, sfr_dict, time_values, metallicity_values=None
):
    """
    General function to handle the selection of SFR values given time values and metallicity values
    """

    # handle value extraction for non-redshift
    time_values_ = time_values
    padded_time_bin_edges = sfr_dict["padded_time_bin_edges"]
    if config["time_type"] != "redshift":
        time_values_ = time_values_.value
        padded_time_bin_edges = padded_time_bin_edges.value

    ################
    # Calculate time indices and determine the SFR values
    time_indices = (
        np.digitize(time_values_, bins=padded_time_bin_edges, right=False) - 1
    )

    # retrieve SFR values
    starformation_values = sfr_dict["padded_starformation_rate_array"][time_indices]

    ###
    # Handle whether we want to select on metallicity too
    if metallicity_values is not None:

        # Get indices for metallicity values
        metallicity_indices = (
            np.digitize(
                metallicity_values,
                bins=sfr_dict["padded_metallicity_bin_edges"],
                right=False,
            )
            - 1
        )

        # Calculate rates
        starformation_values = sfr_dict[
            "padded_metallicity_weighted_starformation_rate_array"
        ][time_indices, metallicity_indices]

    return starformation_values


def calculate_origin_time_array(config, data_dict, convolution_time_bin_center):
    """
    Function to calculate the origin time array
    """

    config["logger"].debug("Calculating origin-time array")

    # if convolution method and SFR is the in lookback time, then we can just subtract
    if config["time_type"] == "lookback_time":
        origin_time_array = (
            np.ones(data_dict["delay_time"].shape) * convolution_time_bin_center
            + data_dict["delay_time"]
        )
        config["logger"].debug(
            "Calculating origin-time array based on lookback_time: {}".format(
                origin_time_array
            )
        )
    elif config["time_type"] == "redshift":
        origin_time_array = calculate_origin_redshift_array(
            config=config,
            convolution_redshift_value=convolution_time_bin_center,
            data_dict=data_dict,
        )
        config["logger"].debug(
            "Calculating origin-time array based on redshift: {}".format(
                origin_time_array
            )
        )
    else:
        raise ValueError("Choice for time-type unknown. {}".format(config["time_type"]))

    return origin_time_array


def calculate_digitized_sfr_rates_for_forward_convolution(
    config, convolution_instruction, sfr_dict, data_dict, time_bin_info_dict
):
    """
    Function to calculate the total starformation occuring in a particular time
    bin

    if metallicity information is not required, this yields a scalar value

    if it is required, this yields a vector with values matching
    `total_star_formation_mass * (dP/dZ_{j})*dZ_{j}` where Z_{j} is the
    metallicity-bin in which the system falls

    TODO: abstract the actual SFR rate sampling
    TODO: move tde docstrings somewhere else
    """

    #########
    # backward convolution
    if convolution_instruction["convolution_direction"] == "backward":
        raise ValueError(
            "Currently backward convolution by sampling with unbinned data is not supported"
        )

    #
    total_star_formation_in_lookback_time_bin = general_sfr_digitise_function(
        config=config,
        sfr_dict=sfr_dict,
        time_values=time_bin_info_dict["bin_center"],
        metallicity_values=(
            data_dict["metallicity"] if "metallicity" in data_dict else None
        ),
    )

    #########

    #
    lookback_time_bin_size = time_bin_info_dict["bin_size"]
    lookback_time_bin_lower_edge = time_bin_info_dict["bin_edge_lower"]

    # multiply by binsize
    if convolution_instruction["multiply_by_sfr_time_binsize"]:
        #
        total_star_formation_in_lookback_time_bin = (
            total_star_formation_in_lookback_time_bin * lookback_time_bin_size
        )

    #
    config["logger"].warning(
        "Lower time bin {} upper time bin {} total mass formed {}".format(
            lookback_time_bin_lower_edge,
            lookback_time_bin_lower_edge + lookback_time_bin_size,
            total_star_formation_in_lookback_time_bin,
        )
    )

    return total_star_formation_in_lookback_time_bin


def calculate_digitized_sfr_rates_binned_data_for_backward_convolution(
    config,
    convolution_instruction,
    convolution_time_bin_center,
    data_dict,
    sfr_dict,
    delay_time_data_bin_info_dict,
):
    """
    Function to handle convolving binned data

    This function performs the following steps:
    - sets up the shifted data-bin edges
    - loops over each left-right edge pair and determines which SFR bins that
      edge-pair spans/overlaps with, the fractional overlap etc. for each
    - left-right edge pair, loops over the overlapping SFR bins

    NOTE: does not support redshift-based convolution
    """

    ############################
    # Unpack and set up

    #
    if config["time_type"] != "lookback_time":
        raise ValueError("Time-type must be `lookback_time`")

    #
    config["logger"].debug(
        "Convolving the binned data for convolution_time_bin_center: {}".format(
            convolution_time_bin_center
        )
    )

    local_system_indices = np.arange(len(data_dict["delay_time_data_bin_index"]))

    ####
    # unpack time bin info
    delay_time_data_bin_edges = delay_time_data_bin_info_dict[
        "delay_time_data_bin_edges"
    ].to(u.yr)
    delay_time_data_bin_sizes = np.diff(delay_time_data_bin_edges)

    config["logger"].info(
        "delay_time_data_bin_edges: {}\ndelay_time_data_bin_sizes: {}".format(
            delay_time_data_bin_edges, delay_time_data_bin_sizes
        )
    )

    ####
    # Get left and right edges
    left_delay_time_data_bin_edges = delay_time_data_bin_edges[:-1]
    right_delay_time_data_bin_edges = delay_time_data_bin_edges[1:]

    config["logger"].info(
        "left_delay_time_data_bin_edges: {}\nright_delay_time_data_bin_edges: {}".format(
            left_delay_time_data_bin_edges, right_delay_time_data_bin_edges
        )
    )

    ####
    # shift time bin data
    shifted_left_delay_time_data_bin_edges = (
        left_delay_time_data_bin_edges + convolution_time_bin_center
    )
    shifted_right_delay_time_data_bin_edges = (
        right_delay_time_data_bin_edges + convolution_time_bin_center
    )

    config["logger"].info(
        "shifted_left_delay_time_data_bin_edges: {}\nshifted_right_delay_time_data_bin_edges: {}".format(
            shifted_left_delay_time_data_bin_edges,
            shifted_right_delay_time_data_bin_edges,
        )
    )

    ####
    # read out sfr info
    sfr_bin_edges = sfr_dict["time_bin_edges"].to(u.yr)
    sfr_bin_sizes = np.diff(sfr_bin_edges)

    config["logger"].info(
        "sfr_bin_edges: {}\nsfr_bin_sizes: {}".format(sfr_bin_edges, sfr_bin_sizes)
    )

    ####
    # Set up empty sfr rates
    sfr_rates = (
        np.zeros(len(data_dict["delay_time_data_bin_index"]))
        * sfr_dict["starformation_rate_array"].unit
    )
    if convolution_instruction["multiply_by_sfr_time_binsize"]:
        sfr_rates = sfr_rates * delay_time_data_bin_sizes.unit

    config["logger"].info("$$$$$$$$$$$$$$$$$$$$$$$$$$")

    ##########
    # Loop over the data time-bin edge pairs
    #  and calculate the fraction of overlap of this edge pair with the SFR bins
    #  and for each SFR bin that they overlap with, calculate the SFR rates
    #  and store these and
    for delay_time_data_bin_i, (
        delay_time_data_bin_size_i,
        shifted_left_delay_time_data_bin_edge,
        shifted_right_delay_time_data_bin_edge,
    ) in enumerate(
        list(
            zip(
                delay_time_data_bin_sizes,
                shifted_left_delay_time_data_bin_edges,
                shifted_right_delay_time_data_bin_edges,
            )
        )
    ):

        #
        config["logger"].info(
            "delay_time_data_bin_i: {}\ndelay_time_data_bin_size_i: {}".format(
                delay_time_data_bin_i, delay_time_data_bin_size_i
            )
        )
        config["logger"].info(
            "shifted_left_delay_time_data_bin_edge: {}\nshifted_right_delay_time_data_bin_edge: {}".format(
                shifted_left_delay_time_data_bin_edge,
                shifted_right_delay_time_data_bin_edge,
            )
        )

        #########
        # check if we extend beyond or below all of the sfr bins
        if shifted_left_delay_time_data_bin_edge >= sfr_bin_edges[-1]:
            config["logger"].warning(
                "left-most delay time bin edge {} extends beyond the rightmost sfr bin edge: {}. skipping current delay time bin and breaking this loop.".format(
                    shifted_left_delay_time_data_bin_edge, sfr_bin_edges[-1]
                )
            )
            break
        if shifted_right_delay_time_data_bin_edge <= sfr_bin_edges[0]:
            config["logger"].warning(
                "right-most delay time bin edge {} extends below the leftmost sfr bin edge: {}. skipping current delay time bin and breaking this loop.".format(
                    shifted_right_delay_time_data_bin_edge, sfr_bin_edges[0]
                )
            )
            break

        #########
        # Calculate/determine the indices of the systems matching the current delay-time data bin index
        matching_delay_time_data_bin_system_indices = local_system_indices[
            data_dict["delay_time_data_bin_index"] == delay_time_data_bin_i
        ]
        config["logger"].info(
            "matching_delay_time_data_bin_system_indices: {}".format(
                matching_delay_time_data_bin_system_indices
            )
        )

        #########
        # Determine bin overlap fractions
        overlap_fractions = calculate_overlap_fractions(
            shifted_left_delay_time_data_bin_edge=shifted_left_delay_time_data_bin_edge,
            shifted_right_delay_time_data_bin_edge=shifted_right_delay_time_data_bin_edge,
            sfr_bin_sizes=sfr_bin_sizes,
            sfr_bin_edges=sfr_bin_edges,
        )

        config["logger"].info(
            "overlap_fractions:\n{}".format(
                "\n\t".join(
                    [
                        "{}: {}".format(key, value)
                        for key, value in overlap_fractions.items()
                    ]
                )
            )
        )

        ###########
        # loop over overlapping sfr bins
        # - Using the overlap-fraction dict information we can loop over the
        #   SFR bins and fetch the rates for the relevant systems
        #
        # TODO: the code below loops over the non-zero overlap bin indices, gets the SFR value and stores it.
        #   This can be done in 1 vector operation.
        config["logger"].info("=========================")

        combined_matching_delay_time_data_bin_sfr_rates = (
            np.zeros(matching_delay_time_data_bin_system_indices.shape)
            * sfr_dict["starformation_rate_array"].unit
            * sfr_bin_sizes.unit
        )

        #
        for overlap_bin_i, sfr_bin_index in enumerate(
            overlap_fractions["non_zero_overlap_with_sfr_bins"]
        ):
            config["logger"].info(
                "sfr_bin_index: {} overlap_bin_i: {}".format(
                    sfr_bin_index, overlap_bin_i
                )
            )

            # Automatic method
            matching_delay_time_data_bin_sfr_rates = general_sfr_digitise_function(
                config=config,
                sfr_dict=sfr_dict,
                time_values=sfr_dict["time_bin_centers"][
                    np.repeat(
                        sfr_bin_index, matching_delay_time_data_bin_system_indices.shape
                    )
                ],
                metallicity_values=(
                    data_dict["metallicity"][
                        matching_delay_time_data_bin_system_indices
                    ]
                    if "metallicity" in data_dict
                    else None
                ),
            )

            #
            config["logger"].info(
                "matching_delay_time_data_bin_sfr_rates: {}".format(
                    matching_delay_time_data_bin_sfr_rates
                )
            )

            #######
            # Weight the rates properly

            # Multiply by the fraction that the data time-bin overlaps
            weighted_matching_delay_time_data_bin_sfr_rates = (
                matching_delay_time_data_bin_sfr_rates
                * overlap_fractions["normalized_combined_overlap_array"][sfr_bin_index]
            )

            config["logger"].info(
                "Weighing the SFR rates with the fraction of overlap of bin: {}: {}".format(
                    overlap_fractions["normalized_combined_overlap_array"][
                        sfr_bin_index
                    ],
                    weighted_matching_delay_time_data_bin_sfr_rates,
                )
            )

            # Multiply by the width of the sfr bin
            weighted_matching_delay_time_data_bin_sfr_rates *= sfr_bin_sizes[
                sfr_bin_index
            ]

            config["logger"].info(
                "Multiplying the SFR rates matching SFR bin size: {}: {}".format(
                    sfr_bin_sizes[sfr_bin_index],
                    weighted_matching_delay_time_data_bin_sfr_rates,
                )
            )

            ########
            # Store the data in the combined array
            combined_matching_delay_time_data_bin_sfr_rates += (
                weighted_matching_delay_time_data_bin_sfr_rates
            )

            config["logger"].info(
                "Added the local rates to the combined rates of this data-time bin: {}".format(
                    combined_matching_delay_time_data_bin_sfr_rates
                )
            )

        #################
        #
        config["logger"].info(
            "Finished looping over overlapping SFR bins. Generated combined_matching_delay_time_data_bin_sfr_rates {}.\n Finalising the rate calculation".format(
                combined_matching_delay_time_data_bin_sfr_rates
            )
        )

        #############
        # Calculate capped time bin size. The right edge of the time-bin can extend beyond the final SFR bin
        capped_delay_time_data_bin_size_i = delay_time_data_bin_size_i
        sum_overlapping_sfr_bin_size = np.sum(
            overlap_fractions["normalized_combined_overlap_array"] * sfr_bin_sizes
        )
        if sum_overlapping_sfr_bin_size < capped_delay_time_data_bin_size_i:
            capped_delay_time_data_bin_size_i = sum_overlapping_sfr_bin_size
            config["logger"].info(
                "Capped delay-time data bin normalisition width to {}.".format(
                    capped_delay_time_data_bin_size_i
                )
            )

        #############
        # re-weight them to make average starformation rate
        combined_matching_delay_time_data_bin_sfr_rates /= (
            capped_delay_time_data_bin_size_i
        )
        config["logger"].info(
            "Divided sfr rates with {} to {}".format(
                capped_delay_time_data_bin_size_i,
                combined_matching_delay_time_data_bin_sfr_rates,
            )
        )

        #########
        # multiply by data time-bin if we to multiply by bin size
        if convolution_instruction["multiply_by_sfr_time_binsize"]:
            combined_matching_delay_time_data_bin_sfr_rates *= (
                capped_delay_time_data_bin_size_i
            )
            config["logger"].info(
                "Multiplying the rates by capped data-time binsize {} to {}".format(
                    capped_delay_time_data_bin_size_i,
                    combined_matching_delay_time_data_bin_sfr_rates,
                )
            )

        #########
        # store data in grand array
        sfr_rates[matching_delay_time_data_bin_system_indices] = (
            combined_matching_delay_time_data_bin_sfr_rates
        )

    config["logger"].debug(
        "Handled convolution of binned data at convolution bin-center {}".format(
            convolution_time_bin_center
        )
    )
    config["logger"].info(
        "Final sfr_rates for convolution bin-center {}: {}".format(
            convolution_time_bin_center, sfr_rates
        )
    )

    return sfr_rates


def calculate_digitized_sfr_rates_non_binned_data_for_backward_convolution(
    config, convolution_instruction, convolution_time_bin_center, data_dict, sfr_dict
):
    """ """

    ###########
    # calculate origin time
    origin_time_array = calculate_origin_time_array(
        config=config,
        data_dict=data_dict,
        convolution_time_bin_center=convolution_time_bin_center,
    )

    #########
    #
    digitised_sfr_rates = general_sfr_digitise_function(
        config=config,
        sfr_dict=sfr_dict,
        time_values=origin_time_array,
        metallicity_values=(
            data_dict["metallicity"] if "metallicity" in data_dict else None
        ),
    )

    ###################
    # Handle multiplication by sfr bin-size
    # TODO: clean and handle implementation
    # TODO: perhaps this can also be put into the sfr calculation function
    if convolution_instruction["multiply_by_sfr_time_binsize"]:

        # get indices
        time_binsize_indices = (
            np.digitize(
                origin_time_array, bins=sfr_dict["padded_time_bin_edges"], right=False
            )
            - 1
        )

        # get time-binsizes
        time_binsizes = sfr_dict["padded_time_bin_sizes"]

        # update sfr_rates
        digitised_sfr_rates = digitised_sfr_rates * time_binsizes[time_binsize_indices]

    return digitised_sfr_rates


def calculate_starformation(  # DH0001
    config, convolution_instruction, data_dict, sfr_dict, time_bin_info_dict
):
    """
    Main function that handles choices for starformation calculation
    """

    if convolution_instruction["convolution_direction"] == "backward":
        # with backward sampling the star formation for binned data needs to perform a weighted averaging over the SFR bins
        if convolution_instruction["contains_binned_data"]:
            starformation = (
                calculate_digitized_sfr_rates_binned_data_for_backward_convolution(
                    config=config,
                    convolution_instruction=convolution_instruction,
                    convolution_time_bin_center=time_bin_info_dict["bin_center"],
                    data_dict=data_dict,
                    sfr_dict=sfr_dict,
                    delay_time_data_bin_info_dict=convolution_instruction[
                        "delay_time_data_bin_info_dict"
                    ],
                )
            )
        # otherwise we just find out the starformation at the exact birth time of the system given the convolution time and the delay time.
        else:
            starformation = (
                calculate_digitized_sfr_rates_non_binned_data_for_backward_convolution(
                    config=config,
                    convolution_instruction=convolution_instruction,
                    convolution_time_bin_center=time_bin_info_dict["bin_center"],
                    data_dict=data_dict,
                    sfr_dict=sfr_dict,
                )
            )
    elif convolution_instruction["convolution_direction"] == "forward":
        # forward sampling just takes the value in the current bin
        starformation = calculate_digitized_sfr_rates_for_forward_convolution(
            config=config,
            convolution_instruction=convolution_instruction,
            sfr_dict=sfr_dict,
            data_dict=data_dict,
            time_bin_info_dict=time_bin_info_dict,
        )
    else:
        raise ValueError("convolution direction not supported")

    return starformation
