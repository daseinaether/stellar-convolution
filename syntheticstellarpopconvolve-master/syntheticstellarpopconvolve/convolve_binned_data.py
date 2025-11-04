"""
Functions to support convolution of binned data.

Binned data is a new version of ensemble-based data, but inflated. Because
that data usually is binned, and the time-bin spans some range rather than a
distinct point in time, we should provide support for that.

Routine that handles calculating the overlap of a time-bin series with the starformation rate bins
- calculates distances from edges
- calculates overlap fraction of sfr bins
- calculates fraction of time

TODO: allow using CDF to re-scale. Data within the bin may not be distributed uniformly per se.
"""

import numpy as np


def calculate_overlap_fractions(
    shifted_left_delay_time_data_bin_edge,
    shifted_right_delay_time_data_bin_edge,
    sfr_bin_sizes,
    sfr_bin_edges,
):
    """
    Function to calculate the overlap

    TODO: consider returning only the information of the bins that overlap.
    """

    assert (
        shifted_left_delay_time_data_bin_edge < shifted_right_delay_time_data_bin_edge
    )

    delay_time_bin_size = (
        shifted_right_delay_time_data_bin_edge - shifted_left_delay_time_data_bin_edge
    )

    ##############
    # calculate distances
    left_distances = sfr_bin_edges[1:] - shifted_left_delay_time_data_bin_edge
    right_distances = shifted_right_delay_time_data_bin_edge - sfr_bin_edges[:-1]

    ##############
    # mask by negatives
    left_distances[left_distances < 0] = 0
    right_distances[right_distances < 0] = 0

    ##############
    # Construct the combined overlap array
    combined_overlap_array = sfr_bin_sizes.astype(float)
    combined_overlap_array[left_distances == 0] = 0
    combined_overlap_array[right_distances == 0] = 0

    # Handle non-zero entries
    left_nonzero_indices = np.nonzero(left_distances)[0]
    right_nonzero_indices = np.nonzero(right_distances)[0]

    # TODO: fix situations in which there are no non-zero entries

    #
    leftmost_nonzero_index = left_nonzero_indices[0]
    rightmost_nonzero_index = right_nonzero_indices[-1]

    # If they are the same, that means they both fall in the same bin
    if leftmost_nonzero_index == rightmost_nonzero_index:
        nonzero_index = leftmost_nonzero_index  # they're the same

        # if both lie in the same bin, then its just the distance between the two data-time bin edges
        combined_overlap_array[nonzero_index] = (
            shifted_right_delay_time_data_bin_edge
            - shifted_left_delay_time_data_bin_edge
        )
    else:
        combined_overlap_array[leftmost_nonzero_index] = left_distances[
            leftmost_nonzero_index
        ]
        combined_overlap_array[rightmost_nonzero_index] = right_distances[
            rightmost_nonzero_index
        ]

    ##############
    # normalize to fraction of the sfr bin
    normalized_combined_overlap_array = combined_overlap_array / sfr_bin_sizes

    ##############
    # get fraction of time-bin
    time_bin_fraction = combined_overlap_array / delay_time_bin_size

    ##############
    # calculate cumulative fraction to allow re-weighting with in-bin expected distribution
    cumulative_time_bin_fraction = np.cumsum(time_bin_fraction)
    change_cumulative_time_bin_fraction = np.diff(
        cumulative_time_bin_fraction
    )  # TODO: this input should be appended with 0 because essentially we're missing that now.

    ##############
    # non-zero sfr-bins. NOTE: these are the indices we should loop over
    non_zero_overlap_with_sfr_bins = np.nonzero(normalized_combined_overlap_array)[0]

    return {
        "combined_overlap_array": combined_overlap_array,
        "normalized_combined_overlap_array": normalized_combined_overlap_array,
        "time_bin_fraction": time_bin_fraction,
        "cumulative_time_bin_fraction": cumulative_time_bin_fraction,
        "change_cumulative_time_bin_fraction": change_cumulative_time_bin_fraction,
        "non_zero_overlap_with_sfr_bins": non_zero_overlap_with_sfr_bins,
    }


if __name__ == "__main__":

    #
    delay_time_data_bin_info = {
        "delay_time_data_bin_edges": np.arange(0, 2, 1),
    }

    #
    shift = 5.6

    sfr_bin_edges = np.arange(0, 20, 5)
    sfr_bin_sizes = np.diff(sfr_bin_edges)

    delay_time_data_bin_edges = delay_time_data_bin_info["delay_time_data_bin_edges"]
    delay_time_data_bin_sizes = np.diff(delay_time_data_bin_edges)

    #
    left_delay_time_data_bin_edges = delay_time_data_bin_edges[:-1]
    right_delay_time_data_bin_edges = delay_time_data_bin_edges[1:]

    shifted_left_delay_time_data_bin_edges = left_delay_time_data_bin_edges + shift
    shifted_right_delay_time_data_bin_edges = right_delay_time_data_bin_edges + shift

    ##########
    # Loop over the data time-bins
    for time_bin_i, (
        time_bin_size_i,
        shifted_left_delay_time_data_bin_edge,
        shifted_right_delay_time_data_bin_edge,
    ) in enumerate(
        list(
            zip(
                delay_time_data_bin_sizes,
                shifted_left_delay_time_data_bin_edges,
                shifted_right_delay_time_data_bin_edges,
            )
        )[:1]
    ):

        #
        overlap_fractions = calculate_overlap_fractions(
            shifted_left_delay_time_data_bin_edge=shifted_left_delay_time_data_bin_edge,
            shifted_right_delay_time_data_bin_edge=shifted_right_delay_time_data_bin_edge,
            sfr_bin_sizes=sfr_bin_sizes,
            sfr_bin_edges=sfr_bin_edges,
        )
