"""Routines for convolution-by-sampling

initial idea with simple situation

sfr [Msun /yr], global
fixed Z

partially grid-like

with a SFR evaluated in lookback time in bins, with t_l,i lookback times and
dt_l,i binsize and edges t_l,i-0.5, t_l,i+0.5

in a given bin we have the total mass formed in stars sfr(t_l=t_l,i) * dt_l,i
= m_tot,i

now, we have some systems of interest (e.g. dwd), gained through pop-synth
simulations. these have a.o. the property normalized yield, i..e number per
formed solar mass

Y_j [Msun]

total number of system j sampled;
Y_j * M_tot,i = N_j

if N_j > 1:
- take X systems where X = floor(N_j)
- N_j-x is then < 1
- take random number from uniform dist, P. if P < N_j-x: accept, else not

then we have a bunch of systems (which can include the same system)
but in that array, assign random lookback time between the bin edges

assign radnom position

- this sampling stategy can be multiprocssed easily (lookback time bins) can
- also easily be extended to include metallicity naturally handles unequal
- yield per systems

Notes:
- this method does not turn things around like the others do. We start at a
given lookback time bin for all systems. We sample a set of systems based on
the total starformation within that lookback time bin, and the normalized
yields of the systems. We then assign a birth lookback time to the systems
(taken randomly between the bin edges)
"""

import numpy as np

from syntheticstellarpopconvolve.post_convolution_hook_routines import (
    handle_post_convolution_function,
)


def convolution_by_sampling_post_convolution_hook_wrapper(
    config,
    sfr_dict,
    data_dict,
    time_bin_info_dict,
    convolution_instruction,
    convolution_results,
    #
    persistent_data=None,
    previous_convolution_results=None,
):
    """
    Function to wrap the post-convolution function call for event-convolution by sampling.

    rules:
    - additional data can be added to the result_dict
    - the number of systems can lower than before the call
    """

    #
    name = "convolution by sampling"

    #
    config["logger"].warning(
        "Handling post-convolution function hook call for {}".format(name)
    )

    #############
    # call hook
    convolution_results = handle_post_convolution_function(
        config=config,
        sfr_dict=sfr_dict,
        data_dict=data_dict,
        time_bin_info_dict=time_bin_info_dict,
        convolution_instruction=convolution_instruction,
        convolution_results=convolution_results,
        name=name,
        #
        persistent_data=persistent_data,
        previous_convolution_results=previous_convolution_results,
    )

    return convolution_results


def select_dict_entries_with_new_indices(sampled_data_dict, new_indices):
    """
    Function to select dict entires with new indices
    """

    sampled_data_dict = {
        data_key: sampled_data_dict[data_key][new_indices]
        for data_key in sampled_data_dict.keys()
        if not data_key == "name"
    }

    return sampled_data_dict


def sample_systems(
    yield_array,
    lookback_time_bin_size,
    lookback_time_bin_lower_edge,
    config,
    convolution_instruction,
):
    """
    General function to handle sampling a set of systems based on
    normalized yields and a total mass of stars formed
    """

    ###########
    #
    config["logger"].warning(
        "Sampling systems. Using yield array {}".format(yield_array)
    )

    ############
    #
    local_indices = np.arange(yield_array.shape[0])

    ############
    # first sample systems that have a should form at least one time, but only
    # the down-rounded number of times
    integer_formations = np.array(np.floor(yield_array), dtype=int)
    integer_sampled_formation_indices = np.repeat(local_indices, integer_formations)

    ############
    # then sample using the remainder (all parts with number between 0 and 1)

    # select the remainder
    fractional_formations = yield_array - integer_formations

    # take a random set to sample the fractional formations
    random_chance = np.random.random(fractional_formations.shape)

    # Sample the indices
    fractional_sampled_formation_indices = local_indices[
        random_chance < fractional_formations
    ]

    ############
    # Combine the sampled indice
    combined_sampled_indices = np.sort(
        np.concatenate(
            [integer_sampled_formation_indices, fractional_sampled_formation_indices]
        )
    )

    ############
    # Construct the payload
    convolution_results = {"sampled_indices": combined_sampled_indices}

    ############
    #
    config["logger"].warning(
        "Sampled {} systems.".format(combined_sampled_indices.shape)
    )

    return convolution_results
