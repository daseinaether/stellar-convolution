"""
Functions to convolve the T0 format with sampling

TODO: move the calculations to the post-convolution hook
TODO: determine which systems that are (at present day) in the lisa frequency range should have interacted through RLOF
TODO: of the systems that are not RLOFing and are within the lisa waveband, store: indices, source.f_orb_now. the rest can be retrieved elsewhere
"""

import copy
import json
import os
import time

import astropy.units as u
import h5py
import numpy as np
import pandas as pd

from syntheticstellarpopconvolve import convolve, default_convolution_config
from syntheticstellarpopconvolve.convolve_stochastically import (
    select_dict_entries_with_new_indices,
)
from syntheticstellarpopconvolve.general_functions import temp_dir


def get_mass_norm():
    return 4476544.539875359 * u.Msun


TMP_DIR = temp_dir("code", "convolve_stochastically", clean_path=True)

np.random.seed(0)


def post_convolution_function(
    config, job_dict, sfr_dict, data_dict, convolution_results, convolution_instruction
):
    """
    Post-convolution function to handle integrating the systems forward in time and finding those that end up in the LISA waveband.

    using local_indices to select everything and using Alexey's distance sampler to handle sampling the distances
    """

    # unpack data
    system_indices = convolution_results["indices"]
    local_indices = np.arange(len(system_indices))

    # sample distances
    dist = np.random.randint(0, 1e6, size=len(local_indices))

    #
    convolution_results["dists"] = dist * u.kpc

    # shuffle and select 2 sets
    shuffled_indices = np.arange(len(system_indices))
    np.random.shuffle(shuffled_indices)

    #
    random_length = np.random.randint(0, high=len(shuffled_indices))

    #
    shuffled_indices_1 = shuffled_indices[:random_length]
    shuffled_indices_2 = shuffled_indices[random_length:]

    # Select all results from set 1
    convolution_result_1 = select_dict_entries_with_new_indices(
        sampled_data_dict=convolution_results,
        new_indices=shuffled_indices_1,
    )
    convolution_result_1["name"] = "set_1"

    # Select all results from set 2
    convolution_result_2 = select_dict_entries_with_new_indices(
        sampled_data_dict=convolution_results,
        new_indices=shuffled_indices_2,
    )
    convolution_result_2["name"] = "set_2"

    # split into two
    convolution_results = [convolution_result_1, convolution_result_2]

    return convolution_results


def prepare_input_data():
    """
    Example function to prepare the input data
    """

    BinCodex_events_filename = (
        "/home/david/Desktop/bincodex_results/example_BinCodex.h5"
    )

    #
    BinCodex_T0_events = pd.read_hdf(
        BinCodex_events_filename,
        "T0",
    )

    ##################
    # update T0 output

    # get mass normalisation
    mass_normalisation_fiducial = get_mass_norm()

    # set normalised yield
    BinCodex_T0_events["normalized_yield"] = 1 / mass_normalisation_fiducial

    # Query the dataset to select the formation of the WDs

    # to check if things start with some number its easier to turn them into strings
    BinCodex_T0_events["str_event"] = BinCodex_T0_events["event"].astype(str)
    BinCodex_T0_events["str_type1"] = BinCodex_T0_events["type1"].astype(str)
    BinCodex_T0_events["str_type2"] = BinCodex_T0_events["type2"].astype(str)

    # first, lets query the type-changing events. Any type-change will do
    wd_binaries = BinCodex_T0_events.query("str_event.str.startswith('1')")

    # The type should change to a WD-type (and the other should already be one)
    wd_binaries = wd_binaries.query("str_type1.str.startswith('2')")
    wd_binaries = wd_binaries.query("str_type2.str.startswith('2')")

    # lets delete the string versions of the columns again
    wd_binaries = wd_binaries.drop(columns=["str_event", "str_type1", "str_type2"])

    # lets also delete the original dataframe
    del BinCodex_T0_events

    return wd_binaries


LIGHTWEIGHT = False

###################
# Read T0 output
start = time.time()
wd_binaries = prepare_input_data()

##################
#

# create file
input_hdf5_filename = os.path.join(TMP_DIR, "input_hdf5.h5")
output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5.h5")
input_hdf5_file = h5py.File(input_hdf5_filename, "w")

# Create groups main
input_hdf5_file.create_group("input_data")
input_hdf5_file.create_group("config")

# add group for events
input_hdf5_file.create_group("input_data/events")

# Write population config to file
input_hdf5_file.create_dataset("config/population", data=json.dumps({}))

# close
input_hdf5_file.close()

# store the data frame in the hdf5file
wd_binaries.to_hdf(input_hdf5_filename, key="input_data/events/stochastic_example")

#
convolution_config = copy.copy(default_convolution_config)
convolution_config["input_filename"] = input_hdf5_filename
convolution_config["output_filename"] = output_hdf5_filename
convolution_config["tmp_dir"] = TMP_DIR
convolution_config["redshift_interpolator_data_output_filename"] = os.path.join(
    TMP_DIR, "interpolator_dict.p"
)
convolution_config["multiply_by_time_binsize"] = False
convolution_config["filter_future_events"] = False


###
# convolution instructions
convolution_config["convolution_instructions"] = [
    {
        "convolution_type": "sample",
        "input_data_name": "stochastic_example",
        "output_data_name": "stochastic_example",
        "filter_future_events": False,
        "post_convolution_function": post_convolution_function,
        "data_column_dict": {
            # required
            "normalized_yield": "normalized_yield",
            "delay_time": {"column_name": "time", "unit": u.Myr},
        },
    },
]

#
convolution_config["time_type"] = "lookback_time"
convolution_config["convolution_lookback_time_bin_edges"] = np.arange(3, 6, 1) * u.Gyr
print(convolution_config["convolution_lookback_time_bin_edges"])
quit()

# construct the sfr-dict (NOTE: this uses absolute SFR, not metallicity dependent)
sfr_dict = {}
sfr_dict["lookback_time_bin_edges"] = (np.arange(0, 10, 1) * u.Gyr).to(u.yr)

#
scale = 1e-5
sfr_dict["starformation_rate_array"] = (
    scale * np.ones(sfr_dict["lookback_time_bin_edges"].shape[0] - 1) * u.Msun / u.yr
)  # example of a constant star-formation rate. this could be anything of course.

# store
convolution_config["SFR_info"] = sfr_dict

# convolve
convolve(config=convolution_config)

print("finished convolution")


# read out content and integrate until today
with h5py.File(convolution_config["output_filename"], "r") as output_hdf5_file:

    print(output_hdf5_file.keys())
    print(output_hdf5_file["output_data"].keys())
    print(output_hdf5_file["output_data/event"].keys())
    print(output_hdf5_file["output_data/event/stochastic_example"].keys())
    print(
        output_hdf5_file[
            "output_data/event/stochastic_example/stochastic_example"
        ].keys()
    )
    # print(
    #     output_hdf5_file[
    #         "output_data/event/stochastic_example/stochastic_example/convolution_results/1500000000.0 yr"
    #     ].keys()
    # )

    print(
        output_hdf5_file[
            "output_data/event/stochastic_example/stochastic_example/convolution_results/set_1"
        ].keys()
    )
    print(
        output_hdf5_file[
            "output_data/event/stochastic_example/stochastic_example/convolution_results/set_2"
        ].keys()
    )

    quit()

    print(
        output_hdf5_file[
            "output_data/event/stochastic_example/stochastic_example/convolution_results"
        ].keys()
    )

    formation_time_bin_keys = list(
        output_hdf5_file[
            "output_data/event/stochastic_example/stochastic_example/convolution_results"
        ].keys()
    )

    ################
    #

    # loop over the formation-time bins
    formation_time_bin_keys = sorted(
        formation_time_bin_keys, key=lambda x: float(x.split(" ")[0])
    )
    for formation_time_bin_key in formation_time_bin_keys:

        # formation_time_bin_key = "3500000000.0 yr"
        print("=================================")
        print(f"formation_time_bin_key: {formation_time_bin_key}")
        print("=================================")

        ###########
        # Read out data

        # convert units
        unit_dict = json.loads(
            output_hdf5_file[
                f"output_data/event/stochastic_example/stochastic_example/convolution_results/{formation_time_bin_key}"
            ].attrs["units"]
        )
        unit_dict = {key: u.Unit(val) for key, val in unit_dict.items()}
        print(unit_dict)


#         indices = output_hdf5_file[
#             f"output_data/event/stochastic_example/stochastic_example/convolution_results/{formation_time_bin_key}/indices"
#         ][()]
#         # print(indices)

#         print(type(indices))
