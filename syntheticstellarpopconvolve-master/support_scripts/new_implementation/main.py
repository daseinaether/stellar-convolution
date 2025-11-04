"""
Function to handle binned data convolution
"""

import copy
import json
import logging
import os

import astropy.units as u
import h5py
import numpy as np
import pandas as pd

from syntheticstellarpopconvolve import (
    convolve,
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.general_functions import (
    calculate_bin_edges,
    generate_boilerplate_outputfile,
    temp_dir,
)

#
TMP_DIR = temp_dir("code", "convolve_stochastically", clean_path=True)


def post_convolution_function(
    config,
    sfr_dict,
    time_bin_info_dict,
    data_dict,
    convolution_results,
    convolution_instruction,
):
    """
    Post-convolution function to handle integrating the systems forward in
    time and finding those that end up in the LISA waveband.

    using local_indices to select everything and using Alexey's distance
    sampler to handle sampling the distances
    """

    print("=========")
    print(time_bin_info_dict)
    print(data_dict)
    print(convolution_results)

    # # unpack data
    # system_indices = convolution_results["indices"]

    # print(system_indices)

    return convolution_results


##############
# create file
output_hdf5_filename = os.path.join(TMP_DIR, "input_hdf5.h5")
generate_boilerplate_outputfile(output_hdf5_filename)

##############
# Create input data
# records = [
#     {"time": 0.25, "value": 10, "probability": 1},
#     {"time": 1.25, "probability": 2, "value": 20},
#     {"time": 2.25, "probability": 3, "value": 30},
#     # {"time": 0.25, "value": 11, "probability": 1.1},
#     # {"time": 3.25, "probability": 4, "value": 40},
# ]

records = [
    {"time": 0.5, "value": 10, "probability": 1},
    {"time": 1.5, "probability": 2, "value": 20},
    {"time": 2.5, "probability": 3, "value": 30},
]


example_dataframe = pd.DataFrame.from_records(records)

sorted_unique_time_centers = np.sort(example_dataframe["time"].unique())
time_bin_edges = calculate_bin_edges(sorted_unique_time_centers)

# store the data frame in the hdf5file
example_dataframe.to_hdf(output_hdf5_filename, key="input_data/binned_example")

##########
#
convolution_config = copy.copy(default_convolution_config)
convolution_config["output_filename"] = output_hdf5_filename
convolution_config["tmp_dir"] = TMP_DIR
convolution_config["multiprocessing"] = False
convolution_config["logger"].setLevel(logging.INFO)
convolution_config["multiply_by_sfr_time_binsize"] = False
convolution_config["multiply_by_convolution_time_binsize"] = False
convolution_config["time_type"] = "lookback_time"
convolution_config["convolution_lookback_time_bin_edges"] = np.arange(0, 3, 1) * u.yr

###
# convolution instructions
convolution_config["convolution_instructions"] = [
    {
        **default_convolution_instruction,
        "convolution_type": "sample",
        "input_data_name": "binned_example",
        "output_data_name": "binned_example",
        "contains_binned_data": True,
        "delay_time_data_bin_info_dict": {
            "delay_time_data_bin_edges": time_bin_edges * u.yr
        },
        "data_column_dict": {
            "normalized_yield": {"column_name": "probability", "unit": u.yr / u.Msun},
            "delay_time": {"column_name": "time", "unit": u.yr},
            "value": "value",
        },
        "post_convolution_function": post_convolution_function,
    },
]

############
# construct the sfr-dict (NOTE: this uses absolute SFR, not metallicity dependent)
sfr_dict = {}
sfr_dict["lookback_time_bin_edges"] = np.arange(0, 10, 1) * u.yr
sfr_dict["starformation_rate_array"] = (
    np.arange(0, len(sfr_dict["lookback_time_bin_edges"]) - 1) ** 2 * u.Msun / u.yr
)
convolution_config["SFR_info"] = sfr_dict

# convolve
convolve(config=convolution_config)

print("finished convolution")


quit()
# read out content and integrate until today
with h5py.File(convolution_config["output_filename"], "r") as output_hdf5_file:
    #
    main_group = "output_data/binned_example/binned_example/convolution_results/"

    ################
    #

    # loop over the formation-time bins
    formation_time_bin_keys = list(output_hdf5_file[main_group].keys())
    formation_time_bin_keys = sorted(
        formation_time_bin_keys, key=lambda x: float(x.split(" ")[0])
    )
    for formation_time_bin_key in formation_time_bin_keys:

        print("=================================")
        print(f"formation_time_bin_key: {formation_time_bin_key}")
        print("=================================")

        ###########
        # Read out data

        # convert units
        unit_dict = json.loads(
            output_hdf5_file[f"{main_group}/{formation_time_bin_key}"].attrs["units"]
        )
        unit_dict = {key: u.Unit(val) for key, val in unit_dict.items()}
        print(unit_dict)

        #
        yield_result = output_hdf5_file[f"{main_group}/{formation_time_bin_key}/yield"][
            ()
        ]

        print(yield_result)
