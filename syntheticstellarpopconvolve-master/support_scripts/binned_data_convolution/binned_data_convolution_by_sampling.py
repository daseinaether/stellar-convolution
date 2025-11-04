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
    sample_around_bin_center,
    temp_dir,
)


def post_convolution_function(
    config,
    sfr_dict,
    data_dict,
    time_bin_info_dict,
    convolution_results,
    convolution_instruction,
):
    """ """

    delay_time_values = convolution_results["delay_time"]
    sampled_delay_time_values = sample_around_bin_center(
        time_bin_edges, delay_time_values.value
    )
    sampled_event_time_values = (
        convolution_results["formation_lookback_times"] - delay_time_values
    )

    print(
        'convolution_results["formation_lookback_times"]',
        convolution_results["formation_lookback_times"],
    )
    print('convolution_results["delay_time"]', convolution_results["delay_time"])
    print("sampled_delay_time_values", sampled_delay_time_values)
    print("sampled_event_time_values", sampled_event_time_values)
    # check which ones are negative

    #########

    log10Teff_bin_edges = np.array([1.5, 2.5, 3.5])
    sampled_log10Teff_values = sample_around_bin_center(
        log10Teff_bin_edges, convolution_results["log10Teff"]
    )
    print('convolution_results["log10Teff"]', convolution_results["log10Teff"])
    print("sampled_log10Teff_values", sampled_log10Teff_values)

    return convolution_results


# records = [
#     {"log10Teff": 3, "time": 2 + 0.25, "value": 10, "probability": 1},
#     {"log10Teff": 3, "time": 2 + 1.25, "probability": 2, "value": 20},
#     {"log10Teff": 3.5, "time": 2 + 2.25, "probability": 3, "value": 30},
#     {"log10Teff": 3, "time": 2 + 0.25, "value": 11, "probability": 1.1},
#     {"log10Teff": 3.5, "time": 2 + 3.25, "probability": 4, "value": 40},
#     {"log10Teff": 3, "time": 2 + 0.25, "value": 10, "probability": 1},
#     {"log10Teff": 3.5, "time": 2 + 1.25, "probability": 2, "value": 20},
#     {"log10Teff": 4, "time": 2 + 2.25, "probability": 3, "value": 30},
#     {"log10Teff": 3.5, "time": 2 + 0.25, "value": 11, "probability": 1.1},
#     {"log10Teff": 4, "time": 2 + 3.25, "probability": 4, "value": 40},
# ]

records = [
    {"time": 0.25, "log10Teff": 3.0, "probability": 10},
    {"time": 1.25, "log10Teff": 0, "probability": 0},
]


example_dataframe = pd.DataFrame.from_records(records)

sorted_unique_time_centers = np.sort(example_dataframe["time"].unique())
time_bin_edges = calculate_bin_edges(sorted_unique_time_centers)

#
TMP_DIR = temp_dir("code", "convolve_stochastically", clean_path=True)

# create file
input_hdf5_filename = os.path.join(TMP_DIR, "input_hdf5.h5")
output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5.h5")
input_hdf5_file = h5py.File(input_hdf5_filename, "w")

# Create groups main
input_hdf5_file.create_group("input_data")
input_hdf5_file.create_group("config")

# close
input_hdf5_file.close()


# store the data frame in the hdf5file
example_dataframe.to_hdf(input_hdf5_filename, key="input_data/binned_example")

#
convolution_config = copy.copy(default_convolution_config)
convolution_config["input_filename"] = input_hdf5_filename
convolution_config["output_filename"] = output_hdf5_filename
convolution_config["tmp_dir"] = TMP_DIR
convolution_config["multiprocessing"] = False
convolution_config["logger"].setLevel(logging.CRITICAL)

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
            # required
            # "normalized_yield": {"column_name": "probability", "unit": u.Msun / u.Msun},
            "normalized_yield": {"column_name": "probability", "unit": 1 / u.Msun},
            # "normalized_yield": "probability",
            "delay_time": {"column_name": "time", "unit": u.yr},
            "log10Teff": "log10Teff",
        },
        "post_convolution_function": post_convolution_function,
        "assign_formation_lookback_time": False,
    },
]

#
convolution_config["time_type"] = "lookback_time"
convolution_config["convolution_lookback_time_bin_edges"] = np.arange(0, 2, 1) * u.yr
# print(convolution_config["normalized_yield_unit"])

# construct the sfr-dict (NOTE: this uses absolute SFR, not metallicity dependent)
sfr_dict = {}
sfr_dict["lookback_time_bin_edges"] = np.arange(0, 2, 1) * u.yr

sfr_dict["starformation_rate_array"] = (
    np.ones(len(sfr_dict["lookback_time_bin_edges"]) - 1) * u.Msun / u.yr
)

# store
convolution_config["SFR_info"] = sfr_dict

#
convolution_config["multiply_by_sfr_time_binsize"] = False
convolution_config["multiply_by_convolution_time_binsize"] = False

# convolve
convolve(config=convolution_config)

print("finished convolution")


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

        # #
        # yield_result = output_hdf5_file[f"{main_group}/{formation_time_bin_key}/yield"][
        #     ()
        # ]

        # print(yield_result)
