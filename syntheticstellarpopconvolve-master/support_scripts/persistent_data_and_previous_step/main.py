"""

"""

import copy
import json
import os

import astropy.units as u
import h5py
import numpy as np
import pandas as pd

from syntheticstellarpopconvolve import convolve, default_convolution_config
from syntheticstellarpopconvolve.general_functions import temp_dir

records = [{"a": 2, "b": 3, "normalized_yield": 1, "time": 1}]
dummy_df = pd.DataFrame.from_records(records)
TMP_DIR = temp_dir("code", "persistent_data", clean_path=True)


def post_convolution_function(
    config,
    sfr_dict,
    data_dict,
    convolution_results,
    convolution_instruction,
    persistent_data,
    previous_convolution_results,
):
    """
    Post-convolution function to handle integrating the systems forward in time and finding those that end up in the LISA waveband.

    using local_indices to select everything and using Alexey's distance sampler to handle sampling the distances
    """

    print(persistent_data)
    print(previous_convolution_results)

    if "test" not in persistent_data:
        persistent_data["test"] = 0
    persistent_data["test"] += 1

    return convolution_results


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
dummy_df.to_hdf(input_hdf5_filename, key="input_data/events/example")


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
convolution_config["multiprocessing"] = False


###
# convolution instructions
convolution_config["convolution_instructions"] = [
    {
        "convolution_type": "integrate",
        "input_data_name": "example",
        "output_data_name": "example",
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
convolution_config["convolution_lookback_time_bin_edges"] = np.arange(0, 6, 1) * u.Gyr


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
    print(output_hdf5_file["output_data/event/example"].keys())
    print(output_hdf5_file["output_data/event/example/example"].keys())
    print(
        output_hdf5_file["output_data/event/example/example/convolution_results"].keys()
    )
    print(
        output_hdf5_file[
            "output_data/event/example/example/convolution_results/0.5 Gyr"
        ].keys()
    )

#     print(
#         output_hdf5_file[
#             "output_data/event/stochastic_example/stochastic_example/convolution_results/set_1"
#         ].keys()
#     )
#     print(
#         output_hdf5_file[
#             "output_data/event/stochastic_example/stochastic_example/convolution_results/set_2"
#         ].keys()
#     )

#     quit()

#     print(
#         output_hdf5_file[
#             "output_data/event/stochastic_example/stochastic_example/convolution_results"
#         ].keys()
#     )

#     formation_time_bin_keys = list(
#         output_hdf5_file[
#             "output_data/event/stochastic_example/stochastic_example/convolution_results"
#         ].keys()
#     )

#     ################
#     #

#     # loop over the formation-time bins
#     formation_time_bin_keys = sorted(
#         formation_time_bin_keys, key=lambda x: float(x.split(" ")[0])
#     )
#     for formation_time_bin_key in formation_time_bin_keys:

#         # formation_time_bin_key = "3500000000.0 yr"
#         print("=================================")
#         print(f"formation_time_bin_key: {formation_time_bin_key}")
#         print("=================================")

#         ###########
#         # Read out data

#         # convert units
#         unit_dict = json.loads(
#             output_hdf5_file[
#                 f"output_data/event/stochastic_example/stochastic_example/convolution_results/{formation_time_bin_key}"
#             ].attrs["units"]
#         )
#         unit_dict = {key: u.Unit(val) for key, val in unit_dict.items()}
#         print(unit_dict)


# #         indices = output_hdf5_file[
# #             f"output_data/event/stochastic_example/stochastic_example/convolution_results/{formation_time_bin_key}/indices"
# #         ][()]
# #         # print(indices)

# #         print(type(indices))
