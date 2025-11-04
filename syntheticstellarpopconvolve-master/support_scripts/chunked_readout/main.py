import logging
import os

import astropy.units as u
import h5py
import numpy as np
import pandas as pd
import pkg_resources

from syntheticstellarpopconvolve import (
    convolve,
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.general_functions import (
    generate_boilerplate_outputfile,
    temp_dir,
)

TMP_DIR = temp_dir(
    "tests",
    "tests_convolution",
    "general_tests",
    clean_path=True,
)

###################
# Set up data
BinCodex_events_filename = pkg_resources.resource_filename(
    "syntheticstellarpopconvolve",
    "example_data/example_BinCodex_dwd.h5",
)


##############
# create file
output_hdf5_filename = os.path.join(TMP_DIR, "input_hdf5.h5")
generate_boilerplate_outputfile(output_hdf5_filename)


#
BinCodex_T0_events = pd.read_hdf(
    BinCodex_events_filename,
    "T0",
)

##################
# update T0 output

# get mass normalisation
mass_normalisation_fiducial = 4476544.539875359 * u.Msun

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

convolution_config = {**default_convolution_config}
convolution_config["output_filename"] = output_hdf5_filename

convolution_config["logger"].setLevel(logging.INFO)

# store the data frame in the hdf5file
wd_binaries.to_hdf(convolution_config["output_filename"], key="input_data/dummy")
# print(len(wd_binaries.index))

# with pd.HDFStore(convolution_config["output_filename"], "r") as store:
#     key = "input_data/dummy"
#     print("Available keys:", store.keys())
#     n_rows = store.get_storer(key).nrows
#     storer = store.get_storer(key)
#     print("Format:", storer.format)  # Should be 'table' for .nrows to work
#     print(n_rows)


# Set up SFR
convolution_config["SFR_info"] = {
    "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.Gyr,
    "starformation_rate_array": np.array([2, 1, 1, 1, 1]) * u.Msun / u.yr,
}

# set up convolution bins
convolution_config["convolution_lookback_time_bin_edges"] = np.array([0, 1]) * u.Gyr

# lookback time convolution only
convolution_config["time_type"] = "lookback_time"

#
convolution_config["output_filename"] = output_hdf5_filename

convolution_config["redshift_interpolator_data_output_filename"] = os.path.join(
    TMP_DIR, "interpolator_dict.p"
)

#
convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp_{}".format("test"))

#
convolution_config["convolution_instructions"] = [
    {
        **default_convolution_instruction,
        "input_data_name": "dummy",
        "output_data_name": "dummy",
        "convolution_type": "integrate",
        "chunked_readout": True,
        "chunk_size": 2500,
        "chunk_total": 2,
        "data_column_dict": {
            # required
            "normalized_yield": "normalized_yield",
            "delay_time": {"column_name": "time", "unit": u.Myr},
        },
        "multiply_by_sfr_time_binsize": True,
    },
]

# convolve
convolve(config=convolution_config)

# read out content and integrate until today
with h5py.File(convolution_config["output_filename"], "r") as output_hdf5_file:

    print(output_hdf5_file["output_data/dummy"].keys())
    print(output_hdf5_file["output_data/dummy/0/dummy/"].keys())
    yield_res = output_hdf5_file[
        "output_data/dummy/0/dummy/convolution_results/0.5 Gyr/yield"
    ][()]
    print(yield_res.shape)
