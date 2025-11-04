import copy
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
    extract_unit_dict,
    generate_boilerplate_outputfile,
    temp_dir,
)

TMP_DIR = temp_dir("examples", "minimal_working_example", clean_path=True)

# Create instance of output
output_hdf5_filename = os.path.join(TMP_DIR, "output_example.h5")
generate_boilerplate_outputfile(output_hdf5_filename)

# SET UP DATA
example_data = {
    "delay_time": np.array([0, 1, 2, 3]),
    "value": np.array([3, 2, 1, 0]),
    "probability": np.array([1, 2, 3, 4]),
}
example_df = pd.DataFrame.from_records(example_data)
example_df.to_hdf(output_hdf5_filename, key="input_data/example")

# Set up global configuration
convolution_config = copy.copy(default_convolution_config)
convolution_config["output_filename"] = output_hdf5_filename

# Set up SFR
convolution_config["SFR_info"] = {
    "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
    "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
}

# set up convolution bin edges
convolution_config["convolution_lookback_time_bin_edges"] = np.array([0, 1]) * u.yr

# Set up the convolution instructions
convolution_config["convolution_instructions"] = [
    {
        **default_convolution_instruction,
        "input_data_name": "example",
        "output_data_name": "example",
        "data_column_dict": {
            "delay_time": "delay_time",
            "normalized_yield": {"column_name": "probability", "unit": 1 / u.Msun},
        },
    }
]

# run convolution
convolve(convolution_config)

# read out results
with h5py.File(convolution_config["output_filename"], "r") as output_hdf5file:
    groupname = "output_data/example/example/convolution_results/0.5 yr/"

    data = output_hdf5file[groupname + "/yield"][()]
    unit_dict = extract_unit_dict(output_hdf5file, groupname)

    print(data)
