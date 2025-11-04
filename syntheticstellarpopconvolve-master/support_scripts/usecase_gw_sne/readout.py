"""
Function to readout the
"""

import os
import shutil

import h5py
import pandas as pd

this_file_dir = os.path.dirname(os.path.abspath(__file__))

input_hdf5_filename = "/home/david/projects/binary_c_root/sspc_convolution/support_scripts/usecase_gw_sne/data/EVENTS_V2.2.2_SEMI_HIGH_RES_SCHNEIDER_WIND_PPISN_NEW_FRYER_DELAYED/dco_convolution_results/convolution_results.h5"
minimized_hdf5_filename = "data/minimized_hdf5.h5"

with h5py.File(input_hdf5_filename, "r") as input_hdf5_file:
    print(input_hdf5_file.keys())
    print(input_hdf5_file["data"].keys())
    print(input_hdf5_file["data/events"].keys())

# load df
df = pd.read_hdf(input_hdf5_filename, key="data/combined_dataframes")
print(df.columns)
# print(len(reread.index))

# slim down
minimized_df = df.query("stellar_type_1==14 & stellar_type_2==14")
print(len(minimized_df.index))

# drop columns
minimized_df["delay_time_in_years"] = (
    minimized_df["merger_time_values_in_years"]
    + minimized_df["formation_time_values_in_years"]
)

drop_cols = [
    "random_seed",
    "uuid",
    "local_index",
    "undergone_CE_with_MS_donor",
    "undergone_CE_with_HG_donor",
    "merger_time_values_in_years",
    "formation_time_values_in_years",
]
existing_cols = list(minimized_df.columns[:])
for drop_col in drop_cols:
    existing_cols.remove(drop_col)

#
minimized_df = minimized_df[existing_cols]
print(minimized_df.columns)

# Randomly sample 7 elements from your dataframe
minimized_df = minimized_df.sample(n=500)
print(len(minimized_df.index))

#
if os.path.isfile(minimized_hdf5_filename):
    os.remove(minimized_hdf5_filename)
minimized_df.to_hdf(minimized_hdf5_filename, "data/combined_dataframes")

# copy to
target_file = os.path.join(
    this_file_dir,
    "../../syntheticstellarpopconvolve/example_data/example_data_usecase_gw.h5",
)
shutil.copy(
    minimized_hdf5_filename,
    os.path.join(
        this_file_dir,
        "../../syntheticstellarpopconvolve/example_data/example_data_usecase_gw.h5",
    ),
)

target_df = pd.read_hdf(target_file, key="data/combined_dataframes")
print(target_df.columns)
