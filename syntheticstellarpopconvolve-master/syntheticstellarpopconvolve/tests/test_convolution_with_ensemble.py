# """
# Testcases for convolution_with_ensemble file

# TODO: test things with redshift
# TODO: test things with multiply SFR histories
# """

# import copy
# import json
# import logging
# import os
# import unittest

# import astropy.units as u
# import h5py
# import numpy as np
# import pkg_resources

# from syntheticstellarpopconvolve import convolve, default_convolution_config
# from syntheticstellarpopconvolve.general_functions import temp_dir

# TMP_DIR = temp_dir(
#     "tests", "tests_convolution", "test_convolution_with_ensemble", clean_path=True
# )


# class test_convolution_with_ensemble(unittest.TestCase):
#     """ """

#     def test_convolution_with_events_with_lookback_time(self):

#         #
#         convolution_config = copy.copy(default_convolution_config)
#         convolution_config["logger"].setLevel(logging.CRITICAL)

#         #
#         input_hdf5_filename = os.path.join(TMP_DIR, "input_hdf5_sfr_only.h5")
#         output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5_sfr_only.h5")

#         ##############
#         # SET UP DATA
#         dummy_ensemble = {
#             "metallicity": {
#                 "0": {
#                     "delay_time": {
#                         "0": {"a": {"1": 1}},
#                         "1": {"a": {"1": 2}},
#                         "2": {"a": {"1": 3}},
#                         "3": {"a": {"1": 4}},
#                     }
#                 }
#             }
#         }

#         #############
#         # create input HDF5 file
#         with h5py.File(input_hdf5_filename, "w") as input_hdf5_file:

#             ######################
#             # Create groups
#             input_hdf5_file.create_group("input_data")
#             input_hdf5_file.create_group("input_data/ensemble")
#             input_hdf5_file.create_group("config")

#             ###############
#             # Readout population settings
#             population_settings_filename = pkg_resources.resource_filename(
#                 "syntheticstellarpopconvolve",
#                 "example_data/example_population_settings.json",
#             )

#             with open(population_settings_filename, "r") as f:
#                 population_settings = json.loads(f.read())

#             # Delete some stuff from the settings
#             del population_settings["population_settings"]["bse_options"]["metallicity"]

#             # Write population config to file
#             input_hdf5_file.create_dataset(
#                 "config/population", data=json.dumps(population_settings)
#             )

#             #
#             input_hdf5_file.create_dataset(
#                 "input_data/ensemble/dummy", data=json.dumps(dummy_ensemble)
#             )

#         # Set up SFR
#         convolution_config["SFR_info"] = {
#             "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
#             "starformation_rate_array": np.array([1, 1, 1, 1, 1])
#             * u.Msun
#             / u.yr
#             / u.Gpc**3,
#         }

#         # set up convolution bins
#         convolution_config["convolution_lookback_time_bin_edges"] = (
#             np.array([0, 1, 2, 3, 4]) * u.yr
#         )

#         # lookback time convolution only
#         convolution_config["time_type"] = "lookback_time"

#         #
#         convolution_config["input_filename"] = input_hdf5_filename
#         convolution_config["output_filename"] = output_hdf5_filename

#         convolution_config["redshift_interpolator_data_output_filename"] = os.path.join(
#             TMP_DIR, "interpolator_dict.p"
#         )

#         #
#         convolution_config["convolution_instructions"] = [
#             {
#                 "input_data_name": "dummy",
#                 "convolution_type": "integrate",
#                 "output_data_name": "dummy",
#                 "data_layer_dict": {
#                     # "delay_time": 3,
#                     "delay_time": {
#                         "multiply_by_binsize": True,
#                         "layer_depth": 3,
#                         # "binsizes": [10, 10, 10, 10]
#                     },
#                     # 'metallicity': 1,
#                 },
#             }
#         ]

#         #
#         convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")

#         #
#         convolve(config=convolution_config)

#         #
#         with h5py.File(output_hdf5_filename, "r") as output_hdf5_file:

#             #
#             arr_ = output_hdf5_file[
#                 "output_data/ensemble/dummy/dummy/convolution_results/0.5 yr"
#             ]["yield"]
#             self.assertTrue(np.array_equal(arr_, np.array([1, 2, 3, 4])))

#             #
#             arr_ = output_hdf5_file[
#                 "output_data/ensemble/dummy/dummy/convolution_results/1.5 yr"
#             ]["yield"][()]
#             self.assertTrue(np.array_equal(arr_, np.array([1, 2, 3, 4])))

#             #
#             arr_ = output_hdf5_file[
#                 "output_data/ensemble/dummy/dummy/convolution_results/2.5 yr"
#             ]["yield"][()]
#             self.assertTrue(np.array_equal(arr_, np.array([1, 2, 3, 0])))

#             #
#             arr_ = output_hdf5_file[
#                 "output_data/ensemble/dummy/dummy/convolution_results/3.5 yr"
#             ]["yield"][()]
#             self.assertTrue(np.array_equal(arr_, np.array([1, 2, 0, 0])))


# if __name__ == "__main__":
#     unittest.main()
