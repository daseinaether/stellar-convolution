# """
# Testcases for convolve_ensembles file
# """

# import copy
# import json
# import os
# import unittest

# import astropy.units as u
# import h5py
# import numpy as np
# import pkg_resources

# from syntheticstellarpopconvolve import default_convolution_config
# from syntheticstellarpopconvolve.check_and_update_convolution_config import (
#     check_and_update_convolution_config,
# )
# from syntheticstellarpopconvolve.convolve_ensembles import (
#     convolve_ensemble_by_integration,
#     ensemble_handle_marginalisation,
#     ensemble_handle_SFR_multiplication,
#     extract_ensemble_data,
#     handle_binsize_multiplication_factor,
# )
# from syntheticstellarpopconvolve.general_functions import temp_dir
# from syntheticstellarpopconvolve.prepare_output_file import prepare_output_file

# TMP_DIR = temp_dir(
#     "tests", "tests_convolution", "tests_convolve_ensembles", clean_path=True
# )


# class test_convolve_ensemble_by_integration(unittest.TestCase):
#     def setUp(self):
#         #
#         input_hdf5_filename = os.path.join(TMP_DIR, "input_hdf5_sfr_only.h5")
#         output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5_sfr_only.h5")

#         ##############
#         # SET UP DATA
#         self.dummy_ensemble = {
#             "metallicity": {
#                 "0": {
#                     "delay_time": {
#                         "0": {"a": {"1": 1}, "b": {"1": 1}},
#                         "1": {"a": {"1": 2}, "b": {"1": 2}},
#                         "2": {"a": {"1": 3}, "b": {"1": 3}},
#                         "3": {"a": {"1": 4}, "b": {"1": 4}},
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
#                 "input_data/dummy", data=json.dumps(self.dummy_ensemble)
#             )

#         #
#         self.convolution_config = copy.copy(default_convolution_config)

#         # Set up SFR
#         self.convolution_config["SFR_info"] = {
#             "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
#             "starformation_rate_array": np.array([1, 1, 1, 1, 1])
#             * u.Msun
#             / u.yr
#             / u.Gpc**3,
#         }

#         # set up convolution binsqq
#         self.convolution_config["convolution_lookback_time_bin_edges"] = (
#             np.array([0, 1, 2, 3, 4]) * u.yr
#         )

#         # lookback time convolution only
#         self.convolution_config["time_type"] = "lookback_time"

#         #
#         self.convolution_config["input_filename"] = input_hdf5_filename
#         self.convolution_config["output_filename"] = output_hdf5_filename

#         self.convolution_config["redshift_interpolator_data_output_filename"] = (
#             os.path.join(TMP_DIR, "interpolator_dict.p")
#         )

#         #
#         self.convolution_config["convolution_instructions"] = [
#             {
#                 "input_data_name": "dummy",
#                 "output_data_name": "dummy",
#                 "convolution_type": "integrate",
#                 "data_layer_dict": {
#                     "delay_time": 3,
#                 },
#             },
#         ]

#         #
#         self.convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")

#         #
#         check_and_update_convolution_config(self.convolution_config)

#         #
#         prepare_output_file(config=self.convolution_config)

#     def test_normal(self):
#         _, data_dict, _ = extract_ensemble_data(
#             config=self.convolution_config,
#             convolution_instruction=self.convolution_config["convolution_instructions"][
#                 0
#             ],
#         )

#         sfr_dict = self.convolution_config["SFR_info"]
#         import logging

#         self.convolution_config["logger"].setLevel(logging.CRITICAL)
#         time_bin_info_dict = {
#             "bin_number": 0,
#             "bin_center": 0.5 * u.yr,
#             "bin_edge_lower": 0,
#             "bin_size": 1 * u.yr,
#             "bin_type": "convolution time",
#             "time_type": self.convolution_config["time_type"],
#         }

#         #
#         result_dict = convolve_ensemble_by_integration(
#             sfr_dict=sfr_dict,
#             time_bin_info_dict=time_bin_info_dict,
#             config=self.convolution_config,
#             convolution_instruction=self.convolution_config["convolution_instructions"][
#                 0
#             ],
#             data_dict=data_dict,
#         )

#         #
#         expected_stripped_ensemble = {
#             "metallicity": {
#                 "0": {
#                     "delay_time": {
#                         "0": {"a": {"1": 0}, "b": {"1": 0}},
#                         "1": {"a": {"1": 0}, "b": {"1": 0}},
#                         "2": {"a": {"1": 0}, "b": {"1": 0}},
#                         "3": {"a": {"1": 0}, "b": {"1": 0}},
#                     }
#                 }
#             }
#         }

#         #
#         self.assertTrue("convolution_results" in result_dict)
#         np.testing.assert_array_equal(
#             result_dict["convolution_results"]["yield"],
#             np.array([1, 1, 2, 2, 3, 3, 4, 4]) * (1.0 / u.yr / u.Gpc**3),
#         )

#         #
#         self.assertTrue("stripped_ensemble" in result_dict["convolution_results"])
#         self.assertTrue(
#             result_dict["convolution_results"]["stripped_ensemble"]
#             == {"dummy": expected_stripped_ensemble}
#         )
#         import logging

#         self.convolution_config["logger"].setLevel(logging.CRITICAL)


# class test_ensemble_handle_SFR_multiplication(unittest.TestCase):
#     def setUp(self):
#         #
#         input_hdf5_filename = os.path.join(TMP_DIR, "input_hdf5_sfr_only.h5")
#         output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5_sfr_only.h5")

#         ##############
#         # SET UP DATA
#         self.dummy_ensemble = {
#             "metallicity": {
#                 "0": {
#                     "delay_time": {
#                         "0": {"a": {"1": 1}, "b": {"1": 1}},
#                         "1": {"a": {"1": 2}, "b": {"1": 2}},
#                         "2": {"a": {"1": 3}, "b": {"1": 3}},
#                         "3": {"a": {"1": 4}, "b": {"1": 4}},
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
#                 "input_data/dummy", data=json.dumps(self.dummy_ensemble)
#             )

#         #
#         self.convolution_config = copy.copy(default_convolution_config)

#         # Set up SFR
#         self.convolution_config["SFR_info"] = {
#             "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.Gyr,
#             "starformation_rate_array": np.array([2, 1, 1, 1, 1])
#             * u.Msun
#             / u.yr
#             / u.Gpc**3,
#         }

#         # set up convolution bins
#         self.convolution_config["convolution_lookback_time_bin_edges"] = (
#             np.array([0, 1, 2, 3, 4]) * u.Gyr
#         )

#         # lookback time convolution only
#         self.convolution_config["time_type"] = "lookback_time"

#         #
#         self.convolution_config["input_filename"] = input_hdf5_filename
#         self.convolution_config["output_filename"] = output_hdf5_filename

#         self.convolution_config["redshift_interpolator_data_output_filename"] = (
#             os.path.join(TMP_DIR, "interpolator_dict.p")
#         )

#         #
#         self.convolution_config["convolution_instructions"] = [
#             {
#                 "input_data_name": "dummy",
#                 "output_data_name": "dummy",
#                 "convolution_type": "integrate",
#                 "data_layer_dict": {
#                     "delay_time": 3,
#                 },
#             },
#         ]

#         #
#         self.convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")

#         #
#         check_and_update_convolution_config(config=self.convolution_config)

#         #
#         prepare_output_file(config=self.convolution_config)

#     def test_ensemble_handle_SFR_multiplication_normal(self):

#         #
#         ensemble = self.dummy_ensemble["metallicity"]["0"]["delay_time"]["0"]

#         #
#         data_dict = {"delay_time": 0}

#         sfr_dict = self.convolution_config["SFR_info"]

#         time_bin_info_dict = {
#             "bin_number": 0,
#             "bin_center": 0.5 * u.Gyr,
#             "bin_edge_lower": 0,
#             "bin_size": 1 * u.Gyr,
#             "bin_type": "convolution time",
#             "time_type": self.convolution_config["time_type"],
#         }

#         ensemble = ensemble_handle_SFR_multiplication(
#             sfr_dict=sfr_dict,
#             time_bin_info_dict=time_bin_info_dict,
#             config=self.convolution_config,
#             convolution_instruction=self.convolution_config["convolution_instructions"][
#                 0
#             ],
#             ensemble=ensemble,
#             data_dict=data_dict,
#             extra_value_dict=None,
#         )

#         unit = u.Msun / u.yr / u.Gpc**3
#         expected_ensemble = {"a": {"1": 2.0 * unit}, "b": {"1": 2.0 * unit}}

#         self.assertTrue(ensemble == expected_ensemble)

#     def test_ensemble_handle_SFR_multiplication_extra_value(self):

#         #
#         ensemble = self.dummy_ensemble["metallicity"]["0"]["delay_time"]["0"]

#         #
#         data_dict = {"delay_time": 0}
#         extra_value_dict = {"time_bin": 3, "metallicity_bin": 4}

#         sfr_dict = self.convolution_config["SFR_info"]

#         time_bin_info_dict = {
#             "bin_number": 0,
#             "bin_center": 0.5 * u.Gyr,
#             "bin_edge_lower": 0,
#             "bin_size": 1 * u.Gyr,
#             "bin_type": "convolution time",
#             "time_type": self.convolution_config["time_type"],
#         }

#         ensemble = ensemble_handle_SFR_multiplication(
#             sfr_dict=sfr_dict,
#             time_bin_info_dict=time_bin_info_dict,
#             config=self.convolution_config,
#             convolution_instruction=self.convolution_config["convolution_instructions"][
#                 0
#             ],
#             ensemble=ensemble,
#             data_dict=data_dict,
#             extra_value_dict=extra_value_dict,
#         )

#         unit = u.Msun / u.yr / u.Gpc**3
#         expected_ensemble = {"a": {"1": 24.0 * unit}, "b": {"1": 24.0 * unit}}

#         self.assertTrue(ensemble == expected_ensemble)


# class test_handle_binsize_multiplication_factor(unittest.TestCase):
#     def test_handle_binsize_multiplication_factor_no_binsize_multiplication(self):
#         convolution_config = copy.copy(default_convolution_config)

#         binsizes, extra_value_dict = handle_binsize_multiplication_factor(
#             config=convolution_config,
#             ensemble={"0.1": 2, "0.2": 2, "0.4": 3},
#             data_layer_dict_entry={"binsizes": [1, 2]},
#             key="0.1",
#             key_i=0,
#             binsizes=None,
#             extra_value_dict={},
#             name="delay_time",
#         )

#         self.assertFalse(extra_value_dict)

#     def test_handle_binsize_multiplication_factor_binsizes_passed(self):
#         convolution_config = copy.copy(default_convolution_config)

#         binsizes, extra_value_dict = handle_binsize_multiplication_factor(
#             config=convolution_config,
#             ensemble={"0.1": 2, "0.2": 2, "0.4": 3},
#             data_layer_dict_entry={"multiply_by_binsize": True},
#             key="0.1",
#             key_i=0,
#             binsizes=[0.1, 0.15, 0.2],
#             extra_value_dict={},
#             name="delay_time",
#         )
#         self.assertTrue(extra_value_dict == {"delay_time_binsize": 0.1})

#     def test_handle_binsize_multiplication_factor_binsizes_calculated(self):
#         convolution_config = copy.copy(default_convolution_config)

#         binsizes, extra_value_dict = handle_binsize_multiplication_factor(
#             config=convolution_config,
#             ensemble={"0.1": 2, "0.2": 2, "0.4": 3},
#             data_layer_dict_entry={"multiply_by_binsize": True},
#             key="0.2",
#             key_i=1,
#             binsizes=None,
#             extra_value_dict={},
#             name="delay_time",
#         )

#         self.assertAlmostEqual(extra_value_dict["delay_time_binsize"], 0.15, 6)


# if __name__ == "__main__":
#     unittest.main()
