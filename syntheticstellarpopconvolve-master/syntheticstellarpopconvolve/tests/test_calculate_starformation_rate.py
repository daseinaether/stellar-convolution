"""
This is the unittest file for the calculate_starformation_rate.py source file

TODO: binned data backward integrate absolute
TODO: binned data backward integrate metallicity-weighted

TODO: forward integrate absolute
TODO: forward integrate metallicity-weighted
"""

import copy
import logging
import os
import unittest

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve import (
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.calculate_starformation_rate import (
    calculate_digitized_sfr_rates_binned_data_for_backward_convolution,
    calculate_digitized_sfr_rates_non_binned_data_for_backward_convolution,
    calculate_origin_time_array,
    general_sfr_digitise_function,
)
from syntheticstellarpopconvolve.check_and_update_sfr_dict import (
    check_and_update_sfr_dict,
)
from syntheticstellarpopconvolve.general_functions import (
    calculate_bin_edges,
    calculate_bincenters,
    temp_dir,
)
from syntheticstellarpopconvolve.prepare_redshift_interpolator import (
    prepare_redshift_interpolator,
)

np.random.seed(0)

TMP_DIR = temp_dir(
    "tests",
    "tests_convolution",
    "tests_calculat_starformation_rate",
    clean_path=True,
)

# class test_calculate_digitized_sfr_rates(unittest.TestCase):
#     def setUp(self):
#         #
#         input_hdf5_filename = os.path.join(TMP_DIR, "input_hdf5_sfr_only.h5")
#         output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5_sfr_only.h5")

#         ##############
#         # SET UP DATA
#         self.dummy_data = {
#             "delay_time": np.array([0, 1, 2, 3]) * u.yr,
#             "probability": np.array([1, 2, 3, 4]),
#         }
#         dummy_df = pd.DataFrame.from_records(self.dummy_data)

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

#         ##############
#         # Store data in pandas
#         dummy_df.to_hdf(input_hdf5_filename, key="input_data/{}".format("dummy"))

#         #
#         self.convolution_config = copy.copy(default_convolution_config)

#         # Set up SFR
#         self.convolution_config["SFR_info"] = {
#             "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * 1e9 * u.yr,
#             "starformation_rate_array": np.array([1, 2, 3, 4, 5])
#             * u.Msun
#             / u.yr
#             / u.Gpc**3,
#         }

#         # set up convolution bins
#         self.convolution_config["convolution_lookback_time_bin_edges"] = (
#             np.array([0, 1, 2, 3, 4]) * 1e9 * u.yr
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
#                 **default_convolution_instruction,
#                 "input_data_name": "dummy",
#                 "output_data_name": "dummy",
#                 "convolution_type": "integrate",
#                 "data_column_dict": {
#                     "delay_time": "delay_time",
#                     "normalized_yield": "probability",
#                 },
#             },
#         ]

#         #
#         self.convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")

#         #
#         check_and_update_convolution_config(self.convolution_config)

#         #
#         check_and_prepare_output_file(config=self.convolution_config)

#     def test_calculate_digitized_sfr_rates_sfr_only(self):

#         digitized_sfr_rates = calculate_digitized_sfr_rates(
#             config=self.convolution_config,
#             convolution_time_bin_center=0.5 * 1e9 * u.yr,
#             data_dict={"delay_time": np.array([-1, 1, 2, 3, 100]) * 1e9 * u.yr},
#             sfr_dict=self.convolution_config["SFR_info"],
#             convolution_instruction=self.convolution_config["convolution_instructions"][
#                 0
#             ],
#         )
#         output_unit = u.Msun / u.yr / u.Gpc**3

#         np.testing.assert_array_equal(
#             digitized_sfr_rates, np.array([0.0, 2.0, 3.0, 4.0, 0.0]) * output_unit
#         )

#     # def test_calculate_digitized_sfr_rates_metallicity(self):

#     #     self.convolution_config["SFR_info"]["metallicity_bin_edges"] = (
#     #         self.convolution_config["SFR_info"]["starformation_array"]
#     #         * np.ones(
#     #             (self.convolution_config["SFR_info"]["starformation_array"].shape[0], 3)
#     #         ).T
#     #     )

#     #     # print(self.convolution_config["SFR_info"]["metallicity_bin_edges"])

#     #     #
#     #     sfr_dict = update_sfr_dict(
#     #         sfr_dict=self.convolution_config["SFR_info"], config=self.convolution_config
#     #     )

#     #     digitized_sfr_rates = calculate_digitized_sfr_rates(
#     #         config=self.convolution_config,
#     #         convolution_time_bin_center=0.5 * 1e9 * u.yr,
#     #         data_dict={"delay_time": np.array([-1, 1, 2, 3, 100]) * 1e9},
#     #         sfr_dict=sfr_dict,
#     #     )

#     #     np.testing.assert_array_equal(
#     #         digitized_sfr_rates, np.array([0.0, 2.0, 3.0, 4.0, 0.0])
#     #     )


class test_general_sfr_digitise_function(unittest.TestCase):
    def test_general_sfr_digitise_function_absolute(self):

        ##############
        dummy_data = {
            "delay_time": np.array([0, 1, 2, 3]) * u.yr,
            "value": np.array([3, 2, 1, 0]),
            "probability": np.array([1, 2, 3, 4]),
        }

        convolution_config = copy.copy(default_convolution_config)
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
        }

        # Set up SFR
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
        }
        convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1]) * u.yr
        )
        convolution_config["logger"].setLevel(logging.CRITICAL)
        convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")
        convolution_config["SFR_info"] = check_and_update_sfr_dict(
            sfr_dict=convolution_config["SFR_info"],
            requires_name=False,
            requires_metallicity_info=False,
            time_type=convolution_config["time_type"],
            config=convolution_config,
        )

        starformation_rate_array = general_sfr_digitise_function(
            config=convolution_config,
            sfr_dict=convolution_config["SFR_info"],
            time_values=dummy_data["delay_time"] + 0.5 * u.yr,
        )

        np.testing.assert_array_almost_equal(
            starformation_rate_array.value, np.array([1, 2, 3, 4])
        )
        self.assertEqual(starformation_rate_array.unit, u.Msun / u.yr)

    def test_general_sfr_digitise_function_metallicity_weighted(self):

        ##############
        dummy_data_with_metallicity = {
            "delay_time": np.array([0, 1, 2, 3, 0, 1, 2, 3]) * u.yr,
            "value": np.array([3, 2, 1, 0, 3, 2, 1, 0]),
            "probability": np.array([1, 2, 3, 4, 1, 2, 3, 4]),
            "metallicity": np.array([0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75]),
        }

        convolution_config = copy.copy(default_convolution_config)
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
        }

        # Set up SFR
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
            "metallicity_bin_edges": np.array([0.0, 0.5, 1.0]),
            "metallicity_distribution_array": np.array(
                [[0.4, 0.6, 0.8, 0.9, 1], [1.6, 1.4, 1.2, 1.1, 1]]
            ).T,
        }

        convolution_config["logger"].setLevel(logging.CRITICAL)
        convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")
        convolution_config["SFR_info"] = check_and_update_sfr_dict(
            sfr_dict=convolution_config["SFR_info"],
            requires_name=False,
            requires_metallicity_info=True,
            time_type=convolution_config["time_type"],
            config=convolution_config,
        )

        starformation_rate_array = general_sfr_digitise_function(
            config=convolution_config,
            sfr_dict=convolution_config["SFR_info"],
            time_values=dummy_data_with_metallicity["delay_time"] + 0.5 * u.yr,
            metallicity_values=dummy_data_with_metallicity["metallicity"],
        )

        np.testing.assert_array_almost_equal(
            starformation_rate_array.value,
            np.array(
                [
                    0.2,
                    0.6,
                    1.2,
                    1.8,
                    0.8,
                    1.4,
                    1.8,
                    2.2,
                ]
            ),
        )
        self.assertEqual(starformation_rate_array.unit, u.Msun / u.yr)


class test_calculate_origin_time_array(unittest.TestCase):
    def test_calculate_origin_time_array_lookback(self):
        convolution_config = copy.copy(default_convolution_config)
        origin_time_array = calculate_origin_time_array(
            config=convolution_config,
            data_dict={"delay_time": np.array([1, 2, 3]) * 1e9 * u.yr},
            convolution_time_bin_center=0.5 * 1e9 * u.yr,
        )

        np.testing.assert_array_equal(
            origin_time_array, np.array([1.5, 2.5, 3.5]) * 1e9 * u.yr
        )

    def test_calculate_origin_time_array_redshift(self):
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["redshift_interpolator_data_output_filename"] = os.path.join(
            TMP_DIR, "interpolator_dict.p"
        )
        convolution_config["time_type"] = "redshift"
        convolution_config = prepare_redshift_interpolator(convolution_config)

        origin_time_array = calculate_origin_time_array(
            config=convolution_config,
            data_dict={"delay_time": np.array([1, 2, 3]) * 1e9 * u.yr},
            convolution_time_bin_center=0.5,  # redshift
        )
        # output_unit = u.Msun/u.yr/u.Gpc**3

        np.testing.assert_array_almost_equal(
            origin_time_array,
            np.array([0.6501032923316669, 0.8336451543045214, 1.0661079791875108]),
        )


class test_calculate_digitized_sfr_rates_non_binned_data_for_backward_convolution(
    unittest.TestCase
):

    def test_calculate_digitized_sfr_rates_non_binned_data_for_backward_convolution_absolute(
        self,
    ):

        ##############
        dummy_data = {
            "delay_time": np.array([0, 1, 2, 3]) * u.yr,
            "value": np.array([3, 2, 1, 0]),
            "probability": np.array([1, 2, 3, 4]),
        }

        # Set up SFR
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
        }
        convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1]) * u.yr
        )
        convolution_config["logger"].setLevel(logging.CRITICAL)
        convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")
        convolution_config["SFR_info"] = check_and_update_sfr_dict(
            sfr_dict=convolution_config["SFR_info"],
            requires_name=False,
            requires_metallicity_info=False,
            time_type=convolution_config["time_type"],
            config=convolution_config,
        )

        #########
        #
        data_dict = dummy_data
        config = convolution_config
        convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "convolution_direction": "backward",
            "convolution_type": "integrate",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
        }
        sfr_dict = convolution_config["SFR_info"]
        convolution_lookback_time_bin_edges = convolution_config[
            "convolution_lookback_time_bin_edges"
        ]
        convolution_time_bin_center = calculate_bincenters(
            convolution_lookback_time_bin_edges
        )[0]

        #
        starformation_rate_array = (
            calculate_digitized_sfr_rates_non_binned_data_for_backward_convolution(
                config=config,
                convolution_instruction=convolution_instruction,
                convolution_time_bin_center=convolution_time_bin_center,
                data_dict=data_dict,
                sfr_dict=sfr_dict,
            )
        )

        np.testing.assert_array_almost_equal(
            starformation_rate_array.value, np.array([1, 2, 3, 4])
        )
        self.assertEqual(starformation_rate_array.unit, u.Msun / u.yr)

    def test_calculate_digitized_sfr_rates_non_binned_data_for_backward_convolution_metallicity_weighted(
        self,
    ):

        ##############
        dummy_data_with_metallicity = {
            "delay_time": np.array([0, 1, 2, 3, 0, 1, 2, 3]) * u.yr,
            "value": np.array([3, 2, 1, 0, 3, 2, 1, 0]),
            "probability": np.array([1, 2, 3, 4, 1, 2, 3, 4]),
            "metallicity": np.array([0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75]),
        }

        # Set up SFR
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
            "metallicity_bin_edges": np.array([0.0, 0.5, 1.0]),
            "metallicity_distribution_array": np.array(
                [[0.4, 0.6, 0.8, 0.9, 1], [1.6, 1.4, 1.2, 1.1, 1]]
            ).T,
        }
        convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1]) * u.yr
        )
        convolution_config["logger"].setLevel(logging.CRITICAL)
        convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")
        convolution_config["SFR_info"] = check_and_update_sfr_dict(
            sfr_dict=convolution_config["SFR_info"],
            requires_name=False,
            requires_metallicity_info=True,
            time_type=convolution_config["time_type"],
            config=convolution_config,
        )

        #########
        #
        data_dict = dummy_data_with_metallicity
        config = convolution_config
        convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "convolution_direction": "backward",
            "convolution_type": "integrate",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
        }
        sfr_dict = convolution_config["SFR_info"]
        convolution_lookback_time_bin_edges = convolution_config[
            "convolution_lookback_time_bin_edges"
        ]
        convolution_time_bin_center = calculate_bincenters(
            convolution_lookback_time_bin_edges
        )[0]

        #
        starformation_rate_array = (
            calculate_digitized_sfr_rates_non_binned_data_for_backward_convolution(
                config=config,
                convolution_instruction=convolution_instruction,
                convolution_time_bin_center=convolution_time_bin_center,
                data_dict=data_dict,
                sfr_dict=sfr_dict,
            )
        )

        np.testing.assert_array_almost_equal(
            starformation_rate_array.value,
            np.array(
                [
                    0.2,
                    0.6,
                    1.2,
                    1.8,
                    0.8,
                    1.4,
                    1.8,
                    2.2,
                ]
            ),
        )
        self.assertEqual(starformation_rate_array.unit, u.Msun / u.yr)


class test_calculate_digitized_sfr_rates_binned_data_for_backward_convolution(
    unittest.TestCase
):

    def test_calculate_digitized_sfr_rates_binned_data_for_backward_convolution_absolute(
        self,
    ):

        ############################
        #

        ##############
        dummy_data = {
            "delay_time": np.array([0, 1, 2, 3]) * u.yr,
            "value": np.array([3, 2, 1, 0]),
            "probability": np.array([1, 2, 3, 4]),
        }

        # Set up SFR
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
        }
        convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1]) * u.yr
        )
        convolution_config["logger"].setLevel(logging.CRITICAL)
        convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")
        convolution_config["SFR_info"] = check_and_update_sfr_dict(
            sfr_dict=convolution_config["SFR_info"],
            requires_name=False,
            requires_metallicity_info=False,
            time_type=convolution_config["time_type"],
            config=convolution_config,
        )

        sorted_unique_time_centers = np.sort(np.unique(dummy_data["delay_time"]))
        time_bin_edges = calculate_bin_edges(sorted_unique_time_centers)

        #########
        #
        data_dict = dummy_data
        config = convolution_config
        convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "contains_binned_data": True,
            "convolution_direction": "backward",
            "convolution_type": "integrate",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
            "delay_time_data_bin_info_dict": {
                "delay_time_data_bin_edges": time_bin_edges
            },
        }
        sfr_dict = convolution_config["SFR_info"]
        convolution_lookback_time_bin_edges = convolution_config[
            "convolution_lookback_time_bin_edges"
        ]
        convolution_time_bin_center = calculate_bincenters(
            convolution_lookback_time_bin_edges
        )[0]

        dummy_data["delay_time_data_bin_index"] = (
            np.digitize(
                dummy_data["delay_time"].to(u.yr),
                convolution_instruction["delay_time_data_bin_info_dict"][
                    "delay_time_data_bin_edges"
                ].to(u.yr),
            )
            - 1
        )

        #
        starformation_rate_array = (
            calculate_digitized_sfr_rates_binned_data_for_backward_convolution(
                config=config,
                convolution_instruction=convolution_instruction,
                convolution_time_bin_center=convolution_time_bin_center,
                data_dict=data_dict,
                sfr_dict=sfr_dict,
                delay_time_data_bin_info_dict=convolution_instruction[
                    "delay_time_data_bin_info_dict"
                ],
            )
        )

        np.testing.assert_array_almost_equal(
            starformation_rate_array.value, np.array([1, 2, 3, 4])
        )
        self.assertEqual(starformation_rate_array.unit, u.Msun / u.yr)

    def test_calculate_digitized_sfr_rates_binned_data_for_backward_convolution_absolute_shifted(
        self,
    ):

        ##############
        dummy_data_shifted = {
            "delay_time": np.array([0.5, 1.5, 2.5, 3.5]) * u.yr,
            "value": np.array([3, 2, 1, 0]),
            "probability": np.array([1, 2, 3, 4]),
        }

        # Set up SFR
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
        }
        convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1]) * u.yr
        )
        convolution_config["logger"].setLevel(logging.CRITICAL)
        convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")
        convolution_config["SFR_info"] = check_and_update_sfr_dict(
            sfr_dict=convolution_config["SFR_info"],
            requires_name=False,
            requires_metallicity_info=False,
            time_type=convolution_config["time_type"],
            config=convolution_config,
        )

        sorted_unique_time_centers = np.sort(
            np.unique(dummy_data_shifted["delay_time"])
        )
        time_bin_edges = calculate_bin_edges(sorted_unique_time_centers)

        #########
        # pack up
        data_dict = dummy_data_shifted
        config = convolution_config
        convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "contains_binned_data": True,
            "convolution_direction": "backward",
            "convolution_type": "integrate",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
            "delay_time_data_bin_info_dict": {
                "delay_time_data_bin_edges": time_bin_edges
            },
        }
        sfr_dict = convolution_config["SFR_info"]
        convolution_lookback_time_bin_edges = convolution_config[
            "convolution_lookback_time_bin_edges"
        ]
        convolution_time_bin_center = calculate_bincenters(
            convolution_lookback_time_bin_edges
        )[0]

        dummy_data_shifted["delay_time_data_bin_index"] = (
            np.digitize(
                dummy_data_shifted["delay_time"].to(u.yr),
                convolution_instruction["delay_time_data_bin_info_dict"][
                    "delay_time_data_bin_edges"
                ].to(u.yr),
            )
            - 1
        )

        #
        starformation_rate_array = (
            calculate_digitized_sfr_rates_binned_data_for_backward_convolution(
                config=config,
                convolution_instruction=convolution_instruction,
                convolution_time_bin_center=convolution_time_bin_center,
                data_dict=data_dict,
                sfr_dict=sfr_dict,
                delay_time_data_bin_info_dict=convolution_instruction[
                    "delay_time_data_bin_info_dict"
                ],
            )
        )

        np.testing.assert_array_almost_equal(
            starformation_rate_array.value,
            np.array([1.5, 2.5, 3.5, 4.5]),
        )
        self.assertEqual(starformation_rate_array.unit, u.Msun / u.yr)

    def test_calculate_digitized_sfr_rates_binned_data_for_backward_convolution_metallicity_weighted(
        self,
    ):

        ############################
        #

        ##############
        dummy_data = {
            "delay_time": np.array([0, 1, 2, 3, 0, 1, 2, 3]) * u.yr,
            "value": np.array([3, 2, 1, 0, 3, 2, 1, 0]),
            "probability": np.array([1, 2, 3, 4, 1, 2, 3, 4]),
            "metallicity": np.array([0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75]),
        }

        # Set up SFR
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
            "metallicity_bin_edges": np.array([0.0, 0.5, 1.0]),
            "metallicity_distribution_array": np.array(
                [[0.4, 0.6, 0.8, 0.9, 1], [1.6, 1.4, 1.2, 1.1, 1]]
            ).T,
        }
        convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1]) * u.yr
        )
        convolution_config["logger"].setLevel(logging.CRITICAL)
        convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")
        convolution_config["SFR_info"] = check_and_update_sfr_dict(
            sfr_dict=convolution_config["SFR_info"],
            requires_name=False,
            requires_metallicity_info=False,
            time_type=convolution_config["time_type"],
            config=convolution_config,
        )

        sorted_unique_time_centers = np.sort(np.unique(dummy_data["delay_time"]))
        time_bin_edges = calculate_bin_edges(sorted_unique_time_centers)

        #########
        #
        data_dict = dummy_data
        config = convolution_config
        convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "contains_binned_data": True,
            "convolution_direction": "backward",
            "convolution_type": "integrate",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
            "delay_time_data_bin_info_dict": {
                "delay_time_data_bin_edges": time_bin_edges
            },
        }
        sfr_dict = convolution_config["SFR_info"]
        convolution_lookback_time_bin_edges = convolution_config[
            "convolution_lookback_time_bin_edges"
        ]
        convolution_time_bin_center = calculate_bincenters(
            convolution_lookback_time_bin_edges
        )[0]

        dummy_data["delay_time_data_bin_index"] = (
            np.digitize(
                dummy_data["delay_time"].to(u.yr),
                convolution_instruction["delay_time_data_bin_info_dict"][
                    "delay_time_data_bin_edges"
                ].to(u.yr),
            )
            - 1
        )

        #
        starformation_rate_array = (
            calculate_digitized_sfr_rates_binned_data_for_backward_convolution(
                config=config,
                convolution_instruction=convolution_instruction,
                convolution_time_bin_center=convolution_time_bin_center,
                data_dict=data_dict,
                sfr_dict=sfr_dict,
                delay_time_data_bin_info_dict=convolution_instruction[
                    "delay_time_data_bin_info_dict"
                ],
            )
        )

        np.testing.assert_array_almost_equal(
            starformation_rate_array.value,
            np.array(
                [
                    0.2,
                    0.6,
                    1.2,
                    1.8,
                    0.8,
                    1.4,
                    1.8,
                    2.2,
                ]
            ),
        )
        self.assertEqual(starformation_rate_array.unit, u.Msun / u.yr)

    def test_calculate_digitized_sfr_rates_binned_data_for_backward_convolution_metallicity_weighted_shifted(
        self,
    ):

        ############################
        #

        ##############
        dummy_data = {
            "delay_time": np.array([0.5, 1.5, 2.5, 3.5, 0.5, 1.5, 2.5, 3.5]) * u.yr,
            "value": np.array([3, 2, 1, 0, 3, 2, 1, 0]),
            "probability": np.array([1, 2, 3, 4, 1, 2, 3, 4]),
            "metallicity": np.array([0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75]),
        }

        # Set up SFR
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
            "metallicity_bin_edges": np.array([0.0, 0.5, 1.0]),
            "metallicity_distribution_array": np.array(
                [[0.4, 0.6, 0.8, 0.9, 1], [1.6, 1.4, 1.2, 1.1, 1]]
            ).T,
        }
        convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1]) * u.yr
        )
        convolution_config["logger"].setLevel(logging.CRITICAL)
        convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")
        convolution_config["SFR_info"] = check_and_update_sfr_dict(
            sfr_dict=convolution_config["SFR_info"],
            requires_name=False,
            requires_metallicity_info=False,
            time_type=convolution_config["time_type"],
            config=convolution_config,
        )

        sorted_unique_time_centers = np.sort(np.unique(dummy_data["delay_time"]))
        time_bin_edges = calculate_bin_edges(sorted_unique_time_centers)

        #########
        #
        data_dict = dummy_data
        config = convolution_config
        convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "contains_binned_data": True,
            "convolution_direction": "backward",
            "convolution_type": "integrate",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
            "delay_time_data_bin_info_dict": {
                "delay_time_data_bin_edges": time_bin_edges
            },
        }
        sfr_dict = convolution_config["SFR_info"]
        convolution_lookback_time_bin_edges = convolution_config[
            "convolution_lookback_time_bin_edges"
        ]
        convolution_time_bin_center = calculate_bincenters(
            convolution_lookback_time_bin_edges
        )[0]

        dummy_data["delay_time_data_bin_index"] = (
            np.digitize(
                dummy_data["delay_time"].to(u.yr),
                convolution_instruction["delay_time_data_bin_info_dict"][
                    "delay_time_data_bin_edges"
                ].to(u.yr),
            )
            - 1
        )

        #
        starformation_rate_array = (
            calculate_digitized_sfr_rates_binned_data_for_backward_convolution(
                config=config,
                convolution_instruction=convolution_instruction,
                convolution_time_bin_center=convolution_time_bin_center,
                data_dict=data_dict,
                sfr_dict=sfr_dict,
                delay_time_data_bin_info_dict=convolution_instruction[
                    "delay_time_data_bin_info_dict"
                ],
            )
        )

        np.testing.assert_array_almost_equal(
            starformation_rate_array.value,
            np.array(
                [
                    0.4,
                    0.9,
                    1.5,
                    2.15,
                    1.1,
                    1.6,
                    2.0,
                    2.35,
                ]
            ),
        )
        self.assertEqual(starformation_rate_array.unit, u.Msun / u.yr)


if __name__ == "__main__":
    unittest.main()
