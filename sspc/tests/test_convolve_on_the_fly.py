"""
Unit tests for convolution on the fly
"""

import copy
import logging
import os
import unittest

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve import (
    convolve,
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.check_and_prepare_output_file import (
    check_and_prepare_output_file,
)
from syntheticstellarpopconvolve.check_and_update_convolution_config import (
    check_and_update_convolution_config,
)
from syntheticstellarpopconvolve.convolve_on_the_fly import (
    convolve_on_the_fly,
    convolve_on_the_fly_post_convolution_hook_wrapper,
    handle_call_on_the_fly_function,
)
from syntheticstellarpopconvolve.general_functions import (
    create_time_bin_info_dict,
    generate_boilerplate_outputfile,
    temp_dir,
)

TMP_DIR = temp_dir(
    "tests", "tests_convolution", "test_convolution_on_the_fly", clean_path=True
)


def dummy_on_the_fly_function(total_star_formation_in_bin):
    """
    on-the-fly function that should work
    """

    return {}


def wrong_arguments_on_the_fly_function():
    """
    on-the-fly function does not have the correct arguments
    """

    return {}


def wrong_return_type_on_the_fly_function(total_star_formation_in_bin):
    """
    On-the-fly function that returns the wrong type of object (None in this case)
    """


def metallicity_required_not_included_on_the_fly_function(total_star_formation_in_bin):
    """
    On-the-fly function that returns the wrong type of object (None in this case)
    """


def metallicity_required_on_the_fly_function(
    total_star_formation_in_bin, metallicity_distribution
):
    """
    On-the-fly function that returns the wrong type of object (None in this case)
    """

    return {}


class test_convolve_on_the_fly(unittest.TestCase):
    """ """

    def setUp(self):
        # setup output file
        output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5.h5")
        generate_boilerplate_outputfile(output_hdf5_filename)

        #
        self.convolution_config = copy.copy(default_convolution_config)
        self.convolution_config["logger"].setLevel(logging.CRITICAL)

        # Set up SFR
        self.convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1, 1, 1]) * u.Msun / u.yr,
        }

        # set up convolution bins
        self.convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1, 2, 3, 4]) * u.yr
        )

        # lookback time convolution only
        self.convolution_config["time_type"] = "lookback_time"

        #
        self.convolution_config["output_filename"] = output_hdf5_filename

        self.convolution_config["redshift_interpolator_data_output_filename"] = (
            os.path.join(TMP_DIR, "interpolator_dict.p")
        )

        #
        self.convolution_config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "binary_c",
                "output_data_name": "BHBH",
                "convolution_type": "on-the-fly",
                "convolution_direction": "forward",
                "data_column_dict": {
                    "delay_time": "delay_time",
                    "normalized_yield": "probability",
                },
                "on_the_fly_function": dummy_on_the_fly_function,
            },
        ]

        #
        self.convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")

        #
        check_and_update_convolution_config(self.convolution_config)

        #
        check_and_prepare_output_file(config=self.convolution_config)

    def test_convolve_on_the_fly_wrong_arguments_on_the_fly_function(self):
        #
        self.convolution_config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "binary_c",
                "output_data_name": "BHBH",
                "convolution_type": "on-the-fly",
                "convolution_direction": "forward",
                "data_column_dict": {
                    "delay_time": "delay_time",
                    "normalized_yield": "probability",
                },
                "on_the_fly_function": wrong_arguments_on_the_fly_function,
            },
        ]

        #
        sfr_dict = self.convolution_config["SFR_info"]

        time_bin_info_dict = {
            "bin_number": 0,
            "bin_center": 0.5 * u.yr,
            "bin_edge_lower": 0,
            "bin_size": 1 * u.yr,
            "bin_type": "starformation time",
            "time_type": self.convolution_config["time_type"],
        }

        with self.assertRaises(ValueError):
            convolve_on_the_fly(
                config=self.convolution_config,
                sfr_dict=sfr_dict,
                convolution_instruction=self.convolution_config[
                    "convolution_instructions"
                ][0],
                time_bin_info_dict=time_bin_info_dict,
            )

    def test_convolve_on_the_fly_normal(self):
        #
        sfr_dict = self.convolution_config["SFR_info"]

        time_bin_info_dict = {
            "bin_number": 0,
            "bin_center": 0.5 * u.yr,
            "bin_edge_lower": 0,
            "bin_size": 1 * u.yr,
            "bin_type": "starformation time",
            "time_type": self.convolution_config["time_type"],
        }

        #
        convolve_on_the_fly(
            config=self.convolution_config,
            sfr_dict=sfr_dict,
            convolution_instruction=self.convolution_config["convolution_instructions"][
                0
            ],
            time_bin_info_dict=time_bin_info_dict,
        )

    def test_convolve_on_the_fly_wrong_return_type_on_the_fly_function(self):

        #
        self.convolution_config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "binary_c",
                "output_data_name": "BHBH",
                "convolution_type": "on-the-fly",
                "convolution_direction": "forward",
                "data_column_dict": {
                    "delay_time": "delay_time",
                    "normalized_yield": "probability",
                },
                "on_the_fly_function": wrong_return_type_on_the_fly_function,
            },
        ]

        #
        sfr_dict = self.convolution_config["SFR_info"]

        time_bin_info_dict = {
            "bin_number": 0,
            "bin_center": 0.5 * u.yr,
            "bin_edge_lower": 0,
            "bin_size": 1 * u.yr,
            "bin_type": "starformation time",
            "time_type": self.convolution_config["time_type"],
        }

        with self.assertRaises(ValueError):
            convolve_on_the_fly(
                config=self.convolution_config,
                sfr_dict=sfr_dict,
                convolution_instruction=self.convolution_config[
                    "convolution_instructions"
                ][0],
                time_bin_info_dict=time_bin_info_dict,
            )

    def test_convolve_on_the_fly_metallicity_required_not_included_on_the_fly_function(
        self,
    ):

        #
        self.convolution_config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "binary_c",
                "output_data_name": "BHBH",
                "convolution_type": "on-the-fly",
                "convolution_direction": "forward",
                "data_column_dict": {
                    "delay_time": "delay_time",
                    "normalized_yield": "probability",
                },
                "on_the_fly_function": metallicity_required_not_included_on_the_fly_function,
            },
        ]

        # Set up SFR
        sfr_dict = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1]) * u.Msun / u.yr,
            "metallicity_bin_edges": np.array([0.01, 0.1, 0.2, 0.3]),
            "metallicity_distribution_array": np.array(
                [[1, 2, 3], [4, 5, 6], [4, 5, 6]]
            ),
        }

        #
        self.convolution_config["SFR_info"] = sfr_dict

        #
        check_and_update_convolution_config(self.convolution_config)

        #
        sfr_dict = self.convolution_config["SFR_info"]

        time_bin_info_dict = {
            "bin_number": 0,
            "bin_center": 0.5 * u.yr,
            "bin_edge_lower": 0,
            "bin_size": 1 * u.yr,
            "bin_type": "starformation time",
            "time_type": self.convolution_config["time_type"],
        }

        with self.assertRaises(ValueError):
            convolve_on_the_fly(
                config=self.convolution_config,
                sfr_dict=sfr_dict,
                convolution_instruction=self.convolution_config[
                    "convolution_instructions"
                ][0],
                time_bin_info_dict=time_bin_info_dict,
            )

    def test_convolve_on_the_fly_metallicity_required_on_the_fly_function(
        self,
    ):

        #
        self.convolution_config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "binary_c",
                "output_data_name": "BHBH",
                "convolution_type": "on-the-fly",
                "convolution_direction": "forward",
                "data_column_dict": {
                    "delay_time": "delay_time",
                    "normalized_yield": "probability",
                },
                "on_the_fly_function": metallicity_required_on_the_fly_function,
            },
        ]

        # Set up SFR
        sfr_dict = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1]) * u.Msun / u.yr,
            "metallicity_bin_edges": np.array([0.01, 0.1, 0.2, 0.3]),
            "metallicity_distribution_array": np.array(
                [[1, 2, 3], [4, 5, 6], [4, 5, 6]]
            ),
        }

        #
        self.convolution_config["SFR_info"] = sfr_dict

        #
        check_and_update_convolution_config(self.convolution_config)

        #
        sfr_dict = self.convolution_config["SFR_info"]

        time_bin_info_dict = {
            "bin_number": 0,
            "bin_center": 0.5 * u.yr,
            "bin_edge_lower": 0,
            "bin_size": 1 * u.yr,
            "bin_type": "starformation time",
            "time_type": self.convolution_config["time_type"],
        }

        convolve_on_the_fly(
            config=self.convolution_config,
            sfr_dict=sfr_dict,
            convolution_instruction=self.convolution_config["convolution_instructions"][
                0
            ],
            time_bin_info_dict=time_bin_info_dict,
        )

    def test_convolve_on_the_fly_convolve(self):

        #
        self.convolution_config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "binary_c",
                "output_data_name": "BHBH",
                "convolution_type": "on-the-fly",
                "convolution_direction": "forward",
                "data_column_dict": {
                    "delay_time": "delay_time",
                    "normalized_yield": "probability",
                },
                "on_the_fly_function": dummy_on_the_fly_function,
            },
        ]

        # Set up SFR
        sfr_dict = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1]) * u.Msun / u.yr,
        }

        #
        self.convolution_config["SFR_info"] = sfr_dict

        convolve(config=self.convolution_config)


def post_convolution_function_no_change(convolution_results):
    """ """

    return convolution_results


class test_convolve_on_the_fly_post_convolution_hook_wrapper(unittest.TestCase):
    def test_convolve_on_the_fly_post_convolution_hook_wrapper(self):

        convolution_results = {"sampled_indices": np.array([1, 2])}

        convolution_instruction = {
            **default_convolution_instruction,
            "convolution_type": "sampling",
            "post_convolution_function": post_convolution_function_no_change,
        }

        time_bin_info_dict = create_time_bin_info_dict(
            config={**default_convolution_config},
            convolution_instruction=convolution_instruction,
            bin_number=1,
            bin_edge_lower=0.5,
            bin_center=0.75,
            bin_size=1,
            bin_type="lookback",
        )

        # any change is allowed
        convolution_results = convolve_on_the_fly_post_convolution_hook_wrapper(
            config={**default_convolution_config},
            convolution_instruction={
                **convolution_instruction,
                "post_convolution_function": post_convolution_function_no_change,
            },
            convolution_results=convolution_results,
            sfr_dict={"lookback_time_edges": np.array([1, 2])},
            time_bin_info_dict=time_bin_info_dict,
        )


def on_the_fly_function():
    pass


def on_the_fly_function_missing_total_star_formation_in_bin():
    pass


def on_the_fly_function_missing_metallicity_distribution(total_star_formation_in_bin):
    pass


def on_the_fly_function_missing_bad_return_type(
    total_star_formation_in_bin, metallicity_distribution
):
    pass


def on_the_fly_function_missing_(total_star_formation_in_bin, metallicity_distribution):
    return {}


class test_handle_call_on_the_fly_function(unittest.TestCase):

    def test_handle_call_on_the_fly_function_no_mass_unit(self):

        convolution_instruction = {
            **default_convolution_instruction,
            "convolution_type": "on_the_fly",
            "on_the_fly_function": on_the_fly_function,
        }

        time_bin_info_dict = create_time_bin_info_dict(
            config={**default_convolution_config},
            convolution_instruction=convolution_instruction,
            bin_number=1,
            bin_edge_lower=0.5 * u.yr,
            bin_center=0.75 * u.yr,
            bin_size=1,
            bin_type="lookback",
        )

        # Set up SFR
        sfr_dict = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1]) * u.Msun / u.yr,
        }

        # will raise issue because unit is not correct of sfr * binsize
        with self.assertRaises(ValueError):
            handle_call_on_the_fly_function(
                config={**default_convolution_config},
                time_bin_info_dict=time_bin_info_dict,
                sfr_dict=sfr_dict,
                convolution_instruction=convolution_instruction,
            )

    def test_handle_call_on_the_fly_function_no_total_star_formation_in_bin(self):

        convolution_instruction = {
            **default_convolution_instruction,
            "convolution_type": "on_the_fly",
            "on_the_fly_function": on_the_fly_function_missing_total_star_formation_in_bin,
        }

        time_bin_info_dict = create_time_bin_info_dict(
            config={**default_convolution_config},
            convolution_instruction=convolution_instruction,
            bin_number=1,
            bin_edge_lower=0.5 * u.yr,
            bin_center=0.75 * u.yr,
            bin_size=1 * u.yr,
            bin_type="lookback",
        )

        # Set up SFR
        sfr_dict = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1]) * u.Msun / u.yr,
        }

        # will raise issue because total_star_formation_in_bin should be an argument in function
        with self.assertRaises(ValueError):

            handle_call_on_the_fly_function(
                config={**default_convolution_config},
                time_bin_info_dict=time_bin_info_dict,
                sfr_dict=sfr_dict,
                convolution_instruction=convolution_instruction,
            )

    def test_handle_call_on_the_fly_function_no_metallicity_distribution(self):

        convolution_instruction = {
            **default_convolution_instruction,
            "convolution_type": "on_the_fly",
            "on_the_fly_function": on_the_fly_function_missing_metallicity_distribution,
        }

        time_bin_info_dict = create_time_bin_info_dict(
            config={**default_convolution_config},
            convolution_instruction=convolution_instruction,
            bin_number=1,
            bin_edge_lower=0.5 * u.yr,
            bin_center=0.75 * u.yr,
            bin_size=1 * u.yr,
            bin_type="lookback",
        )

        # Set up SFR
        sfr_dict = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1]) * u.Msun / u.yr,
            "metallicity_distribution_array": np.array([[1, 1]]),
        }

        # will raise issue because 'metallicity_distribution' should be an argument in function
        with self.assertRaises(ValueError):
            handle_call_on_the_fly_function(
                config={**default_convolution_config},
                time_bin_info_dict=time_bin_info_dict,
                sfr_dict=sfr_dict,
                convolution_instruction=convolution_instruction,
            )

    def test_handle_call_on_the_fly_function_bad_return_type(self):

        convolution_instruction = {
            **default_convolution_instruction,
            "convolution_type": "on_the_fly",
            "on_the_fly_function": on_the_fly_function_missing_bad_return_type,
        }

        time_bin_info_dict = create_time_bin_info_dict(
            config={**default_convolution_config},
            convolution_instruction=convolution_instruction,
            bin_number=1,
            bin_edge_lower=0.5 * u.yr,
            bin_center=0.75 * u.yr,
            bin_size=1 * u.yr,
            bin_type="lookback",
        )

        # Set up SFR
        sfr_dict = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1]) * u.Msun / u.yr,
            "metallicity_distribution_array": np.array([[1, 1]]),
        }

        # will raise issue because returned object is not a dictionary
        with self.assertRaises(ValueError):
            handle_call_on_the_fly_function(
                config={**default_convolution_config},
                time_bin_info_dict=time_bin_info_dict,
                sfr_dict=sfr_dict,
                convolution_instruction=convolution_instruction,
            )

    def test_handle_call_on_the_fly_function_(self):

        convolution_instruction = {
            **default_convolution_instruction,
            "convolution_type": "on_the_fly",
            "on_the_fly_function": on_the_fly_function_missing_,
        }

        time_bin_info_dict = create_time_bin_info_dict(
            config={**default_convolution_config},
            convolution_instruction=convolution_instruction,
            bin_number=1,
            bin_edge_lower=0.5 * u.yr,
            bin_center=0.75 * u.yr,
            bin_size=1 * u.yr,
            bin_type="lookback",
        )

        # Set up SFR
        sfr_dict = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1]) * u.Msun / u.yr,
            "metallicity_distribution_array": np.array([[1, 1]]),
        }

        # will raise issue because returned object is not a dictionary
        handle_call_on_the_fly_function(
            config={**default_convolution_config},
            time_bin_info_dict=time_bin_info_dict,
            sfr_dict=sfr_dict,
            convolution_instruction=convolution_instruction,
        )

    # ######
    # # Get quantities
    # bin_number = time_bin_info_dict["bin_number"]
    # bin_size = time_bin_info_dict["bin_size"]
    # bin_lower_edge = time_bin_info_dict["bin_edge_lower"]
    # sfr = sfr_dict["starformation_rate_array"][bin_number]
    # total_star_formation_in_bin = sfr * bin_size

    # # Check if the total star formation is a mass-type value
    # if not is_mass_unit(total_star_formation_in_bin):
    #     raise ValueError(
    #         "The total star formation in current bin ({}) is not of a mass-type unit. Something wrong with either the sfr ({}) or the time-bin size ({})".format(
    #             total_star_formation_in_bin, sfr, bin_size
    #         )
    #     )

    # def test_handle_post_convolution_function_no_convolution_results(self):

    #     convolution_instruction = {
    #         **default_convolution_instruction,
    #         "post_convolution_function": post_convolution_function_no_convolution_result_argument,
    #     }

    #     with self.assertRaises(ValueError):
    #         handle_post_convolution_function(
    #             config={**default_convolution_config},
    #             sfr_dict={},
    #             data_dict={},
    #             time_bin_info_dict={},
    #             convolution_instruction=convolution_instruction,
    #             convolution_results={},
    #             name="no_convolution_results",
    #         )

    # def test_handle_post_convolution_function_wrong_convolution_results_type(self):

    #     convolution_instruction = {
    #         **default_convolution_instruction,
    #         "post_convolution_function": post_convolution_function_wrong_convolution_results_type,
    #     }

    #     with self.assertRaises(ValueError):
    #         handle_post_convolution_function(
    #             config={**default_convolution_config},
    #             sfr_dict={},
    #             data_dict={},
    #             time_bin_info_dict={},
    #             convolution_instruction=convolution_instruction,
    #             convolution_results={},
    #             name="wrong_convolution_results_type",
    #         )

    # def test_handle_post_convolution_function_wrong_elements_in_convolution_results_list(
    #     self,
    # ):

    #     convolution_instruction = {
    #         **default_convolution_instruction,
    #         "post_convolution_function": post_convolution_function_wrong_elements_in_convolution_results_list,
    #     }

    #     with self.assertRaises(ValueError):
    #         handle_post_convolution_function(
    #             config={**default_convolution_config},
    #             sfr_dict={},
    #             data_dict={},
    #             time_bin_info_dict={},
    #             convolution_instruction=convolution_instruction,
    #             convolution_results={},
    #             name="wrong_elements_in_convolution_results_list",
    #         )

    # def test_handle_post_convolution_function_no_name_in_convolution_results_lists(
    #     self,
    # ):

    #     convolution_instruction = {
    #         **default_convolution_instruction,
    #         "post_convolution_function": post_convolution_function_no_name_in_convolution_results_lists,
    #     }

    #     with self.assertRaises(ValueError):
    #         handle_post_convolution_function(
    #             config={**default_convolution_config},
    #             sfr_dict={},
    #             data_dict={},
    #             time_bin_info_dict={},
    #             convolution_instruction=convolution_instruction,
    #             convolution_results={},
    #             name="no_name_in_convolution_results_lists",
    #         )


if __name__ == "__main__":
    unittest.main()
