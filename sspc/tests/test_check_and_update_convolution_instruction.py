"""
This is the unittest file for the check_and_update_convolution_instruction.py source file
"""

import unittest

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve import (
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.check_and_update_convolution_instruction import (
    check_and_update_convolution_instructions,
    check_convolution_instruction,
    check_delay_time_data_bin_info_dict,
    check_metallicity,
)


class test_check_delay_time_data_bin_info_dict(unittest.TestCase):
    def test_check_delay_time_data_bin_info_dict_missing_key(self):
        delay_time_data_bin_info_dict = {}

        with self.assertRaises(ValueError):
            check_delay_time_data_bin_info_dict(delay_time_data_bin_info_dict)

    def test_check_delay_time_data_bin_info_dict_wrong_unit(self):
        delay_time_data_bin_info_dict = {
            "delay_time_data_bin_edges": np.array([0, 1] * u.Msun)
        }

        with self.assertRaises(ValueError):
            check_delay_time_data_bin_info_dict(delay_time_data_bin_info_dict)

    # if "delay_time_data_bin_edges" not in delay_time_data_bin_info_dict:
    #     raise ValueError(
    #         "`delay_time_data_bin_edges` is required in the delay_time_data_bin_info_dict when convolving binned data"
    #     )

    # if not is_time_unit(delay_time_data_bin_info_dict["delay_time_data_bin_edges"]):
    #     raise ValueError("Please express 'delay_time_data_bin_edges' in units of time")


class test_check_metallicity(unittest.TestCase):
    def setUp(self):

        self.convolution_instruction_with_metallicity = {
            "data_column_dict": {"metallicity": "Fe/H"}
        }
        self.convolution_instruction_missing_metallicity = {
            "data_column_dict": {"no_metallicity_key": "some_value"}
        }

        self.config = default_convolution_config

    def test_check_metallicity_with_metallicity(self):
        config = {**self.config, "SFR_info": {"metallicity_": 1}}

        check_metallicity(
            convolution_config=config,
            convolution_instruction=self.convolution_instruction_with_metallicity,
        )
        # No exception should be raised

    def test_check_metallicity_missing_metallicity(self):
        config = {**self.config, "SFR_info": {"metallicity_": 1}}
        with self.assertRaises(ValueError):
            check_metallicity(
                convolution_config=config,
                convolution_instruction=self.convolution_instruction_missing_metallicity,
            )


class test_check_convolution_instruction(unittest.TestCase):
    def setUp(self):
        self.event_convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "event_data",
            "output_data_name": "output_event_data",
            "convolution_type": "integrate",
            "data_column_dict": {"delay_time": "delay", "normalized_yield": "rate"},
        }

        self.config = default_convolution_config

    def test_check_convolution_instruction_event_type(self):
        check_convolution_instruction(
            convolution_instruction=self.event_convolution_instruction,
            convolution_config=self.config,
        )
        # No exception should be raised

    def test_check_convolution_instruction_missing_event_required_key(self):
        event_convolution_instruction_missing_key = {
            **default_convolution_instruction,
            "input_data_name": "event_data",
            "data_column_dict": {"normalized_yield": "rate"},
        }
        with self.assertRaises(ValueError):
            check_convolution_instruction(
                convolution_instruction=event_convolution_instruction_missing_key,
                convolution_config=self.config,
            )

    def test_check_convolution_instruction_event_missing_metallicity(self):
        event_convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "event_data",
            "output_data_name": "output_event_data",
            "convolution_type": "integrate",
            "data_column_dict": {"delay_time": "delay", "normalized_yield": "rate"},
        }

        config = {**self.config, "SFR_info": {"metallicity_": 1}}
        with self.assertRaises(ValueError):
            check_convolution_instruction(
                convolution_instruction=event_convolution_instruction,
                convolution_config=config,
            )


class test_check_and_update_convolution_instructions(unittest.TestCase):

    def test_check_and_update_convolution_instructions_no_convolution_instructions(
        self,
    ):

        config = default_convolution_config

        with self.assertRaises(ValueError):
            check_and_update_convolution_instructions(convolution_config=config)


if __name__ == "__main__":
    unittest.main()
