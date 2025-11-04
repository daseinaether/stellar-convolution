"""
Testcases for check_convolution_config file

TODO: rebuild this
"""

import copy
import os
import unittest

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve import (
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.check_and_update_convolution_config import (
    check_and_update_convolution_config,
    check_convolution_config,
)
from syntheticstellarpopconvolve.general_functions import temp_dir

TMP_DIR = temp_dir(
    "tests",
    "tests_convolution",
    "tests_check_and_update_convolution_config",
    clean_path=True,
)


class test_check_convolution_config(unittest.TestCase):
    def test_check_convolution_config_with_valid_input(self):
        config_with_valid_convolution = {
            "check_convolution_config": True,
            "logger": default_convolution_config["logger"],
            "time_type": "lookback_time",
            "convolution_lookback_time_bin_edges": np.array([1, 2]) * u.yr,
            "convolution_instructions": [
                {
                    **default_convolution_instruction,
                    "input_data_name": "event_data",
                    "output_data_name": "output_event_data",
                    "convolution_type": "integrate",
                    "data_column_dict": {
                        "delay_time": "delay",
                        "normalized_yield": "rate",
                        "metallicity": "metallicity",
                    },
                },
            ],
            "SFR_info": [
                {
                    "name": "test",
                    "lookback_time_bin_edges": np.array([0, 1, 2, 3]) * 1e9 * u.yr,
                    "starformation_rate_array": np.array([1, 2, 3]) * u.Msun / u.yr,
                    "metallicity_bin_edges": np.array([0.1, 0.2, 0.3, 0.4]),
                    "metallicity_distribution_array": np.array(
                        [[0.5, 0.6, 0.7], [0.5, 0.6, 0.7], [0.5, 0.6, 0.7]]
                    ),
                }
            ],
        }

        check_convolution_config(config_with_valid_convolution)
        # No exception should be raised

    def test_check_convolution_config_missing_convolution_instruction(self):

        config_with_missing_convolution_instruction = {
            "check_convolution_config": True,
            "logger": default_convolution_config["logger"],
            "time_type": "redshift",
            "convolution_instructions": [],  # Missing convolution instruction
            "SFR_info": [
                {
                    "lookback_time_bin_edges": [0, 1, 2, 3],
                    "starformation_rate_array": [1, 2, 3] * u.Msun / u.yr,
                    "metallicity_bin_edges": [0.1, 0.2, 0.3],
                    "metallicity_distribution_array": [
                        [0.5, 0.6, 0.7],
                        [0.5, 0.6, 0.7],
                        [0.5, 0.6, 0.7],
                    ],
                }
            ],
        }

        with self.assertRaises(ValueError):
            check_convolution_config(config_with_missing_convolution_instruction)

    def test_check_convolution_config_missing_SFR_info(self):
        #
        config_with_missing_SFR_info = {
            "check_convolution_config": True,
            "logger": default_convolution_config["logger"],
            "time_type": "redshift",
            "convolution_instructions": [
                {
                    "convolution_type": "integrate",
                    "input_data_name": "event_data",
                    "output_data_name": "output_event_data",
                    "data_column_dict": {
                        "delay_time": "delay",
                        "normalized_yield": "rate",
                    },
                }
            ],
            "SFR_info": [],  # Missing SFR info
        }

        with self.assertRaises(ValueError):
            check_convolution_config(config_with_missing_SFR_info)


class test_update_convolution_config(unittest.TestCase):
    """ """

    def test_update_convolution_config(self):

        config = copy.copy(default_convolution_config)
        config["redshift_interpolator_data_output_filename"] = os.path.join(
            TMP_DIR, "interpolator_dict.p"
        )
        #
        config["input_filename"] = os.path.join(TMP_DIR, "input_hdf5_sfr_only.h5")
        config["output_filename"] = os.path.join(TMP_DIR, "output_hdf5_sfr_only.h5")
        config["time_type"] = "redshift"

        # Set up SFR
        config["SFR_info"] = {
            "redshift_bin_edges": np.array([0, 1, 2, 3, 4, 5]),
            "starformation_rate_array": np.array([1, 1, 1, 1, 1])
            * u.Msun
            / u.yr
            / u.Gpc**3,
        }

        #
        config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "dummy",
                "output_data_name": "dummy",
                "convolution_type": "integrate",
                "data_column_dict": {
                    "delay_time": "delay_time",
                    "normalized_yield": "probability",
                },
            },
        ]

        with open(config["input_filename"], "a") as _:
            pass

        config["convolution_redshift_bin_edges"] = np.array([1, 2])

        #
        check_and_update_convolution_config(config=config)

        #
        self.assertEqual(config["convolution_time_bin_centers"].tolist(), [1.5])


if __name__ == "__main__":
    unittest.main()
