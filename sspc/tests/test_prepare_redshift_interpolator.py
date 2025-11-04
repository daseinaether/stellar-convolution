"""
Testcases for prepare_redshift_interpolator file
"""

import copy
import os
import pickle
import unittest

import astropy.units as u
import numpy as np
from scipy import interpolate

from syntheticstellarpopconvolve import default_convolution_config
from syntheticstellarpopconvolve.cosmology_utils import redshift_to_lookback_time
from syntheticstellarpopconvolve.general_functions import temp_dir
from syntheticstellarpopconvolve.prepare_redshift_interpolator import (
    create_interpolation_datasets,
    load_interpolation_data,
    prepare_redshift_interpolator,
)

TMP_DIR = temp_dir(
    "tests", "tests_convolution", "tests_prepare_redshift_interpolator", clean_path=True
)


class test_prepare_redshift_interpolator(unittest.TestCase):
    def setUp(self):
        self.convolution_config = copy.copy(default_convolution_config)

        # Set up SFR
        self.convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]),
            "starformation_rate_array": np.array([1, 1, 1, 1, 1])
            * u.Msun
            / u.yr
            / u.Gpc**3,
        }

        # set up convolution bins
        self.convolution_config["convolution_time_bin_edges"] = np.array(
            [0, 1, 2, 3, 4]
        )

        # lookback time convolution only
        self.convolution_config["time_type"] = "redshift"

        self.convolution_config["redshift_interpolator_use_log"] = True

        self.convolution_config["redshift_interpolator_data_output_filename"] = (
            os.path.join(TMP_DIR, "interpolator_dict.p")
        )

    def test_prepare_redshift_interpolator_general(self):
        config = prepare_redshift_interpolator(self.convolution_config)
        redshift_interpolator_dict = config["interpolators"]

        #
        self.assertTrue(
            "redshift_to_lookback_time_interpolator" in redshift_interpolator_dict
        )
        self.assertTrue(
            "lookback_time_to_redshift_interpolator" in redshift_interpolator_dict
        )

        #
        example = interpolate.interp1d(
            [1, 2],
            [3, 4],
            bounds_error=False,
            fill_value=0,
        )

        #
        self.assertTrue(
            isinstance(
                redshift_interpolator_dict["redshift_to_lookback_time_interpolator"],
                type(example),
            )
        )
        self.assertTrue(
            isinstance(
                redshift_interpolator_dict["lookback_time_to_redshift_interpolator"],
                type(example),
            )
        )

        #
        lookback_time_at_redshift_one = redshift_to_lookback_time(
            redshift=1, cosmology=self.convolution_config["cosmology"]
        )

        interpolated_redshift = redshift_interpolator_dict[
            "redshift_to_lookback_time_interpolator"
        ]([1])
        interpolated_lookback_time = redshift_interpolator_dict[
            "lookback_time_to_redshift_interpolator"
        ]([interpolated_redshift])

        self.assertAlmostEqual(
            lookback_time_at_redshift_one.value,
            interpolated_redshift[0],
            6,
            "lookback time at redshift=1 incorrect",
        )

        self.assertAlmostEqual(
            1, interpolated_lookback_time[0], 6, "lookback time at redshift=1 incorrect"
        )

    def test_prepare_redshift_interpolator_not_load(self):
        self.convolution_config["time_type"] = "lookback_time"

        config = prepare_redshift_interpolator(self.convolution_config)
        self.assertTrue("interpolators" not in config)


class test_create_interpolation_datasets(unittest.TestCase):
    def setUp(self):
        self.convolution_config = copy.copy(default_convolution_config)
        self.convolution_config["redshift_interpolator_data_output_filename"] = (
            os.path.join(TMP_DIR, "interpolator_dict.p")
        )

    def test_create_interpolation_datasets(self):
        create_interpolation_datasets(self.convolution_config)

        #
        self.assertTrue(
            os.path.isfile(
                self.convolution_config["redshift_interpolator_data_output_filename"]
            )
        )

        # Reload dict
        with open(
            self.convolution_config["redshift_interpolator_data_output_filename"], "rb"
        ) as f:
            interpolation_data_dict = pickle.load(f)

        self.assertTrue("redshift_data" in interpolation_data_dict)
        self.assertTrue("lookback_time_data" in interpolation_data_dict)
        self.assertTrue("min_redshift" in interpolation_data_dict)
        self.assertTrue("max_redshift" in interpolation_data_dict)
        self.assertTrue("redshift_stepsize" in interpolation_data_dict)
        self.assertTrue("interpolate_log" in interpolation_data_dict)


class test_load_interpolation_data(unittest.TestCase):
    def setUp(self):
        self.convolution_config = copy.copy(default_convolution_config)

        # Set up SFR
        self.convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]),
            "starformation_rate_array": np.array([1, 1, 1, 1, 1])
            * u.Msun
            / u.yr
            / u.Gpc**3,
        }

        # set up convolution bins
        self.convolution_config["convolution_time_bin_edges"] = np.array(
            [0, 1, 2, 3, 4]
        )

        # lookback time convolution only
        self.convolution_config["time_type"] = "redshift"

        self.convolution_config["redshift_interpolator_use_log"] = True

        self.convolution_config["redshift_interpolator_data_output_filename"] = (
            os.path.join(TMP_DIR, "interpolator_dict.p")
        )

    def test_general(self):
        redshift_interpolator_dict = load_interpolation_data(self.convolution_config)

        #
        self.assertTrue(
            "redshift_to_lookback_time_interpolator" in redshift_interpolator_dict
        )
        self.assertTrue(
            "lookback_time_to_redshift_interpolator" in redshift_interpolator_dict
        )


if __name__ == "__main__":
    unittest.main()
