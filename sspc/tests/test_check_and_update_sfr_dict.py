"""
TODO: there is some issue here for situations where the metallicity distribution is a square.
"""

import logging
import unittest

import astropy.units as u
import numpy as np
from astropy.cosmology import Planck13 as cosmo  # Planck 2013

from syntheticstellarpopconvolve import default_convolution_config
from syntheticstellarpopconvolve.check_and_update_sfr_dict import (
    check_sfr_dict,
    pad_sfr_dict,
    update_sfr_dict,
)
from syntheticstellarpopconvolve.general_functions import temp_dir

TMP_DIR = temp_dir(
    "tests",
    "tests_convolution",
    "tests_check_and_update_sfr_dict",
    clean_path=True,
)


class test_pad_sfr_dict(unittest.TestCase):
    def setUp(self):

        logger = logging.getLogger(__name__)
        FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s ] %(asctime)s: %(message)s"
        logging.basicConfig(format=FORMAT)
        logger.setLevel(logging.INFO)

        self.config = {
            "logger": logger,
            "time_type": "lookback_time",
            "cosmology": cosmo,
        }  # Example config
        self.sfr_dict = {
            "lookback_time_bin_edges": np.array([1, 2, 3, 4]),
            "redshift_bin_edges": np.array([0.1, 0.2, 0.3, 0.4]),
            "starformation_rate_array": np.array([10, 20, 30]),
            "metallicity_bin_edges": np.array([0.01, 0.1, 0.2]),
            "metallicity_distribution_array": np.array([[1, 2, 3], [4, 5, 6]]).T,
        }

    def test_pad_sfr_dict_lookback_time(self):
        padded_sfr_dict = pad_sfr_dict(self.config, self.sfr_dict)
        self.assertTrue("padded_lookback_time_bin_edges" in padded_sfr_dict)
        self.assertTrue(
            np.array_equal(
                padded_sfr_dict["padded_lookback_time_bin_edges"],
                np.array([1 - 1e13, 1, 2, 3, 4, 4 + 1e13]),
            )
        )
        self.assertTrue("padded_starformation_rate_array" in padded_sfr_dict)
        self.assertTrue(
            np.array_equal(
                padded_sfr_dict["padded_starformation_rate_array"],
                np.array([0, 10, 20, 30, 0]),
            )
        )

    def test_pad_sfr_dict_redshift(self):
        self.config["time_type"] = "redshift"
        padded_sfr_dict = pad_sfr_dict(self.config, self.sfr_dict)
        self.assertTrue("padded_redshift_bin_edges" in padded_sfr_dict)
        self.assertTrue(
            np.array_equal(
                padded_sfr_dict["padded_redshift_bin_edges"],
                np.array([0.1 - 1e13, 0.1, 0.2, 0.3, 0.4, 0.4 + 1e13]),
            )
        )
        self.assertTrue("padded_starformation_rate_array" in padded_sfr_dict)
        self.assertTrue(
            np.array_equal(
                padded_sfr_dict["padded_starformation_rate_array"],
                np.array([0, 10, 20, 30, 0]),
            )
        )

    def test_pad_sfr_dict_metallicity(self):
        self.sfr_dict = update_sfr_dict(sfr_dict=self.sfr_dict, config=self.config)
        padded_sfr_dict = pad_sfr_dict(self.config, self.sfr_dict)
        self.assertTrue("padded_metallicity_bin_edges" in padded_sfr_dict)

        self.assertTrue(
            np.array_equal(
                padded_sfr_dict["padded_metallicity_bin_edges"],
                np.array([-1e-20, 0.01, 0.1, 0.2, 2.0]),
            )
        )

        #
        self.assertTrue("padded_metallicity_distribution_array" in padded_sfr_dict)
        expected_array = np.array(
            [[0, 0, 0, 0, 0], [0, 1, 2, 3, 0], [0, 4, 5, 6, 0], [0, 0, 0, 0, 0]]
        ).T
        self.assertTrue(
            np.array_equal(
                padded_sfr_dict["padded_metallicity_distribution_array"],
                expected_array,
            )
        )

        #
        self.assertTrue(
            "padded_metallicity_weighted_starformation_rate_array" in padded_sfr_dict
        )
        expected_array = (
            np.array(
                [[0, 0, 0, 0, 0], [0, 1, 2, 3, 0], [0, 4, 5, 6, 0], [0, 0, 0, 0, 0]]
            )
            * np.array([0, 10, 20, 30, 0])[np.newaxis, :]
            * np.diff(padded_sfr_dict["padded_metallicity_bin_edges"])[:, np.newaxis]
        ).T
        self.assertTrue(
            np.array_equal(
                padded_sfr_dict["padded_metallicity_weighted_starformation_rate_array"],
                expected_array,
            )
        )

    def test_pad_sfr_dict_metallicity_square(self):
        sfr_dict = {
            "lookback_time_bin_edges": np.array([1, 2, 3, 4]),
            "starformation_rate_array": np.array([10, 20, 30]),
            "metallicity_bin_edges": np.array([0.01, 0.1, 0.2, 0.3]),
            "metallicity_distribution_array": np.array(
                [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
            ).T,
        }

        sfr_dict = update_sfr_dict(sfr_dict=sfr_dict, config=self.config)

        # TODO: add explicit check for what the content should be

    def test_pad_sfr_dict_metallicity_non_square(self):
        sfr_dict = {
            "lookback_time_bin_edges": np.array([1, 2, 3, 4]),
            "starformation_rate_array": np.array([10, 20, 30]),
            "metallicity_bin_edges": np.array([0.01, 0.1, 0.2]),
            # "metallicity_distribution_array": np.array([[1, 2, 3], [4, 5, 6]]),
            "metallicity_distribution_array": np.array([[1, 2], [3, 4], [5, 6]]),
        }

        sfr_dict = update_sfr_dict(sfr_dict=sfr_dict, config=self.config)

        # TODO: add explicit check for what the content should be


class test_check_sfr_dict(unittest.TestCase):
    def setUp(self):
        self.sfr_dict = {
            "name": "test_sfr_dict",
            "lookback_time_bin_edges": np.array([1, 2, 3]) * 1e9 * u.yr,
            "starformation_rate_array": np.array([10, 20]) * u.Msun / u.yr,
            "metallicity_bin_edges": np.array([0.01, 0.1, 0.2, 0.3]),
            "metallicity_distribution_array": np.array([[1, 2, 3], [4, 5, 6]]),
        }

        self.config = default_convolution_config

    def test_check_sfr_dict_with_name(self):
        requires_name = True
        requires_metallicity_info = True
        time_type = "lookback_time"

        # No exception should be raised
        check_sfr_dict(
            sfr_dict=self.sfr_dict,
            requires_name=requires_name,
            requires_metallicity_info=requires_metallicity_info,
            time_type=time_type,
            config=self.config,
        )

    def test_check_sfr_dict_without_name(self):
        requires_name = True
        requires_metallicity_info = True
        time_type = "lookback_time"

        del self.sfr_dict["name"]  # Removing the name key

        with self.assertRaises(ValueError):
            check_sfr_dict(
                sfr_dict=self.sfr_dict,
                requires_name=requires_name,
                requires_metallicity_info=requires_metallicity_info,
                time_type=time_type,
                config=self.config,
            )

    def test_check_sfr_dict_without_metallicity_info(self):
        requires_name = True
        requires_metallicity_info = True
        time_type = "lookback_time"

        del self.sfr_dict[
            "metallicity_bin_edges"
        ]  # Removing the metallicity_bin_edges key

        with self.assertRaises(ValueError):
            check_sfr_dict(
                sfr_dict=self.sfr_dict,
                requires_name=requires_name,
                requires_metallicity_info=requires_metallicity_info,
                time_type=time_type,
                config=self.config,
            )

    def test_check_sfr_dict_without_time_type_info(self):
        requires_name = True
        requires_metallicity_info = True
        time_type = "lookback_time"

        del self.sfr_dict[
            "lookback_time_bin_edges"
        ]  # Removing the lookback_time_bin_edges key

        with self.assertRaises(ValueError):
            check_sfr_dict(
                sfr_dict=self.sfr_dict,
                requires_name=requires_name,
                requires_metallicity_info=requires_metallicity_info,
                time_type=time_type,
                config=self.config,
            )

    def test_check_sfr_dict_without_lookback_time_unit(self):
        requires_name = True
        requires_metallicity_info = True
        time_type = "lookback_time"

        self.sfr_dict["lookback_time_bin_edges"] = np.array([1, 2, 3]) * 1e9

        with self.assertRaises(ValueError):
            check_sfr_dict(
                sfr_dict=self.sfr_dict,
                requires_name=requires_name,
                requires_metallicity_info=requires_metallicity_info,
                time_type=time_type,
                config=self.config,
            )

    def test_check_sfr_dict_redshift_wrong_bin_edges(self):
        requires_name = True
        requires_metallicity_info = True
        time_type = "redshift"

        with self.assertRaises(ValueError):
            check_sfr_dict(
                sfr_dict=self.sfr_dict,
                requires_name=requires_name,
                requires_metallicity_info=requires_metallicity_info,
                time_type=time_type,
                config=self.config,
            )

    def test_check_sfr_dict_redshift_correct_bin_edges(self):
        self.sfr_dict["redshift_bin_edges"] = np.array([0, 1, 2])

        requires_name = True
        requires_metallicity_info = True
        time_type = "redshift"

        check_sfr_dict(
            sfr_dict=self.sfr_dict,
            requires_name=requires_name,
            requires_metallicity_info=requires_metallicity_info,
            time_type=time_type,
            config=self.config,
        )


class test_update_sfr_dict(unittest.TestCase):
    def test_update_sfr_dict_lookback_time(self):
        config = {"logger": logging.getLogger(__name__), "time_type": "lookback_time"}

        sfr_dict = {
            "lookback_time_bin_edges": np.array([1, 2, 3, 4]),
        }

        sfr_dict = update_sfr_dict(sfr_dict, config)

        # Expected output
        expected_bin_centers = np.array([1.5, 2.5, 3.5])

        # Assertions
        assert (
            "lookback_time_bin_centers" in sfr_dict
        ), "Key 'lookback_time_bin_centers' is missing in the updated dictionary"
        assert (
            "time_bin_centers" in sfr_dict
        ), "Key 'time_bin_centers' is missing in the updated dictionary"

        np.testing.assert_array_almost_equal(
            sfr_dict["lookback_time_bin_centers"], expected_bin_centers, decimal=10
        )
        np.testing.assert_array_almost_equal(
            sfr_dict["time_bin_centers"], expected_bin_centers, decimal=10
        )

    def test_update_sfr_dict_redshift(self):
        config = {
            "logger": logging.getLogger(__name__),
            "time_type": "redshift",
            "cosmology": cosmo,
        }

        sfr_dict = {
            "redshift_bin_edges": np.array([0.1, 0.2, 0.3, 0.4]),
        }

        sfr_dict = update_sfr_dict(sfr_dict, config)

        # Expected output
        expected_bin_centers = np.array([0.15, 0.25, 0.35])

        # Assertions
        assert (
            "redshift_bin_centers" in sfr_dict
        ), "Key 'redshift_bin_centers' is missing in the updated dictionary"
        assert (
            "time_bin_centers" in sfr_dict
        ), "Key 'time_bin_centers' is missing in the updated dictionary"

        np.testing.assert_array_almost_equal(
            sfr_dict["redshift_bin_centers"], expected_bin_centers, decimal=10
        )
        np.testing.assert_array_almost_equal(
            sfr_dict["time_bin_centers"], expected_bin_centers, decimal=10
        )

    def test_update_sfr_dict_metallicity_info(self):
        config = {"logger": logging.getLogger(__name__), "time_type": "lookback_time"}

        sfr_dict = {
            "lookback_time_bin_edges": np.array([1, 2, 3, 4]),
            "starformation_rate_array": np.array([10, 20, 30]),
            "metallicity_bin_edges": np.array([0.01, 0.1, 0.2]),
            "metallicity_distribution_array": np.array([[1, 2, 3], [4, 5, 6]]).T,
        }

        # When calling the function
        sfr_dict = update_sfr_dict(sfr_dict, config)

        # Expected outputs
        expected_metallicity_bin_centers = np.array([0.055, 0.15])
        expected_metallicity_bin_sizes = np.array([0.09, 0.1])
        expected_weighted_sfr = np.array(
            [
                [10 * 0.09 * 1, 10 * 0.1 * 4],  # First row
                [20 * 0.09 * 2, 20 * 0.1 * 5],  # Second row
                [30 * 0.09 * 3, 30 * 0.1 * 6],  # Third row
            ]
        )

        # Assertions
        assert (
            "include_metallicity_info" in sfr_dict
        ), "Key 'include_metallicity_info' is missing in the updated dictionary"
        assert (
            sfr_dict["include_metallicity_info"] is True
        ), "'include_metallicity_info' should be True when 'metallicity_distribution_array' is present"

        assert (
            "metallicity_bin_centers" in sfr_dict
        ), "Key 'metallicity_bin_centers' is missing"
        np.testing.assert_array_almost_equal(
            sfr_dict["metallicity_bin_centers"],
            expected_metallicity_bin_centers,
            err_msg="metallicity_bin_centers values are incorrect",
        )

        assert (
            "metallicity_bin_sizes" in sfr_dict
        ), "Key 'metallicity_bin_sizes' is missing"
        np.testing.assert_array_almost_equal(
            sfr_dict["metallicity_bin_sizes"],
            expected_metallicity_bin_sizes,
            err_msg="metallicity_bin_sizes values are incorrect",
        )

        assert (
            "metallicity_weighted_starformation_rate_array" in sfr_dict
        ), "Key 'metallicity_weighted_starformation_rate_array' is missing"
        np.testing.assert_array_almost_equal(
            sfr_dict["metallicity_weighted_starformation_rate_array"],
            expected_weighted_sfr,
            err_msg="metallicity_weighted_starformation_rate_array values are incorrect",
        )

    def test_update_sfr_dict_store_redshift_shell_info(self):
        config = {
            "logger": logging.getLogger(__name__),
            "time_type": "redshift",
            "cosmology": cosmo,
        }
        sfr_dict = {
            "redshift_bin_edges": np.array([0.1, 0.2, 0.3, 0.4]),
        }
        sfr_dict = update_sfr_dict(sfr_dict, config)

        assert (
            "redshift_shell_volume_dict" in sfr_dict.keys()
        ), "'redshift_shell_volume_dict' is not present in the sfr dictionary"


if __name__ == "__main__":
    unittest.main()
