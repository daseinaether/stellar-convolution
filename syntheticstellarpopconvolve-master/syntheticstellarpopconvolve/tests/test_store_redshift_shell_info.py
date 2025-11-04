"""
Testcases for store_redshift_shell_info file
"""

import copy
import unittest

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve import default_convolution_config
from syntheticstellarpopconvolve.general_functions import temp_dir
from syntheticstellarpopconvolve.store_redshift_shell_info import (
    create_shell_volume_dict,
    store_redshift_shell_info,
)

TMP_DIR = temp_dir(
    "tests", "tests_convolution", "tests_store_redshift_shell_info", clean_path=True
)


class test_store_redshift_shell_info(unittest.TestCase):
    def test_store_redshift_shell_info_redshift(self):

        #
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["time_type"] = "redshift"

        # Set up SFR
        convolution_config["SFR_info"] = {
            "redshift_bin_edges": np.array([0, 1, 2, 3, 4, 5]),
            "starformation_rate_array": np.array([1, 1, 1, 1, 1])
            * u.Msun
            / u.yr
            / u.Gpc**3,
        }

        #
        convolution_config = store_redshift_shell_info(
            convolution_config, sfr_dict=convolution_config["SFR_info"]
        )

        #
        self.assertTrue("redshift_shell_volume_dict" in convolution_config.keys())
        self.assertTrue(len(convolution_config["redshift_shell_volume_dict"]) > 0)
        self.assertTrue(0.5 in convolution_config["redshift_shell_volume_dict"].keys())
        self.assertTrue(
            "shell_volume" in convolution_config["redshift_shell_volume_dict"][0.5]
        )
        self.assertTrue(
            "lower_edge_shell_redshift"
            in convolution_config["redshift_shell_volume_dict"][0.5]
        )
        self.assertTrue(
            "upper_edge_shell_redshift"
            in convolution_config["redshift_shell_volume_dict"][0.5]
        )
        self.assertTrue(
            "lower_edge_shell_lookback_time"
            in convolution_config["redshift_shell_volume_dict"][0.5]
        )
        self.assertTrue(
            "upper_edge_shell_lookback_time"
            in convolution_config["redshift_shell_volume_dict"][0.5]
        )
        self.assertTrue(
            "center_shell" in convolution_config["redshift_shell_volume_dict"][0.5]
        )
        self.assertTrue(
            "delta_shell_lookback_time"
            in convolution_config["redshift_shell_volume_dict"][0.5]
        )
        self.assertTrue(
            "delta_shell_redshift"
            in convolution_config["redshift_shell_volume_dict"][0.5]
        )

    def test_store_redshift_shell_info_no_redshift(self):

        #
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["time_type"] = "lookback_time"

        # Set up SFR
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]),
            "starformation_rate_array": np.array([1, 1, 1, 1, 1])
            * u.Msun
            / u.yr
            / u.Gpc**3,
        }

        #
        convolution_config = store_redshift_shell_info(
            convolution_config, sfr_dict=convolution_config["SFR_info"]
        )

        self.assertTrue("redshift_shell_volume_dict" not in convolution_config.keys())


class test_create_shell_volume_dict(unittest.TestCase):
    def test_create_shell_volume_dict(self):

        #
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["time_type"] = "redshift"

        #
        shell_volume_dict = create_shell_volume_dict(
            redshift_bin_edges=np.array([0, 1, 2, 3, 4, 5]),
            config=convolution_config,
        )

        #
        self.assertTrue(0.5 in shell_volume_dict.keys())
        self.assertTrue(1.5 in shell_volume_dict.keys())
        self.assertTrue(2.5 in shell_volume_dict.keys())
        self.assertTrue(3.5 in shell_volume_dict.keys())
        self.assertTrue(4.5 in shell_volume_dict.keys())
        self.assertTrue("shell_volume" in shell_volume_dict[0.5])
        self.assertTrue("lower_edge_shell_redshift" in shell_volume_dict[0.5])
        self.assertTrue("upper_edge_shell_redshift" in shell_volume_dict[0.5])
        self.assertTrue("lower_edge_shell_lookback_time" in shell_volume_dict[0.5])
        self.assertTrue("upper_edge_shell_lookback_time" in shell_volume_dict[0.5])
        self.assertTrue("center_shell" in shell_volume_dict[0.5])
        self.assertTrue("delta_shell_lookback_time" in shell_volume_dict[0.5])


# def store_redshift_shell_info(config, sfr_dict):
#     """
#     Function to add the redshift shell info dict to the hdf5 file
#     """

#     if config["time_type"] == "redshift":
#         ##################
#         # Create shell volume dict
#         redshift_shell_volume_dict = create_shell_volume_dict(
#             redshift_bin_edges=sfr_dict["redshift_bin_edges"],
#             config=config,
#         )
#         sfr_dict["redshift_shell_volume_dict"] = redshift_shell_volume_dict

#     return sfr_dict


if __name__ == "__main__":
    unittest.main()
