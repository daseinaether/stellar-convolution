"""
Tests for backward convolution of binned data.
"""

import copy
import logging
import os
import unittest

import astropy.units as u
import h5py
import numpy as np
import pandas as pd

from syntheticstellarpopconvolve import (
    convolve,
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.general_functions import (
    extract_unit_dict,
    generate_boilerplate_outputfile,
    temp_dir,
)

np.random.seed(0)

TMP_DIR = temp_dir(
    "tests",
    "tests_convolution",
    "tests_convolve_nonbinned_data_with_backward_convolution",
    clean_path=True,
)


def integrate_post_convolution_function(
    config,
    sfr_dict,
    data_dict,
    convolution_results,
    convolution_instruction,
):
    convolution_results["yield"] = convolution_results["yield"] * 0

    return convolution_results


def sample_post_convolution_function(
    config,
    sfr_dict,
    data_dict,
    convolution_results,
    convolution_instruction,
):

    convolution_results["sampled_value"] = data_dict["value"][
        convolution_results["sampled_indices"]
    ]

    return convolution_results


class test_convolve_nonbinned_data_with_backward_convolution(unittest.TestCase):
    """
    TODO: make a more complicated post convolution hook function test
    """

    def setUp(self):
        #
        output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5_sfr_only.h5")
        generate_boilerplate_outputfile(output_hdf5_filename)

        ##############
        # SET UP DATA
        self.dummy_data = {
            "delay_time": np.array([0, 1, 2, 3]),
            "value": np.array([3, 2, 1, 0]),
            "probability": np.array([1, 2, 3, 4]),
        }
        dummy_df = pd.DataFrame.from_records(self.dummy_data)

        ##############
        # Store data in pandas
        dummy_df.to_hdf(output_hdf5_filename, key="input_data/{}".format("dummy"))

        ###################
        #
        self.convolution_config = copy.copy(default_convolution_config)
        self.convolution_config["logger"].setLevel(logging.CRITICAL)
        self.convolution_config["output_filename"] = output_hdf5_filename
        self.convolution_config["redshift_interpolator_data_output_filename"] = (
            os.path.join(TMP_DIR, "interpolator_dict.p")
        )
        self.convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")

        # Set up SFR
        self.convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1, 1, 1])
            * u.Msun
            / u.yr
            / u.Gpc**3,
        }

        # set up convolution bins
        self.convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1]) * u.yr
        )

    def test_convolve_nonbinned_data_with_backward_convolution_integrate(self):
        #
        normal_convolution_instructions = {
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
        self.convolution_config["convolution_instructions"] = [
            normal_convolution_instructions
        ]

        #
        convolve(self.convolution_config)

        with h5py.File(
            self.convolution_config["output_filename"], "r"
        ) as output_hdf5file:
            groupname = "output_data/dummy/dummy/convolution_results/0.5 yr/"

            data = output_hdf5file[groupname + "/yield"][()]
            unit_dict = extract_unit_dict(output_hdf5file, groupname)

            np.testing.assert_array_equal(
                data * unit_dict["yield"],
                np.array([1, 2, 3, 4.0]) * (1.0 / u.yr / u.Gpc**3),
            )

    def test_convolve_nonbinned_data_with_backward_convolution_integrate_post_convolution_simple(
        self,
    ):

        #
        convolution_instructions = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "convolution_type": "integrate",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
            "post_convolution_function": integrate_post_convolution_function,
        }

        #
        self.convolution_config["convolution_instructions"] = [convolution_instructions]

        #
        convolve(self.convolution_config)

        with h5py.File(
            self.convolution_config["output_filename"], "r"
        ) as output_hdf5file:
            groupname = "output_data/dummy/dummy/convolution_results/0.5 yr/"

            data = output_hdf5file[groupname + "/yield"][()]
            unit_dict = extract_unit_dict(output_hdf5file, groupname)

            np.testing.assert_array_equal(
                data * unit_dict["yield"],
                np.zeros(self.dummy_data["probability"].shape)
                * (1.0 / u.yr / u.Gpc**3),
            )

    def test_convolve_nonbinned_data_with_backward_convolution_sample_fractional(self):

        #
        convolution_instructions = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "convolution_type": "sample",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
            "multiply_by_sfr_time_binsize": True,
        }

        # Set up SFR
        self.convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1.1, 2.2, 3.3, 4.4, 5.5])
            * u.Msun
            / u.yr,
        }

        self.convolution_config["convolution_instructions"] = [convolution_instructions]

        #
        convolve(self.convolution_config)

        with h5py.File(
            self.convolution_config["output_filename"], "r"
        ) as output_hdf5file:
            groupname = "output_data/dummy/dummy/convolution_results/0.5 yr/"

            data = output_hdf5file[groupname + "/sampled_indices"][()]

            unique, counts = np.unique(data, return_counts=True)

            np.testing.assert_array_equal(unique, np.array([0, 1, 2, 3]))
            np.testing.assert_array_equal(counts, np.array([1, 4, 10, 18]))

    def test_convolve_nonbinned_data_with_backward_convolution_sample_integer(self):

        #
        convolution_instructions = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "convolution_type": "sample",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
            "multiply_by_sfr_time_binsize": True,
        }

        # Set up SFR
        self.convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
        }

        self.convolution_config["convolution_instructions"] = [convolution_instructions]

        #
        convolve(self.convolution_config)

        with h5py.File(
            self.convolution_config["output_filename"], "r"
        ) as output_hdf5file:
            groupname = "output_data/dummy/dummy/convolution_results/0.5 yr/"

            data = output_hdf5file[groupname + "/sampled_indices"][()]
            unique, counts = np.unique(data, return_counts=True)

            np.testing.assert_array_equal(unique, np.array([0, 1, 2, 3]))
            np.testing.assert_array_equal(counts, np.array([1, 4, 9, 16]))

    def test_convolve_nonbinned_data_with_backward_convolution_sample_integer_post_convolution_simple(
        self,
    ):

        #
        convolution_instructions = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "convolution_type": "sample",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
                "value": "value",
            },
            "multiply_by_sfr_time_binsize": True,
            "post_convolution_function": sample_post_convolution_function,
        }

        # Set up SFR
        self.convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 2, 3, 4, 5]) * u.Msun / u.yr,
        }

        self.convolution_config["convolution_instructions"] = [convolution_instructions]

        #
        convolve(self.convolution_config)

        with h5py.File(
            self.convolution_config["output_filename"], "r"
        ) as output_hdf5file:
            groupname = "output_data/dummy/dummy/convolution_results/0.5 yr/"

            data = output_hdf5file[groupname + "/sampled_value"][()]

            unique, counts = np.unique(data, return_counts=True)

            np.testing.assert_array_equal(unique, np.array([0, 1, 2, 3]))
            np.testing.assert_array_equal(counts, np.array([16, 9, 4, 1]))


if __name__ == "__main__":
    unittest.main()
