"""
"""

import logging
import unittest

import astropy.units as u
import h5py
import numpy as np
import pandas as pd
import pkg_resources

from syntheticstellarpopconvolve import convolve, default_convolution_instruction
from syntheticstellarpopconvolve.convolution_by_sampling import (
    select_dict_entries_with_new_indices,
)
from syntheticstellarpopconvolve.general_functions import temp_dir
from syntheticstellarpopconvolve.tests.utils import Boilerplate

TMP_DIR = temp_dir(
    "tests",
    "tests_convolution",
    "general_tests",
    clean_path=True,
)

np.random.seed(0)


##################
# Non unit-tests but sanity checks
#


####
# Define
def postprocessing_multiple_dicts(
    config, sfr_dict, data_dict, convolution_results, convolution_instruction
):
    """
    Post-convolution function to handle integrating the systems forward in time and finding those that end up in the LISA waveband.

    using local_indices to select everything and using Alexey's distance sampler to handle sampling the distances
    """

    # unpack data
    system_indices = convolution_results["sampled_indices"]
    local_indices = np.arange(len(system_indices))

    # add distances
    convolution_results["dists"] = (
        np.random.randint(0, 1e6, size=len(local_indices)) * u.kpc
    )

    # shuffle and select 2 sets
    shuffled_indices = np.arange(len(system_indices))
    np.random.shuffle(shuffled_indices)
    random_length = np.random.randint(0, high=len(shuffled_indices))

    # Select all results from set 1
    convolution_result_1 = select_dict_entries_with_new_indices(
        sampled_data_dict=convolution_results,
        new_indices=shuffled_indices[:random_length],
    )

    # Select all results from set 2
    convolution_result_2 = select_dict_entries_with_new_indices(
        sampled_data_dict=convolution_results,
        new_indices=shuffled_indices[random_length:],
    )

    # split into two
    convolution_results = [convolution_result_1, convolution_result_2]

    return convolution_results


def postprocessing_multiple_dicts_with_name(
    config, sfr_dict, data_dict, convolution_results, convolution_instruction
):
    convolution_results = postprocessing_multiple_dicts(
        config=config,
        sfr_dict=sfr_dict,
        data_dict=data_dict,
        convolution_results=convolution_results,
        convolution_instruction=convolution_instruction,
    )

    convolution_results[0]["name"] = "set_1"
    convolution_results[1]["name"] = "set_2"

    return convolution_results


def postprocessing_multiple_dicts_without_name(
    config, sfr_dict, data_dict, convolution_results, convolution_instruction
):
    convolution_results = postprocessing_multiple_dicts(
        config=config,
        sfr_dict=sfr_dict,
        data_dict=data_dict,
        convolution_results=convolution_results,
        convolution_instruction=convolution_instruction,
    )

    return convolution_results


class test_postprocessing(unittest.TestCase, Boilerplate):

    def setUp(self):
        self.setup(
            name="test_postprocessing",
            tmp_dir=TMP_DIR,
            add_population_settings=False,
            sfr_unit=u.Msun / u.yr,
        )
        self.convolution_config["logger"].setLevel(logging.CRITICAL)

        ###################
        # Set up data
        BinCodex_events_filename = pkg_resources.resource_filename(
            "syntheticstellarpopconvolve",
            "example_data/example_BinCodex_dwd.h5",
        )

        #
        BinCodex_T0_events = pd.read_hdf(
            BinCodex_events_filename,
            "T0",
        )

        ##################
        # update T0 output

        # get mass normalisation
        mass_normalisation_fiducial = 4476544.539875359 * u.Msun

        # set normalised yield
        BinCodex_T0_events["normalized_yield"] = 1 / mass_normalisation_fiducial

        # Query the dataset to select the formation of the WDs

        # to check if things start with some number its easier to turn them into strings
        BinCodex_T0_events["str_event"] = BinCodex_T0_events["event"].astype(str)
        BinCodex_T0_events["str_type1"] = BinCodex_T0_events["type1"].astype(str)
        BinCodex_T0_events["str_type2"] = BinCodex_T0_events["type2"].astype(str)

        # first, lets query the type-changing events. Any type-change will do
        wd_binaries = BinCodex_T0_events.query("str_event.str.startswith('1')")

        # The type should change to a WD-type (and the other should already be one)
        wd_binaries = wd_binaries.query("str_type1.str.startswith('2')")
        wd_binaries = wd_binaries.query("str_type2.str.startswith('2')")

        # lets delete the string versions of the columns again
        wd_binaries = wd_binaries.drop(columns=["str_event", "str_type1", "str_type2"])

        # lets also delete the original dataframe
        del BinCodex_T0_events

        # store the data frame in the hdf5file
        wd_binaries.to_hdf(
            self.convolution_config["output_filename"], key="input_data/dummy"
        )

    def test_postprocessing_multiple_dictionaries_with_name(self):

        #
        self.convolution_config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "dummy",
                "output_data_name": "dummy",
                "convolution_type": "sample",
                "data_column_dict": {
                    # required
                    "normalized_yield": "normalized_yield",
                    "delay_time": {"column_name": "time", "unit": u.Myr},
                },
                "post_convolution_function": postprocessing_multiple_dicts_with_name,
                "multiply_by_sfr_time_binsize": True,
            },
        ]

        # convolve
        convolve(config=self.convolution_config)

        # read out content and integrate until today
        with h5py.File(
            self.convolution_config["output_filename"], "r"
        ) as output_hdf5_file:

            self.assertTrue(
                "set_1"
                in output_hdf5_file[
                    "output_data/dummy/dummy/convolution_results/"
                ].keys()
            )
            self.assertTrue(
                "set_2"
                in output_hdf5_file[
                    "output_data/dummy/dummy/convolution_results/"
                ].keys()
            )

            #
            indices_1 = output_hdf5_file[
                "output_data/dummy/dummy/convolution_results/set_1/0.5 Gyr/sampled_indices"
            ][()]
            self.assertTrue(len(indices_1) == 2502)

            #
            indices_2 = output_hdf5_file[
                "output_data/dummy/dummy/convolution_results/set_2/0.5 Gyr/sampled_indices"
            ][()]
            self.assertTrue(len(indices_2) == 1139)

    def test_postprocessing_multiple_dictionaries_without_name(self):

        #
        self.convolution_config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "dummy",
                "output_data_name": "dummy",
                "convolution_type": "sample",
                "data_column_dict": {
                    # required
                    "normalized_yield": "normalized_yield",
                    "delay_time": {"column_name": "time", "unit": u.Myr},
                },
                "post_convolution_function": postprocessing_multiple_dicts_without_name,
                "filter_future_events": False,
            },
        ]

        # Check if ValueError is raised
        with self.assertRaises(ValueError):
            # convolve
            convolve(config=self.convolution_config)


if __name__ == "__main__":
    unittest.main()
