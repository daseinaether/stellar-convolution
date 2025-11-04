"""
Testcases for convolution_with_events file

TODO: test things with redshift
TODO: test things with multiply SFR histories
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
    generate_boilerplate_outputfile,
    temp_dir,
)

TMP_DIR = temp_dir(
    "tests", "tests_convolution", "test_convolution_with_events", clean_path=True
)


class test_convolution_with_events(unittest.TestCase):
    """ """

    def test_convolution_with_events_with_lookback_time(self):
        """ """

        #
        output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5_sfr_only.h5")
        generate_boilerplate_outputfile(output_hdf5_filename)

        ##############
        # SET UP DATA
        dummy_data = {
            "delay_time": np.array([0, 1, 2, 3]),
            "probability": np.array([1, 2, 3, 4]),
        }
        dummy_df = pd.DataFrame.from_records(dummy_data)

        ##############
        # Store data in pandas
        dummy_df.to_hdf(output_hdf5_filename, key="input_data/{}".format("dummy"))

        #
        convolution_config = copy.copy(default_convolution_config)
        convolution_config["logger"].setLevel(logging.CRITICAL)

        # Set up SFR
        convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1, 1, 1])
            * u.Msun
            / u.yr
            / u.Gpc**3,
        }

        # set up convolution bins
        convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1, 2, 3, 4]) * u.yr
        )

        # lookback time convolution only
        convolution_config["time_type"] = "lookback_time"

        #
        convolution_config["output_filename"] = output_hdf5_filename

        convolution_config["redshift_interpolator_data_output_filename"] = os.path.join(
            TMP_DIR, "interpolator_dict.p"
        )

        #
        convolution_config["convolution_instructions"] = [
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

        #
        convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")

        #
        convolve(config=convolution_config)

        #
        with h5py.File(output_hdf5_filename, "r") as output_hdf5_file:

            #
            arr_ = output_hdf5_file[
                "output_data/dummy/dummy/convolution_results/0.5 yr/yield"
            ][()]
            self.assertTrue(np.array_equal(arr_, np.array([1, 2, 3, 4])))

            #
            arr_ = output_hdf5_file[
                "output_data/dummy/dummy/convolution_results/1.5 yr/yield"
            ][()]
            self.assertTrue(np.array_equal(arr_, np.array([1, 2, 3, 4])))

            #
            arr_ = output_hdf5_file[
                "output_data/dummy/dummy/convolution_results/2.5 yr/yield"
            ][()]
            self.assertTrue(np.array_equal(arr_, np.array([1, 2, 3, 0])))

            #
            arr_ = output_hdf5_file[
                "output_data/dummy/dummy/convolution_results/3.5 yr/yield"
            ][()]
            self.assertTrue(np.array_equal(arr_, np.array([1, 2, 0, 0])))

            #
            self.assertTrue("SFR_info" in output_hdf5_file["output_data"].attrs.keys())


if __name__ == "__main__":
    unittest.main()
