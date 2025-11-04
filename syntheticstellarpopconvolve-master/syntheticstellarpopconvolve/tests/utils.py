"""
Util functions for the tests
"""

import copy
import os

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve import default_convolution_config
from syntheticstellarpopconvolve.general_functions import (
    generate_boilerplate_outputfile,
)


class Boilerplate:

    def __init__(self):
        """ """

    def setup(self, name, tmp_dir, add_population_settings, sfr_unit):
        #
        output_hdf5_filename = os.path.join(tmp_dir, "output_hdf5_{}.h5".format(name))
        generate_boilerplate_outputfile(output_hdf5_filename)

        #
        self.convolution_config = copy.copy(default_convolution_config)

        # Set up SFR
        self.convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.Gyr,
            "starformation_rate_array": np.array([2, 1, 1, 1, 1]) * sfr_unit,
        }

        # set up convolution bins
        self.convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1, 2, 3, 4]) * u.Gyr
        )

        # lookback time convolution only
        self.convolution_config["time_type"] = "lookback_time"

        #
        self.convolution_config["output_filename"] = output_hdf5_filename

        self.convolution_config["redshift_interpolator_data_output_filename"] = (
            os.path.join(tmp_dir, "interpolator_dict.p")
        )

        #
        self.convolution_config["tmp_dir"] = os.path.join(
            tmp_dir, "tmp_{}".format(name)
        )
