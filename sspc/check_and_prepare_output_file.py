"""
Function to copy the input file and
"""

import json

import h5py

from syntheticstellarpopconvolve.general_functions import JsonCustomEncoder


def check_and_prepare_output_file(config):
    """
    Function to prepare the output file, create some initial groups and store the configuration
    """

    #
    config["logger"].debug("Checking and preparing convolution output file")

    #
    output_file = h5py.File(config["output_filename"], "r")

    # check if there is data in the file
    if "input_data" not in output_file.keys():
        raise ValueError("Please provide a 'input_data' group in the output hdf5file.")

    # check if there is a config group and
    if "config" not in output_file.keys():
        raise ValueError("Please provide a 'config' group in the output hdf5file.")

    output_file.close()

    # Store config.
    with h5py.File(config["output_filename"], "a") as output_hdf5file:

        if "convolution" not in output_hdf5file["config"].keys():

            # Store convolution configuration in
            output_hdf5file["config"].create_dataset(
                "convolution", data=json.dumps(config, cls=JsonCustomEncoder)
            )

        else:
            config["logger"].warning(
                "tried to store config in output file, but was already present. "
            )
