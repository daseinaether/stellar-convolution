"""
File containing the default values and the validations for the configuration of the convolution
"""

import logging
import os
from typing import Callable

import astropy.units as u
import numpy as np
import voluptuous as vol
from astropy.cosmology import Planck13 as cosmo  # Planck 2013

ALLOWED_NUMERICAL_TYPES = (int, float, complex, np.number)
dimensionless_unit = u.m / u.m

##########
# Logger configuration
logger = logging.getLogger(__name__)
FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s ] %(asctime)s: %(message)s"
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.CRITICAL)


#################
# Validation routines
def list_of_dicts_validation(value):
    if isinstance(value, list):
        for el in value:
            if not isinstance(el, dict):
                raise ValueError(
                    "All entries in the list should be dictionary-type objects"
                )
    else:
        raise ValueError("Input has to either be a list or a dict")


def dict_or_list_of_dicts_validation(value):
    if not isinstance(value, (list, dict)):
        raise ValueError("Input has to either be a list or a dict")
    if isinstance(value, list):
        list_of_dicts_validation(value)


def unit_validation(value):
    # allow 'astropy.units.core.Unit, 'astropy.units.quantity.Quantity', 'astropy.units.core.CompositeUnit'
    if not isinstance(value, (type(u.yr), type(1 / u.Msun), type(u.Msun**-1))):
        raise ValueError("Input has to be a astropy-unit object")


def logger_validation(value):
    if not isinstance(value, type(logger)):
        raise ValueError("Input has to be a logging-type object")


def callable_validation(value):
    if not isinstance(value, Callable):
        raise ValueError("Input has to be a callable")


def callable_or_none_validation(value):
    if value is not None:
        if not isinstance(value, Callable):
            raise ValueError("Input has to be a callable")


def array_validation(value):
    if not isinstance(value, type(np.array([]))):
        raise ValueError("Input has to be a numpy array")


def existing_path_validation(value):
    if isinstance(value, str):
        if not os.path.isfile(value):
            raise ValueError("File doesnt exist")
    else:
        raise ValueError("Please provide a string-based input")


#
boolean_int_validation = vol.All(vol.Range(max=1), vol.Boolean())
float_or_int = vol.Or(float, int)

############################
#
default_convolution_config_dict = {
    ###################
    # Convolution configuration
    "time_type": {
        "value": "lookback_time",
        "description": "Time-type used in convolution. Can be either 'redshift' or 'lookback_time'",
        "validation": vol.All(
            str,
            vol.In(["redshift", "lookback_time"]),
        ),
    },
    # ###################
    # # Starformation related
    "SFR_info": {
        "value": {},
        "description": "dictionary containing the starformation rate info. Can also be a list of dictionaries.",
        "validation": dict_or_list_of_dicts_validation,
    },
    # Convolution time bins
    "convolution_lookback_time_bin_edges": {
        "value": None,
        "description": "Lookback-time bin-edges used in convolution.",  # TODO: update this if we do things with units
        "validation": array_validation,
    },
    "convolution_redshift_bin_edges": {
        "value": None,
        "description": "Redshift bin-edges used in convolution.",  # TODO: update this if we do things with units
        "validation": array_validation,
    },
    ###################
    # Redshift interpolator settings
    "redshift_interpolator_force_rebuild": {
        "value": False,
        "description": "Whether to force rebuild the redshift interpolator.",  # TODO: expand explanation
        "validation": boolean_int_validation,
    },
    "redshift_interpolator_rebuild_when_settings_mismatch": {
        "value": True,
        "description": "Whether to rebuild the redshift interpolator when the config of the existing one don't match with the current config.",  # TODO: expand explanation
        "validation": boolean_int_validation,
    },
    "redshift_interpolator_stepsize": {
        "value": 0.001,
        "description": "Stepsize for the redshift interpolation.",  # TODO: expand explanation. Also consider if this is the best way
        "validation": float_or_int,
    },
    "redshift_interpolator_min_redshift_if_log": {
        "value": 1e-5,
        "description": "Minimum redshift for the redshift interpolator if using log spacing.",  # TODO: expand explanation. Also consider if this is the best way
        "validation": float_or_int,
    },
    "redshift_interpolator_min_redshift": {
        "value": 0,
        "description": "Minimum redshift for the redshift interpolator",  # TODO: expand explanation. Also consider if this is the best way
        "validation": float_or_int,
    },
    "redshift_interpolator_max_redshift": {
        "value": 50,
        "description": "Minimum redshift for the redshift interpolator",  # TODO: expand explanation. Also consider if this is the best way
        "validation": float_or_int,
    },
    "redshift_interpolator_use_log": {
        "value": True,
        "description": "Whether to interpolate in log redshift.",  # TODO: expand explanation
        "validation": boolean_int_validation,
    },
    "redshift_interpolator_data_output_filename": {
        # "value": os.path.join(
        #     os.path.dirname(os.path.abspath(__file__)), "interpolator_data_dict.p"
        # ),
        "value": None,
        "description": "Filename for the redshift interpolator object.",  # TODO: expand explanation. Also consider if this is the best way
        "validation": str,
    },
    ###################
    # Multiprocessing settings
    "multiprocessing": {
        "value": True,
        "description": "Flag whether to enable multiprocessing. True for multiprocessing, which allows faster convolution but does not allow the use of the previous convolution results and the persistent data. False for sequential convolution, which is slower but previous convolution results and the persistent data is available here.",
        "validation": int,
    },
    "num_cores": {
        "value": 1,
        "description": "Number of cores to use to do the convolution",
        "validation": int,
    },
    "max_job_queue_size": {
        "value": 8,
        "description": "Max number of jobs in the multiprocessing queue for the convolution.",
        "validation": int,
    },
    ###################
    # custom convolution
    "custom_convolution_function": {
        "value": None,
        "description": "",  # TODO: expand
        "validation": callable_or_none_validation,
    },
    "custom_data_extraction_function": {
        "value": None,
        "description": "",  # TODO: expand
        "validation": callable_or_none_validation,
    },
    "include_custom_rates": {
        "value": False,
        "description": "Whether to include custom, user-specified, rates. See 'custom_rates_function'.",  # TODO: expand
        "validation": boolean_int_validation,
    },
    "custom_rates_function": {
        "value": None,
        "description": "Custom rate function used in the convolution.",  # TODO: expand
        "validation": callable_or_none_validation,
    },
    ###################
    # Logger
    "logger": {
        "value": logger,
        "description": "Logger object.",
        "validation": logger_validation,
    },
    ###################
    # Output processing settings
    "write_to_hdf5": {
        "value": True,
        "description": "Whether to write the pickle-files from the convolution back to the main hdf5 file",
        "validation": boolean_int_validation,
    },
    "remove_pickle_files": {
        "value": True,
        "description": "Flag whether to remove all the pickle files after writing them to the main hdf5 file",
        "validation": boolean_int_validation,
    },
    ###################
    # unsorted
    #
    "tmp_dir": {
        "value": "/tmp/sspc",
        "description": "Target directory for the tmp files.",  # TODO: expand explanation. Also consider if this is the best way
        "validation": str,
    },
    "cosmology": {
        "value": cosmo,
        "description": "Astropy cosmology used throughout the code. ",  # TODO: expand explanation
        # "validation": # TODO: add validation
    },
    "convolution_instructions": {
        "value": [{}],
        "description": "List of instructions for the convolution. ",  # TODO: expand explanation
        "validation": list_of_dicts_validation,  # TODO: lets also allow just 1 convolution instruction as a dictionary
        # "validation": # NOTE: validation handled with custom function
    },
    "output_filename": {
        "value": "",
        "description": "Full path to output hdf5 filename. This should point to a file that already contains the input data that will be used to do the convolution with.",
        "validation": str,
    },
    "check_convolution_config": {
        "value": True,
        "description": "Flag whether to validate the configuration dictionary before running the convolution code.",
        "validation": boolean_int_validation,
    },
    "delay_time_default_unit": {
        "value": u.yr,
        "description": "Default unit used for the delay-time data. NOTE: this can be overridden in data_dict column or layer entries.",
        "validation": unit_validation,
    },
    "default_normalized_yield_unit": {
        "value": 1.0 / u.Msun,
        "description": "Default unit used for the normalized-yield data.",
        "validation": unit_validation,
    },
}

# extract only values
default_convolution_config = {
    key: value["value"] for key, value in default_convolution_config_dict.items()
}

# extract only descriptions
default_convolution_config_descriptions = {
    key: value["description"] for key, value in default_convolution_config_dict.items()
}
