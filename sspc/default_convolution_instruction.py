"""
File containing the default convolution-instruction.
"""

import astropy.units as u
import numpy as np
import voluptuous as vol

from syntheticstellarpopconvolve.default_convolution_config import (
    boolean_int_validation,
    callable_or_none_validation,
)

ALLOWED_NUMERICAL_TYPES = (int, float, complex, np.number)
dimensionless_unit = u.m / u.m

# any relevant option needs to be in here
default_convolution_instruction_dict = {
    ########################
    # Unsorted
    "chunked_readout": {
        "value": False,
        "description": "Flag to read the input data in chunks. See `chunk_size`.",
        "validation": boolean_int_validation,
    },
    "chunk_size": {
        "value": 0,
        "description": "chunk size for the data readout.",
        "validation": int,
    },
    "chunk_total": {
        "value": 0,
        "description": "Total number of chunks to be considered. Should be an integer rounded up calculated as ceil(<total entries in input dataframe>/<chunk size>).",
        "validation": int,
    },
    "multiply_by_convolution_time_binsize": {
        "value": False,
        "description": "Flag to multiply the convolution results by the convolution time-bin size. Not supported when time_type=='redshift'.",
        "validation": boolean_int_validation,
    },
    "multiply_by_sfr_time_binsize": {
        "value": False,
        "description": "Flag to multiply the convolution results by the starformation rate time bin size. Not supported when time_type=='redshift'.",
        "validation": boolean_int_validation,
    },
    "convolution_direction": {
        "value": "backward",
        "description": "Choice of convolution direction. 'backward' convolves the data such that every event occurs at the current time by looking what the starformation rate is for each given delay time. 'forward' convolution generates each system at the same time and looks at when events happen afterwards based on their delay times. Note: neither option is supported in all choices of `convolution_type` and `contains_binned_data`.",
        "validation": vol.All(
            str,
            vol.In(["backward", "forward"]),
        ),
    },
    "reverse_convolution": {
        "value": False,
        "description": "Flag to reverse the convolution direction. If True, we start with the bin furthest back in time and work to. Useful in combination with `convolution_config['multiprocessing']=False`, and the `previous_convolution_results` and `persistant_data` objects.",
        "validation": boolean_int_validation,
    },
    "contains_binned_data": {
        "value": False,
        "description": "Flag to indicate whether the input data is binned (in time). If so, the user should provide additional information",
        "validation": boolean_int_validation,
    },
    "data_time_bin_info": {
        "value": {},
        "description": "Dictionary containing data time-bin information when convolving binned data.",
        "validation": dict,
    },
    ########################
    #
    "input_data_name": {
        "value": "input_data",
        "description": "Name of to the current input dataset. Will be used to extract the data from the provided input-hdf5 file (expected in /input_data/<input data name>/), and will be used in the output-data path. ",
        "validation": str,
    },
    "output_data_name": {
        "value": "output_data",
        "description": "Name assigned to the current output dataset. Will be used in the output-data path. Can be useful when running a convolution with the same input-data but e.g. with a different post-convolution function.",
        "validation": str,
    },
    "convolution_type": {
        "value": "integrate",
        "description": "Method of convolution. The three choices are as follows. 'integrate': Convolution by integration uses backward convolution to multiply the normalized_yield of the systems with the starformation rate at the time the system would be born, given the delay time and the target event time. This is particularly useful when you are just interested in the (transient) event. 'sample': Convolution by sampling uses forward convolution to 'sample' systems according to the yield () and assigns an event-time based on the delay time and the birth time. This is particularly useful if you want to post-process systems after they are born/the event occurs, like integrating the orbit of double compact object forward in time due to gravitational wave radiation (see LISA project example in `examples/notebook_usecases`). 'on-the-fly': Convolution by simulating the systems on-the-fly. Requires a method that uses the total mass in star formation, and optionally the metallicity distribution, combined with a population synthesis code, to evolve systems on the fly.",
        "validation": vol.All(
            str,
            vol.In(["integrate", "sample", "on-the-fly"]),
        ),
    },
    "data_column_dict": {
        "value": {},
        "description": "Dictionary containing the mapping between the names of the columns in the pandas dataframe of the input data and the names as they are used in the convolution framework, as `{<framework_data_name>: <pandas_column_name>}`. Mappings for TODO are required. Any extra columns that are provided are accessible by the post-convolution function. Entries in this dictionary can be either names or dictionaries themselves, allowing more advanced functionality. See examples/notebook_convolution_advanced",
        "validation": dict,  # TODO: make more specific: dict with string values or dict with dict values
    },
    "post_convolution_function": {
        "value": None,
        "description": "Function that performs post-convolution operations on the convolved data, like applying detection probability weights, further integration of systems or just general filtering of the data. Different `convolution_type`s allows for different modifications of the data, in that convolution of `event-based` data by 'sampling' allows TODO. The arguments of this function should be chosen from: 'config', 'time_value', 'convolution_instruction', 'data_dict' and the contents of 'post_convolution_function_extra_parameters'. For more explanation about this function see the convolution notebook.",
        "validation": callable_or_none_validation,
    },
    "post_convolution_function_extra_parameters": {
        "value": {},
        "description": "Dictionary containing additional arguments that can be accessed by the extra_weights_function. Note: using this can often be avoided by using functools.partial to fix certain input parameters.",
        "validation": dict,
    },
    # "assign_formation_lookback_time": {
    #     "value": True,
    #     "description": "Flag to indicate to assign formation-lookback times during convolution-by-sampling. If false, assigning event-lookback times and filtering future events will not be possible, and `assign_event_lookback_time` and `filter_future_events` are ignored. Based on the current convolution bin and delay time of the events. Not recommended when using binned data.",
    #     "validation": boolean_int_validation,
    # },
    # "assign_event_lookback_time": {
    #     "value": True,
    #     "description": "Flag to indicate to assign event-lookback times during convolution-by-sampling. If false, filtering future events will not be possible, and `filter_future_events` is ignored. Based on the current convolution bin and delay time of the events. Not recommended when using binned data.",
    #     "validation": boolean_int_validation,
    # },
    # "filter_future_events": {
    #     "value": True,
    #     "description": "Flag to control filtering out future events during convolution-by-sampling. See also `assign_event_lookback_time`. Not recommended when using binned data.",
    #     "validation": boolean_int_validation,
    # },
}


# extract only values
default_convolution_instruction = {
    key: value["value"] for key, value in default_convolution_instruction_dict.items()
}

# extract only descriptions
default_convolution_instruction_descriptions = {
    key: value["description"]
    for key, value in default_convolution_instruction_dict.items()
}
