"""
Some general functions related to the convolution codebase.

Mostly unsorted, likely better placed in together with related functionality.
"""

import functools
import json
import logging
import os
import shutil
import tempfile
import warnings
from inspect import isfunction

import astropy.units as u
import h5py
import numpy as np
import pandas as pd
import psutil
from astropy.cosmology import Planck13 as cosmo  # Planck 2013
from scipy import interpolate

logger = logging.getLogger(__name__)

dimensionless_unit = u.m / u.m


def maybe_strip_scaled_dimensionless(q):
    if not isinstance(q, u.Quantity):
        return q
    try:
        scale = q.to_value(u.dimensionless_unscaled)
        warnings.warn(
            f"Detected dimensionless-but-scaled quantity: {q.unit}. Stripping units."
        )
        return scale
    except u.UnitConversionError:
        return q


def print_hdf5_structure(f, subkey=None, detailed=True):
    if detailed:

        def _print_tree(name, obj):  # DH0001
            if isinstance(obj, h5py.Dataset):
                print(f"{name}: Dataset, shape={obj.shape}, dtype={obj.dtype}")
            elif isinstance(obj, h5py.Group):
                print(f"{name}: Group")

    else:

        def _print_tree(name, obj):  # DH0001
            print(subkey + "/" + name)

    if subkey is not None:
        f[subkey].visititems(_print_tree)
    else:
        f.visititems(_print_tree)


def generate_data_dict(config, convolution_instruction):
    """
    Function to generate the data dict.
    """

    # on the fly sampling generates its own data
    if "convolution_type" in convolution_instruction:
        if convolution_instruction["convolution_type"] == "on-the-fly":
            return config, {}, convolution_instruction

    #
    config["logger"].debug(
        "Generating data_dict using the extractor function {}".format(
            extract_data.__name__,
        )
    )

    # otherwise extract
    config, data_dict, convolution_instruction = extract_data(
        config=config, convolution_instruction=convolution_instruction
    )

    return config, data_dict, convolution_instruction


def extract_data(config, convolution_instruction):
    """
    Function to extract the data from the correct table and store the information in the correct column.

    Only extracts what is required by the data column dict
    """

    #
    data_dict = {}

    if convolution_instruction["chunked_readout"]:

        chunk_number = convolution_instruction["chunk_number"]
        chunksize = convolution_instruction["chunk_size"]  # Number of rows per chunk

        # Calculate row range
        start = chunk_number * chunksize
        stop = start + chunksize

        #
        df = pd.read_hdf(
            config["output_filename"],
            "/input_data/{}".format(convolution_instruction["input_data_name"]),
            start=start,
            stop=stop,
        )
    else:
        #
        df = pd.read_hdf(
            config["output_filename"],
            "/input_data/{}".format(convolution_instruction["input_data_name"]),
        )

    data_column_dict = convolution_instruction["data_column_dict"]

    # add all the columns to the data dictionary. This automatically handles the correct additional columns for the extra weights function
    for column in data_column_dict.keys():
        config["logger"].debug(
            "Extracting {} as the {} data".format(data_column_dict[column], column)
        )

        # if its a string we just assume its the column name
        if isinstance(data_column_dict[column], str):
            data_dict[column] = df[data_column_dict[column]].to_numpy()

            #################
            # Handle unit for delay-time
            if column == "delay_time":
                data_dict[column] = (
                    data_dict[column] * config["delay_time_default_unit"]
                )

        elif isinstance(data_column_dict[column], dict):
            if "column_name" not in data_column_dict[column]:
                raise ValueError(
                    "Please provide the input-data column name through the 'column_name' key."
                )

            # extract data with the explicit column name entry
            data = df[data_column_dict[column]["column_name"]].to_numpy()

            #################
            # Handle conversion
            data = handle_custom_scaling_or_conversion(
                config=config,
                data_layer_or_column_dict_entry=data_column_dict[column],
                value=data,
            )

            # Store
            data_dict[column] = data

            #################
            # Handle unit for delay-time
            # TODO: this should just take whatever unit is provided
            if column == "delay_time":
                if "unit" in data_column_dict[column].keys():
                    unit = data_column_dict[column]["unit"]
                else:
                    unit = config["delay_time_default_unit"]

                #
                data_dict[column] = data_dict[column] * unit
        else:
            raise ValueError("input type not supported.")

    ##########
    # If we have binned data we should addd the delay time bin indices to the
    if convolution_instruction["contains_binned_data"]:
        data_dict["delay_time_data_bin_index"] = (
            np.digitize(
                data_dict["delay_time"].to(u.yr),
                convolution_instruction["delay_time_data_bin_info_dict"][
                    "delay_time_data_bin_edges"
                ].to(u.yr),
            )
            - 1
        )

    #
    return config, data_dict, convolution_instruction


def sample_around_bin_center(bin_edges, values):
    """
    Basic function to handle sampling around bincenter given bin edges and values.

    Note: this does not handle values that fall outside of the bins well
    """

    bin_widths = np.diff(bin_edges)

    indices = np.digitize(values, bin_edges) - 1

    # get random values and scale
    random_arr = np.random.random(indices.shape) - 0.5
    random_arr = random_arr * bin_widths[indices]

    # Add to values
    sampled_values = values + random_arr

    return sampled_values


def create_job_dict(
    config, sfr_dict, data_dict, convolution_instruction, time_bin_info_dict, bin_number
):
    """
    Function to create the job dict
    """

    # Set up job dict
    job_dict = {
        "job_number": bin_number,
        "time_bin_info_dict": time_bin_info_dict,
        "sfr_dict": sfr_dict,
        "convolution_instruction": convolution_instruction,
        "data_dict": data_dict,
        "output_dir": get_tmp_dir(
            config=config,
            convolution_instruction=convolution_instruction,
            sfr_dict=sfr_dict,
        ),
    }

    return job_dict


def create_time_bin_info_dict(
    config,
    convolution_instruction,
    bin_number,
    bin_center,
    bin_edge_lower,
    bin_size,
    bin_type,
):
    """
    Function to set up the time bin info dict
    """

    time_bin_info_dict = {
        "bin_number": bin_number,
        "bin_center": bin_center,
        "bin_edge_lower": bin_edge_lower,
        "bin_size": bin_size,
        "bin_type": bin_type,
        "time_type": config["time_type"],
        "reverse_bin_order": convolution_instruction["reverse_convolution"],
        "convolution_direction": convolution_instruction["convolution_direction"],
    }

    return time_bin_info_dict


def get_physical_dimensions(unit):
    """Return the physical dimensions of a unit in sorted [M][L][T] notation using SI base units."""
    # Decompose into SI base units
    decomposed = unit.decompose(bases=u.si.bases)

    # Extract base units and their powers
    si_dimensions = {
        str(base.physical_type): power
        for base, power in zip(decomposed.bases, decomposed.powers)
    }

    # Define standard SI dimension notation
    notation_map = {
        "mass": "M",
        "length": "L",
        "time": "T",
        "current": "I",
        "temperature": "Θ",
        "amount of substance": "N",
        "luminous intensity": "J",
    }

    # Convert to notation, ensuring unknown types don't appear
    sorted_notation = [
        f"[{notation_map[ptype]}^{power}]" if power != 1 else f"[{notation_map[ptype]}]"
        for ptype, power in sorted(si_dimensions.items())  # Sorting works correctly now
        if ptype in notation_map  # Ignore unknown types
    ]

    return "".join(sorted_notation) or "[dimensionless]"


def generate_boilerplate_outputfile(outputfile_name):
    """
    Function to generate a boilerplate output file structure so the user does not have to worry about including the correct groups.
    """

    # create file
    output_hdf5_file = h5py.File(outputfile_name, "w")

    # Create groups main
    output_hdf5_file.create_group("input_data")
    output_hdf5_file.create_group("config")

    # close
    output_hdf5_file.close()


def extract_unit_dict(output_hdf5_file, key):
    """
    FUnction to extract the unit dict from the hdf5 file
    """

    # convert units
    unit_dict = json.loads(output_hdf5_file[key].attrs["units"])
    unit_dict = {key: u.Unit(val) for key, val in unit_dict.items()}

    return unit_dict


def get_username():
    """
    Function to get the username of the user that spawned the current process
    """

    return psutil.Process().username()


class JsonCustomEncoder(json.JSONEncoder):
    """Support for data types that JSON default encoder
    does not do.

    This includes:

        * Numpy array or number
        * Complex number
        * Set
        * Bytes
        * astropy.UnitBase
        * astropy.Quantity

    Examples
    --------
    >>> import json
    >>> import numpy as np
    >>> from astropy.utils.misc import JsonCustomEncoder
    >>> json.dumps(np.arange(3), cls=JsonCustomEncoder)
    '[0, 1, 2]'

    copied from astropy and extended

    """

    def default(self, obj):  # DH0001
        import numpy as np
        from astropy import units as u

        if isinstance(obj, u.Quantity):
            return dict(value=obj.value, unit=obj.unit.to_string())
        if isinstance(obj, (np.number, np.ndarray)):
            return obj.tolist()
        elif isinstance(obj, complex):
            return [obj.real, obj.imag]
        elif isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, bytes):  # pragma: py3
            return obj.decode()
        elif isinstance(obj, interpolate.interp1d):
            return str(obj)
        elif isinstance(obj, (u.UnitBase, u.FunctionUnitBase)):
            if obj == u.dimensionless_unscaled:
                obj = "dimensionless_unit"
                return str(obj)
            else:
                return obj.to_string()
        elif isinstance(obj, type(logger)):
            return str(obj)
        elif isinstance(obj, type(cosmo)):
            return str(obj)
        elif isinstance(obj, type(cosmo)):
            return str(obj)
        elif isfunction(obj):
            return str(obj)

        return json.JSONEncoder.default(self, obj)


####
# General functions
def custom_json_serializer(obj):  # DH0001
    """
    Custom serialiser for binary_c to use when functions are present in the dictionary
    that we want to export.

    Function objects will be turned into str representations of themselves

    Args:
        obj: The object that might not be serialisable

    Returns:
        Either string representation of object if the object is a function, or the object itself
    """

    if isinstance(obj, u.Quantity):
        return obj.value
    elif isinstance(obj, interpolate.interp1d):
        return str(obj)
    # elif isinstance(o, t)
    return obj


def verbose_print(  # DH0001
    message: str, verbosity: int, minimal_verbosity: int
) -> None:
    """
    Function that decides whether to print a message based on the current verbosity
    and its minimum verbosity

    if verbosity is equal or higher than the minimum, then we print

    Args:
        message: message to print
        verbosity: current verbosity level
        minimal_verbosity: threshold verbosity above which to print
    """

    if verbosity >= minimal_verbosity:
        print(message)


def vb(message, verbosity, minimal_verbosity):  # DH0001
    """
    Shorthand for verbose_print
    """

    verbose_print(message, verbosity, minimal_verbosity)


def calculate_bincenters(array, convert="linear"):
    """
    Function to calculate bincenters

    TODO: allow other conversions
    """

    if convert == "linear":
        bincenters = (array[1:] + array[:-1]) / 2
    else:
        raise ValueError(f"convert choice {convert} is unknown")

    return bincenters


def calculate_bin_edges(arr):
    """
    Function to calculate the edge values given a bunch of centers
    """

    #
    diff = np.diff(arr)
    edge_values = (arr[1:] + arr[:-1]) / 2

    #
    edge_values = np.insert(
        edge_values,
        0,
        edge_values[0] - diff[0],
        axis=0,
    )

    #
    edge_values = np.insert(
        edge_values,
        edge_values.shape[0],
        edge_values[-1] + diff[-1],
        axis=0,
    )

    return edge_values


def pad_function(array, left_val, right_val, relative_to_edge_val, axis=0):
    """
    Function to pad an array
    """

    # copy
    padded_array = array[:]

    # check if there are units involved
    try:
        unit = padded_array.unit

        left_val = left_val * unit
        right_val = right_val * unit
    except AttributeError:
        pass

    #
    if relative_to_edge_val:
        #
        padded_array = np.insert(
            padded_array,
            0,
            padded_array[axis] + left_val,
            axis=axis,
        )

        #
        padded_array = np.insert(
            padded_array,
            padded_array.shape[axis],
            padded_array[-1] + right_val,
            axis=axis,
        )
    else:

        #
        padded_array = np.insert(
            padded_array,
            0,
            left_val,
            axis=axis,
        )

        #
        padded_array = np.insert(
            padded_array,
            padded_array.shape[axis],
            right_val,
            axis=axis,
        )

    return padded_array


def generate_group_name(convolution_instruction, sfr_dict):
    """
    Function to generate the group name. Also provides layers
    """

    #
    elements = []

    if sfr_dict is None:
        sfr_dict = {}

    #
    if sfr_dict.get("name", None) is not None:
        elements.append(sfr_dict["name"])

    # store input name
    elements.append(convolution_instruction["input_data_name"])

    # store chunkname
    if (
        convolution_instruction["chunked_readout"]
        and "chunk_number" in convolution_instruction
    ):
        elements.append(str(convolution_instruction["chunk_number"]))

    # store output name
    elements.append(convolution_instruction["output_data_name"])

    # construct groupname
    groupname = "/".join(elements)

    return groupname, elements


def get_tmp_dir(config, convolution_instruction, sfr_dict=None):
    """
    Function to get tmp dir
    """

    #
    groupname, _ = generate_group_name(
        convolution_instruction=convolution_instruction, sfr_dict=sfr_dict
    )

    #
    tmp_dir = os.path.join(config["tmp_dir"], groupname)

    return tmp_dir


def handle_custom_scaling_or_conversion(config, data_layer_or_column_dict_entry, value):
    """
    Function that handles multiplying the key of the ensemble with some value or with some function
    """

    ###########
    # Handle logic of multiple steps
    if ("conversion_factor" in data_layer_or_column_dict_entry.keys()) and (
        "conversion_function" in data_layer_or_column_dict_entry.keys()
    ):
        raise ValueError(
            "We currently do not support both a conversion factor and a conversion function"
        )

    ###########
    # convert data by a function
    if "conversion_factor" in data_layer_or_column_dict_entry.keys():
        value = value * data_layer_or_column_dict_entry["conversion_factor"]

        #
        config["logger"].debug(
            "Applying conversion factor {} on data column {}".format(
                data_layer_or_column_dict_entry["conversion_factor"],
                data_layer_or_column_dict_entry,
            )
        )

    ###########
    # convert data by a function
    if "conversion_function" in data_layer_or_column_dict_entry.keys():
        value = data_layer_or_column_dict_entry["conversion_function"](value)

        #
        config["logger"].debug(
            "Applying conversion function {} on data column {}".format(
                data_layer_or_column_dict_entry["conversion_function"].__name__,
                data_layer_or_column_dict_entry,
            )
        )

    return value


def temp_dir(*child_dirs: str, clean_path=False) -> str:
    """
    Function to create directory within the TMP directory of the file system, starting with `/<TMP>/binary_c_python-<username>`

    Makes use of os.makedirs exist_ok which requires python 3.2+

    Args:
        *child_dirs: str input where each next input will be a child of the previous full_path. e.g. ``temp_dir('tests', 'grid')`` will become ``'/tmp/binary_c_python-<username>/tests/grid'``
        *clean_path (optional): Boolean to make sure that the directory is cleaned if it exists
    Returns:
        the path of a sub directory called binary_c_python in the TMP of the file system
    """

    tmp_dir = tempfile.gettempdir()
    username = get_username()
    full_path = os.path.join(tmp_dir, "sspc-{}".format(username))

    # loop over the other paths if there are any:
    if child_dirs:
        for extra_dir in child_dirs:
            full_path = os.path.join(full_path, extra_dir)

    # Check if we need to clean the path
    if clean_path and os.path.isdir(full_path):
        shutil.rmtree(full_path)

    #
    os.makedirs(full_path, exist_ok=True)

    return full_path


def check_required(config, required_list):
    """
    Function to check if the keys in the required_list are present in the convolution_instruction dict
    """

    for key in required_list:
        if key not in config.keys():
            raise ValueError(
                "{} is required in the convolution_instruction".format(key)
            )


def is_time_unit(parameter):
    """
    Function to check if a parameter has time-units
    """

    try:
        parameter.to(u.yr)
        return True
    except u.core.UnitConversionError:
        return False
    except AttributeError:
        return False


def is_mass_unit(parameter):
    """
    Function to check if a parameter has time-units
    """

    try:
        parameter.to(u.kg)
        return True
    except u.core.UnitConversionError:
        return False
    except AttributeError:
        return False


def has_unit(parameter, fail_on_dimensionless=True):
    """
    Function to check if a parameter has any unit assigned to it
    """

    try:
        unit = parameter.unit

        if fail_on_dimensionless:
            dimensionless_unit = u.m / u.m
            if unit == dimensionless_unit:
                return False
        return True
    except:
        return False


has_unit_dimensionless_okay = functools.partial(has_unit, fail_on_dimensionless=False)


def get_normalized_yield_unit(config, convolution_instruction):
    """
    Function to get the normalized yield unit either from config or from convolution_instruction
    """

    #
    normalized_yield_unit = config["default_normalized_yield_unit"]

    if "normalized_yield" not in convolution_instruction["data_column_dict"]:
        raise ValueError(
            "'normalized_yield' should be provided in the 'data_column_dict'"
        )

    if isinstance(
        convolution_instruction["data_column_dict"]["normalized_yield"], dict
    ):
        if "unit" in convolution_instruction["data_column_dict"]["normalized_yield"]:
            normalized_yield_unit = convolution_instruction["data_column_dict"][
                "normalized_yield"
            ]["unit"]

    return normalized_yield_unit
