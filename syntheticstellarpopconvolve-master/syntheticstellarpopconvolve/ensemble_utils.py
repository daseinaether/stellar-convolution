"""
File containing utility functions related to ensemble-based data
"""

# import bz2
# import collections
# import copy
# import gzip
# import sys
import time
from collections import OrderedDict

# import astropy.units as u
# import msgpack
import numpy as np
import pandas as pd

from syntheticstellarpopconvolve.default_convolution_config import (
    ALLOWED_NUMERICAL_TYPES,
)

# import simplejson
# from halo import Halo


# from syntheticstellarpopconvolve.general_functions import (
#     calculate_bin_edges,
#     handle_custom_scaling_or_conversion,
#     has_unit_dimensionless_okay,
# )

# def ensemble_compression(filename):  # DH0001
#     """
#     Return the compression type of the ensemble file, based on its filename extension.
#     """

#     if filename.endswith(".bz2"):
#         return "bzip2"
#     if filename.endswith(".gz"):
#         return "gzip"
#     return None


# def open_ensemble(filename, encoding="utf-8"):  # DH0001
#     """
#     Function to open an ensemble at filename for reading and decompression if required.
#     """

#     compression = ensemble_compression(filename)
#     if ensemble_file_type(filename) == "msgpack":
#         flags = "rb"
#     else:
#         flags = "rt"
#     if compression == "bzip2":
#         file_object = bz2.open(filename, flags, encoding=encoding)
#     elif compression == "gzip":
#         file_object = gzip.open(filename, flags, encoding=encoding)
#     else:
#         file_object = open(filename, flags, encoding=encoding)
#     return file_object


# def keys_to_floats(input_dict: dict) -> dict:  # DH0001
#     """
#     Function to convert all the keys of the dictionary to float to float

#     we need to convert keys to floats:
#         this is ~ a factor 10 faster than David's ``recursive_change_key_to_float`` routine, probably because this version only does the float conversion, nothing else.

#     Args:
#         input_dict: dict of which we want to turn all the keys to float types if possible

#     Returns:
#         new_dict: dict of which the keys have been turned to float types where possible
#     """

#     # this adopts the type correctly *and* is fast
#     new_dict = type(input_dict)()

#     for k, v in input_dict.items():
#         # convert key to a float, if we can
#         # otherwise leave as is
#         try:
#             newkey = float(k)
#         except ValueError:
#             newkey = k

#         # act on value(s)
#         if isinstance(v, list):
#             # list data
#             new_dict[newkey] = [
#                 (
#                     keys_to_floats(item)
#                     if isinstance(item, collections.abc.Mapping)
#                     else item
#                 )
#                 for item in v
#             ]
#         elif isinstance(v, collections.abc.Mapping):
#             # dict, ordereddict, etc. data
#             new_dict[newkey] = keys_to_floats(v)
#         else:
#             # assume all other data are scalars
#             new_dict[newkey] = v

#     return new_dict


# def ensemble_file_type(filename):  # DH0001
#     """
#     Returns the file type of an ensemble file.
#     """

#     if ".json" in filename:
#         filetype = "JSON"
#     elif ".msgpack" in filename:
#         filetype = "msgpack"
#     else:
#         filetype = None
#     return filetype


# def load_ensemble(  # DH0001
#     filename,
#     convert_float_keys=True,
#     select_keys=None,
#     timing=False,
#     flush=False,
#     quiet=False,
# ):  # DH0001
#     """
#     Function to load an ensemeble file, even if it is compressed,
#     and return its contents to as a Python dictionary.

#     Args:
#         convert_float_keys : if True, converts strings to floats.
#         select_keys : a list of keys to be selected from the ensemble.
#     """

#     # open the file

#     # load with some info to the terminal
#     if not quiet:
#         print("Loading JSON...", flush=flush)

#     # open the ensemble and get the file type
#     file_object = open_ensemble(filename)
#     filetype = ensemble_file_type(filename)

#     if not filetype or not file_object:
#         print(
#             "Unknown filetype : your ensemble should be saved either as JSON or msgpack data.",
#             flush=flush,
#         )
#         sys.exit()

#     if quiet:
#         tstart = time.time()
#         if filetype == "JSON":
#             data = simplejson.load(file_object)
#             file_object.close()
#         elif filetype == "msgpack":
#             data = msgpack.load(file_object, object_hook=_hook)  # noqa: F821
#             file_object.close()
#         if timing:
#             print(
#                 "\n\nTook {} s to load the data\n\n".format(time.time() - tstart),
#                 flush=True,
#             )
#     else:
#         with Halo(text="Loading", interval=250, spinner="moon", color="yellow"):
#             tstart = time.time()
#             _loaded = False

#             def _hook(obj):  # DH0001
#                 """
#                 Hook to load ensemble
#                 """

#                 nonlocal _loaded
#                 if not _loaded:
#                     _loaded = True
#                     print(
#                         "\nLoaded {} data, now putting in a dictionary".format(
#                             filetype
#                         ),
#                         flush=True,
#                     )
#                 return obj

#             if filetype == "JSON":
#                 # orjson promises to be fast, but it doesn't seem to be
#                 # and fails on "Infinity"... oops
#                 # data = orjson.loads(file_object.read())

#                 # simplejson is faster than standard json and "just works"
#                 # on the big Moe set in 37s
#                 if not quiet:
#                     data = simplejson.load(file_object, object_hook=_hook)
#                 else:
#                     data = simplejson.load(file_object)
#                 file_object.close()

#                 # standard json module
#                 # on the big Moe set takes 42s
#                 # data = json.load(file_object,
#                 #                 object_hook=_hook)
#             elif filetype == "msgpack":
#                 data = msgpack.load(file_object, object_hook=_hook)
#                 file_object.close()

#             if timing:
#                 print(
#                     "\n\nTook {} s to load the data\n\n".format(time.time() - tstart),
#                     flush=True,
#                 )

#     # strip non-selected keys, if a list is given in select_keys
#     if select_keys:
#         keys = list(data["ensemble"].keys())
#         for key in keys:
#             if key not in select_keys:
#                 del data["ensemble"][key]

#     # perhaps convert floats?
#     tstart = time.time()
#     if convert_float_keys:
#         # timings are for 100 iterations on the big Moe data set
#         # data = format_ensemble_results(data) # 213s
#         # data = recursive_change_key_to_float(data) # 61s
#         data = keys_to_floats(data)  # 6.94s

#         if timing:
#             print(
#                 "\n\nTook {} s to convert floats\n\n".format(time.time() - tstart),
#                 flush=True,
#             )

#     # return data
#     return data


# class AutoVivificationDict(dict):  # DH0001
#     """
#     Implementation of perl's autovivification feature, by overriding the
#     get item and the __iadd__ operator (https://docs.python.org/3/reference/datamodel.html?highlight=iadd#object.__iadd__)

#     This allows to set values within a subdict that might not exist yet:

#     Example:
#         newdict = {}
#         newdict['example']['mass'] += 10
#         print(newdict)
#         >>> {'example': {'mass': 10}}
#     """

#     def __getitem__(self, item):  # DH0001
#         """
#         Getitem function for the autovivication dict
#         """

#         try:
#             return dict.__getitem__(self, item)
#         except KeyError:
#             value = self[item] = type(self)()
#             return value

#     def __iadd__(self, other):  # DH0001
#         """
#         iadd function (handling the +=) for the autovivication dict.
#         """

#         # if a value does not exist, assume it is 0.0
#         try:
#             self += other
#         except:
#             self = other
#         return self


# def merge_dicts(
#     dict_1: dict,
#     dict_2: dict,
#     use_ordereddict=True,
#     allow_matching_key_type_mismatch=True,
# ) -> dict:
#     """
#     Function to merge two dictionaries in a custom way. Taken from binarycpython

#     Behaviour:

#     When dict keys are only present in one of either:
#         - we just add the content to the new dict

#     When dict keys are present in both, we decide based on the value types how to combine them:
#         - dictionaries will be merged by calling recursively calling this function again
#         - numbers will be added
#         - (opt) lists will be appended
#         - booleans are merged with logical OR
#         - identical strings are just set to the string
#         - non-identical strings are concatenated
#         - NoneTypes are set to None
#         - In the case that the instances do not match: for now I will raise an error

#     Args:
#         dict_1: first dictionary
#         dict_2: second dictionary

#     Returns:
#         Merged dictionary

#     """

#     # Set up new dict
#     if use_ordereddict:
#         new_dict = collections.OrderedDict()
#     else:
#         new_dict = {}

#     ##################
#     #
#     keys_1 = dict_1.keys()
#     keys_2 = dict_2.keys()

#     ##################
#     # Find overlapping keys of both dicts
#     overlapping_keys = set(keys_1).intersection(set(keys_2))

#     # Find the keys that are unique
#     unique_to_dict_1 = set(keys_1).difference(set(keys_2))
#     unique_to_dict_2 = set(keys_2).difference(set(keys_1))

#     ##################
#     # Add the unique keys to the new dict
#     for key in unique_to_dict_1:
#         # If these items are numerical or string, then just put them in
#         if isinstance(dict_1[key], ALLOWED_NUMERICAL_TYPES + (str,)):
#             new_dict[key] = dict_1[key]
#         # Else, to be safe we should deepcopy them
#         else:
#             copy_dict = dict_1[key]
#             new_dict[key] = copy_dict

#     for key in unique_to_dict_2:
#         # If these items are numerical or string, then just put them in
#         if isinstance(dict_2[key], ALLOWED_NUMERICAL_TYPES + (str,)):
#             new_dict[key] = dict_2[key]
#         # Else, to be safe we should deepcopy them
#         else:
#             copy_dict = dict_2[key]
#             new_dict[key] = copy_dict

#     ##################
#     # Go over the common keys:
#     for key in overlapping_keys:

#         ##################
#         # If they keys are not the same, it depends on their type whether we still deal with them at all, or just raise an error
#         if not isinstance(dict_1[key], type(dict_2[key])):

#             ##################
#             # Exceptions: numbers can be added
#             if isinstance(dict_1[key], ALLOWED_NUMERICAL_TYPES) and isinstance(
#                 dict_2[key], ALLOWED_NUMERICAL_TYPES
#             ):
#                 new_dict[key] = dict_1[key] + dict_2[key]

#             ##################
#             # Exceptions: versions of dicts can be merged
#             elif isinstance(
#                 dict_1[key], (dict, collections.OrderedDict, type(AutoVivificationDict))
#             ) and isinstance(
#                 dict_2[key], (dict, collections.OrderedDict, type(AutoVivificationDict))
#             ):
#                 new_dict[key] = merge_dicts(
#                     dict_1[key],
#                     dict_2[key],
#                     use_ordereddict=use_ordereddict,
#                     allow_matching_key_type_mismatch=allow_matching_key_type_mismatch,
#                 )

#             ##################
#             #
#             if not allow_matching_key_type_mismatch:
#                 print(
#                     "Error key: {} value: {} type: {} and key: {} value: {} type: {} are not of the same type and cannot be merged".format(
#                         key,
#                         dict_1[key],
#                         type(dict_1[key]),
#                         key,
#                         dict_2[key],
#                         type(dict_2[key]),
#                     )
#                 )
#                 raise ValueError

#             ##################
#             # one key is None, just use the other
#             elif dict_1[key] is None:
#                 try:
#                     new_dict[key] = dict_2[key]
#                 except:
#                     msg = f"{key}: Failed to set from {dict_2[key]} when other key was of NoneType "
#                     raise ValueError(msg)

#             elif dict_1[key] is None:
#                 try:
#                     new_dict[key] = dict_1[key]
#                 except:
#                     msg = f"{key}: Failed to set from {dict_1[key]} when other key was of NoneType "
#                     raise ValueError(msg)

#             # string-int clash : convert both to ints and save
#             elif (
#                 isinstance(dict_1[key], str)
#                 and isinstance(dict_2[key], int)
#                 or isinstance(dict_1[key], int)
#                 and isinstance(dict_2[key], str)
#             ):
#                 try:
#                     new_dict[key] = int(dict_1[key]) + int(dict_2[key])
#                 except ValueError as e:
#                     msg = "{}: Failed to convert string (either '{}' or '{}') to an int".format(
#                         key, dict_1[key], dict_2[key]
#                     )
#                     raise ValueError(msg) from e

#             # string-float clash : convert both to floats and save
#             elif (
#                 isinstance(dict_1[key], str)
#                 and isinstance(dict_2[key], float)
#                 or isinstance(dict_1[key], float)
#                 and isinstance(dict_2[key], str)
#             ):
#                 try:
#                     new_dict[key] = float(dict_1[key]) + float(dict_2[key])
#                 except ValueError as e:
#                     msg = "{}: Failed to convert string (either '{}' or '{}') to an float".format(
#                         key, dict_1[key], dict_2[key]
#                     )
#                     raise ValueError(msg) from e

#             # If the above cases have not dealt with it, then we should raise an error
#             else:
#                 msg = "merge_dicts error: key: {key} value: {value1} type: {type1} and key: {key} value: {value2} type: {type2} are not of the same type and cannot be merged".format(
#                     key=key,
#                     value1=dict_1[key],
#                     type1=type(dict_1[key]),
#                     value2=dict_2[key],
#                     type2=type(dict_2[key]),
#                 )
#                 raise ValueError(msg)

#         # Here the keys are the same type
#         # Here we check for the cases that we want to explicitly catch. Ints will be added,
#         # floats will be added, lists will be appended (though that might change) and dicts will be
#         # dealt with by calling this function again.
#         else:
#             # ints
#             # Booleans (has to be the type Bool, not just a 0 or 1)
#             if isinstance(dict_1[key], bool) and isinstance(dict_2[key], bool):
#                 new_dict[key] = dict_1[key] or dict_2[key]

#             elif isinstance(dict_1[key], int) and isinstance(dict_2[key], int):
#                 new_dict[key] = dict_1[key] + dict_2[key]

#             elif isinstance(dict_1[key], np.int64) and isinstance(
#                 dict_2[key], np.int64
#             ):
#                 new_dict[key] = dict_1[key] + dict_2[key]

#             # floats
#             elif isinstance(dict_1[key], float) and isinstance(dict_2[key], float):
#                 new_dict[key] = dict_1[key] + dict_2[key]

#             # lists
#             elif isinstance(dict_1[key], list) and isinstance(dict_2[key], list):
#                 new_dict[key] = dict_1[key] + dict_2[key]

#             # Astropy quantities (using a dummy type representing the numpy array)
#             elif isinstance(dict_1[key], type(np.array([1]) * u.m)) and isinstance(
#                 dict_2[key], type(np.array([1]) * u.m)
#             ):
#                 new_dict[key] = dict_1[key] + dict_2[key]

#             # dicts
#             elif isinstance(dict_1[key], dict) and isinstance(dict_2[key], dict):
#                 new_dict[key] = merge_dicts(
#                     dict_1[key],
#                     dict_2[key],
#                     use_ordereddict=use_ordereddict,
#                     allow_matching_key_type_mismatch=allow_matching_key_type_mismatch,
#                 )

#             # strings
#             elif isinstance(dict_1[key], str) and isinstance(dict_2[key], str):
#                 if dict_1[key] == dict_2[key]:
#                     # same strings
#                     new_dict[key] = dict_1[key]
#                 else:
#                     # different strings: just concatenate them
#                     new_dict[key] = dict_1[key] + dict_2[key]

#             # None types
#             elif dict_1[key] is None and dict_2[key] is None:
#                 new_dict[key] = None

#             else:
#                 msg = "Object types {}: {} ({}), {} ({}) not supported.".format(
#                     key,
#                     dict_1[key],
#                     type(dict_1[key]),
#                     dict_2[key],
#                     type(dict_2[key]),
#                 )
#                 raise ValueError(msg)

#     #
#     return new_dict


# def multiply_ensemble(ensemble, factor):
#     """
#     Function to multiply the endpoints for
#     """

#     for key in ensemble.keys():
#         if isinstance(ensemble[key], dict):
#             multiply_ensemble(ensemble=ensemble[key], factor=factor)
#         elif isinstance(ensemble[key], (int, float)):
#             ensemble[key] = ensemble[key] * factor


# def check_if_value_layer(keys):
#     """
#     Function to determine whether the supplied set of keys consists of value-keys or name-keys
#     """

#     is_value_layer = True
#     for key in keys:
#         try:
#             float(key)
#         except ValueError:
#             is_value_layer = False
#     return is_value_layer


# def get_layer_iterable(ensemble, is_value_layer):
#     """
#     Function to get the layer iterable
#     """

#     ########
#     # determine iterable
#     if is_value_layer:
#         iterable = sorted(ensemble.keys(), key=lambda x: float(x))
#     else:
#         iterable = ensemble.keys()
#     return iterable


# def check_if_value_layer_and_get_layer_iterable(ensemble):
#     """
#     Function to handle checking whether the layer is a data layer and then build up the iterable
#     """

#     #
#     is_value_layer = check_if_value_layer(ensemble.keys())
#     iterable = get_layer_iterable(ensemble, is_value_layer)

#     return is_value_layer, iterable


# def get_ensemble_binsizes(config, ensemble, data_layer_dict_entry):
#     """
#     Function to calculate the binsizes, taking into account transformations of the numbers.
#     """

#     #######
#     # TODO: make this optional
#     calculate_edges_before_transformations = True

#     #######
#     # if binsizes are provided then just use those
#     if "binsizes" in data_layer_dict_entry:
#         return data_layer_dict_entry["binsizes"]

#     #######
#     # loop over sorted (float) version of the ensemble keys:
#     sorted_ensemble_keys = list(sorted(ensemble.keys(), key=lambda x: float(x)))

#     #
#     values = np.array([float(el) for el in sorted_ensemble_keys])

#     # Transform edges
#     if calculate_edges_before_transformations:

#         # calculate edge values
#         edge_values = calculate_bin_edges(values)

#         # perform transformations
#         transformed_edge_values = np.array(
#             [
#                 handle_custom_scaling_or_conversion(
#                     config=config,
#                     data_layer_or_column_dict_entry=data_layer_dict_entry,
#                     value=edge_value,
#                 )
#                 for edge_value in edge_values
#             ]
#         )
#     else:
#         # perform transformations
#         transformed_values = np.array(
#             [
#                 handle_custom_scaling_or_conversion(
#                     config=config,
#                     data_layer_or_column_dict_entry=data_layer_dict_entry,
#                     value=value,
#                 )
#                 for value in values
#             ]
#         )

#         # calculate edge values
#         transformed_edge_values = calculate_bin_edges(transformed_values)

#     # calculate binsizes
#     binsizes = np.diff(transformed_edge_values)

#     #
#     config["logger"].debug(
#         "Values: {} transformed_edge_values: {} binsizes: {}".format(
#             values, transformed_edge_values, binsizes
#         )
#     )

#     return binsizes


# ################
# # data-layer dict functionality
# def get_data_layer_dict_values(data_layer_dict):
#     """
#     function to extract data layer values
#     """

#     data_layer_values = []

#     for key in data_layer_dict:
#         if isinstance(data_layer_dict[key], int):
#             data_layer_values.append(data_layer_dict[key])
#         elif isinstance(data_layer_dict[key], dict):
#             data_layer_values.append(data_layer_dict[key]["layer_depth"])
#         else:
#             raise ValueError("input type not supported.")

#     return data_layer_values


# def get_deepest_data_layer_depth(data_layer_dict):
#     """
#     Function to get the deepest data layer depth
#     """

#     data_layer_values = get_data_layer_dict_values(data_layer_dict=data_layer_dict)

#     return max(data_layer_values)


# def shift_layers_dict(data_layer_dict, shift_value):
#     """
#     Function to shift the layer depths
#     """

#     new_layer_depth_dict = copy.copy(data_layer_dict)

#     for key in data_layer_dict.keys():
#         shift_data_layer(
#             data_layer_dict=new_layer_depth_dict, key=key, shift=shift_value
#         )

#     return new_layer_depth_dict


# def shift_data_layer(data_layer_dict, key, shift):
#     """
#     Function to shift the data layer with a certain value
#     """

#     if isinstance(data_layer_dict[key], int):
#         data_layer_dict[key] += shift
#     elif isinstance(data_layer_dict[key], dict):
#         data_layer_dict[key]["layer_depth"] += shift
#     else:
#         raise ValueError("input type not supported.")


# def invert_data_layer_dict(data_layer_dict):
#     """
#     Function to invert the data layer dictionary.

#     This function does not truly swap key and values because if entries in the
#     data_layer_dict contain dictionaries we will extract the layer depth from
#     them instead of making the dictionary the key.
#     """

#     inverted_data_layer_dict = {}

#     for key in data_layer_dict:
#         if isinstance(data_layer_dict[key], int):
#             inverted_data_layer_dict[data_layer_dict[key]] = key
#         elif isinstance(data_layer_dict[key], dict):
#             inverted_data_layer_dict[data_layer_dict[key]["layer_depth"]] = key
#         else:
#             raise ValueError("input type not supported.")

#     return inverted_data_layer_dict


# ##################
# # other
# def shift_layers_list(layer_list, shift_value):
#     """
#     Function to shift the entries in a list by a certain value
#     """

#     for el_i, _ in enumerate(layer_list):
#         layer_list[el_i] += shift_value

#     return layer_list


# ################
# # endpoints functionality
# def extract_endpoints(ensemble, endpoint_list=None, found_units=False):
#     """
#     Function to strip the endpoints
#     """

#     if endpoint_list is None:
#         endpoint_list = []

#     # Handle recursive
#     if isinstance(ensemble, (dict, OrderedDict)):
#         for key in ensemble.keys():
#             if isinstance(ensemble[key], (dict, OrderedDict)):
#                 endpoint_list, unit_endpoint_list = extract_endpoints(
#                     ensemble[key],
#                     endpoint_list=list(endpoint_list),
#                     found_units=found_units,
#                 )
#             elif has_unit_dimensionless_okay(ensemble[key]):
#                 found_units = True
#                 endpoint_list.append(ensemble[key])
#             elif isinstance(ensemble[key], (int, float)):
#                 endpoint_list.append(ensemble[key])

#         return endpoint_list, found_units

#     return np.array(endpoint_list), found_units


# def attach_endpoints(ensemble, endpoint_array, counter=0, depth=0):
#     """
#     Function to attach endpoints
#     """

#     #
#     if isinstance(ensemble, (dict, OrderedDict)):

#         #
#         for key in ensemble.keys():
#             if isinstance(ensemble[key], (dict, OrderedDict)):
#                 counter = attach_endpoints(
#                     ensemble[key],
#                     endpoint_array=endpoint_array,
#                     counter=counter,
#                     depth=depth + 1,
#                 )
#             elif isinstance(ensemble[key], (int, float)):
#                 ensemble[key] = endpoint_array[counter]
#                 counter += 1
#             elif has_unit_dimensionless_okay(ensemble[key]):
#                 ensemble[key] = endpoint_array[counter]
#                 counter += 1
#             if ensemble[key] is None:
#                 ensemble[key] = endpoint_array[counter]
#                 counter += 1

#         # if we're back all the way up and we've not exhausted the array then something is wrong
#         if depth == 0:
#             if len(endpoint_array) != counter:
#                 raise ValueError(
#                     "Somehow we have not reached the end of the endpoint array"
#                 )
#         return counter

#     raise ValueError("Unsupported type: {}".format(type(ensemble)))


# def set_endpoints(ensemble, value):
#     """
#     Function to set the endpoints to a fixed value
#     """

#     if isinstance(ensemble, (dict, OrderedDict)):
#         for key in ensemble.keys():
#             if isinstance(ensemble[key], (dict, OrderedDict)):
#                 set_endpoints(ensemble=ensemble[key], value=value)
#             elif isinstance(ensemble[key], (int, float)):
#                 ensemble[key] = value
#             elif has_unit_dimensionless_okay(ensemble[key]):
#                 ensemble[key] = value
#             elif ensemble[key] is None:
#                 ensemble[key] = value

#         if len(ensemble.keys()) == 0:
#             raise ValueError("Encountered empty ensemble")

#     else:
#         raise ValueError("Unsupported type: {}".format(type(ensemble)))


# def strip_ensemble_endpoints(ensemble):
#     """
#     Function to strip the ensemble from its endpoints
#     """

#     # extract endpoints
#     endpoints, found_units = extract_endpoints(ensemble=ensemble)

#     # set original ensemble endpoints to 0
#     set_endpoints(ensemble=ensemble, value=0)

#     return ensemble, endpoints, found_units


# ################
# def get_depth_ensemble_first_endpoint(ensemble, depth=0):
#     """
#     Function to get the maximum depth in an ensemble. Stops at the first branch
#     that does not contain a key or nested dict (i.e. this function does not crawl through the entire ensemble)
#     """

#     if isinstance(ensemble, (dict)):
#         for key in ensemble.keys():
#             return get_depth_ensemble_first_endpoint(ensemble[key], depth=depth + 1)
#     else:
#         return depth


# def get_max_depth_ensemble(ensemble):
#     """
#     Function to find the maximum depth in a nested dictionary
#     """

#     def _find_max_depth_helper(ensemble, depth):  # DH0001
#         if not isinstance(ensemble, dict):
#             return depth
#         max_depth = depth
#         for value in ensemble.values():
#             max_depth = max(max_depth, _find_max_depth_helper(value, depth + 1))
#         return max_depth

#     return _find_max_depth_helper(ensemble, 0)


# def get_depth_ensemble_all_endpoints(ensemble, depth=0, endpoint_depths=None):
#     """
#     Function to get the maximum depth in an ensemble. Stops at the first branch
#     that does not contain a key or nested dict (i.e. this function does not crawl through the entire ensemble)
#     """

#     if endpoint_depths is None:
#         endpoint_depths = []

#     if isinstance(ensemble, (dict)):
#         for key in ensemble.keys():
#             get_depth_ensemble_all_endpoints(
#                 ensemble[key], depth=depth + 1, endpoint_depths=endpoint_depths
#             )
#     else:
#         endpoint_depths.append(depth)
#     return endpoint_depths


# ################
# # Functions to get the ensemble structure.
# def _get_ensemble_structure(ensemble, structure_dict, max_depth, depth=0):
#     """
#     Recursive function acompanying the "get_ensemble_structure" function
#     """

#     #
#     if depth < max_depth:
#         for key in ensemble.keys():
#             # add to structure
#             if key not in structure_dict[depth]:
#                 structure_dict[depth].append(key)

#             # check the rest
#             structure_dict = _get_ensemble_structure(
#                 ensemble=ensemble[key],
#                 structure_dict=structure_dict,
#                 max_depth=max_depth,
#                 depth=depth + 1,
#             )

#     return structure_dict


# def get_ensemble_structure(ensemble, named_layer_list=None):
#     """
#     Function to generate the ensemble structure.

#     Optionally, if a named layer list is provided, this function will raise an error if named layers contain more than one value
#     """

#     # check if there are endpoints of varying depth
#     depth_all_endpoints = get_depth_ensemble_all_endpoints(ensemble)

#     if len(np.unique(depth_all_endpoints)) > 1:
#         raise ValueError("This ensemble has endpoints of varying depth. abort")

#     #
#     max_depth = depth_all_endpoints[0]

#     # setup structure dict
#     structure_dict = {el: [] for el in range(max_depth)}

#     # generate structure
#     structure_dict = _get_ensemble_structure(
#         ensemble=ensemble, structure_dict=structure_dict, max_depth=max_depth
#     )

#     # check content if necessary
#     if named_layer_list is not None:
#         for named_layer_depth in named_layer_list:
#             if len(structure_dict[named_layer_depth]) > 1:
#                 raise ValueError("Multiple values at named layer depth")

#     return structure_dict


# ##############
# # marginaliser functions


# def ensemble_marginalise_layer(ensemble, marginalisation_depth, depth=0):
#     """
#     Function to marginalise a layer of the ensemble
#     """

#     if marginalisation_depth == 0:
#         raise ValueError("We can't currently marginalise with marginalisation_depth=-1")

#     if depth + 1 == marginalisation_depth:
#         # merge subdicts
#         for key in ensemble.keys():
#             merged_subdicts = {}
#             for subdict_key in ensemble[key].keys():
#                 merged_subdicts = merge_dicts(
#                     merged_subdicts,
#                     ensemble[key][subdict_key],
#                     use_ordereddict=False,
#                     allow_matching_key_type_mismatch=False,
#                 )
#             ensemble[key] = merged_subdicts
#     else:
#         for key in ensemble.keys():
#             ensemble[key] = ensemble_marginalise_layer(
#                 ensemble=ensemble[key],
#                 depth=depth + 1,
#                 marginalisation_depth=marginalisation_depth,
#             )

#     return ensemble


###################
# Inflation functions:


# def inflate_ensemble(ensemble_data):
#     """
#     Function to inflate an ensemble, taking all the values for each datalayer and making a rectangular grid for it

#     The first value should be a namelayer
#     """

#     parameter_name = list(ensemble_data.keys())[0]

#     next_layer_keys = list(ensemble_data[parameter_name].keys())

#     first_of_next_next_layer = ensemble_data[parameter_name][next_layer_keys[0]]
#     is_final_layer = not isinstance(first_of_next_next_layer, (dict, OrderedDict))

#     # if this is the final layer, then handle this layer with the dedicated function
#     if is_final_layer:
#         flattened_data = flatten_data_ensemble1d(ensemble_data, parameter_name).T
#         return flattened_data

#     # If its not the final layer, we should call this function again and return the result
#     combined_array = None
#     for valuekey in next_layer_keys:
#         # Check if its the final layer
#         next_next_layer = ensemble_data[parameter_name][valuekey]

#         # combine result with the current value key
#         res = inflate_ensemble(next_next_layer)

#         # look at how many rows we have
#         rows = res.shape[0]
#         cols = res.shape[-1]

#         # Get a column with the current valuekeys
#         new_column = np.array([float(valuekey)] * rows)

#         # Create a new array with the shape that can fit the new column
#         new_array = np.zeros((rows, cols + 1))

#         # Set the new column as the first one and the rest in the rest
#         new_array[:, 0] = new_column
#         new_array[:, 1:] = res

#         # Combine the arrays
#         if combined_array is None:
#             combined_array = new_array
#         else:
#             combined_array = np.append(combined_array, new_array, axis=0)

#     # Return the results
#     return combined_array


def flatten_data_ensemble1d(input_dict, named_subkey_1=None):
    """
    Functon to get the subkey and its associated value from a dictionary.

    This dictionary can't be deeper than 2 levels (i.e. it has to be the last named key.

    This requires the input dictionary to be structured as follows:

    {
        named_subkey_1:
        {
            numerical subkeys 1:
                value,
            ..
        }
    }
    """

    #
    if named_subkey_1 is None:
        _input_dict = input_dict
    else:
        _input_dict = input_dict[named_subkey_1]

    #
    data = []
    for key_1 in sorted(_input_dict):
        value = _input_dict[key_1]
        data.append([key_1, float(value)])
    data = np.array(data).T

    return data


def inflate_ensemble_with_lists_and_named_layers(ensemble_data):
    """
    Function to inflate an ensemble by using dataframes. taking all the values for each datalayer and making a rectangular grid for it

    this function assumes named layers are present between each value layer
    """

    parameter_name = list(ensemble_data.keys())[0]
    if isinstance(ensemble_data[parameter_name], ALLOWED_NUMERICAL_TYPES):
        return [[ensemble_data[parameter_name]]]

    next_layer_keys = list(ensemble_data[parameter_name].keys())
    first_of_next_next_layer = ensemble_data[parameter_name][next_layer_keys[0]]
    is_final_layer = not isinstance(first_of_next_next_layer, (dict, OrderedDict))

    # if this is the final layer, then handle this layer with the dedicated function
    if is_final_layer:
        flattened_data = flatten_data_ensemble1d(
            ensemble_data, named_subkey_1=parameter_name
        )

        flattened_data_list = [[], []]
        flattened_data_list[0] = flattened_data[0].tolist()
        flattened_data_list[1] = flattened_data[1].astype(float).tolist()

        return flattened_data_list

    # If its not the final layer, we should call this function again and return the result
    combined_list = None
    for valuekey in next_layer_keys:
        next_next_layer = ensemble_data[parameter_name][valuekey]

        # combine result with the current value key
        res = inflate_ensemble_with_lists_and_named_layers(next_next_layer)

        if combined_list is None:
            combined_list = [[] for _ in range(len(res) + 1)]
        valuekey_list = [valuekey for _ in range(len(res[0]))]

        # Create add the current key to the combined_list
        combined_list[0] = combined_list[0] + valuekey_list

        # add the results of the previous inflation to the results
        for res_list_i, res_list in enumerate(res):
            combined_list[res_list_i + 1] = combined_list[res_list_i + 1] + res_list

    # Return the results
    return combined_list


def inflate_ensemble_with_lists_without_named_layers(ensemble_data):
    """
    Function to inflate an ensemble by using dataframes. taking all the values for each datalayer and making a rectangular grid for it

    this function assumes named layers are present between each value layer
    """

    # parameter_name = list(ensemble_data.keys())[0]
    next_layer_keys = list(ensemble_data.keys())
    first_of_next_layer = ensemble_data[next_layer_keys[0]]
    is_final_layer = not isinstance(first_of_next_layer, (dict, OrderedDict))

    # if this is the final layer, then handle this layer with the dedicated function
    if is_final_layer:
        flattened_data = flatten_data_ensemble1d(ensemble_data)

        flattened_data_list = [[], []]
        flattened_data_list[0] = flattened_data[0].tolist()
        flattened_data_list[1] = flattened_data[1].astype(float).tolist()

        return flattened_data_list

    # If its not the final layer, we should call this function again and return the result
    combined_list = None
    for valuekey in next_layer_keys:
        next_layer = ensemble_data[valuekey]

        # combine result with the current value key
        res = inflate_ensemble_with_lists_without_named_layers(next_layer)

        if combined_list is None:
            combined_list = [[] for _ in range(len(res) + 1)]
        valuekey_list = [valuekey for _ in range(len(res[0]))]

        # Create add the current key to the combined_list
        combined_list[0] = combined_list[0] + valuekey_list

        # add the results of the previous inflation to the results
        for res_list_i, res_list in enumerate(res):
            combined_list[res_list_i + 1] = combined_list[res_list_i + 1] + res_list

    # Return the results
    return combined_list


def find_columnames_recursively(ensemble_data, columnnames=None):
    """
    Function to find all the column names recursively

    This function should only be called on ensemble datasets that do not have differnt tree structures in them

    The first layer should be a namelayer
    """

    # get the column name
    if columnnames is None:
        new_columnnames = [list(ensemble_data.keys())[0]]
    else:
        new_columnnames = columnnames + [list(ensemble_data.keys())[0]]

    # Check if we are in the lowest layer
    next_layer_keys = list(ensemble_data[new_columnnames[-1]].keys())
    next_next_layer = ensemble_data[new_columnnames[-1]][next_layer_keys[0]]

    # Call itself or return if
    if isinstance(next_next_layer, (dict, OrderedDict)):
        return find_columnames_recursively(next_next_layer, columnnames=new_columnnames)

    return new_columnnames


def convert_ensemble_to_dataframe(
    ensemble_data, verbose=False, contains_named_layers=True, columnames=None
):
    """
    Function to inflate an ensemble, which will transform an ensemble (i.e. a nested histogram),
    into a rectangular representation of the same data in the form of a pandas dataframe.

    This will increase the size of the data by about a factor of 2 (i think),
    but will make certain operations much easier as we can extract a particular column as a numpy array.


    If the dataframe has named-layers, this can be handled by setting `contains_named_layers=True`. If it does not contain named layers, you can provide the column names manually through the columnnames parameter

    The final layer is assumed called 'probability' (used in grid-based pop-synth), but it may be that that layer contains e.g. normalized yield instead.
    """

    # Convert to dataframe
    start = time.time()

    if contains_named_layers:
        try:
            columnames = find_columnames_recursively(ensemble_data)
            data_list = inflate_ensemble_with_lists_and_named_layers(ensemble_data)

        except:
            columnames = find_columnames_recursively({"ensemble": ensemble_data})
            data_list = inflate_ensemble_with_lists_and_named_layers(
                {"ensemble": ensemble_data}
            )

        #
    else:
        if columnames is None:
            raise ValueError(
                "If the data does not contain named layers you should provide the column names yourself"
            )
        #
        data_list = inflate_ensemble_with_lists_without_named_layers(ensemble_data)

    ##########
    #
    df = pd.DataFrame(data_list)
    df = df.transpose()

    try:
        df.columns = columnames + ["probability"]
    except:
        df.columns = columnames[1:] + ["probability"]

    if verbose:
        stop = time.time()
        print("Converting ensemble data to dataframe took {}s".format(stop - start))

    return df


if __name__ == "__main__":
    # example_ensemble_data = {
    #     "a": {
    #         "5": {
    #             "b": {
    #                 "1": 0.5,
    #                 "2": 0.5,
    #             }
    #         },
    #         "6": {
    #             "b": {
    #                 "3": 0.5,
    #                 "4": 0.5,
    #             }
    #         },
    #         "7": {
    #             "b": {
    #                 "8": 0.5,
    #                 "9": 0.5,
    #             }
    #         },
    #     }
    # }

    # df = convert_ensemble_to_dataframe(
    #     example_ensemble_data, contains_named_layers=True
    # )
    # print(df)

    import json

    import pkg_resources

    # from syntheticstellarpopconvolve.ensemble_utils import inflate_ensemble_with_lists_and_named_layers, convert_ensemble_to_dataframe
    # load the data
    example_ensemble_filename = pkg_resources.resource_filename(
        "syntheticstellarpopconvolve", "example_data/example_ensemble.json"
    )
    with open(example_ensemble_filename, "r") as f_ensemble:
        ensemble = json.loads(f_ensemble.read())
    # print(ensemble['ensemble']['Xyield'])

    # Xyield = ensemble['ensemble']

    inflated_ensemble = convert_ensemble_to_dataframe(
        ensemble_data=ensemble["ensemble"]["Xyield"],
        verbose=False,
        contains_named_layers=True,
    )
    print(inflated_ensemble)

    print(inflated_ensemble["probability"].to_numpy())
