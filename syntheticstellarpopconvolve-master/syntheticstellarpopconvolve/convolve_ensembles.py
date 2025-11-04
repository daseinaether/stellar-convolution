# """
# Ensemble convolution functions

# TODO: put all the general ensemble files into the ensemble file
# TODO: allow the yield-rate and extra weights function to pass units to the ensemble and when
# stripping the end-points we should make sure to handle that properly.
# TODO: move the ensemble utils functions to a different place maybe.
# """

# import json

# import h5py
# import numpy as np

# from syntheticstellarpopconvolve.default_convolution_config import (
#     ALLOWED_NUMERICAL_TYPES,
# )
# from syntheticstellarpopconvolve.ensemble_utils import (
#     attach_endpoints,
#     check_if_value_layer_and_get_layer_iterable,
#     ensemble_marginalise_layer,
#     get_data_layer_dict_values,
#     get_deepest_data_layer_depth,
#     get_ensemble_binsizes,
#     invert_data_layer_dict,
#     multiply_ensemble,
#     shift_data_layer,
#     shift_layers_dict,
#     shift_layers_list,
#     strip_ensemble_endpoints,
# )
# from syntheticstellarpopconvolve.general_functions import (
#     calculate_digitized_sfr_rates,
#     handle_custom_scaling_or_conversion,
# )
# from syntheticstellarpopconvolve.post_convolution_hook_routines import (
#     handle_post_convolution_function,
# )


# def convolve_ensemble_integration_post_convolution_hook_wrapper(
#     config,
#     sfr_dict,
#     data_dict,
#     time_bin_info_dict,
#     convolution_instruction,
#     ensemble,
#     #
#     persistent_data=None,
#     previous_convolution_results=None,
# ):
#     """
#     Function to wrap the post-convolution function call for ensemble-convolution by integration.

#     rules:
#     - additional data can be added to the result_dict
#     - the number of systems can lower than before the call

#     NOTE: the result_dict is expected to be updated in-place.
#     TODO: perhaps we should inflate the ensemble and pass it as the data for the data dict here.
#     """

#     #
#     name = "convolve-ensemble by integration"

#     #
#     config["logger"].warning(
#         "Handling post-convolution function hook call for {}".format(name)
#     )

#     #############
#     # pre-call setup
#     stripped_ensemble, stripped_endpoints, found_units = strip_ensemble_endpoints(
#         ensemble=ensemble
#     )

#     convolution_results = {"yield": stripped_endpoints}

#     #
#     num_systems_before = len(convolution_results[list(convolution_results.keys())[0]])

#     #############
#     # call hook
#     handle_post_convolution_function(
#         config=config,
#         sfr_dict=sfr_dict,
#         data_dict=data_dict,
#         time_bin_info_dict=time_bin_info_dict,
#         convolution_instruction=convolution_instruction,
#         convolution_results=convolution_results,
#         name=name,
#         #
#         persistent_data=persistent_data,
#         previous_convolution_results=previous_convolution_results,
#     )

#     #############
#     # check output

#     # Abort when its a list
#     if isinstance(convolution_results, list):
#         raise ValueError(
#             "Currently returning multiple convolution_results after calling the post-convolution function is not supported for convolution of ensemble-based results."
#         )

#     # if its a dict, check the length
#     num_systems_after = len(convolution_results[list(convolution_results.keys())[0]])

#     #
#     if num_systems_before != num_systems_after:
#         raise ValueError(
#             "{} post-convolution function has changed the number of systems stored in the output dict. Due to current data structure decisions this is not supported currently. Please make sure that the number of systems before and after calling this function stays equal.".format(
#                 name
#             )
#         )

#     #
#     num_output_entries = len(convolution_results.keys())

#     if num_output_entries > 1:
#         raise ValueError(
#             "{} post-convolution function has added additional entries to the result dictionary. Due tot current data structure decisions this is not supported currently. Please make sure that the number the result dictionary only contains one entry".format(
#                 name
#             )
#         )

#     ################
#     # re-attach the updated data to the ensemble
#     # TODO: allow the result dict to contain more data than just 1 number.
#     attach_endpoints(ensemble=ensemble, endpoint_array=convolution_results["yield"])

#     return ensemble


# def handle_binsize_multiplication_factor(
#     config,
#     ensemble,
#     data_layer_dict_entry,
#     key,
#     key_i,
#     binsizes,
#     name,
#     extra_value_dict,
# ):
#     """
#     Function to handle the calculation of the binsize multiplication
#     """

#     ##########
#     # check if we want to multiply this by the binsizes of the ensemble
#     if data_layer_dict_entry.get("multiply_by_binsize", False):
#         # Determine binsizes

#         # Calculate the binsizes (otherwise assume its already calculated)
#         if binsizes is None:
#             binsizes = get_ensemble_binsizes(
#                 config=config,
#                 ensemble=ensemble,
#                 data_layer_dict_entry=data_layer_dict_entry,
#             )

#             #
#             config["logger"].debug("Calculated binsizes: {}".format(binsizes))

#         # select current binsize
#         binsize = binsizes[key_i]

#         #
#         config["logger"].debug(
#             "Selected binsize {} for {} ({})".format(binsize, key, key_i)
#         )

#         # # TODO: check if the binsize extends beyond the last starformation
#         # binsize = restrict_binsize()
#         # config["logger"].debug("Restricted binsize to {}".format(binsize))

#         #
#         extra_value_dict["{}_binsize".format(name)] = binsize

#     return binsizes, extra_value_dict


# def ensemble_handle_marginalisation(
#     config, ensemble, convolution_instruction, is_pre_conv
# ):
#     """
#     Function to handle ensemble marginalisation pre and post convolution.

#     Inspects the marginalisation layer dict and marginalises deepest-first

#     If a data dict is provided and the 'is_pre_conv' flag is set to true, layers that match those in the data dict get skipped
#     """

#     ########
#     # if there is no marginalisation provided, we don't really have to do anything
#     if "marginalisation_list" not in convolution_instruction.keys():
#         return config, ensemble, convolution_instruction

#     #
#     config["logger"].debug(
#         "Handling ensemble marginalisation (is_pre_conv: {})".format(is_pre_conv)
#     )

#     ########
#     # find out which layers we can actually remove (to avoid removing the actual data we need to convolve.)
#     to_remove_layers = []
#     for layer in convolution_instruction["marginalisation_list"]:
#         # skip if the current layer exists in the data layer dict and we are doing pre-convolution
#         if is_pre_conv:
#             if layer in get_data_layer_dict_values(
#                 data_layer_dict=convolution_instruction["data_layer_dict"]
#             ):
#                 continue
#         to_remove_layers.append(layer)

#     #########
#     # Store which layers we should keep
#     to_keep_layers = list(
#         set(convolution_instruction["marginalisation_list"]) - set(to_remove_layers)
#     )
#     convolution_instruction["marginalisation_list"] = to_keep_layers

#     # sort list so we do deepest first
#     to_remove_layers = sorted(to_remove_layers, reverse=True)

#     #########
#     # Process layers
#     for layer in to_remove_layers:
#         #
#         config["logger"].debug("Marginalising layer {}".format(layer))

#         # perform marginalisation
#         ensemble = ensemble_marginalise_layer(
#             ensemble=ensemble, marginalisation_depth=layer
#         )

#         # update any data layer dict entry with an updated depth
#         if is_pre_conv:
#             for key, value in convolution_instruction["data_layer_dict"].items():
#                 if value > layer:
#                     shift_data_layer(
#                         data_layer_dict=convolution_instruction["data_layer_dict"],
#                         key=key,
#                         shift=-1,
#                     )

#         # Also update the part of the marginalisation list that will not currently be handled
#         for index, value in enumerate(convolution_instruction["marginalisation_list"]):
#             if value > layer:
#                 convolution_instruction["marginalisation_list"][index] = value - 1

#     return config, ensemble, convolution_instruction


# ###########
# # Main ensemble convolution functions
# def extract_ensemble_data(config, convolution_instruction):
#     """
#     Function to extract the ensemble-type data
#     """

#     #
#     config["logger"].debug(
#         "Extracting ensemble data {}".format(convolution_instruction["input_data_name"])
#     )

#     data_dict = {}
#     with h5py.File(config["output_filename"], "r") as output_hdf5file:
#         data_dict["ensemble_data"] = {
#             convolution_instruction["input_data_name"]: json.loads(
#                 output_hdf5file[
#                     "input_data/ensemble/{}".format(
#                         convolution_instruction["input_data_name"]
#                     )
#                 ][()]
#             # )
#         }

#     ##########
#     # we need to shift the layers because the way the data gets stored is by
#     # adding a layer at the start
#     convolution_instruction["data_layer_dict"] = shift_layers_dict(
#         convolution_instruction["data_layer_dict"], shift_value=1
#     )
#     if "marginalisation_list" in convolution_instruction:
#         convolution_instruction["marginalisation_list"] = shift_layers_list(
#             convolution_instruction["marginalisation_list"], shift_value=1
#         )

#     ##########
#     # handle pre-convolution ensemble marginalisation
#     (
#         config,
#         data_dict["ensemble_data"],
#         convolution_instruction,
#     ) = ensemble_handle_marginalisation(
#         config=config,
#         ensemble=data_dict["ensemble_data"],
#         convolution_instruction=convolution_instruction,
#         is_pre_conv=True,
#     )

#     return config, data_dict, convolution_instruction


# def ensemble_handle_SFR_multiplication(
#     config,
#     convolution_instruction,
#     sfr_dict,
#     time_bin_info_dict,
#     ensemble,
#     data_dict,
#     extra_value_dict=None,
#     #
#     persistent_data=None,
#     previous_convolution_results=None,
# ):
#     """
#     Function to handle multiplying the provided ensemble with a.
#     """

#     #
#     if extra_value_dict is None:
#         extra_value_dict = {}
#     extra_value = np.prod(list(extra_value_dict.values()))

#     #
#     config["logger"].debug(
#         "Convolving ensemble data with SFR rate at convolution bin {} ({}) with data_dict: {} and multiplying by {} ({})".format(
#             time_bin_info_dict["bin_number"],
#             time_bin_info_dict["bin_center"],
#             data_dict,
#             extra_value,
#             extra_value_dict,
#         )
#     )

#     # TODO: put below in a dedicated function

#     ##############
#     # to re-use the event-based functionality we should cast all the data in the data_dict into numpy arrays
#     new_data_dict = {}
#     for key, value in data_dict.items():
#         try:
#             new_data_dict[key] = np.array([value])
#         except TypeError:
#             new_data_dict[key] = np.array([value.value]) * value.unit
#     data_dict = new_data_dict

#     #############
#     digitized_sfr_rates = calculate_digitized_sfr_rates(
#         config=config,
#         convolution_time_bin_center=time_bin_info_dict["bin_center"],
#         data_dict=data_dict,
#         sfr_dict=sfr_dict,
#     )

#     # Multiply ensemble with SFR and extra value
#     multiply_ensemble(ensemble=ensemble, factor=digitized_sfr_rates[0] * extra_value)

#     ################
#     # Handle post-convolution
#     ensemble = convolve_ensemble_integration_post_convolution_hook_wrapper(
#         config=config,
#         sfr_dict=sfr_dict,
#         time_bin_info_dict=time_bin_info_dict,
#         data_dict=data_dict,
#         convolution_instruction=convolution_instruction,
#         ensemble=ensemble,
#         #
#         persistent_data=persistent_data,
#         previous_convolution_results=previous_convolution_results,
#     )

#     return ensemble


# def ensemble_convolve_ensemble(
#     config,
#     sfr_dict,
#     convolution_instruction,
#     time_bin_info_dict,
#     ensemble,
#     depth=0,
#     data_dict=None,
#     extra_value_dict=None,
#     #
#     persistent_data=None,
#     previous_convolution_results=None,
# ):
#     """
#     Recursive function that handles convolving the ensemble.
#     """

#     ########
#     #
#     if data_dict is None:
#         data_dict = {}
#     if extra_value_dict is None:
#         extra_value_dict = {}
#     binsizes = None

#     ########
#     # unpack some data
#     data_layer_dict = convolution_instruction["data_layer_dict"]
#     inverted_data_layer_dict = convolution_instruction["inverted_data_layer_dict"]
#     deepest_data_layer_depth = convolution_instruction["deepest_data_layer_depth"]
#     data_layer_values = convolution_instruction["data_layer_values"]

#     #
#     config["logger"].debug(
#         "Convolving ensemble at depth {} using data_dict: {} data_layer_dict: {}, deepest_data_layer_depth: {}, data_layer_values: {}".format(
#             depth,
#             data_dict,
#             data_layer_dict,
#             deepest_data_layer_depth,
#             data_layer_values,
#         )
#     )

#     ########
#     #
#     if isinstance(ensemble, dict):

#         ########
#         # check if we are in a value-layer
#         is_value_layer, layer_iterable = check_if_value_layer_and_get_layer_iterable(
#             ensemble
#         )
#         config["logger"].debug("Iterable: {}".format(layer_iterable))

#         ########
#         # Go over the keys
#         for key_i, key in enumerate(layer_iterable):
#             config["logger"].debug("Current layer key: {} ({})".format(key, key_i))

#             #################
#             # check if the layer is one of the target depths to pick up data for the data dict and to find
#             if depth in data_layer_values:
#                 #
#                 name = inverted_data_layer_dict[depth]

#                 config["logger"].debug(
#                     "The layer ({}) is that of data layer {}. Storing value {} in data_dict under {}".format(
#                         depth,
#                         name,
#                         key,
#                         name,
#                     )
#                 )

#                 ###########
#                 # if we are in a data layer but also this layer is not a value layer then there is an issue
#                 if not is_value_layer:
#                     raise ValueError(
#                         "The current layer has not been picked up as a value layer, but is configured to be a data-layer. Please check your data-layer dict."
#                     )

#                 ###########
#                 # if its a integer we just assume the current layer we hit does not have to be converted (other than to float)
#                 if isinstance(data_layer_dict[name], int):
#                     value = float(key)

#                     #################
#                     # Handle unit for delay-time
#                     if name == "delay_time":
#                         value = value * config["delay_time_default_unit"]

#                 # if its a dictionary, we have more options: convert with factor, convert with function, calculate binsize
#                 elif isinstance(data_layer_dict[name], dict):

#                     ########
#                     #
#                     value = float(key)

#                     ########
#                     # multiply by factor or apply function on value
#                     value = handle_custom_scaling_or_conversion(
#                         config=config,
#                         data_layer_or_column_dict_entry=data_layer_dict[name],
#                         value=value,
#                     )

#                     #################
#                     # Handle unit for delay-time
#                     if name == "delay_time":
#                         if "delay_time_unit" in data_layer_dict[name].keys():
#                             value = value * data_layer_dict[name]["delay_time_unit"]
#                         else:
#                             value = value * config["delay_time_default_unit"]

#                     ########
#                     # Determine binsize multiplication factor
#                     binsizes, extra_value_dict = handle_binsize_multiplication_factor(
#                         config=config,
#                         ensemble=ensemble,
#                         data_layer_dict_entry=data_layer_dict[name],
#                         name=name,
#                         key=key,
#                         key_i=key_i,
#                         binsizes=binsizes,
#                         extra_value_dict=extra_value_dict,
#                     )
#                 #
#                 else:
#                     raise ValueError("input type not supported.")

#                 # store
#                 data_dict[inverted_data_layer_dict[depth]] = value

#                 #
#                 config["logger"].debug("data dict: {}".format(data_dict))

#             #################
#             # multiply with starformation if we reached a depth that is deeper than any data layer.
#             if depth >= deepest_data_layer_depth:
#                 # multiplication with SFR-related things here
#                 ensemble[key] = ensemble_handle_SFR_multiplication(
#                     config=config,
#                     sfr_dict=sfr_dict,
#                     convolution_instruction=convolution_instruction,
#                     time_bin_info_dict=time_bin_info_dict,
#                     ensemble=ensemble[key],
#                     data_dict=data_dict,
#                     extra_value_dict=extra_value_dict,
#                     #
#                     persistent_data=persistent_data,
#                     previous_convolution_results=previous_convolution_results,
#                 )
#             else:
#                 # call self with increased depth
#                 ensemble[key] = ensemble_convolve_ensemble(
#                     sfr_dict=sfr_dict,
#                     config=config,
#                     time_bin_info_dict=time_bin_info_dict,
#                     convolution_instruction=convolution_instruction,
#                     ensemble=ensemble[key],
#                     depth=depth + 1,
#                     data_dict=data_dict,
#                     extra_value_dict=extra_value_dict,
#                     #
#                     persistent_data=persistent_data,
#                     previous_convolution_results=previous_convolution_results,
#                 )

#     elif isinstance(ensemble, ALLOWED_NUMERICAL_TYPES):
#         raise ValueError(
#             "Arrived at a layer in the ensemble that is of numerical type (depth={}), likely the endpoint. This should not happen".format(
#                 depth
#             )
#         )

#     return ensemble


# def extract_units_from_endpoints(endpoints):
#     """
#     The endpoint array possibly contains values with units.
#     """

#     units = []
#     unique_units = []

#     for endpoint_i, endpoint in enumerate(endpoints):
#         unit = endpoint.unit
#         units.append(unit)
#         if unit not in unique_units:
#             unique_units.append(unit)
#             if len(unique_units) > 1:
#                 raise ValueError(
#                     "Multiple units in the same endpoint array not supported currently"
#                 )
#         endpoints[endpoint_i] = endpoints[endpoint_i].value

#     # multiply the whole array with the single, first unit. That should be the only unit cause otherwise there would be an error
#     return np.array(endpoints) * units[0]


# def convolve_ensemble_by_integration(
#     time_bin_info_dict,
#     config,
#     convolution_instruction,
#     data_dict,
#     sfr_dict,
#     #
#     persistent_data=None,
#     previous_convolution_results=None,
# ):
#     """
#     Function for the multiprocessing worker to convolve ensemble-based data.

#     There are some requirements for the ensemble structure, but the nested
#     dictionary that is passed by data_dict needs to contain at least time or
#     log10_time

#     Moreover, the end-point nodes are expected to contain the
#     quantity-per-unit-mass. In that way we do not have to rely on extracting
#     that from the meta-data and stuff

#     Note: ensemble convolution only supports convolution by integration at this point.
#     """

#     #
#     config["logger"].debug(
#         "Convolving ensemble-based data {} for bin_center {}".format(
#             convolution_instruction["input_data_name"], time_bin_info_dict["bin_center"]
#         )
#     )

#     # ##########
#     # #
#     # config["logger"].debug("Worker {}".format(job_dict['worker_ID']) +
#     #     "Ensemble convolution by integration: {} bin center: {}: Calculating {} rates".format(
#     #         bin_type,
#     #         bin_center,
#     #         convolution_instruction["input_data_name"],
#     #     )
#     # )

#     # pre-convolution preparation
#     ensemble = data_dict["ensemble_data"]
#     data_layer_dict = convolution_instruction["data_layer_dict"]

#     # check if we want to supply a fixed metallicity
#     data_dict = {}
#     if "metallicity_value" in convolution_instruction:
#         data_dict["metallicity"] = convolution_instruction["metallicity_value"]

#     # add some extra things to the convolution instruction
#     # TODO this can be placed elsewhere?
#     # TODO: what the difference between max_depth and deepest_data_layer_depth?
#     convolution_instruction["deepest_data_layer_depth"] = get_deepest_data_layer_depth(
#         data_layer_dict=data_layer_dict
#     )
#     convolution_instruction["inverted_data_layer_dict"] = invert_data_layer_dict(
#         data_layer_dict=data_layer_dict
#     )
#     convolution_instruction["data_layer_values"] = get_data_layer_dict_values(
#         data_layer_dict=data_layer_dict
#     )

#     # convolution
#     ensemble = ensemble_convolve_ensemble(
#         ensemble=ensemble,
#         sfr_dict=sfr_dict,
#         convolution_instruction=convolution_instruction,
#         config=config,
#         time_bin_info_dict=time_bin_info_dict,
#         data_dict=data_dict,
#     )

#     # marginalisation
#     config, ensemble, convolution_instruction = ensemble_handle_marginalisation(
#         config=config,
#         ensemble=ensemble,
#         convolution_instruction=convolution_instruction,
#         is_pre_conv=False,
#     )

#     # detach endpoints from ensemble
#     stripped_ensemble, stripped_endpoints, found_units = strip_ensemble_endpoints(
#         ensemble=ensemble
#     )

#     # extract units from array (or rather, make it an array with a unit, instead of a array of values with units)
#     # if found_units:
#     stripped_endpoints = extract_units_from_endpoints(endpoints=stripped_endpoints)

#     # put back the units
#     stripped_endpoints = stripped_endpoints * config["normalized_yield_unit"]

#     #
#     convolution_result = {"yield": stripped_endpoints}

#     if time_bin_info_dict["bin_number"] == 0:
#         convolution_result["stripped_ensemble"] = stripped_ensemble

#     return {"convolution_results": convolution_result}
