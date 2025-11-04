"""
Functions for on-the-fly convolutions.

This is mostly experimental, but in short:
- based on the total mass formed into stars in a given target convolution bin, and potentially the metallicity distribution,
- the user can provide a call to a population-synthesis code that evolves a population on the fly.

The user is responsible for all the conversions and reweighting here.

what do we need:
- on-the-fly function to convolve

"""

from syntheticstellarpopconvolve.general_functions import is_mass_unit
from syntheticstellarpopconvolve.post_convolution_hook_routines import (
    extract_arguments,
    handle_post_convolution_function,
)


def convolve_on_the_fly_post_convolution_hook_wrapper(
    config,
    sfr_dict,
    convolution_instruction,
    time_bin_info_dict,
    convolution_results,
    #
    persistent_data=None,
    previous_convolution_results=None,
):
    """
    Function to wrap the post-convolution function call for event-convolution by integration.

    rules:
    - additional data can be added to the convolution_results
    - the number of systems has to be equal to before the post-convolution function.
    """

    #
    name = "convolve on-the-fly"

    #
    config["logger"].warning(
        "Handling post-convolution function hook call for {}".format(name)
    )

    #############
    # call hook
    convolution_results = handle_post_convolution_function(
        config=config,
        sfr_dict=sfr_dict,
        data_dict={},
        time_bin_info_dict=time_bin_info_dict,
        convolution_instruction=convolution_instruction,
        convolution_results=convolution_results,
        name=name,
        #
        persistent_data=persistent_data,
        previous_convolution_results=previous_convolution_results,
    )

    return convolution_results


def handle_call_on_the_fly_function(
    config, time_bin_info_dict, sfr_dict, convolution_instruction
):
    """
    Function to call an external evolution code to perform on-the-fly evolution
    """

    ######
    # Get quantities
    bin_number = time_bin_info_dict["bin_number"]
    bin_size = time_bin_info_dict["bin_size"]
    bin_lower_edge = time_bin_info_dict["bin_edge_lower"]
    sfr = sfr_dict["starformation_rate_array"][bin_number]
    total_star_formation_in_bin = sfr * bin_size

    # Check if the total star formation is a mass-type value
    if not is_mass_unit(total_star_formation_in_bin):
        raise ValueError(
            "The total star formation in current bin ({}) is not of a mass-type unit. Something wrong with either the sfr ({}) or the time-bin size ({})".format(
                total_star_formation_in_bin, sfr, bin_size
            )
        )

    # ##################
    # # check whether the yield is dimensionless

    # # it has to be dimensionless, otherwise its not really a count.
    # # force into cgs (basically to ensure that Gyr/yr is seen as dimensionless with a scale)
    # if has_unit(yield_array.cgs, fail_on_dimensionless=True):
    #     raise ValueError(
    #         "Combined formation yield (unit: {}. dimension: {}) has to be dimensionless for convolution by sampling. The total star formation in bin (unit: {}. dimension: {}) times the normalized yield (unit: {}. dimension: {}) should not have a unit anymore.".format(
    #             yield_array.unit.to_string(),
    #             get_physical_dimensions(yield_array.unit),
    #             starformation.unit.to_string(),
    #             get_physical_dimensions(starformation.unit),
    #             normalized_yield_unit.unit.to_string(),
    #             get_physical_dimensions(normalized_yield_unit.unit),
    #         )
    #     )

    #
    config["logger"].warning(
        "Lower time bin {} upper time bin {} total mass formed {}".format(
            bin_lower_edge,
            bin_lower_edge + bin_size,
            total_star_formation_in_bin,
        )
    )

    #
    metallicity_distribution = (
        sfr_dict["metallicity_distribution_array"][:, bin_number]
        if "metallicity_distribution_array" in sfr_dict
        else None
    )

    ######
    # Call user-provided function including sfr info
    on_the_fly_function = convolution_instruction.get("on_the_fly_function", None)

    if on_the_fly_function is None:
        raise ValueError(
            "Can't perform on-the-fly convolution if no `on_the_fly_function` is provided. Please add a function to the `on_the_fly_function` field in the `convolution_instruction` dict"
        )

    # Construct what parameters are available for the extra function
    available_parameters = {
        # Standard info
        "config": config,
        "sfr_dict": sfr_dict,
        "time_bin_info_dict": time_bin_info_dict,
        "convolution_instruction": convolution_instruction,
        # Explicit info
        "total_star_formation_in_bin": total_star_formation_in_bin,
        "metallicity_distribution": metallicity_distribution,
        **convolution_instruction.get("on_the_fly_function_extra_parameters", {}),
    }

    # Extract the correct things from the available parameters
    on_the_fly_function_args = extract_arguments(
        func=on_the_fly_function,
        arg_dict=available_parameters,
    )

    # Enforce that certain arguments are present:
    if "total_star_formation_in_bin" not in on_the_fly_function_args:
        raise ValueError(
            "`total_star_formation_in_bin` is a required argument in the `on_the_fly_function` call."
        )

    if "metallicity_distribution_array" in sfr_dict:
        if "metallicity_distribution" not in on_the_fly_function_args:
            raise ValueError(
                "`metallicity_distribution` is a required argument in the `on_the_fly_function` call when including metallicity information in the starformation rate dict"
            )

    #
    config["logger"].debug(
        "Handling `on_the_fly_function` function call using function {} and arguments {}".format(
            convolution_instruction["on_the_fly_function"].__name__,
            on_the_fly_function,
        )
    )

    # Call on-the-fly function
    convolution_results = on_the_fly_function(**on_the_fly_function_args)

    # check result type
    if not isinstance(convolution_results, dict):
        raise ValueError(
            "The object-type returned by `on_the_fly_function` ({}) must be of a dictionary type.".format(
                type(convolution_results)
            )
        )

    return convolution_results


def convolve_on_the_fly(
    config,
    sfr_dict,
    convolution_instruction,
    time_bin_info_dict,
    #
    persistent_data=None,
    previous_convolution_results=None,
):
    """ """

    #
    config["logger"].warning(
        "Performing on-the-fly convolution at bin-center {}".format(
            time_bin_info_dict["bin_center"]
        )
    )

    ######
    # Handle calling
    convolution_results = handle_call_on_the_fly_function(
        config=config,
        time_bin_info_dict=time_bin_info_dict,
        sfr_dict=sfr_dict,
        convolution_instruction=convolution_instruction,
    )

    ######
    # Handle post-convolution function
    convolution_results = convolve_on_the_fly_post_convolution_hook_wrapper(
        config=config,
        sfr_dict=sfr_dict,
        convolution_instruction=convolution_instruction,
        time_bin_info_dict=time_bin_info_dict,
        convolution_results=convolution_results,
        #
        persistent_data=persistent_data,
        previous_convolution_results=previous_convolution_results,
    )

    return {"convolution_results": convolution_results}
