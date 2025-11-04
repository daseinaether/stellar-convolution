"""
File containing methods to support post-convolution hook functionality
"""

import inspect


def extract_arguments(func, arg_dict):
    """
    Function that extracts the entries in 'arg_dict' that are arguments to the function 'func'
    """

    # get various arg types
    signature = inspect.signature(func)
    all_args = inspect.getfullargspec(func).args
    args_with_defaults = [
        k
        for k, v in signature.parameters.items()
        if v.default is not inspect.Parameter.empty
    ]
    args_without_defaults = [arg for arg in all_args if arg not in args_with_defaults]

    for arg in args_without_defaults:
        if arg not in arg_dict.keys():
            raise ValueError(
                "Argument '{}' of postconvolution function is not part of the available information. Please only choose arguments from the following list: '{}'".format(
                    arg, list(arg_dict.keys())
                )
            )

    # construct args
    args = {arg: arg_dict[arg] for arg in args_without_defaults}

    # check if kwonlyargs are also passed along
    args_for_args_with_defaults = {
        arg: arg_dict[arg] for arg in args_with_defaults if arg in arg_dict.keys()
    }

    # combine args
    combined_args = {**args, **args_for_args_with_defaults}

    return combined_args


def handle_post_convolution_function(
    config,
    sfr_dict,
    data_dict,
    time_bin_info_dict,
    convolution_instruction,
    convolution_results,
    name,
    #
    persistent_data=None,
    previous_convolution_results=None,
):
    """
    Function to handle post-convolution function call.

    An example of a post-convolution call is integrating systems to
    present-day time with LegWork and filtering out systems that do
    not fall within the LISA frequency range or that have merged by
    the present-day.

    Another example is to integrate systems through a gravitational
    potential based on the sampled position and a certain integration
    time.
    """

    post_convolution_function = convolution_instruction.get(
        "post_convolution_function", None
    )

    if post_convolution_function is not None:

        # Construct what parameters are available for the extra function
        available_parameters = {
            "config": config,
            "sfr_dict": sfr_dict,
            "data_dict": data_dict,
            "convolution_results": convolution_results,
            "time_bin_info_dict": time_bin_info_dict,
            "convolution_instruction": convolution_instruction,
            "persistent_data": persistent_data,
            "previous_convolution_results": previous_convolution_results,
            **convolution_instruction.get(
                "post_convolution_function_extra_parameters", {}
            ),
        }

        # Make sure we extract the correct things from the available parameters
        post_convolution_function_args = extract_arguments(
            func=post_convolution_function,
            arg_dict=available_parameters,
        )

        # Enforce that certain arguments are present:
        if "convolution_results" not in post_convolution_function_args:
            raise ValueError(
                "`convolution_results` is a required argument in the `post_convolution_function` call."
            )

        #
        config["logger"].debug(
            "Handling '{}' post-convolution function call using function {} and arguments {}".format(
                name,
                convolution_instruction["post_convolution_function"].__name__,
                post_convolution_function_args,
            )
        )

        # Call post-convolution function
        convolution_results = post_convolution_function(
            **post_convolution_function_args
        )

        ################
        # Check shape/type of results

        # check if the result is a list
        if isinstance(convolution_results, list):
            for convolution_result in convolution_results:
                # check if the elements are dicts
                if not isinstance(convolution_result, dict):
                    raise ValueError(
                        "The result dict object must be a dictionary type object after the post-convolution call. It's now a {}-type object".format(
                            type(convolution_result)
                        )
                    )

                # check if a name is provided to the convolution result
                if "name" not in convolution_result.keys():
                    raise ValueError(
                        "When returning multiple result-dicts, the result-dicts need to have a 'name' entry to identify and store them correctly. Please provide one."
                    )

        # Otherwise check if the convolution_results is a dict object
        elif not isinstance(convolution_results, dict):
            raise ValueError(
                "The result dict object must be a dictionary type object after the post-convolution call. It's now a {}-type object. Please ensure that you return a dictionary.".format(
                    type(convolution_results)
                )
            )

    return convolution_results
