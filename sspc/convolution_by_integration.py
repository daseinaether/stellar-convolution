# """
# Functions to convolve events
# """

from syntheticstellarpopconvolve.post_convolution_hook_routines import (
    handle_post_convolution_function,
)


def convolution_by_integration_post_convolution_hook_wrapper(
    config,
    sfr_dict,
    data_dict,
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
    name = "convolution by integration"

    #
    config["logger"].warning(
        "Handling post-convolution function hook call for {}".format(name)
    )

    #############
    # pre-call setup
    num_systems_before = len(convolution_results[list(convolution_results.keys())[0]])

    #############
    # call hook
    convolution_results = handle_post_convolution_function(
        config=config,
        sfr_dict=sfr_dict,
        data_dict=data_dict,
        time_bin_info_dict=time_bin_info_dict,
        convolution_instruction=convolution_instruction,
        convolution_results=convolution_results,
        name=name,
        #
        persistent_data=persistent_data,
        previous_convolution_results=previous_convolution_results,
    )

    #############
    # check output
    if isinstance(convolution_results, list):
        for convolution_result in convolution_results:

            #############
            # check output
            num_systems_after = len(
                convolution_result[list(convolution_result.keys())[0]]
            )

            #
            if num_systems_before != num_systems_after:
                raise ValueError(
                    "post-convolution function for event-convolution by integration has changed the number of systems stored in the output dict. Due to current data structure decisions this is not supported. Please make sure that the number of systems before and after calling this function stays equal."
                )
    else:
        #############
        # check output
        num_systems_after = len(
            convolution_results[list(convolution_results.keys())[0]]
        )

        #
        if num_systems_before != num_systems_after:
            raise ValueError(
                "post-convolution function for event-convolution by integration has changed the number of systems stored in the output dict. Due to current data structure decisions this is not supported. Please make sure that the number of systems before and after calling this function stays equal."
            )

    return convolution_results
