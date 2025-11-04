"""
Functions to check and update the convolution instructions
"""

import voluptuous as vol

from syntheticstellarpopconvolve.default_convolution_instruction import (
    default_convolution_instruction_dict,
)
from syntheticstellarpopconvolve.general_functions import check_required, is_time_unit


def check_metallicity(convolution_config, convolution_instruction):
    """
    Function to check the metallicity
    """

    #
    if isinstance(convolution_config["SFR_info"], dict):
        requires_metallicity_data = any(
            [key.startswith("metallicity") for key in convolution_config["SFR_info"]]
        )
    elif isinstance(convolution_config["SFR_info"], list):
        requires_metallicity_data = any(
            [
                any([key.startswith("metallicity") for key in sfr_dict])
                for sfr_dict in convolution_config["SFR_info"]
            ]
        )

    if requires_metallicity_data:
        if "metallicity" not in convolution_instruction["data_column_dict"].keys():
            raise ValueError(
                "If metallicity information is provided in the starformation dictionary, then you need to provide metallicity info in the 'data_column_dict' too."
            )


def check_delay_time_data_bin_info_dict(delay_time_data_bin_info_dict):
    """
    Function to check data time bin info dict
    """

    if "delay_time_data_bin_edges" not in delay_time_data_bin_info_dict:
        raise ValueError(
            "`delay_time_data_bin_edges` is required in the delay_time_data_bin_info_dict when convolving binned data"
        )

    if not is_time_unit(delay_time_data_bin_info_dict["delay_time_data_bin_edges"]):
        raise ValueError("Please express 'delay_time_data_bin_edges' in units of time")


def check_convolution_instruction(convolution_instruction, convolution_config):
    """
    Function to check convolution instructions
    """

    ##########
    # from the main dictionary, create a validation scheme
    validation_dict = {
        key: value["validation"]
        for key, value in default_convolution_instruction_dict.items()
        if "validation" in value
    }
    validation_schema = vol.Schema(validation_dict, extra=vol.ALLOW_EXTRA)

    ##########
    # do the basic validation
    for parameter, parameter_dict in convolution_config.items():

        ##########
        # Custom rules. we can decide to skip checking the input on some occasions

        #
        validation_schema({parameter: parameter_dict})

    #######
    # required for all
    check_required(
        config=convolution_instruction,
        required_list=["input_data_name", "output_data_name"],
    )

    ###################
    # checks for particular types of configurations
    if convolution_instruction["convolution_type"] == "integrate":

        check_required(
            config=convolution_instruction,
            required_list=[
                "data_column_dict",
            ],
        )

        #
        check_required(
            config=convolution_instruction["data_column_dict"],
            required_list=[
                "normalized_yield",
            ],
        )

        # check how metallicity is treated
        check_metallicity(
            convolution_config=convolution_config,
            convolution_instruction=convolution_instruction,
        )

        check_required(
            config=convolution_instruction,
            required_list=[
                "contains_binned_data",
            ],
        )

        #
        if convolution_instruction["contains_binned_data"]:
            check_required(
                config=convolution_instruction,
                required_list=[
                    "delay_time_data_bin_info_dict",
                ],
            )

            check_delay_time_data_bin_info_dict(
                delay_time_data_bin_info_dict=convolution_instruction[
                    "delay_time_data_bin_info_dict"
                ],
            )

        else:
            check_required(
                config=convolution_instruction["data_column_dict"],
                required_list=[
                    "delay_time",
                ],
            )

    elif convolution_instruction["convolution_type"] == "sample":

        check_required(
            config=convolution_instruction,
            required_list=[
                "data_column_dict",
            ],
        )

        #
        check_required(
            config=convolution_instruction["data_column_dict"],
            required_list=[
                "normalized_yield",
            ],
        )

    elif convolution_instruction["convolution_type"] == "on-the-fly":

        check_required(
            config=convolution_instruction,
            required_list=["on_the_fly_function"],
        )

    else:
        raise ValueError(
            "convolution type {} unsupported".format(
                convolution_instruction["convolution_type"]
            )
        )


def check_and_update_convolution_instruction(  # DH0001
    convolution_instruction, convolution_config
):
    """
    Function to check convolution instructions
    """

    # check
    check_convolution_instruction(
        convolution_instruction=convolution_instruction,
        convolution_config=convolution_config,
    )

    # TODO: add call to update convolution instruction


def check_and_update_convolution_instructions(convolution_config):
    """
    Main function to check the convolution instructions.
    """

    if convolution_config["convolution_instructions"]:
        for convolution_instruction in convolution_config["convolution_instructions"]:
            check_and_update_convolution_instruction(
                convolution_instruction=convolution_instruction,
                convolution_config=convolution_config,
            )
    else:
        raise ValueError("Please provide at least one convolution intruction")
