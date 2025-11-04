import copy

import astropy.units as u

dimensionless_unit = u.m / u.m


def build_description_table(table_name, parameter_list, description_dict):  # DH0001
    """
    Function to create a table containing the description of the options
    """

    #
    indent = "   "

    # Get parameter list and parse descriptions
    parameter_list_with_descriptions = [
        [
            parameter,
            parse_description(
                description_dict=description_dict[parameter], add_validation=False
            ),
        ]
        for parameter in parameter_list
    ]

    # Construct parameter list
    rst_table = """
.. list-table:: {}
{}:widths: 25, 75
{}:header-rows: 1
""".format(
        table_name, indent, indent
    )

    #
    rst_table += "\n"
    rst_table += indent + "* - Option\n"
    rst_table += indent + "  - Description\n"

    for parameter_el in parameter_list_with_descriptions:
        rst_table += indent + "* - {}\n".format(parameter_el[0])
        rst_table += indent + "  - {}\n".format(parameter_el[1])

    return rst_table


def parse_description(  # DH0001
    description_dict, add_unit=True, add_value=True, add_validation=True
):
    """
    Function to parse the description for a given parameter
    """

    # Make a local copy
    description_dict = copy.copy(description_dict)

    ############
    # Add description
    description_string = "Description:\n   "

    # Clean description text
    description_text = description_dict["description"].strip()

    if description_text:
        description_text = description_text[0].capitalize() + description_text[1:]
        if description_text[-1] != ".":
            description_text = description_text + "."
    description_string += description_text

    ##############
    # Add unit (in latex)
    if add_unit and "unit" in description_dict:
        if description_dict["unit"] != dimensionless_unit:
            description_string = description_string + "\n\nUnit: [{}].".format(
                description_dict["unit"].to_string("latex_inline")
            )

    ##############
    # Add default value
    if add_value and "value" in description_dict:
        # Clean
        if isinstance(description_dict["value"], str) and (
            "/home" in description_dict["value"]
        ):
            description_dict["value"] = "example path"

        # Write
        description_string = description_string + "\n\nDefault value:\n   {}".format(
            description_dict["value"]
        )

    ##############
    # Add validation
    if add_validation and "validation" in description_dict:
        # Write
        description_string = description_string + "\n\nValidation:\n   {}".format(
            description_dict["validation"]
        )

    # Check if there are newlines, and replace them by newlines with indent
    description_string = description_string.replace("\n", "\n       ")

    return description_string


#############
# Utilities to build the description table
def write_convolution_config_and_instruction_documentation_to_rst_file(  # DH0001
    convolution_config_defaults_dict,
    convolution_instruction_defaults_dict,
    output_file: str,
) -> None:
    """
    Function that writes the descriptions of the grid options to an rst file

    Args:
        output_file: target file where the grid options descriptions are written to
    """

    ###############
    # Check input
    if not output_file.endswith(".rst"):
        msg = "Filename doesn't end with .rst, please provide a proper filename"
        raise ValueError(msg)

    ############
    # Set up intro
    page_text = ""
    title = "Convolution options"
    page_text += title + "\n"
    page_text += "=" * len(title) + "\n\n"
    page_text += "The following page contains documentation on the convolution-config options and the convolution-instructions options."
    page_text += "\n\n"
    page_text += "The convolution-config dictionary provides the `global` configuration of the convolution code. This dictionary has to be passed to the main entrypoint function as `convolve(config=convolution_config)`."
    page_text += "\n\n"
    page_text += "The convolution-instruction dictionary provides configuration on a per-convolution basis, allowing the main convolution function to perform a series of convolutions with different configurations. This dictionary has to be passed to the convolution-config function as `convolution_config['convolution_instructions'] = [{convolution_instruction_dict_1, ...}]`."
    page_text += "\n\n"

    ###########
    # set up documentation for the convolution-config

    #
    convolution_config_descriptions_dict = {}
    for key, value in convolution_config_defaults_dict.items():
        convolution_config_descriptions_dict[key] = {}
        convolution_config_descriptions_dict[key]["description"] = value["description"]
        convolution_config_descriptions_dict[key]["value"] = value["value"]

        if "validation" in value:
            convolution_config_descriptions_dict[key]["validation"] = value[
                "validation"
            ]

    #
    convolution_config_description_text = ""
    convolution_config_description_title = "Convolution-config options"
    convolution_config_description_text = convolution_config_description_title + "\n"
    convolution_config_description_text += (
        "-" * len(convolution_config_description_title) + "\n\n"
    )
    # convolution_config_description_text += "In this section we list the public options for the population code. These are meant to be changed by the user.\n"
    convolution_config_description_text += build_description_table(
        table_name="Convolution-config options",
        parameter_list=sorted(list(convolution_config_defaults_dict.keys())),
        description_dict=convolution_config_descriptions_dict,
    )
    page_text += convolution_config_description_text
    page_text += "\n\n"

    ###########
    # set up documentation for the convolution-instruction

    #
    convolution_instruction_descriptions_dict = {}
    for key, value in convolution_instruction_defaults_dict.items():
        convolution_instruction_descriptions_dict[key] = {}
        convolution_instruction_descriptions_dict[key]["description"] = value[
            "description"
        ]
        convolution_instruction_descriptions_dict[key]["value"] = value["value"]

        if "validation" in value:
            convolution_instruction_descriptions_dict[key]["validation"] = value[
                "validation"
            ]

    #
    convolution_instruction_description_text = ""
    convolution_instruction_description_title = "Convolution-instruction options"
    convolution_instruction_description_text = (
        convolution_instruction_description_title + "\n"
    )
    convolution_instruction_description_text += (
        "-" * len(convolution_instruction_description_title) + "\n\n"
    )
    # convolution_config_description_text += "In this section we list the public options for the population code. These are meant to be changed by the user.\n"
    convolution_instruction_description_text += build_description_table(
        table_name="Convolution-instruction options",
        parameter_list=sorted(list(convolution_instruction_defaults_dict.keys())),
        description_dict=convolution_instruction_descriptions_dict,
    )
    page_text += convolution_instruction_description_text
    page_text += "\n\n"

    ###############
    # write to file
    with open(output_file, "w") as f:
        f.write(page_text)
