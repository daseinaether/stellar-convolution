"""
Configuration file for the Sphinx documentation builder.

This file only contains a selection of the most common options. For a full
list see the documentation:
    https://www.sphinx-doc.org/en/master/usage/configuration.html
    https://brendanhasz.github.io/2019/01/05/sphinx.html
    https://www.sphinx-doc.org/en/1.5/ext/example_google.html


https://stackoverflow.com/questions/22256995/restructuredtext-in-sphinx-and-docstrings-single-vs-double-back-quotes-or-back

TODO: https://stackoverflow.com/questions/64443832/sorting-table-with-rst-and-python-sphinx-in-html-output implement datatables

"""

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import datetime
import os
import sys

from git import Repo

from syntheticstellarpopconvolve.default_convolution_config import (
    default_convolution_config_dict,
)
from syntheticstellarpopconvolve.default_convolution_instruction import (
    default_convolution_instruction_dict,
)
from syntheticstellarpopconvolve.docs import (
    write_convolution_config_and_instruction_documentation_to_rst_file,
)


def write_custom_footer():
    """
    Function to write the custom footer to the template file
    """

    #
    branch_infix = "/-/tree/"
    # commit_infix = "/-/commit/"

    ############
    # Construct binary_c-python git information
    base_dir = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )
    local_repo = Repo(path=base_dir)

    sspc_git_branch_name = str(local_repo.active_branch.name)
    # sspc_git_revision = str(local_repo.active_branch.commit)

    #
    sspc_git_root = "https://gitlab.com/dhendriks/syntheticstellarpopconvolve"
    sspc_branch_url = sspc_git_root + branch_infix + sspc_git_branch_name
    # sspc_commit_url = (
    #     sspc_git_root + commit_infix + sspc_git_revision.replace('"', "").split(":")[-1]
    # )

    ############
    # Construct footer text

    string = """
<br><br>
Generated on Synthetic Stellar Pop Convolve branch {sspc_git_branch_name}: <a href="{sspc_branch_url}">git branch url</a>.
""".format(
        sspc_git_branch_name=sspc_git_branch_name,
        sspc_branch_url=sspc_branch_url,
    )

    # if not os.getenv("READTHEDOCS", False):
    #     string += ' and <a href="{sspc_commit_url}">git commit url</a>'.format(
    #         sspc_commit_url=sspc_commit_url,
    #     )

    # Set up template
    output_text = """
{{% extends '!footer.html' %}}

{{%- block extrafooter %}}
{footerstring}
{{% endblock %}}
"""

    # format
    formatted_text = output_text.format(footerstring=string).strip()

    # Write to file
    with open("_templates/footer.html", "w") as outfile_filehandle:
        outfile_filehandle.write(formatted_text)


# #
# def patched_m2r2_setup(app):
#     """
#     Function to handle the markdown parsing better
#     """

#     try:
#         return current_m2r2_setup(app)
#     except AttributeError:
#         app.add_source_suffix(".md", "markdown")
#         app.add_source_parser(m2r2.M2RParser)
#     return dict(
#         version=m2r2.__version__,
#         parallel_read_safe=True,
#         parallel_write_safe=True,
#     )


# Include paths for python code
sys.path.insert(0, os.path.abspath("."))

# include paths for c code
cautodoc_root = os.path.abspath("../../")

# -- Project information -----------------------------------------------------

project = "Synthetic Stellar Pop Convolve (SSPC)"
copyright = "{}, David Hendriks".format(datetime.datetime.now().year)
author = "David Hendriks"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "hawkmoth",
    # "m2r2",
    "myst_parser",
    "sphinx_rtd_theme",
    "sphinx_autodoc_typehints",  # https://mypy.readthedocs.io/en/stable/cheat_sheet_py3.html
    "nbsphinx",
    "sphinx_math_dollar",
    "sphinx.ext.mathjax",
]

mathjax3_config = {
    "tex": {
        "inlineMath": [["\\(", "\\)"]],
        "displayMath": [["\\[", "\\]"]],
    }
}

# Napoleon settings
napoleon_google_docstring = (
    True  # https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
)
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# source_suffix = [".rst", ".md"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

# html_theme = "alabaster"
html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    "css/custom.css",
]


#######
# logo
html_logo = "_static/logo/sspc_logo.png"


# """Patching m2r2"""
# current_m2r2_setup = m2r2.setup

# #
# m2r2.setup = patched_m2r2_setup

print("Generating convolution_documentation.rst")
write_convolution_config_and_instruction_documentation_to_rst_file(
    convolution_config_defaults_dict=default_convolution_config_dict,
    convolution_instruction_defaults_dict=default_convolution_instruction_dict,
    output_file="convolution_documentation.rst",
)
print("Done")

# Generate a custom footer
if not os.getenv("READTHEDOCS", False):
    print("Generating custom footer")
    write_custom_footer()
    print("Done")
