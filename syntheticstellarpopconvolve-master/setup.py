"""
Setup script for Synthetic stellar pop convolve (SSPC)
"""

import os
from distutils.core import setup

#
this_file = os.path.abspath(__file__)
this_file_dir = os.path.dirname(this_file)

############################################################
# Defining functionality
############################################################


# Functions
def version():
    """
    opens VERSION and returns version number
    """

    # with open("VERSION") as file:
    #     return file.read().strip()

    about = {}
    with open("syntheticstellarpopconvolve/_version.py") as f:
        exec(f.read(), about)

    return about["__version__"]


def readme():
    """
    Opens readme file and returns content
    """

    with open("README.md") as file:
        return file.read()


def license():
    """
    Opens license file and returns the content
    """

    with open("LICENSE.md") as file:
        return file.read()


def requirements(directory):
    """
    Opens requirements.txt and returns content as a list
    """

    requirements_file = os.path.join(directory, "requirements.txt")

    # Read out file and construct list
    requirements_list = []
    with open(requirements_file) as f:
        for el in f.readlines():
            requirements_list.append(el.strip())

    return requirements_list


############################################################
# Main setup function call
############################################################

setup(
    name="syntheticstellarpopconvolve",
    version=version(),
    description="synthetic stellar pop convolve is a python package to convolve output of stellar population synthesis codes with star formation rates.",
    author="David Hendriks",
    author_email="davidhendriks93@gmail.com",
    long_description=readme(),
    long_description_content_type="text/markdown",
    url="https://gitlab.com/dhendriks/syntheticstellarpopconvolve",
    project_urls={
        "Documentation": "https://synthetic-stellar-pop-convolve.readthedocs.io/en/latest/",
        "Source": "https://gitlab.com/dhendriks/syntheticstellarpopconvolve",
        "Bug Tracker": "https://gitlab.com/dhendriks/syntheticstellarpopconvolve/-/issues",
    },
    license="gpl",
    keywords=[
        "astrophysics",
        "stellar evolution",
        "population synthesis",
        "starformation rate history",
        "convolution",
    ],  # Keywords that define your package best
    packages=[
        "syntheticstellarpopconvolve",
        "syntheticstellarpopconvolve.tests",
        "syntheticstellarpopconvolve.usecase_notebook_utils",
        "syntheticstellarpopconvolve.example_data",
    ],
    # package_data={
    #     "syntheticstellarpopconvolve": [
    #         "example_data/*.dat",
    #         "example_data/*.json",
    #     ],
    #     # "convolution": [
    #     #     "example_data",
    #     #     "example_data",
    #     # ],
    # },
    # setup_requires=['pbr'],
    # pbr=True,
    include_package_data=True,
    install_requires=requirements(this_file_dir),
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: C",
        "Topic :: Education",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    # cmdclass={"build": CustomBuildCommand},
)
