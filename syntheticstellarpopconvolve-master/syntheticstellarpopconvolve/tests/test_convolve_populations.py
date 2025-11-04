"""
Testcases for convolve_populations file
"""

# import copy
# import os
import unittest

# from syntheticstellarpopconvolve import (
#     default_convolution_config,
#     default_convolution_instruction,
# )
# from syntheticstellarpopconvolve.check_and_prepare_output_file import (
#     check_and_prepare_output_file,
# )
# from syntheticstellarpopconvolve.check_and_update_convolution_config import (
#     check_and_update_convolution_config,
# )
from syntheticstellarpopconvolve.general_functions import (  # generate_boilerplate_outputfile,
    temp_dir,
)

# import astropy.units as u
# import numpy as np
# import pandas as pd


TMP_DIR = temp_dir(
    "tests", "tests_convolution", "tests_convolve_populations", clean_path=True
)


if __name__ == "__main__":
    unittest.main()
