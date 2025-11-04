"""
Testcases for convolve file
"""

import unittest

from syntheticstellarpopconvolve.general_functions import temp_dir

TMP_DIR = temp_dir("tests", "tests_convolution", "tests_convolve", clean_path=True)

if __name__ == "__main__":
    unittest.main()
