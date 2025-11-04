"""
Testcases for convolution_sfr_distributions file
"""

import unittest

from syntheticstellarpopconvolve.general_functions import temp_dir

TMP_DIR = temp_dir(
    "tests",
    "tests_convolution",
    "tests_starformation_rate_distributions",
    clean_path=True,
)

if __name__ == "__main__":
    unittest.main()
