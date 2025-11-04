"""
Testcases for check_and_prepare_output_file file
"""

import copy
import os
import unittest

import h5py

from syntheticstellarpopconvolve import default_convolution_config
from syntheticstellarpopconvolve.check_and_prepare_output_file import (
    check_and_prepare_output_file,
)
from syntheticstellarpopconvolve.general_functions import temp_dir

TMP_DIR = temp_dir(
    "tests", "tests_convolution", "test_check_and_prepare_output_file", clean_path=True
)


class test_check_and_prepare_output_file(unittest.TestCase):
    """ """

    def setUp(self):
        """
        create input files
        """

        self.config_working = copy.copy(default_convolution_config)
        self.config_working["input_filename"] = os.path.join(
            TMP_DIR, "working_file.hdf5"
        )
        self.config_working["output_filename"] = os.path.join(
            TMP_DIR, "output_file.hdf5"
        )

        #
        self.working_hdf5_file = h5py.File(self.config_working["output_filename"], "w")
        self.working_hdf5_file.create_group("input_data")
        self.working_hdf5_file.create_group("config")
        self.working_hdf5_file.close()

    def test_check_and_prepare_output_file(self):

        # run preparing of output file
        check_and_prepare_output_file(self.config_working)

        #
        self.assertTrue(os.path.isfile(self.config_working["output_filename"]))

        with h5py.File(self.config_working["output_filename"], "r") as output_hdf5file:

            self.assertTrue("input_data" in output_hdf5file.keys())
            self.assertTrue("config" in output_hdf5file.keys())
            self.assertTrue("convolution" in output_hdf5file["config"].keys())


# TODO: integrate the stuff below
# """
# Testcases for check_input_file
# """

# import copy
# import os
# import sys
# import unittest

# import h5py

# from syntheticstellarpopconvolve import default_convolution_config
# from syntheticstellarpopconvolve.check_input_file import check_input_file
# from syntheticstellarpopconvolve.general_functions import temp_dir

# TMP_DIR = temp_dir(
#     "tests", "tests_convolution", "tests_check_input_file", clean_path=True
# )


# def assertMayRaise(self, exception, expr, **kwargs):
#     if exception is None:
#         try:
#             expr(kwargs)
#         except:
#             info = sys.exc_info()
#             self.fail("%s raised" % repr(info[0]))
#     else:
#         self.assertRaises(exception, expr, **kwargs)


# unittest.TestCase.assertMayRaise = assertMayRaise


# class test_check_input_file(unittest.TestCase):
#     """ """

#     def setUp(self):
#         """
#         create input files
#         """

#         self.config_working = copy.copy(default_convolution_config)
#         self.config_working["input_filename"] = os.path.join(
#             TMP_DIR, "working_file.hdf5"
#         )

#         self.config_not_working = copy.copy(default_convolution_config)
#         self.config_not_working["input_filename"] = os.path.join(
#             TMP_DIR, "not_working_file.hdf5"
#         )

#         #
#         self.working_hdf5_file = h5py.File(self.config_working["input_filename"], "w")
#         self.working_hdf5_file.create_group("input_data")
#         self.working_hdf5_file.create_group("config")
#         self.working_hdf5_file.create_group("config/population")
#         self.working_hdf5_file.close()

#         #
#         self.not_working_hdf5_file = h5py.File(
#             self.config_not_working["input_filename"], "w"
#         )
#         self.not_working_hdf5_file.close()

#     def test_check_input_file_not_working(self):
#         self.assertRaises(ValueError, check_input_file, self.config_not_working)

#     def test_check_input_file_working(self):
#         check_input_file(self.config_working)


# if __name__ == "__main__":
#     unittest.main()


if __name__ == "__main__":
    unittest.main()
