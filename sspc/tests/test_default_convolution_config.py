"""
Testcases for default_convolution_config.
"""

import logging
import os
import unittest

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve.default_convolution_config import (
    array_validation,
    callable_or_none_validation,
    callable_validation,
    dict_or_list_of_dicts_validation,
    existing_path_validation,
    list_of_dicts_validation,
    logger_validation,
    unit_validation,
)
from syntheticstellarpopconvolve.general_functions import temp_dir

TMP_DIR = temp_dir(
    "tests", "tests_convolution", "tests_default_convolution_config", clean_path=True
)


def assertDoesNotRaise(self, exception, expr, *args, **kwargs):

    # Use try-except block to catch exceptions
    try:
        expr(*args, **kwargs)
    except exception as e:
        self.fail(f"Unexpected exception raised: {e}")


unittest.TestCase.assertDoesNotRaise = assertDoesNotRaise


class test_list_of_dicts_validation(unittest.TestCase):
    def test_list_of_dicts_validation(self):

        # should raise a valueerror if the object we pass is not a logging object
        self.assertRaises(ValueError, list_of_dicts_validation, [[]])

    def test_list_of_dicts_validation_no_raise(self):
        self.assertDoesNotRaise(ValueError, list_of_dicts_validation, [{}])


class test_dict_or_list_of_dicts_validation(unittest.TestCase):
    def test_dict_or_list_of_dicts_validation(self):

        # should raise a valueerror if the object we pass is not a logging object
        self.assertRaises(ValueError, dict_or_list_of_dicts_validation, [[]])

    def test_dict_or_list_of_dicts_validation_no_raise_dict(self):
        self.assertDoesNotRaise(ValueError, dict_or_list_of_dicts_validation, {})

    def test_dict_or_list_of_dicts_validation_no_raise_list_of_dicts(self):
        self.assertDoesNotRaise(ValueError, dict_or_list_of_dicts_validation, [{}])


class test_unit_validation(unittest.TestCase):
    def test_unit_validation(self):

        # should raise a valueerror if the object we pass is not a logging object
        self.assertRaises(ValueError, unit_validation, 1)

    def test_unit_validation_no_raise(self):
        self.assertDoesNotRaise(ValueError, unit_validation, u.yr)


class test_logger_validation(unittest.TestCase):
    def test_logger_validation(self):

        # should raise a valueerror if the object we pass is not a logging object
        self.assertRaises(ValueError, logger_validation, {})

    def test_logger_validation_no_raise(self):
        self.assertDoesNotRaise(
            ValueError, logger_validation, logging.getLogger(__name__)
        )


class test_callable_validation(unittest.TestCase):
    """ """

    def test_callable_validation(self):
        # should raise a valueerror if the object we pass is not a logging object
        self.assertRaises(ValueError, callable_validation, {})

    def test_callable_validation_no_raise(self):
        self.assertDoesNotRaise(ValueError, callable_validation, callable_validation)


class test_callable_or_none_validation(unittest.TestCase):
    def test_callable_or_none_validation(self):
        # should raise a valueerror if the object we pass is not a logging object
        self.assertRaises(ValueError, callable_or_none_validation, {})

    def test_callable_or_none_validation_no_raise(self):
        self.assertDoesNotRaise(
            ValueError, callable_or_none_validation, callable_validation
        )

        self.assertDoesNotRaise(ValueError, callable_or_none_validation, None)


class test_array_validation(unittest.TestCase):
    def test_array_validation(self):
        # should raise a valueerror if the object we pass is not a logging object
        self.assertRaises(ValueError, array_validation, {})

    def test_array_validation_no_raise(self):
        self.assertDoesNotRaise(ValueError, array_validation, np.array([]))


class test_existing_path_validation(unittest.TestCase):
    def test_existing_path_validation(self):
        # should raise a valueerror if the object we pass is not a logging object
        self.assertRaises(ValueError, existing_path_validation, {})

    def test_existing_path_validation_no_raise(self):
        targetfile = os.path.join(TMP_DIR, "test.txt")
        with open(targetfile, "w") as f:
            f.write("hello")

        self.assertDoesNotRaise(ValueError, existing_path_validation, targetfile)


if __name__ == "__main__":
    unittest.main()
