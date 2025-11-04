"""
Testcases for post_convolution_hook_routines file
"""

import unittest

import numpy as np

from syntheticstellarpopconvolve import (
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.general_functions import temp_dir
from syntheticstellarpopconvolve.post_convolution_hook_routines import (  # handle_extra_weights_function,
    extract_arguments,
    handle_post_convolution_function,
)

np.random.seed(0)


TMP_DIR = temp_dir(
    "tests",
    "tests_convolution",
    "tests_post_convolution_hook_routines",
    clean_path=True,
)


class test_extract_arguments(unittest.TestCase):
    def test_extract_arguments_1_extra(self):
        def funca(a, b):
            pass

        args = extract_arguments(funca, {"a": 2, "b": 3, "c": 4})
        self.assertEqual(args, {"a": 2, "b": 3})

    def test_extract_arguments_exact(self):
        def funca(a, b):
            pass

        args = extract_arguments(funca, {"a": 2, "b": 3})
        self.assertEqual(args, {"a": 2, "b": 3})

    def test_extract_arguments_default_args_only(self):
        def funcb(a, b, c=3):
            pass

        args = extract_arguments(funcb, {"a": 2, "b": 3})
        self.assertEqual(args, {"a": 2, "b": 3})

    def test_extract_arguments_default_all(self):
        def funcb(a, b, c=3):
            pass

        args = extract_arguments(funcb, {"a": 2, "b": 3, "c": 4})
        self.assertEqual(args, {"a": 2, "b": 3, "c": 4})

    def test_extract_arguments_default_1_extra(self):
        def funcb(a, b, c=3):
            pass

        args = extract_arguments(funcb, {"a": 2, "b": 3, "c": 4, "d": 5})
        self.assertEqual(args, {"a": 2, "b": 3, "c": 4})

    def test_extract_arguments_missing(self):
        def funca(a, b):
            pass

        with self.assertRaises(ValueError):
            extract_arguments(funca, {"a": 2})


def post_convolution_function_no_convolution_result_argument():
    pass


def post_convolution_function_wrong_convolution_results_type(convolution_results):
    return ()


def post_convolution_function_wrong_elements_in_convolution_results_list(
    convolution_results,
):
    return [()]


def post_convolution_function_no_name_in_convolution_results_lists(convolution_results):
    return [{}]


class test_handle_post_convolution_function(unittest.TestCase):
    def test_handle_post_convolution_function_no_convolution_results(self):

        convolution_instruction = {
            **default_convolution_instruction,
            "post_convolution_function": post_convolution_function_no_convolution_result_argument,
        }

        with self.assertRaises(ValueError):
            handle_post_convolution_function(
                config={**default_convolution_config},
                sfr_dict={},
                data_dict={},
                time_bin_info_dict={},
                convolution_instruction=convolution_instruction,
                convolution_results={},
                name="no_convolution_results",
            )

    def test_handle_post_convolution_function_wrong_convolution_results_type(self):

        convolution_instruction = {
            **default_convolution_instruction,
            "post_convolution_function": post_convolution_function_wrong_convolution_results_type,
        }

        with self.assertRaises(ValueError):
            handle_post_convolution_function(
                config={**default_convolution_config},
                sfr_dict={},
                data_dict={},
                time_bin_info_dict={},
                convolution_instruction=convolution_instruction,
                convolution_results={},
                name="wrong_convolution_results_type",
            )

    def test_handle_post_convolution_function_wrong_elements_in_convolution_results_list(
        self,
    ):

        convolution_instruction = {
            **default_convolution_instruction,
            "post_convolution_function": post_convolution_function_wrong_elements_in_convolution_results_list,
        }

        with self.assertRaises(ValueError):
            handle_post_convolution_function(
                config={**default_convolution_config},
                sfr_dict={},
                data_dict={},
                time_bin_info_dict={},
                convolution_instruction=convolution_instruction,
                convolution_results={},
                name="wrong_elements_in_convolution_results_list",
            )

    def test_handle_post_convolution_function_no_name_in_convolution_results_lists(
        self,
    ):

        convolution_instruction = {
            **default_convolution_instruction,
            "post_convolution_function": post_convolution_function_no_name_in_convolution_results_lists,
        }

        with self.assertRaises(ValueError):
            handle_post_convolution_function(
                config={**default_convolution_config},
                sfr_dict={},
                data_dict={},
                time_bin_info_dict={},
                convolution_instruction=convolution_instruction,
                convolution_results={},
                name="no_name_in_convolution_results_lists",
            )


if __name__ == "__main__":
    unittest.main()
