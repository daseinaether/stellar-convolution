"""
Testcases for ensemble_utils.py file
"""

import unittest

import numpy as np
import pandas as pd

from syntheticstellarpopconvolve.ensemble_utils import (
    convert_ensemble_to_dataframe,
    find_columnames_recursively,
    flatten_data_ensemble1d,
    inflate_ensemble_with_lists_and_named_layers,
    inflate_ensemble_with_lists_without_named_layers,
)
from syntheticstellarpopconvolve.general_functions import temp_dir

TMP_DIR = temp_dir(
    "tests", "tests_convolution", "tests_ensemble_utils", clean_path=True
)


example_ensemble_coherent = {
    "a": {
        "1": {
            "b": {
                "1": {"c": {"1": 1, "2": 2}},
                "2": {
                    "c": {
                        "2": 3,
                        "3": 4,
                    }
                },
            }
        },
        "2": {
            "b": {
                "2": {
                    "c": {
                        "3": 5,
                        "4": 6,
                    }
                },
                "3": {"c": {"4": 7, "5": 8}},
            }
        },
    }
}


example_ensemble_small = {
    "a": {
        "1": {
            "b": {
                "1": 1,
                "2": 2,
            }
        },
        "2": {
            "b": {
                "2": 3,
                "3": 4,
            }
        },
    }
}


class test_find_columnames_recursively(unittest.TestCase):
    def test_find_columnames_recursively(self):
        columnnames = find_columnames_recursively(example_ensemble_coherent)
        expected_columnnames = ["a", "b", "c"]
        self.assertEqual(columnnames, expected_columnnames)


class test_convert_ensemble_to_dataframe(unittest.TestCase):
    def test_convert_ensemble_to_dataframe(self):
        df = convert_ensemble_to_dataframe(example_ensemble_small)
        df["a"] = df["a"].astype(int)
        df["b"] = df["b"].astype(int)
        df["probability"] = df["probability"].astype(float)

        records = [
            {"a": 1, "b": 1, "probability": 1.0},
            {"a": 1, "b": 2, "probability": 2.0},
            {"a": 2, "b": 2, "probability": 3.0},
            {"a": 2, "b": 3, "probability": 4.0},
        ]

        expected_df = pd.DataFrame.from_records(records)
        self.assertTrue(df.equals(expected_df))

    def test_convert_ensemble_to_dataframe_structure(self):
        df = convert_ensemble_to_dataframe({"test": example_ensemble_small})

        df["a"] = df["a"].astype(int)
        df["b"] = df["b"].astype(int)
        df["probability"] = df["probability"].astype(float)

        records = [
            {"ensemble": "test", "a": 1, "b": 1, "probability": 1.0},
            {"ensemble": "test", "a": 1, "b": 2, "probability": 2.0},
            {"ensemble": "test", "a": 2, "b": 2, "probability": 3.0},
            {"ensemble": "test", "a": 2, "b": 3, "probability": 4.0},
        ]

        expected_df = pd.DataFrame.from_records(records)
        self.assertTrue(df.equals(expected_df))


class test_flatten_data_ensemble1d(unittest.TestCase):
    def test_flatten_data_ensemble1d_no_named_key(self):
        input_file = {1: 1, 2: 2}

        flattened_data_ensemble1d = flatten_data_ensemble1d(input_file)

        np.testing.assert_array_almost_equal(
            flattened_data_ensemble1d, np.array([[1.0, 2.0], [1.0, 2.0]])
        )

    def test_flatten_data_ensemble1d(self):
        input_file = {"a": {1: 1, 2: 2}}

        flattened_data_ensemble1d = flatten_data_ensemble1d(input_file, "a")

        np.testing.assert_array_almost_equal(
            flattened_data_ensemble1d, np.array([[1.0, 2.0], [1.0, 2.0]])
        )


class test_inflate_ensemble_with_lists_and_named_layers(unittest.TestCase):
    def setUp(self):
        self.example_ensemble_data = {
            "a": {
                "5": {
                    "b": {
                        "1": 0.5,
                        "2": 0.5,
                    }
                },
                "6": {
                    "b": {
                        "3": 0.5,
                        "4": 0.5,
                    }
                },
                "7": {
                    "b": {
                        "8": 0.5,
                        "9": 0.5,
                    }
                },
            }
        }

        # self.example_ensemble_data_2 = {
        #     "ensemble_data": {
        #         "a": {
        #             "5": {
        #                 "b": {
        #                     "1": 0.5,
        #                     "2": 0.5,
        #                 }
        #             },
        #             "6": {
        #                 "b": {
        #                     "3": 0.5,
        #                     "4": 0.5,
        #                 }
        #             },
        #             "7": {
        #                 "b": {
        #                     "8": 0.5,
        #                     "9": 0.5,
        #                 }
        #             },
        #         }
        #     }
        # }

    def test_named_layers(self):
        inflated_ensemble = inflate_ensemble_with_lists_and_named_layers(
            self.example_ensemble_data
        )
        expected_result = [
            ["5", "5", "6", "6", "7", "7"],
            ["1", "2", "3", "4", "8", "9"],
            [0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
        ]
        self.assertEqual(expected_result, inflated_ensemble)

    # def test_named_layers(self):
    #     # Test with named layers
    #     df = convert_ensemble_to_dataframe(
    #         self.example_ensemble_data, contains_named_layers=True
    #     )
    #     self.assertEqual(list(df.columns), ["a", "b", "probability"])
    #     self.assertEqual(len(df), 6)

    # def test_unnamed_layers_with_columnames(self):
    #     # Test without named layers but with provided column names
    #     columnames = ["custom1", "a", "custom2", "b"]
    #     df = convert_ensemble_to_dataframe(
    #         self.example_ensemble_data,
    #         contains_named_layers=False,
    #         columnames=columnames,
    #     )
    #     self.assertEqual(list(df.columns), columnames + ["probability"])
    #     self.assertEqual(len(df), 6)

    # def test_unnamed_layers_without_columnames(self):
    #     # Test without named layers and without column names (should raise ValueError)
    #     with self.assertRaises(ValueError):
    #         convert_ensemble_to_dataframe(
    #             self.example_ensemble_data, contains_named_layers=False
    #         )

    # def test_named_layers_2(self):
    #     # Test with named layers
    #     df = convert_ensemble_to_dataframe(
    #         self.example_ensemble_data_2, contains_named_layers=True
    #     )
    #     self.assertEqual(list(df.columns), ["a", "b", "probability"])
    #     print(df)
    #     self.assertEqual(len(df), 6)

    # def test_unnamed_layers_with_columnames_2(self):
    #     # Test without named layers but with provided column names
    #     columnames = ["ensemble", "custom1", "a", "custom2", "b"]
    #     df = convert_ensemble_to_dataframe(
    #         self.example_ensemble_data_2, contains_named_layers=False, columnames=columnames
    #     )
    #     self.assertEqual(list(df.columns), columnames+["probability"])
    #     self.assertEqual(len(df), 6)

    # def test_unnamed_layers_without_columnames_2(self):
    #     # Test without named layers and without column names (should raise ValueError)
    #     with self.assertRaises(ValueError):
    #         convert_ensemble_to_dataframe(self.example_ensemble_data_2, contains_named_layers=False)


class test_inflate_ensemble_with_lists_without_named_layers(unittest.TestCase):
    def setUp(self):
        self.example_ensemble_data = {
            5: {1: 0.5, 2: 0.5},
            6: {3: 0.5, 4: 0.5},
        }

        self.example_ensemble_data_2 = {
            1: {5: {1: 0.5, 2: 0.5}, 6: {3: 0.5, 4: 0.5}},
            2: {5: {1: 0.5, 2: 0.5}, 6: {3: 0.5, 4: 0.5}},
        }

    def test_inflate_ensemble_with_lists_without_named_layers(self):
        inflated_ensemble = inflate_ensemble_with_lists_without_named_layers(
            self.example_ensemble_data
        )
        expected_result = [[5, 5, 6, 6], [1.0, 2.0, 3.0, 4.0], [0.5, 0.5, 0.5, 0.5]]
        self.assertEqual(expected_result, inflated_ensemble)

    def test_inflate_ensemble_with_lists_without_named_layers_2(self):
        inflated_ensemble = inflate_ensemble_with_lists_without_named_layers(
            self.example_ensemble_data_2
        )
        expected_result = [
            [1, 1, 1, 1, 2, 2, 2, 2],
            [5, 5, 6, 6, 5, 5, 6, 6],
            [1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0],
            [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
        ]
        self.assertEqual(expected_result, inflated_ensemble)

    # def test_unnamed_layers_with_columnames(self):
    #     # Test without named layers but with provided column names
    #     columnames = ["a", "b"]
    #     df = convert_ensemble_to_dataframe(
    #         self.example_ensemble_data,
    #         contains_named_layers=False,
    #         columnames=columnames,
    #     )
    #     self.assertEqual(list(df.columns), columnames + ["probability"])
    #     self.assertEqual(len(df), 6)

    # def test_named_layers_without_columnames(self):
    #     # Test named layers
    #     with self.assertRaises(ValueError):
    #         convert_ensemble_to_dataframe(
    #             self.example_ensemble_data, contains_named_layers=False
    #         )

    # def test_unnamed_layers_with_columnames_2(self):
    #     # Test without named layers but with provided column names
    #     columnames = ["a", "b", "c"]
    #     df = convert_ensemble_to_dataframe(
    #         self.example_ensemble_data_2,
    #         contains_named_layers=False,
    #         columnames=columnames,
    #     )
    #     self.assertEqual(list(df.columns), columnames + ["probability"])
    #     self.assertEqual(len(df), 18)

    # def test_unnamed_layers_without_columnames_2(self):
    #     # Test without named layers and without column names (should raise ValueError)
    #     with self.assertRaises(ValueError):
    #         convert_ensemble_to_dataframe(
    #             self.example_ensemble_data_2, contains_named_layers=False
    #         )


# class test__get_ensemble_structure(unittest.TestCase):
#     def test_basic_structure(self):
#         ensemble = {"a": {"b": {}, "c": {}}, "d": {"e": {}}}
#         expected_structure = {0: ["a", "d"], 1: ["b", "c", "e"]}
#         self.assertEqual(
#             _get_ensemble_structure(ensemble, {0: [], 1: []}, 2), expected_structure
#         )

#     def test_max_depth_reached(self):
#         ensemble = {"a": {"b": {}}, "c": {"d": {}}}
#         expected_structure = {0: ["a", "c"]}
#         self.assertEqual(
#             _get_ensemble_structure(ensemble, {0: []}, 1), expected_structure
#         )

#     def test_empty_ensemble(self):
#         ensemble = {}
#         expected_structure = {0: []}
#         self.assertEqual(
#             _get_ensemble_structure(ensemble, {0: []}, 1), expected_structure
#         )

#     def test_nested_empty_ensemble(self):
#         ensemble = {"a": {"b": {"c": {}}}}
#         expected_structure = {0: ["a"], 1: ["b"], 2: ["c"]}
#         self.assertEqual(
#             _get_ensemble_structure(ensemble, {0: [], 1: [], 2: []}, 3),
#             expected_structure,
#         )


# class test_get_ensemble_binsizes(unittest.TestCase):
#     def test_get_ensemble_binsizes_predetermined_binsizes(self):
#         binsizes = get_ensemble_binsizes(
#             config={},
#             ensemble={"0.1": 2, "0.2": 2},
#             data_layer_dict_entry={"binsizes": [1, 2]},
#         )

#         self.assertTrue(binsizes == [1, 2])

#     def test_get_ensemble_binsizes_no_scaling(self):

#         convolution_config = copy.copy(default_convolution_config)

#         binsizes = get_ensemble_binsizes(
#             config=convolution_config,
#             ensemble={"0.1": 2, "0.2": 2, "0.4": 3},
#             data_layer_dict_entry={},
#         )

#         #
#         np.testing.assert_array_almost_equal(binsizes, np.array([0.1, 0.15, 0.2]))

#     def test_get_ensemble_binsizes_factor_scaling(self):

#         convolution_config = copy.copy(default_convolution_config)

#         binsizes = get_ensemble_binsizes(
#             config=convolution_config,
#             ensemble={"0.1": 2, "0.2": 2, "0.4": 3},
#             data_layer_dict_entry={"conversion_factor": 2},
#         )

#         #
#         np.testing.assert_array_almost_equal(binsizes, np.array([0.2, 0.3, 0.4]))

#     def test_get_ensemble_binsizes_factor_function(self):

#         convolution_config = copy.copy(default_convolution_config)

#         binsizes = get_ensemble_binsizes(
#             config=convolution_config,
#             ensemble={"0.1": 2, "0.2": 2, "0.4": 3},
#             data_layer_dict_entry={"conversion_function": lambda x: 10**x},
#         )

#         #
#         np.testing.assert_array_almost_equal(
#             binsizes, np.array([0.29051909, 0.58272477, 1.16701535])
#         )


# class test_ensemble_marginalise_layer(unittest.TestCase):
#     def test_ensemble_marginalise_layer(self):
#         dummy_ensemble = {
#             "metallicity": {
#                 "0": {
#                     "delay_time": {
#                         "0": {"a": {"1": 1}, "b": {"1": 1}},
#                         "1": {"a": {"1": 2}, "b": {"1": 2}},
#                         "2": {"a": {"1": 3}, "b": {"1": 3}},
#                         "3": {"a": {"1": 4}, "b": {"1": 4}},
#                     }
#                 }
#             }
#         }

#         target_ensemble = {
#             "metallicity": {
#                 "0": {
#                     "delay_time": {
#                         "0": {"1": 2},
#                         "1": {"1": 4},
#                         "2": {"1": 6},
#                         "3": {"1": 8},
#                     }
#                 }
#             }
#         }

#         #
#         marginalized_ensemble = ensemble_marginalise_layer(
#             ensemble=dummy_ensemble, marginalisation_depth=4
#         )
#         self.assertTrue(target_ensemble == marginalized_ensemble)


# class test_get_deepest_data_layer_depth(unittest.TestCase):
#     def setUp(self):
#         self.data_layer_dict = {
#             "key1": 1,
#             "key2": {"layer_depth": 2},
#             "key3": 3,
#             "key4": {"layer_depth": 4},
#         }

#     def test_get_deepest_data_layer_depth(self):
#         deepest_depth = get_deepest_data_layer_depth(self.data_layer_dict)
#         self.assertEqual(
#             deepest_depth, 4
#         )  # Deepest data layer depth in the provided dictionary

#     def test_get_deepest_data_layer_depth_empty_dict(self):
#         empty_dict = {}
#         with self.assertRaises(ValueError):
#             get_deepest_data_layer_depth(empty_dict)

#     def test_get_deepest_data_layer_depth_invalid_input(self):
#         invalid_data_layer_dict = {"key1": [1, 2, 3]}
#         with self.assertRaises(ValueError):
#             get_deepest_data_layer_depth(invalid_data_layer_dict)


# class test_get_data_layer_dict_values(unittest.TestCase):
#     def setUp(self):
#         self.data_layer_dict = {
#             "key1": 1,
#             "key2": {"layer_depth": 2},
#             "key3": 3,
#             "key4": {"layer_depth": 4},
#         }

#     def test_get_data_layer_dict_values(self):
#         data_layer_values = get_data_layer_dict_values(self.data_layer_dict)
#         expected_values = [1, 2, 3, 4]
#         self.assertListEqual(data_layer_values, expected_values)

#     def test_get_data_layer_dict_values_invalid_input(self):
#         invalid_data_layer_dict = {"key1": [1, 2, 3]}
#         with self.assertRaises(ValueError):
#             get_data_layer_dict_values(invalid_data_layer_dict)


# class test_multiply_ensemble(unittest.TestCase):
#     def setUp(self):
#         self.nested_dict = {
#             "a": {"b": {"c": 1, "d": 2}, "e": 3},
#             "f": {"g": {"h": {"i": {"j": {"k": {"l": 10}}}}}},
#         }

#     def test_multiply_ensemble(self):
#         factor = 2
#         multiplied_dict = {
#             "a": {"b": {"c": 2, "d": 4}, "e": 6},
#             "f": {"g": {"h": {"i": {"j": {"k": {"l": 20}}}}}},
#         }
#         multiply_ensemble(self.nested_dict, factor)
#         self.assertEqual(self.nested_dict, multiplied_dict)


# class test_invert_data_layer_dict(unittest.TestCase):
#     def setUp(self):
#         self.data_layer_dict = {
#             "key1": 1,
#             "key2": {"layer_depth": 2},
#             "key3": 3,
#             "key4": {"layer_depth": 4},
#         }

#     def test_invert_data_layer_dict(self):
#         inverted_dict = invert_data_layer_dict(self.data_layer_dict)
#         expected_inverted_dict = {1: "key1", 2: "key2", 3: "key3", 4: "key4"}
#         self.assertDictEqual(inverted_dict, expected_inverted_dict)

#     def test_invert_data_layer_dict_invalid_input(self):
#         invalid_data_layer_dict = {"key1": [1, 2, 3]}
#         with self.assertRaises(ValueError):
#             invert_data_layer_dict(invalid_data_layer_dict)


# class test_get_depth_ensemble_first_endpoint(unittest.TestCase):
#     def setUp(self):
#         self.nested_dict = {
#             "a": {"b": {"c": 1, "d": 2}, "e": 3},
#             "f": {"g": {"h": {"i": {"j": {"k": {"l": 10}}}}}},
#         }

#     def test_get_depth_ensemble_first_endpoint(self):
#         depth = get_depth_ensemble_first_endpoint(self.nested_dict)
#         self.assertEqual(
#             depth, 3
#         )  # Depth of the first endpoint in the nested dictionary


# class test_get_max_depth_ensemble(unittest.TestCase):
#     def setUp(self):
#         self.nested_dict = {
#             "a": {"b": {"c": 1, "d": 2}, "e": 3},
#             "f": {"g": {"h": {"i": {"j": {"k": {"l": 10}}}}}},
#         }

#     def test_find_max_depth(self):
#         max_depth = get_max_depth_ensemble(self.nested_dict)
#         self.assertEqual(max_depth, 7)  # Maximum depth is 7 in this nested dictionary


# class test_get_depth_ensemble_all_endpoints(unittest.TestCase):
#     def setUp(self):
#         self.nested_dict = {
#             "a": {"b": {"c": 1, "d": 2}, "e": 3},
#             "f": {"g": {"h": {"i": {"j": {"k": {"l": 10}}}}}},
#         }

#     def test_get_depth_ensemble_all_endpoints(self):
#         endpoint_depths = get_depth_ensemble_all_endpoints(self.nested_dict)
#         expected_depths = [
#             3,
#             3,
#             2,
#             7,
#         ]  # Depths of all endpoints in the nested dictionary
#         self.assertEqual(endpoint_depths, expected_depths)


# class test_shift_layers_dict(unittest.TestCase):
#     def setUp(self):
#         self.data_layer_dict = {
#             "layer1": 1,
#             "layer2": {"layer_depth": 2},
#             "layer3": {"layer_depth": 3},
#         }

#     def test_shift_layers_dict_positive_shift(self):
#         shift_value = 5
#         expected_result = {
#             "layer1": 6,
#             "layer2": {"layer_depth": 7},
#             "layer3": {"layer_depth": 8},
#         }
#         self.assertEqual(
#             shift_layers_dict(self.data_layer_dict, shift_value), expected_result
#         )

#     def test_shift_layers_dict_negative_shift(self):
#         shift_value = -2
#         expected_result = {
#             "layer1": -1,
#             "layer2": {"layer_depth": 0},
#             "layer3": {"layer_depth": 1},
#         }
#         self.assertEqual(
#             shift_layers_dict(self.data_layer_dict, shift_value), expected_result
#         )

#     def test_shift_layers_dict_zero_shift(self):
#         shift_value = 0
#         expected_result = {
#             "layer1": 1,
#             "layer2": {"layer_depth": 2},
#             "layer3": {"layer_depth": 3},
#         }
#         self.assertEqual(
#             shift_layers_dict(self.data_layer_dict, shift_value), expected_result
#         )

#     def test_shift_layers_dict_empty_dict(self):
#         empty_dict = {}
#         shift_value = 5
#         expected_result = {}
#         self.assertEqual(shift_layers_dict(empty_dict, shift_value), expected_result)

#     def test_shift_layers_dict_unsupported(self):
#         shift_value = 1
#         data_layer_dict = {
#             "layer1": 1,
#             "layer2": {"layer_depth": 2},
#             "layer3": "unsupported",
#         }
#         with self.assertRaises(ValueError):
#             shift_layers_dict(data_layer_dict, shift_value)


# class test_shift_data_layer(unittest.TestCase):
#     def test_shift_data_layer_int_value(self):
#         data_layer_dict = {"layer1": 1}
#         key = "layer1"
#         shift = 5
#         expected_result = {"layer1": 6}
#         shift_data_layer(data_layer_dict, key, shift)
#         self.assertEqual(data_layer_dict, expected_result)

#     def test_shift_data_layer_dict_value(self):
#         data_layer_dict = {"layer2": {"layer_depth": 2}}
#         key = "layer2"
#         shift = -2
#         expected_result = {"layer2": {"layer_depth": 0}}
#         shift_data_layer(data_layer_dict, key, shift)
#         self.assertEqual(data_layer_dict, expected_result)

#     def test_shift_data_layer_unsupported_input(self):
#         data_layer_dict = {"layer3": "unsupported"}
#         key = "layer3"
#         shift = 5
#         with self.assertRaises(ValueError):
#             shift_data_layer(data_layer_dict, key, shift)


# class test_shift_layers_list(unittest.TestCase):
#     def test_shift_layers_list_positive_shift(self):
#         layer_list = [1, 2, 3, 4]
#         shift_value = 5
#         expected_result = [6, 7, 8, 9]
#         self.assertEqual(shift_layers_list(layer_list, shift_value), expected_result)

#     def test_shift_layers_list_negative_shift(self):
#         layer_list = [1, 2, 3, 4]
#         shift_value = -2
#         expected_result = [-1, 0, 1, 2]
#         self.assertEqual(shift_layers_list(layer_list, shift_value), expected_result)

#     def test_shift_layers_list_zero_shift(self):
#         layer_list = [1, 2, 3, 4]
#         shift_value = 0
#         expected_result = [1, 2, 3, 4]
#         self.assertEqual(shift_layers_list(layer_list, shift_value), expected_result)

#     def test_shift_layers_list_empty_list(self):
#         layer_list = []
#         shift_value = 5
#         expected_result = []
#         self.assertEqual(shift_layers_list(layer_list, shift_value), expected_result)


# class test_strip_ensemble_endpoints(unittest.TestCase):
#     # TODO: add with units
#     def setUp(self):
#         self.nested_dict = {"a": {"b": {"c": 1, "d": 2}, "e": 3}, "f": {"g": 4}, "h": 5}

#     def test_strip_ensemble_endpoints(self):
#         expected_endpoints = [1, 2, 3, 4, 5]
#         ensemble, endpoints, found_units = strip_ensemble_endpoints(self.nested_dict)
#         self.assertEqual(endpoints, expected_endpoints)

#         # Ensure original ensemble endpoints are set to 0
#         extracted_endpoints, found_units = extract_endpoints(ensemble)
#         self.assertTrue(all(endpoint == 0 for endpoint in extracted_endpoints))


# class test_check_if_value_layer(unittest.TestCase):
#     def setUp(self):
#         self.name_keys = {"name1": 1, "name2": 2, "name3": 3}
#         self.value_keys = {"1": 1, "2": 2, "3": 3}
#         self.mixed_keys = {"1": 1, "name2": 2, "3": 3}

#     def test_value_layer(self):
#         self.assertTrue(check_if_value_layer(self.value_keys.keys()))

#     def test_name_layer(self):
#         self.assertFalse(check_if_value_layer(self.name_keys.keys()))

#     def test_mixed_layer(self):
#         self.assertFalse(check_if_value_layer(self.mixed_keys.keys()))


# class test_get_layer_iterable(unittest.TestCase):
#     def setUp(self):
#         self.name_keys = {"name1": 1, "name2": 2, "name3": 3}
#         self.value_keys = {"1": 1, "2": 2, "3": 3}
#         self.unsorted_value_keys = {"1": 1, "-2": 2, "3": 3}
#         self.mixed_keys = {"1": 1, "name2": 2, "3": 3}

#     def test_get_layer_iterable_value_layer(self):
#         self.assertEqual(
#             list(get_layer_iterable(self.value_keys, True)), ["1", "2", "3"]
#         )

#     def test_get_layer_iterable_unsorted_value_layer(self):
#         self.assertEqual(
#             list(get_layer_iterable(self.unsorted_value_keys, True)), ["-2", "1", "3"]
#         )

#     def test_get_layer_iterable_name_layer(self):
#         self.assertEqual(
#             list(get_layer_iterable(self.name_keys, False)), ["name1", "name2", "name3"]
#         )

#     def test_get_layer_iterable_mixed_layer(self):
#         self.assertEqual(
#             list(get_layer_iterable(self.mixed_keys, False)), ["1", "name2", "3"]
#         )


# class test_check_if_value_layer_and_get_layer_iterable(unittest.TestCase):
#     def setUp(self):
#         self.name_keys = {"name1": 1, "name2": 2, "name3": 3}
#         self.value_keys = {"1": 1, "2": 2, "3": 3}
#         self.unsorted_value_keys = {"1": 1, "-2": 2, "3": 3}
#         self.mixed_keys = {"1": 1, "name2": 2, "3": 3}

#     def test_value_layer_and_get_layer_iterable(self):
#         is_value_layer, iterable = check_if_value_layer_and_get_layer_iterable(
#             self.value_keys
#         )
#         self.assertTrue(is_value_layer)
#         self.assertEqual(list(iterable), ["1", "2", "3"])

#     def test_unsorted_value_layer_and_get_layer_iterable(self):
#         is_value_layer, iterable = check_if_value_layer_and_get_layer_iterable(
#             self.unsorted_value_keys
#         )
#         self.assertTrue(is_value_layer)
#         self.assertEqual(list(iterable), ["-2", "1", "3"])

#     def test_name_layer_and_get_layer_iterable(self):
#         is_value_layer, iterable = check_if_value_layer_and_get_layer_iterable(
#             self.name_keys
#         )
#         self.assertFalse(is_value_layer)
#         self.assertEqual(list(iterable), ["name1", "name2", "name3"])

#     def test_mixed_layer_and_get_layer_iterable(self):
#         is_value_layer, iterable = check_if_value_layer_and_get_layer_iterable(
#             self.mixed_keys
#         )
#         self.assertFalse(is_value_layer)
#         self.assertEqual(list(iterable), ["1", "name2", "3"])


# class test_extract_endpoints(unittest.TestCase):
#     def setUp(self):
#         self.nested_dict = {"a": {"b": {"c": 1, "d": 2}, "e": 3}, "f": {"g": 4}, "h": 5}
#         self.flat_dict = {"a": 1, "b": 2, "c": 3}
#         self.empty_dict = {}

#     def test_extract_endpoints_nested_dict(self):
#         expected_endpoints = [1, 2, 3, 4, 5]
#         endpoints, found_units = extract_endpoints(self.nested_dict)
#         self.assertEqual(endpoints, expected_endpoints)

#     def test_extract_endpoints_flat_dict(self):
#         expected_endpoints = [1, 2, 3]
#         endpoints, found_units = extract_endpoints(self.flat_dict)
#         self.assertEqual(endpoints, expected_endpoints)

#     def test_extract_endpoints_empty_dict(self):
#         endpoints, found_units = extract_endpoints(self.empty_dict)
#         self.assertEqual(endpoints, [])


# class test_attach_endpoints(unittest.TestCase):
#     def setUp(self):
#         self.nested_dict = {
#             "a": {"b": {"c": None, "d": None}, "e": None},
#             "f": {"g": None},
#             "h": None,
#         }
#         self.endpoint_array = [1, 2, 3, 4, 5]

#     def test_attach_endpoints_nested_dict(self):
#         expected_dict = {"a": {"b": {"c": 1, "d": 2}, "e": 3}, "f": {"g": 4}, "h": 5}
#         attach_endpoints(self.nested_dict, self.endpoint_array)
#         self.assertEqual(self.nested_dict, expected_dict)

#     def test_attach_endpoints_empty_dict(self):
#         empty_dict = {}
#         with self.assertRaises(ValueError):
#             attach_endpoints(empty_dict, self.endpoint_array)


# class test_set_endpoints(unittest.TestCase):
#     def setUp(self):
#         self.nested_dict = {
#             "a": {"b": {"c": None, "d": None}, "e": None},
#             "f": {"g": None},
#             "h": None,
#         }
#         self.value = 100

#     def test_set_endpoints_nested_dict(self):
#         expected_dict = {
#             "a": {"b": {"c": 100, "d": 100}, "e": 100},
#             "f": {"g": 100},
#             "h": 100,
#         }
#         set_endpoints(self.nested_dict, self.value)
#         self.assertEqual(self.nested_dict, expected_dict)

#     def test_set_endpoints_empty_dict(self):
#         empty_dict = {}
#         with self.assertRaises(ValueError):
#             set_endpoints(empty_dict, self.value)


# class test_get_ensemble_structure(unittest.TestCase):
#     def setUp(self):
#         # Define some sample ensembles for testing
#         self.ensemble_single_depth = {"layer1": {"node1": 1, "node2": 2, "node3": 3}}

#         self.ensemble_multiple_depth = {
#             "layer1": {
#                 "node1": {"subnode1": 11, "subnode2": 12},
#                 "node2": {"subnode1": 21, "subnode2": 22},
#             }
#         }

#         self.ensemble_named_layer = {
#             "layer1": {
#                 "node1": {"subnode1": 11, "subnode2": 12},
#                 "node2": {"subnode1": 21, "subnode2": 22},
#             },
#             "layer2": {
#                 "node1": {"subnode1": 31, "subnode2": 32},
#                 "node2": {"subnode1": 41, "subnode2": 42},
#             },
#         }

#     def test_get_ensemble_structure_single_depth(self):
#         expected_structure = {0: ["layer1"], 1: ["node1", "node2", "node3"]}
#         self.assertEqual(
#             get_ensemble_structure(self.ensemble_single_depth), expected_structure
#         )

#     def test_get_ensemble_structure_multiple_depth(self):
#         expected_structure = {
#             0: ["layer1"],
#             1: ["node1", "node2"],
#             2: ["subnode1", "subnode2"],
#         }
#         self.assertEqual(
#             get_ensemble_structure(self.ensemble_multiple_depth), expected_structure
#         )

#     def test_get_ensemble_structure_named_layer(self):
#         # Test with named layer list provided

#         # Test with named layer containing multiple values
#         with self.assertRaises(ValueError):
#             get_ensemble_structure(self.ensemble_named_layer, named_layer_list=[1])


if __name__ == "__main__":
    unittest.main()
