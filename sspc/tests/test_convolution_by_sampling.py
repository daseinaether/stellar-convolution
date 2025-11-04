"""
This is the unittest file for the convolution_by_sampling.py source file
"""

import unittest

import astropy.units as u
import numpy as np

from syntheticstellarpopconvolve import (
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.convolution_by_sampling import (
    convolution_by_sampling_post_convolution_hook_wrapper,
    sample_systems,
    select_dict_entries_with_new_indices,
)
from syntheticstellarpopconvolve.general_functions import create_time_bin_info_dict

np.random.seed(0)


class test_select_dict_entries_with_new_indices(unittest.TestCase):
    def test_select_dict_entries_with_new_indices(self):

        sampled_data_dict = {
            "values": np.array([3, 2, 1]),
            "sampled_indices": np.array([1, 2, 3]),
        }

        new_indices = np.array([0, 1, 1, 2, 2, 2])

        new_data = select_dict_entries_with_new_indices(
            sampled_data_dict=sampled_data_dict, new_indices=new_indices
        )

        np.testing.assert_array_equal(new_data["values"], np.array([3, 2, 2, 1, 1, 1]))

        np.testing.assert_array_equal(
            new_data["sampled_indices"], np.array([1, 2, 2, 3, 3, 3])
        )


class test_sample_systems(unittest.TestCase):
    def test_sample_systems(self):

        yield_array = np.array([10, 5, 2, 1, 0.5])
        lookback_time_bin_size = 1 * u.Gyr
        lookback_time_bin_lower_edge = 1 * u.Gyr
        config = default_convolution_config
        convolution_instruction = default_convolution_instruction

        convolution_results = sample_systems(
            yield_array=yield_array,
            lookback_time_bin_size=lookback_time_bin_size,
            lookback_time_bin_lower_edge=lookback_time_bin_lower_edge,
            config=config,
            convolution_instruction=convolution_instruction,
        )

        sampled_indices = convolution_results["sampled_indices"]

        unique, counts = np.unique(sampled_indices, return_counts=True)

        np.testing.assert_array_equal(unique, np.array([0, 1, 2, 3, 4]))
        np.testing.assert_array_equal(counts, np.array([10, 5, 2, 1, 1]))


def post_convolution_function_no_change(convolution_results):
    """ """

    return convolution_results


def post_convolution_function_sliced(convolution_results):
    """ """

    return select_dict_entries_with_new_indices(convolution_results, np.array([1]))


class test_convolution_by_sampling_post_convolution_hook_wrapper(unittest.TestCase):
    def test_convolution_by_sampling_post_convolution_hook_wrapper(self):

        convolution_results = {"sampled_indices": np.array([1, 2])}

        convolution_instruction = {
            **default_convolution_instruction,
            "convolution_type": "sampling",
            "post_convolution_function": post_convolution_function_no_change,
        }

        time_bin_info_dict = create_time_bin_info_dict(
            config=default_convolution_config,
            convolution_instruction=convolution_instruction,
            bin_number=1,
            bin_edge_lower=0.5,
            bin_center=0.75,
            bin_size=1,
            bin_type="lookback",
        )

        # any change is allowed
        convolution_results = convolution_by_sampling_post_convolution_hook_wrapper(
            config={**default_convolution_config},
            convolution_instruction={
                **convolution_instruction,
                "post_convolution_function": post_convolution_function_no_change,
            },
            convolution_results=convolution_results,
            sfr_dict={"lookback_time_edges": np.array([1, 2])},
            data_dict={"values": np.array([5, 2])},
            time_bin_info_dict=time_bin_info_dict,
        )

        convolution_results = convolution_by_sampling_post_convolution_hook_wrapper(
            config={**default_convolution_config},
            convolution_instruction={
                **convolution_instruction,
                "post_convolution_function": post_convolution_function_sliced,
            },
            convolution_results=convolution_results,
            sfr_dict={"lookback_time_edges": np.array([1, 2])},
            data_dict={"values": np.array([5, 2])},
            time_bin_info_dict=time_bin_info_dict,
        )


if __name__ == "__main__":
    unittest.main()
