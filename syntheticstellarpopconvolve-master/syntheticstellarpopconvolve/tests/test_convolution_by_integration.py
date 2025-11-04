"""
This is the unittest file for the convolution_by_integration.py source file
"""

import unittest

import numpy as np

from syntheticstellarpopconvolve import (
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.convolution_by_integration import (
    convolution_by_integration_post_convolution_hook_wrapper,
)
from syntheticstellarpopconvolve.convolution_by_sampling import (
    select_dict_entries_with_new_indices,
)
from syntheticstellarpopconvolve.general_functions import create_time_bin_info_dict


def post_convolution_function_no_change(convolution_results):
    """ """

    return convolution_results


def post_convolution_function_sliced(convolution_results):
    """ """

    return select_dict_entries_with_new_indices(convolution_results, np.array([1]))


class test_convolution_by_integration_post_convolution_hook_wrapper(unittest.TestCase):
    def test_convolution_by_integration_post_convolution_hook_wrapper(self):

        convolution_results = {"yield": np.array([1, 2])}

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
        convolution_results = convolution_by_integration_post_convolution_hook_wrapper(
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

    def test_convolution_by_interation_post_convolution_hook_wrapper_changed_length(
        self,
    ):

        convolution_results = {"yield": np.array([1, 2])}

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

        #
        with self.assertRaises(ValueError):
            convolution_results = (
                convolution_by_integration_post_convolution_hook_wrapper(
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
            )


if __name__ == "__main__":
    unittest.main()
