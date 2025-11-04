"""
This is the unittest file for the convolve_binned_data.py source file
"""

import unittest

import numpy as np

from syntheticstellarpopconvolve.convolve_binned_data import calculate_overlap_fractions


class test_calculate_overlap_fractions(unittest.TestCase):
    """
    TODO: make SFR bins of larger width to get the normalized value right
    """

    def test_calculate_overlap_fractions_one_bin_overlap(self):

        #
        delay_time_data_bin_info = {
            "delay_time_data_bin_edges": np.arange(0, 2, 1),
        }

        #
        shift = 5.6

        sfr_bin_edges = np.arange(0, 20, 5)
        sfr_bin_sizes = np.diff(sfr_bin_edges)

        delay_time_data_bin_edges = delay_time_data_bin_info[
            "delay_time_data_bin_edges"
        ]

        #
        left_delay_time_data_bin_edges = delay_time_data_bin_edges[:-1]
        right_delay_time_data_bin_edges = delay_time_data_bin_edges[1:]

        shifted_left_delay_time_data_bin_edges = left_delay_time_data_bin_edges + shift
        shifted_right_delay_time_data_bin_edges = (
            right_delay_time_data_bin_edges + shift
        )

        #
        overlap_fractions = calculate_overlap_fractions(
            shifted_left_delay_time_data_bin_edge=shifted_left_delay_time_data_bin_edges[
                0
            ],
            shifted_right_delay_time_data_bin_edge=shifted_right_delay_time_data_bin_edges[
                0
            ],
            sfr_bin_sizes=sfr_bin_sizes,
            sfr_bin_edges=sfr_bin_edges,
        )

        expected_overlap_fractions = {
            "combined_overlap_array": np.array([0.0, 1.0, 0.0]),
            "normalized_combined_overlap_array": np.array([0.0, 0.2, 0.0]),
            "time_bin_fraction": np.array([0.0, 1.0, 0.0]),
            "cumulative_time_bin_fraction": np.array([0.0, 1.0, 1.0]),
            "change_cumulative_time_bin_fraction": np.array([1.0, 0.0]),
            "non_zero_overlap_with_sfr_bins": np.array([1]),
        }

        for key in expected_overlap_fractions.keys():
            self.assertTrue(
                np.array_equal(overlap_fractions[key], expected_overlap_fractions[key])
            )

    def test_calculate_overlap_fractions_multi_bin_overlap(self):

        #
        delay_time_data_bin_info = {
            "delay_time_data_bin_edges": np.arange(0, 5, 4),
        }

        #
        shift = 0.25

        sfr_bin_edges = np.arange(0, 6, 1)
        sfr_bin_sizes = np.diff(sfr_bin_edges)

        delay_time_data_bin_edges = delay_time_data_bin_info[
            "delay_time_data_bin_edges"
        ]

        #
        left_delay_time_data_bin_edges = delay_time_data_bin_edges[:-1]
        right_delay_time_data_bin_edges = delay_time_data_bin_edges[1:]

        shifted_left_delay_time_data_bin_edges = left_delay_time_data_bin_edges + shift
        shifted_right_delay_time_data_bin_edges = (
            right_delay_time_data_bin_edges + shift
        )

        #
        overlap_fractions = calculate_overlap_fractions(
            shifted_left_delay_time_data_bin_edge=shifted_left_delay_time_data_bin_edges[
                0
            ],
            shifted_right_delay_time_data_bin_edge=shifted_right_delay_time_data_bin_edges[
                0
            ],
            sfr_bin_sizes=sfr_bin_sizes,
            sfr_bin_edges=sfr_bin_edges,
        )

        expected_overlap_fractions = {
            "combined_overlap_array": np.array([0.75, 1.0, 1.0, 1.0, 0.25]),
            "normalized_combined_overlap_array": np.array([0.75, 1.0, 1.0, 1.0, 0.25]),
            "time_bin_fraction": np.array([0.1875, 0.25, 0.25, 0.25, 0.0625]),
            "cumulative_time_bin_fraction": np.array(
                [0.1875, 0.4375, 0.6875, 0.9375, 1.0]
            ),
            "change_cumulative_time_bin_fraction": np.array([0.25, 0.25, 0.25, 0.0625]),
            "non_zero_overlap_with_sfr_bins": np.array([0, 1, 2, 3, 4]),
        }

        for key in expected_overlap_fractions.keys():
            self.assertTrue(
                np.array_equal(overlap_fractions[key], expected_overlap_fractions[key])
            )


if __name__ == "__main__":
    unittest.main()
