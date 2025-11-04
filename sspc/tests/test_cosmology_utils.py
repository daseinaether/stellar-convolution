"""
Testcases for convolve_ensembles file
"""

import unittest

from astropy.cosmology import Planck13 as cosmo  # Planck 2013

from syntheticstellarpopconvolve.cosmology_utils import (
    age_of_universe_to_redshift,
    lookback_time_to_redshift,
    redshift_to_age_of_universe,
    redshift_to_lookback_time,
)
from syntheticstellarpopconvolve.general_functions import temp_dir

TMP_DIR = temp_dir(
    "tests", "tests_convolution", "tests_cosmology_utils", clean_path=True
)


class test_age_of_universe_to_redshift(unittest.TestCase):
    """ """

    def test_age_of_universe_to_redshift(self):

        redshift = age_of_universe_to_redshift(age_of_universe=1, cosmology=cosmo)

        self.assertAlmostEqual(
            5.676612273980656, redshift, 6, "redshift at age of universe=1Gyr incorrect"
        )


class test_lookback_time_to_redshift(unittest.TestCase):
    """ """

    def test_lookback_time_to_redshift(self):
        redshift = lookback_time_to_redshift(lookback_time=1, cosmology=cosmo)

        self.assertAlmostEqual(
            0.07304795128154057, redshift, 6, "redshift at lookback time=1Gyr incorrect"
        )

    def test_lookback_time_to_redshift_at_zero(self):

        redshift = lookback_time_to_redshift(lookback_time=0, cosmology=cosmo)

        self.assertAlmostEqual(
            6.930802008253847e-08,
            redshift,
            8,
            "redshift at lookback time=1Gyr incorrect",
        )


class test_redshift_to_lookback_time(unittest.TestCase):
    """ """

    def test_redshift_to_lookback_time(self):
        lookback_time = redshift_to_lookback_time(redshift=1, cosmology=cosmo)

        self.assertAlmostEqual(
            7.9331257328432265,
            lookback_time.value,
            6,
            "lookback time at redshift=1 incorrect",
        )


class test_redshift_to_age_of_universe(unittest.TestCase):
    """ """

    def test_redshift_to_age_of_universe(self):
        age_of_universe = redshift_to_age_of_universe(redshift=1, cosmology=cosmo)

        self.assertAlmostEqual(
            5.863165023612126,
            age_of_universe.value,
            6,
            "age of universe at redshift=1 incorrect",
        )


# #
# ##################
# # Functions to test the interpolation
# def test_redshift_to_lookback_time(config, size_test_sample, test_log=False):
#     """
#     Redshift to lookback time interpolation test
#     """

#     # Load the interpolators:
#     interpolators_dict = load_interpolation_data(config)

#     # Set up the interpolator
#     redshift_to_lookback_time_interpolator = interpolators_dict[
#         "redshift_to_lookback_time_interpolator"
#     ]

#     # Create test sample redshift
#     if test_log:
#         random_redshift_sample = 10 ** np.random.uniform(
#             np.log10(
#                 config["min_interpolation_redshift"]
#                 + config["min_redshift_change_if_log_sampling"]
#             ),
#             np.log10(config["max_interpolation_redshift"]),
#             size_test_sample,
#         )
#     else:
#         random_redshift_sample = np.random.uniform(
#             config["min_interpolation_redshift"],
#             config["max_interpolation_redshift"],
#             size_test_sample,
#         )
#     # print(random_redshift_sample)

#     # Transform the redshift to time to compare to
#     start_real_time = time.time()
#     true_lookback_time_results = (
#         config["cosmo"].age(0).value - config["cosmo"].age(random_redshift_sample).value
#     )
#     # print(true_lookback_time_results)
#     stop_real_time = time.time()
#     print(
#         "Took {} seconds to generate the real data for {} samples".format(
#             stop_real_time - start_real_time, size_test_sample
#         )
#     )

#     # Create time values by interpolation
#     start_interpolation_time = time.time()
#     interpolated_lookback_time_results = redshift_to_lookback_time_interpolator(
#         random_redshift_sample
#     )
#     # print(redshift_to_time_interpolated_time_sample)
#     stop_interpolation_time = time.time()
#     print(
#         "Took {} seconds to interpolate the data for {} samples".format(
#             stop_interpolation_time - start_interpolation_time, size_test_sample
#         )
#     )
#     print(
#         "Generating redshift to lookback time via interpolation is {:.2e} times faster than normal function".format(
#             (stop_real_time - start_real_time)
#             / (stop_interpolation_time - start_interpolation_time)
#         )
#     )

#     # Calculat the abs fractional errors
#     fractional_error = np.abs(
#         (true_lookback_time_results - interpolated_lookback_time_results)
#         / true_lookback_time_results
#     )
#     # print(fractional_error)

#     result_dict = {
#         "redshift_sample": random_redshift_sample,
#         "true_lookback_time_results": true_lookback_time_results,
#         "interpolated_lookback_time_results": interpolated_lookback_time_results,
#         "fractional_error": fractional_error,
#         "sample_size": size_test_sample,
#     }

#     return result_dict

# def test_lookback_time_to_redshift(config, size_test_sample, test_log=False):
#     """
#     Lookback time to redshift test
#     """

#     #
#     min_lookback_time = redshift_to_lookback_time(config["min_interpolation_redshift"])
#     max_lookback_time = redshift_to_lookback_time(config["max_interpolation_redshift"])

#     # Load the interpolators:
#     interpolators_dict = load_interpolation_data(config)

#     # Set up the interpolator
#     lookback_time_to_redshift_interpolator = interpolators_dict[
#         "lookback_time_to_redshift_interpolator"
#     ]

#     # Create test sample of lookback times
#     if test_log:
#         if min_lookback_time == 0:
#             min_lookback_time += redshift_to_lookback_time(
#                 config["min_redshift_change_if_log_sampling"]
#             )
#         random_lookback_times_sample = 10 ** np.random.uniform(
#             np.log10(min_lookback_time.value),
#             np.log10(max_lookback_time.value),
#             size_test_sample,
#         )
#     else:
#         random_lookback_times_sample = np.random.uniform(
#             min_lookback_time.value, max_lookback_time.value, size_test_sample
#         )
#     # print(random_lookback_times_sample)

#     # Transform the redshift to time to compare to
#     start_real_time = time.time()
#     true_redshift_results = np.array(
#         [
#             lookback_time_to_redshift(lookback_time)
#             for lookback_time in random_lookback_times_sample
#         ]
#     )
#     # print(true_redshift_results)
#     stop_real_time = time.time()
#     print(
#         "Took {} seconds to generate the real data for {} samples".format(
#             stop_real_time - start_real_time, size_test_sample
#         )
#     )

#     # Create redshift values by interpolation
#     start_interpolation_time = time.time()
#     interpolated_redshift_results = lookback_time_to_redshift_interpolator(
#         random_lookback_times_sample
#     )
#     # print(interpolated_redshift_results)
#     stop_interpolation_time = time.time()
#     print(
#         "Took {} seconds to interpolate the data for {} samples".format(
#             stop_interpolation_time - start_interpolation_time, size_test_sample
#         )
#     )
#     print(
#         "Generating lookback time to redshift via interpolation is {:.2e} times faster than normal function".format(
#             (stop_real_time - start_real_time)
#             / (stop_interpolation_time - start_interpolation_time)
#         )
#     )

#     # Calculate the abs fractional errors
#     fractional_error = np.abs(
#         (true_redshift_results - interpolated_redshift_results) / true_redshift_results
#     )
#     # print(fractional_error)

#     result_dict = {
#         "lookback_times_sample": random_lookback_times_sample,
#         "true_redshift_results": true_redshift_results,
#         "interpolated_redshift_results": interpolated_redshift_results,
#         "fractional_error": fractional_error,
#         "sample_size": size_test_sample,
#     }

#     return result_dict


if __name__ == "__main__":
    unittest.main()
