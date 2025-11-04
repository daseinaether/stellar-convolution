"""
Testcases for general_functions file
"""

import copy
import json
import logging
import os
import tempfile
import unittest

import astropy.units as u
import h5py
import numpy as np
import pandas as pd

from syntheticstellarpopconvolve import (
    default_convolution_config,
    default_convolution_instruction,
)
from syntheticstellarpopconvolve.check_and_prepare_output_file import (
    check_and_prepare_output_file,
)
from syntheticstellarpopconvolve.check_and_update_convolution_config import (
    check_and_update_convolution_config,
)
from syntheticstellarpopconvolve.general_functions import (  # calculate_digitized_sfr_rates,; calculate_origin_time_array,
    JsonCustomEncoder,
    calculate_bin_edges,
    calculate_bincenters,
    check_required,
    create_job_dict,
    create_time_bin_info_dict,
    extract_data,
    extract_unit_dict,
    generate_boilerplate_outputfile,
    generate_data_dict,
    generate_group_name,
    get_normalized_yield_unit,
    get_physical_dimensions,
    get_tmp_dir,
    get_username,
    handle_custom_scaling_or_conversion,
    has_unit,
    is_mass_unit,
    is_time_unit,
    pad_function,
    sample_around_bin_center,
    temp_dir,
)

TMP_DIR = temp_dir(
    "tests", "tests_convolution", "tests_general_functions", clean_path=True
)

np.random.seed(0)


class test_generate_data_dict(unittest.TestCase):
    def test_generate_data_dict_events(self):

        #
        output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5_sfr_only.h5")
        generate_boilerplate_outputfile(output_hdf5_filename)

        ##############
        # SET UP DATA
        self.dummy_data = {
            "delay_time": np.array([0, 1, 2, 3]),
            "probability": np.array([1, 2, 3, 4]),
        }
        dummy_df = pd.DataFrame.from_records(self.dummy_data)

        ##############
        # Store data in pandas
        dummy_df.to_hdf(output_hdf5_filename, key="input_data/{}".format("dummy"))

        #
        self.convolution_config = copy.copy(default_convolution_config)

        # Set up SFR
        self.convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]) * u.yr,
            "starformation_rate_array": np.array([1, 1, 1, 1, 1])
            * u.Msun
            / u.yr
            / u.Gpc**3,
        }

        # set up convolution bins
        self.convolution_config["convolution_lookback_time_bin_edges"] = (
            np.array([0, 1, 2, 3, 4]) * u.yr
        )

        # lookback time convolution only
        self.convolution_config["time_type"] = "lookback_time"

        #
        self.convolution_config["output_filename"] = output_hdf5_filename

        self.convolution_config["redshift_interpolator_data_output_filename"] = (
            os.path.join(TMP_DIR, "interpolator_dict.p")
        )

        #
        self.convolution_config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "dummy",
                "output_data_name": "dummy",
                "convolution_type": "integrate",
                "data_column_dict": {
                    "delay_time": "delay_time",
                    "normalized_yield": "probability",
                },
            },
        ]

        #
        self.convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")

        #
        check_and_prepare_output_file(config=self.convolution_config)

        #
        check_and_update_convolution_config(self.convolution_config)

        #
        normal_convolution_instructions = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
        }

        _, data_dict, _ = generate_data_dict(
            config=self.convolution_config,
            convolution_instruction=normal_convolution_instructions,
        )

        #
        np.testing.assert_array_equal(
            data_dict["delay_time"], self.dummy_data["delay_time"] * u.yr
        )


class test_extract_data(unittest.TestCase):
    """ """

    def setUp(self):
        #
        output_hdf5_filename = os.path.join(TMP_DIR, "output_hdf5_sfr_only.h5")
        generate_boilerplate_outputfile(output_hdf5_filename)

        ##############
        # SET UP DATA
        self.dummy_data = {
            "delay_time": np.array([0, 1, 2, 3]),
            "probability": np.array([1, 2, 3, 4]),
        }
        dummy_df = pd.DataFrame.from_records(self.dummy_data)

        ##############
        # Store data in pandas
        dummy_df.to_hdf(output_hdf5_filename, key="input_data/{}".format("dummy"))

        #
        self.convolution_config = copy.copy(default_convolution_config)

        # Set up SFR
        self.convolution_config["SFR_info"] = {
            "lookback_time_bin_edges": np.array([0, 1, 2, 3, 4, 5]),
            "starformation_rate_array": np.array([1, 1, 1, 1, 1])
            * u.Msun
            / u.yr
            / u.Gpc**3,
        }

        # set up convolution bins
        self.convolution_config["convolution_time_bin_edges"] = np.array(
            [0, 1, 2, 3, 4]
        )

        # lookback time convolution only
        self.convolution_config["time_type"] = "lookback_time"

        #
        self.convolution_config["output_filename"] = output_hdf5_filename

        self.convolution_config["redshift_interpolator_data_output_filename"] = (
            os.path.join(TMP_DIR, "interpolator_dict.p")
        )

        #
        self.convolution_config["convolution_instructions"] = [
            {
                **default_convolution_instruction,
                "input_data_name": "dummy",
                "output_data_name": "dummy",
                "data_column_dict": {
                    "delay_time": "delay_time",
                    "normalized_yield": "probability",
                },
            },
        ]

        #
        self.convolution_config["tmp_dir"] = os.path.join(TMP_DIR, "tmp")

        #
        check_and_prepare_output_file(config=self.convolution_config)

    def test_extract_data_normal(self):
        #
        normal_convolution_instructions = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "data_column_dict": {
                "delay_time": "delay_time",
                "normalized_yield": "probability",
            },
        }

        #
        _, data_dict, _ = extract_data(
            config=self.convolution_config,
            convolution_instruction=normal_convolution_instructions,
        )

        #
        np.testing.assert_array_equal(
            data_dict["delay_time"], self.dummy_data["delay_time"] * u.yr
        )

    def test_extract_data_factor_multiply(self):
        factor_convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "data_column_dict": {
                "delay_time": {"column_name": "delay_time", "conversion_factor": 2},
                "normalized_yield": "probability",
            },
        }

        #
        _, data_dict, _ = extract_data(
            config=self.convolution_config,
            convolution_instruction=factor_convolution_instruction,
        )

        #
        np.testing.assert_array_equal(
            data_dict["delay_time"], 2 * self.dummy_data["delay_time"] * u.yr
        )

    def test_extract_data_function_multiply(self):
        ###########
        # function multiplying
        function_convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "dummy",
            "output_data_name": "dummy",
            "data_column_dict": {
                "delay_time": {
                    "column_name": "delay_time",
                    "conversion_function": lambda x: x**2,
                },
                "normalized_yield": "probability",
            },
        }

        #
        _, data_dict, _ = extract_data(
            config=self.convolution_config,
            convolution_instruction=function_convolution_instruction,
        )

        #
        np.testing.assert_array_equal(
            data_dict["delay_time"], (self.dummy_data["delay_time"] ** 2) * u.yr
        )

    def test_extract_data_not_existing(self):
        ###########
        # Non existent
        faulty_convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "dummy2",
            "output_data_name": "dummy",
            "data_column_dict": {
                "delay_time": {
                    "column_name": "delay_time",
                    "conversion_function": lambda x: x**2,
                },
                "normalized_yield": "probability",
            },
        }

        with self.assertRaises(KeyError):

            #
            _, data_dict, _ = extract_data(
                config=self.convolution_config,
                convolution_instruction=faulty_convolution_instruction,
            )


class test_create_time_bin_info_dict(unittest.TestCase):
    def test_create_time_bin_info_dict(self):
        """ """

        time_bin_info_dict = create_time_bin_info_dict(
            config=default_convolution_config,
            convolution_instruction=default_convolution_instruction,
            bin_number=1,
            bin_edge_lower=0.5,
            bin_center=0.75,
            bin_size=1,
            bin_type="lookback",
        )

        self.assertEqual(
            time_bin_info_dict["time_type"], default_convolution_config["time_type"]
        )
        self.assertEqual(
            time_bin_info_dict["convolution_direction"],
            default_convolution_instruction["convolution_direction"],
        )
        self.assertEqual(
            time_bin_info_dict["reverse_bin_order"],
            default_convolution_instruction["reverse_convolution"],
        )


class test_sample_around_bin_center(unittest.TestCase):
    def test_sample_around_bin_center_basic(self):
        """ """

        bin_edges = np.array([0, 1])
        values = np.repeat(0.5, 10)

        sampled_values = sample_around_bin_center(bin_edges, values)

        expected_values = [
            0.5488135,
            0.71518937,
            0.60276338,
            0.54488318,
            0.4236548,
            0.64589411,
            0.43758721,
            0.891773,
            0.96366276,
            0.38344152,
        ]

        np.testing.assert_array_almost_equal(sampled_values, expected_values)

    def test_sample_around_bin_center_shifted_multiple_values(self):
        """ """

        bin_edges = np.array([2.5, 7.5, 12.5])
        values = np.concatenate([np.repeat(5, 10), np.repeat(10, 10)])

        sampled_values = sample_around_bin_center(bin_edges, values)

        expected_values = [
            6.45862519,
            5.1444746,
            5.34022281,
            7.12798319,
            2.85518029,
            2.9356465,
            2.60109199,
            6.66309923,
            6.39078375,
            6.85006074,
            12.39309171,
            11.49579282,
            9.80739681,
            11.40264588,
            8.09137213,
            10.69960511,
            8.21676644,
            12.22334459,
            10.10924161,
            9.5733097,
        ]

        np.testing.assert_array_almost_equal(sampled_values, expected_values)


class test_create_job_dict(unittest.TestCase):

    def test_create_job_dict(self):
        """ """

        job_dict = create_job_dict(
            config=default_convolution_config,
            convolution_instruction=default_convolution_instruction,
            time_bin_info_dict={"bin_number": 1},
            bin_number=1,
            sfr_dict={"lookback_time_edges": np.array([1, 2])},
            data_dict={"dummy": np.array([1, 2])},
        )

        self.assertTrue("job_number" in job_dict.keys())

        self.assertTrue("output_dir" in job_dict.keys())
        self.assertEqual(job_dict["output_dir"], "/tmp/sspc/input_data/output_data")


class test_get_physical_dimensions(unittest.TestCase):
    """ """

    def test_get_physical_dimensions_dimensionless(self):

        dimensionless_unit = u.m / u.m

        physical_dimensions = get_physical_dimensions(dimensionless_unit)

        self.assertEqual(physical_dimensions, "[dimensionless]")

    def test_get_physical_dimensions_starformation_rate(self):
        starformation_rate = u.Msun / u.yr

        physical_dimensions = get_physical_dimensions(starformation_rate)

        self.assertEqual(physical_dimensions, "[M][T^-1]")

    def test_get_physical_dimensions_mixed_base(self):
        mixed_base = u.cm / u.Msun

        physical_dimensions = get_physical_dimensions(mixed_base)

        self.assertEqual(physical_dimensions, "[L][M^-1]")


class test_generate_boilerplate_outputfile(unittest.TestCase):
    """ """

    def test_generate_boilerplate_outputfile(self):

        output_filename = os.path.join(
            TMP_DIR, "test_generate_boilerplate_outputfile.hdf5"
        )
        generate_boilerplate_outputfile(output_filename)

        # check if it is a file
        self.assertTrue(os.path.isfile(output_filename))

        with h5py.File(output_filename, "r") as output_hdf5_file:

            self.assertTrue("input_data" in output_hdf5_file.keys())
            self.assertTrue("config" in output_hdf5_file.keys())


class test_get_normalized_yield_unit(unittest.TestCase):
    """ """

    def test_get_normalized_yield_unit_no_normalized_yield(self):

        #
        tmp_convolution_config = copy.copy(default_convolution_config)
        tmp_convolution_instruction = copy.copy(default_convolution_instruction)

        with self.assertRaises(ValueError):
            get_normalized_yield_unit(
                tmp_convolution_config, tmp_convolution_instruction
            )

    def test_get_normalized_yield_unit_default(self):

        #
        tmp_convolution_config = copy.copy(default_convolution_config)
        tmp_convolution_instruction = copy.copy(default_convolution_instruction)
        tmp_convolution_instruction["data_column_dict"] = {
            "normalized_yield": "normalized_yield"
        }

        #
        unit = get_normalized_yield_unit(
            tmp_convolution_config, tmp_convolution_instruction
        )
        expected_unit = 1 * 1 / u.Msun

        #
        self.assertEqual(unit, expected_unit)

    def test_get_normalized_yield_unit_custom(self):

        #
        tmp_convolution_config = copy.copy(default_convolution_config)
        tmp_convolution_instruction = copy.copy(default_convolution_instruction)
        tmp_convolution_instruction["data_column_dict"] = {
            "normalized_yield": {"name": "normalized_yield", "unit": u.Msun}
        }

        #
        unit = get_normalized_yield_unit(
            tmp_convolution_config, tmp_convolution_instruction
        )
        expected_unit = u.Msun

        #
        self.assertEqual(unit, expected_unit)


class test_extract_unit_dict(unittest.TestCase):
    """ """

    def test_extract_unit_dict(self):

        tmp_output_filename = os.path.join(TMP_DIR, "test_extract_unit_dict.h5py")
        groupname = "test_extract_unit_dict"

        #######
        # Store data
        unit_dict = {"a": u.Msun}

        with h5py.File(tmp_output_filename, "a") as output_hdf5file:
            output_hdf5file.create_group(groupname)

            output_hdf5file["test_extract_unit_dict"].attrs["units"] = json.dumps(
                unit_dict, cls=JsonCustomEncoder
            )

        #######
        # Store data
        with h5py.File(tmp_output_filename, "r") as output_hdf5file:
            read_unit_dict = extract_unit_dict(output_hdf5file, groupname)

            self.assertTrue(read_unit_dict == unit_dict)


class test_is_mass_unit(unittest.TestCase):
    """ """

    def test_is_mass_unit(self):
        mass_unit_value = 1 * u.g

        self.assertTrue(is_mass_unit(mass_unit_value))

    def test_is_not_mass_unit(self):
        no_unit_value = 1
        self.assertFalse(is_mass_unit(no_unit_value))

    def test_is_unit_but_not_mass_unit(self):
        wrong_unit_value = 1 * u.yr
        self.assertFalse(is_mass_unit(wrong_unit_value))


class test_is_time_unit(unittest.TestCase):
    """ """

    def test_is_time_unit(self):
        time_unit_value = 1 * u.s

        self.assertTrue(is_time_unit(time_unit_value))

    def test_is_not_time_unit(self):
        no_unit_value = 1
        self.assertFalse(is_time_unit(no_unit_value))

    def test_is_unit_but_not_time_unit(self):
        wrong_unit_value = 1 * u.m
        self.assertFalse(is_time_unit(wrong_unit_value))


class test_has_unit(unittest.TestCase):
    """ """

    def test_unit(self):
        unit_value = 1 * u.m

        self.assertTrue(has_unit(unit_value))

    def test_no_unit(self):
        no_unit_value = 1

        self.assertFalse(has_unit(no_unit_value))

    def test_dimensionless_unit(self):

        dimensionless_unit = u.m / u.m

        dimensionless_value = 1 * dimensionless_unit

        self.assertTrue(has_unit(dimensionless_value, fail_on_dimensionless=False))

        self.assertFalse(has_unit(dimensionless_value, fail_on_dimensionless=True))


class test_get_username(unittest.TestCase):
    """ """

    def test_get_username(self):
        username = get_username()

        # should be a string
        self.assertTrue(isinstance(username, str))

        # should be of some lenght
        self.assertTrue(len(username) > 0)


class test_temp_dir(unittest.TestCase):
    """
    Unittests for temp_dir
    """

    def test_create_temp_dir(self):
        """
        Test making a temp directory and comparing that to what it should be
        """

        #
        username = get_username()
        general_temp_dir = tempfile.gettempdir()

        # Get username
        username = get_username()
        sspc_temp_dir = os.path.join(temp_dir())

        #
        self.assertTrue(
            os.path.isdir(os.path.join(general_temp_dir, "sspc-{}".format(username)))
        )
        self.assertTrue(
            os.path.join(general_temp_dir, "sspc-{}".format(username)) == sspc_temp_dir
        )


class test_check_required(unittest.TestCase):
    def setUp(self):
        self.config = {
            "input_shape": (32, 32, 3),
            "output_shape": (10,),
            "learning_rate": 0.001,
        }

    def test_check_required_all_present(self):
        required_list = ["input_shape", "output_shape", "learning_rate"]
        check_required(self.config, required_list)
        # No exception should be raised

    def test_check_required_missing_key(self):
        required_list = ["input_shape", "output_shape", "learning_rate", "batch_size"]
        with self.assertRaises(ValueError):
            check_required(self.config, required_list)

    def test_check_required_empty_list(self):
        required_list = []
        check_required(self.config, required_list)
        # No exception should be raised


class test_handle_custom_scaling_or_conversion(unittest.TestCase):
    def setUp(self):
        self.data_layer_dict = {
            "factor": {"layer_depth": 2, "conversion_factor": 2},
            "function": {"layer_depth": 4, "conversion_function": lambda x: x**2},
            "both": {
                "layer_depth": 4,
                "conversion_function": lambda x: x**2,
                "conversion_factor": 2,
            },
        }

        self.logger = logging.getLogger(__name__)
        FORMAT = "[%(filename)s:%(lineno)s - %(funcName)20s ] %(asctime)s: %(message)s"
        logging.basicConfig(format=FORMAT)
        self.logger.setLevel(logging.INFO)

    def test_factor_array(self):
        array = handle_custom_scaling_or_conversion(
            config={"logger": self.logger},
            data_layer_or_column_dict_entry=self.data_layer_dict["factor"],
            value=np.array([1, 2]),
        )
        np.testing.assert_array_equal(array, np.array([2, 4]))

    def test_function_array(self):
        array = handle_custom_scaling_or_conversion(
            config={"logger": self.logger},
            data_layer_or_column_dict_entry=self.data_layer_dict["function"],
            value=np.array([1, 2]),
        )
        np.testing.assert_array_equal(array, np.array([1, 4]))

    def test_factor_scalar(self):
        value = handle_custom_scaling_or_conversion(
            config={"logger": self.logger},
            data_layer_or_column_dict_entry=self.data_layer_dict["factor"],
            value=1,
        )
        self.assertEqual(value, 2)

    def test_function_scalar(self):
        value = handle_custom_scaling_or_conversion(
            config={"logger": self.logger},
            data_layer_or_column_dict_entry=self.data_layer_dict["function"],
            value=2,
        )
        self.assertEqual(value, 4)

    def test_get_deepest_data_layer_depth_both(self):
        with self.assertRaises(ValueError):
            handle_custom_scaling_or_conversion(
                config={"logger": self.logger},
                data_layer_or_column_dict_entry=self.data_layer_dict["both"],
                value=2,
            )


class test_calculate_bincenters(unittest.TestCase):
    def test_calculate_bincenters_linear(self):
        array = np.array([1.0, 2, 3, 4, 5])
        expected_bincenters = np.array([1.5, 2.5, 3.5, 4.5])
        bincenters = calculate_bincenters(array, convert="linear")
        np.testing.assert_array_equal(bincenters, expected_bincenters)


class test_calculate_bin_edges(unittest.TestCase):
    def test_calculate_bin_edges(self):
        arr = np.array([1.0, 2, 3, 4, 5])
        expected_bin_edges = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
        bin_edges = calculate_bin_edges(arr)
        np.testing.assert_array_equal(bin_edges, expected_bin_edges)


class test_pad_function(unittest.TestCase):
    def test_pad_function_relative_to_edge_val_axis_0(self):
        array = np.array([1.0, 2, 3, 4, 5])
        left_val = -0.5
        right_val = 0.5
        relative_to_edge_val = True
        expected_padded_array = np.array([0.5, 1, 2, 3, 4, 5, 5.5])
        padded_array = pad_function(
            array, left_val, right_val, relative_to_edge_val, axis=0
        )
        np.testing.assert_array_equal(padded_array, expected_padded_array)

    def test_pad_function_absolute_axis_1(self):
        array = np.array([[1.0, 2, 3], [4, 5, 6]])
        left_val = 0
        right_val = 0
        relative_to_edge_val = False
        expected_padded_array = np.array([[0, 1, 2, 3, 0], [0, 4, 5, 6, 0]])
        padded_array = pad_function(
            array, left_val, right_val, relative_to_edge_val, axis=1
        )
        np.testing.assert_array_equal(padded_array, expected_padded_array)


class test_generate_group_name(unittest.TestCase):
    def setUp(self):
        self.convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "input_image",
            "output_data_name": "output_image",
        }
        self.sfr_dict = {"name": "test_group"}

    def test_generate_group_name_with_sfr(self):
        groupname, elements = generate_group_name(
            self.convolution_instruction, self.sfr_dict
        )
        expected_groupname = "test_group/input_image/output_image"
        expected_elements = ["test_group", "input_image", "output_image"]
        self.assertEqual(groupname, expected_groupname)
        self.assertListEqual(elements, expected_elements)

    def test_generate_group_name_without_sfr(self):
        groupname, elements = generate_group_name(self.convolution_instruction, {})
        expected_groupname = "input_image/output_image"
        expected_elements = ["input_image", "output_image"]
        self.assertEqual(groupname, expected_groupname)
        self.assertListEqual(elements, expected_elements)


class test_get_tmp_dir(unittest.TestCase):
    def setUp(self):
        self.convolution_instruction = {
            **default_convolution_instruction,
            "input_data_name": "input_image",
            "output_data_name": "output_image",
        }

    def test_get_tmp_dir(self):
        tmp_dir = get_tmp_dir(
            config={"tmp_dir": TMP_DIR},
            convolution_instruction=self.convolution_instruction,
        )
        self.assertEqual(tmp_dir, os.path.join(TMP_DIR, "input_image/output_image"))


if __name__ == "__main__":
    unittest.main()
