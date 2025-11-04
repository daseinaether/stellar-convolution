"""
Main convolution test script
"""

# pylint: disable=W0611
# flake8: noqa
import unittest

from syntheticstellarpopconvolve.tests.general_tests import test_postprocessing
from syntheticstellarpopconvolve.tests.test_calculate_birth_redshift_array import (
    test_calculate_origin_redshift_array,
)
from syntheticstellarpopconvolve.tests.test_calculate_starformation_rate import (
    test_calculate_digitized_sfr_rates_binned_data_for_backward_convolution,
    test_calculate_digitized_sfr_rates_non_binned_data_for_backward_convolution,
    test_calculate_origin_time_array,
    test_general_sfr_digitise_function,
)
from syntheticstellarpopconvolve.tests.test_check_and_prepare_output_file import (
    test_check_and_prepare_output_file,
)
from syntheticstellarpopconvolve.tests.test_check_and_update_convolution_config import (
    test_check_convolution_config,
    test_update_convolution_config,
)
from syntheticstellarpopconvolve.tests.test_check_and_update_convolution_instruction import (
    test_check_and_update_convolution_instructions,
    test_check_convolution_instruction,
    test_check_delay_time_data_bin_info_dict,
    test_check_metallicity,
)
from syntheticstellarpopconvolve.tests.test_check_and_update_sfr_dict import (
    test_check_sfr_dict,
    test_pad_sfr_dict,
    test_update_sfr_dict,
)
from syntheticstellarpopconvolve.tests.test_convolution_by_integration import (
    test_convolution_by_integration_post_convolution_hook_wrapper,
)
from syntheticstellarpopconvolve.tests.test_convolution_by_sampling import (
    test_convolution_by_sampling_post_convolution_hook_wrapper,
    test_sample_systems,
    test_select_dict_entries_with_new_indices,
)
from syntheticstellarpopconvolve.tests.test_convolution_with_events import (
    test_convolution_with_events,
)
from syntheticstellarpopconvolve.tests.test_convolve_binned_data import (
    test_calculate_overlap_fractions,
)
from syntheticstellarpopconvolve.tests.test_convolve_binned_data_with_backward_convolution import (
    test_convolve_binned_data_with_backward_convolution,
)
from syntheticstellarpopconvolve.tests.test_convolve_nonbinned_data_with_backward_convolution import (
    test_convolve_nonbinned_data_with_backward_convolution,
)
from syntheticstellarpopconvolve.tests.test_convolve_on_the_fly import (
    test_convolve_on_the_fly,
    test_convolve_on_the_fly_post_convolution_hook_wrapper,
    test_handle_call_on_the_fly_function,
)
from syntheticstellarpopconvolve.tests.test_cosmology_utils import (
    test_age_of_universe_to_redshift,
    test_lookback_time_to_redshift,
    test_redshift_to_age_of_universe,
    test_redshift_to_lookback_time,
)
from syntheticstellarpopconvolve.tests.test_default_convolution_config import (
    test_array_validation,
    test_callable_or_none_validation,
    test_callable_validation,
    test_dict_or_list_of_dicts_validation,
    test_existing_path_validation,
    test_list_of_dicts_validation,
    test_logger_validation,
    test_unit_validation,
)
from syntheticstellarpopconvolve.tests.test_ensemble_utils import (
    test_convert_ensemble_to_dataframe,
    test_find_columnames_recursively,
    test_flatten_data_ensemble1d,
    test_inflate_ensemble_with_lists_and_named_layers,
    test_inflate_ensemble_with_lists_without_named_layers,
)
from syntheticstellarpopconvolve.tests.test_general_functions import (
    test_calculate_bin_edges,
    test_calculate_bincenters,
    test_check_required,
    test_create_job_dict,
    test_create_time_bin_info_dict,
    test_extract_data,
    test_extract_unit_dict,
    test_generate_boilerplate_outputfile,
    test_generate_data_dict,
    test_generate_group_name,
    test_get_normalized_yield_unit,
    test_get_physical_dimensions,
    test_get_tmp_dir,
    test_get_username,
    test_handle_custom_scaling_or_conversion,
    test_has_unit,
    test_is_mass_unit,
    test_is_time_unit,
    test_pad_function,
    test_sample_around_bin_center,
    test_temp_dir,
)
from syntheticstellarpopconvolve.tests.test_post_convolution_hook_routines import (
    test_extract_arguments,
    test_handle_post_convolution_function,
)
from syntheticstellarpopconvolve.tests.test_prepare_redshift_interpolator import (
    test_create_interpolation_datasets,
    test_load_interpolation_data,
    test_prepare_redshift_interpolator,
)
from syntheticstellarpopconvolve.tests.test_store_redshift_shell_info import (
    test_create_shell_volume_dict,
    test_store_redshift_shell_info,
)

# from syntheticstellarpopconvolve.tests.test_check_and_update_convolution_config import (
#     test_check_convolution_config,
#     test_check_convolution_instruction,
#     test_check_metallicity,
#     test_check_required,
# )


if __name__ == "__main__":
    unittest.main()
