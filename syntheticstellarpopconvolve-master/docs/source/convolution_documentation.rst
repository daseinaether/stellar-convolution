Convolution options
===================

The following page contains documentation on the convolution-config options and the convolution-instructions options.

The convolution-config dictionary provides the `global` configuration of the convolution code. This dictionary has to be passed to the main entrypoint function as `convolve(config=convolution_config)`.

The convolution-instruction dictionary provides configuration on a per-convolution basis, allowing the main convolution function to perform a series of convolutions with different configurations. This dictionary has to be passed to the convolution-config function as `convolution_config['convolution_instructions'] = [{convolution_instruction_dict_1, ...}]`.

Convolution-config options
--------------------------


.. list-table:: Convolution-config options
   :widths: 25, 75
   :header-rows: 1

   * - Option
     - Description
   * - SFR_info
     - Description:
          Dictionary containing the starformation rate info. Can also be a list of dictionaries.
       
       Default value:
          {}
   * - check_convolution_config
     - Description:
          Flag whether to validate the configuration dictionary before running the convolution code.
       
       Default value:
          True
   * - convolution_instructions
     - Description:
          List of instructions for the convolution.
       
       Default value:
          [{}]
   * - convolution_lookback_time_bin_edges
     - Description:
          Lookback-time bin-edges used in convolution.
       
       Default value:
          None
   * - convolution_redshift_bin_edges
     - Description:
          Redshift bin-edges used in convolution.
       
       Default value:
          None
   * - cosmology
     - Description:
          Astropy cosmology used throughout the code.
       
       Default value:
          FlatLambdaCDM(name="Planck13", H0=67.77 km / (Mpc s), Om0=0.30712, Tcmb0=2.7255 K, Neff=3.046, m_nu=[0.   0.   0.06] eV, Ob0=0.048252)
   * - custom_convolution_function
     - Description:
          
       
       Default value:
          None
   * - custom_data_extraction_function
     - Description:
          
       
       Default value:
          None
   * - custom_rates_function
     - Description:
          Custom rate function used in the convolution.
       
       Default value:
          None
   * - default_normalized_yield_unit
     - Description:
          Default unit used for the normalized-yield data.
       
       Default value:
          1.0 1 / solMass
   * - delay_time_default_unit
     - Description:
          Default unit used for the delay-time data. NOTE: this can be overridden in data_dict column or layer entries.
       
       Default value:
          yr
   * - include_custom_rates
     - Description:
          Whether to include custom, user-specified, rates. See 'custom_rates_function'.
       
       Default value:
          False
   * - logger
     - Description:
          Logger object.
       
       Default value:
          <Logger syntheticstellarpopconvolve.default_convolution_config (CRITICAL)>
   * - max_job_queue_size
     - Description:
          Max number of jobs in the multiprocessing queue for the convolution.
       
       Default value:
          8
   * - multiprocessing
     - Description:
          Flag whether to enable multiprocessing. True for multiprocessing, which allows faster convolution but does not allow the use of the previous convolution results and the persistent data. False for sequential convolution, which is slower but previous convolution results and the persistent data is available here.
       
       Default value:
          True
   * - num_cores
     - Description:
          Number of cores to use to do the convolution.
       
       Default value:
          1
   * - output_filename
     - Description:
          Full path to output hdf5 filename. This should point to a file that already contains the input data that will be used to do the convolution with.
       
       Default value:
          
   * - redshift_interpolator_data_output_filename
     - Description:
          Filename for the redshift interpolator object.
       
       Default value:
          None
   * - redshift_interpolator_force_rebuild
     - Description:
          Whether to force rebuild the redshift interpolator.
       
       Default value:
          False
   * - redshift_interpolator_max_redshift
     - Description:
          Minimum redshift for the redshift interpolator.
       
       Default value:
          50
   * - redshift_interpolator_min_redshift
     - Description:
          Minimum redshift for the redshift interpolator.
       
       Default value:
          0
   * - redshift_interpolator_min_redshift_if_log
     - Description:
          Minimum redshift for the redshift interpolator if using log spacing.
       
       Default value:
          1e-05
   * - redshift_interpolator_rebuild_when_settings_mismatch
     - Description:
          Whether to rebuild the redshift interpolator when the config of the existing one don't match with the current config.
       
       Default value:
          True
   * - redshift_interpolator_stepsize
     - Description:
          Stepsize for the redshift interpolation.
       
       Default value:
          0.001
   * - redshift_interpolator_use_log
     - Description:
          Whether to interpolate in log redshift.
       
       Default value:
          True
   * - remove_pickle_files
     - Description:
          Flag whether to remove all the pickle files after writing them to the main hdf5 file.
       
       Default value:
          True
   * - time_type
     - Description:
          Time-type used in convolution. Can be either 'redshift' or 'lookback_time'.
       
       Default value:
          lookback_time
   * - tmp_dir
     - Description:
          Target directory for the tmp files.
       
       Default value:
          /tmp/sspc
   * - write_to_hdf5
     - Description:
          Whether to write the pickle-files from the convolution back to the main hdf5 file.
       
       Default value:
          True


Convolution-instruction options
-------------------------------


.. list-table:: Convolution-instruction options
   :widths: 25, 75
   :header-rows: 1

   * - Option
     - Description
   * - chunk_size
     - Description:
          Chunk size for the data readout.
       
       Default value:
          0
   * - chunk_total
     - Description:
          Total number of chunks to be considered. Should be an integer rounded up calculated as ceil(<total entries in input dataframe>/<chunk size>).
       
       Default value:
          0
   * - chunked_readout
     - Description:
          Flag to read the input data in chunks. See `chunk_size`.
       
       Default value:
          False
   * - contains_binned_data
     - Description:
          Flag to indicate whether the input data is binned (in time). If so, the user should provide additional information.
       
       Default value:
          False
   * - convolution_direction
     - Description:
          Choice of convolution direction. 'backward' convolves the data such that every event occurs at the current time by looking what the starformation rate is for each given delay time. 'forward' convolution generates each system at the same time and looks at when events happen afterwards based on their delay times. Note: neither option is supported in all choices of `convolution_type` and `contains_binned_data`.
       
       Default value:
          backward
   * - convolution_type
     - Description:
          Method of convolution. The three choices are as follows. 'integrate': Convolution by integration uses backward convolution to multiply the normalized_yield of the systems with the starformation rate at the time the system would be born, given the delay time and the target event time. This is particularly useful when you are just interested in the (transient) event. 'sample': Convolution by sampling uses forward convolution to 'sample' systems according to the yield () and assigns an event-time based on the delay time and the birth time. This is particularly useful if you want to post-process systems after they are born/the event occurs, like integrating the orbit of double compact object forward in time due to gravitational wave radiation (see LISA project example in `examples/notebook_usecases`). 'on-the-fly': Convolution by simulating the systems on-the-fly. Requires a method that uses the total mass in star formation, and optionally the metallicity distribution, combined with a population synthesis code, to evolve systems on the fly.
       
       Default value:
          integrate
   * - data_column_dict
     - Description:
          Dictionary containing the mapping between the names of the columns in the pandas dataframe of the input data and the names as they are used in the convolution framework, as `{<framework_data_name>: <pandas_column_name>}`. Mappings for TODO are required. Any extra columns that are provided are accessible by the post-convolution function. Entries in this dictionary can be either names or dictionaries themselves, allowing more advanced functionality. See examples/notebook_convolution_advanced.
       
       Default value:
          {}
   * - data_time_bin_info
     - Description:
          Dictionary containing data time-bin information when convolving binned data.
       
       Default value:
          {}
   * - input_data_name
     - Description:
          Name of to the current input dataset. Will be used to extract the data from the provided input-hdf5 file (expected in /input_data/<input data name>/), and will be used in the output-data path.
       
       Default value:
          input_data
   * - multiply_by_convolution_time_binsize
     - Description:
          Flag to multiply the convolution results by the convolution time-bin size. Not supported when time_type=='redshift'.
       
       Default value:
          False
   * - multiply_by_sfr_time_binsize
     - Description:
          Flag to multiply the convolution results by the starformation rate time bin size. Not supported when time_type=='redshift'.
       
       Default value:
          False
   * - output_data_name
     - Description:
          Name assigned to the current output dataset. Will be used in the output-data path. Can be useful when running a convolution with the same input-data but e.g. with a different post-convolution function.
       
       Default value:
          output_data
   * - post_convolution_function
     - Description:
          Function that performs post-convolution operations on the convolved data, like applying detection probability weights, further integration of systems or just general filtering of the data. Different `convolution_type`s allows for different modifications of the data, in that convolution of `event-based` data by 'sampling' allows TODO. The arguments of this function should be chosen from: 'config', 'time_value', 'convolution_instruction', 'data_dict' and the contents of 'post_convolution_function_extra_parameters'. For more explanation about this function see the convolution notebook.
       
       Default value:
          None
   * - post_convolution_function_extra_parameters
     - Description:
          Dictionary containing additional arguments that can be accessed by the extra_weights_function. Note: using this can often be avoided by using functools.partial to fix certain input parameters.
       
       Default value:
          {}
   * - reverse_convolution
     - Description:
          Flag to reverse the convolution direction. If True, we start with the bin furthest back in time and work to. Useful in combination with `convolution_config['multiprocessing']=False`, and the `previous_convolution_results` and `persistant_data` objects.
       
       Default value:
          False


