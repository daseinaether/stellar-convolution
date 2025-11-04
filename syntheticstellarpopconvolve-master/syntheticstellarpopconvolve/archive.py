# def create_bins_from_centers(centers):
#     """
#     Function to create a set of bin edges from a set of bin centers. Assumes the two endpoints have the same binwidth as their neighbours
#     """

#     # Create bin edges minus the outer two
#     bin_edges = (centers[1:] + centers[:-1]) / 2

#     # Add to left
#     bin_edges = np.append(np.array(bin_edges[0] - np.diff(centers)[0]), bin_edges)

#     # Add to right
#     bin_edges = np.append(bin_edges, np.array(bin_edges[-1] + np.diff(centers)[-1]))

#     return bin_edges


# def create_extended_time_bins(config):
#     """
#     Function to create extended time bins so we
#     """

#     # Extend the time bins to have the digitize put the stuff that falls out of the bounds in a column with 0 SFR
#     extended_time_bins = copy.copy(config["time_bins"])
#     extended_time_bins = np.insert(
#         extended_time_bins,
#         0,
#         extended_time_bins[0] - np.diff(config["time_bins"])[0],
#     )
#     extended_time_bins = np.insert(
#         extended_time_bins, extended_time_bins.size, extended_time_bins[-1] + 1e100
#     )

#     return extended_time_bins


# def create_metallicity_redshift_dataframe(
#     amt_z_values, amt_metallicity_values, z_values, metallicities, config_dict_cosmology
# ):
#     """
#     Function to create the metallicity redshift dataframe
#     """

#     solar_value = config_dict_cosmology["solar_value"]
#     min_value_metalprob = config_dict_cosmology["min_value_metalprob"]

#     #
#     all_redshifts = np.zeros(amt_z_values * amt_metallicity_values)
#     all_metallicity_values = np.zeros(amt_z_values * amt_metallicity_values)
#     all_probabilities = np.zeros(amt_z_values * amt_metallicity_values)

#     #
#     for i, redshift in enumerate(z_values):
#         metallicity_distribution_values = metallicity_distribution(
#             metallicities, redshift, config_dict_cosmology
#         )
#         normed_metallicity_distribution_values = (
#             metallicity_distribution_values / np.sum(metallicity_distribution_values)
#         )

#         normed_metallicity_distribution_values[
#             normed_metallicity_distribution_values < min_value_metalprob
#         ] = 0
#         renormed_normed_metallicity_distribution_values = (
#             normed_metallicity_distribution_values
#             / np.sum(normed_metallicity_distribution_values)
#         )

#         all_redshifts[
#             i * amt_metallicity_values : (i + 1) * amt_metallicity_values
#         ] = redshift
#         all_metallicity_values[
#             i * amt_metallicity_values : (i + 1) * amt_metallicity_values
#         ] = (metallicities / solar_value)
#         all_probabilities[
#             i * amt_metallicity_values : (i + 1) * amt_metallicity_values
#         ] = renormed_normed_metallicity_distribution_values

#     # Put the metallicity evolution in a dataframe
#     new_dataframe = pd.DataFrame(
#         data=np.array([all_redshifts, all_metallicity_values, all_probabilities]).T,
#         columns=["redshift", "metallicity", "probability"],
#     )

#     return new_dataframe


# def create_metallicity_lookback_dataframe(
#     amt_lookback_values,
#     amt_metallicity_values,
#     lookback_values,
#     metallicities,
#     config_dict_cosmology,
# ):
#     """
#     Function to create the metallicity redshift dataframe
#     """

#     solar_value = config_dict_cosmology["solar_value"]
#     min_value_metalprob = config_dict_cosmology["min_value_metalprob"]

#     #
#     all_lookbacks = np.zeros(amt_lookback_values * amt_metallicity_values)
#     all_metallicity_values = np.zeros(amt_lookback_values * amt_metallicity_values)
#     all_probabilities = np.zeros(amt_lookback_values * amt_metallicity_values)

#     for i, lookback in enumerate(lookback_values):

#         # Turn it around
#         universe_age = cosmo.age(0).value - lookback

#         # calc redshift
#         redshift = astropy.cosmology.z_at_value(cosmo.age, float(universe_age) * u.Gyr)

#         # Calculate rest of the values
#         metallicity_distribution_values = metallicity_distribution(
#             metallicities, redshift, config_dict_cosmology
#         )
#         normed_metallicity_distribution_values = (
#             metallicity_distribution_values / np.sum(metallicity_distribution_values)
#         )

#         normed_metallicity_distribution_values[
#             normed_metallicity_distribution_values < min_value_metalprob
#         ] = 0
#         renormed_normed_metallicity_distribution_values = (
#             normed_metallicity_distribution_values
#             / np.sum(normed_metallicity_distribution_values)
#         )

#         all_lookbacks[
#             i * amt_metallicity_values : (i + 1) * amt_metallicity_values
#         ] = lookback
#         all_metallicity_values[
#             i * amt_metallicity_values : (i + 1) * amt_metallicity_values
#         ] = (metallicities / solar_value)
#         all_probabilities[
#             i * amt_metallicity_values : (i + 1) * amt_metallicity_values
#         ] = renormed_normed_metallicity_distribution_values

#     # Put the metallicity evolution in a dataframe
#     new_dataframe = pd.DataFrame(
#         data=np.array([all_lookbacks, all_metallicity_values, all_probabilities]).T,
#         columns=["lookback", "metallicity", "probability"],
#     )

#     return new_dataframe


# """
# Star formation rate and metallicity distribution functions along with some other cosmology functions
# """

# import numpy as np
# from grav_waves.gw_analysis.functions.cosmology_functions import starformation_rate
# from scipy.stats import norm as NormDist
# import pandas as pd
# from david_phd_functions.plotting import custom_mpl_settings
# from grav_waves.settings import cosmo

# custom_mpl_settings.load_mpl_rc()

# import astropy

# # def madau_dickinson_sfr(z, a, b, c, d):
# #     """
# #     Function from Neijsel et al 2019 (https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.3740N/abstract)
# #     """

# #     rate = a * np.power(1 + z, b) / (1 + np.power((1 + z) / c, d))

# #     return rate * (u.Msun / (u.yr * (u.Mpc**3)))


# def starformation_rate_dummy(constant=1):
#     """
#     Dummy starformation rate function. returns <constant> * (u.Msun/(u.yr * (u.Mpc**3)))
#     """

#     return constant * (u.Msun / (u.yr * (u.Gpc**3)))


# def starformation_rate(z, cosmology_configuration, verbosity=0):
#     """
#     Function to determine the star formation rate.

#     Input:
#         z: redshift value
#         cosmology_configuration: dict containing the choice of star formation rate function ('star_formation_rate_function') and a set of arguments ('star_formation_rate_args')

#     Output:
#         starformation rate in u.Msun/(u.yr * (u.Gpc**3))
#     """

#     if cosmology_configuration["star_formation_rate_function"] == "madau_dickinson_sfr":
#         return starformation_rate_coen19(
#             z, **cosmology_configuration["star_formation_rate_args"]
#         ).to(u.Msun / (u.yr * (u.Gpc**3)))
#     elif cosmology_configuration["star_formation_rate_function"] == "dummy":
#         return starformation_rate_dummy().to(u.Msun / (u.yr * (u.Gpc**3)))
#     else:
#         raise ValueError("Unknown SFR function")

# def metallicity_distribution(metallicity, redshift, cosmology_configuration):
#     """
#     Function that calculates the metallicity probability distribution for a given metallicity and redshift

#     Will call the chosen function with the arguments (metallicity, redshift, **cosmology_configuration['metallicity_distribution_args'])

#     Input:
#         metallicity: metallicty (non-log) value
#         redshift: current redshift
#         cosmology_configuration: dict containing information on the choice of prescription: 'metallicity_distribution_function' and a set of arguments for that function ('metallicity_distribution_args')

#     Output:
#         metallicity probability distribution: dP/dZ (will always be in in linear Z, not in log Z)
#     """

#     if cosmology_configuration["metallicity_distribution_function"] == "coen19":
#         return metallicity_distribution_coen19(
#             metallicity,
#             redshift,
#             **cosmology_configuration["metallicity_distribution_args"],
#         )
#     elif cosmology_configuration["metallicity_distribution_function"] == "dummy":
#         return metallicity_distribution_dummy()
#     elif cosmology_configuration["metallicity_distribution_function"] == "vanSon21":
#         # TODO: implement Z-dist lieke
#         raise NotImplementedError(
#             "van Son et al 2021 metallicity distribution is not implemented yet"
#         )
#     else:
#         raise ValueError("metallicity distribution option not known")
#
# def normalize_metallicities_distribution(
#     redshift_value,
#     cosmology_configuration,
#     min_log10_Z=-12,
#     max_log10_Z=0,
#     steps_log10_Z=10000,
# ):
#     """
#     Function that calculates the normalisation constant for the metallicities_distribution at a given redshift
#     """

#     #
#     stepsize_log10_Z = (max_log10_Z - min_log10_Z) / steps_log10_Z

#     # Do numpy:
#     # start_numpy = time.time()
#     logz_array = np.linspace(min_log10_Z, max_log10_Z, steps_log10_Z)
#     z_array = 10**logz_array

#     dPdZ_array = metallicity_distribution(
#         metallicity=z_array,
#         redshift=redshift_value,
#         cosmology_configuration=cosmology_configuration,
#     )
#     dPdlog10Z_array = dPdZ_array * np.log(10) * z_array

#     P_array = dPdlog10Z_array * stepsize_log10_Z

#     total_numpy = np.sum(P_array)
#     # end_numpy = time.time()

#     #     print(
#     #         """
#     # normalize_metallicities_distribution:
#     #         redshift_value: {}
#     #         min_log10_Z: {}
#     #         max_log10_Z: {}
#     #         steps_log10_Z: {}

#     #         loop method: totalP: {} (took {}s)
#     #         numpy method: totalP: {} (took {}s)
#     # """.format(
#     #         redshift_value,
#     #         min_log10_Z,
#     #         max_log10_Z,
#     #         steps_log10_Z,
#     #         total, end_loop-start_loop,
#     #         total_numpy, end_numpy-start_numpy,
#     #         )
#     #     )

#     return total_numpy


##############
# Old archive code

#     starformation_array = job_dict["sfr_dict"]["starformation_array"]

#     # TODO: i think i want to change this entirely:
#     # - Copy the ensemble structure locally.
#     # - do 'convolution' by just multiplying the entire subtree with a number
#     # - retain structure like the initial ensemble
#     # - extract ONLY the endpoint data as an array
#     # - store that in the file
#     # - use specialized functions to 'zip back' those endpoints onto the original ensemble.

#     #
#     extended_time_bins = create_extended_time_bins(config=config)

#     #
#     maxtime = redshift_to_lookback_time(
#         config["redshift_first_SFR"], cosmology=config["cosmology"]
#     )

#     # data = load_ensemble(ensemble_file)
#     ensemble_data = data_dict["ensemble_data"]
#     bse_options = config["population_settings"]["population_settings"]["bse_options"]

#     convolved_ensemble = {}

#     # Handle scalar ensemble convolution
#     if "scalars" in ensemble_data.keys():
#         convolved_ensemble["scalars"] = convolve_scalar_ensemble_data(
#             scalar_ensemble=ensemble_data["scalars"],
#             time_value=time_value,
#             job_dict=job_dict,
#             config=config,
#             convolution_instruction=convolution_instruction,
#         )


# def convolve_general_ensemble_data():
#     """
#     Sort of template for ensemble data convolution
#     """

#     # TODO: Copy original
#     # TODO: make another copy with all values set to 0
#     # TODO: Loop over metallicity if necessary
#     # TODO: Loop over time
#     # TODO: calculate SFR (time and metallicity if necessary)
#     # TODO: calculate timestep
#     # TODO: loop over all sub-ensembles to multiply each by either just SFR or SFR * dt
#     # TODO: merge with total
#     # strip endpoints to array
#     # if i == 0: also include stripped ensemble in output


# def convolve_scalar_ensemble_data(
#     scalar_ensemble, time_value, job_dict, config, convolution_instruction
# ):
#     """
#     Dedicated function to handle the convolution of ensemble data
#     """

#     #
#     bse_options = config["population_settings"]["population_settings"]["bse_options"]

#     # Get event types
#     eventtypes = [
#         s for s in scalar_ensemble if s.startswith("SN") or s.startswith("MERGER")
#     ]

#     # Get star types
#     startypes = [
#         s
#         for s in scalar_ensemble
#         if s.startswith("SPECTRAL_TYPE_")
#         or s.startswith("NUMBER_OF")
#         or s.startswith("TEST")
#     ]

#     #################
#     # Loop over metallicity
#     for metallicity in scalar_ensemble["metallicity"]:
#         print(metallicity)

#         print(scalar_ensemble["metallicity"][metallicity].keys())
#         print(eventtypes, startypes)

#     print(scalar_ensemble.keys())

#     # Set up times
#     logtimestep = bse_options["ensemble_logdt"]
#     minlogtime = bse_options["ensemble_startlogtime"]
#     maxlogtime = -1e30
#     for event in eventtypes:
#         for time in scalar_ensemble[event]:
#             maxlogtime = max(maxlogtime, time)
#     print(f"time {minlogtime} to {maxlogtime} logtimestep {logtimestep}")

#     quit()

#     #
#     rate = AutoVivificationDict()
#     number = AutoVivificationDict()
#     Mtot = 0
#     dtyearstot = 0

#     # loop over time
#     logtime = minlogtime
#     while logtime <= maxlogtime:

#         ###########
#         # TODO: package this in a function and handle the redshift conversion here as well

#         # Convert from log time to normal time
#         dlogt = logtime - logtimestep  # log time bin "width"
#         t = 10.0 ** (logtime + 6.0)  # linear time in yr
#         dtyears = t - 10.0 ** (dlogt + 6.0)  # linear time bin "width" in yr

#         # shift by time
#         t = t + data_dict["time_value"]

#         # correction factor if we exceed maxtime
#         if t > maxtime:
#             dtyears *= (t - maxtime) / dtyears
#         dtyearstot += dtyears

#         # Get indices for birth redshift
#         digitized_time_indices = (
#             np.digitize(np.array([t]), bins=extended_time_bins, right=False) - 1
#         )

#         # TODO: currently lets just take the starformation rate. we'll generalize this later to include metallicity
#         dM = starformation_array[digitized_time_indices]

#         print(dM)

#         # in this dtyears we formed dM (Msun) mass into stars
#         # dM = dtyears * SFR(t)
#         Mtot += dM

#         print(
#             f"logt={logtime:15g} t={t*1e-6:15g} Myr dlogt={dlogt:15g} dtyears={dtyears*1e-6:15g} Myr"
#         )

#         def convolve_ensemble_rates(scalar_ensemble, eventtypes, rate):
#             """ """
#             # rates: we have to take into account the binsize in time
#             for event in eventtypes:
#                 if logtime in scalar_ensemble[event]:
#                     # 1e-6 converts from per Myr to per year
#                     drate = 1e-6 * scalar_ensemble[event][logtime]

#                     rate[event] += drate * dM

#         def convolve_ensemble_number(scalar_ensemble, startype):
#             # numbers of stars
#             for startype in startypes:
#                 if logtime in scalar_ensemble[startype]:
#                     # binary_c's ensemble gives us the dp per unit mass into
#                     # stars in each bin. There is no time involved here,
#                     # so we don't need to convert from Myr to yr.
#                     dpdM = scalar_ensemble[startype][logtime]

#                     # Nstars = probability * SFR
#                     #
#                     # dpdM is the probability per unit mass
#                     # in this timestep : convert to dN
#                     dN = dpdM * dM

#                     # convert tests to 1/unit time
#                     if startype == "TEST_COUNT" or startype == "TEST_RATE":
#                         dN /= maxtime * 1e-6

#                     # hence the number of this type of star
#                     number[startype] += dN

#         # increase the time
#         logtime = float("%g" % (logtime + logtimestep))
#         # and loop

#     # ensemble
#     return ensemble
