# """ DH0001_file

# Script to make population statistics and their
# history for a given star-formation rate, e.g.
# for the Milky Way for which a SFR(t) is given.

# Args:
#     1) the ensemble data file
#     2) output file (if not given, uses stdout)
# """

# # import numba; from numba import njit
# import datetime
# import math
# import sys
# import time

# import numexpr as ne
# import numpy as np
# import py_rinterpolate
# import pyfma
# import scipy
# from binarycpython.utils.dicts import keys_to_floats
# from binarycpython.utils.ensemble import load_ensemble

# ############################################################

# starttime = 0.0
# dt = 1e6  # linear timestep
# present_day = 13.7e9  # Galactic age, now (yr)
# maxtime = 0.1 * present_day  # present_day # age of the Milky Way (yr)
# ntimebins = int(1 + maxtime / dt)
# event_types = ("SN", "MERGER", "MERGER_INTO")
# star_types = ("SPECTRAL_TYPE", "NUMBER_OF", "TEST", "MASS_IN")
# observations = {
#     "MASS_IN_STARS": "4.6e10 to 6.43e10 Msun",
#     "NUMBER_OF_STARS": "1e11 to 4e11",
# }


# ############################################################
# # functions
# #
# def obs(x):
#     if x in observations:
#         return f"Observed: {observations[x]:25s}"
#     else:
#         return ""


# def SFR(time):
#     ############################################################
#     # Star-formation rate as a function of time (years)
#     # since the birth of the milky way, based on Gaia DR2
#     # (Mor et al. 2019)
#     #
#     # Input: time since birth of the Galaxy / years
#     ############################################################

#     # time now
#     time_Gyr = time * 1e-9

#     # hence age in Gyr
#     age_Gyr = (present_day - time) * 1e-9

#     # data only goes back 10 Gyr
#     age_Gyr = min(10.0, age_Gyr)

#     # rough fit to Mor et al. (2019) : Msun/Gyr/pc^2
#     sfr = (0.7946) * math.exp((0.2641) * age_Gyr) + (7.3643) * math.exp(
#         -((age_Gyr - (2.5566)) ** 2) / (3.036)
#     )

#     # convert to Msun/year assuming (as Mor+ do)
#     # 1Msun/year now == 1.58Msun/Gyr/pc^2 (data point for now)
#     sfr *= 1.0 / 1.58

#     # no young stars
#     # print(f"Time {time:g} -> time_Gyr {time_Gyr:g} -> age_Gyr {age_Gyr:g} -> sfr {sfr}")

#     return sfr


# def times(ensemble):
#     # return a sorted list of times in the ensemble
#     times = {}
#     for scalar in ensemble["scalars"]:
#         for logt in ensemble["scalars"][scalar]:
#             times[float(logt)] = 1
#     return sorted(times.keys())


# def ensemble_to_image(ensemble):
#     ############################################################
#     # convert an ensemble's scalars to a 2D array
#     # in which column 0 is the time
#     #
#     # returns the converted image, and a tuple
#     # of the number of columns and rows
#     ############################################################
#     _times = times(ensemble)
#     scalarkeys = list(ensemble["scalars"].keys())
#     ncols = 1 + len(scalarkeys)  # +1 for log time
#     nrows = len(_times)
#     N = np.zeros((nrows, ncols))
#     nrows = N.shape[0]
#     ncols = N.shape[1]
#     for nrow, _time in enumerate(_times):
#         # column time_column is the time
#         N[nrow, 0] = _time
#         # all others are data
#         for ncol, scalar in enumerate(scalarkeys, start=1):
#             if ncol < ncols:
#                 N[nrow, ncol] = ensemble["scalars"][scalar].get(_time, 0.0)

#     return N, (ncols, nrows)


# def timestamp():
#     presentDate = datetime.datetime.now()
#     return datetime.datetime.timestamp(presentDate)


# def logimage_to_linearimage(log_image):
#     # convert a log-time ensemble image to a linear-time ensemble image
#     # (without any interpolation)
#     linear_image = np.copy(log_image)
#     # times : log -> linear
#     linear_image[:, 0] = ne.evaluate("10.0**x", local_dict={"x": linear_image[:, 0]})
#     # other values are per unit linear time, so can remain the same
#     return linear_image


# def regularize_image(irregular_image, wanted_times):
#     ############################################################
#     # Convert an image at non-fixed times to
#     # an image at the times we want, i.e.
#     # at "regular" times.
#     #
#     # Input:
#     # 1) an irregularly-timed image
#     # 2) a list of wanted times
#     #
#     # Returns the regularized 2D "image" [[],[],...]
#     #
#     ############################################################
#     nrows = irregular_image.shape[0]
#     ncols = irregular_image.shape[1]
#     n_wanted_times = len(wanted_times)

#     rinterpolator = py_rinterpolate.Rinterpolate(
#         table=np.ndarray.flatten(irregular_image, order="C"),
#         nparams=1,  # time is only parameter
#         ndata=ncols - 1,  # all other columns are data
#         nlines=nrows,  # manually set this!
#         verbosity=0,
#     )
#     regular_image = np.zeros((n_wanted_times, ncols))

#     for ntime, time in enumerate(wanted_times):
#         result = rinterpolator.interpolate([time])
#         regular_image[ntime, 0] = time
#         regular_image[ntime, 1:] = result

#     del rinterpolator  # garbage cleanup
#     return regular_image.reshape((n_wanted_times, ncols))


# # Set up file paths
# ensemble_file = str(sys.argv[1])
# outfile = str(sys.argv[2]) if len(sys.argv) > 2 else "/dev/stdout"
# print(f"Ensemble data from {ensemble_file}")
# print(f"Output results to {outfile}")

# # load the ensemble and check it has data
# data = load_ensemble(ensemble_file, allow_nan=True, filter_nones=True)
# ensemble = data["ensemble"]
# metadata = data["metadata"]
# bse_options = metadata["settings"]["population_settings"]["bse_options"]

# if False:
#     # dummy data
#     bse_options["ensemble_logtimes"] = 0
#     ensemble["scalars"] = {
#         "TEST": {
#             "0.0": 0.0,
#             "1.0": 2.0,
#             "2.0": 2.0,
#             "3.0": 0.0,
#             "4.0": 0.0,
#             "5.0": 0.0,
#             "6.0": 0.0,
#             "7.0": 0.0,
#             "8.0": 0.0,
#             "9.0": 0.0,
#             "10.0": 0.0,
#         }
#     }

# # convert keys to floats (if possible)
# ensemble = keys_to_floats(ensemble)

# # load logtime ensemble into a 2D image
# N, (ncols, nrows) = ensemble_to_image(ensemble)

# # convert to linear times?
# if bse_options.get("ensemble_logtimes", 0):
#     N = logimage_to_linearimage(N)

# # convert time from Myr to yr
# N[:, 0] *= 1e6

# # make a list of times we want to output
# wanted_times = np.arange(0.0, maxtime + dt, dt)

# # regularize the image to the wanted times
# N = regularize_image(N, wanted_times)

# # number of header and data items
# header = ("Time", "dt", "dM")  # header columns (for logging)
# nheader = len(header)  # used often
# ndata = ncols - 1  # number of data columns

# # compute future mass ejected as a function of time in an array
# dM = np.zeros((ntimebins + 1))
# for i, tnow in enumerate(wanted_times):
#     dM[i] = dt * SFR(tnow)

# # make Galactic history
# Galactic_history = np.zeros((ntimebins + 1, nheader + ndata))
# j = len(wanted_times)

# ## try jit?
# # @njit(parallel=True)
# # def jit_fma(a,b,c):
# #    return a + b * c

# tcpu0 = timestamp()
# dfdt0 = None
# fprev = 0.0
# for i, tnow in enumerate(wanted_times):
#     # logging
#     tcpu = timestamp()
#     f = float(i) / float(len(wanted_times))
#     if tcpu > tcpu0 + 1.0:
#         if dfdt0 == None:
#             dfdt0 = f / (tcpu - tcpu0)
#         elif f > fprev + 1e-2:
#             n = 1.05
#             fcomplete = 1.0 - (1.0 - f) ** (n + 1)
#             print(f"t={tnow: 15g} y : {100.0*fcomplete:5.2f}%", end="\r")
#             fprev = f

#     # mass added in stars in this timestep
#     _dM = dM[i]
#     # _dM = 1.0 if i == 0 else 0.0 # testing

#     # first three cols are: time, dt, dM (for output, e.g. plotting)
#     Galactic_history[i, 0:2] = tnow
#     Galactic_history[i, 1] = dt
#     Galactic_history[i, 2] = _dM

#     # rest are counts/rates: convolve these
#     if _dM > 0.0:
#         # convolve (this is the slow step!)

#         # 9.30user 1.15system 0:09.25elapsed 112%CPU
#         # Galactic_history[i:i+j,nheader:] += _dM * N[0:j, 1:]

#         # 8.77user 1.11system 0:08.67elapsed 113%CPU
#         # np.add(Galactic_history[i:i+j,nheader:],
#         #       _dM * N[0:j, 1:],
#         #       out=Galactic_history[i:i+j,nheader:])

#         # single evaulate call:
#         # 13.90user 1.81system 0:05.62elapsed 279%CPU
#         #
#         # with re-evaluate:
#         # 13.00user 1.17system 0:04.83elapsed 293%CPU
#         #
#         # definitely the best multi-core performance
#         # if more overall CPU time and doesn't use all
#         # the cores all the time
#         #
#         ldict = {"x": Galactic_history[i : i + j, nheader:], "y": N[0:j, 1:]}
#         # 46.07user 1.82system 0:13.89elapsed 344%CPU
#         Galactic_history[i : i + j, nheader:] = (
#             ne.evaluate("x + _dM * y", local_dict=ldict)
#             if i == 0
#             else ne.re_evaluate(local_dict=ldict)
#         )

#         # jit uses more CPU time than numexp but overall runs
#         # more slowly! (overheads?)
#         # 106.84user 1.36system 0:16.09elapsed 672%CPU
#         # Galactic_history[i:i+j,nheader:] = jit_fma(Galactic_history[i:i+j,nheader:],
#         #                                           _dM,
#         #                                           N[0:j, 1:])

#         # (much) slower! also drops the final line (?)
#         # 17.11user 1.12system 0:17.03elapsed 107%CPU
#         # Galactic_history[i:i+j,nheader:] = pyfma.fma(_dM, N[0:j, 1:],  Galactic_history[i+j,nheader:])
#     j -= 1

# # remove final row
# Galactic_history = Galactic_history[:-1]

# # output
# if outfile == "/dev/stdout":
#     prestring = "OUTPUT "
#     f = sys.stdout
# else:
#     prestring = ""
#     f = open(outfile, "w")
# f.write(
#     f"{prestring}{' '.join(header)} {' '.join(list(ensemble['scalars'].keys()))}\n{prestring}"
# )
# np.savetxt(f, Galactic_history, fmt="%g", newline=f"\n{prestring}")
# f.write("\n")
# f.close()
# sys.exit()

# ############################################################
# # TODO!
# ############################################################

# if False:
#     if not "scalars" in ensemble:
#         print("No scalars found in the ensemble: have your stars run for long enough?")
#         exit()

#     # SNIa_types = ('SNIA_He','SNIA_DD','SNIA_CHAND','SNIA_He_Coal','SNIA_CHAND_Coal','SNIA_Hybrid_HeCOWD','SNIA_Hybrid_HeCOWD_subluminous')

#     eventtypes = [s for s in ensemble["scalars"] if s.startswith(event_types)]

#     startypes = [s for s in ensemble["scalars"] if s.startswith(star_types)]

#     print(eventtypes)

#     logtimestep = bse_options["ensemble_logdt"]
#     minlogtime = bse_options["ensemble_startlogtime"]
#     maxlogtime = -1e30
#     for event in eventtypes:
#         for time in ensemble["scalars"][event]:
#             maxlogtime = max(maxlogtime, time)

#     print(f"time {minlogtime} to {maxlogtime} logtimestep {logtimestep}")

#     rate = AutoVivificationDict()
#     number = AutoVivificationDict()
#     Mtot_formed = 0
#     dtyearstot = 10.0 ** (minlogtime + 6.0 - logtimestep)

#     # loop over time
#     logtime = minlogtime
#     while logtime <= maxlogtime:
#         ############################################################
#         # time bins are from tlow=logt-logtimestep to thigh=logt
#         ############################################################
#         logt = logtime + 6.0  # log time in years
#         t = 10.0**logt  # linear time in yr
#         tlow = 10.0 ** (logt - logtimestep)
#         thigh = t
#         thalf = 0.5 * (tlow + thigh)
#         dtyears = thigh - tlow

#         if tlow > maxtime:
#             # none of the bin is before the present day
#             ftime = 0.0
#         elif thigh > maxtime:
#             # some of the bin is before the present day
#             ftime = (maxtime - tlow) / dtyears
#         else:
#             # all of the bin is before the present day
#             ftime = 1.0
#         dtyears *= ftime

#         # total years simulated
#         dtyearstot += dtyears

#         if ftime > 0.0 and ftime <= 1.0:
#             print(
#                 f"At log(time/Myr) {logtime:12g} log(time/yr) {logt:12g} time {t:12g}, time bin from {tlow:12g} to {thigh:12g}, dtyears = {dtyears:12g}, tot {dtyearstot:12g} {dtyearstot/maxtime:12g}, ftime = {ftime:g}"
#             )

#             # in this dtyears we formed dM (Msun) mass into stars
#             dM = dtyears * SFR(thalf)
#             Mtot_formed += dM

#             # rates
#             for event in eventtypes:
#                 if logtime in ensemble["scalars"][event]:
#                     # 1e-6 converts from per Myr to per year
#                     drate = 1e-6 * ensemble["scalars"][event][logtime]

#                     rate[event] += drate * dM

#             # numbers of stars
#             for startype in startypes:
#                 if logtime in ensemble["scalars"][startype]:
#                     # binary_c's ensemble gives us the dprob/dM i.e.
#                     # whatever per unit mass into stars in each bin.
#                     #
#                     # There is no time involved here,
#                     # so we don't need to convert from Myr to yr.
#                     dpdM = ensemble["scalars"][startype][logtime]

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

#     print(f"Simulated {dtyearstot:g} years of evolution (should be {maxtime:g})\n")

#     if "qcrit_nuclear_burning" in bse_options:
#         print(f"\nQcrit algorithm : {bse_options['qcrit_nuclear_burning']}\n")

#     for opt in bse_options:
#         print(f"{opt} = {bse_options[opt]}")
#     print()

#     # rates
#     print("events")
#     for prefix in event_types:
#         total = 0
#         subtypes = sorted(
#             [
#                 s
#                 for s in rate
#                 if s.startswith(prefix)
#                 and not s.endswith(
#                     ("LUMINOSITY", "POWER", "DARK", "BRIGHT", "MOMENTUM")
#                 )
#             ]
#         )
#         for event_type in subtypes:
#             if prefix == "MERGER_INTO" or not "_INTO_" in event_type:
#                 total += rate[event_type]

#         for event_type in subtypes:
#             if prefix == "MERGER_INTO" or not "_INTO_" in event_type:
#                 if rate[event_type] == 0.0:
#                     print(
#                         f"{event_type:30s} {rate[event_type]: 15g} = {100*rate[event_type]/total: 8.2f}% i.e. never {obs(event_type)}"
#                     )
#                 else:
#                     print(
#                         f"{event_type:30s} {rate[event_type]: 15g} = {100*rate[event_type]/total: 6.2f}% i.e. every {1.0/rate[event_type]: 15g} years or {rate[event_type]/Mtot_formed: 15g} per Msun per year {obs(event_type)}"
#                     )

#         print()

#     # number of stars
#     print("stars")
#     print(star_types)
#     print(startypes)
#     for prefix in star_types:
#         total = 0
#         subtypes = sorted([s for s in number if s.startswith(prefix)])
#         for startype in subtypes:
#             total += number[startype]

#         for startype in subtypes:
#             if number[startype] == 0.0:
#                 print(
#                     f"{startype:30s} {number[startype]: 15g} = {100*number[startype]/total: 8.2f}% i.e. never {obs(startype)}"
#                 )
#             else:
#                 print(f"{startype:30s} {number[startype]: 15g} {obs(startype)}")

#         print()

#     print(f"Formed {Mtot_formed:g} Msun in stars in {dtyearstot:g} yr of simulation")
#     print(
#         f"average stellar mass {Mtot_formed/number['NUMBER_OF_STARS']:g} should be ~ 0.5Msun"
#     )
