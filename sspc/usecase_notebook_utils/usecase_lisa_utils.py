"""
Copy of galactic stellar density distribution files written by Alexey for the LISA UCB project

Created on Wed Jun  5 18:05:24 2024

@author: alexey
"""

import sys
import time

import astropy.constants as const
import astropy.units as u
import numpy as np
import pandas as pd
import scipy as sp
from scipy.interpolate import CubicSpline, interp1d

# Units are kpc, Gyr


def zCDFInv(Xiz, Hz):  # DH0001
    zCoord = -Hz * np.log(1 - Xiz)
    return zCoord


def RCDFInv(Xir, Hr):  # DH0001
    # Get the parameters for the inverse CDF
    def RCD(R):  # DH0001
        Res = (1 - np.exp(-R / Hr)) - (R / Hr) * np.exp(-R / Hr) - Xir
        return Res

    Sol = sp.optimize.root_scalar(RCD, bracket=(0.0001 * Hr, 20 * Hr))
    if Sol.converged:
        R = Sol.root
    else:
        print("The radial solution did not converge")
        sys.exit()
    return R


def Sample1D(Hr, Hz):  # DH0001
    RRand = np.random.uniform()
    ZRand = np.random.uniform()
    ZSign = np.sign(np.random.uniform() - 0.5)
    R = RCDFInv(RRand, Hr)
    Z = zCDFInv(ZRand, Hz) * ZSign
    Th = 2.0 * np.pi * np.random.uniform()
    X = R * np.cos(Th)
    Y = R * np.sin(Th)
    Age = np.random.uniform(0, 12)

    XRel = X - 8.0
    YRel = Y
    ZRel = Z

    RRel = np.sqrt(XRel**2 + YRel**2 + ZRel**2)

    ResDict = {
        "Age": Age,
        "Xkpc": X,
        "Ykpc": Y,
        "Zkpc": Z,
        "Rkpc": R,
        "Th": Th,
        "XRelkpc": XRel,
        "YRelkpc": YRel,
        "ZRelkpc": ZRel,
        "RRelkpc": RRel,
    }

    return ResDict


def Sample1DPop(NBin, Hr, Hz):  # DH0001
    RRandSet = np.random.uniform(0, 1, NBin)
    ZRandSet = np.random.uniform(0, 1, NBin)
    ZSignSet = np.sign(np.random.uniform(0, 2, NBin) - 1)
    RSet = np.asarray([RCDFInv(RRandSet[i], Hr) for i in range(NBin)], dtype=float)
    ZSet = (
        np.asarray([zCDFInv(ZRandSet[i], Hz) for i in range(NBin)], dtype=float)
        * ZSignSet
    )
    ThSet = np.random.uniform(0, 2.0 * np.pi, NBin)
    XSet = RSet * np.cos(ThSet)
    YSet = RSet * np.sin(ThSet)
    AgeSet = np.random.uniform(0, 12, NBin)
    IDSet = np.arange(NBin) + 1

    ResDict = {
        "ID": IDSet,
        "Ages": AgeSet,
        "Xkpc": XSet,
        "Ykpc": YSet,
        "Zkpc": ZSet,
        "Rkpc": RSet,
        "Th": ThSet,
    }
    ResDF = pd.DataFrame(ResDict)
    return ResDF


ExportTable = False

if ExportTable:
    NBin = 10**4
    Hr = 4
    Hz = 0.5

    Res = Sample1DPop(NBin, Hr, Hz)

    Res.to_csv("./GalTest.csv", index=False)


def sample_distances_simple(NBin):  # DH0001
    """
    Simple distance sampler using the
    """

    Hr = 4
    Hz = 0.5

    Res = Sample1DPop(NBin, Hr, Hz)

    del Res["Ages"]
    del Res["ID"]
    del Res["Th"]
    del Res["Rkpc"]

    galcen_distance = 8.122  # * u*kpc

    Res["Xkpc_rel"] = Res["Xkpc"] - galcen_distance
    Res["Ykpc_rel"] = Res["Ykpc"]
    Res["Zkpc_rel"] = Res["Zkpc"]

    Res["Distance_to_sun"] = np.sqrt(
        (Res["Xkpc_rel"] ** 2) + (Res["Ykpc_rel"] ** 2) + (Res["Zkpc_rel"])
    )

    return Res["Distance_to_sun"].to_numpy() * u.kpc


###########################
#
###########################


def get_bin_frac_ratio(IC_model, binary_fraction=0.5):  # DH0001

    # these are the hard coded ratios based on initial conditions sampling
    # tests done by K. Breivik using binary fractions from 0.1-1.0
    ratio_dict = {
        "m2_min_05": [
            5.497065159607499,
            4.028385607219321,
            3.1320443249502383,
            2.541669993010094,
            2.105948317768266,
            1.786314805728416,
            1.5164506548687628,
            1.313816074228385,
            1.141178996587137,
            0.9968331286953663,
            0.8730241509546721,
            0.7696547388705937,
            0.6788737987994921,
            0.5998311488995663,
            0.5308139080430776,
            0.4699942084336498,
            0.41210469109378484,
            0.36061186819185237,
            0.31590620434046357,
            0.27425387234835646,
            0.23543985197379194,
            0.19987998623485503,
            0.16926018076106297,
            0.13969418254623414,
            0.11143821417288673,
            0.08563401307212755,
            0.06239655229658728,
            0.03973086411370632,
            0.019317445609880482,
            0.0,
        ],
        "log_uniform_porb": [
            5.886862859310299,
            4.3420507204072525,
            3.405093627693288,
            2.7466102648326007,
            2.27052928301355,
            1.9273925451530112,
            1.6405699036728623,
            1.4078089101945124,
            1.2331941437437164,
            1.069552690381338,
            0.9456643329362575,
            0.8327939389648276,
            0.7359840090556246,
            0.6480801205088771,
            0.5722387172813348,
            0.5039632125801372,
            0.4440067853899282,
            0.3913613364579611,
            0.3402636518883438,
            0.29843507961701765,
            0.2543196815715736,
            0.21741447251907284,
            0.18304208420964435,
            0.15107685764473266,
            0.12041807829958329,
            0.09269327057199453,
            0.06744511227698262,
            0.04274378038546733,
            0.02137711323393566,
            0.0,
        ],
        "ecc_uniform": [
            5.889211849314321,
            4.361889407239361,
            3.4020828542097057,
            2.746002579297502,
            2.2673692464054183,
            1.9141940732871088,
            1.6405959517007267,
            1.4094329827706373,
            1.230470758764182,
            1.072042568070819,
            0.943431196191181,
            0.8333833246264916,
            0.7309333535548196,
            0.6475817013002746,
            0.574417702589979,
            0.5056945660498336,
            0.44488117789149545,
            0.3884654642073608,
            0.3388900220256573,
            0.2962626999212119,
            0.2557025971124346,
            0.21685490011138606,
            0.18305848181583317,
            0.15006625511716296,
            0.12090585737779058,
            0.09221860126908878,
            0.06728132687768605,
            0.043606002397485674,
            0.021239106748386388,
            0.0,
        ],
        "qmin_01": [
            5.823149995795932,
            4.257871067861341,
            3.3350623602831586,
            2.701648948025014,
            2.241763868198779,
            1.8890779951715042,
            1.6099148365641083,
            1.3953172844589194,
            1.2110201528032853,
            1.0592751609913709,
            0.9225325052920952,
            0.8191978467059325,
            0.7180939238557436,
            0.6326983982513528,
            0.5640129919280746,
            0.49581103805511545,
            0.43644888054132813,
            0.3834832365660745,
            0.3356910316544453,
            0.29248525094604144,
            0.2511267079306599,
            0.21240840773299557,
            0.1776591271238454,
            0.14765505334148785,
            0.11783814018942763,
            0.09097678374483624,
            0.06650626075146698,
            0.04210676407248778,
            0.0203586397778618,
            0.0,
        ],
        "fiducial": [
            5.932296012899532,
            4.357823767133379,
            3.394145723221763,
            2.743531290685849,
            2.261347038645783,
            1.9173522895783706,
            1.6465752814658685,
            1.409849902360307,
            1.232375306106761,
            1.0745464651773011,
            0.9431907171044254,
            0.8331178361129032,
            0.7317070852105148,
            0.6506914750851542,
            0.5707015794152772,
            0.5048166662383198,
            0.44379738107168043,
            0.3913719830983486,
            0.33791478738462305,
            0.29476240868778614,
            0.2551125640746029,
            0.21752011832877077,
            0.18147403979103155,
            0.14971671211685017,
            0.1209730374015313,
            0.09372046390290077,
            0.06743095018359142,
            0.043359831467175966,
            0.02054866077700454,
            0.0,
        ],
        "ecc_thermal": [
            5.9200035744176125,
            4.341656072119606,
            3.4004867453530436,
            2.747854677340573,
            2.2726296004773694,
            1.9140533603622882,
            1.6388693872876943,
            1.4154780593615193,
            1.2294748612513091,
            1.0760766054593665,
            0.9442728392343459,
            0.8321788961090938,
            0.7330076399499459,
            0.6473674909796342,
            0.5733259021802497,
            0.5031744591982074,
            0.44384269602503534,
            0.38875209683437956,
            0.3429980155840034,
            0.2955812055263655,
            0.25491188680860966,
            0.21846399667917327,
            0.18130416084339843,
            0.14918121603752904,
            0.12215748474240215,
            0.09270863580218562,
            0.06716028032728252,
            0.04318939158189473,
            0.020791845074920087,
            0.0,
        ],
    }
    binfracs = np.linspace(0.1, 1.0, 30)

    # select the list of ratios based on the initial conditions model
    ratio = ratio_dict[IC_model]

    # set up a spline to get the ratio for any binfrac
    r_spline = CubicSpline(binfracs, ratio)

    return r_spline(binary_fraction)


def get_mass_norm(IC_model, binary_fraction=0.5):  # DH0001
    """selects the mass normalization for the
    initial conditions sample set based on
    the IC_model name and a binary fraction

    Parameters
    ----------
    IC_model : `str`
        initial conditions model chosen from:
            ecc_uniform, ecc_thermal, porb_log_uniform, m2_min_05, qmin_01, fiducial

    Returns
    -------
    mass_norm : `float`
        the total ZAMS mass of the initial stellar population
        including single and binary stars
    """

    mass_binaries = {
        "ecc_uniform": 2720671.1164002735,
        "ecc_thermal": 2700943.07050043,
        "porb_log_uniform": 2713046.6197530716,
        "m2_min_05": 2905718.830512573,
        "qmin_01": 5510313.245766795,
        "fiducial": 2697557.2681495477,
    }

    # get the ratio of singles to binaries for the selected binary fraction
    ratio = get_bin_frac_ratio(IC_model, binary_fraction=binary_fraction)
    mass_total = mass_binaries[IC_model] * (1 + ratio)

    return mass_total


######################
#
######################


def get_period(semimajor_axis, m1, m2):  # DH0001
    """
    function to get the periods of the systems
    """

    mtot = m1 + m2
    p2 = ((semimajor_axis**3) * 4 * np.pi**2) / (const.G * mtot)
    p = np.sqrt(p2)

    return p.to(u.yr)


def precompute_radial_cdf(Hr, num_points=1000, R_max=20):  # DH0001
    """
    Precompute the radial CDF and its inverse for efficient sampling.
    """
    # Generate a fine grid of R values
    R_values = np.linspace(0, R_max * Hr, num_points)
    CDF_values = (1 - np.exp(-R_values / Hr)) - (R_values / Hr) * np.exp(-R_values / Hr)

    # Ensure the CDF is monotonic
    CDF_values[0] = 0.0  # Explicitly set the minimum CDF value
    CDF_values[-1] = 1.0  # Explicitly set the maximum CDF value (asymptotic)

    # Create an interpolation function for the inverse CDF
    inverse_cdf = interp1d(
        CDF_values, R_values, bounds_error=False, fill_value="extrapolate"
    )
    return inverse_cdf


def RCDFInv_interpolated(Xir, inverse_cdf):  # DH0001
    """
    Use the precomputed inverse CDF to sample R values.
    """
    return inverse_cdf(Xir)


def Sample1DPop_interpolated(NBin, Hr, Hz, inverse_cdf):  # DH0001
    """
    Sample 1D population using precomputed inverse CDF.
    """
    RRandSet = np.random.uniform(0, 1, NBin)
    ZRandSet = np.random.uniform(0, 1, NBin)
    ZSignSet = np.random.choice([-1, 1], NBin)

    # Use interpolation for RCDFInv
    RSet = RCDFInv_interpolated(RRandSet, inverse_cdf)
    ZSet = zCDFInv(ZRandSet, Hz) * ZSignSet
    ThSet = np.random.uniform(0, 2.0 * np.pi, NBin)

    # Compute Cartesian coordinates
    XSet = RSet * np.cos(ThSet)
    YSet = RSet * np.sin(ThSet)

    return np.stack((XSet, YSet, ZSet), axis=1)


def sample_distances_interpolated(NBin, Hr, inverse_cdf):  # DH0001
    """
    Main function to sample distances using interpolation for the radial CDF.
    """

    Hz = 0.5  # Vertical scale height
    galcen_distance = 8.122  # kpc

    # Sample positions
    positions = Sample1DPop_interpolated(NBin, Hr, Hz, inverse_cdf)

    # Compute relative coordinates and distances
    XSet_rel = positions[:, 0] - galcen_distance
    YSet_rel = positions[:, 1]
    ZSet_rel = positions[:, 2]

    distances = np.sqrt(XSet_rel**2 + YSet_rel**2 + ZSet_rel**2)

    return distances * u.kpc

def sample_distances_and_angles_interpolated(NBin, Hr, inverse_cdf):  # DH0001
    """
    Main function to sample distances using interpolation for the radial CDF.
    """

    Hz = 0.5  # Vertical scale height
    galcen_distance = 8.122  # kpc

    # Sample positions
    positions = Sample1DPop_interpolated(NBin, Hr, Hz, inverse_cdf)

    # Compute relative coordinates and distances
    XSet_rel = positions[:, 0] - galcen_distance
    YSet_rel = positions[:, 1]
    ZSet_rel = positions[:, 2]

    # Compute distance from the Sun
    d = np.sqrt(XSet_rel**2 + YSet_rel**2 + ZSet_rel**2) * u.kpc  # Distance

    # Compute Galactic longitude and latitude
    l = np.arctan2(y, x_sun_rel) * u.rad  # Galactic longitude
    b = np.arcsin(z / d) * u.rad  # Galactic latitude

    # Convert to degrees
    l = l.to(u.deg)
    b = b.to(u.deg)

    # Create SkyCoord in Galactic frame
    galactic_coords = SkyCoord(l=l, b=b, distance=d, frame="galactic")

    # Convert to ICRS (RA/Dec)
    icrs_coords = galactic_coords.transform_to("icrs")

    # Extract RA and Dec
    ra = icrs_coords.ra.deg
    dec = icrs_coords.dec.deg

    return d, ra, dec



if __name__ == "__main__":
    import functools

    Hr = 4  # Radial scale length

    # Precompute the inverse CDF for radial distribution
    inverse_cdf = precompute_radial_cdf(Hr)

    bound_sample_distances_interpolated = functools.partial(
        sample_distances_interpolated, inverse_cdf=inverse_cdf
    )

    num_system_indices = 100
    convolution_results = {}
    num_repeats = 5

    #
    range_sample_sizes = [
        10,
        100,
        1000,
        10000,
        100000,
        1000000,
        10000000,
        # 100000000,
    ]

    #
    for sample_size in range_sample_sizes:
        # sample distances and add those to the result dict
        start_dist = time.time()
        for _ in range(num_repeats):

            dist = bound_sample_distances_interpolated(NBin=sample_size)
            convolution_results["dists"] = dist
        stop_dist = time.time()

        print(
            "Distance-sampling {} systems on average took {}s".format(
                sample_size, (stop_dist - start_dist) / num_repeats
            )
        )
