"""
Simple test to set up
"""

import agama
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

den = agama.Density(
    dict(
        type="Disk", surfaceDensity=8.95679e08, scaleRadius=2.5, scaleHeight=0.3
    ),  # thin disc from McMillan 2017
    dict(
        type="Disk", surfaceDensity=1.83444e08, scaleRadius=3.0, scaleHeight=0.9
    ),  # thick disc
    dict(
        type="Spheroid",
        densityNorm=9.8351e10,
        axisRatioZ=0.5,
        gamma=0,
        beta=1.8,
        scaleRadius=0.075,
        outerCutoffRadius=2.1,
    ),  # bulge (if you need it)
)
positions = den.sample(1000)[0]

# # Print first few positions
# print(stars[:5])


# Sun's Galactic coordinates
SUN_X = 8.2 * u.kpc  # kpc (distance from Galactic center)
SUN_Y = 0.0 * u.kpc  # kpc
SUN_Z = 0.0 * u.kpc  # kpc

# # Example Cartesian positions of stars (assumed to be in a Galactic frame)
# x = np.array([10, 12, 15])  # kpc (from Galactic Center)
# y = np.array([2, -5, 3])    # kpc
# z = np.array([0.5, -1, 2])  # kpc

x = positions[:, 0] * u.kpc
y = positions[:, 1] * u.kpc
z = positions[:, 2] * u.kpc

# Adjust x-coordinate to be Sun-relative
x_sun_rel = x - SUN_X

# Compute distance from the Sun
d = np.sqrt(x_sun_rel**2 + y**2 + z**2)  # Distance

# Compute Galactic longitude and latitude
l = np.arctan2(y, x_sun_rel)  # Galactic longitude
b = np.arcsin(z / d)  # Galactic latitude

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

# print("RA:", ra)
# print("Dec:", dec)

magnitude = 21 * np.ones(dec.shape)
# print('magnitude: ', magnitude)


# configure magnitude

# calculate detection probability

import os
from pathlib import Path

ding = (
    Path(os.getenv("GAIAUNLIMITED_DATADIR", "~/.gaiaunlimited")).expanduser().resolve()
)


from gaiaunlimited import DR3SelectionFunction  # Initialize Gaia DR3 Selection Function

sf = DR3SelectionFunction()

# Compute detection probabilities
p_detect = sf.query(coords=icrs_coords, gmag=magnitude)

print(np.min(p_detect))
print(np.max(p_detect))
