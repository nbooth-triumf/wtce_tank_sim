import csv
import numpy as np
from tank_sim.Projectiles import Projectiles
from tank_sim.Histogram import Histogram
from tank_sim.Cylinder import Cylinder

path_header = "C:/Users/booth/PycharmProjects/wtce_tank_sim/"

# Define knowns (from 2021.3.31 WTCE Proposal pg 30)
tank_height_m = 4.0      # metres
tank_diam_m = 4.1        # metres

# Import csv of cos(theta) and create pdfs
description_cos_th = "Cos(th) hist for d=5"
file_cos_th = "data/dist_cos_th_d=5.csv"
histCosTh = Histogram(description_cos_th, file_cos_th)
histCosTh.create_normalized()
histCosTh.create_cumulative()

"""
# Import csv of phi and create pdfs
description_phi = "Phi dist for d=5"
file_phi = "Data/dist_phi_d=5.csv"
histPhi = Histogram(description_phi, file_phi)
histPhi.create_normalize()
histPhi.create_cumulative()
"""

# Create Projectiles object of the extracted histograms
description_projectile = "Direction Information for simulated Photons"
num_photons = 1 * 10**4     # one million photons as default
myPhotons = Projectiles(description_projectile, num_photons)

# Generate arbitrary number of photons using Monte Carlo method
myPhotons.generate_direction_coords(pdf_cos_th=histCosTh.cumulative)
photon_trajectories = myPhotons.coords_list

# Create cylinder object of simulation tank and place an mPMT at a random height within it
myTank = Cylinder(tank_height_m, tank_diam_m)
mpmt_height_m = 0   # tank_height_m * np.random.uniform(-0.5, 0.5)

# Determine impact coordinates for each photon trajectory
myTank.generate_impact_coords(photon_trajectories, mpmt_height_m)

# Create graphic
graphic_file = path_header + "tank_sim"
myTank.create_graphics(graphic_file)

# Create Intensity distribution
intensity_file = path_header + "tank_intensity"
myTank.create_intensity(intensity_file)
