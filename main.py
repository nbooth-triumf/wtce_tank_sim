import numpy as np
import os
from tank_sim.Projectiles import Projectiles
from tank_sim.Histogram import Histogram
from tank_sim.Cylinder import Cylinder
from tank_sim.Graphics import Graphics

# Directions
make_projectiles = True
make_cylinder = True
open_object = ''
save_object = ''

# Define parameters of experiment being simulated (from 2021.3.31 WTCE Proposal pg 30 and others)
tank_height_m = 4.0         # metres
tank_diam_m = 4.1           # metres
pmt_eff_area = 65 * 10**-4  # metres^2
pmt_efficiency = 0.2        # assumed
dynamic_min = 0             # photons
dynamic_max = 50            # photons

# Define knowns (from data and experiments)
num_photons = 1 * 10**6     # one million photons as default
cos_th_min = 0.5            # dimensionless
cos_th_max = 1.0            # dimensionless
phi_min = 0                 # radians
phi_max = 2*np.pi           # radians

# Define cos(th) directory
data_group = "22_07_08"
pipe_label = "Pipe A"
working_folder = "C:/Users/booth/PycharmProjects/wtce_tank_sim/"
pickle_folder = working_folder + "pickle/"
path_header = working_folder + data_group + "/"
description_cos_th = "Cos(th) data from " + data_group
cos_th_folder = path_header + "data/cos_th_dist/"

# Initialize
description_projectile = "Direction Information for simulated Photons"
myPhotons = Projectiles(description_projectile, num_photons, cos_th_min, cos_th_max, phi_min, phi_max)
myTank = Cylinder(tank_height_m, tank_diam_m, pmt_eff_area)

""""""

if open_object == 'projectiles':
    # Open pickle file
    myPhotons = Projectiles.open_projectiles(pickle_folder + "projectiles.p")
    print("photon pickle opened successfully")

elif open_object == 'cylinder':
    # Open pickle file
    myTank = Cylinder.open_cylinder(pickle_folder + "cylinder.p")
    print("tank pickle opened successfully")
else:
    print("Will not open any pickle files.")
    # if not opening anything nor making anything
    if not make_projectiles and not make_cylinder:
        print("Directions incompatible with procedures. Please try different directions.")

""""""

if make_projectiles:
    # Combine all files in directory into full histogram
    histCosTh = Histogram(description_cos_th, cos_th_folder)
    histCosTh.combine_files()

    # Create normalized and cumulative pdfs
    histCosTh.create_normalized()
    histCosTh.create_cumulative()

    # Plot created distributions
    cos_th_range = np.linspace(cos_th_min, cos_th_max, histCosTh.num_bins, endpoint=False)
    cos_th_range = [round(value, 2) for value in cos_th_range]
    histCosTh.plot_dists(cos_th_range, 'cos(theta)', cos_th_folder + 'Histograms')

    """
    # Import csv of phi and create pdfs
    description_phi = "Phi dist for d=5"
    file_phi = "Data/dist_phi_d=5.csv"
    histPhi = Histogram(description_phi, file_phi)
    histPhi.create_normalize()
    histPhi.create_cumulative()
    """

    # Generate arbitrary number of photons using Monte Carlo method
    print("Time to generate some photons!")
    myPhotons.generate_direction_coords(pdf_cos_th=histCosTh.normalized)

    if save_object == 'projectiles':
        # Empty pickle folder to avoid memory issues
        for f in os.listdir(pickle_folder):
            os.remove(os.path.join(pickle_folder, f))

        # Save photons object
        myPhotons.save_projectiles(pickle_folder + "projectiles.p")
        print("photon pickle saved successfully")

""""""

if make_cylinder:
    # Place detector
    mpmt_height_m = tank_height_m * 0  # * np.random.uniform(-0.5, 0.5)

    # Determine impact coordinates for each photon trajectory
    myTank.generate_impact_coords(myPhotons.coords_list, mpmt_height_m)
    myTank.organize_data(stepsize=10**-2, pmt_eff=pmt_efficiency)

    if save_object == 'cylinder':
        # Empty pickle folder to avoid memory issues
        for f in os.listdir(pickle_folder):
            os.remove(os.path.join(pickle_folder, f))

        # Save pickle
        myTank.save_cylinder(pickle_folder + "cylinder.p")
        print("tank pickle saved successfully")

""""""

# Create graphics object and create figures
myGraphics = Graphics(myTank)
graphic_file = path_header + "tank_sim"
myGraphics.create_graphics(graphic_file, data_group, pipe_label)

# Create Intensity distribution
intensity_file = path_header + "tank_intensity"
myGraphics.create_intensity(intensity_file, data_group, pipe_label, log=False)
myGraphics.create_intensity(intensity_file, data_group, pipe_label, log=True)
