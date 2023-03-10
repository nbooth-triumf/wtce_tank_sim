import numpy as np
import csv
import os


class Projectiles(object):
    """
    Base class for list of projectiles and their direction coordinates
    """

    def __init__(self, description, num_projectiles, cos_th_min, cos_th_max, phi_min, phi_max):
        """
        Constructor
        """
        self.description = description      # Description of dataset
        self.num_proj = num_projectiles     # Number of projectiles to be given direction

        self.cos_th_min = cos_th_min        # Define cos_th range
        self.cos_th_max = cos_th_max
        self.phi_min = phi_min              # Define phi range
        self.phi_max = phi_max

        self.coords_list = []  # Initialize data structure of direction coordinates

    """"""

    def generate_direction_coords(self, pdf_cos_th, pdf_phi=[0]):
        count = 0  # track number of photon pairs

        while count < self.num_proj:
            # Determine random cos(theta)
            cos_th_rand = self.monte_carlo(pdf_cos_th,
                                           rand_var_min=self.cos_th_min, rand_var_max=self.cos_th_max)

            # Determine phi from uniform distribution unless otherwise specified
            if pdf_phi != [0]:
                phi_rand = self.monte_carlo(pdf_phi, rand_var_min=self.phi_min, rand_var_max=self.phi_max)
            else:
                phi_rand = np.random.uniform(self.phi_min, self.phi_max)

            # Combine and store
            direction_pair = [cos_th_rand, phi_rand]
            self.coords_list.append(direction_pair)
            count = count + 1

    """"""

    def monte_carlo(self, cumul_pdf, rand_var_min, rand_var_max):
        # Extract bin count and bin width from pdf
        bin_count = len(cumul_pdf)
        bin_width = (rand_var_max - rand_var_min) / bin_count

        # Generate random value between 0 and 1, uniformly
        value = np.random.uniform(0, 1.0)

        # Locate bin number where value is included in cumulative
        correct_bin = 0  # Initialize
        while correct_bin < bin_count:
            if value <= cumul_pdf[correct_bin]:
                # This is the correct bin
                break
            else:
                # Continue to next bin
                correct_bin = correct_bin + 1

        # Account for edge cases
        if correct_bin == 0:
            bin_to_left = 0
            counts_to_left = 0
        else:
            bin_to_left = correct_bin - 1
            counts_to_left = cumul_pdf[bin_to_left]

        # Determine exact random variable value relative to where value is between
        # cumulative probability density values
        cumulative_range = cumul_pdf[correct_bin] - counts_to_left
        relative_value = value - counts_to_left
        ratio = relative_value / cumulative_range
        bin_offset = ratio * bin_width
        rand_var = rand_var_min + bin_width * correct_bin + bin_offset

        return rand_var
