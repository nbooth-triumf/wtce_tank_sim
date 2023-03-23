import numpy as np
import pickle


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

    def monte_carlo(self, norm_pdf, rand_var_min, rand_var_max):
        # Extract bin count, bin width, and normalized_max from pdf
        bin_count = len(norm_pdf)
        bin_width = (rand_var_max - rand_var_min) / bin_count
        norm_max = max(norm_pdf)
        f_big = 1.2*norm_max
        point_is_rejected = True

        # Generate rand_var from pdf and constraints
        while point_is_rejected:
            # Acceptance-Rejection method
            rand1 = np.random.uniform(0, 1)
            x_rand = rand_var_min + (rand_var_max - rand_var_min)*rand1
            rand2 = np.random.uniform(0, 1)
            f_rand = f_big*rand2

            # Determine bin and compare to pdf
            bin_rand = int(np.floor((x_rand - rand_var_min) / bin_width))
            f_x_rand = norm_pdf[bin_rand]
            if f_rand < f_x_rand:
                rand_var = x_rand
                point_is_rejected = False

        return rand_var

    # Save a projectiles object for later runs of code
    def save_projectiles(self, file_name):
        # open and close file around performing action
        with open(file_name, "wb") as f:
            pickle.dump(self, f)

    """"""

    # Open a projectiles object that was saved from a previous run
    @classmethod
    def open_projectiles(cls, file_name):
        # open and close file around performing action
        with open(file_name, "rb") as f:
            return pickle.load(f)

    """"""
