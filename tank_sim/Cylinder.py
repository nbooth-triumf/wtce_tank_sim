import numpy as np
import matplotlib.pyplot as plt


class Cylinder(object):
    """
    Base class for defining tank region and corresponding geometry
    """

    def __init__(self, height_m, diameter_m):
        """
        Constructor - HEIGHT AND DIAMETER IN METRES
        """
        self.height = height_m * 1000       # Tank height [=] mm
        self.z_limit = self.height / 2      # Assuming origin in centre of tank, min/max value of z
        self.diameter = diameter_m * 1000   # Tank diameter [=] mm
        self.radius = self.diameter / 2     # Tank radius [=] mm

        self.num_photons = None             # Initialize number of projectiles
        self.wall_impact_coords = []        # Initialize wall impact coords of [gamma, impact_height]
        self.lid_impact_coords = []         # Initialize lid_impact_coords of [r, alpha]
        self.base_impact_coords = []        # Initialize base_impact_coords of [r, alpha]

        self.lid_radii = []                 # Initialize relevant data structures for plotting
        self.lid_alpha = []                 # on individual axes for lid, wall, and base
        self.wall_alpha = []
        self.wall_height = []
        self.base_radii = []
        self.base_alpha = []

        self.lid_r = None
        self.lid_a = None
        self.wall_a = None
        self.wall_z = None
        self.base_r = None
        self.base_a = None

        self.lid_counts = None
        self.wall_counts = None
        self.base_counts = None
        self.counts_max = None

        self.troubleshooting_info = []      # Data structure to save relevant values if error occurs

    """"""

    def generate_impact_coords(self, directional_coords, source_height_m):
        # convert to mm
        source_height = source_height_m * 1000      # [=] mm
        self.num_photons = len(directional_coords)

        for n in range(self.num_photons):
            # Extract data
            [cos_th, phi] = directional_coords[n]
            theta_init = np.arccos(cos_th)

            # Account for air / plastic / water interfaces via Snell's Law
            # Refractive index of plastic is not needed due to algebra
            n_air = 1.0003
            n_water = 1.333
            n_ratio = n_air / n_water
            theta = np.arcsin(n_ratio * np.sin(theta_init))

            # Define unit vectors in the source (prime) reference frame from given theta and phi
            x_prime_hat = np.sin(theta)*np.sin(phi)
            y_prime_hat = np.sin(theta)*np.cos(phi)
            z_prime_hat = np.cos(theta)

            # Rotate and shift unit vector from source frame to tank frame
            x_hat = y_prime_hat
            y_hat = z_prime_hat
            z_hat = x_prime_hat

            # Project lightray onto the XY plane and determine polynomial y = m*x + b
            # Determine intersection of polynomial and circle x**2 + y*2 = R**2, where R is the known radius

            # Account for limits
            if theta == 0 or np.abs(phi) == np.pi/2:
                x_impact = 0
                y_impact = self.radius
            # General
            else:
                m_xy = y_hat / x_hat
                x_impact = 2 * m_xy * self.radius / (m_xy**2 + 1)
                y_impact = m_xy * x_impact - self.radius

            # Troubleshooting - comment out by default
            """
            r_math = np.sqrt(x_impact**2 + y_impact**2)
            r_actual = self.radius
            if np.abs(r_math - r_actual) > 0.001:
                print("impact coordinate match is incorrect")
                break
            """

            # Use similar triangles and unit vectors to obtain z
            delta_z_y = (y_impact + self.radius) * (z_hat / y_hat)
            # Account for limits in x_hat (y_hat is never zero)
            if x_hat == 0.0:
                delta_z_x = 0.0
            else:
                delta_z_x = x_impact * (z_hat / x_hat)

            if np.abs(delta_z_y - delta_z_x) > 0.001:
                # If these do not match, the geometry/vectors are incorrect
                print("impact coordinate math is incorrect")
                break
            else:
                # If they do match, save as the correct value
                delta_z = delta_z_x
            z_impact = source_height + delta_z

            # Determine if height of impact point is within tank
            if np.abs(z_impact) > self.z_limit:
                # Impact point is outside tank.
                # Use similar triangles to determine x_cap and y_cap at z = z_limit
                x_cap = x_impact * self.z_limit / np.abs(z_impact)
                y_cap = (y_impact + self.radius) * self.z_limit / np.abs(z_impact) - self.radius

                # Convert [x_cap, y_cap, z_limit] to [r, alpha, z_limit]
                r = np.sqrt(x_cap**2 + y_cap**2)
                alpha = np.arctan2(y_cap, x_cap) - np.pi/2  # Rotate from x-axis centric to y-axis centric

                # Store impact coords in correct data structure
                if z_impact > 0:
                    self.lid_impact_coords.append([r, alpha, self.z_limit])
                else:
                    self.base_impact_coords.append([r, alpha, self.z_limit])

            else:
                # Impact point is inside tank
                # Convert [x, y, z] to [r, alpha, z] and confirm r = self.radius
                r = np.sqrt(x_impact**2 + y_impact**2)
                alpha = np.arctan2(x_impact, y_impact)

                # Sanity check
                if np.abs(r - self.radius) < 0.001:     # [=] mm
                    self.wall_impact_coords.append([r, alpha, z_impact])
                else:
                    print("Error encountered. Wall impact not at wall.")
                    self.troubleshooting_info.append([cos_th, phi, x_impact, y_impact, z_impact])

    """"""

    def create_graphics(self, file_name, show=False):
        # Prepare data
        self.isolate_axes()     # Extract each individual axes as a separate dataset
        self.make_meshes()      # Create 2D arrays in polar, cart, polar (for lid, wall, base, respectively)

        # Plot
        scatter_name = file_name + "_scatter.jpg"
        self.make_scatter_plot(scatter_name, show)
        heatmap_name = file_name + "_heatmap.jpg"
        self.make_heatmap(heatmap_name, show)
        """
        plt.imshow(self.wall_counts)
        ax = plt.gca()
        ax.set_ylim(ax.get_ylim()[::-1])
        plt.colorbar()
        plt.savefig(heatmap_name)
        """
    """"""

    def isolate_axes(self):
        # Put data into plottable form by separating each axes into a distinct data structure
        # Lid
        for n in range(len(self.lid_impact_coords)):
            # z is constant, extract r and alpha coords
            [radius, alpha, z] = self.lid_impact_coords[n]
            self.lid_radii.append(radius)
            self.lid_alpha.append(alpha)

        # Wall
        for n in range(len(self.wall_impact_coords)):
            # r is constant, extract alpha and z coords
            [radius, alpha, z] = self.wall_impact_coords[n]
            self.wall_alpha.append(alpha)
            self.wall_height.append(z)

        # Base
        for n in range(len(self.base_impact_coords)):
            # z is constant, extract r and alpha coords
            [radius, alpha, z] = self.base_impact_coords[n]
            self.base_radii.append(radius)
            self.base_alpha.append(alpha)

    """"""

    def make_meshes(self, num_elements=100):
        # Initialize span of each dimension
        r_span = np.linspace(0, self.radius, num_elements)
        r_bin_width = self.radius / num_elements
        alpha_elements = 3 * num_elements       # to scale wall x-axis appropriately
        alpha_span = np.linspace(-np.pi, np.pi, alpha_elements)
        alpha_bin_width = 2*np.pi / alpha_elements
        z_span = np.linspace(-self.height/2, self.height/2, num_elements)
        z_bin_width = self.height / num_elements

        # Mesh
        self.lid_r, self.lid_a = np.meshgrid(r_span, alpha_span)
        self.wall_a, self.wall_z = np.meshgrid(alpha_span, z_span)
        self.base_r, self.base_a = np.meshgrid(r_span, alpha_span)

        # Propagate count values into meshes
        # Lid
        self.lid_counts = np.zeros((alpha_elements, num_elements))
        for n in range(len(self.lid_impact_coords)):
            # r_bin
            r_coord = self.lid_radii[n]
            r_bin = int(np.floor(r_coord / r_bin_width))
            if r_bin == num_elements:
                r_bin = num_elements - 1            # Account for edges

            # alpha_bin
            alpha_coord = self.lid_alpha[n]
            alpha_bin = int(np.floor((alpha_coord + np.pi) / alpha_bin_width))
            if alpha_bin == alpha_elements:
                alpha_bin = alpha_elements - 1        # Account for edges

            # add photon to [alpha_bin, r_bin]
            self.lid_counts[alpha_bin, r_bin] = self.lid_counts[alpha_bin, r_bin] + 1

        # Wall
        self.wall_counts = np.zeros((num_elements, alpha_elements))
        for n in range(len(self.wall_impact_coords)):
            # alpha_bin
            alpha_coord = self.wall_alpha[n]
            alpha_bin = int(np.floor((alpha_coord + np.pi) / alpha_bin_width))
            if alpha_bin == alpha_elements:
                alpha_bin = alpha_elements - 1        # Account for edges

            # z_bin
            z_coord = self.wall_height[n]
            z_bin = int(np.floor((z_coord + self.height/2) / z_bin_width))
            if z_bin == num_elements:
                z_bin = num_elements - 1        # Account for edges

            # add photon to [z_bin, alpha_bin]
            self.wall_counts[z_bin, alpha_bin] = self.wall_counts[z_bin, alpha_bin] + 1

        # Base
        self.base_counts = np.zeros((alpha_elements, num_elements))
        for n in range(len(self.base_impact_coords)):
            # r_bin
            r_coord = self.base_radii[n]
            r_bin = int(np.floor(r_coord / r_bin_width))
            if r_bin == num_elements:
                r_bin = num_elements - 1  # Account for edges

            # alpha_bin
            alpha_coord = self.base_alpha[n]
            alpha_bin = int(np.floor((alpha_coord + np.pi) / alpha_bin_width))
            if alpha_bin == alpha_elements:
                alpha_bin = alpha_elements - 1  # Account for edges

            # add photon to [alpha_bin, r_bin]
            self.base_counts[alpha_bin, r_bin] = self.base_counts[alpha_bin, r_bin] + 1

    """"""

    def make_scatter_plot(self, file_name, show):
        # Initialize whole figure
        fig = plt.figure(figsize=(10, 10))

        # Create lid subplot
        lid = fig.add_subplot(3, 1, 1, projection='polar')
        lid.set_theta_zero_location("S")
        lid.scatter(self.lid_alpha, self.lid_radii)
        self.convert_polar_xticks_to_radians(lid)
        lid.set_rticks(np.arange(0, self.radius, 500))
        lid.set_rlabel_position(135)

        # Create wall subplot
        wall = fig.add_subplot(3, 1, 2)
        wall.scatter(self.wall_alpha, self.wall_height)
        plt.xlim(-np.pi, np.pi)
        wall.set_xticks(np.arange(-np.pi, np.pi, np.pi/4))
        self.convert_polar_xticks_to_radians(wall)
        plt.ylim(-self.z_limit, self.z_limit)
        plt.grid()

        # Create base subplot
        base = fig.add_subplot(3, 1, 3, projection='polar')
        base.set_theta_zero_location("N")
        base.set_theta_direction(-1)
        base.scatter(self.base_alpha, self.base_radii)
        self.convert_polar_xticks_to_radians(base)
        base.set_rticks(np.arange(0, self.radius, 500))
        base.set_rlabel_position(135)

        fig.tight_layout()

        if show:
            fig.show()
        fig.savefig(file_name)
        plt.close()

    """"""

    def make_heatmap(self, file_name, show):
        # Find maximum counts in the simulation
        lid_max = np.amax(self.lid_counts)
        wall_max = np.amax(self.wall_counts)
        base_max = np.amax(self.base_counts)
        self.counts_max = max(lid_max, wall_max, base_max)

        # Initialize whole figure
        fig = plt.figure(figsize=(12, 12))

        lid = plt.subplot2grid((3, 3), (0, 1), projection='polar')
        wall = plt.subplot2grid((3, 3), (1, 0), colspan=3)
        base = plt.subplot2grid((3, 3), (2, 1), projection='polar')

        # Create lid subplot
        lid.set_theta_zero_location("S")
        lid.pcolormesh(self.lid_a, self.lid_r, self.lid_counts, vmin=0, vmax=self.counts_max, cmap='viridis')
        self.convert_polar_xticks_to_radians(lid)
        lid.set_rticks(np.arange(0, self.radius, 500))
        lid.tick_params(axis='y', colors='orange')
        lid.set_rlabel_position(215)
        lid.grid()

        # Create wall subplot
        mid = wall.imshow(self.wall_counts, vmin=0, vmax=self.counts_max, cmap='viridis')
        wall.set_ylim(wall.get_ylim()[::-1])
        wall.grid()

        # Create base subplot
        base.set_theta_zero_location("N")
        base.set_theta_direction(-1)
        base.pcolormesh(self.base_a, self.base_r, self.base_counts, vmin=0, vmax=self.counts_max, cmap='viridis')
        self.convert_polar_xticks_to_radians(base)
        base.set_rticks(np.arange(0, self.radius, 500))
        base.tick_params(axis='y', colors='orange')
        base.set_rlabel_position(165)
        base.grid()

        cbar_ax = fig.add_axes([0.8, 0.05, 0.05, 0.25])
        fig.colorbar(mid, cax=cbar_ax)

        plt.tight_layout()
        if show:
            plt.show()
        plt.savefig(file_name)
        plt.close()

    """"""

    def format_radians_label(self, float_in):
        # Converts a float value in radians into a
        # string representation of that float
        string_out = str((float_in / np.pi)) + "π"

        return string_out

    """"""

    def convert_polar_xticks_to_radians(self, ax):
        # Converts x-tick labels from degrees to radians

        # Get the x-tick positions (returns in radians)
        label_positions = ax.get_xticks()

        # Convert to a list since we want to change the type of the elements
        labels = list(label_positions)

        # Convert from [0, 2*np.pi] to [-np.pi, np.pi]
        for n in range(len(labels)):
            val = labels[n]
            if val > np.pi:
                labels[n] = val - 2*np.pi

        # Format each label
        labels = [self.format_radians_label(angle) for angle in labels]

        # Keep xtick locations the same but change labels to new labels
        ax.set_xticks(label_positions)
        ax.set_xticklabels(labels)

    """"""

    def create_intensity(self, file_name, show=False):
        # Initialize
        intensities = []
        num_bins = int(self.counts_max) + 1
        tick_freq = round(num_bins / 10)
        if tick_freq < 1:
            tick_freq = 1

        # Extract data
        for n in range(len(self.lid_counts)):
            for m in range(len(self.lid_counts[n])):
                counts = int(self.lid_counts[n][m])
                if counts != 0:     intensities.append(counts)

        for n in range(len(self.wall_counts)):
            for m in range(len(self.wall_counts[n])):
                counts = int(self.wall_counts[n][m])
                if counts != 0:     intensities.append(counts)

        for n in range(len(self.base_counts)):
            for m in range(len(self.base_counts[n])):
                counts = int(self.base_counts[n][m])
                if counts != 0:     intensities.append(counts)

        # Organize data
        num_inten = len(intensities)
        frac_bins = []
        log_steps = []
        exp = 0
        log_steps.append(exp)       # len(edges) must be len(values) + 1
        cum_frac = 0
        while cum_frac < 0.99:
            cumulative_count = 0
            log_limit = 10**exp
            for n in range(num_inten):
                if intensities[n] < log_limit:
                    cumulative_count = cumulative_count + 1
            cum_frac = cumulative_count / num_inten
            frac_bins.append(cum_frac)
            exp = exp + 1
            log_steps.append(exp)

        # Plot
        fig = plt.figure()
        plt.stairs(frac_bins, log_steps)
        plt.title('Cumulative Intensity Histogram')
        plt.xlabel('Log10 of Number of Photons')
        plt.ylabel('Number of Bins (Fraction of total)')
        plt.grid()
        if show:
            fig.show()
        fig.savefig(file_name)
        plt.close(fig)

    """"""
