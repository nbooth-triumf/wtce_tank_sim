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

        self.troubleshooting_info = []      # Data structure to save relevant values if error occurs

    """"""

    def generate_impact_coords(self, directional_coords, source_height_m):
        # convert to mm
        source_height = source_height_m * 1000      # [=] mm
        self.num_photons = len(directional_coords)

        for n in range(self.num_photons):
            # Extract data
            [cos_th, phi] = directional_coords[n]
            theta = np.arccos(cos_th)

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
            delta_z_x = x_impact * (z_hat / x_hat)
            delta_z_y = (y_impact + self.radius) * (z_hat / y_hat)
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
        self.make_scatter_plot(file_name, show)
        self.make_heatmap(file_name, show)

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
        alpha_span = np.linspace(0, 2*np.pi, num_elements)
        z_span = np.linspace(-self.height/2, self.height/2, num_elements)

        # Mesh
        self.lid_r, self.lid_a = np.meshgrid(r_span, alpha_span)
        self.wall_a, self.wall_z = np.meshgrid(alpha_span, z_span)
        self.base_r, self.base_a = np.meshgrid(r_span, alpha_span)

    """"""

    def make_scatter_plot(self, file_name, show):
        # Initialize whole figure
        fig = plt.figure(figsize=(10, 10))

        # Create lid subplot
        lid = fig.add_subplot(3, 1, 1, projection='polar')
        lid.set_theta_zero_location("S")
        lid.scatter(self.lid_alpha, self.lid_radii)
        self.convert_polar_xticks_to_radians(lid)

        # Create wall subplot
        wall = fig.add_subplot(3, 1, 2)
        wall.scatter(self.wall_alpha, self.wall_height)
        plt.xlim(-np.pi, np.pi)
        plt.ylim(-self.z_limit, self.z_limit)
        plt.grid()

        # Create base subplot
        base = fig.add_subplot(3, 1, 3, projection='polar')
        base.set_theta_zero_location("N")
        base.scatter(self.base_alpha, self.base_radii)
        self.convert_polar_xticks_to_radians(base)
        fig.tight_layout()

        if show:
            fig.show()
        fig.savefig(file_name)
        plt.close()

    """"""

    def make_heatmap(self, file_name, show):
        print("hewwo world")

    """"""

    def format_radians_label(self, float_in):
        # Converts a float value in radians into a
        # string representation of that float
        string_out = str(float_in / np.pi) + "π"

        return string_out

    """"""

    def convert_polar_xticks_to_radians(self, ax):
        # Converts x-tick labels from degrees to radians

        # Get the x-tick positions (returns in radians)
        label_positions = ax.get_xticks()

        # Convert to a list since we want to change the type of the elements
        labels = list(label_positions)

        # Format each label
        labels = [self.format_radians_label(angle) for angle in labels]

        # Keep xtick locations the same but change labels to new labels
        ax.set_xticks(label_positions)
        ax.set_xticklabels(labels)

    """"""
