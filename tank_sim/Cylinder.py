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
                alpha = np.arctan2(y_cap, x_cap)

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

    def create_graphic(self, file_name, show=False):
        # Put data into plottable form
        # Lid
        lid_radii = []
        lid_alpha = []
        for n in range(len(self.lid_impact_coords)):
            [radius, alpha, z] = self.lid_impact_coords[n]
            lid_radii.append(radius)
            lid_alpha.append(alpha)

        # Wall
        wall_alpha = []
        wall_height = []
        for n in range(len(self.wall_impact_coords)):
            [radius, alpha, z] = self.wall_impact_coords[n]
            wall_alpha.append(alpha)
            wall_height.append(z)

        # Base
        base_radii = []
        base_alpha = []
        for n in range(len(self.base_impact_coords)):
            [radius, alpha, z] = self.base_impact_coords[n]
            base_radii.append(radius)
            base_alpha.append(alpha)

        # Plot
        fig = plt.figure()
        lid = fig.add_subplot(3, 1, 1, projection='polar')
        lid.scatter(lid_alpha, lid_radii)
        wall = fig.add_subplot(3, 1, 2)
        wall.scatter(wall_alpha, wall_height)
        plt.xlim(-np.pi, np.pi)
        plt.ylim(-self.z_limit, self.z_limit)
        plt.grid()
        base = fig.add_subplot(3, 1, 3, projection='polar')
        base.scatter(base_alpha, base_radii)

        if show:
            fig.show()
        fig.savefig(file_name)
        plt.close()

    """"""

    def brute_force_coords(self, directional_coords, source_height_m):
        # convert to mm
        source_height = source_height_m * 1000  # [=] mm
        self.num_photons = len(directional_coords)

        # iterate through each photon coordinates
        for n in range(self.num_photons):
            # Extract data
            [cos_th, phi] = directional_coords[n]
            source_theta = np.arccos(cos_th)

            # Determine quadrant of impact from phi. It is easier to determine the impact location
            # in quadrant 1 then assign the correct signage than to track signage as we go.
            source_phi = phi        # Initialize
            horizontal = 0          # Initialize
            vertical = 0            # Initialize
            if 0 <= phi < np.pi/2:
                # Stay in quadrant 1
                horizontal = 1
                vertical = 1
                source_phi = source_phi
            elif np.pi/2 <= phi < np.pi:
                # Move from quadrant 2 to 1
                horizontal = -1
                vertical = 1
                source_phi = np.pi - source_phi
            elif -np.pi/2 <= phi < 0:
                # Move from quadrant 3 to 1
                horizontal = 1
                vertical = -1
                source_phi = np.abs(source_phi)
            elif -np.pi <= phi < -np.pi/2:
                # Move from quadrant 4 to 1
                horizontal = -1
                vertical = -1
                source_phi = np.pi - np.abs(source_phi)
            else:
                print("Invalid phi angle.")

            """  
            Point D is directly opposite the source O, separated by the diameter of the tank. Following
            theta and phi, the lightray (OP) from source O will hit the tangent plane of point D at point P,
            creating a right triangle OBP with point B on D's tangent plane. This lightray intersects with
            the tank cylinder at an point A, creating another right triangle OCA with point C directly
            above the source O and at the same height of intersection point A.
            
            Point A can therefore be defined by a height h, the distance between O and A along the z axis,
            and an angle gamma, the number of radians between the radius to point A and the diameter OD.
            """

            # Find angle gamma, the radians of angleAZE, where E is directly above D and across from C.
            segmentDP = self.diameter * np.tan(source_theta)
            angleBDP = source_phi
            segmentBD = segmentDP * np.cos(angleBDP)
            angleBOD = np.arctan(segmentBD / self.diameter)
            angleAZE = 2 * angleBOD                             # gamma

            # Find impact height h_impact, the length of semgentOC
            segmentBP = segmentDP * np.sin(angleBDP)
            segmentOB = np.sqrt(segmentBD**2 + self.diameter**2)
            anglePOB = np.arctan(segmentBP / segmentOB)
            angleCAO = anglePOB                                 # Properties of parallel lines
            angleAZC = np.pi - angleAZE
            segmentCA = self.diameter * np.sin(angleAZC / 2)    # Chord based trigonometry
            segmentOC = segmentCA * np.tan(angleCAO)            # h_impact

            # Apply signage to gamma and h_impact
            gamma = horizontal * angleAZE
            h_impact = vertical * segmentOC
            impact_height = source_height + h_impact  # Impact height relative to bottom of tank

            # Determine if impact point A is within tank dimensions
            if 0 <= impact_height < self.height:
                # Impact is on tank wall -> save collision coordinates
                self.wall_impact_coords.append([impact_height, gamma])
            else:
                # Assume impact happens at point L on endcap. Triangle OCA is similar to new triangle OKL,
                # where K is the edge of the tank directly below O. Remove signage as before.
                if impact_height < 0:
                    segmentOK = source_height
                else:
                    segmentOK = self.height - source_height

                segmentKL = segmentCA * (segmentOK / segmentOC)

                # Angle between segmentKL and the diameter of the endcap is equal to angleBOD, so
                # components can be found
                segmentKL_1 = segmentKL * np.sin(angleBOD)
                segmentKL_2 = segmentKL * np.cos(angleBOD)

                # Point L is therefore offset from the z axis by segmentKL_1 in the x' direction and
                # by segmentKL_2 - tank_radius in the y' direction
                r_x = segmentKL_1
                r_y = segmentKL_2 - self.radius

                # Determine impact coords [r, alpha] on endcaps
                r = np.sqrt(r_x**2 + r_y**2)
                alpha = np.arctan(r_y/r_x)

                # Append impact coords to correct data structure
                if impact_height < 0:
                    self.base_impact_coords.append([r, alpha])
                else:
                    self.lid_impact_coords.append([r, alpha])

    """"""
