import numpy as np


class Cylinder(object):
    """
    Base class for defining tank region and corresponding geometry
    """

    def __init__(self, height_m, diameter_m):
        """
        Constructor - HEIGHT AND DIAMETER IN METRES
        """
        self.height = height_m * 1000       # Tank height [=] mm
        self.diameter = diameter_m * 1000   # Tank diameter [=] mm
        self.radius = self.diameter / 2     # Tank radius [=] mm

        self.num_photons = None             # Initialize number of projectiles
        self.wall_impact_coords = []        # Initialize wall impact coords of [gamma, impact_height]
        self.lid_impact_coords = []         # Initialize lid_impact_coords of [r, alpha]
        self.base_impact_coords = []        # Initialize base_impact_coords of [r, alpha]

    """"""

    def generate_impact_coords(self, directional_coords, source_height_m):
        # convert to mm
        source_height = source_height_m * 1000      # [=] mm
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
