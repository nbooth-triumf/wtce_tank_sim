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
