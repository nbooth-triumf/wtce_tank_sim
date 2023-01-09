import numpy as np
import csv
import os


class Dataset(object):
    """
    Base class for histograms of the angular distribution of light source
    """

    def __init__(self, description, file_path, sig_figs=3):
        """
        Constructor
        """
        self.description = description      # Description of dataset
        self.sig_figs = sig_figs            # Number of decimal points for rounding

        self.file_path = file_path                              # Full file path of image location
        self.file_name = os.path.splitext(self.file_path)[0]    # Name of file without file extension
        self.file_type = os.path.splitext(self.file_path)[1]    # Type of file (file extension only)

        self.data = []          # Extracted histogram
        self.normalized = []    # Normalized so area under curve = 1
        self.cumulative = []    # Separate histogram of cumulative probability

        self.num_bins = None    # Number of bins in histogram

    """"""

    def open_data_csv(self):
        # Open csv and save each line as a data point
        with open(self.file_path, 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            for row in csv_reader:
                self.data.append(row[0])

        # Convert numbers in list to integers
        self.data = list(map(int, self.data))

        # Extract knowns
        self.num_bins = len(self.data)

    """"""

    def normalize(self):
        # Sum up all bins in dataset
        summed_total = sum(self.data)

        # Determine amount of total in each bin and save
        for n in range(self.num_bins):
            value = self.data[n]
            norm_value = round(value / summed_total, self.sig_figs)
            self.normalized.append(norm_value)

    """"""

    def cumulative(self):
        # Iterate through each bin
        for n in range(self.num_bins):
            # Determine total of current bin and all bins to left
            sum_to_bin = sum(self.normalized[0:(n + 1)])
            # Round to avoid residual errors
            sum_to_bin = round(sum_to_bin, self.sig_figs)
            # Save
            self.cumulative.append(sum_to_bin)
