import csv
import os
import numpy as np
import matplotlib.pyplot as plt


class Histogram(object):
    """
    Base class for histograms of the angular distribution of light source
    """

    def __init__(self, description, data_folder, sig_figs=3):
        """
        Constructor
        """
        self.description = description      # cos(th) or phi with separation distance
        self.directory = data_folder        # Folder where histograms can be found
        self.sig_figs = sig_figs            # Number of decimal points for rounding

        self.num_files = None               # Number of files in directory
        self.file_paths = []                # File paths for each file in directory
        self.file_names = []                # File names for each file
        self.file_types = []                # File extension for each file
        self.file_data = []                 # Data structure of all datasets from all files
        self.file_data_start = []           # First full, nonzero bin in a file's data

        self.full_hist = []                 # Combined histogram of all files in directory
        self.normalized = []                # Normalized so area under curve = 1
        self.cumulative = []                # Separate histogram of cumulative probability

        self.num_bins = None                # Number of bins in histogram

        self.open_all_files()

    """"""

    def open_all_files(self):
        # Iterate over all files in given directory and save the file paths
        for file_path in os.listdir(self.directory):
            f = os.path.join(self.directory, file_path)
            # checking if it is a file
            if os.path.isfile(f):
                self.file_paths.append(f)
        self.num_files = len(self.file_paths)

        # Extract data sets from files
        for file_path in self.file_paths:
            # Extract and save relevant identification information
            file_name = os.path.splitext(file_path)[0]
            self.file_names.append(file_name)
            file_type = os.path.splitext(file_path)[1]
            self.file_types.append(file_type)

            # Open file and identify first full bin
            [data, num_bins] = self.open_csv(file_path)
            self.file_data.append(data)
            first_full_bin = self.find_first_bin(data)
            self.file_data_start.append(first_full_bin)

            # Confirm num_bins is constant for each file
            if self.num_bins is None:
                self.num_bins = num_bins
            elif self.num_bins != num_bins:
                print('Error: Data files are not compatible.')

    """"""

    def combine_files(self):
        # Combine datasets using found starting bin numbers
        self.full_hist = np.zeros(self.num_bins)

        for n in range(self.num_files):
            # Define Start
            start_index = self.file_data_start[n]
            current_index = start_index

            # Define Stop
            if n == self.num_files-1:
                stop_index = self.num_bins
            else:
                stop_index = self.file_data_start[n + 1]

            # Iterate to stop
            while current_index < stop_index:
                self.full_hist[current_index] = self.file_data[n][current_index]
                current_index = current_index + 1

            # Adjust histogram values to make distribution continuous at stop
            # at this point, current_index = stop_index
            if n == self.num_files - 1:
                # Last file does need to be adjusted
                continue
            else:
                # Adjust file_data[n] bins according to the normalized offset
                current_stop_data = self.file_data[n][current_index]
                next_start_data = self.file_data[n+1][current_index]
                for m in range(start_index, stop_index):
                    # Normalize current data point then scale to next data point
                    unadjusted = self.file_data[n][m]
                    norm_val = unadjusted / current_stop_data
                    adjusted = int(round(norm_val * next_start_data, self.sig_figs))
                    self.full_hist[m] = adjusted

    """"""

    def open_csv(self, file_path):
        # Initialize
        file_data = []

        # Open csv and save each line as a data point
        with open(file_path, 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            for row in csv_reader:
                file_data.append(row[0])

        # Convert numbers in list to integers and extract number of bins
        file_data = list(map(int, file_data))
        file_bins = len(file_data)

        return [file_data, file_bins]

    """"""

    def find_first_bin(self, dataset):
        # Initialize
        num_bins = len(dataset)
        current_index = 0
        first_full_bin = None

        # Iterate to second to last bin
        while current_index < (num_bins - 1):
            # Extract
            current_data = dataset[current_index]

            # If bin is empty or negative, proceed to next bin
            if current_data <= 0:
                current_index = current_index + 1

            # If bin is positive nonzero, investigate further
            else:
                # Compare current bin to next bin
                next_data = dataset[current_index + 1]
                comparison = (next_data - current_data) / current_data

                # If change from current bin to next bin is much greater than the current bin value,
                # the current bin is not full.
                if comparison > 1:
                    current_index = current_index + 1
                # otherwise, save it as first full bin
                else:
                    first_full_bin = current_index
                    break

        # Return
        if first_full_bin is None:
            # no full bin found before final bin; current_index is equal to num_bins - 1
            if dataset[current_index] > 0:
                # Final bin is first full bin
                return current_index
            else:
                # Dataset is empty. Alert user.
                print('Dataset has no full bins. Element 0 is returned.')
                return 0
        else:
            # If first_full_bin is defined, return it.
            return first_full_bin

    """"""

    def create_normalized(self):
        # Sum up all bins in dataset
        summed_total = sum(self.full_hist)

        # Determine amount of total in each bin and save
        for n in range(self.num_bins):
            value = self.full_hist[n]
            norm_value = round(value / summed_total, self.sig_figs)
            self.normalized.append(norm_value)

    """"""

    def create_cumulative(self):
        # Iterate through each bin
        for n in range(self.num_bins):
            # Determine total of current bin and all bins to left
            sum_to = n + 1
            sum_to_bin = sum(self.normalized[0:sum_to])

            # Round to avoid residual errors
            sum_to_bin = round(sum_to_bin, self.sig_figs)

            # Save
            self.cumulative.append(sum_to_bin)

        # Force cumulative to end at 1.0
        if self.cumulative[self.num_bins - 1] > 1.0:
            self.cumulative[self.num_bins - 1] = 1.0

    """"""

    def plot_dists(self, bin_range, x_label, path_header):
        # Plot all distributions
        # combined full histogram
        label = 'full_hist'
        self.plot_histogram(bin_range, self.full_hist, label, x_label, path_header+'/'+label+'.png')

        # Normalized
        label = 'normalized_dist'
        self.plot_histogram(bin_range, self.normalized, label, x_label, path_header+'/'+label+'.png')

        # Cumulative
        label = 'cumulative_dist'
        self.plot_histogram(bin_range, self.cumulative, label, x_label, path_header+'/'+label+'.png')

    """"""

    def plot_histogram(self, bin_range, distribution, title, x_label, file_name, show=False):
        # Define bins
        xmin = bin_range[0]
        xstep = bin_range[1] - bin_range[0]
        xmax = bin_range[-1] + xstep
        bins_plot = np.linspace(xmin, xmax, self.num_bins+1)

        # Plot
        fig = plt.figure()
        plt.stairs(distribution, bins_plot)
        plt.title(title)
        plt.xlim((xmin, xmax))
        plt.xlabel(x_label)
        plt.grid()
        if show:
            fig.show()
        fig.savefig(file_name)

    """"""
