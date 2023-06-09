import matplotlib.pyplot as plt
import numpy as np


class Graphics(object):
    """
    Base class for plotting all figures
    """

    def __init__(self, testCylinder):
        """
        Constructor - HEIGHT AND DIAMETER IN METRES
        """

        myCylinder = testCylinder
        self.num_photons = myCylinder.num_photons
        self.radius = myCylinder.radius
        self.z_limit = myCylinder.z_limit

        self.lid_areas = myCylinder.lid_pmt_areas
        self.wall_areas = myCylinder.wall_pmt_areas
        self.base_areas = myCylinder.base_pmt_areas

        self.lid_counts = myCylinder.lid_counts
        self.wall_counts = myCylinder.wall_counts
        self.base_counts = myCylinder.base_counts
        self.lid_pmt_counts = myCylinder.lid_pmt_counts
        self.wall_pmt_counts = myCylinder.wall_pmt_counts
        self.base_pmt_counts = myCylinder.base_pmt_counts
        self.lid_log_counts = myCylinder.lid_log_counts
        self.wall_log_counts = myCylinder.wall_log_counts
        self.base_log_counts = myCylinder.base_log_counts

        self.counts_max = myCylinder.counts_max
        self.counts_pmt_max = myCylinder.counts_pmt_max
        self.counts_log_max = myCylinder.counts_log_max

        self.lid_alpha = myCylinder.lid_alpha
        self.lid_radii = myCylinder.lid_radii
        self.wall_alpha = myCylinder.wall_alpha
        self.wall_height = myCylinder.wall_height
        self.base_alpha = myCylinder.base_alpha
        self.base_radii = myCylinder.base_radii

        self.lid_a = myCylinder.lid_a
        self.lid_r = myCylinder.lid_r
        self.base_a = myCylinder.base_a
        self.base_r = myCylinder.base_r

        self.intensities = []
        self.log_intensities = []
        self.weights = []
        self.log_weights = []

    """"""

    def create_graphics(self, file_name, data_description, show=False):
        # Plot scatter plot of simulation
        scatter_name = file_name + "_scatter"
        #  self.make_scatter_plot(scatter_name, show)

        # Standardize colour bars and create heatmaps for raw, adjusted, and log data
        uniform = False  # For troubleshooting
        v_max = max(self.counts_max, self.counts_pmt_max, self.counts_log_max)
        heatmap_name = file_name + "_heatmap"
        if uniform:
            self.make_heatmap(heatmap_name, v_max, data_description, label='raw', show=show)
            self.make_heatmap(heatmap_name, v_max, data_description, label='adjusted', show=show)
            self.make_heatmap(heatmap_name, v_max, data_description, label='log_scale', show=show)
        else:
            self.make_heatmap(heatmap_name, self.counts_max, data_description, label='raw', show=show)
            self.make_heatmap(heatmap_name, self.counts_pmt_max, data_description, label='adjusted', show=show)
            self.make_heatmap(heatmap_name, self.counts_log_max, data_description, label='log_scale', show=show)

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
        wall.set_xticks(np.arange(-np.pi, np.pi, np.pi / 4))
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

    def make_heatmap(self, file_name, v_max, data_description, label, show=False):
        # Initialize to keep Python from yelling about warnings
        lid_to_plot = None
        wall_to_plot = None
        base_to_plot = None
        v_max_plot = v_max

        # Define data to plot
        if label == 'raw':
            lid_to_plot = self.lid_counts
            wall_to_plot = self.wall_counts
            base_to_plot = self.base_counts
        elif label == 'adjusted':
            # Convert to counts per square centimetre
            lid_to_plot = self.lid_pmt_counts / 10**4
            wall_to_plot = self.wall_pmt_counts / 10**4
            base_to_plot = self.base_pmt_counts / 10**4
            v_max_plot = v_max / 10**4
        elif label == 'log_scale':
            lid_to_plot = self.lid_log_counts
            wall_to_plot = self.wall_log_counts
            base_to_plot = self.base_log_counts
        else:
            print("Plot label not recognized.")

        # Initialize whole figure
        fig = plt.figure(figsize=(12, 12))
        plot_title = 'Heatmap of Data Simulation using ' + label + ' data from ' + data_description
        plt.title(plot_title)

        lid = plt.subplot2grid((3, 3), (0, 1), projection='polar')
        wall = plt.subplot2grid((3, 3), (1, 0), colspan=3)
        base = plt.subplot2grid((3, 3), (2, 1), projection='polar')

        # Create lid subplot
        lid.set_theta_zero_location("S")
        lid.pcolormesh(self.lid_a, self.lid_r, lid_to_plot, vmin=0, vmax=v_max_plot, cmap='viridis')
        self.convert_polar_xticks_to_radians(lid)
        lid.set_rticks(np.arange(0, self.radius, 500))
        lid.tick_params(axis='y', colors='orange')
        lid.set_rlabel_position(215)
        lid.grid()

        # Create wall subplot
        mid = wall.imshow(wall_to_plot, vmin=0, vmax=v_max_plot, cmap='viridis')
        wall.set_ylim(wall.get_ylim()[::-1])
        wall.grid()

        # Create base subplot
        base.set_theta_zero_location("N")
        base.set_theta_direction(-1)
        base.pcolormesh(self.base_a, self.base_r, base_to_plot, vmin=0, vmax=v_max_plot, cmap='viridis')
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
        plt.savefig(file_name + "_" + label)
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
                labels[n] = val - 2 * np.pi

        # Format each label
        labels = [self.format_radians_label(angle) for angle in labels]

        # Keep xtick locations the same but change labels to new labels
        ax.set_xticks(label_positions)
        ax.set_xticklabels(labels)

    """"""

    def create_intensity(self, file_name, data_group, pipe_label, log, show=False):
        if log:
            self.log_intensity(file_name + "_log", data_group, pipe_label, show)
        else:
            self.pmt_intensity(file_name, data_group, pipe_label, show)

    """"""

    def pmt_intensity(self, file_name, data_group, pipe_label, show=False):
        # Extract data
        for n in range(len(self.lid_pmt_counts)):
            for m in range(len(self.lid_pmt_counts[n])):
                # convert to counts per square centimetre
                counts = int(self.lid_pmt_counts[n][m]) / 10**4
                area = self.lid_areas[n][m] * 10**4
                if counts != 0:
                    self.intensities.append(counts)
                    self.weights.append(area)

        for n in range(len(self.wall_pmt_counts)):
            for m in range(len(self.wall_pmt_counts[n])):
                # convert to counts per square centimetre
                counts = int(self.wall_pmt_counts[n][m]) / 10**4
                area = self.wall_areas[n][m] * 10**4
                if counts != 0:
                    self.intensities.append(counts)
                    self.weights.append(area)

        for n in range(len(self.base_pmt_counts)):
            for m in range(len(self.base_pmt_counts[n])):
                # convert to counts per square centimetre
                counts = int(self.base_pmt_counts[n][m]) / 10**4
                area = self.base_areas[n][m] * 10**4
                if counts != 0:
                    self.intensities.append(counts)
                    self.weights.append(area)

        # Plot
        fig = plt.figure()
        plt.hist(x=self.intensities, weights=self.weights)
        plt.xlim(left=0)
        plt.ylim(bottom=0)
        plt.title('Intensity Histogram for ' + pipe_label)
        plt.xlabel('Number of Photons per cm')
        plt.ylabel('Number of Bins')
        plt.grid()
        if show:
            fig.show()
        to_save = file_name + '_' + data_group
        fig.savefig(to_save)
        plt.close(fig)

    """"""

    def log_intensity(self, file_name, data_group, pipe_label, show=False):
        # Extract data
        for n in range(len(self.lid_log_counts)):
            for m in range(len(self.lid_log_counts[n])):
                log_counts = self.lid_log_counts[n][m]
                bin_area = self.lid_areas[n][m]
                counts = (log_counts / bin_area) / 10**4
                area = bin_area * 10**4
                if counts != 0:
                    self.log_intensities.append(counts)
                    self.log_weights.append(area)

        for n in range(len(self.wall_log_counts)):
            for m in range(len(self.wall_log_counts[n])):
                log_counts = self.wall_log_counts[n][m]
                bin_area = self.wall_areas[n][m]
                counts = (log_counts / bin_area) / 10**4
                area = bin_area * 10**4
                if counts != 0:
                    self.log_intensities.append(counts)
                    self.log_weights.append(area)

        for n in range(len(self.base_log_counts)):
            for m in range(len(self.base_log_counts[n])):
                log_counts = self.base_log_counts[n][m]
                bin_area = self.base_areas[n][m]
                counts = (log_counts / bin_area) / 10**4
                area = bin_area * 10**4
                if counts != 0:
                    self.log_intensities.append(counts)
                    self.log_weights.append(area)

        # Bins
        log_step = 0.1
        max_log = int(np.ceil(np.max(self.log_intensities)))
        bins = np.arange(0, max_log, log_step)

        # Plot
        fig = plt.figure()
        plt.hist(x=self.log_intensities, bins=bins, weights=self.log_weights)
        plt.xlim([0, max_log])
        plt.ylim(bottom=0)
        plt.title('Log10 Intensity Histogram for ' + pipe_label)
        plt.xlabel('Log10 of Number of Photons per cm')
        plt.ylabel('Number of Bins')
        plt.grid()
        if show:
            fig.show()
        to_save = file_name + '_' + data_group
        fig.savefig(to_save)
        plt.close(fig)

    """"""
