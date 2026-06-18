import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from pathlib import Path

#Set Paths for analysis
data_paths = ["C:/Users/mbl-4/Desktop/DataAnalysis_MBL_Neurobiology_2026_OBTeam/2026_MBL_OB_Project/gc_mus_EPSC/Exports/260609_0_char.csv",
              #"E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260611_0_char.csv",
              #"E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260611_3_char.csv",
              #"E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260612_6_char.csv",
            "C:/Users/mbl-4/Desktop/DataAnalysis_MBL_Neurobiology_2026_OBTeam/2026_MBL_OB_Project/Exports/260611_3_char.csv"]
              #"E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260613_2_char.csv",
              #"E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260606_1_char.csv"] #should be char file
              # "E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260609_0_char.csv"

#Initiate list to hold pandas dataframes
data = []
#Load in data
for path in data_paths:
    df = pd.read_csv(path, encoding='cp1252')
    data.append(df)


def plot_dotplot(data_1, dot_color, plot_title, y_axis_label):
    """
    Function to make a dot plot with mean (± SEM)

    :param data_1: List of datapoints for dataset 1
    :param dot_color: Hexcode color for first group
    :param plot_title: Title of plot as a string
    :param y_axis_label: y-axis label for graph

    Returns: matplotlib figure object
    """

    # dot edge properties
    dot_edge_color = 'black'
    dot_edge_width = 1

    # Font sizes
    title_fontsize = 18
    axis_label_fontsize = 16

    # Figure size (width, height in inches)
    fig_width = 4
    fig_height = 6

    # ==============================================================================
    # CREATE THE PLOT
    # ==============================================================================

    # Calculate means for bar heights
    data_1_mean = np.mean(data_1)

    # Calculate SEM
    #data_1_sem = np.std(data_1, ddof=1) / np.sqrt(len(data_1))

    # Set up the figure with high DPI for publication quality
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)

    # X-axis positions
    x_pos = 0

    # Add individual data points with jitter
    np.random.seed(42)  # For reproducible jitter

    # Plot dataset 1 data points
    x_jitter_1 = np.random.normal(0, 0.05, size=len(data_1))
    ax.scatter(
        np.full(len(data_1), x_pos),
        data_1,
        color='#FFFFFF',
        s=50,
        alpha=0.6,
        edgecolors='black',
        linewidths=0.5,
        zorder=3
    )

    # Plot mean as a horizontal line
    ax.plot(
        [x_pos - 0.03, x_pos + 0.03],
        [data_1_mean, data_1_mean],
        color=dot_color,
        linewidth=3,
        solid_capstyle='butt',
        zorder=4
    )

    # Set labels and title
    ax.set_ylabel(y_axis_label, fontsize=axis_label_fontsize, fontweight='bold')
    ax.set_title(plot_title, fontsize=title_fontsize, fontweight='bold', pad=20)

    # Set x-axis ticks and labels
    ax.set_xticks([])
    ax.set_xticklabels([])

    # Set y-axis tick label size
    ax.tick_params(axis='y', labelsize=12)

    # Remove top and right spines (borders)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Make left and bottom spines thicker
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

    # Determine appropriate y-axis limits
    all_data = np.concatenate([data_1])
    y_min = 0
    y_max = max(0, np.max(all_data), data_1_mean)

    # Add some padding (10% of the range)
    y_range = y_max - y_min
    y_min = 0
    y_max += 0.1 * y_range

    ax.set_ylim(y_min, y_max)
    ax.set_xlim(-0.1, 0.1)
    # Adjust layout to prevent label cutoff
    fig.tight_layout()

    # Return the figure object
    return fig


#############################
# Plot Cell resistance
#############################

cell_R = []

# Extract data
for d in data:
    # If each file has one resistance value, take that value directly
    cell_R_val = d['Cell Resistance B (Mohms)'].iloc[1]
    cell_R.append(cell_R_val)

# Plot titles and labels
cell_R_title = 'Cell Resistance'
cell_R_label = 'Resistance (Mohms)'

# Line Color
line1_color = "#459B3D"

# Make dot plot with mean
cell_R_plot = plot_dotplot(cell_R, line1_color, cell_R_title, cell_R_label)

#############################
# Plot Access resistance
#############################

access_R = []

# Extract data
for d in data:
    access_R_val = d['Access Resistance B (Mohms)'].iloc[1]
    access_R.append(access_R_val)

# Plot titles and labels
access_R_title = 'Access Resistance'
access_R_label = 'Access Resistance (Mohms)'

# Line Color
line1_color = "#459B3D"

# Make dot plot with mean
access_R_plot = plot_dotplot(access_R, line1_color, access_R_title, access_R_label)
#############################
# Plot Cell Capacitance
#############################

cell_C = []

# Extract data
for d in data:
    cell_C_val = d['Capacitance B (pF)'].iloc[1]
    cell_C.append(cell_C_val)

# Plot titles and labels
cell_C_title = 'Cell Capacitance'
cell_C_label = 'Capacitance (pF)'

# Line Color
line1_color = "#459B3D"

# Make dot plot with mean
cell_C_plot = plot_dotplot(cell_C, line1_color, cell_C_title, cell_C_label)

################
# Export Plots #
################

path = input('Enter file path for figure export: ')
file_path = path.strip().strip('"').strip("'")

export_dir = Path(file_path)
export_dir.mkdir(parents=True, exist_ok=True)

groupname = input('Enter qualifier for group, i.e. GC for granular cells, all for all cells, etc.: ')

export_path = export_dir / f"{groupname}_VC_cell_resistance.pdf"
cell_R_plot.savefig(export_path)

export_path = export_dir / f"{groupname}_VC_access_resistance.pdf"
access_R_plot.savefig(export_path)

export_path = export_dir / f"{groupname}_VC_cell_capacitance.pdf"
cell_C_plot.savefig(export_path)