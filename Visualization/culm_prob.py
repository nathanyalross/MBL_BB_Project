from pathlib import Path
from plot_func import plot_cumulative_probability 
import pandas as pd


#Set Paths for analysis
raw_data = ['E:/MBL_neurobiology/bulb_baddies/Exports/260610_1.csv']

#Load in data
data = []
for path in raw_data:
    df = pd.read_csv(path)
    #Append to list and transpose dataframe so columns are data type and rows are sweep number
    data.append(df)

#Function to plot culmulative Probability for Amplitude (Inspect plot_func for more info)
amp_fig,ax,amp_results = plot_cumulative_probability(data, 'Amplitude (pA)')

#Function to plot culmulative Probability for Decay Time
decay_fig,ax,decay_results = plot_cumulative_probability(data, 'Decay Tau')

#Function to plot culmulative Probability for Charge (Inspect plot_func for more info)
charge_fig,ax,charge_results = plot_cumulative_probability(data, 'Charge (pC)')

#Function to plot culmulative Probability for Decay Time (Inspect plot_func for more info)
rt_fig,ax,rt_results = plot_cumulative_probability(data, 'Rise time (ms)')

################
# Export Plots #
################

#specify export path
path = input('Enter file path for figure export: ')
file_path = path.strip().strip('"').strip("'")

export_dir = Path(file_path)
export_dir.mkdir(parents=True, exist_ok=True)

groupname = input('Enter qualifier for group, i.e. SST for SST neurons, all for all cells, etc.: ')

export_path = export_dir / f"{groupname}_amp_culm_prob.pdf"
amp_fig.savefig(export_path)

export_path = export_dir / f"{groupname}_decay_culm_prob.pdf"
decay_fig.savefig(export_path)

export_path = export_dir / f"{groupname}_charge_culm_prob.pdf"
charge_fig.savefig(export_path)

export_path = export_dir / f"{groupname}_rise_culm_prob.pdf"
rt_fig.savefig(export_path)