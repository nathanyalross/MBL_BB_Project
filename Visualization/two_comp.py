import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from pathlib import Path
from plot_func import plot_bar_two_datasets 
from plot_func import plot_line

#Set Paths for analysis
data_paths = ["E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260606_1_summary.csv",
              #"E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260609_0_summary.csv", #weird noise
              #"E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260611_0_summary.csv", #weird noise
              #"E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260611_3_summary.csv", #weird noise
              "E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260612_6_summary.csv",
              "E:/MBL_neurobiology/bulb_baddies/gc_mus_EPSC/Exports/260613_1_summary.csv"
              ] #should be summary file

statistics = pd.DataFrame(columns = ['Test Statistic', 'P-value', 'Test Used'])


#Initiate list to hold pandas dataframes
data = []
#Load in data
for path in data_paths:
    df = pd.read_csv(path, index_col=False)
    #set index as the sweep column
    df = df.set_index('Unnamed: 0')
    #Append to list and transpose dataframe so columns are data type and rows are sweep number
    data.append(df.T)

########################
# Plot Freq Timeseries #
########################

#First, load in and plot data
bl_freqs = []
drug_freqs = []
#Extract BL vs drug data
freq_data= pd.DataFrame()
i=0
for d in data:
    #Select data in freq column and append to dataframe
    freq = d['Event frequency (Hz)'].to_list()
    #Find average frequency for first 5 minutes
    #bl_freq_val = np.mean(freq[2:4]) original
    bl_freq_val = np.mean(freq[0:4])
    bl_freqs.append(bl_freq_val)
    #Find average frequency for minute 7-12 (drug application)
    #drug_freq_val = np.mean(freq[6:11]) original
    drug_freq_val = np.mean(freq[4:9]) #drug applied at sweep 4
    drug_freqs.append(drug_freq_val)
    #Append values to dataframe for graphing
    freq_data[f'Cell_{i}']=freq
    i=i+1


#Plot titles and labels
freq_title = 'Event Frequency'
freq_y_label = 'Frequency (Hz)'

#Line Colors (color names or hex codes)
line1_color = "#459B3D"
line2_color = "#BED2BC" 

freq_plot = plot_line(freq_data, line1_color, freq_title, freq_y_label)

#Statistical comparisons baseline vs drug
freq_title = 'Frequency'
freq_x_label = ['BL', 'Muscarine']

freq_bar, test_stat, p_val, test = plot_bar_two_datasets(bl_freqs, line1_color,  drug_freqs, line2_color, freq_x_label, freq_title, freq_y_label)
statistics.loc['Frequency'] = [test_stat, p_val, test]

########################
# Plot Charge Timeseries #
########################

# #First, load in and plot data
# bl_charge = []
# drug_charge = []
# #Extract GFP+ data
# charge_data= pd.DataFrame()
# i=0
# for d in data:
#     #Select data in charge column and append to dataframe
#     charge = d['Charge (pC)'].to_list()
#     #Find average charge for first 5 minutes
#     bl_charge_val = np.mean(charge[0:4])
#     bl_charge.append(bl_charge_val)
#     #Find average charge for minute 7-12 (drug application)
#     drug_charge_val = np.mean(charge[4:9])
#     drug_charge.append(drug_charge_val)
#     #Append values to dataframe for graphing
#     charge_data[f'Cell_{i}']=charge
#     i=i+1
#
# #Plot titles and labels
# charge_title = 'Event Charge'
# charge_y_label = 'Normalized Charge'
#
# #Line Colors (color names or hex codes)
# line1_color = "#459B3D"
# line2_color = "#BED2BC"
#
# charge_plot = plot_line(charge_data, line1_color, charge_title, charge_y_label)
#
# #Statistical comparisons
# charge_title = 'Charge'
# charge_x_label = ['BL', 'Muscarine']
#
# charge_bar, test_stat, p_val, test = plot_bar_two_datasets(bl_charge, line1_color, drug_charge, line2_color, charge_x_label, charge_title, charge_y_label)
# statistics.loc['Charge'] = [test_stat, p_val, test]

#############################
# Plot Amplitude Timeseries #
#############################

#First, load in and plot data
bl_amp = []
drug_amp = []
#Extract data
amp_data= pd.DataFrame()
i=0
for d in data:
    #Select data in amplitude column and append to dataframe
    amp = d['Amplitude (pA)'].to_list()
    #Find average amplitude for first 5 minutes
    bl_amp_val = np.nanmean(amp[0:4])
    bl_amp.append(bl_amp_val)
    #Find average amplitude for minute 7-12 (drug application)
    drug_amp_val = np.nanmean(amp[4:9])
    drug_amp.append(drug_amp_val)
    #Append values to dataframe for graphing
    amp_data[f'Cell_{i}']=amp
    i=i+1

#Plot titles and labels
amp_title = 'Event Amplitude'
amp_y_label = 'Amplitude (pA)'

#Line Colors (color names or hex codes)
line1_color = "#459B3D"
line2_color = "#BED2BC" 

amplitude_plot = plot_line(amp_data, line1_color, amp_title, amp_y_label)

#Statistical comparisons baseline vs drug
amp_title = 'Amplitude'
amp_x_label = ['BL', 'Muscarine']

amp_bar, test_stat, p_val, test = plot_bar_two_datasets(bl_amp, line1_color, drug_amp, line2_color, amp_x_label, amp_title, amp_y_label)
statistics.loc['Amplitude'] = [test_stat, p_val, test]

#############################
# Plot Rise Time Timeseries #
#############################

#First, load in and plot data
bl_rise = []
drug_rise = []
#Extract data
rise_data= pd.DataFrame()
i=0
for d in data:
    #Select data in rise time column and append to dataframe
    rise = d['Rise time (ms)'].to_list()
    #Find average rise time for first 5 minutes
    bl_rise_val = np.nanmean(rise[0:4])
    bl_rise.append(bl_rise_val)
    #Find average rise time for minute 7-12 (drug application)
    drug_rise_val = np.nanmean(rise[4:9])
    drug_rise.append(drug_rise_val)
    #Append values to dataframe for graphing
    rise_data[f'Cell_{i}']=rise
    i=i+1

#Plot titles and labels
rise_title = 'Event Rise Time'
rise_y_label = 'Rise Time (ms)'

#Line Colors (color names or hex codes)
line1_color = "#459B3D"
line2_color = "#BED2BC" 

rise_time_plot = plot_line(rise_data, line1_color, rise_title, rise_y_label)

#Statistical comparisons 
rise_title = 'Rise Time'
rise_x_label = ['BL', 'Muscarine']

rise_bar, test_stat, p_val, test = plot_bar_two_datasets(bl_rise, line1_color,  drug_rise, line2_color, rise_x_label, rise_title, rise_y_label)
statistics.loc['Rise Time'] = [test_stat, p_val, test]

#############################
# Plot Decay Time Timeseries #
#############################

#First, load in and plot data
bl_decay = []
drug_decay = []
#Extract data
decay_data= pd.DataFrame()
i=0
for d in data:
    #Select data in rise time column and append to dataframe
    decay = d['Decay tau (ms)'].to_list()
    #Find average rise time for first 5 minutes
    bl_decay_val = np.nanmean(decay[0:4])
    bl_decay.append(bl_decay_val)
    #Find average decay time for minute 7-12 (drug application)
    drug_decay_val = np.nanmean(decay[4:9])
    drug_decay.append(drug_decay_val)
    #Append values to dataframe for graphing
    decay_data[f'Cell_{i}']=decay
    i=i+1

#Plot titles and labels
decay_title = 'Event Decay Time'
decay_y_label = 'Decay Time (ms)'

#Line Colors (color names or hex codes)
line1_color = "#459B3D"
line2_color = "#BED2BC" 

decay_time_plot = plot_line(decay_data, line1_color, decay_title, decay_y_label)

#Statistical comparisons 
decay_title = 'Decay Time'
decay_x_label = ['BL', 'Muscarine']

decay_bar, test_stat, p_val, test = plot_bar_two_datasets(bl_decay, line1_color,  drug_decay, line2_color, decay_x_label, decay_title, decay_y_label)
statistics.loc['Decay Time'] = [test_stat, p_val, test]


################
# Export Plots #
################

path = input('Enter file path for figure export: ')
file_path = path.strip().strip('"').strip("'")

export_dir = Path(file_path)
export_dir.mkdir(parents=True, exist_ok=True)

groupname = input('Enter qualifier for group, i.e. SST for SST neurons, all for all cells, etc.: ')

export_path = export_dir / f"{groupname}_VC_Frequency_line.pdf"
freq_plot.savefig(export_path)
export_path = export_dir / f"{groupname}_VC_Frequency_bar.pdf"
freq_bar.savefig(export_path)

# export_path = export_dir / f"{groupname}_VC_Charge_line.pdf"
# charge_plot.savefig(export_path)
# export_path = export_dir / f"{groupname}_VC_Charge_bar.pdf"
# charge_bar.savefig(export_path)

export_path = export_dir / f"{groupname}_VC_Amplitude_line.pdf"
amplitude_plot.savefig(export_path)
export_path = export_dir / f"{groupname}_VC_Amplitude_bar.pdf"
amp_bar.savefig(export_path)

export_path = export_dir / f"{groupname}_VC_Rise_Time_line.pdf"
rise_time_plot.savefig(export_path)
export_path = export_dir / f"{groupname}_VC_Rise_Time_bar.pdf"
rise_bar.savefig(export_path)

export_path = export_dir / f"{groupname}_VC_Decay_Time_line.pdf"
decay_time_plot.savefig(export_path)
export_path = export_dir / f"{groupname}_VC_Decay_Time_bar.pdf"
decay_bar.savefig(export_path)

export_path = export_dir /f"{groupname}_stats.csv"
statistics.to_csv(export_path, index=True)