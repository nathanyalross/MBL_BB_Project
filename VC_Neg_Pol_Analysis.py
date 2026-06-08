import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from pathlib import Path
from scipy import signal
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid  # modern alternative to np.trapz
import os
import json
from datetime import datetime
from typing import List, Dict, Any

#Initiate list of paramaters for analysis
params = []

def meta_analysis(path: str, input_list, params: List[str]) -> Dict[str, Any]:
    """
    Function to create a meta-analysis json file that includes all input dataframe names and analyses ran

    Will check path to see if a meta analysis has already been done, if a json file is in the path
    it will append to that file

    Args: 
        path: path containing processed data
        input_list: name of input data
        params: analysis parameters for data

    Returns:
        meta_json: exports json file to path that includes meta data, overwriting existing one if needed
    """
    
    # Define the meta-analysis file path
    meta_file_path = os.path.join(path, "meta_analysis.json")
    
    # Initialize or load existing meta-analysis data
    # If the json file already exists:
    if os.path.exists(meta_file_path):
        #try opening the existing json file
        try:
            with open(meta_file_path, 'r') as f:
                meta_data = json.load(f)
            print(f"Loaded existing meta-analysis file from {meta_file_path}")
        #If opening fails, just report an error and create new json file
        except (json.JSONDecodeError, IOError) as e:
            print(f"Error loading existing meta-analysis file: {e}")
            print("Creating new meta-analysis file...")
            meta_data = {}
    #If there is no existing json file then a new file will be created
    else:
        meta_data = {}
        print("Creating new meta-analysis file...")
    
    # Initialize structure if it doesn't exist
    if "meta_analysis" not in meta_data:
        meta_data["meta_analysis"] = {
            "created_date": datetime.now().isoformat(),
            "last_updated": datetime.now().isoformat(),
            "sessions": []
        }
    
    # Update last_updated timestamp
    meta_data["meta_analysis"]["last_updated"] = datetime.now().isoformat()
    
    # Create current session data
    current_session = {
        #Will identify which session this is in order
        "session_id": len(meta_data["meta_analysis"]["sessions"]) + 1,
        #Identifies time that the session started
        "timestamp": datetime.now().isoformat(),
        #Lists input dataframes
        "Input Data": input_list,
        #Lists the analyses that were performed
        "Parameters Used": params
    }
    
    # Check if this exact session already exists (avoid duplicates)
    session_exists = False
    #Iterate through all of the sessions in the metadata
    for session in meta_data["meta_analysis"]["sessions"]:
        #If all of the session data matches a past session, then it will update the timestamp for identical session
        if (session.get("Input Data") == input_list and 
            session.get("Parameters Used") == params):
            session_exists = True
            print("Identical session found. Updating timestamp...")
            #Update timestamp of identical session (re-analyzed)
            session["timestamp"] = datetime.now().isoformat()
            break
    
    #If the current session does not exist then add the current session to metadata as new session
    if not session_exists:
        meta_data["meta_analysis"]["sessions"].append(current_session)
        print(f"Added new session with {len(input_list)} datasets")
    
    # Add summary statistics
    all_datasets = set()
    
    #Look through all sessions in file and find all input dataframes and sessions used
    for session in meta_data["meta_analysis"]["sessions"]:
        all_datasets.update(session.get("Input Data", []))
    
    #Create summary data
    meta_data["meta_analysis"]["summary"] = {
        "total_unique_datasets": len(all_datasets),
        "total_sessions": len(meta_data["meta_analysis"]["sessions"]),
        "all_datasets_used": sorted(list(all_datasets)),
    }
    
    # Ensure directory exists
    os.makedirs(path, exist_ok=True)
    
    # Write the updated meta-analysis file
    try:
        with open(meta_file_path, 'w') as f:
            json.dump(meta_data, f, indent=2, sort_keys=True)
        print(f"Meta-analysis file saved to: {meta_file_path}")
    except IOError as e:
        print(f"Error saving meta-analysis file: {e}")
        return None
    
    return meta_data

def explore_h5(file_path):
    """Print the data organization structure"""
    datasets = []
    with h5py.File(file_path, 'r') as f:
        def collect_datasets(name, obj):
            if isinstance(obj, h5py.Dataset):
                datasets.append(name)
                print(f"Dataset: {name}")
                print(f"  Shape: {obj.shape}")
                print(f"  Dtype: {obj.dtype}")

        print("File Structure and Datasets:")
        f.visititems(collect_datasets)

    return datasets

def apply_lowpass_filter(data, cutoff_hz=1500, sampling_rate=20000, order=4):
    """Apply low-pass Butterworth filter to data"""
    nyquist = sampling_rate / 2
    normalized_cutoff = cutoff_hz / nyquist
    b, a = signal.butter(order, normalized_cutoff, btype='low')
    return signal.filtfilt(b, a, data)

def median_filter_baseline(data, window_size=1000):
    """Calculate baseline using median filter to avoid event contamination"""
    # Ensure window_size is odd
    if window_size % 2 == 0:
        window_size += 1

    # Make sure window_size is not larger than data length
    window_size = min(window_size, len(data))
    if window_size % 2 == 0:
        window_size -= 1

    # Apply median filter to get baseline estimate
    baseline = signal.medfilt(data, kernel_size=window_size)
    return baseline

def detect_events(data, time_axis, threshold_factor, min_amplitude,
                  min_rise_time, max_rise_time, min_decay_time):
    """
    Detect using median-based detection

    Parameters:
    - data: current trace in pA
    - time_axis: time axis in seconds
    - threshold_factor: multiplier for noise-based threshold
    - min_amplitude: minimum event amplitude in pA
    - min_rise_time: minimum rise time in seconds
    - max_rise_time: maximum rise time in seconds
    - min_decay_time: minimum decay time in seconds
    """

    # Calculate baseline and noise
    baseline = median_filter_baseline(data, window_size=1001)  # Use odd number either 1001 or 501
    noise_trace = baseline - data #Inverted for negative polarity
    noise_std = np.std(noise_trace)

    # Set detection threshold (positive for negative events)
    threshold = threshold_factor * noise_std

    print(f"Noise std: {noise_std:.2f} pA")
    print(f"Detection threshold: {threshold:.2f} pA")

    # Calculate sampling rate and distances
    dt = time_axis[1] - time_axis[0] if len(time_axis) > 1 else 1 / 10000
    sampling_rate = 1 / dt

    # Convert time-based parameters to sample distances
    min_distance_samples = int(0.005 * sampling_rate)  # 5ms (0.005) instead of 10ms (0.010)
    min_width_samples = int(0.0005 * sampling_rate)  # 0.5ms (0.0005) instead of 1ms (0.001)

    # Find peaks above threshold
    peaks, properties = signal.find_peaks(
        noise_trace,
        height=threshold,
        distance=max(1, min_distance_samples),
        width=max(1, min_width_samples)
    )

    print(f"Initial peaks found: {len(peaks)}")

    events = []

    for peak_idx in peaks:
        try:
            # Calculate event properties
            event_props = analyze_event(data, time_axis, peak_idx, baseline[peak_idx])

            # Filter events based on criteria
            #if the calculated amplitude is greater than or equal to minimum amplitude
            if (event_props['amplitude'] >= min_amplitude and
                    # And rise time is between min and max requirement
                    min_rise_time <= event_props['rise_time'] <= max_rise_time and
                    #And decay time is greater than the minimum
                    event_props['decay_tau'] >= min_decay_time):
                events.append(event_props)

        except Exception as e:
            print(f"Error analyzing event at index {peak_idx}: {e}")
            continue
        
    #print(event_props)

    print(f"Events passing criteria: {len(events)}")

    return events, threshold, baseline

def analyze_event(data, time_axis, peak_idx, baseline_value):
    """Analyze event properties using working debug logic without printing."""

    # Determine time units
    # Find sampling rate - if there is no time axis defaults to 20 kHz
    dt = time_axis[1] - time_axis[0] if len(time_axis) > 1 else 1 / 10000
    # Returns T/F value if your sampling rate is greater than 20 kHz (or 0.00005)
    time_in_ms = dt > 0.001
    # Sampling rate for analysis - if time_in_ms is true then it will default to dividing sampling rate by 1000. otherwise, keep original sampling rate
    dt_seconds = dt / 1000 if time_in_ms else dt

    # Define analysis window
    #Window before peak index set to larger value either 10 samples or 0.02/sampling rate (if 20 kHz, then this is 400)
    window_pre = max(10, int(0.020 / dt_seconds))
    #Window after peak index set to larger value either 50 samples or 0.1/sr (if 20kHz then 2000)
    window_post = max(50, int(0.100 / dt_seconds))
    #Start index is peak idx-pre window as long as it isn't negative
    start_idx = max(0, peak_idx - window_pre)
    #End index same logic as above
    end_idx = min(len(data), peak_idx + window_post)

    # Baseline: use data from 5-15 ms *after* the peak
    pre_start = min(len(data), peak_idx + int(0.005 / dt_seconds))
    pre_end = min(len(data), peak_idx + int(0.015 / dt_seconds))
    event_baseline = np.mean(data[pre_start:pre_end]) if pre_end > pre_start else baseline_value

    # Peak properties: calculate amplitude as a positive value
    peak_time = time_axis[peak_idx]
    peak_amplitude = event_baseline - data[peak_idx]

    # Onset: 10% of peak, relative to baseline
    onset_threshold = event_baseline - 0.1 * peak_amplitude
    onset_idx = next((i for i in range(peak_idx, start_idx - 1, -1)
                      if data[i] >= onset_threshold), max(0, peak_idx - int(0.005 / dt_seconds)))

    # Rise time: 10% to 90%
    rise_90_threshold = event_baseline + 0.9 * peak_amplitude
    rise_90_idx = next((i for i in range(onset_idx, peak_idx + 1)
                        if data[i] >= rise_90_threshold), peak_idx)
    
    #Rise time is set as the maximum between 0.0001 (0.1 ms) or difference between 90% and 10%, whichever larger
    rise_time = max(0.0001, time_axis[rise_90_idx] - time_axis[onset_idx])

    # Decay tau
    decay_start_idx = peak_idx
    #Decay end is set as min between end of sweep and peak idx + (0.1/sampling rate, 2000 for 20 kHz)
    decay_end_idx = min(len(data), peak_idx + int(0.100 / dt_seconds))
    decay_tau = 0.020
    if decay_end_idx > decay_start_idx + 20:
        try:
            #Timeseries data between peak and decay end window
            decay_data = data[decay_start_idx:decay_end_idx] - event_baseline
            decay_data = event_baseline - data[decay_start_idx:decay_end_idx]  # Invert data for fit
            #Time of decay pulled (length of decay start index)
            decay_time = time_axis[decay_start_idx:decay_end_idx] - time_axis[decay_start_idx]
            #If time is in ms as determined earlier, then divide decay time by 1000 to get seconds
            if time_in_ms:
                decay_time = decay_time / 1000
            ### If the decay timeseries is greter than 10% of amplitude (flawed due to not taking into account baseline?)
            valid = decay_data > (0.1 * peak_amplitude)
            if np.sum(valid) > 10:
                popt, _ = curve_fit(
                    lambda t, A, tau: A * np.exp(-t / tau),
                    decay_time[valid],
                    decay_data[valid],
                    p0=[peak_amplitude, 0.020],
                    bounds=([0, 0.001], [5 * peak_amplitude, 0.200]),
                    maxfev=1000
                )
                if 0.001 <= popt[1] <= 0.200:
                    decay_tau = popt[1]
        except Exception:
            pass

    # Charge calculation
    charge_end_idx = next(
        #Finding range, either end of sweep or .08/sample rate (1600 for 20 kHz)
        (i for i in range(peak_idx, min(len(data), peak_idx + int(0.080 / dt_seconds)))
         #If the value at current index relative to baseline is less than 0.1 times absolute value of peak amplitude
         if abs(data[i] - event_baseline) <= 0.1 * abs(peak_amplitude)),
        # Then find the min between end of data and peak + 0.05/sr (1,000 for 20 kHz), selecting your charge end index
        min(len(data), peak_idx + int(0.050 / dt_seconds))
    )

    #Set values for charge, expected charge, and duration
    charge = 0.0
    expected_charge = 0.0
    duration_s = 0.0

    #If the onset index is less than end index, then proceed
    if onset_idx < charge_end_idx:
        #Charge current is set as the data values for the window with baseline subtraction
        charge_current = data[onset_idx:charge_end_idx] - event_baseline
        charge_current = event_baseline - data[onset_idx:charge_end_idx]  # Invert for charge calculation
        if peak_amplitude > 0:
            charge_current = np.maximum(charge_current, 0)
            # If greater than 0, then the charge_current is set as max between current and 0 for each index value
        else:
            charge_current = np.abs(np.minimum(charge_current, 0))
            #If less than zero than the charge current for each index is the absolute value or zero if turns positive

        #If sampling is in ms then convert time values to ms
        charge_time = time_axis[onset_idx:charge_end_idx]
        if time_in_ms:
            charge_time = charge_time / 1000

        #If there is a charge time to analyze:
        if len(charge_time) > 1:
            #Duration of charge window
            duration_s = charge_time[-1] - charge_time[0]
            # Find mean of time window
            mean_current = np.mean(charge_current)
            # find mean value times duration
            rectangular = mean_current * duration_s
            # Multiply peak amp value by the duration of window and divide by two. Why though?
            triangular = (peak_amplitude * duration_s) / 2
            #Calculate AUC using trapezoid - window values and the timespan of the window
            trapz = trapezoid(charge_current, charge_time)
            #This trapz output is the charge
            charge = trapz
            #Expected charge is simply the triangular calculation - for what purpose?
            expected_charge = triangular

            # Fallbacks
            if abs(charge - triangular) > 2 * triangular:
                charge = triangular
            elif abs(charge - rectangular) > max(0.5, 0.3 * max(abs(charge), abs(rectangular))):
                charge = rectangular

    return {
        'peak_time': peak_time,
        'peak_idx': peak_idx,
        'onset_time': time_axis[onset_idx] if onset_idx < len(time_axis) else peak_time,
        'amplitude': peak_amplitude,
        'rise_time': rise_time,
        'decay_tau': decay_tau,
        'charge': charge,
        'expected_charge': expected_charge,
        'duration': duration_s,
        'baseline': event_baseline
    }

def calculate_interevent_intervals(events):
    """Calculate interevent intervals"""
    if len(events) < 2:
        return []

    intervals = []
    for i in range(1, len(events)):
        interval = events[i]['peak_time'] - events[i - 1]['peak_time']
        intervals.append(interval)

    return intervals

def plot_VC(file_path):
    """Plot multiple sweeps from the VC dataset concatenated sequentially in time"""
    # ==============================================================================
    # MODIFY THESE VALUES FOR YOUR ANALYSIS
    # ==============================================================================
    sweep_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]  # Change this to plot different sweeps 
    # ==============================================================================

    #Following loop will dynamically find the VC_Cont dataset
    with h5py.File(file_path, 'r') as f:
        # Search for VC_Continuous dataset in Data group
        dataset_name = None
        for key in f['Data'].keys():
            if re.match(r'R\d+_S\d+_VC_Continuous', key):
                dataset_name = f'Data/{key}'
                break
    
    #dataset_name = 'Data/R2_S1_VC_Continuous'  #If you need to manually set the dataset name set it here

    with h5py.File(file_path, 'r') as f:
        data = f[dataset_name][:]

        # Check the actual number of sweeps available
        num_sweeps = data.shape[1]
        print(f"Dataset shape: {data.shape}")
        print(f"Number of sweeps available: {num_sweeps}")

        # Filter sweep_numbers to only include available sweeps
        available_sweeps = [s for s in sweep_numbers if s <= num_sweeps]
        if not available_sweeps:
            print(f"No valid sweeps found. Available sweeps: 1 to {num_sweeps}")
            return []

        print(f"Plotting sweeps: {available_sweeps}")

        # Create figure for concatenated traces
        fig, ax = plt.subplots(figsize=(15, 6))

        # Prepare arrays to concatenate all sweeps
        all_time = []
        all_current = []

        # Calculate time axis for one sweep
        num_samples = data.shape[0]
        # Based on the data shape (14000 samples), assuming 10 kHz sampling rate
        time_axis = np.arange(num_samples) * (1 / 10000)  # 10 kHz = 1/10000 s per sample
        sweep_duration_s = time_axis[-1]  # Duration of one sweep in seconds

        for i, sweep_number in enumerate(available_sweeps):
            # Convert sweep number to index (1-indexed to 0-indexed)
            sweep_index = sweep_number - 1

            # Get current data for this sweep
            current_pA = data[:, sweep_index] * 1e12  # Convert to pA
            current_pA = apply_lowpass_filter(current_pA, cutoff_hz=1500)  # Apply 5 Hz filter

            # Create time axis for this sweep, offset by previous sweeps
            time_offset = i * sweep_duration_s  # Each sweep starts after the previous one
            time_sweep = time_axis + time_offset

            # Store data
            all_time.extend(time_sweep)
            all_current.extend(current_pA)

        # Convert to numpy arrays
        all_time = np.array(all_time)
        all_current = np.array(all_current)

        # Plot the concatenated trace
        ax.plot(all_time, all_current, color='black', linewidth=0.5)
        # Define scale bar sizes
        time_scale_bar = 60  # In seconds
        current_scale_bar = 100  # In pA
        
        # Get data ranges
        time_range = all_time.max() - all_time.min()
        current_range = all_current.max() - all_current.min()
        
        # Position scale bars in bottom-right area (as percentages of data range)
        # Horizontal position: 85% from left edge
        # Vertical position: 10% from bottom edge
        scale_bar_x = all_time.min() + 0.85 * time_range
        scale_bar_y = all_current.min() + 0.10 * current_range
        
        # Draw horizontal scale bar (time)
        ax.plot([scale_bar_x, scale_bar_x + time_scale_bar], 
                [scale_bar_y, scale_bar_y], 
                color='black', linewidth=3)
        
        # Draw vertical scale bar (current)
        ax.plot([scale_bar_x, scale_bar_x], 
                [scale_bar_y, scale_bar_y + current_scale_bar], 
                color='black', linewidth=3)
        
        # Add labels centered along scale bars
        # Time label (centered below horizontal bar)
        ax.text(scale_bar_x + time_scale_bar/2, 
                scale_bar_y - 0.05 * current_range,
                '60 s', 
                ha='center', va='top', fontsize=10, color='black')
        
        # Current label (centered along vertical bar, rotated)
        ax.text(scale_bar_x - 0.02 * time_range, 
                scale_bar_y + current_scale_bar/2,
                '100 pA', 
                ha='center', va='bottom', fontsize=10, color='black', rotation=90)
        # ===================================================================

        plt.title('')
        plt.xlabel("",font='Arial',fontsize=16, fontweight='bold')
        plt.ylabel("",font='Arial',fontsize=16, fontweight='bold')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.tight_layout()
        plt.show()

        return fig, all_time, all_current

def plot_event_detection(file_path):
#New Plot Event Detection function that will plot and save information for each sweep individually
    """Plot VC data with event detection and analysis"""

    # ==============================================================================
    # MODIFY THESE VALUES FOR YOUR ANALYSIS
    # ==============================================================================
    sweep_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

    # Detection parameters - keep same as plot_event_overlay
    threshold_factor = 4.0  # Threshold = threshold_factor * noise_std
    min_amplitude = 6  # Minimum event amplitude in pA
    min_rise_time = 0.0001 # Minimum rise time in seconds (0.001 equals 1ms)
    max_rise_time = 0.01 # Maximum rise time in seconds (0.1 equals 100ms)
    min_decay_time = 0.0001 # Minimum decay time in seconds (1ms)

    params.append([f'Threshold factor:{threshold_factor}', f'Minimum amplitude:{min_amplitude}', f'Minimum Rise Time:{min_rise_time}',
                  f'Maximum Rise Time:{max_rise_time}', f'Minimum Decay Time:{min_decay_time}'])

    # Analysis window parameters
    analysis_duration = 60  # Duration of analysis window (seconds)
    # ==============================================================================

    #Following loop will process data for each sweep
    with h5py.File(file_path, 'r') as f:
        # Search for VC_Continuous dataset in Data group
        dataset_name = None
        for key in f['Data'].keys():
            if re.match(r'R\d+_S\d+_VC_Continuous', key):
                dataset_name = f'Data/{key}'
                break

        #dataset_name = 'Data/R10_S1_VC_Continuous'   #If you need to set dataset name manually do it here
        
        print(f"Using dataset: {dataset_name}")
        data = f[dataset_name][:]

        # Check the actual number of sweeps available
        num_sweeps = data.shape[1]
        print(f"Dataset shape: {data.shape}")
        print(f"Number of sweeps available: {num_sweeps}")

        # Filter sweep_numbers to only include available sweeps
        available_sweeps = [s for s in sweep_numbers if s <= num_sweeps]
        if not available_sweeps:
            print(f"No valid sweeps found. Available sweeps: 1 to {num_sweeps}")
            return []

        print(f"Processing sweeps for event detection: {available_sweeps}")

        # Prepare arrays to concatenate all sweeps
        all_time = []
        all_current = []
        all_events = []

        # Calculate time axis for one sweep
        num_samples = data.shape[0]
        sampling_rate = 10000  # 10 kHz
        time_axis = np.arange(num_samples) / sampling_rate
        sweep_duration_s = time_axis[-1]

        # Process each sweep
        for i, sweep_number in enumerate(available_sweeps):
            sweep_index = sweep_number - 1
            current_pA = data[:, sweep_index] * 1e12  # Convert to pA
            current_pA = apply_lowpass_filter(current_pA, cutoff_hz=5000)  # Apply 5 Hz filter

            # Create time axis for this sweep, offset by previous sweeps
            time_offset = i * sweep_duration_s
            time_sweep = time_axis + time_offset

            # Detect events in this sweep
            print(f"\nProcessing sweep {sweep_number}...")
            events, threshold, baseline = detect_events(
                current_pA, time_axis, threshold_factor, min_amplitude,
                min_rise_time, max_rise_time, min_decay_time
            )

            # Adjust event times for concatenated display
            for event in events:
                event['peak_time'] += time_offset
                event['onset_time'] += time_offset

            all_events.extend(events)
            all_time.extend(time_sweep)
            all_current.extend(current_pA)

        # Convert to numpy arrays
        all_time = np.array(all_time)
        all_current = np.array(all_current)

        start_times= list(range(0, int(len(all_time)/sampling_rate), analysis_duration))

        #List to store dictionaries for summary stats
        sweep_data = []
        #List to store event data for each sweep
        event_data = []

        # Bin data and iterate through on a 60 second basis
        for x, start_time in enumerate(start_times):
            # Create figure for analysis window only
            plt.figure(figsize=(15, 6))

            # Analysis window (5 seconds)
            analysis_end_time = start_time + analysis_duration

            # Find indices for analysis window
            analysis_mask = (all_time >= start_time) & (all_time <= analysis_end_time)
            analysis_time = all_time[analysis_mask]
            analysis_current = all_current[analysis_mask]

            if len(analysis_time) > 0:
                plt.plot(analysis_time, analysis_current, color='black', linewidth=1)

                # Mark events in analysis window
                analysis_events = [event for event in all_events
                                if start_time <= event['peak_time'] <= analysis_end_time]

                if analysis_events:
                    event_times = [event['peak_time'] for event in analysis_events]
                    event_amplitudes = []

                    for event_time in event_times:
                        idx = np.argmin(np.abs(analysis_time - event_time))
                        event_amplitudes.append(analysis_current[idx])

                    plt.scatter(event_times, event_amplitudes, color='red', s=50, zorder=5)

                    # Add event numbers
                    for i, (t, a) in enumerate(zip(event_times, event_amplitudes)):
                        plt.annotate(f'{i + 1}', (t, a), xytext=(5, 10),
                                    textcoords='offset points', fontsize=8, color='red')

                plt.title(
                    f'Event Detection Analysis Window: {start_time}-{analysis_end_time}s ({len(analysis_events)} events)')
                plt.xlabel("Time (s)")
                plt.ylabel("Current (pA)")
                plt.grid(False)
                plt.xlim(start_time, analysis_end_time)

                # Print event properties for analysis window
                if analysis_events:
                    print(f"\n{'=' * 80}")
                    print(f"EVENT PROPERTIES - Analysis Window ({start_time}-{analysis_end_time}s)")
                    print(f"{'=' * 80}")

                    # Calculate interevent intervals
                    intervals = calculate_interevent_intervals(analysis_events)

                    for i, event in enumerate(analysis_events):
                        # Create a dictionary for event data
                        event_dict = {
                            'Sweep Number': x + 1,
                            'Time (s)': event['peak_time'],
                            'Amplitude (pA)': event['amplitude'],
                            'Rise time (ms)': event['rise_time'],
                            'Decay Tau': event['decay_tau'],
                            'Charge (pC)': event['charge'],
                            'Charge expected (pC)': event['expected_charge'] * 1,
                            'Duration of event': event['duration'],
                        }

                        event_data.append(event_dict)

                    # Summary statistics
                    amplitudes = [event['amplitude'] for event in analysis_events]
                    rise_times = [event['rise_time'] * 1000 for event in analysis_events]
                    decay_taus = [event['decay_tau'] * 1000 for event in analysis_events if
                                not np.isnan(event['decay_tau'])]
                    charges = [event['charge'] * 1 for event in analysis_events]
                    expected_charges = [event['expected_charge'] * 1 for event in analysis_events]
                    event_dur = [event['duration'] * 1 for event in analysis_events]

                    # Create a dictionary for this sweep's data
                    sweep_dict = {
                        'Number of Events': len(analysis_events),
                        'Event frequency (Hz)': len(analysis_events) / analysis_duration,
                        'Amplitude (pA)': np.mean(amplitudes),
                        'SD of Amplitude (pA)': np.std(amplitudes),
                        'Rise time (ms)': np.mean(rise_times),
                        'SD of Rise time (ms)': np.std(rise_times),
                        'Charge (pC)': np.mean(charges),
                        'SD Charge (pC)': np.std(charges),
                        'Charge expected (pC)': np.mean(expected_charges),
                        'SD Charge expected (pC)': np.std(expected_charges),
                        'Event Durations (s)': np.mean(event_dur),
                        'SD Event Durations (s)': np.std(event_dur),
                        'SD Noise (pA)': threshold
                    }
                    
                    if decay_taus:
                        sweep_dict['Decay tau (ms)'] = np.mean(decay_taus)
                        sweep_dict['SD Decay tau (ms)'] = np.std(decay_taus)
                    
                    if intervals:
                        sweep_dict['Interevent Interval (ms)'] = np.mean(intervals) * 1000
                        sweep_dict['SD Interevent Interval (ms)'] = np.std(intervals) * 1000
                    
                    sweep_data.append(sweep_dict)
                elif len(analysis_events)==0:
                    sweep_dict = {
                        'Number of Events': 0,
                        'Event frequency (Hz)': 0,
                        'Amplitude (pA)': 'NA',
                        'SD of Amplitude (pA)': 'NA',
                        'Rise time (ms)': 'NA',
                        'SD of Rise time (ms)': 'NA',
                        'Charge (pC)': 0,
                        'SD Charge (pC)': 0,
                        'Charge expected (pC)': 0,
                        'Duration of event': 'NA',
                        'SD Noise (pA)': 'NA'
                    }
                    sweep_data.append(sweep_dict)

                    event_dict = {
                        'Sweep Number': x,
                        'Time (s)': 'NA',
                        'Amplitude (pA)': 0,
                        'Rise time (ms)': 0,
                        'Decay Tau': 0,
                        'Charge (pC)': 0,
                        'Charge expected (pC)': 0,
                        'Duration of event': 0
                    }
                    event_data.append(event_dict)

            else:
                plt.text(0.5, 0.5, 'No data in analysis window',
                        transform=plt.gca().transAxes, ha='center', va='center')
                plt.title('Analysis Window: No Data')

            plt.tight_layout()
            plt.show()

        sum_df = pd.DataFrame(sweep_data).T
        sum_df.columns = [f'Sweep {i+1}' for i in range(len(sweep_data))]
        event_df = pd.DataFrame(event_data)
        return sum_df,event_df

def plot_event_overlay(file_path):
    """Plot all detected events overlaid with average waveform using improved peak alignment"""

    # ==============================================================================
    # MODIFY THESE VALUES FOR YOUR ANALYSIS (same as plot_event_detection)
    # ==============================================================================
    sweep_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

    #Following loop will dynamically find the VC_Cont dataset
    with h5py.File(file_path, 'r') as f:
        # Search for VC_Continuous dataset in Data group
        dataset_name = None
        for key in f['Data'].keys():
            if re.match(r'R\d+_S\d+_VC_Continuous', key):
                dataset_name = f'Data/{key}'
                break

    #dataset_name = 'Data/R10_S1_VC_Continuous'   #If you need to manually set the dataset name do that here

    # Detection parameters - keep same as plot_event_detection
    threshold_factor = 4.0  # Threshold = threshold_factor * noise_std
    min_amplitude = 6  # Minimum event amplitude in pA
    min_rise_time = 0.0001 # Minimum rise time in seconds (0.001 equals 1ms)
    max_rise_time = 0.01 # Maximum rise time in seconds (0.1 equals 100ms)
    min_decay_time = 0.001 # Minimum decay time in seconds (1ms)
    max_amplitude = 2000  # Maximum event amplitude in pA
    max_decay_time = 0.20

    params.append([f'Maximum Amplitude:{max_amplitude}', f'Maximum Decay Time:{max_decay_time}'])

    # Analysis window parameters
    analysis_start_time = 0.0
    analysis_duration = 60

    # Event overlay parameters
    pre_event_time = 0.025  # 25 ms before event peak
    post_event_time = 0.300  # 300 ms after event peak

    # NEW: Peak alignment parameters
    peak_search_window = 0.0004  # 4 ms window (2 ms before and after detected peak)
    # ==============================================================================

    with h5py.File(file_path, 'r') as f:
        data = f[dataset_name][:]

        # Process data same as plot_event_detection
        available_sweeps = [s for s in sweep_numbers if s <= data.shape[1]]
        all_time = []
        all_current = []
        all_events = []

        num_samples = data.shape[0]
        sampling_rate = 10000 #10 kHz
        time_axis = np.arange(num_samples) / sampling_rate
        sweep_duration_s = time_axis[-1]

        for i, sweep_number in enumerate(available_sweeps):
            sweep_index = sweep_number - 1
            current_pA = data[:, sweep_index] * 1e12
            current_pA = apply_lowpass_filter(current_pA, cutoff_hz=9999) #Different than filter applied earlier (1500)

            time_offset = i * sweep_duration_s
            time_sweep = time_axis + time_offset

            events, threshold, baseline = detect_events(
                current_pA, time_axis, threshold_factor, min_amplitude,
                min_rise_time, max_rise_time, min_decay_time
            )

            # Filter events with max parameters and adjust timing
            filtered_events = []
            for event in events:
                if (event['amplitude'] <= max_amplitude and
                        event['decay_tau'] <= max_decay_time):
                    event['peak_time'] += time_offset
                    event['onset_time'] += time_offset
                    filtered_events.append(event)

            all_events.extend(filtered_events)
            all_time.extend(time_sweep)
            all_current.extend(current_pA)

        all_time = np.array(all_time)
        all_current = np.array(all_current)

        # Filter events in analysis window
        analysis_end_time = analysis_start_time + analysis_duration
        analysis_events = [event for event in all_events
                           if analysis_start_time <= event['peak_time'] <= analysis_end_time]

        if not analysis_events:
            print("No events found in analysis window")
            return

        # Extract individual event traces with improved peak alignment
        dt = all_time[1] - all_time[0]
        pre_samples = int(pre_event_time / dt)
        post_samples = int(post_event_time / dt)
        peak_search_samples = int(peak_search_window / (2 * dt))  # samples for ±0.5ms window

        event_traces = []
        event_time_axis = np.arange(-pre_samples, post_samples) * dt * 1000  # Convert to ms

        plt.figure(figsize=(4, 6))

        for i, event in enumerate(analysis_events):
            # Find approximate event peak index in concatenated data
            approx_peak_idx = np.argmin(np.abs(all_time - event['peak_time']))

            # Define search window around the detected peak
            search_start = max(0, approx_peak_idx - peak_search_samples)
            search_end = min(len(all_current), approx_peak_idx + peak_search_samples)

            # Cross-correlation alignment (first event becomes template)
            if i == 0:  # Store first valid event as template
                template_start = max(0, approx_peak_idx - int(0.005 / dt))  # 5ms before peak
                template_end = min(len(all_current), approx_peak_idx + int(0.015 / dt))  # 15ms after peak
                if template_end > template_start:
                    template = all_current[template_start:template_end] - event['baseline']

            # Cross-correlate current event with template
            if 'template' in locals() and len(template) > 10:
                # Extract longer segment for correlation
                corr_start = max(0, approx_peak_idx - int(0.010 / dt))
                corr_end = min(len(all_current), approx_peak_idx + int(0.020 / dt))
                corr_data = all_current[corr_start:corr_end] - event['baseline']

                if len(corr_data) >= len(template):
                    correlation = np.correlate(corr_data, template, mode='valid')
                    if len(correlation) > 0:
                        best_match_idx = np.argmax(correlation)
                        true_peak_idx = corr_start + best_match_idx + len(template) // 2
                    else:
                        true_peak_idx = approx_peak_idx
                else:
                    true_peak_idx = approx_peak_idx
            else:
                true_peak_idx = approx_peak_idx

            # Extract trace around the refined peak
            start_idx = max(0, true_peak_idx - pre_samples)
            end_idx = min(len(all_current), true_peak_idx + post_samples)

            if end_idx - start_idx == len(event_time_axis):
                event_trace = all_current[start_idx:end_idx] - event['baseline']
                event_traces.append(event_trace)

                # Plot individual event in gray
                plt.plot(event_time_axis, event_trace, color='gray', alpha=0.6, linewidth=0.6)

        # Calculate and plot average waveform
        if event_traces:
            avg_trace = np.mean(event_traces, axis=0)
            plt.plot(event_time_axis, avg_trace, color='black', linewidth=1.2, label=f'Average (n={len(event_traces)})')

            #plt.axvline(x=0, color='red', linestyle='--', alpha=0.5, label='Peak')
            plt.xlabel('Time relative to peak (ms)')
            plt.ylabel('Current (pA)')
            plt.title(f'Event Overlay - {len(event_traces)} events (Improved Alignment)')
            plt.legend()
            plt.grid(False)
            plt.tight_layout()
            plt.show()

            print(f"Plotted {len(event_traces)} events with improved peak alignment")
            print(f"Used {peak_search_window * 1000:.1f} ms search window for peak refinement")

        else:
            print("No valid event traces extracted")
        
# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

#Altered approach to manually enter path
path = input('Enter file path for analysis: ')
file_path = path.strip("'")
path_object= Path(file_path)
filename = path_object.stem

export_path = path_object.parent
export_path = str(export_path) + '/Exports'
export_dir = Path(export_path)
export_dir.mkdir(parents=True, exist_ok=True)

# Explore file structure
print("Exploring file structure...")
datasets = explore_h5(file_path)

print(f"\n{'=' * 60}")
print("PLOTTING ORIGINAL VC SWEEP")
print(f"{'=' * 60}")

# Plot the original concatenated trace
trace, time_data, current_data = plot_VC(file_path)
export_path = export_dir / f"{filename}_VC_trace.pdf"
trace.savefig(export_path)

print(f"\n{'=' * 60}")
print("PLOTTING EVENT DETECTION ANALYSIS")
print(f"{'=' * 60}")

# Plot event detection analysis and export summary statistics
sum_data, event_data = plot_event_detection(file_path)
export_path = export_dir / f"{filename}_summary.csv"
sum_data.to_csv(export_path, index=True)
export_path = export_dir / f"{filename}.csv"
event_data.to_csv(export_path, index=True)

print(f"\n{'=' * 60}")
print("PLOTTING EVENT OVERLAY")
print(f"{'=' * 60}")

# Plot event overlay
plot_event_overlay(file_path)

#Export metadata
meta_analysis(export_dir, str(filename), params)