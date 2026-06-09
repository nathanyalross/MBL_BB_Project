import matplotlib.pyplot as plt
import h5py
import numpy as np
import csv
import os
import re
from pathlib import Path
import os

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

def plot_VC(file_path):
    """Plot specific sweeps from the VC2 dataset and perform measurements"""
    # ==============================================================================
    # MODIFY THESE VALUES FOR YOUR ANALYSIS
    # ==============================================================================

    #Manually specify which run to analyze in command prompt
    run = input('Specify Initial IV Curve Run, number only: ')

    with h5py.File(file_path, 'r') as f:
        # Search for VC_Continuous dataset in Data group
        dataset_name = None
        for key in f['Data'].keys():
            if re.match(f'R{run}_S1_VC_IV_fast', key):
            #if re.match(f'R{run}+_S\d+_VC_IV_Marcela', key):
                dataset_name = f'Data/{key}'
                print(dataset_name)
                break

    #If you need to manually enter dataset name set that here
    #dataset_name = 'Data/R4_S1_VC_IV_Marcela'  #If you need to manually set the dataset name set it here

    sweep_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13]  # Change this to plot different sweeps
    # ==============================================================================

    with h5py.File(file_path, 'r') as f:
        data = f[dataset_name][:]

        # Store results for all sweeps
        all_results = []

        for sweep_number in sweep_numbers:
            # Convert sweep number to index (1-indexed to 0-indexed)
            sweep_index = sweep_number - 1

            # Create time axis (20 kHz = 50 μs per sample)
            time_axis = np.arange(12000) * 50e-6  # 50 μs per sample
            time_ms = time_axis * 1000  # Convert to milliseconds 
            current_pA = data[:, sweep_index] * 1e12  # Convert to pA

            # ==============================================================================
            # MEASUREMENTS
            # ==============================================================================

            # 1. Minimum value between 50-51 ms
            mask_50_51 = (time_ms >= 50) & (time_ms <= 51)
            min_value_50_51 = np.min(current_pA[mask_50_51])

            # 2. Baseline (first 50 ms) vs 140-145 ms difference
            mask_baseline = time_ms <= 50
            mask_140_145 = (time_ms >= 140) & (time_ms <= 145)
            baseline_avg = np.mean(current_pA[mask_baseline])
            avg_140_145 = np.mean(current_pA[mask_140_145])
            difference_baseline_140_145 = avg_140_145 - baseline_avg

            # 3. AUC calculation using 50-80 ms average level as reference
            # Use horizontal line at 50-80 ms average level
            mask_50_80 = (time_ms >= 50) & (time_ms <= 80)
            time_projection_50_80 = time_ms[mask_50_80]

            # Horizontal reference line at 140-145 ms average level
            projected_line_50_80 = np.full_like(time_projection_50_80, avg_140_145)

            # Calculate AUC between actual data and 140-145 ms average level (50-80 ms)
            actual_data_projection_50_80 = current_pA[mask_50_80]

            # Only calculate AUC for values BELOW the 140-145 ms average level
            difference_from_140_145 = actual_data_projection_50_80 - projected_line_50_80
            below_140_145_mask = difference_from_140_145 < 0
            difference_from_140_145[~below_140_145_mask] = 0  # Set values above 140-145 ms level to 0

            auc_projection_50_80 = np.trapezoid(difference_from_140_145, time_projection_50_80)


            # 4. AUC calculation using baseline (0-50 ms) average level as reference
            # Use horizontal line at baseline average level
            mask_151_249 = (time_ms >= 151) & (time_ms <= 249)
            time_projection_151_249 = time_ms[mask_151_249]

            # Horizontal reference line at baseline average level
            projected_line_151_249 = np.full_like(time_projection_151_249, baseline_avg)

            # Calculate AUC between actual data and baseline average level (150-250 ms)
            actual_data_projection_151_249 = current_pA[mask_151_249]

            # Only calculate AUC for values ABOVE the baseline
            difference_from_baseline = actual_data_projection_151_249 - projected_line_151_249
            above_baseline_mask = difference_from_baseline > 0
            difference_from_baseline[~above_baseline_mask] = 0  # Set values below baseline to 0

            auc_projection_151_249 = np.trapezoid(difference_from_baseline, time_projection_151_249)

            # Store results for this sweep
            sweep_results = {
                'sweep_number': sweep_number,
                'min_value_50_51': min_value_50_51,
                'baseline_avg': baseline_avg,
                'avg_140_145': avg_140_145,
                'difference_baseline_140_145': difference_baseline_140_145,
                'auc_projection_50_80': auc_projection_50_80,
                'auc_projection_150_250': auc_projection_151_249
            }
            all_results.append(sweep_results)

            # ==============================================================================
            # PRINT MEASUREMENTS FOR THIS SWEEP
            # ==============================================================================
            print(f"\n{'=' * 50}")
            print(f"MEASUREMENTS FOR VC2 SWEEP {sweep_number}")
            print(f"{'=' * 50}")
            print(f"1. Minimum value between 50-51 ms: {min_value_50_51:.2f} pA")
            print(f"2. Baseline average (0-50 ms): {baseline_avg:.2f} pA")
            print(f"   Average 140-145 ms: {avg_140_145:.2f} pA")
            print(f"   Difference (140-145 ms - baseline): {difference_baseline_140_145:.2f} pA")
            print(f"3. AUC relative to 140-145 ms level (50-80 ms): {auc_projection_50_80:.2f} pA·ms")
            print(f"4. AUC relative to baseline level (150-246 ms): {auc_projection_151_249:.2f} pA·ms")
            print(f"{'=' * 50}")

            # Plot the sweep with measurement annotations
            plt.figure(figsize=(6,6))
            plt.plot(time_ms, current_pA, 'black', linewidth=0.5, label='Data')

            # Add shaded regions only for the measured areas (not full Y-axis)
            y_min, y_max = plt.ylim()

            # Shade the specific measurement regions
            plt.fill_between(time_ms, current_pA, 0, where=mask_50_51, alpha=0.3, color='red',
                             label='Min region (50-51 ms)')
            plt.fill_between(time_ms, current_pA, baseline_avg, where=mask_baseline, alpha=0.3, color='blue',
                             label='Baseline (0-50 ms)')
            plt.fill_between(time_ms, current_pA, avg_140_145, where=mask_140_145, alpha=0.3, color='green',
                             label='End region (140-145 ms)')

            # Add AUC shading relative to 140-145 ms average level (50-80 ms) - ONLY BELOW
            below_140_145_in_region = (current_pA < avg_140_145) & mask_50_80
            plt.fill_between(time_ms, current_pA, avg_140_145, where=below_140_145_in_region, alpha=0.3, color='purple',
                             label='AUC relative to 140-145 ms level (50-80 ms) - BELOW only')

            # Add AUC shading relative to baseline average level (150-250 ms) - ONLY ABOVE baseline
            above_baseline_in_region = (current_pA > baseline_avg) & mask_151_249
            plt.fill_between(time_ms, current_pA, baseline_avg, where=above_baseline_in_region, alpha=0.3,
                             color='orange',
                             label='AUC relative to baseline level (150-250 ms) - ABOVE only')

            # Add horizontal reference lines for averages
            plt.axhline(y=baseline_avg, color='blue', linestyle=':', alpha=0.7)
            plt.axhline(y=avg_140_145, color='green', linestyle=':', alpha=0.7)

            plt.title(f'VC BEGIN - Sweep {sweep_number}')
            plt.xlabel("Time (ms)")
            plt.ylabel("Current (pA)")
            plt.legend()
            plt.grid(False)
            plt.tight_layout()
            #Option to show plot
            #plt.show()

        # ==============================================================================
        # SUMMARY TABLE FOR ALL SWEEPS
        # ==============================================================================
        print(f"\n{'=' * 95}")
        print("SUMMARY TABLE FOR ALL VC2 SWEEPS")
        print(f"{'=' * 95}")
        print(f"{'Sweep':<6} {'Min 50-51ms':<12} {'Baseline':<12} {'Avg 140-145ms':<14} {'Difference':<12} {'AUC 50-80ms':<12} {'AUC 150-250ms':<14}")
        print(f"{'No.':<6} {'(pA)':<12} {'(pA)':<12} {'(pA)':<14} {'(pA)':<12} {'(pA·ms)':<12} {'(pA·ms)':<14}")
        print(f"{'-' * 95}")

        for result in all_results:
            print(f"{result['sweep_number']:<6} {result['min_value_50_51']:<12.2f} {result['baseline_avg']:<12.2f} {result['avg_140_145']:<14.2f} {result['difference_baseline_140_145']:<12.2f} {result['auc_projection_50_80']:<12.2f} {result['auc_projection_150_250']:<14.2f}")

        print(f"{'=' * 95}")

        return all_results

def plot_VC2(file_path):
    """Plot specific sweeps from the VC2 dataset and perform measurements"""
    # ==============================================================================
    # MODIFY THESE VALUES FOR YOUR ANALYSIS
    # ==============================================================================

    #Manually specify which run to analyze in command prompt
    run = input('Specify Final IV Curve Run, number only: ')

    with h5py.File(file_path, 'r') as f:
        # Search for VC_Continuous dataset in Data group
        dataset_name = None
        for key in f['Data'].keys():
            if re.match(f'R{run}_S1_VC_IV_fast', key):
            #if re.match(f'R{run}+_S\d+_VC_IV_Marcela', key):
                dataset_name = f'Data/{key}'
                break

    #dataset_name = 'Data/R9_S2_VC_IV_Marcela'  #If you need to manually set the dataset name set it here

    sweep_numbers = [1,2,3,4,5,6,7,8,9,10,11,12,13]  # Change this to plot different sweeps
    # ==============================================================================

    with h5py.File(file_path, 'r') as f:
        data = f[dataset_name][:]

        # Store results for all sweeps
        all_results = []

        for sweep_number in sweep_numbers:
            # Convert sweep number to index (1-indexed to 0-indexed)
            sweep_index = sweep_number - 1

            # Create time axis (20 kHz = 50 μs per sample)
            time_axis = np.arange(12000) * 50e-6  # 50 μs per sample
            time_ms = time_axis * 1000  # Convert to milliseconds
            current_pA = data[:, sweep_index] * 1e12  # Convert to pA

            # ==============================================================================
            # MEASUREMENTS
            # ==============================================================================

            # 1. Minimum value between 50-51 ms
            mask_50_51 = (time_ms >= 50) & (time_ms <= 51)
            min_value_50_51 = np.min(current_pA[mask_50_51])

            # 2. Baseline (first 50 ms) vs 140-145 ms difference
            mask_baseline = time_ms <= 50
            mask_140_145 = (time_ms >= 140) & (time_ms <= 145)
            baseline_avg = np.mean(current_pA[mask_baseline])
            avg_140_145 = np.mean(current_pA[mask_140_145])
            difference_baseline_140_145 = avg_140_145 - baseline_avg

            # 3. AUC calculation using 140-145 ms average level as reference
            # Use horizontal line at 140-145 ms average level
            mask_50_80 = (time_ms >= 50) & (time_ms <= 80)
            time_projection_50_80 = time_ms[mask_50_80]

            # Horizontal reference line at 140-145 ms average level
            projected_line_50_80 = np.full_like(time_projection_50_80, avg_140_145)

            # Calculate AUC between actual data and 140-145 ms average level (50-80 ms)
            actual_data_projection_50_80 = current_pA[mask_50_80]

            # Only calculate AUC for values BELOW the 140-145 ms average level
            difference_from_140_145 = actual_data_projection_50_80 - projected_line_50_80
            below_140_145_mask = difference_from_140_145 < 0
            difference_from_140_145[~below_140_145_mask] = 0  # Set values above 140-145 ms level to 0

            auc_projection_50_80 = np.trapezoid(difference_from_140_145, time_projection_50_80)

            # 4. AUC calculation using baseline (0-50 ms) average level as reference
            # Use horizontal line at baseline average level
            mask_151_249 = (time_ms >= 151) & (time_ms <= 249)
            time_projection_151_249 = time_ms[mask_151_249]

            # Horizontal reference line at baseline average level
            projected_line_151_249 = np.full_like(time_projection_151_249, baseline_avg)

            # Calculate AUC between actual data and baseline average level (150-250 ms)
            actual_data_projection_151_249 = current_pA[mask_151_249]

            # Only calculate AUC for values ABOVE the baseline
            difference_from_baseline = actual_data_projection_151_249 - projected_line_151_249
            above_baseline_mask = difference_from_baseline > 0
            difference_from_baseline[~above_baseline_mask] = 0  # Set values below baseline to 0

            auc_projection_151_249 = np.trapezoid(difference_from_baseline, time_projection_151_249)

            # Store results for this sweep
            sweep_results = {
                'sweep_number': sweep_number,
                'min_value_50_51': min_value_50_51,
                'baseline_avg': baseline_avg,
                'avg_140_145': avg_140_145,
                'difference_baseline_140_145': difference_baseline_140_145,
                'auc_projection_50_80': auc_projection_50_80,
                'auc_projection_150_250': auc_projection_151_249
            }
            all_results.append(sweep_results)

            # ==============================================================================
            # PRINT MEASUREMENTS FOR THIS SWEEP
            # ==============================================================================
            print(f"\n{'=' * 50}")
            print(f"MEASUREMENTS FOR VC2 SWEEP {sweep_number}")
            print(f"{'=' * 50}")
            print(f"1. Minimum value between 50-51 ms: {min_value_50_51:.2f} pA")
            print(f"2. Baseline average (0-50 ms): {baseline_avg:.2f} pA")
            print(f"   Average 140-145 ms: {avg_140_145:.2f} pA")
            print(f"   Difference (140-145 ms - baseline): {difference_baseline_140_145:.2f} pA")
            print(f"3. AUC relative to 140-145 ms level (50-80 ms): {auc_projection_50_80:.2f} pA·ms")
            print(f"4. AUC relative to baseline level (150-246 ms): {auc_projection_151_249:.2f} pA·ms")
            print(f"{'=' * 50}")

            # Plot the sweep with measurement annotations
            plt.figure(figsize=(6, 6))
            plt.plot(time_ms, current_pA, 'black', linewidth=0.5, label='Data')

            # Add shaded regions only for the measured areas (not full Y-axis)
            y_min, y_max = plt.ylim()

            # Shade the specific measurement regions
            plt.fill_between(time_ms, current_pA, 0, where=mask_50_51, alpha=0.3, color='red',
                             label='Min region (50-51 ms)')
            plt.fill_between(time_ms, current_pA, baseline_avg, where=mask_baseline, alpha=0.3, color='blue',
                             label='Baseline (0-50 ms)')
            plt.fill_between(time_ms, current_pA, avg_140_145, where=mask_140_145, alpha=0.3, color='green',
                             label='End region (140-145 ms)')

            # Add AUC shading relative to 140-145 ms average level (50-80 ms) - ONLY BELOW
            below_140_145_in_region = (current_pA < avg_140_145) & mask_50_80
            plt.fill_between(time_ms, current_pA, avg_140_145, where=below_140_145_in_region, alpha=0.3, color='purple',
                             label='AUC relative to 140-145 ms level (50-80 ms) - BELOW only')

            # Add AUC shading relative to baseline average level (150-250 ms) - ONLY ABOVE baseline
            above_baseline_in_region = (current_pA > baseline_avg) & mask_151_249
            plt.fill_between(time_ms, current_pA, baseline_avg, where=above_baseline_in_region, alpha=0.3,
                             color='orange',
                             label='AUC relative to baseline level (150-250 ms) - ABOVE only')

            # Add horizontal reference lines for averages
            plt.axhline(y=baseline_avg, color='blue', linestyle=':', alpha=0.7)
            plt.axhline(y=avg_140_145, color='green', linestyle=':', alpha=0.7)

            plt.title(f'VC END - Sweep {sweep_number}')
            plt.xlabel("Time (ms)")
            plt.ylabel("Current (pA)")
            plt.legend()
            plt.grid(False)
            plt.tight_layout()
            #Option to show plot
            #plt.show()

        # ==============================================================================
        # SUMMARY TABLE FOR ALL SWEEPS
        # ==============================================================================
        print(f"\n{'=' * 95}")
        print("SUMMARY TABLE FOR ALL VC2 SWEEPS")
        print(f"{'=' * 95}")
        print(f"{'Sweep':<6} {'Min 50-51ms':<12} {'Baseline':<12} {'Avg 140-145ms':<14} {'Difference':<12} {'AUC 50-80ms':<12} {'AUC 150-250ms':<14}")
        print(f"{'No.':<6} {'(pA)':<12} {'(pA)':<12} {'(pA)':<14} {'(pA)':<12} {'(pA·ms)':<12} {'(pA·ms)':<14}")
        print(f"{'-' * 95}")

        for result in all_results:
            print(f"{result['sweep_number']:<6} {result['min_value_50_51']:<12.2f} {result['baseline_avg']:<12.2f} {result['avg_140_145']:<14.2f} {result['difference_baseline_140_145']:<12.2f} {result['auc_projection_50_80']:<12.2f} {result['auc_projection_150_250']:<14.2f}")

        print(f"{'=' * 95}")

        return all_results

def save_results_to_csv(file_path, vc_results, vc2_results, output_filename=None, file_location=None):
    """Save all results to a CSV file with horizontal structure"""

    # Generate output filename if not provided
    if output_filename is None:
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_filename = f"{base_name}_results.csv"

    # Add file location if provided
    if file_location is not None:
        # Create directory if it doesn't exist
        os.makedirs(file_location, exist_ok=True)
        output_filename = os.path.join(file_location, output_filename)

    # Get the base filename for the header
    file_name = os.path.splitext(os.path.basename(file_path))[0]

    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)

        # Create the complete header row with everything
        header_row = [file_name]

        # Add VC BEGIN section
        header_row.extend([
            'VC BEGIN',
            'Sweep No.',
            'Min 50-51ms (pA)',
            'Baseline (pA)',
            'Avg 140-145ms (pA)',
            'Difference (pA)',
            'AUC 50-80ms (pA·ms)',
            'AUC 150-250ms (pA·ms)'
        ])

        # Add VC END section
        header_row.extend([
            'VC END',
            'Sweep No.',
            'Min 50-51ms (pA)',
            'Baseline (pA)',
            'Avg 140-145ms (pA)',
            'Difference (pA)',
            'AUC 50-80ms (pA·ms)',
            'AUC 150-250ms (pA·ms)'
        ])

        # Add new columns before CC section
        header_row.extend([
            'Access Resistance B (Mohms)',
            'Cell Resistance B (Mohms)',
            'Capacitance B (pF)',
            'Access Resistance E (Mohms)',
            'Cell Resistance E (Mohms)',
            'Capacitance E (pF)'
        ])

        writer.writerow(header_row)

        # Find the maximum number of sweeps across all datasets
        max_sweeps = max(len(vc_results), len(vc2_results))

        # Write data rows (starting from row 2)
        for i in range(max_sweeps):
            data_row = ['']  # Empty first column

            # VC BEGIN data
            data_row.append('')  # Empty for section header column
            if i < len(vc_results):
                result = vc_results[i]
                data_row.extend([
                    result['sweep_number'],
                    f"{result['min_value_50_51']:.2f}",
                    f"{result['baseline_avg']:.2f}",
                    f"{result['avg_140_145']:.2f}",
                    f"{result['difference_baseline_140_145']:.2f}",
                    f"{result['auc_projection_50_80']:.2f}",
                    f"{result['auc_projection_150_250']:.2f}"
                ])
            else:
                data_row.extend(['', '', '', '', '', '', ''])  # Empty columns

            # VC END data
            data_row.append('')  # Empty for section header column
            if i < len(vc2_results):
                result = vc2_results[i]
                data_row.extend([
                    result['sweep_number'],
                    f"{result['min_value_50_51']:.2f}",
                    f"{result['baseline_avg']:.2f}",
                    f"{result['avg_140_145']:.2f}",
                    f"{result['difference_baseline_140_145']:.2f}",
                    f"{result['auc_projection_50_80']:.2f}",
                    f"{result['auc_projection_150_250']:.2f}"
                ])
            else:
                data_row.extend(['', '', '', '', '', '', ''])  # Empty columns

            if i < len(vc_results):
                # Access resistance calculation: -20/D_value*1000 for first row, -10/D_value*1000 for second, etc.
                d_value = vc_results[i]['min_value_50_51']  # This is column D
                if i == 0:
                    access_resistance = -20 / d_value * 1000
                elif i == 1:
                    access_resistance = -10 / d_value * 1000
                elif i == 2:
                    access_resistance = 0
                elif i == 3:
                    access_resistance = 10 / d_value * 1000
                elif i == 4:
                    access_resistance = 20 / d_value * 1000
                else:
                    access_resistance = ''  # Or continue pattern as needed

                data_row.extend([
                    f"{access_resistance:.2f}" if isinstance(access_resistance, (int, float)) else access_resistance,
                ])
            else:
                data_row.extend([''])  # Empty columns

            if i < len(vc_results):
                # Access resistance calculation: -20/E_value*1000 for first row, -10/E_value*1000 for second, etc.
                g_value = vc_results[i]['difference_baseline_140_145']  # This is column G
                if i == 0:
                    cell_resistance = -20 / g_value * 1000
                elif i == 1:
                    cell_resistance = -10 / g_value * 1000
                elif i == 2:
                    cell_resistance = 0
                elif i == 3:
                    cell_resistance = 10 / g_value * 1000
                elif i == 4:
                    cell_resistance = 20 / g_value * 1000
                else:
                    cell_resistance = ''  # Or continue pattern as needed

                data_row.extend([
                    f"{cell_resistance:.2f}" if isinstance(cell_resistance, (int, float)) else cell_resistance,
                ])
            else:
                data_row.extend([''])  # Empty columns


            if i < len(vc_results):
                # Access resistance calculation: -20/E_value*1000 for first row, -10/E_value*1000 for second, etc.
                h_value = vc_results[i]['auc_projection_50_80']  # This is column h
                if i == 0:
                    capacitance = h_value / -20
                elif i == 1:
                    capacitance = h_value / -10
                elif i == 2:
                    capacitance = 0
                elif i == 3:
                    capacitance = h_value / 10
                elif i == 4:
                    capacitance = h_value / 20
                else:
                    capacitance = ''  # Or continue pattern as needed

                data_row.extend([
                    f"{capacitance:.2f}" if isinstance(capacitance,
                                                           (int, float)) else capacitance,
                ])
            else:
                data_row.extend([''])  # Empty columns

            if i < len(vc2_results):
                # Access resistance calculation: -20/D_value*1000 for first row, -10/D_value*1000 for second, etc.
                d_value = vc2_results[i]['min_value_50_51']  # This is column D
                if i == 0:
                    access_resistance = -20 / d_value * 1000
                elif i == 1:
                    access_resistance = -10 / d_value * 1000
                elif i == 2:
                    access_resistance = 0
                elif i == 3:
                    access_resistance = 10 / d_value * 1000
                elif i == 4:
                    access_resistance = 20 / d_value * 1000
                else:
                    access_resistance = ''  # Or continue pattern as needed

                data_row.extend([
                    f"{access_resistance:.2f}" if isinstance(access_resistance, (int, float)) else access_resistance,
                ])
            else:
                data_row.extend([''])  # Empty columns

            # Add new column data before CC section
            if i < len(vc2_results):
                # Access resistance calculation: -20/E_value*1000 for first row, -10/E_value*1000 for second, etc.
                g_value = vc2_results[i]['difference_baseline_140_145']  # This is column G
                if i == 0:
                    cell_resistance = -20 / g_value * 1000
                elif i == 1:
                    cell_resistance = -10 / g_value * 1000
                elif i == 2:
                    cell_resistance = 0
                elif i == 3:
                    cell_resistance = 10 / g_value * 1000
                elif i == 4:
                    cell_resistance = 20 / g_value * 1000
                else:
                    cell_resistance = ''  # Or continue pattern as needed

                data_row.extend([
                    f"{cell_resistance:.2f}" if isinstance(cell_resistance, (int, float)) else cell_resistance,
                ])
            else:
                data_row.extend([''])  # Empty columns


            # Add new column data before CC section
            if i < len(vc2_results):
                # Access resistance calculation: -20/E_value*1000 for first row, -10/E_value*1000 for second, etc.
                h_value = vc2_results[i]['auc_projection_50_80']  # This is column h
                if i == 0:
                    capacitance = h_value / -20
                elif i == 1:
                    capacitance = h_value / -10
                elif i == 2:
                    capacitance = 0
                elif i == 3:
                    capacitance = h_value / 10
                elif i == 4:
                    capacitance = h_value / 20
                else:
                    capacitance = ''  # Or continue pattern as needed

                data_row.extend([
                    f"{capacitance:.2f}" if isinstance(capacitance,
                                                           (int, float)) else capacitance,
                ])
            else:
                data_row.extend([''])  # Empty columns

            writer.writerow(data_row)

        # Add some empty rows before action potential times
        writer.writerow([])
        writer.writerow([])
        writer.writerow([])
        writer.writerow([])
        writer.writerow([])

        # Action Potential Times Section
        writer.writerow(['Action Potential Times (s)'])
        writer.writerow(['Sweep No.', 'AP Times'])

    print(f"\nResults saved to: {output_filename}")
    return output_filename

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

path = input('Enter file path for analysis: ')
file_path = path.strip().strip('"').strip("'")

path_object = Path(file_path)
filename = path_object.stem

export_dir = path_object.parent / "Exports"
export_dir.mkdir(parents=True, exist_ok=True)

# Explore file structure
print("Exploring file structure...")
datasets = explore_h5(file_path)

print(f"\n{'=' * 60}")
print("PLOTTING VC SWEEP")
print(f"{'=' * 60}")

# Plot the specified sweep from VC dataset and store results
vc_results = plot_VC(file_path)

print(f"\n{'=' * 60}")
print("PLOTTING VC2 SWEEP")
print(f"{'=' * 60}")

# Plot the specified sweep from VC2 dataset and store results
vc2_results = plot_VC2(file_path)

print(f"\n{'=' * 60}")
print("PLOTTING CC SWEEP")
print(f"{'=' * 60}")


# Save with custom filename and location
base_filename = os.path.splitext(os.path.basename(file_path))[0]
output_name = f"{base_filename}_char.csv"
save_results_to_csv(file_path, vc_results, vc2_results,
                    output_filename=output_name, file_location=export_dir)
