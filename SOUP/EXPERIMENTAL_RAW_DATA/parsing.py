import os
import pandas as pd
import numpy as np
import re
from pathlib import Path

def parse_filename(filename):
    # Purpose: extracts the pH and starting Fe/CN concentrations encoded in a raw data filename.
    # Inputs are: filename, a raw data CSV filename, e.g. 'A1_0p6mMCN_0p02mMFe_pH6.csv'.
    # Outputs are: a dictionary {'pH': int or None, 'Fe_mM': float or None, 'CN_mM': float or None}, with each value parsed out of the filename (a value that can't be found is None).
    # Called by: process_all_data(), in this file.
    # The .csv extension is removed here.
    name = filename.replace('.csv', '')

    ph_match = re.search(r'pH(\d+)', name)
    pH = int(ph_match.group(1)) if ph_match else None

    # The Fe concentration is extracted here, e.g. '0p02mMFe' -> 0.02.
    fe_match = re.search(r'(\d+p\d+)mMFe', name)
    if fe_match:
        Fe_mM = float(fe_match.group(1).replace('p', '.'))
    else:
        Fe_mM = None

    # The CN concentration is extracted here, e.g. '0p6mMCN' -> 0.6.
    cn_match = re.search(r'(\d+p\d+)mMCN', name)
    if cn_match:
        CN_mM = float(cn_match.group(1).replace('p', '.'))
    else:
        CN_mM = None

    return {'pH': pH, 'Fe_mM': Fe_mM, 'CN_mM': CN_mM}

def absorbance_to_ferrocyanide(absorbance):
    # Purpose: converts a measured absorbance reading to a ferrocyanide concentration, using the calibration curve.
    # Inputs are: absorbance, the measured absorbance at 217 nm.
    # Outputs are: the corresponding ferrocyanide concentration, in mM, using the calibration A_217nm = 22.6876 * [Fe(CN)6^4-] + 0.15298 (rearranged to solve for concentration).
    # Called by: process_all_data(), in this file.
    return (absorbance - 0.15298) / 22.6876

def process_all_data(base_path):
    # Purpose: walks every raw CSV file across both experiment folders and combines them into one dataset with metadata attached.
    # Inputs are: base_path, the folder containing the two experiment-type subfolders (CN_conc_pH_dependence and Fe_conc_pH_dependence), each with ROUND1/ROUND2/ROUND3 subfolders of raw CSV files.
    # Outputs are: a single combined DataFrame with one row per (time, absorbance) reading across every file found, with metadata columns (pH, starting Fe/CN concentrations, round, experiment type, filename) attached - or None if no data was found.
    # Called by: directly from the entry point at the bottom of this file (if __name__ == "__main__").
    all_data = []

    folders = ['CN_conc_pH_dependence', 'Fe_conc_pH_dependence']

    for folder in folders:
        folder_path = Path(base_path) / folder

        if not folder_path.exists():
            print(f"Warning: {folder_path} does not exist, skipping...")
            continue

        for round_num in [1, 2, 3]:
            round_folder = folder_path / f'ROUND{round_num}'

            if not round_folder.exists():
                print(f"Warning: {round_folder} does not exist, skipping...")
                continue

            for csv_file in round_folder.glob('*.csv'):
                try:
                    params = parse_filename(csv_file.name)

                    df = pd.read_csv(csv_file, header=None, names=['time_min', 'absorbance'])

                    df['time_s'] = df['time_min'] * 60

                    # The absorbance is converted to ferrocyanide concentration here, first in mM, then converted to M.
                    FeCN6_mM = absorbance_to_ferrocyanide(df['absorbance'])
                    df['FeCN6_M'] = FeCN6_mM / 1000

                    # The metadata columns are added here (converting mM to M).
                    df['pH'] = params['pH']
                    df['Fe_init_M'] = params['Fe_mM'] / 1000 if params['Fe_mM'] is not None else None
                    df['HCN_init_M'] = params['CN_mM'] / 1000 if params['CN_mM'] is not None else None
                    df['round'] = round_num
                    df['experiment_type'] = folder
                    df['filename'] = csv_file.name

                    df = df[['round', 'experiment_type', 'filename', 'pH',
                            'Fe_init_M', 'HCN_init_M', 'time_s', 'FeCN6_M']]

                    all_data.append(df)
                    print(f"Processed: {csv_file.name}")

                except Exception as e:
                    print(f"Error processing {csv_file.name}: {str(e)}")

    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        return combined_df
    else:
        print("No data was processed!")
        return None

def average_kinetics_data(combined_df, output_file, target_pH=9):
    # Purpose: filters the combined dataset to one pH, removes the zero-HCN controls, and averages across repeated rounds to produce mean/std/SEM per condition and time point.
    # Inputs are: combined_df, the combined DataFrame produced by process_all_data(); output_file, the path to write the averaged CSV file to; target_pH, which pH to filter for (default 9).
    # Outputs are: a DataFrame of the averaged data (also written to output_file), or None if there is no data at target_pH, or none left once the zero-HCN controls are removed. This filters combined_df down to target_pH, drops the zero-HCN control experiments, then groups by experimental condition and exact time point, averaging across the repeated rounds to get a mean, standard deviation, and count for each group.
    # Called by: directly from the entry point at the bottom of this file (if __name__ == "__main__").

    print(f"Total rows in combined dataset: {len(combined_df)}")
    print(f"pH values present: {sorted(combined_df['pH'].unique())}")
    print(f"Rounds present: {sorted(combined_df['round'].unique())}")

    df_filtered = combined_df[combined_df['pH'] == target_pH].copy()
    print(f"\nRows after filtering to pH {target_pH}: {len(df_filtered)}")

    if len(df_filtered) == 0:
        print(f"ERROR: No data found for pH {target_pH}!")
        return None

    # The zero-HCN control experiments are removed here.
    df_filtered = df_filtered[df_filtered['HCN_init_M'] > 0].copy()
    print(f"Rows after removing 0 HCN controls: {len(df_filtered)}")

    if len(df_filtered) == 0:
        print(f"ERROR: No data remaining after removing 0 HCN!")
        return None

    print(f"\nExperiment types at pH {target_pH}:")
    print(df_filtered['experiment_type'].value_counts())

    print(f"\nUnique Fe concentrations (M): {sorted(df_filtered['Fe_init_M'].unique())}")
    print(f"Unique HCN concentrations (M): {sorted(df_filtered['HCN_init_M'].unique())}")

    # Every experimental condition and exact time point is grouped here, averaging across 'round' for the repeated experiments.
    grouping_cols = ['experiment_type', 'pH', 'Fe_init_M', 'HCN_init_M', 'time_s']

    print("\nCalculating averages across rounds...")

    averaged = df_filtered.groupby(grouping_cols).agg({
        'FeCN6_M': ['mean', 'std', 'count']
    }).reset_index()

    averaged.columns = ['experiment_type', 'pH', 'Fe_init_M', 'HCN_init_M',
                        'time_s', 'FeCN6_M_mean', 'FeCN6_M_std', 'n_replicates']

    averaged = averaged[['experiment_type', 'pH', 'Fe_init_M', 'HCN_init_M',
                        'time_s', 'FeCN6_M_mean', 'FeCN6_M_std', 'n_replicates']]

    averaged = averaged.sort_values(['experiment_type', 'Fe_init_M', 'HCN_init_M', 'time_s'])

    print(f"\nAveraged data contains {len(averaged)} rows")
    print(f"Number of replicates per condition: {averaged['n_replicates'].unique()}")

    # The standard error of the mean is calculated here.
    averaged['FeCN6_M_sem'] = averaged['FeCN6_M_std'] / np.sqrt(averaged['n_replicates'])

    averaged.to_csv(output_file, index=False)
    print(f"\nAveraged data saved to: {output_file}")

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"\nConditions at pH {target_pH}:")
    condition_summary = averaged.groupby(['experiment_type', 'Fe_init_M', 'HCN_init_M']).size()
    for idx, count in condition_summary.items():
        exp_type, fe, hcn = idx
        print(f"  {exp_type}: Fe={fe:.4f} M, HCN={hcn if pd.notna(hcn) else 'N/A'} M -> {count} time points")

    return averaged

if __name__ == "__main__":
    # The base path is set here.
    base_path = ""  
    target_pH = 9

    # Step 1: every raw CSV is compiled into one combined dataset, with metadata attached.
    print("Starting data processing...")
    combined_data = process_all_data(base_path)

    if combined_data is not None:
        combined_output_file = "combined_kinetics_data.csv"
        combined_data.to_csv(combined_output_file, index=False)
        print(f"\nData processing complete!")
        print(f"Total rows: {len(combined_data)}")
        print(f"Output saved to: {combined_output_file}")

        print("\nData summary:")
        print(combined_data.groupby(['round', 'experiment_type', 'pH']).size())

        # Step 2: the combined dataset is filtered to target_pH and averaged across the repeated rounds.
        print("\n" + "="*60)
        print(f"Averaging across rounds at pH {target_pH}...")
        print("="*60)
        averaged_output_file = f"averaged_kinetics_pH{target_pH}.csv"
        averaged_data = average_kinetics_data(combined_data, averaged_output_file, target_pH=target_pH)

        if averaged_data is not None:
            print("\n" + "="*60)
            print("Processing complete!")
            print("="*60)
            print("\nFirst few rows of averaged data:")
            print(averaged_data.head(20))

            print("\nColumn descriptions:")
            print("  - time_s: Time in seconds")
            print("  - FeCN6_M_mean: Average [Fe(CN)6^4-] concentration in M")
            print("  - FeCN6_M_std: Standard deviation across rounds")
            print("  - FeCN6_M_sem: Standard error of the mean (SEM)")
            print("  - n_replicates: Number of rounds averaged (should be 3)")
    else:
        print("No data to save.")
