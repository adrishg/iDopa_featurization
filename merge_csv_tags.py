# -*- coding: utf-8 -*-
import pandas as pd
import os
import argparse
import sys

# --- Configuration Defaults (can be overridden by command-line arguments) ---
DEFAULT_TAG_COLUMN_NAME = 'Tag'  # The expected name for your tag column in both CSVs
DEFAULT_COLUMNS_TO_MERGE_FROM_SECONDARY = []

# --- Helper Function for User Confirmation ---
def confirm_action(prompt):
    while True:
        response = input(prompt + " (y/n): ").strip().lower()
        if response == 'y':
            return True
        elif response == 'n':
            return False
        else:
            print("Invalid input. Please enter 'y' or 'n'.")

# --- Normalization Function ---
def normalize_tag(tag):
    """Convert tag to lowercase, replace underscores with spaces, and strip whitespace."""
    return str(tag).replace('_', ' ').strip().lower()

# --- Script Logic ---
def merge_csv_files(primary_csv_path, secondary_csv_path, output_csv_path, ref_column_name, columns_to_merge):
    if not os.path.exists(primary_csv_path):
        print(f"Error: Primary CSV file not found at '{primary_csv_path}'. Please check the path.")
        return
    if not os.path.exists(secondary_csv_path):
        print(f"Error: Secondary CSV file not found at '{secondary_csv_path}'. Please check the path.")
        return

    try:
        df_primary = pd.read_csv(primary_csv_path, header=0)
        df_secondary = pd.read_csv(secondary_csv_path, header=0)

        print(f"Primary CSV loaded. Rows: {len(df_primary)}, Columns: {df_primary.columns.tolist()}")
        print(f"Secondary CSV loaded. Rows: {len(df_secondary)}, Columns: {df_secondary.columns.tolist()}")

        primary_tag_col_to_use = ref_column_name
        secondary_tag_col_to_use = ref_column_name

        if ref_column_name is None:
            if df_primary.columns[0] == DEFAULT_TAG_COLUMN_NAME:
                primary_tag_col_to_use = DEFAULT_TAG_COLUMN_NAME
                print(f"Primary CSV: Found default tag column '{DEFAULT_TAG_COLUMN_NAME}' as the first column.")
            else:
                first_col_header_primary = df_primary.columns[0]
                if confirm_action(f"Primary CSV: The first column header is '{first_col_header_primary}'. Is this your tag column?"):
                    primary_tag_col_to_use = first_col_header_primary
                elif DEFAULT_TAG_COLUMN_NAME in df_primary.columns:
                    if confirm_action(f"Primary CSV: '{DEFAULT_TAG_COLUMN_NAME}' found elsewhere. Do you want to use it?"):
                        primary_tag_col_to_use = DEFAULT_TAG_COLUMN_NAME
                    else:
                        print("Please re-run the script and specify the correct tag column using --ref_column.")
                        return
                else:
                    print("Could not determine primary tag column. Exiting.")
                    return

            if df_secondary.columns[0] == DEFAULT_TAG_COLUMN_NAME:
                secondary_tag_col_to_use = DEFAULT_TAG_COLUMN_NAME
                print(f"Secondary CSV: Found default tag column '{DEFAULT_TAG_COLUMN_NAME}' as the first column.")
            else:
                first_col_header_secondary = df_secondary.columns[0]
                if confirm_action(f"Secondary CSV: The first column header is '{first_col_header_secondary}'. Is this your tag column?"):
                    secondary_tag_col_to_use = first_col_header_secondary
                elif DEFAULT_TAG_COLUMN_NAME in df_secondary.columns:
                    if confirm_action(f"Secondary CSV: '{DEFAULT_TAG_COLUMN_NAME}' found elsewhere. Do you want to use it?"):
                        secondary_tag_col_to_use = DEFAULT_TAG_COLUMN_NAME
                    else:
                        print("Please re-run the script and specify the correct tag column using --ref_column.")
                        return
                else:
                    print("Could not determine secondary tag column. Exiting.")
                    return

        if primary_tag_col_to_use not in df_primary.columns:
            raise ValueError(f"Primary tag column '{primary_tag_col_to_use}' not found in primary CSV.")
        if secondary_tag_col_to_use not in df_secondary.columns:
            raise ValueError(f"Secondary tag column '{secondary_tag_col_to_use}' not found in secondary CSV.")

        print(f"Using '{primary_tag_col_to_use}' from primary CSV and '{secondary_tag_col_to_use}' from secondary CSV for matching.")

        for col_name in columns_to_merge:
            if col_name not in df_secondary.columns:
                print(f"Warning: Column '{col_name}' from 'columns_to_merge' not found in the secondary CSV. Skipping this column.")
                continue
            df_primary[f'Merged_{col_name}'] = ''

        df_secondary['__normalized_tag__'] = df_secondary[secondary_tag_col_to_use].apply(normalize_tag)

        for index, primary_row in df_primary.iterrows():
            primary_tag = str(primary_row[primary_tag_col_to_use]).strip()
            normalized_primary_tag = normalize_tag(primary_tag)

            matching_secondary_rows = df_secondary[df_secondary['__normalized_tag__'] == normalized_primary_tag]

            if not matching_secondary_rows.empty:
                secondary_data = matching_secondary_rows.iloc[0]
                for col_name in columns_to_merge:
                    if col_name in secondary_data:
                        df_primary.at[index, f'Merged_{col_name}'] = secondary_data[col_name]

        df_secondary.drop(columns='__normalized_tag__', inplace=True)

        df_primary.to_csv(output_csv_path, index=False)
        print(f"\nSuccessfully merged data and saved to '{output_csv_path}'.")

    except Exception as e:
        print(f"An error occurred during merging: {e}")
        sys.exit(1)

# --- Main Execution ---
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge two CSV files based on a common tag column, allowing for flexible tag matching.")
    parser.add_argument('--primary_csv', type=str, required=True, help="Path to the primary CSV file.")
    parser.add_argument('--secondary_csv', type=str, required=True, help="Path to the secondary CSV file.")
    parser.add_argument('--output_csv', type=str, default='merged_output.csv', help="Output path for the merged CSV.")
    parser.add_argument('--ref_column', type=str, help=f"Optional: Tag column to use for both CSVs. If not provided, attempts to use '{DEFAULT_TAG_COLUMN_NAME}'.")
    parser.add_argument('--columns_to_merge', nargs='+', required=True, help="List of column headers from secondary CSV to merge.")

    args = parser.parse_args()

    merge_csv_files(args.primary_csv, args.secondary_csv, args.output_csv, args.ref_column, args.columns_to_merge)

