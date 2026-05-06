# -*- coding: utf-8 -*-
import pandas as pd
import os
import argparse
import sys

def normalize_tag(tag):
    return str(tag).replace('_', ' ').strip().lower()

def merge_csv_files(primary_csv_path, secondary_csv_path, output_csv_path, ref_column_name, columns_to_merge):
    if not os.path.exists(primary_csv_path):
        print(f"Error: Primary CSV not found at '{primary_csv_path}'.")
        return
    if not os.path.exists(secondary_csv_path):
        print(f"Error: Secondary CSV not found at '{secondary_csv_path}'.")
        return

    df_primary = pd.read_csv(primary_csv_path)
    df_secondary = pd.read_csv(secondary_csv_path)

    tag_col = ref_column_name or df_primary.columns[0]

    if tag_col not in df_primary.columns:
        raise ValueError(f"Tag column '{tag_col}' not found in primary CSV.")
    if tag_col not in df_secondary.columns:
        raise ValueError(f"Tag column '{tag_col}' not found in secondary CSV.")

    valid_cols = [c for c in columns_to_merge if c in df_secondary.columns]
    for c in set(columns_to_merge) - set(valid_cols):
        print(f"Warning: Column '{c}' not found in secondary CSV. Skipping.")

    df_primary['__norm__'] = df_primary[tag_col].apply(normalize_tag)
    df_secondary['__norm__'] = df_secondary[tag_col].apply(normalize_tag)

    df_merged = df_primary.merge(
        df_secondary[['__norm__'] + valid_cols],
        on='__norm__',
        how='left'
    ).drop(columns='__norm__')

    df_merged.to_csv(output_csv_path, index=False)
    print(f"Merged data saved to '{output_csv_path}'.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge two CSV files on a common tag column.")
    parser.add_argument('--primary_csv', type=str, required=True)
    parser.add_argument('--secondary_csv', type=str, required=True)
    parser.add_argument('--output_csv', type=str, default='merged_output.csv')
    parser.add_argument('--ref_column', type=str, default=None)
    parser.add_argument('--columns_to_merge', nargs='+', required=True)
    args = parser.parse_args()

    merge_csv_files(args.primary_csv, args.secondary_csv, args.output_csv, args.ref_column, args.columns_to_merge)
