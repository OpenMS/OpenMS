#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "pandas",
#     "pyarrow",
#     "numpy",
# ]
# ///
"""
Script to create a test mzParquet file for testing the OpenMS MzParquetFile reader.
This creates a synthetic MS dataset with MS1 and MS2 spectra following the mzParquet schema.
"""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np
import argparse
from pathlib import Path

def create_test_mzparquet(filename, num_ms1_spectra=100, num_ms2_spectra=500, peaks_per_spectrum=50):
    """
    Create a test mzParquet file with synthetic MS data.
    
    IMPORTANT: Ensures that all peaks for a single spectrum are within the same row group.
    This is a key requirement for the OpenMS implementation.
    
    Parameters:
    -----------
    filename : str
        Output filename for the parquet file
    num_ms1_spectra : int
        Number of MS1 spectra to generate
    num_ms2_spectra : int
        Number of MS2 spectra to generate
    peaks_per_spectrum : int
        Average number of peaks per spectrum
    """
    
    print(f"Creating test mzParquet file: {filename}")
    print(f"MS1 spectra: {num_ms1_spectra}, MS2 spectra: {num_ms2_spectra}")
    print(f"Average peaks per spectrum: {peaks_per_spectrum}")
    
    # Group data by row groups to ensure spectra don't span multiple row groups
    row_groups = []
    current_row_group = []
    max_rows_per_group = 10000
    current_row_count = 0
    
    scan_id = 1
    
    # Generate MS1 spectra
    for ms1_idx in range(num_ms1_spectra):
        rt = 10.0 + ms1_idx * 2.0  # RT from 10 to 210 seconds
        num_peaks = np.random.poisson(peaks_per_spectrum)
        
        # Check if this spectrum would fit in current row group
        if current_row_count + num_peaks > max_rows_per_group and current_row_group:
            # Finish current row group and start new one
            row_groups.append(current_row_group)
            current_row_group = []
            current_row_count = 0
        
        # Generate peaks for this spectrum
        mz_values = np.sort(np.random.uniform(100, 2000, num_peaks))
        intensities = np.random.exponential(10000, num_peaks).astype(np.uint32)
        
        spectrum_rows = []
        for mz, intensity in zip(mz_values, intensities):
            spectrum_rows.append({
                'scan': scan_id,
                'level': 1,
                'rt': rt,
                'mz': float(mz),
                'intensity': int(intensity),
                'collision_energy': None,
                'ion_mobility': np.random.uniform(0.5, 1.5) if np.random.random() < 0.3 else None,
                'isolation_lower': None,
                'isolation_upper': None,
                'precursor_scan': None,
                'precursor_mz': None,
                'precursor_charge': None
            })
        
        # Add all peaks of this spectrum to current row group
        current_row_group.extend(spectrum_rows)
        current_row_count += len(spectrum_rows)
        scan_id += 1
        
        # Generate MS2 spectra for some MS1 precursors
        if ms1_idx < num_ms2_spectra:
            num_ms2_for_this_ms1 = np.random.poisson(2) + 1  # 1-5 MS2 per MS1
            
            for ms2_idx in range(num_ms2_for_this_ms1):
                if scan_id > num_ms1_spectra + num_ms2_spectra:
                    break
                
                # Pick a random precursor from the MS1 peaks
                if len(mz_values) > 0:
                    precursor_mz = np.random.choice(mz_values[mz_values > 400])  # Only peptide-like m/z
                    precursor_charge = np.random.choice([2, 3, 4], p=[0.6, 0.3, 0.1])
                    collision_energy = np.random.uniform(25, 35)
                    
                    # MS2 RT slightly after MS1
                    ms2_rt = rt + 0.1 + ms2_idx * 0.05
                    
                    # Generate MS2 peaks
                    num_ms2_peaks = np.random.poisson(peaks_per_spectrum // 2)
                    
                    # Check if this MS2 spectrum would fit in current row group
                    if current_row_count + num_ms2_peaks > max_rows_per_group and current_row_group:
                        # Finish current row group and start new one
                        row_groups.append(current_row_group)
                        current_row_group = []
                        current_row_count = 0
                    
                    ms2_mz_values = np.sort(np.random.uniform(50, precursor_mz, num_ms2_peaks))
                    ms2_intensities = np.random.exponential(5000, num_ms2_peaks).astype(np.uint32)
                    
                    # Isolation window
                    isolation_width = 2.0
                    isolation_lower = precursor_mz - isolation_width/2
                    isolation_upper = precursor_mz + isolation_width/2
                    
                    ms2_spectrum_rows = []
                    for mz, intensity in zip(ms2_mz_values, ms2_intensities):
                        ms2_spectrum_rows.append({
                            'scan': scan_id,
                            'level': 2,
                            'rt': ms2_rt,
                            'mz': float(mz),
                            'intensity': int(intensity),
                            'collision_energy': collision_energy,
                            'ion_mobility': np.random.uniform(0.5, 1.5) if np.random.random() < 0.2 else None,
                            'isolation_lower': isolation_lower,
                            'isolation_upper': isolation_upper,
                            'precursor_scan': scan_id - num_ms2_for_this_ms1 - ms2_idx,  # Reference to MS1
                            'precursor_mz': precursor_mz,
                            'precursor_charge': precursor_charge
                        })
                    
                    # Add all peaks of this MS2 spectrum to current row group
                    current_row_group.extend(ms2_spectrum_rows)
                    current_row_count += len(ms2_spectrum_rows)
                    scan_id += 1
    
    # Add the last row group
    if current_row_group:
        row_groups.append(current_row_group)
    
    print(f"Created {len(row_groups)} row groups")
    print(f"Row group sizes: {[len(rg) for rg in row_groups]}")
    
    # Flatten all row groups into a single DataFrame
    all_data = []
    for rg in row_groups:
        all_data.extend(rg)
    
    df = pd.DataFrame(all_data)
    
    print(f"Generated {len(df)} total peaks across {df['scan'].nunique()} spectra")
    print(f"MS levels: {sorted(df['level'].unique())}")
    print(f"RT range: {df['rt'].min():.1f} - {df['rt'].max():.1f} seconds")
    print(f"m/z range: {df['mz'].min():.1f} - {df['mz'].max():.1f}")
    
    # Verify no spectrum spans multiple row groups
    verify_spectrum_integrity(row_groups)
    
    # Define schema matching the OpenMS implementation
    schema = pa.schema([
        ('scan', pa.uint32()),
        ('level', pa.uint32()),
        ('rt', pa.float32()),
        ('mz', pa.float32()),
        ('intensity', pa.uint32()),
        ('collision_energy', pa.float32()),
        ('ion_mobility', pa.float32()),
        ('isolation_lower', pa.float32()),
        ('isolation_upper', pa.float32()),
        ('precursor_scan', pa.uint32()),
        ('precursor_mz', pa.float32()),
        ('precursor_charge', pa.uint32())
    ])
    
    # Convert row groups to Arrow Tables and write with explicit row group boundaries
    tables = []
    for rg_data in row_groups:
        rg_df = pd.DataFrame(rg_data)
        tables.append(pa.Table.from_pandas(rg_df, schema=schema))
    
    # Write with explicit row group control
    write_parquet_with_rowgroups(tables, filename)
    
    print(f"Successfully created {filename}")
    
    # Print some statistics
    parquet_file = pq.ParquetFile(filename)
    print(f"File size: {Path(filename).stat().st_size / 1024:.1f} KB")
    print(f"Number of row groups: {parquet_file.num_row_groups}")
    print(f"Rows per row group: {[parquet_file.metadata.row_group(i).num_rows for i in range(parquet_file.num_row_groups)]}")

def verify_spectrum_integrity(row_groups):
    """Verify that no spectrum spans multiple row groups."""
    print("Verifying spectrum integrity...")
    
    all_scans_by_group = []
    for i, rg in enumerate(row_groups):
        scans_in_group = set(row['scan'] for row in rg)
        all_scans_by_group.append(scans_in_group)
        print(f"Row group {i}: {len(scans_in_group)} unique scans")
    
    # Check for overlapping scans between row groups
    for i in range(len(all_scans_by_group)):
        for j in range(i + 1, len(all_scans_by_group)):
            overlap = all_scans_by_group[i].intersection(all_scans_by_group[j])
            if overlap:
                raise ValueError(f"ERROR: Scans {overlap} appear in both row group {i} and {j}!")
    
    print("✓ Verified: No spectrum spans multiple row groups")

def write_parquet_with_rowgroups(tables, filename):
    """Write multiple Arrow tables as separate row groups."""
    
    # Create a writer that allows us to control row groups
    schema = tables[0].schema
    
    with pq.ParquetWriter(filename, schema, compression='snappy') as writer:
        for table in tables:
            # Write each table as a separate row group
            writer.write_table(table)
    
    print(f"Written {len(tables)} row groups to {filename}")

def main():
    parser = argparse.ArgumentParser(description='Create test mzParquet files for OpenMS testing')
    parser.add_argument('output', help='Output parquet filename')
    parser.add_argument('--ms1-spectra', type=int, default=100, help='Number of MS1 spectra (default: 100)')
    parser.add_argument('--ms2-spectra', type=int, default=500, help='Number of MS2 spectra (default: 500)')
    parser.add_argument('--peaks-per-spectrum', type=int, default=50, help='Average peaks per spectrum (default: 50)')
    parser.add_argument('--large', action='store_true', help='Create large test file (10k MS1, 50k MS2)')
    
    args = parser.parse_args()
    
    if args.large:
        create_test_mzparquet(
            args.output, 
            num_ms1_spectra=10000, 
            num_ms2_spectra=50000, 
            peaks_per_spectrum=75
        )
    else:
        create_test_mzparquet(
            args.output,
            num_ms1_spectra=args.ms1_spectra,
            num_ms2_spectra=args.ms2_spectra, 
            peaks_per_spectrum=args.peaks_per_spectrum
        )

if __name__ == '__main__':
    main()
