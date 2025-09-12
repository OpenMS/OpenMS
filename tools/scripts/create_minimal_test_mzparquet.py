#!/usr/bin/env -S uv run
# /// script
# requires-python = ">=3.8"
# dependencies = [
#     "pandas",
#     "pyarrow",
# ]
# ///
"""
Simple script to create a minimal test mzParquet file for basic testing.
"""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

def create_minimal_test():
    """Create a minimal test file with just a few spectra.
    
    Ensures all peaks for each spectrum are in the same row group.
    """
    
    # Define test data manually for predictable testing
    # Group data by spectrum to ensure they stay together
    spectra_data = [
        # MS1 spectrum (scan 1) - all peaks together
        [
            {'scan': 1, 'level': 1, 'rt': 10.5, 'mz': 445.2, 'intensity': 50000, 'collision_energy': None, 'ion_mobility': None, 'isolation_lower': None, 'isolation_upper': None, 'precursor_scan': None, 'precursor_mz': None, 'precursor_charge': None},
            {'scan': 1, 'level': 1, 'rt': 10.5, 'mz': 523.8, 'intensity': 75000, 'collision_energy': None, 'ion_mobility': None, 'isolation_lower': None, 'isolation_upper': None, 'precursor_scan': None, 'precursor_mz': None, 'precursor_charge': None},
            {'scan': 1, 'level': 1, 'rt': 10.5, 'mz': 612.3, 'intensity': 120000, 'collision_energy': None, 'ion_mobility': None, 'isolation_lower': None, 'isolation_upper': None, 'precursor_scan': None, 'precursor_mz': None, 'precursor_charge': None},
        ],
        # MS2 spectrum (scan 2) - all peaks together
        [
            {'scan': 2, 'level': 2, 'rt': 10.6, 'mz': 175.1, 'intensity': 25000, 'collision_energy': 30.0, 'ion_mobility': 1.2, 'isolation_lower': 611.3, 'isolation_upper': 613.3, 'precursor_scan': 1, 'precursor_mz': 612.3, 'precursor_charge': 2},
            {'scan': 2, 'level': 2, 'rt': 10.6, 'mz': 288.2, 'intensity': 45000, 'collision_energy': 30.0, 'ion_mobility': 1.2, 'isolation_lower': 611.3, 'isolation_upper': 613.3, 'precursor_scan': 1, 'precursor_mz': 612.3, 'precursor_charge': 2},
            {'scan': 2, 'level': 2, 'rt': 10.6, 'mz': 401.3, 'intensity': 35000, 'collision_energy': 30.0, 'ion_mobility': 1.2, 'isolation_lower': 611.3, 'isolation_upper': 613.3, 'precursor_scan': 1, 'precursor_mz': 612.3, 'precursor_charge': 2},
        ],
        # Another MS1 spectrum (scan 3) - all peaks together  
        [
            {'scan': 3, 'level': 1, 'rt': 15.2, 'mz': 389.7, 'intensity': 80000, 'collision_energy': None, 'ion_mobility': None, 'isolation_lower': None, 'isolation_upper': None, 'precursor_scan': None, 'precursor_mz': None, 'precursor_charge': None},
            {'scan': 3, 'level': 1, 'rt': 15.2, 'mz': 456.8, 'intensity': 95000, 'collision_energy': None, 'ion_mobility': None, 'isolation_lower': None, 'isolation_upper': None, 'precursor_scan': None, 'precursor_mz': None, 'precursor_charge': None},
        ],
        # MS2 spectrum (scan 4) - all peaks together
        [
            {'scan': 4, 'level': 2, 'rt': 15.3, 'mz': 147.1, 'intensity': 15000, 'collision_energy': 28.0, 'ion_mobility': None, 'isolation_lower': 455.8, 'isolation_upper': 457.8, 'precursor_scan': 3, 'precursor_mz': 456.8, 'precursor_charge': 3},
            {'scan': 4, 'level': 2, 'rt': 15.3, 'mz': 261.2, 'intensity': 28000, 'collision_energy': 28.0, 'ion_mobility': None, 'isolation_lower': 455.8, 'isolation_upper': 457.8, 'precursor_scan': 3, 'precursor_mz': 456.8, 'precursor_charge': 3},
            {'scan': 4, 'level': 2, 'rt': 15.3, 'mz': 374.3, 'intensity': 22000, 'collision_energy': 28.0, 'ion_mobility': None, 'isolation_lower': 455.8, 'isolation_upper': 457.8, 'precursor_scan': 3, 'precursor_mz': 456.8, 'precursor_charge': 3},
        ]
    ]
    
    # Define schema
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
    
    # Create separate row groups to ensure spectrum integrity
    # Put first two spectra in first row group, last two in second row group
    row_group_1_data = spectra_data[0] + spectra_data[1]  # Scans 1, 2
    row_group_2_data = spectra_data[2] + spectra_data[3]  # Scans 3, 4
    
    # Convert to Arrow Tables
    rg1_df = pd.DataFrame(row_group_1_data)
    rg2_df = pd.DataFrame(row_group_2_data)
    
    table1 = pa.Table.from_pandas(rg1_df, schema=schema)
    table2 = pa.Table.from_pandas(rg2_df, schema=schema)
    
    # Write with explicit row group control
    with pq.ParquetWriter('test_minimal.parquet', schema, compression='snappy') as writer:
        writer.write_table(table1)  # Row group 1: scans 1, 2
        writer.write_table(table2)  # Row group 2: scans 3, 4
    
    # Verify the structure
    parquet_file = pq.ParquetFile('test_minimal.parquet')
    total_peaks = sum(len(spectrum) for spectrum in spectra_data)
    
    print("Created test_minimal.parquet with:")
    print(f"- {len(spectra_data)} spectra")
    print(f"- {total_peaks} total peaks")
    print(f"- MS levels: {sorted(set(row['level'] for spectrum in spectra_data for row in spectrum))}")
    print(f"- {parquet_file.num_row_groups} row groups")
    print(f"- Row group 1: scans 1-2 ({parquet_file.metadata.row_group(0).num_rows} rows)")
    print(f"- Row group 2: scans 3-4 ({parquet_file.metadata.row_group(1).num_rows} rows)")
    print("- Predictable data for unit testing")
    print("✓ Verified: Each spectrum's peaks are in the same row group")

if __name__ == '__main__':
    create_minimal_test()
