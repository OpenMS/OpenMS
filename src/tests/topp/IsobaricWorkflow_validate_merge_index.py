#!/usr/bin/env python3
"""
Validation script for IsobaricWorkflow id_merge_index metadata.

This script checks that:
1. All PeptideIdentifications in the consensusXML have an id_merge_index metavalue
2. The id_merge_index values are consecutive (0, 1, 2, ... for N input files)
3. Map indices in ConsensusFeatures are consecutive

Usage:
    python IsobaricWorkflow_validate_merge_index.py <consensusXML_file> <expected_num_files>
"""

import sys
import xml.etree.ElementTree as ET


def validate_consensus_xml(filename, expected_num_files):
    """Validate id_merge_index in consensusXML file."""
    tree = ET.parse(filename)
    root = tree.getroot()
    
    # Track id_merge_index values found
    merge_indices_found = set()
    features_without_merge_index = 0
    total_features = 0
    
    # Track map indices for consecutive check
    map_indices = set()
    
    # Check each consensus feature
    for feature in root.findall('.//consensusElement'):
        total_features += 1
        
        # Check map indices
        for element in feature.findall('.//element'):
            map_attrib = element.get('map')
            if map_attrib is not None:
                map_indices.add(int(map_attrib))
        
        # Check PeptideIdentification
        for pep_id in feature.findall('.//PeptideIdentification'):
            # Look for id_merge_index in UserParam
            found_merge_index = False
            for user_param in pep_id.findall('.//UserParam'):
                if user_param.get('name') == 'id_merge_index':
                    merge_index = int(user_param.get('value'))
                    merge_indices_found.add(merge_index)
                    found_merge_index = True
                    break
            
            if not found_merge_index:
                features_without_merge_index += 1
    
    # Validate results
    errors = []
    
    if features_without_merge_index > 0:
        errors.append(f"Found {features_without_merge_index} PeptideIdentifications without id_merge_index")
    
    if len(merge_indices_found) != expected_num_files:
        errors.append(f"Expected {expected_num_files} unique id_merge_index values, found {len(merge_indices_found)}: {sorted(merge_indices_found)}")
    
    # Check if merge indices are consecutive starting from 0
    expected_indices = set(range(expected_num_files))
    if merge_indices_found != expected_indices:
        errors.append(f"id_merge_index values are not consecutive. Expected {expected_indices}, found {merge_indices_found}")
    
    # Check if map indices are consecutive
    if map_indices:
        min_map = min(map_indices)
        max_map = max(map_indices)
        expected_map_indices = set(range(min_map, max_map + 1))
        if map_indices != expected_map_indices:
            errors.append(f"Map indices are not consecutive. Expected {expected_map_indices}, found gaps in {map_indices}")
    
    # Print results
    print(f"Validated {total_features} consensus features")
    print(f"Found id_merge_index values: {sorted(merge_indices_found)}")
    print(f"Map indices range: {min(map_indices) if map_indices else 'N/A'} to {max(map_indices) if map_indices else 'N/A'}")
    
    if errors:
        print("\nERRORS:")
        for error in errors:
            print(f"  - {error}")
        return 1
    else:
        print("\n✓ All checks passed!")
        return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    
    consensus_file = sys.argv[1]
    expected_files = int(sys.argv[2])
    
    sys.exit(validate_consensus_xml(consensus_file, expected_files))
