#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
"""
IDMassAccuracy - Calculates a distribution of the mass error from given mass spectra and IDs.

This tool computes mass error distributions by comparing theoretical masses from peptide
identifications against experimental masses from mass spectra. It generates statistics
for both precursor and fragment ion mass errors.
"""

import argparse
import sys
from dataclasses import dataclass
from typing import List, Tuple, Optional

import pyopenms as pms


@dataclass
class MassDifference:
    """Stores mass difference information for a single peak comparison."""
    exp_mz: float
    theo_mz: float
    charge: int
    intensity: float


def get_mass_difference(exp_mz: float, theo_mz: float, use_ppm: bool) -> float:
    """
    Calculate mass difference between experimental and theoretical m/z values.

    Args:
        exp_mz: Experimental m/z value
        theo_mz: Theoretical m/z value
        use_ppm: If True, return error in ppm; otherwise in Da

    Returns:
        Mass difference in Da or ppm
    """
    diff = exp_mz - theo_mz
    if use_ppm:
        return diff / theo_mz * 1e6
    return diff


def calculate_statistics(values: List[float]) -> Tuple[float, float, float]:
    """
    Calculate mean, median, and standard deviation of a list of values.

    Args:
        values: List of numerical values

    Returns:
        Tuple of (mean, median, standard_deviation)
    """
    if not values:
        return 0.0, 0.0, 0.0

    n = len(values)
    mean = sum(values) / n

    sorted_vals = sorted(values)
    if n % 2 == 0:
        median = (sorted_vals[n // 2 - 1] + sorted_vals[n // 2]) / 2
    else:
        median = sorted_vals[n // 2]

    variance = sum((x - mean) ** 2 for x in values) / n if n > 0 else 0
    std_dev = variance ** 0.5

    return mean, median, std_dev


def fit_gaussian(data: List[float], num_bins: int = 100) -> Optional[Tuple[float, float, float]]:
    """
    Fit a Gaussian distribution to the data using histogram-based approach.

    Args:
        data: List of values to fit
        num_bins: Number of histogram bins

    Returns:
        Tuple of (amplitude, mean, sigma) or None if fitting fails
    """
    if len(data) < 3:
        return None

    try:
        # Create histogram
        min_val, max_val = min(data), max(data)
        if min_val == max_val:
            return None

        bin_width = (max_val - min_val) / num_bins
        histogram = [0] * num_bins

        for val in data:
            bin_idx = min(int((val - min_val) / bin_width), num_bins - 1)
            histogram[bin_idx] += 1

        # Find peak for initial guess
        max_count = max(histogram)
        peak_idx = histogram.index(max_count)
        peak_x = min_val + (peak_idx + 0.5) * bin_width

        # Use pyopenms GaussFitter if available
        gf = pms.GaussFitter()
        points = []
        for i in range(num_bins):
            x = min_val + (i + 0.5) * bin_width
            y = float(histogram[i])
            if y > 0:
                points.append(pms.DPosition2(x, y))

        if len(points) < 3:
            return None

        # Set initial parameters (A=amplitude, x0=mean, sigma=width)
        init_result = pms.GaussFitResult(float(max_count), peak_x, (max_val - min_val) / 4.0)
        gf.setInitialParameters(init_result)

        result = gf.fit(points)
        return result.A, result.x0, result.sigma

    except Exception:
        # Fallback to simple statistics-based estimate
        mean, _, std = calculate_statistics(data)
        if std == 0:
            return None
        amplitude = len(data) / (std * (2 * 3.14159) ** 0.5)
        return amplitude, mean, std


def process_precursor_masses(
    exp: pms.MSExperiment,
    peptide_ids: List[pms.PeptideIdentification],
    precursor_error_ppm: bool
) -> List[MassDifference]:
    """
    Calculate precursor mass differences.

    Args:
        exp: MS experiment containing spectra
        peptide_ids: List of peptide identifications
        precursor_error_ppm: If True, calculate error in ppm

    Returns:
        List of MassDifference objects for precursors
    """
    mass_diffs = []

    for pep_id in peptide_ids:
        hits = pep_id.getHits()
        if not hits:
            continue

        # Get experimental precursor m/z
        exp_mz = pep_id.getMZ()
        if exp_mz <= 0:
            continue

        # Get best hit
        best_hit = hits[0]
        seq = best_hit.getSequence()
        charge = best_hit.getCharge()
        if charge == 0:
            charge = pep_id.getMetaValue("spectrum_reference") if pep_id.metaValueExists("spectrum_reference") else 1

        # Calculate theoretical m/z
        theo_mass = seq.getMonoWeight(pms.Residue.ResidueType.Full, charge)
        theo_mz = theo_mass / charge if charge > 0 else theo_mass

        # Find corresponding spectrum to get intensity
        intensity = 1.0  # Default intensity
        rt = pep_id.getRT()
        for spec in exp:
            if abs(spec.getRT() - rt) < 0.01:
                precursors = spec.getPrecursors()
                if precursors:
                    intensity = precursors[0].getIntensity()
                break

        mass_diffs.append(MassDifference(
            exp_mz=exp_mz,
            theo_mz=theo_mz,
            charge=charge,
            intensity=intensity
        ))

    return mass_diffs


def process_fragment_masses(
    exp: pms.MSExperiment,
    peptide_ids: List[pms.PeptideIdentification],
    fragment_mass_tolerance: float,
    fragment_error_ppm: bool
) -> List[MassDifference]:
    """
    Calculate fragment ion mass differences.

    Args:
        exp: MS experiment containing spectra
        peptide_ids: List of peptide identifications
        fragment_mass_tolerance: Tolerance for spectrum alignment
        fragment_error_ppm: If True, calculate error in ppm

    Returns:
        List of MassDifference objects for fragment ions
    """
    mass_diffs = []

    # Set up theoretical spectrum generator
    tsg = pms.TheoreticalSpectrumGenerator()
    tsg_params = tsg.getParameters()
    tsg_params.setValue("add_metainfo", "true")
    tsg_params.setValue("add_b_ions", "true")
    tsg_params.setValue("add_y_ions", "true")
    tsg_params.setValue("add_a_ions", "false")
    tsg_params.setValue("add_losses", "false")
    tsg.setParameters(tsg_params)

    # Set up spectrum aligner
    aligner = pms.SpectrumAlignment()
    aligner_params = aligner.getParameters()
    aligner_params.setValue("tolerance", fragment_mass_tolerance)
    aligner_params.setValue("is_relative_tolerance", "false")
    aligner.setParameters(aligner_params)

    # Create lookup for spectra by RT
    spec_lookup = {}
    for i, spec in enumerate(exp):
        if spec.getMSLevel() == 2:
            spec_lookup[round(spec.getRT(), 4)] = i

    # Process each identification
    for pep_id in peptide_ids:
        hits = pep_id.getHits()
        if not hits:
            continue

        best_hit = hits[0]
        seq = best_hit.getSequence()
        charge = best_hit.getCharge()
        if charge == 0:
            charge = 1

        # Find corresponding experimental spectrum
        rt = pep_id.getRT()
        rt_key = round(rt, 4)

        exp_spec = None
        # Find closest spectrum
        min_diff = float('inf')
        for spec in exp:
            if spec.getMSLevel() == 2:
                diff = abs(spec.getRT() - rt)
                if diff < min_diff:
                    min_diff = diff
                    exp_spec = spec

        if exp_spec is None or len(exp_spec) == 0:
            continue

        # Generate theoretical spectrum
        theo_spec = pms.MSSpectrum()
        try:
            tsg.getSpectrum(theo_spec, seq, 1, min(charge, 2))
        except Exception:
            continue

        if len(theo_spec) == 0:
            continue

        # Sort spectra
        exp_spec.sortByPosition()
        theo_spec.sortByPosition()

        # Align spectra
        alignment = []
        try:
            aligner.getSpectrumAlignment(alignment, theo_spec, exp_spec)
        except Exception:
            continue

        # Collect mass differences from aligned peaks
        for theo_idx, exp_idx in alignment:
            theo_peak = theo_spec[theo_idx]
            exp_peak = exp_spec[exp_idx]

            mass_diffs.append(MassDifference(
                exp_mz=exp_peak.getMZ(),
                theo_mz=theo_peak.getMZ(),
                charge=1,  # Fragment ions typically charge 1 or 2
                intensity=exp_peak.getIntensity()
            ))

    return mass_diffs


def write_mass_differences(
    filename: str,
    mass_diffs: List[MassDifference],
    use_ppm: bool
) -> None:
    """
    Write mass differences to a TSV file.

    Args:
        filename: Output file path
        mass_diffs: List of mass difference data
        use_ppm: If True, output in ppm; otherwise in Da
    """
    with open(filename, 'w') as f:
        # Write header
        unit = "ppm" if use_ppm else "Da"
        f.write(f"exp_mz\ttheo_mz\tcharge\tintensity\tmass_error_{unit}\n")

        # Write data
        for md in mass_diffs:
            error = get_mass_difference(md.exp_mz, md.theo_mz, use_ppm)
            f.write(f"{md.exp_mz:.6f}\t{md.theo_mz:.6f}\t{md.charge}\t{md.intensity:.2f}\t{error:.6f}\n")


def write_statistics(
    filename: str,
    mass_diffs: List[MassDifference],
    use_ppm: bool,
    num_bins: int,
    data_type: str
) -> None:
    """
    Write statistics and Gaussian fit to a file.

    Args:
        filename: Output file path
        mass_diffs: List of mass difference data
        use_ppm: If True, calculate in ppm
        num_bins: Number of histogram bins
        data_type: Description of data (e.g., "precursor" or "fragment")
    """
    if not mass_diffs:
        return

    errors = [get_mass_difference(md.exp_mz, md.theo_mz, use_ppm) for md in mass_diffs]
    mean, median, std = calculate_statistics(errors)
    unit = "ppm" if use_ppm else "Da"

    with open(filename, 'a') as f:
        f.write(f"\n# Statistics for {data_type} mass errors ({unit}):\n")
        f.write(f"# Number of data points: {len(errors)}\n")
        f.write(f"# Mean: {mean:.6f} {unit}\n")
        f.write(f"# Median: {median:.6f} {unit}\n")
        f.write(f"# Standard deviation: {std:.6f} {unit}\n")

        # Gaussian fit
        fit_result = fit_gaussian(errors, num_bins)
        if fit_result:
            amp, x0, sigma = fit_result
            f.write(f"# Gaussian fit - Amplitude: {amp:.4f}, Mean: {x0:.6f}, Sigma: {sigma:.6f}\n")


def id_mass_accuracy(
    in_files: List[str],
    id_files: List[str],
    precursor_out: Optional[str],
    fragment_out: Optional[str],
    precursor_error_ppm: bool,
    fragment_error_ppm: bool,
    fragment_mass_tolerance: float,
    num_bins: int
) -> None:
    """
    Main processing function for ID mass accuracy calculation.

    Args:
        in_files: List of input mzML files
        id_files: List of input idXML files
        precursor_out: Output file for precursor mass errors
        fragment_out: Output file for fragment mass errors
        precursor_error_ppm: Calculate precursor error in ppm
        fragment_error_ppm: Calculate fragment error in ppm
        fragment_mass_tolerance: Tolerance for fragment alignment
        num_bins: Number of histogram bins
    """
    all_precursor_diffs = []
    all_fragment_diffs = []

    # Process each pair of files
    for in_file, id_file in zip(in_files, id_files):
        print(f"Processing {in_file} with {id_file}...")

        # Load mzML file
        exp = pms.MSExperiment()
        pms.MzMLFile().load(in_file, exp)

        # Normalize spectra
        normalizer = pms.Normalizer()
        for spec in exp:
            if spec.getMSLevel() == 2:
                normalizer.filterSpectrum(spec)

        # Load identifications
        protein_ids = []
        peptide_ids = []
        pms.IdXMLFile().load(id_file, protein_ids, peptide_ids)

        print(f"  Loaded {len(exp)} spectra and {len(peptide_ids)} peptide identifications")

        # Calculate precursor mass differences
        if precursor_out:
            precursor_diffs = process_precursor_masses(exp, peptide_ids, precursor_error_ppm)
            all_precursor_diffs.extend(precursor_diffs)
            print(f"  Found {len(precursor_diffs)} precursor mass comparisons")

        # Calculate fragment mass differences
        if fragment_out:
            fragment_diffs = process_fragment_masses(
                exp, peptide_ids, fragment_mass_tolerance, fragment_error_ppm
            )
            all_fragment_diffs.extend(fragment_diffs)
            print(f"  Found {len(fragment_diffs)} fragment mass comparisons")

    # Write outputs
    if precursor_out and all_precursor_diffs:
        write_mass_differences(precursor_out, all_precursor_diffs, precursor_error_ppm)
        write_statistics(precursor_out, all_precursor_diffs, precursor_error_ppm, num_bins, "precursor")
        print(f"Wrote {len(all_precursor_diffs)} precursor mass errors to {precursor_out}")

    if fragment_out and all_fragment_diffs:
        write_mass_differences(fragment_out, all_fragment_diffs, fragment_error_ppm)
        write_statistics(fragment_out, all_fragment_diffs, fragment_error_ppm, num_bins, "fragment")
        print(f"Wrote {len(all_fragment_diffs)} fragment mass errors to {fragment_out}")

    # Print summary statistics
    if all_precursor_diffs:
        errors = [get_mass_difference(md.exp_mz, md.theo_mz, precursor_error_ppm) for md in all_precursor_diffs]
        mean, median, std = calculate_statistics(errors)
        unit = "ppm" if precursor_error_ppm else "Da"
        print(f"\nPrecursor mass error statistics ({unit}):")
        print(f"  Mean: {mean:.6f}, Median: {median:.6f}, Std: {std:.6f}")

    if all_fragment_diffs:
        errors = [get_mass_difference(md.exp_mz, md.theo_mz, fragment_error_ppm) for md in all_fragment_diffs]
        mean, median, std = calculate_statistics(errors)
        unit = "ppm" if fragment_error_ppm else "Da"
        print(f"\nFragment mass error statistics ({unit}):")
        print(f"  Mean: {mean:.6f}, Median: {median:.6f}, Std: {std:.6f}")


def main():
    """Main entry point for the IDMassAccuracy tool."""
    parser = argparse.ArgumentParser(
        description="IDMassAccuracy - Calculates mass error distributions from spectra and IDs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python IDMassAccuracy.py -in sample.mzML -id sample.idXML -precursor_out precursor_errors.tsv -fragment_out fragment_errors.tsv

This tool calculates the mass error distribution by comparing:
  - Precursor masses: experimental vs. theoretical masses from peptide sequences
  - Fragment masses: experimental peaks vs. theoretical fragment ions (b/y ions)

Output files contain tab-separated values with mass errors and optional Gaussian fit statistics.
        """
    )

    # Input files
    parser.add_argument(
        "-in",
        action="store",
        nargs="+",
        type=str,
        dest="in_files",
        required=True,
        metavar="<file>",
        help="Input mzML file(s) containing spectra"
    )

    parser.add_argument(
        "-id_in",
        action="store",
        nargs="+",
        type=str,
        dest="id_files",
        required=True,
        metavar="<file>",
        help="Input idXML file(s) containing peptide identifications"
    )

    # Output files
    parser.add_argument(
        "-precursor_out",
        action="store",
        type=str,
        metavar="<file>",
        help="Output file for precursor mass errors (TSV format)"
    )

    parser.add_argument(
        "-fragment_out",
        action="store",
        type=str,
        metavar="<file>",
        help="Output file for fragment mass errors (TSV format)"
    )

    # Error calculation options
    parser.add_argument(
        "-precursor_error_ppm",
        action="store_true",
        default=False,
        help="Calculate precursor mass error in ppm (default: Da)"
    )

    parser.add_argument(
        "-fragment_error_ppm",
        action="store_true",
        default=False,
        help="Calculate fragment mass error in ppm (default: Da)"
    )

    # Algorithm parameters
    parser.add_argument(
        "-fragment_mass_tolerance",
        action="store",
        type=float,
        default=0.5,
        metavar="<value>",
        help="Fragment mass tolerance for spectrum alignment in Da (default: 0.5)"
    )

    parser.add_argument(
        "-number_of_bins",
        action="store",
        type=int,
        default=100,
        metavar="<n>",
        help="Number of bins for histogram/Gaussian fitting (default: 100)"
    )

    args = parser.parse_args()

    # Validate inputs
    if len(args.in_files) != len(args.id_files):
        parser.error("Number of mzML files must match number of idXML files")

    if not args.precursor_out and not args.fragment_out:
        parser.error("At least one output file (-precursor_out or -fragment_out) must be specified")

    # Validate file existence
    for f in args.in_files:
        try:
            with open(f, 'r'):
                pass
        except FileNotFoundError:
            parser.error(f"Input file not found: {f}")

    for f in args.id_files:
        try:
            with open(f, 'r'):
                pass
        except FileNotFoundError:
            parser.error(f"ID file not found: {f}")

    # Run the tool
    id_mass_accuracy(
        in_files=args.in_files,
        id_files=args.id_files,
        precursor_out=args.precursor_out,
        fragment_out=args.fragment_out,
        precursor_error_ppm=args.precursor_error_ppm,
        fragment_error_ppm=args.fragment_error_ppm,
        fragment_mass_tolerance=args.fragment_mass_tolerance,
        num_bins=args.number_of_bins
    )


if __name__ == "__main__":
    main()
