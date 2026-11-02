/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <string>

namespace MzPeak::Schema::PSI {

/**
 * The type of array encoded into a Parquet column, semantically,
 * denoted using a CURIE from the PSI-MS controlled vocabulary for a
 * child of MS:1000513 (binary data array).  If no existing array
 * fits, use MS:1000786 (non-standard array).
 */
enum class ArrayType {
  /// MS:1000514
  ///
  /// A data array of m/z values (floating point).
  Mz,

  /// MS:1000515
  ///
  /// A data array of intensity values (floating point).
  Intensity,

  /// MS:1000516
  ///
  /// A data array of charge values (32-bit integers).
  Charge,

  /// MS:1000517
  ///
  /// A data array of signal-to-noise values (floating point).
  SignalToNoise,

  /// MS:1000595
  ///
  /// A data array of relative time offset values from a reference
  /// time (floating point).
  RelativeTimeOffset,

  /// MS:1000617
  ///
  /// A data array of electromagnetic radiation wavelength values
  /// (floating point).
  ElectromagneticRadiation,

  /// MS:1000820
  ///
  /// A data array of flow rate measurements (floating point).
  FlowRate,

  /// MS:1000821
  ///
  /// A data array of pressure measurements (floating point).
  Pressure,

  /// MS:1000822
  ///
  /// A data array of temperature measurements (floating point).
  Temperature,

  /// MS:1002478
  ///
  /// Array of mean charge values where the mean charge is calculated
  /// as a weighted mean of the charges of individual peaks that are
  /// aggregated into a processed spectrum (32-bit floating point).
  MeanCharge,

  /// MS:1002529
  ///
  /// A data array of resolution values (floating point).
  Resolution,

  /// MS:1002530
  ///
  /// A data array of signal baseline values (the signal in the
  /// absence of analytes) (floating point).
  SignalBaseline,

  /// MS:1002742
  ///
  /// A data array of noise values (floating point).
  Noise,

  /// MS:1002743
  ///
  /// A data array of parallel, independent m/z values for a sampling
  /// of noise across a spectrum (typically much smaller than
  /// MS:1000514, the m/z array) (floating point).
  SampledNoiseMz,

  /// MS:1002744
  ///
  /// A data array of intensity values for the amplitude of noise
  /// variation superposed on the baseline (MS:1002745) across a
  /// spectrum (for use with MS:1002743, sampled noise m/z array)
  /// (floating point).
  SampledNoiseIntensity,

  /// MS:1002745
  ///
  /// A data array of baseline intensity values (the intensity in the
  /// absence of analytes) for a sampling of noise across a spectrum
  /// (for use with MS:1002743, sampled noise m/z array) (floating
  /// point).
  SampledNoiseBaseline,

  /// MS:1002893
  ///
  /// Abstract array of ion mobility data values. A more specific
  /// child term concept should be specified in data files to make
  /// precise the nature of the data being provided.
  IonMobility,

  /// MS:1002816 — mean ion mobility array.
  MeanIonMobility,

  /// MS:1003006 — mean inverse reduced ion mobility array.
  MeanInverseReducedIonMobility,

  /// MS:1003007 — raw ion mobility array.
  RawIonMobility,

  /// MS:1003008 — raw inverse reduced ion mobility array.
  RawInverseReducedIonMobility,

  /// MS:1002477 — mean ion mobility drift time array.
  MeanIonMobilityDriftTime,

  /// MS:1003153 — raw ion mobility drift time array.
  RawIonMobilityDriftTime,

  /// MS:1003143
  ///
  /// A data array of mass values (floating point).
  Mass,

  /// MS:1003157
  ///
  /// Array of m/z values representing the lower bound m/z of the
  /// quadrupole position at each point in the spectrum.
  ScanningQuadrupolePositionLowerBoundMz,

  /// MS:1003158
  ///
  /// Array of m/z values representing the upper bound m/z of the
  /// quadrupole position at each point in the spectrum.
  ScanningQuadrupolePositionUpperBoundMz,

  /// MS:1000786
  ///
  /// Non-standard array.
  NonStandard,
};

/**
 * Return true when @p t is any flavour of ion-mobility array.
 *
 * The abstract term MS:1002893 has several concrete children, and converters
 * emit the concrete ones — MS:1002816 in particular.  Matching only the
 * abstract term therefore selects nothing and the mobility dimension silently
 * disappears rather than failing.
 */
bool is_ion_mobility(ArrayType);

/**
 * Convert an ArrayType to a string.
 */
std::string array_type_to_string(ArrayType);

/**
 * Parse an ArrayType from a string.
 */
ArrayType array_type_from_string(std::string_view);

} // namespace MzPeak::Schema::PSI
