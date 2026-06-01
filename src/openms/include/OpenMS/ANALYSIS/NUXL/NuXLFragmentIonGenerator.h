// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentAdductDefinition.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentAnnotationHelper.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <vector>
#include <map>
#include <set>
#include <iostream>

namespace OpenMS
{  

/**
  @brief Augments a peptide-cross-link fragment spectrum with NuXL-specific ions that @ref TheoreticalSpectrumGenerator does not produce.

  All entry points are static and append to caller-provided @c PeakSpectrum + parallel
  @c IntegerDataArray (charge) + @c StringDataArray (annotation) arrays. Every appended
  peak has intensity @c 1.0f and charge @c 1 unless a method documents otherwise. None of
  the methods sort the output — call @c sortByPosition afterwards if order matters (the
  composite @ref generatePartialLossSpectrum entry point does this on the final result).

  Ion kinds covered:
    - RNA marker ions (protonated nitrogenous base and its loss variants, e.g. @c U-H2O).
    - Shifted residue-specific immonium ions (Y, W, F, H, C, P, L/I, K, M) carrying a
      configurable mass shift.
    - Non-classical Lys-related immonium variants (@c iK(C5H10N1), @c iK(C6H13N2O)) that
      empirically appear alongside the classical K immonium.
    - Precursor peaks of the peptide + complete-RNA adduct for every charge @c <= the
      precursor charge (used mainly for ETD-type spectra).

  @ingroup Analysis_ID
*/
class OPENMS_DLLAPI NuXLFragmentIonGenerator
{
  public:
  /// Annotation prefix used for marker ions appended by @ref addMS2MarkerIons (e.g. @c "MI:U-H2O").
  static constexpr const char* ANNOTATIONS_MARKER_ION_PREFIX = "MI:";

  /**
    @brief Append singly-charged RNA marker ions (protonated nitrogenous bases and their loss variants) to @p spectrum.

    For each entry in @p marker_ions, appends one peak at @c m.mass @c + @c PROTON_MASS_U
    with intensity @c 1.0f, charge @c 1, and annotation
    @c ANNOTATIONS_MARKER_ION_PREFIX @c + @c m.name (e.g. @c "MI:U-H2O").
    Peaks are appended unconditionally — no collision check against existing entries.

    @param[in]     marker_ions         Marker-ion definitions (name + monoisotopic mass).
    @param[in,out] spectrum            Spectrum to extend; new peaks are appended.
    @param[in,out] spectrum_charge     Parallel charge array; one @c 1 appended per marker ion.
    @param[in,out] spectrum_annotation Parallel annotation array; one prefixed name appended per marker ion.
  */
  static void addMS2MarkerIons(
    const std::vector<NuXLFragmentAdductDefinition>& marker_ions,
    PeakSpectrum& spectrum,
    PeakSpectrum::IntegerDataArray& spectrum_charge,
    PeakSpectrum::StringDataArray& spectrum_annotation);

  /**
    @brief Append shifted residue-specific immonium ions to @p partial_loss_spectrum for amino acids found in @p unmodified_sequence.

    For every supported residue present in @p unmodified_sequence, appends an immonium ion
    at the residue's nominal immonium m/z @em plus @p fragment_shift_mass. Residues covered:
    Y, W, F, H, C, P, L (matches I as well), K, and M. K additionally emits two non-classical
    variants (@c iK(C5H10N1) and @c iK(C6H13N2O)). M additionally emits a methionine-related
    @c CH5S fragment.

    All appended peaks have charge @c 1 and intensity @c 1.0f. Annotations are formatted by
    @ref NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon for the standard immoniums and
    as @c "iK(C5H10N1){shift}" / @c "iK(C6H13N2O){shift}" for the K variants (with @c {shift}
    standing in for the @p fragment_shift_name suffix). No collision
    check against existing entries.

    @param[in]     unmodified_sequence                  Peptide sequence (one-letter codes); only residues present in this string contribute peaks.
    @param[in]     fragment_shift_name                  Annotation suffix (e.g. @c "U-H2O") appended to the K variants and used by @ref NuXLFragmentAnnotationHelper.
    @param[in]     fragment_shift_mass                  Mass shift added to every appended immonium m/z.
    @param[in,out] partial_loss_spectrum                Spectrum to extend.
    @param[in,out] partial_loss_spectrum_charge         Parallel charge array.
    @param[in,out] partial_loss_spectrum_annotation     Parallel annotation array.
  */
  static void addShiftedImmoniumIons(
    const String & unmodified_sequence,
    const String & fragment_shift_name,
    const double fragment_shift_mass,
    PeakSpectrum & partial_loss_spectrum,
    PeakSpectrum::IntegerDataArray& partial_loss_spectrum_charge,
    PeakSpectrum::StringDataArray& partial_loss_spectrum_annotation);

  /**
    @brief Compose a full partial-loss fragment spectrum for one peptide / NA-adduct pair.

    For each @p partial_loss_modification entry @c f:
      -# Calls @ref addShiftedImmoniumIons (charge-1 immonium series for @c f.mass).
      -# For each charge @c z in @c [1, min(precursor_charge, 3)], copies peaks from the
         charge-specific @p patial_loss_template_z{1,2,3} into the output and shifts their
         m/z by @c f.mass @c / @c z. For @c z @c == @c 2 and @c z @c == @c 3 the template is
         filtered down to peaks whose template-charge annotation matches @c z (lower-charge
         peaks already present in the higher-charge templates are skipped).
      -# Charges @c >= @c 4 are not emitted (the loop @c breaks).
      -# Appends @c " " @c + @c f.name to each shifted-peak annotation.

    After the per-modification loop, calls @ref addPrecursorWithCompleteRNA_ for every
    charge in @c [1, precursor_charge] (mainly for ETD-type spectra). Finally,
    @p partial_loss_spectrum is sorted by m/z via @c sortByPosition.

    The output spectrum's @c IntegerDataArrays and @c StringDataArrays are resized to one
    array each on entry.

    @param[in]  unmodified_sequence                       Peptide sequence used for immonium-ion generation.
    @param[in]  fixed_and_variable_modified_peptide_weight Mass of the peptide carrying the fixed and variable modifications (used to compute the precursor-with-complete-RNA m/z).
    @param[in]  precursor_rna_adduct                      String label of the complete RNA adduct (placed inside the @c "[M+nH+...]" precursor annotation).
    @param[in]  precursor_rna_mass                        Monoisotopic mass of the complete RNA adduct.
    @param[in]  precursor_charge                          Upper bound on the emitted fragment / precursor charge states.
    @param[in]  partial_loss_modification                 Definitions of MS2 fragment adducts whose shifts the spectrum should reflect.
    @param[in]  patial_loss_template_z1                  Charge-1 template peaks to shift and append.
    @param[in]  patial_loss_template_z2                  Charge-2 template peaks (filtered to charge-2 entries during the shift step).
    @param[in]  patial_loss_template_z3                  Charge-3 template peaks (filtered to charge-3 entries during the shift step).
    @param[out] partial_loss_spectrum                     Receives the composed spectrum (data arrays resized to one each).
  */
  static void generatePartialLossSpectrum(const String& unmodified_sequence,
                                    const double& fixed_and_variable_modified_peptide_weight,
                                    const String& precursor_rna_adduct,
                                    const double& precursor_rna_mass,
                                    const int& precursor_charge,
                                    const std::vector<NuXLFragmentAdductDefinition>& partial_loss_modification,
                                    const PeakSpectrum& patial_loss_template_z1,
                                    const PeakSpectrum& patial_loss_template_z2,
                                    const PeakSpectrum& patial_loss_template_z3,
                                    PeakSpectrum& partial_loss_spectrum);

  /**
    @brief Append the @c [M+nH+complete-RNA] precursor peak to @p partial_loss_spectrum for the requested charge.

    Computes the m/z as @c (peptide_weight @c + @c precursor_rna_mass @c + @c charge @c * @c PROTON_MASS_U) @c / @c charge
    and appends one peak only when no existing entry is within @c 1e-4 m/z (collision check
    via @c findNearest). Annotation format:
      - @c "[M+H+{adduct}]" for charge @c 1,
      - @c "[M+nH+{adduct}]" with @c n explicit for charge @c >= @c 2,
    where @c {adduct} stands in for @p precursor_rna_adduct.
    Intensity @c 1.0f; the charge data array receives the actual @p charge value (not always @c 1).

    @param[in]     fixed_and_variable_modified_peptide_weight Mass of the modified peptide.
    @param[in]     precursor_rna_adduct                       Adduct label embedded in the annotation.
    @param[in]     precursor_rna_mass                         Monoisotopic mass of the complete RNA adduct.
    @param[in]     charge                                     Charge state to emit.
    @param[in,out] partial_loss_spectrum                      Spectrum to extend.
    @param[in,out] partial_loss_spectrum_charge               Parallel charge array.
    @param[in,out] partial_loss_spectrum_annotation           Parallel annotation array.
  */
  static void addPrecursorWithCompleteRNA_(const double fixed_and_variable_modified_peptide_weight,
                                    const String & precursor_rna_adduct,
                                    const double precursor_rna_mass,
                                    const int charge,
                                    PeakSpectrum & partial_loss_spectrum,
                                    MSSpectrum::IntegerDataArray & partial_loss_spectrum_charge,
                                    MSSpectrum::StringDataArray & partial_loss_spectrum_annotation);

  /**
    @brief Append two non-classical Lys-related immonium ions (@c iK(C5H10N1) and @c iK(C6H13N2O)) when @p unmodified_sequence contains @c K.

    No-op when @p unmodified_sequence has no @c K. For each variant, the peak is appended
    only when no existing entry is within @c 1e-4 m/z (collision check). Both peaks have
    charge @c 1 and intensity @c 1.0f. Unlike @ref addShiftedImmoniumIons, these are not
    mass-shifted — they capture variants observed empirically alongside the classical
    K immonium.

    @param[in]     unmodified_sequence Peptide sequence; only @c K-containing sequences contribute peaks.
    @param[in,out] spectrum            Spectrum to extend.
    @param[in,out] spectrum_charge     Parallel charge array.
    @param[in,out] spectrum_annotation Parallel annotation array.
  */
  static void addSpecialLysImmonumIons(const String& unmodified_sequence,
                                    PeakSpectrum &spectrum,
                                    PeakSpectrum::IntegerDataArray &spectrum_charge,
                                    PeakSpectrum::StringDataArray &spectrum_annotation);
  };
}
