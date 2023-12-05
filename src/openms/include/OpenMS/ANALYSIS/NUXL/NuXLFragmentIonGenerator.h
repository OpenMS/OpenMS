// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

// helper class that adds special ions not covered by TheoreticalSpectrumGenerator
class OPENMS_DLLAPI NuXLFragmentIonGenerator
{
  public:
  // prefix used to denote marker ions in fragment  annotations
  static constexpr const char* ANNOTATIONS_MARKER_ION_PREFIX = "MI:";

  // add RNA-marker ions of charge 1
  // this includes the protonated nitrogenous base and all shifts (e.g., U-H2O, U'-H20, ...)
  static void addMS2MarkerIons(
    const std::vector<NuXLFragmentAdductDefinition>& marker_ions,
    PeakSpectrum& spectrum,
    PeakSpectrum::IntegerDataArray& spectrum_charge,
    PeakSpectrum::StringDataArray& spectrum_annotation);

  static void addShiftedImmoniumIons(
    const String & unmodified_sequence,
    const String & fragment_shift_name,
    const double fragment_shift_mass,
    PeakSpectrum & partial_loss_spectrum,
    PeakSpectrum::IntegerDataArray& partial_loss_spectrum_charge,
    PeakSpectrum::StringDataArray& partial_loss_spectrum_annotation);

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
  static void addPrecursorWithCompleteRNA_(const double fixed_and_variable_modified_peptide_weight,
                                    const String & precursor_rna_adduct,
                                    const double precursor_rna_mass,
                                    const int charge,
                                    PeakSpectrum & partial_loss_spectrum,
                                    MSSpectrum::IntegerDataArray & partial_loss_spectrum_charge,
                                    MSSpectrum::StringDataArray & partial_loss_spectrum_annotation);

  static void addSpecialLysImmonumIons(const String& unmodified_sequence,
                                    PeakSpectrum &spectrum,
                                    PeakSpectrum::IntegerDataArray &spectrum_charge, 
                                    PeakSpectrum::StringDataArray &spectrum_annotation);
  };
}
