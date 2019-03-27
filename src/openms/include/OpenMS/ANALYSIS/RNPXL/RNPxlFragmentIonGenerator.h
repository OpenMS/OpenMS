// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/RNPXL/RNPxlFragmentAdductDefinition.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlFragmentAnnotationHelper.h>
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
class OPENMS_DLLAPI RNPxlFragmentIonGenerator
{
  public:
  // prefix used to denote marker ions in fragment  annotations
  static constexpr const char* ANNOTATIONS_MARKER_ION_PREFIX = "MI:";

  // add RNA-marker ions of charge 1
  // this includes the protonated nitrogenous base and all shifts (e.g., U-H2O, U'-H20, ...)
  static void addMS2MarkerIons(
    const std::vector<RNPxlFragmentAdductDefinition>& marker_ions,
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
                                    const std::vector<RNPxlFragmentAdductDefinition>& partial_loss_modification,
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

