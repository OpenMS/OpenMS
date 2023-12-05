// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentIonGenerator.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentAnnotationHelper.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{
void NuXLFragmentIonGenerator::addMS2MarkerIons(
  const std::vector<NuXLFragmentAdductDefinition> &marker_ions, 
  PeakSpectrum &spectrum,
  PeakSpectrum::IntegerDataArray &spectrum_charge, 
  PeakSpectrum::StringDataArray &spectrum_annotation)
{
  for (auto const & m : marker_ions)
  {
    const double mz = m.mass + Constants::PROTON_MASS_U;

    spectrum.emplace_back(mz, 1.0);
    spectrum_charge.emplace_back(1);
    spectrum_annotation.emplace_back(NuXLFragmentIonGenerator::ANNOTATIONS_MARKER_ION_PREFIX + m.name);  // add name (e.g., MI:U-H2O)
  }
}

void NuXLFragmentIonGenerator::addSpecialLysImmonumIons(
  const String& unmodified_sequence,
  PeakSpectrum &spectrum,
  PeakSpectrum::IntegerDataArray &spectrum_charge, 
  PeakSpectrum::StringDataArray &spectrum_annotation)
{
   if (unmodified_sequence.has('K'))
   {
      const double immonium_ion2_mz = EmpiricalFormula("C5H10N1").getMonoWeight(); 
      // only add special ios if there is not already a peak
      if (spectrum.findNearest(immonium_ion2_mz, 1e-4) == -1)
      {  
        spectrum.emplace_back(immonium_ion2_mz, 1.0);
        spectrum_charge.emplace_back(1);
        spectrum_annotation.emplace_back(String("iK(C5H10N1)"));
      }

      // usually only observed without shift (A. Stuetzer)
      const double immonium_ion3_mz = EmpiricalFormula("C6H13N2O").getMonoWeight(); 
      // only add special ios if there is not already a peak
      if (spectrum.findNearest(immonium_ion3_mz, 1e-4) == -1)
      {  
        spectrum.emplace_back(immonium_ion3_mz, 1.0);
        spectrum_charge.emplace_back(1);
        spectrum_annotation.emplace_back(String("iK(C6H13N2O)"));
      }

    }
}


void NuXLFragmentIonGenerator::addShiftedImmoniumIons(const String &unmodified_sequence,
                                                                    const String &fragment_shift_name,
                                                                    const double fragment_shift_mass,
                                                                    PeakSpectrum &partial_loss_spectrum,
                                                                    PeakSpectrum::IntegerDataArray &partial_loss_spectrum_charge,
                                                                    PeakSpectrum::StringDataArray &partial_loss_spectrum_annotation) 
{
  if (unmodified_sequence.hasSubstring("Y"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C8H10NO").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('Y', fragment_shift_name));
  }

  if (unmodified_sequence.hasSubstring("W"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C10H11N2").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('W', fragment_shift_name));
  }

  if (unmodified_sequence.hasSubstring("F"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C8H10N").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('F', fragment_shift_name));
  }

  if (unmodified_sequence.hasSubstring("H"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C5H8N3").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('H', fragment_shift_name));
  }

  if (unmodified_sequence.hasSubstring("C"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C2H6NS").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('C', fragment_shift_name));
  }

  if (unmodified_sequence.hasSubstring("P"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C4H8N").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('P', fragment_shift_name));
  }

  if (unmodified_sequence.hasSubstring("L") || unmodified_sequence.hasSubstring("I"))
  {
    const double immonium_ion_mz = EmpiricalFormula("C5H12N").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('L', fragment_shift_name));
  }

  if (unmodified_sequence.hasSubstring("K"))
  {
    // classical immonium ion
    const double immonium_ion_mz = EmpiricalFormula("C5H13N2").getMonoWeight() + fragment_shift_mass;
    partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('K', fragment_shift_name));

    // TODO: check if only DNA specific and if also other shifts are observed
    // according to A. Stuetzer mainly observed with C‘-NH3 (94.0167 Da)
    const double immonium_ion2_mz = EmpiricalFormula("C5H10N1").getMonoWeight()  + fragment_shift_mass; 
    partial_loss_spectrum.emplace_back(immonium_ion2_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(String("iK(C5H10N1)" + fragment_shift_name));

    // usually only observed without shift (A. Stuetzer)
    const double immonium_ion3_mz = EmpiricalFormula("C6H13N2O").getMonoWeight()  + fragment_shift_mass; 
    partial_loss_spectrum.emplace_back(immonium_ion3_mz, 1.0);
    partial_loss_spectrum_charge.emplace_back(1);
    partial_loss_spectrum_annotation.emplace_back(String("iK(C6H13N2O)" + fragment_shift_name));
  }

  if (unmodified_sequence.hasSubstring("M"))
  {
    {
      const double immonium_ion_mz = 104.05285 + fragment_shift_mass;
      partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
      partial_loss_spectrum_charge.emplace_back(1);
      partial_loss_spectrum_annotation.emplace_back(NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('M', fragment_shift_name));
    }

    {
      const double immonium_ion_mz = EmpiricalFormula("CH5S").getMonoWeight() + fragment_shift_mass; // methionine related fragment
      partial_loss_spectrum.emplace_back(immonium_ion_mz, 1.0);
      partial_loss_spectrum_charge.emplace_back(1);
      partial_loss_spectrum_annotation.emplace_back(NuXLFragmentAnnotationHelper::getAnnotatedImmoniumIon('M', fragment_shift_name));
    }
  }
}

  /* 
  * Add peaks with shifts induced by the RNA/DNA:
  *   - Precursor with complete NA-oligo for charge 1..z
  *   - Partial shifts (without complete precursor adduct)
  *     - Add shifted immonium ions for charge 1 only
  *     and create shifted shifted b,y,a ions + precursors for charge 1..z (adding the unshifted version and performing the shift)
  *     based on the total_loss_spectrum provide to the method
  */
void NuXLFragmentIonGenerator::generatePartialLossSpectrum(const String &unmodified_sequence,
                                                                         const double &fixed_and_variable_modified_peptide_weight,
                                                                         const String &precursor_rna_adduct,
                                                                         const double &precursor_rna_mass,
                                                                         const int &precursor_charge,
                                                                         const std::vector<NuXLFragmentAdductDefinition> &partial_loss_modification,
                                                                         const PeakSpectrum& partial_loss_template_z1,
                                                                         const PeakSpectrum& partial_loss_template_z2,
                                                                         const PeakSpectrum& partial_loss_template_z3,
                                                                         PeakSpectrum &partial_loss_spectrum)
{
  partial_loss_spectrum.getIntegerDataArrays().resize(1);
  PeakSpectrum::IntegerDataArray& partial_loss_spectrum_charge = partial_loss_spectrum.getIntegerDataArrays()[0];

  partial_loss_spectrum.getStringDataArrays().resize(1); // annotation
  PeakSpectrum::StringDataArray& partial_loss_spectrum_annotation = partial_loss_spectrum.getStringDataArrays()[0];

  // for all observable MS2 adducts ...
  for (Size i = 0; i != partial_loss_modification.size(); ++i)
  {
    // get name and mass of fragment adduct
    const String& fragment_shift_name = partial_loss_modification[i].name; // e.g. U-H2O
    const double fragment_shift_mass = partial_loss_modification[i].mass;

    // ADD: shifted immonium ion peaks of charge 1 (if the amino acid is present in the sequence)
    NuXLFragmentIonGenerator::addShiftedImmoniumIons(
      unmodified_sequence,
      fragment_shift_name,
      fragment_shift_mass,
      partial_loss_spectrum,
      partial_loss_spectrum_charge,
      partial_loss_spectrum_annotation);

    // annotate generated a-,b-,y-ions with fragment shift name
    PeakSpectrum shifted_series_peaks;
    shifted_series_peaks.getStringDataArrays().resize(1); // annotation
    shifted_series_peaks.getIntegerDataArrays().resize(1); // charge
    PeakSpectrum::StringDataArray& shifted_series_annotations = shifted_series_peaks.getStringDataArrays()[0];
    PeakSpectrum::IntegerDataArray& shifted_series_charges = shifted_series_peaks.getIntegerDataArrays()[0];

    // For every charge state
    for (int z = 1; z <= precursor_charge; ++z)
    {
      // 1. add shifted peaks 
      if (z == 1)
      {
        for (Size i = 0; i != partial_loss_template_z1.size(); ++i) 
        { 
          Peak1D p = partial_loss_template_z1[i];
          p.setMZ(p.getMZ() + fragment_shift_mass);         
          shifted_series_peaks.push_back(p);
          shifted_series_annotations.push_back(partial_loss_template_z1.getStringDataArrays()[0][i]);
          shifted_series_charges.push_back(1);
        } 
      }
      else if (z == 2)
      {
        for (Size i = 0; i != partial_loss_template_z2.size(); ++i) 
        { 
          // currently, also contains z=1 precursor peaks which we aleardy added before
          if (partial_loss_template_z2.getIntegerDataArrays()[0][i] == 2)
          { 
            Peak1D p = partial_loss_template_z2[i];
            p.setMZ(p.getMZ() + fragment_shift_mass / 2.0);         
            shifted_series_peaks.push_back(p);
            shifted_series_annotations.push_back(partial_loss_template_z2.getStringDataArrays()[0][i]);
            shifted_series_charges.push_back(2);
          } 
        } 
      }
      else if (z == 3)
      {
        for (Size i = 0; i != partial_loss_template_z3.size(); ++i) 
        { 
          // currently, also contains z=1 and 2 precursor peaks which we aleardy added before
          if (partial_loss_template_z3.getIntegerDataArrays()[0][i] == 3)
          { 
            Peak1D p = partial_loss_template_z3[i];
            p.setMZ(p.getMZ() + fragment_shift_mass / 3.0);         
            shifted_series_peaks.push_back(p);
            shifted_series_annotations.push_back(partial_loss_template_z3.getStringDataArrays()[0][i]);
            shifted_series_charges.push_back(3);
          } 
        } 
      }
      else // don't consider fragment ions with charge >= 4 
      { 
        break; 
      }    
    }

    // 2. add fragment shift name to annotation of shifted peaks
    for (Size j = 0; j != shifted_series_annotations.size(); ++j)
    {
      shifted_series_annotations[j] += " " + fragment_shift_name;
    }

    // append shifted and annotated ion series to partial loss spectrum
    partial_loss_spectrum.insert(partial_loss_spectrum.end(), 
      shifted_series_peaks.begin(), shifted_series_peaks.end());
    // std::move strings during insert
    partial_loss_spectrum_annotation.insert(
      partial_loss_spectrum_annotation.end(),
      make_move_iterator(shifted_series_annotations.begin()),
      make_move_iterator(shifted_series_annotations.end())
    );
    partial_loss_spectrum.getIntegerDataArrays()[0].insert(
      partial_loss_spectrum_charge.end(),
      shifted_series_charges.begin(),
      shifted_series_charges.end()
    );
  }

  // ADD: (mainly for ETD) MS2 precursor peaks of the MS1 adduct (total RNA) carrying peptide for all z <= precursor charge
  for (int charge = 1; charge <= static_cast<int>(precursor_charge); ++charge)
  {
    addPrecursorWithCompleteRNA_(fixed_and_variable_modified_peptide_weight,
                                 precursor_rna_adduct,
                                 precursor_rna_mass,
                                 charge,
                                 partial_loss_spectrum,
                                 partial_loss_spectrum_charge,
                                 partial_loss_spectrum_annotation);
  }

  partial_loss_spectrum.sortByPosition();
}

void NuXLFragmentIonGenerator::addPrecursorWithCompleteRNA_(
  const double fixed_and_variable_modified_peptide_weight, 
  const String &precursor_rna_adduct,
  const double precursor_rna_mass, 
  const int charge, 
  PeakSpectrum &partial_loss_spectrum,
  MSSpectrum::IntegerDataArray &partial_loss_spectrum_charge,
  MSSpectrum::StringDataArray &partial_loss_spectrum_annotation)
{
  const double xl_mz = (fixed_and_variable_modified_peptide_weight + precursor_rna_mass +
                  static_cast<double>(charge) * Constants::PROTON_MASS_U)
                 / static_cast<double>(charge);

  // only add special ions if there is not already a peak
  if (partial_loss_spectrum.findNearest(xl_mz, 1e-4) == -1)
  {  
    partial_loss_spectrum.push_back(Peak1D(xl_mz, 1.0));
    partial_loss_spectrum_charge.push_back(charge);
    if (charge > 1)
    {
      partial_loss_spectrum_annotation.push_back(String("[M+") 
        + String(charge) + "H+" + precursor_rna_adduct + "]");
    } 
    else
    {
      partial_loss_spectrum_annotation.push_back(String("[M+H+") 
        + precursor_rna_adduct + "]");
    }  
  }
}

}

