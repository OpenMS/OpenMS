// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentIonGenerator.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLAnnotateAndLocate.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLFragmentAnnotationHelper.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLConstants.h>

#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Macros.h> // for OPENMS_PRECONDITION
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

//#define DEBUG_OpenNuXL

namespace OpenMS
{

  // create total loss spectrum using new_param as template
  static PeakSpectrum createTotalLossSpectrumForAnnotations(const AASequence& fixed_and_variable_modified_peptide, size_t precursor_charge, Param new_param)
  {
    PeakSpectrum total_loss_spectrum;
    TheoreticalSpectrumGenerator tmp_generator;
    new_param.setValue("add_all_precursor_charges", "true");
    new_param.setValue("add_abundant_immonium_ions", "true");
    new_param.setValue("add_losses", "true");
    new_param.setValue("add_term_losses", "true");
    new_param.setValue("add_a_ions", "true");
    new_param.setValue("add_internal_fragments", "true");
    tmp_generator.setParameters(new_param);
    tmp_generator.getSpectrum(total_loss_spectrum, fixed_and_variable_modified_peptide, 1, precursor_charge);

    const String& unmodified_sequence = fixed_and_variable_modified_peptide.toUnmodifiedString();
    const bool contains_Methionine = unmodified_sequence.has('M');

    if (contains_Methionine) // add mainly DEB + NM related precursor losses 
    {
      static const double M_star_pc_loss = EmpiricalFormula("CH4S").getMonoWeight(); // methionine related loss on precursor (see OpenNuXL for scoring related code)
      for (size_t charge = 1; charge <= precursor_charge; ++charge)
      {
        String ion_name = (charge == 1) ? "[M+H]-CH4S" : "[M+" + String(charge) + "H]-CH4S";              
        total_loss_spectrum.getStringDataArrays()[0].push_back(ion_name);
        total_loss_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX].push_back(charge);      
        double mono_pos = fixed_and_variable_modified_peptide.getMonoWeight(Residue::Full, charge) - M_star_pc_loss; // precursor peak
        total_loss_spectrum.emplace_back(mono_pos / (double)charge, 1.0);
      }
    }
    // add special immonium ions
    NuXLFragmentIonGenerator::addSpecialLysImmonumIons(
      unmodified_sequence,
      total_loss_spectrum, 
      total_loss_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
      total_loss_spectrum.getStringDataArrays()[0]);
    total_loss_spectrum.sortByPosition(); // need to resort after adding special immonium ions
    return total_loss_spectrum;
  }

  using MapIonIndexToFragmentAnnotation = map<Size, vector<NuXLFragmentAnnotationHelper::FragmentAnnotationDetail_> >;

  // ion centric (e.g. b and y-ion) spectrum annotation for unshifted ions (will later be merged with shifted) 
  static vector<PeptideHit::PeakAnnotation> createIonCentricFragmentAnnotationsForUnshiftedIons(
    const PeakSpectrum& total_loss_spectrum, 
    const PeakSpectrum& exp_spectrum, 
    const vector<pair<Size, Size>>& alignment, 
    set<Size>& peak_is_annotated,
    vector<PeptideHit::PeakAnnotation>& annotated_precursor_ions,
    MapIonIndexToFragmentAnnotation& unshifted_b_ions, 
    MapIonIndexToFragmentAnnotation& unshifted_y_ions, 
    MapIonIndexToFragmentAnnotation& unshifted_a_ions,
    vector<PeptideHit::PeakAnnotation>& unshifted_loss_ions, 
    vector<PeptideHit::PeakAnnotation>& annotated_immonium_ions
    )
  {
    const PeakSpectrum::StringDataArray& total_loss_annotations = total_loss_spectrum.getStringDataArrays()[0];
    const PeakSpectrum::IntegerDataArray& total_loss_charges = total_loss_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX];

    // create total loss annotations
    for (auto const & aligned : alignment)
    {
      // information on the experimental fragment in the alignment
      const Size& fragment_index = aligned.second;
      const Peak1D& fragment = exp_spectrum[fragment_index];
      const double fragment_intensity = fragment.getIntensity(); // in percent (%)
      const double fragment_mz = fragment.getMZ();
      

      const String& ion_name = total_loss_annotations[aligned.first];
      const int charge = total_loss_charges[aligned.first];

      OPENMS_PRECONDITION(exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX][fragment_index] == charge, "Charges in alignment must match.");

      // define which ion names are annotated
      if (ion_name[0] == 'y')
      {
        Size loss_first = ion_name.find_first_of('-'); // start of loss
        Size charge_pos = ion_name.find_first_of('+'); // charge indicator at end
        const bool ion_has_neutral_loss = (loss_first != string::npos);

        if (ion_has_neutral_loss) // ion with neutral loss e.g. water
        {
          PeptideHit::PeakAnnotation fa;
          fa.mz = fragment_mz;
          fa.intensity = fragment_intensity;
          fa.charge = charge;
          fa.annotation = ion_name;
          unshifted_loss_ions.push_back(fa);
          peak_is_annotated.insert(aligned.second);
        }
        else // no neutral loss
        {
          String ion_nr_string = ion_name.substr(1, charge_pos - 1);
          Size ion_number = (Size)ion_nr_string.toInt();
          NuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
          unshifted_y_ions[ion_number].push_back(d);
          #ifdef DEBUG_OpenNuXL
            const AASequence& peptide_sequence = fixed_and_variable_modified_peptide.getSuffix(ion_number);
            OPENMS_LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
          #endif
          peak_is_annotated.insert(aligned.second);
        }
      }
      else if (ion_name[0] == 'b')
      {
        Size loss_first = ion_name.find_first_of('-'); // start of loss
        Size charge_pos = ion_name.find_first_of('+'); // charge indicator at end
        const bool ion_has_neutral_loss = (loss_first != string::npos);

        if (ion_has_neutral_loss)
        {
          PeptideHit::PeakAnnotation fa;
          fa.mz = fragment_mz;
          fa.intensity = fragment_intensity;
          fa.charge = charge;
          fa.annotation = ion_name;
          unshifted_loss_ions.push_back(fa);
          peak_is_annotated.insert(aligned.second);
        }
        else
        {
          String ion_nr_string = ion_name.substr(1, charge_pos - 1);
          Size ion_number = (Size)ion_nr_string.toInt();
          #ifdef DEBUG_OpenNuXL
            const AASequence& peptide_sequence = aas.getPrefix(ion_number);
            OPENMS_LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
          #endif
          NuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
          unshifted_b_ions[ion_number].push_back(d);
          peak_is_annotated.insert(aligned.second);
        }
      }
      else if (ion_name[0] == 'a')
      {
        Size loss_first = ion_name.find_first_of('-'); // start of loss
        Size charge_pos = ion_name.find_first_of('+'); // charge indicator at end
        const bool ion_has_neutral_loss = (loss_first != string::npos);

        if (ion_has_neutral_loss)
        {
          PeptideHit::PeakAnnotation fa;
          fa.mz = fragment_mz;
          fa.intensity = fragment_intensity;
          fa.charge = charge;
          fa.annotation = ion_name;
          unshifted_loss_ions.push_back(fa);
          peak_is_annotated.insert(aligned.second);
        }
        else
        {
          String ion_nr_string = ion_name.substr(1, charge_pos - 1);
          auto ion_number = (Size)ion_nr_string.toInt();
          #ifdef DEBUG_OpenNuXL
            const AASequence& peptide_sequence = aas.getPrefix(ion_number);
            OPENMS_LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << peptide_sequence.toString() << " intensity: " << fragment_intensity << endl;
          #endif
          NuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d("", charge, fragment_mz, fragment_intensity);
          unshifted_a_ions[ion_number].push_back(d);
          peak_is_annotated.insert(aligned.second);
        }
      }
      else if (ion_name.hasPrefix("[M+")) // precursor ion
      {
        PeptideHit::PeakAnnotation fa;
        fa.mz = fragment_mz;
        fa.intensity = fragment_intensity;
        fa.charge = charge;
        fa.annotation = ion_name;
        annotated_precursor_ions.push_back(fa);
        peak_is_annotated.insert(aligned.second);
      }
      else if (ion_name.hasPrefix("i")) // immonium ion
      {
        PeptideHit::PeakAnnotation fa;
        fa.mz = fragment_mz;
        fa.intensity = fragment_intensity;
        fa.charge = charge;
        fa.annotation = ion_name;
        annotated_immonium_ions.push_back(fa);
        peak_is_annotated.insert(aligned.second);
      }
      else if (isupper(ion_name[0])) // internal ions
      {
        PeptideHit::PeakAnnotation fa;
        fa.mz = fragment_mz;
        fa.intensity = fragment_intensity;
        fa.charge = charge;
        fa.annotation = ion_name;
        annotated_immonium_ions.push_back(fa);  //TODO: add to annotated_internal_fragment_ions or rename vector
        peak_is_annotated.insert(aligned.second);
      }
    }

    // generate fragment annotation strings for unshifted ions
    vector<PeptideHit::PeakAnnotation> fas;
    if (!unshifted_b_ions.empty())
    {
      const vector<PeptideHit::PeakAnnotation>& fas_tmp = NuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToPHFA("b", unshifted_b_ions);
      fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
    }
    if (!unshifted_y_ions.empty())
    {
      const vector<PeptideHit::PeakAnnotation>& fas_tmp = NuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToPHFA("y", unshifted_y_ions);
      fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
    }
    if (!unshifted_a_ions.empty())
    {
      const vector<PeptideHit::PeakAnnotation>& fas_tmp = NuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToPHFA("a", unshifted_a_ions);
      fas.insert(fas.end(), fas_tmp.begin(), fas_tmp.end());
    }
    if (!annotated_immonium_ions.empty())
    {
      fas.insert(fas.end(), annotated_immonium_ions.begin(), annotated_immonium_ions.end());          
    }
    if (!unshifted_loss_ions.empty())
    {
      fas.insert(fas.end(), unshifted_loss_ions.begin(), unshifted_loss_ions.end());          
    }
    return fas;
  }

  // static
  void NuXLAnnotateAndLocate::annotateAndLocate_(
    const PeakMap& exp, 
    vector<vector<NuXLAnnotatedHit>>& annotated_hits,
    const NuXLModificationMassesResult& mm,
    const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications, 
    const ModifiedPeptideGenerator::MapToResidueType& variable_modifications, 
    Size max_variable_mods_per_peptide, 
    double fragment_mass_tolerance, 
    bool fragment_mass_tolerance_unit_ppm, 
    const NuXLParameterParsing::PrecursorsToMS2Adducts & all_feasible_adducts)   
  {
    TheoreticalSpectrumGenerator partial_loss_spectrum_generator;
    auto param = partial_loss_spectrum_generator.getParameters();
    param.setValue("add_first_prefix_ion", "true");
    param.setValue("add_abundant_immonium_ions", "false"); // we add them manually for charge 1
    param.setValue("add_precursor_peaks", "true");
    param.setValue("add_all_precursor_charges", "false"); // we add them manually for every charge
    param.setValue("add_metainfo", "true");
    param.setValue("add_a_ions", "true");
    param.setValue("add_b_ions", "true");
    param.setValue("add_c_ions", "false");
    param.setValue("add_x_ions", "false");
    param.setValue("add_y_ions", "true");
    param.setValue("add_z_ions", "false");
    param.setValue("add_internal_fragments", "true"); 
    partial_loss_spectrum_generator.setParameters(param);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (SignedSize scan_index = 0; scan_index < (SignedSize)annotated_hits.size(); ++scan_index)
    {
      if (annotated_hits[scan_index].empty()) { continue; }

      const PeakSpectrum & exp_spectrum = exp[scan_index];
      const Size & precursor_charge = exp_spectrum.getPrecursors()[0].getCharge();

      for (auto & a : annotated_hits[scan_index])
      {
        // get unmodified string
        const String unmodified_sequence = a.sequence.getString();

        // initialize result fields
        a.best_localization = unmodified_sequence;
        a.best_localization_score = 0;

        AASequence aas(AASequence::fromString(unmodified_sequence));

        // reapply modifications (because for memory reasons we only stored the index and recreation is fast)
        vector<AASequence> all_modified_peptides;
        ModifiedPeptideGenerator::applyFixedModifications(fixed_modifications, aas);
        ModifiedPeptideGenerator::applyVariableModifications(variable_modifications, aas, max_variable_mods_per_peptide, all_modified_peptides);

        // sequence with modifications - note: reannotated version requires much more memory heavy AASequence object
        const AASequence& fixed_and_variable_modified_peptide = all_modified_peptides[a.peptide_mod_index];
        const double fixed_and_variable_modified_peptide_weight = fixed_and_variable_modified_peptide.getMonoWeight();

        // determine NA on precursor from index in map
        auto mod_combinations_it = mm.mod_combinations.cbegin();
        std::advance(mod_combinations_it, a.NA_mod_index);
        const String precursor_na_adduct = *mod_combinations_it->second.begin(); // TODO: check if it is enough to consider only first precursor adduct ????????????????????????????????????????????????????????
        const double precursor_na_mass = EmpiricalFormula(mod_combinations_it->first).getMonoWeight();

        // generate total loss spectrum for the fixed and variable modified peptide (without NAs) (using the settings for partial loss generation)
        // but as we also add the abundant immonium ions for charge 1 and precursor ions for all charges to get a more complete annotation
        // (these have previously not been used in the scoring of the total loss spectrum)
        PeakSpectrum total_loss_spectrum = createTotalLossSpectrumForAnnotations(fixed_and_variable_modified_peptide, precursor_charge, partial_loss_spectrum_generator.getParameters()); // use same parameters

        // first annotate total loss peaks (these give no information where the actual shift occured)
        #ifdef DEBUG_OpenNuXL
          OPENMS_LOG_DEBUG << "Annotating ion (total loss spectrum): " << fixed_and_variable_modified_peptide.toString()  << endl;
        #endif
        vector<pair<Size, Size>> alignment;

        // align spectra (only allow matching charges)
        DataArrays::FloatDataArray ppm_error_array; // not needed here but filled by alignment
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(alignment, 
          fragment_mass_tolerance, 
          fragment_mass_tolerance_unit_ppm, 
          total_loss_spectrum, 
          exp_spectrum, 
          total_loss_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX], 
          exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX], 
          ppm_error_array);

        // fill annotated spectrum information
        set<Size> peak_is_annotated;  // experimental peak index
        vector<PeptideHit::PeakAnnotation> annotated_precursor_ions; // also used for shifted ones

        MapIonIndexToFragmentAnnotation unshifted_b_ions, unshifted_y_ions, unshifted_a_ions;
        vector<PeptideHit::PeakAnnotation> unshifted_loss_ions, annotated_immonium_ions;

        auto fas = createIonCentricFragmentAnnotationsForUnshiftedIons(total_loss_spectrum, exp_spectrum, alignment, peak_is_annotated, 
          annotated_precursor_ions,
          unshifted_b_ions, 
          unshifted_y_ions, 
          unshifted_a_ions,
          unshifted_loss_ions, 
          annotated_immonium_ions        
          );
            
        // we don't localize on non-cross-links (only annotate)
        if (precursor_na_adduct == "none") 
        { 
          a.fragment_annotations = fas;
          continue; 
        }

        // ion centric (e.g. b and y-ion) spectrum annotation that records all shifts of specific ions (e.g. y5, y5 + U, y5 + C3O)

        // generate all partial loss spectra (excluding the complete loss spectrum) merged into one spectrum
        // 1. get all possible NA fragment shifts in the MS2 (based on the precursor RNA/DNA)
        OPENMS_LOG_DEBUG << "Precursor NA adduct: "  << precursor_na_adduct << endl;

        const vector<NucleotideToFeasibleFragmentAdducts>& feasible_MS2_adducts = all_feasible_adducts.at(precursor_na_adduct).feasible_adducts;

        if (feasible_MS2_adducts.empty()) { continue; } // should not be the case - check case of no nucleotide but base fragment ?

        // 2. retrieve the (nucleotide specific) fragment adducts for the cross-linked nucleotide (annotated in main search)
        auto nt_to_adducts = std::find_if(feasible_MS2_adducts.begin(), feasible_MS2_adducts.end(),
          [&a](NucleotideToFeasibleFragmentAdducts const & item)
          {
            return (item.first == a.cross_linked_nucleotide);
          });

        OPENMS_POSTCONDITION(nt_to_adducts != feasible_MS2_adducts.end(), "Nucleotide not found in mapping to feasible adducts.")

        const vector<NuXLFragmentAdductDefinition>& partial_loss_modification = nt_to_adducts->second;

        // get marker ions (these are not specific to the cross-linked nucleotide but also depend on the whole oligo bound to the precursor)
        const vector<NuXLFragmentAdductDefinition>& marker_ions = all_feasible_adducts.at(precursor_na_adduct).marker_ions;
        OPENMS_LOG_DEBUG << "Marker ions used for this Precursor NA adduct: "  << endl;
        for (auto & fa : marker_ions)
        {
          OPENMS_LOG_DEBUG << fa.name << " " << fa.mass << endl;
        }


        PeakSpectrum partial_loss_spectrum;

        {
          PeakSpectrum partial_loss_template_z1, 
                      partial_loss_template_z2, 
                      partial_loss_template_z3;
       
          partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z1, fixed_and_variable_modified_peptide, 1, 1);
          partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z2, fixed_and_variable_modified_peptide, 2, 2);
          partial_loss_spectrum_generator.getSpectrum(partial_loss_template_z3, fixed_and_variable_modified_peptide, 3, 3);
          NuXLFragmentIonGenerator::generatePartialLossSpectrum(unmodified_sequence,
                                      fixed_and_variable_modified_peptide_weight,
                                      precursor_na_adduct,
                                      precursor_na_mass,
                                      precursor_charge,
                                      partial_loss_modification,
                                      partial_loss_template_z1,
                                      partial_loss_template_z2,
                                      partial_loss_template_z3,
                                      partial_loss_spectrum);
        }

        // add shifted marker ions
        NuXLFragmentIonGenerator::addMS2MarkerIons(
          marker_ions,
          partial_loss_spectrum,
          partial_loss_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX],
          partial_loss_spectrum.getStringDataArrays()[0]);

        partial_loss_spectrum.sortByPosition(); // need to resort after adding marker ions
        
        // ion centric (e.g. b and y-ion) spectrum annotation that records all shifts of specific ions (e.g. y5, y5 + U, y5 + C3O)
        MapIonIndexToFragmentAnnotation shifted_b_ions, shifted_y_ions, shifted_a_ions;
        vector<PeptideHit::PeakAnnotation> shifted_immonium_ions, annotated_marker_ions;

        vector<double> sites_sum_score(aas.size(), 0);

        /////////////////
        // Align partial-loss-spectrum to the experimental measured one
        alignment.clear();
        ppm_error_array.clear();

        // align spectra (only allow matching charges)
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(alignment, 
          fragment_mass_tolerance, 
          fragment_mass_tolerance_unit_ppm, 
          partial_loss_spectrum, 
          exp_spectrum, 
          partial_loss_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX], 
          exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX], 
          ppm_error_array);

        #ifdef OPENUXL_DEBUG
        for (size_t i = 0; i != exp_spectrum.size(); ++i)
        {
          OPENMS_LOG_DEBUG << "exp: " << exp_spectrum[i].getMZ() << "\t" << exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX][i] << endl;
        }
        #endif

        const PeakSpectrum::StringDataArray& partial_loss_annotations = partial_loss_spectrum.getStringDataArrays()[0];
        const PeakSpectrum::IntegerDataArray& partial_loss_charges = partial_loss_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX];

        if (alignment.empty())
        {
          a.fragment_annotations = fas;
          continue;
        }

        /* uncomment to write all annotations to a file(only makes sense if a single spectrum is searched)
        MSExperiment tmp_exp;
        tmp_exp.addSpectrum(total_loss_spectrum);
        tmp_exp.addSpectrum(partial_loss_spectrum);
        MzMLFile().store("theoretical_loss_spectrum.mzML", tmp_exp);
        */

        for (auto pair_it = alignment.begin(); pair_it != alignment.end(); ++pair_it)
        {
          // only annotate experimental peaks with shift - i.e. do not annotated complete loss peaks again
          if (peak_is_annotated.find(pair_it->second) != peak_is_annotated.end()) { continue; }

          // information on the experimental fragment in the alignment
          const Size & fragment_index = pair_it->second;
          const Peak1D & fragment = exp_spectrum[fragment_index];
          const double & fragment_intensity = fragment.getIntensity(); // in percent (%)
          const double & fragment_mz = fragment.getMZ();
          const int & fragment_charge = exp_spectrum.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX][fragment_index];
          
          #ifdef DEBUG_OpenNuXL
            OPENMS_LOG_DEBUG << "fragment_mz:" << fragment_mz << " fragment_charge:" << fragment_charge << endl; 
          #endif

          String ion_name = partial_loss_annotations[pair_it->first];
          const int charge = partial_loss_charges[pair_it->first];

          #ifdef DEBUG_OpenNuXL
            OPENMS_LOG_DEBUG << "theo_name:" << ion_name  << " theo_charge:" << charge << endl; 
          #endif
          vector<String> f;

          ion_name.split(' ', f);  // e.g. "y3 C3O" or just "y2"
          String fragment_shift_name;
          if (f.size() == 2) { fragment_shift_name = f[1]; }

          String fragment_ion_name = f[0]; // e.g. y3

          #ifdef DEBUG_OpenNuXL
            OPENMS_LOG_DEBUG << "Annotating ion: " << ion_name << " at position: " << fragment_mz << " " << " intensity: " << fragment_intensity << endl;
          #endif

          // define which ion names are annotated
          if (fragment_ion_name.hasPrefix("y"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("y", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            NuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
            if (ion_number > 1) // trypsin doesn't cut at cross-linked amino acid
            {
              shifted_y_ions[ion_number].push_back(d);
            }
          }
          else if (fragment_ion_name.hasPrefix("b"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("b", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            NuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
            shifted_b_ions[ion_number].push_back(d);
          }
          else if (fragment_ion_name.hasPrefix("a"))
          {
            String ion_nr_string = fragment_ion_name;
            ion_nr_string.substitute("a", "");
            ion_nr_string.substitute("+", ""); // remove one or multiple '+'
            auto ion_number = (Size)ion_nr_string.toInt();

            NuXLFragmentAnnotationHelper::FragmentAnnotationDetail_ d(fragment_shift_name, charge, fragment_mz, fragment_intensity);
            shifted_a_ions[ion_number].push_back(d);
          }
          else if (ion_name.hasPrefix(NuXLFragmentIonGenerator::ANNOTATIONS_MARKER_ION_PREFIX))
          {
            OPENMS_LOG_DEBUG << "Marker ion aligned: " << ion_name << " fragment_mz: " << fragment_mz << " fragment_charge: " << fragment_charge << endl;
            if (fragment_charge == 1)
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = 1;
              fa.annotation = ion_name;
              annotated_marker_ions.push_back(fa);
            }
            else
            {
              OPENMS_LOG_ERROR << "Unexpected marker ion charge." << endl;
            }            
          }
          else if (ion_name.hasPrefix("i"))
          {
            OPENMS_LOG_DEBUG << "Immonium ion aligned: " << ion_name << " fragment_mz: " << fragment_mz << " fragment_charge: " << fragment_charge << endl;            
            if (fragment_charge == 1)
            {
              PeptideHit::PeakAnnotation fa;
              fa.mz = fragment_mz;
              fa.intensity = fragment_intensity;
              fa.charge = 1;
              fa.annotation = ion_name;
              shifted_immonium_ions.push_back(fa);
            }
            else
            {
              OPENMS_LOG_ERROR << "Unexpected immonium ion charge." << endl;
            }
          }
          else if (ion_name.hasPrefix("[M+"))
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = charge;
            fa.annotation = ion_name;  
            annotated_precursor_ions.push_back(fa);
          }
          else if (isupper(ion_name[0])) // shifted internal ions
          {
            PeptideHit::PeakAnnotation fa;
            fa.mz = fragment_mz;
            fa.intensity = fragment_intensity;
            fa.charge = charge;
            String with_plus = ion_name;
            with_plus.substitute(' ', '+');
            fa.annotation = with_plus; // turn "PEPT U-H2O" into "PEPT+U-H20"
            shifted_immonium_ions.push_back(fa);  //TODO: add to shifted_internal_fragment_ions or rename vector
          }
        }

        // track shifts in n- and c-term ladders (in AAs coordinates)
        // n_shifts and c_shifts will contain the summed intensities over all observed shifts at that position
        // the distinction allows to easily detect prefix and suffix ladders in the next step
        vector<double> n_shifts(sites_sum_score.size(), 0); // vector index 0 == ion index 1
        vector<double> c_shifts(sites_sum_score.size(), 0);

        for (Size i = 0; i != n_shifts.size(); ++i)
        {
          if (shifted_b_ions.find(i + 1) == shifted_b_ions.end()) { continue; }
          for (auto& k : shifted_b_ions[i + 1]) { n_shifts[i] += k.intensity; }
        }

        for (Size i = 0; i != n_shifts.size(); ++i)
        {
          if (shifted_a_ions.find(i + 1) == shifted_a_ions.end()) { continue; }
          for (auto& k : shifted_a_ions[i + 1]) { n_shifts[i] += k.intensity; }
        }

        for (Size i = 0; i != c_shifts.size(); ++i)
        {
          const Size ion_index = c_shifts.size() - i;
          if (shifted_y_ions.find(ion_index) == shifted_y_ions.end()) { continue; }
          for (auto& k : shifted_y_ions[ion_index]) { c_shifts[i] += k.intensity; }
        }

        vector<double> n_noshifts(sites_sum_score.size(), 0);
        vector<double> c_noshifts(sites_sum_score.size(), 0);
        for (Size i = 0; i != n_noshifts.size(); ++i)
        {
          if (unshifted_b_ions.find(i + 1) == unshifted_b_ions.end()) { continue; }
          for (auto& k : unshifted_b_ions[i + 1]) { n_noshifts[i] += k.intensity; }
        }

        for (Size i = 0; i != n_noshifts.size(); ++i)
        {
          if (unshifted_a_ions.find(i + 1) == unshifted_a_ions.end()) { continue; }
          for (auto& k : unshifted_a_ions[i + 1]) { n_noshifts[i] += k.intensity; }
        }

        for (Size i = 0; i != c_noshifts.size(); ++i)
        {
          const Size ion_index = c_noshifts.size() - i;
          if (unshifted_y_ions.find(ion_index) == unshifted_y_ions.end()) { continue; }
          for (auto& k : unshifted_y_ions[ion_index]) { c_noshifts[i] += k.intensity; }
        }

#ifdef DEBUG_OpenNuXL
        cout << "n:";
        for (auto& k : n_shifts) cout << k << " ";
        cout << endl;
        cout << "c:";
        for (auto& k : c_shifts) cout << k << " ";
        cout << endl;
        cout << "n0:";
        for (auto& k : n_noshifts) cout << k << " ";
        cout << endl;
        cout << "c0:";
        for (auto& k : c_noshifts) cout << k << " ";
        cout << endl;
#endif

        // Rules implemented:
        // 1. if cross-link on AA, then the prefix or suffix ending at this AA must be shifted
        // 2. if the previous AA in the prefix / suffix had a stronger shifted signal, then the current on is not the correct one
        // 3. if the current AA is cross-linked, then the previous AA is not cross-linked and we should observe an unshifted prefix / suffix ion
        // Sum up all intensities of shifted prefix / suffix ions
        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          sites_sum_score[i] = 0.0;
          if (n_shifts[i] == 0 && c_shifts[i] == 0) { continue; } // no shifts? no cross-link at this AA

          if (n_shifts[i] > 0)
          {
            // Rules apply only for a3,b3 and higher ions (because we rarely observe a1,b1 ions we can't check for Rule 3)
            if (i >= 2 && n_shifts[i - 1] > n_shifts[i]) continue; // Stronger signal from shifted AA before the current one? Then skip it.
            if (i >= 2 && n_noshifts[i - 1] == 0) continue; // continue if unshifted AA is missing before (left of) the shifted one.
            // sum up all intensities from this position and all longer prefixes that also carry the NA
            for (Size j = i; j != sites_sum_score.size(); ++j) { sites_sum_score[i] += n_shifts[j]; }
          }

          if (c_shifts[i] > 0)
          {
            // Rules apply only for y3 and higher ions (because we rarely observe y1 ions we can't check for Rule 3)
            if (i < c_shifts.size()-2 && c_shifts[i + 1] > c_shifts[i]) continue; // AA after has higher intensity and also shifted? Then skip it.
            if (i < c_noshifts.size()-2 && c_noshifts[i + 1] == 0) continue; // continue if unshifted AA is missing before (right of) the shifted one.
            // sum up all intensities from this position and all longer suffixes that also carry the NA
            for (int j = i; j >= 0; --j) { sites_sum_score[i] += c_shifts[j]; }
          }
        }
#ifdef DEBUG_OpenNuXL
        cout << "site sum score (shifted a/b/y-ions):";
        for (auto& k : sites_sum_score) cout << k << " ";
        cout << endl;
#endif

        #ifdef DEBUG_OpenNuXL
          OPENMS_LOG_DEBUG << "Localisation based on immonium ions: ";
        #endif
        String aas_unmodified = aas.toUnmodifiedString();
        for (Size i = 0; i != aas_unmodified.size(); ++i)
        {
          String origin = String(aas_unmodified[i]);

          for (auto& a : shifted_immonium_ions)
          {
            // compare origin (the AA) of immonium ion to current AA
            if (a.annotation[0] == 'i' && a.annotation[1] == aas_unmodified[i])
            {

              #ifdef DEBUG_OpenNuXL
                OPENMS_LOG_DEBUG << "\n" << a.annotation << " " << "\n";
              #endif
              sites_sum_score[i] += a.intensity;
            }
          }
        }
#ifdef DEBUG_OpenNuXL
        cout << "site sum score (shifted a/b/y-ions & immonium ions):";
        for (auto& k : sites_sum_score) cout << k << " ";
        cout << endl;
#endif

        String best_localization = unmodified_sequence;
        int best_localization_position = -1; // UNKNOWN
        double best_localization_score = 0;
        String localization_scores;
        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          if (sites_sum_score[i] > best_localization_score) { best_localization_score = sites_sum_score[i]; }
        }

        for (Size i = 0; i != sites_sum_score.size(); ++i)
        {
          #ifdef DEBUG_OpenNuXL
            OPENMS_LOG_DEBUG << String::number(100.0 * sites_sum_score[i], 2);
          #endif

          if (i != 0) localization_scores += ',';
          if (sites_sum_score[i] > 0 )
          {
            localization_scores += String::number(100.0 * sites_sum_score[i], 2);
          }
          else
          {
            localization_scores += "0";
          }

          if (best_localization_score > 0.0 && sites_sum_score[i] >= best_localization_score - 1e-6)
          {
            best_localization[i] = tolower(best_localization[i]);
            best_localization_position = i; // Note: check if there are situations where multiple have the same score
          }
        }
        #ifdef DEBUG_OpenNuXL
          OPENMS_LOG_DEBUG << endl;
        #endif

        // create annotation strings for shifted fragment ions
        NuXLFragmentAnnotationHelper::addShiftedPeakFragmentAnnotation_(shifted_b_ions,
                                          shifted_y_ions,
                                          shifted_a_ions,
                                          shifted_immonium_ions,
                                          annotated_marker_ions,
                                          annotated_precursor_ions,
                                          fas);

        // store score of best localization(s)
        a.localization_scores = localization_scores;
        a.best_localization = best_localization;
        a.best_localization_score = best_localization_score;
        a.best_localization_position = best_localization_position;
        a.fragment_annotations = fas;

        #ifdef DEBUG_OpenNuXL1
          OPENMS_LOG_DEBUG << "Ion centric annotation: " << endl;
          OPENMS_LOG_DEBUG << "unshifted b ions: " << endl;
          OPENMS_LOG_DEBUG << NuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("b", unshifted_b_ions) << endl;
          OPENMS_LOG_DEBUG << "unshifted y ions: " << endl;
          OPENMS_LOG_DEBUG << NuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("y", unshifted_y_ions) << endl;
          OPENMS_LOG_DEBUG << "unshifted a ions: " << endl;
          OPENMS_LOG_DEBUG << NuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("a", unshifted_a_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted b ions: " << endl;
          OPENMS_LOG_DEBUG << NuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("b", shifted_b_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted y ions: " << endl;
          OPENMS_LOG_DEBUG << NuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("y", shifted_y_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted a ions: " << endl;
          OPENMS_LOG_DEBUG << NuXLFragmentAnnotationHelper::fragmentAnnotationDetailsToString("a", shifted_a_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted immonium ions: " << endl;
          OPENMS_LOG_DEBUG << NuXLFragmentAnnotationHelper::shiftedIonsToString(shifted_immonium_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted marker ions: " << endl;
          OPENMS_LOG_DEBUG << NuXLFragmentAnnotationHelper::shiftedIonsToString(annotated_marker_ions) << endl;
          OPENMS_LOG_DEBUG << "shifted precursor ions: " << endl;
          OPENMS_LOG_DEBUG << NuXLFragmentAnnotationHelper::shiftedIonsToString(annotated_precursor_ions) << endl;
          OPENMS_LOG_DEBUG << "Localization scores: ";
          OPENMS_LOG_DEBUG << localization_scores << endl;
          OPENMS_LOG_DEBUG << "Localisation based on ion series and immonium ions of all observed fragments: ";
          OPENMS_LOG_DEBUG << best_localization << endl;
        #endif
      }
    }
  }
} // namespace

