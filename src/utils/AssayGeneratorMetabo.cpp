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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMAssay.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/FILTERING/CALIBRATION/PrecursorCorrection.h>

#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h> 
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h> 
#include <OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>

#include <OpenMS/KERNEL/RangeUtils.h>

#include <algorithm>
#include <map> //insert

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
  @page UTILS_AssayGeneratorMetabo AssayGeneratorMetabo

  @brief Generates an assay library using DDA data (Metabolomics)

    <CENTER>
      <table>
          <tr>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
              <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ AssayGeneratorMetabo \f$ \longrightarrow \f$</td>
              <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderMetabo </td>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref OpenSWATH pipeline </td>
          </tr>
          <tr>
              <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref Utils_AccurateMassSearch </td>
          </tr>
      </table>
  </CENTER>

  Build an assay library from DDA data (MS and MS/MS) (mzML).
  Please provide a list of features found in the data (featureXML).

  Features can be detected using the FeatureFinderMetabo (FFM) and identifcation information
  can be added using the AccurateMassSearch feautreXML output.

  If the FFM featureXML is used the "use_known_unknowns" flag is used automatically.

  Internal procedure AssayGeneratorMetabo: \n
  1. Input mzML and featureXML \n
  2. Reannotate precursor mz and intensity \n
  3. Filter feature by number of masstraces \n
  4. Assign precursors to specific feature (FeatureMapping) \n
  5. Extract feature meta information (if possible) \n
  6. Find MS2 spectrum with highest intensity precursor for one feature \n
  7. Dependent on the method use the MS2 with the highest intensity precursor or a consensus spectrum
     for the transition calculation \n
  8. Calculate thresholds (maximum and minimum intensity for transition peak) \n
  9. Extract and write transitions (tsv, traml) \n

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SiriusAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_SiriusAdapter.html
 */

/// @cond TOPPCLASSES

class TOPPAssayGeneratorMetabo :
  public TOPPBase,
  private TransitionTSVFile
{
public:
  TOPPAssayGeneratorMetabo() :
    TOPPBase("AssayGeneratorMetabo", "Assay library generation from DDA data (Metabolomics)", false)
    {}

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("executable", "<executable>", "", "sirius executable e.g. sirius", false, false, ListUtils::create<String>("skipexists"));

    registerInputFileList_("in", "<file(s)>", StringList(), "MzML input file(s) used for assay library generation");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFileList_("in_id", "<file(s)>", StringList(), "FeatureXML input file(s) containing id information (e.g. accurate mass search)");
    setValidFormats_("in_id", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out", "<file>", "", "Assay library output file");
    setValidFormats_("out", ListUtils::create<String>("tsv,traML,pqp"));

    // TODO: add method fragment annotation here? swith the other to heuristcs? s
    registerStringOption_("method", "<choice>", "highest_intensity", "Method used for assay library construction",false);
    setValidStrings_("method", ListUtils::create<String>("highest_intensity,consensus_spectrum"));

    registerFlag_("use_fragment_annotation", "Use sirius fragment annotation", false);
    registerFlag_("use_exact_mass", "Use exact mass for Fragment Annotation", false);
    registerFlag_("exclude_ms2_precursor", "Excludes precursor in ms2 from transition list", false);
    
    // preprocessing
    // TODO: add ppm for precursor mz distance? where is it used 
    registerDoubleOption_("precursor_mz_distance", "<num>", 0.0001, "Max m/z distance of the precursor entries of two spectra to be merged in [Da].", false);
    registerDoubleOption_("precursor_recalibration_window", "<num>", 0.1, "Tolerance window for precursor selection (Annotation of precursor mz and intensity)", false, true);
    registerStringOption_("precursor_recalibration_window_unit", "<choice>", "Da", "Unit of the precursor_mz_tolerance_annotation", false, true);
    setValidStrings_("precursor_recalibration_window_unit", ListUtils::create<String>("Da,ppm"));
    registerDoubleOption_("precursor_rt_tolerance", "<num>", 5, "Tolerance window (left and right) for precursor selection [seconds]", false);
    registerFlag_("use_known_unknowns", "Use features without identification information", false);

    // transition extraction 
    registerIntOption_("min_transitions", "<int>", 3, "minimal number of transitions", false);
    registerIntOption_("max_transitions", "<int>", 6, "maximal number of transitions", false);
    registerDoubleOption_("cosine_similarity_threshold", "<num>", 0.98, "Threshold for cosine similarity of MS2 spectras of same precursor used for consensus spectrum", false);
    registerDoubleOption_("transition_threshold", "<num>", 10, "Further transitions need at least x% of the maximum intensity (default 10%)", false);

    // TODO: what are the parameters for? add all comments - see desisotoper
    registerTOPPSubsection_("deisotoping", "deisotoping");
    registerFlag_("deisotoping:use_deisotoper", "Use deisotper and its parameters", false);
    registerDoubleOption_("deisotoping:fragment_tolerance", "<num>", 1, "Tolerance used to match isotopic peaks (Deisotoping)", false);
    registerStringOption_("deisotoping:fragment_unit", "<choice>", "ppm", "Unit of the fragment tolerance", false);
    setValidStrings_("deisotoping:fragment_unit", ListUtils::create<String>("ppm,Da"));
    registerIntOption_("deisotoping:min_charge", "<num>", 1, "set min charge ", false);
    setMinInt_("deisotoping:min_charge", 1);
    registerIntOption_("deisotoping:max_charge", "<num>", 1, "set max charge", false);
    setMinInt_("deisotoping:max_charge", 1);
    registerIntOption_("deisotoping:min_isopeaks", "<num>", 2, "set min charge ", false);
    setMinInt_("deisotoping:min_isopeaks", 2);
    registerIntOption_("deisotoping:max_isopeaks", "<num>", 3, "set max charge", false);
    setMinInt_("deisotoping:max_isopeaks", 3);
    registerFlag_("deisotoping:keep_only_deisotoped", "Use deisotper and its parameters", false);
    registerFlag_("deisotoping:annotate_charge", "Use deisotper and its parameters", false);

    // sirius 
    registerFullParam_(SiriusAdapterAlgorithm().getDefaults());
  }

  // datastructure used for preprocessing
  // potential transitions for one specific compound
  struct PotentialTransitions
  {
      double precursor_mz;
      double precursor_rt;
      double precursor_int;
      int precursor_charge;
      double transition_quality_score; // here precursor intensity will be used at first, till a better scoring is available
      String compound_name;
      String compound_adduct;
      TargetedExperiment::Compound potential_cmp;
      vector<ReactionMonitoringTransition> potential_rmts;
  };

  struct compare
  {
  	string key;
  	compare(string const &i): key(i){}
  
  	bool operator()(string const &i)
  	{
  		return (i == key);
  	}
  };

  static bool compoundNameLess_(const TOPPAssayGeneratorMetabo::PotentialTransitions& i, const TOPPAssayGeneratorMetabo::PotentialTransitions& j)
  {
    if (i.compound_name < j.compound_name)
    {
      return (i.compound_name < j.compound_name);
    }
    if (i.compound_name == j.compound_name)
    {
      if (i.compound_adduct == j.compound_adduct)
      {
        return (i.precursor_int > j.precursor_int);
      }
      else
      {
        return (i.compound_adduct < j.compound_adduct);
      }
    }
    else
    {
      return false;
    }
  }
 
  // use std::unique with this function to retain the unique element with the highest precursor intensity
  static bool compoundUnique_(const TOPPAssayGeneratorMetabo::PotentialTransitions& i, const TOPPAssayGeneratorMetabo::PotentialTransitions& j)
  {
    if (i.compound_name == j.compound_name && i.compound_adduct == j.compound_adduct && i.precursor_int > j.precursor_int)
    {
      return true;
    }
    else 
    {
      return false;
    }
  }

  static bool intensityLess_(Peak1D a, Peak1D b)
  {
      return (a.getIntensity() < b.getIntensity());
  }

  static bool extractAndCompareScanIndexLess_(const String& i, const String& j)
  {
    return (SiriusMzTabWriter::extract_scan_index(i) < SiriusMzTabWriter::extract_scan_index(j));
  }

  // method to extract a potential transistions based on the ms/ms based of the highest intensity precursor or a consensus spectrum
  vector<PotentialTransitions> extractPotentialTransitions(const PeakMap& spectra,
                                                           const map<BaseFeature const *, vector<size_t>>& feature_ms2_spectra_map,
                                                           const double& precursor_rt_tol,
                                                           const double& precursor_mz_distance,
                                                           const double& cosine_sim_threshold,
                                                           const double& transition_threshold,
                                                           const bool& method_consensus_spectrum,
                                                           const bool& exclude_ms2_precursor,
                                                           const unsigned int& file_counter)
  {
    int transition_group_counter = 0;
    vector<PotentialTransitions> v_pts;

    for (auto it = feature_ms2_spectra_map.begin();
              it != feature_ms2_spectra_map.end();
              ++it)
      {
        TargetedExperiment::Compound cmp;
        cmp.clearMetaInfo();
        vector<ReactionMonitoringTransition> v_rmt;

        String description("UNKNOWN"), sumformula("UNKNOWN"), adduct("UNKNOWN");
        StringList v_description, v_sumformula, v_adduct;

        double feature_rt;
        const BaseFeature *min_distance_feature = it->first;
        feature_rt = min_distance_feature->getRT();

        // extract metadata from featureXML
        if (!(min_distance_feature->getPeptideIdentifications().empty()) &&
            !(min_distance_feature->getPeptideIdentifications()[0].getHits().empty()))
        {
          // accurate mass search may provide multiple possible Hits
          // for heuristics use the identification with the smallest mz error (ppm)
          double min_id_mz_error = std::numeric_limits<double>::max();
          for (unsigned int j = 0; j != min_distance_feature->getPeptideIdentifications()[0].getHits().size(); ++j)
          {
            double current_id_mz_error = min_distance_feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("mz_error_ppm");
            // compare the absolute error absolute error
            if (abs(current_id_mz_error) < min_id_mz_error)
            {
              description = min_distance_feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("description");
              sumformula = min_distance_feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("chemical_formula");
              adduct = min_distance_feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("modifications");

              // change format of description [name] to name
              description.erase(remove_if(begin(description),
                                          end(description),
                                          [](char c) { return c == '[' || c == ']'; }), end(description));

              // change format of adduct information M+H;1+ -> [M+H]1+
              String adduct_prefix = adduct.prefix(';').trim();
              String adduct_suffix = adduct.suffix(';').trim();
              adduct = "[" + adduct_prefix + "]" + adduct_suffix;

              min_id_mz_error = abs(current_id_mz_error);
            }
          }
          // use the identification with the lowest mass deviation
          v_description.push_back(description);
          v_sumformula.push_back(sumformula);
          v_adduct.push_back(adduct);
        }

        double highest_precursor_mz = 0.0;
        float highest_precursor_int = 0.0;
        int highest_precursor_charge = 0;
        MSSpectrum highest_precursor_int_spectrum;
        MSSpectrum transition_spectrum;
        String native_id;

        // find precursor/spectrum with highest intensity precursor
        vector<size_t> index = it->second;

        for (auto index_it = index.begin();
                  index_it != index.end();
                  ++index_it)
        {
          const MSSpectrum &spectrum = spectra[*index_it];

          // check if MS2 spectrum is empty
          if (spectrum.empty())
          {
            LOG_WARN << "Empty MSMS spectrum was provided. Please hava a look. Index: " << *index_it << std::endl;
            continue;
          }

          const vector<Precursor> &precursor = spectrum.getPrecursors();

          // get m/z and intensity of precursor
          // only spectra with precursors are in the map, therefore no need to check for their presence
          double precursor_mz = precursor[0].getMZ();
          float precursor_int = precursor[0].getIntensity();
          int precursor_charge = precursor[0].getCharge();

          native_id = spectrum.getNativeID();

          // spectrum with highest intensity precursor
          if (precursor_int > highest_precursor_int)
          {
            highest_precursor_int = precursor_int;
            highest_precursor_mz = precursor_mz;
            highest_precursor_int_spectrum = spectrum;
            highest_precursor_charge = precursor_charge;
          }
          transition_spectrum = highest_precursor_int_spectrum;
        }

        // if only one MS2 is available and the consensus method is used - jump right to the transition list calculation
        // fallback: highest intensity precursor
        if (method_consensus_spectrum && index.size() >= 2)
        {
          // transform to binned spectra
          vector<BinnedSpectrum> binned;
          vector<MSSpectrum> similar_spectra;
          MSExperiment exp;
          const BinnedSpectrum binned_highest_int(highest_precursor_int_spectrum,
                                                  BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES,
                                                  false,
                                                  1,
                                                  BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

          // calculation of contrast angle (cosine simiarity)
          for (auto index_it = index.begin();
                    index_it != index.end();
                    ++index_it)
          {
            const MSSpectrum &spectrum = spectra[*index_it];
            const BinnedSpectrum binned_spectrum(spectrum,
                                                 BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES,
                                                 false,
                                                 1,
                                                 BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

            BinnedSpectralContrastAngle bspa;
            double cosine_sim = bspa(binned_highest_int, binned_spectrum);

            if (cosine_sim > cosine_sim_threshold)
            {
              similar_spectra.push_back(spectrum);
              exp.addSpectrum(spectrum);
            }
          }

          // calculate consensus spectrum
          exp.sortSpectra();
          SpectraMerger merger;
          Param p;
          p.insert("", SpectraMerger().getDefaults());
          p.setValue("precursor_method:mz_tolerance", precursor_mz_distance);
          p.setValue("precursor_method:rt_tolerance", precursor_rt_tol * 2);
          merger.setParameters(p);

          // all MS spectra should have the same precursor
          merger.mergeSpectraPrecursors(exp);

          // check if all precursors have been merged if not use highest intensity precursor
          if (exp.getSpectra().size() < 2)
          {
            transition_spectrum = exp.getSpectra()[0];
          }
        }

        // TODO: Why would transition spectrum have a size of 0 ? (Empty MS2 spectrum)
        // check if transition spectrum is empty
        if (transition_spectrum.empty())
        {
          LOG_WARN << "Empty transition spectrum was provided." << std::endl;
          continue;
        }

        // transition calculations
        // calculate max intensity peak and threshold
        float max_int = 0.0;
        float min_int = std::numeric_limits<float>::max();

        // sort intensity in MS2 spectrum to extract transitions
        transition_spectrum.sortByIntensity(true);

        // filter out the precursors if they are in the ms2 spectrum;
        if (exclude_ms2_precursor)
        {
          for (auto spec_it = transition_spectrum.begin();
                    spec_it != transition_spectrum.end();
                    ++spec_it)
          {
            if (transition_spectrum.getPrecursors()[0].getMZ() == spec_it->getMZ())
            {
              transition_spectrum.erase(spec_it);
              break;
            }
          }
        }

        // find max and min intensity peak
        max_int = max_element(transition_spectrum.begin(),transition_spectrum.end(), intensityLess_)->getIntensity();
        min_int = min_element(transition_spectrum.begin(),transition_spectrum.end(), intensityLess_)->getIntensity();

        // no peaks or all peaks have same intensity (single peak / noise)
        if (min_int >= max_int)
        {
          continue;
        }

        vector<TargetedExperimentHelper::RetentionTime> v_cmp_rt;
        TargetedExperimentHelper::RetentionTime cmp_rt;
        cmp_rt.setRT(feature_rt);
        v_cmp_rt.push_back(cmp_rt);
        cmp.rts = v_cmp_rt;

        if (description == "UNKNOWN")
        {
          cmp.id = String(transition_group_counter) + "_" + description + "_" + file_counter;
          cmp.setMetaValue("CompoundName", description);
        }
        else
        {
          description = ListUtils::concatenate(v_description, ",");
          cmp.id = String(transition_group_counter) + "_" + description + "_" + file_counter;
          cmp.setMetaValue("CompoundName", description);
        }
        cmp.smiles_string = "NA";
        if (sumformula == "UNKNOWN")
        {
          cmp.molecular_formula = sumformula;
        }
        else
        {
          sumformula = ListUtils::concatenate(v_sumformula, ",");
          // sumformula = v_sumformula[0];
          cmp.molecular_formula = sumformula;
        }
        if (adduct == "UNKNOWN")
        {
          cmp.setMetaValue("Adducts", adduct);
        }
        else
        {
          // only one adduct for each descirption using the lowest mass error
          // adduct = v_adduct[0];
          adduct = ListUtils::concatenate(v_adduct, ",");
          cmp.setMetaValue("Adducts", adduct);
        }

        // threshold should be at x % of the maximum intensity
        // hard minimal threshold of min_int * 1.1
        float threshold_transition = max_int * (transition_threshold / 100);
        float threshold_noise = min_int * 1.1;

        int transition_counter = 0;
        // here ms2 spectra information is used
        for (auto spec_it = transition_spectrum.begin();
                  spec_it != transition_spectrum.end();
                  ++spec_it)
        {
          ReactionMonitoringTransition rmt;
          rmt.clearMetaInfo();

          float current_int = spec_it->getIntensity();
          double current_mz = spec_it->getMZ();

          // write row for each transition
          // current int has to be higher than transition thresold and should not be smaller than threshold noise
          if (current_int > threshold_transition && current_int > threshold_noise)
          {
            float rel_int = current_int / max_int;

            rmt.setPrecursorMZ(highest_precursor_mz);
            rmt.setProductMZ(current_mz);
            rmt.setLibraryIntensity(rel_int);

            if (description == "UNKNOWN")
            {
              rmt.setCompoundRef(String(transition_group_counter) + "_" + description + "_" + file_counter);
              rmt.setNativeID(String(transition_group_counter) + "_" + String(transition_counter) + "_" + description + "_" + file_counter);
            }
            else
            {
              description = ListUtils::concatenate(v_description, ",");
              rmt.setCompoundRef(String(transition_group_counter) + "_" + description + "_" + file_counter);
              rmt.setNativeID(String(transition_group_counter) + "_" + String(transition_counter) + "_" + description + "_" + file_counter);
            }
            v_rmt.push_back(rmt);
            transition_counter += 1;
          }
        }
        transition_group_counter += 1;
        PotentialTransitions pts;
        pts.precursor_mz = highest_precursor_mz;
        pts.precursor_rt = feature_rt;
        pts.precursor_int = highest_precursor_int;
        pts.precursor_charge = highest_precursor_charge;
        //pts.transition_quality =
        pts.compound_name = description;
        pts.compound_adduct = adduct;
        pts.potential_cmp = cmp;
        pts.potential_rmts = v_rmt;
        v_pts.push_back(pts);
      }
      return v_pts;
  }

  // method to extract a potential transitions based on the ms/ms based of the highest intensity precursor with fragment annotation using SIRIUS
  vector<PotentialTransitions> extractPotentialTransitionsFragmentAnnotation(const vector< pair <SiriusMSFile::CompoundInfo, MSSpectrum> >& v_cmp_spec,
                                                                             const double& transition_threshold,
                                                                             const bool& use_exact_mass,
                                                                             const bool& exclude_ms2_precursor,
                                                                             const unsigned int& file_counter)
  {
    int transition_group_counter = 0;
    vector<PotentialTransitions> v_pts;

    for (auto it = v_cmp_spec.begin();
              it != v_cmp_spec.end();
              ++it)
      {
        // check if annotated spectrum exists
        MSSpectrum transition_spectrum;
        transition_spectrum = it->second;
        if (transition_spectrum.empty()) { continue; }

        TargetedExperiment::Compound cmp;
        cmp.clearMetaInfo();
        vector<ReactionMonitoringTransition> v_rmt;

        String description("UNKNOWN"), sumformula("UNKNOWN"), adduct("UNKNOWN");

        double feature_rt;
        feature_rt = it->first.rt;
        description = it->first.des;
        int charge = it->first.charge;

        // use annotated metadata
        sumformula = it->second.getMetaValue("annotated_sumformula");
        adduct = it->second.getMetaValue("annotated_adduct");

        // transition calculations
        // calculate max intensity peak and threshold
        float max_int = 0.0;
        float min_int = std::numeric_limits<float>::max();

        // sort intensity in MS2 spectrum to extract transitions
        transition_spectrum.sortByIntensity(true);

        // have to remove ms2 precursor peak before min/max
        double exact_mass_precursor = 0.0;
        for (auto spec_it = transition_spectrum.begin();
                  spec_it != transition_spectrum.end();
                  ++spec_it)
        {
          int spec_index = spec_it - transition_spectrum.begin();
          OpenMS::DataArrays::StringDataArray explanation_array = transition_spectrum.getStringDataArrays()[0];
          if (explanation_array.getName() != "explanation")
          {
            LOG_WARN << "Fragment explanation was not found. Please check if you annotation works properly." << std::endl;
          }
          else
          {
            // precursor in fragment annotation has the same sumformula as MS1 Precursor
            if (explanation_array[spec_index] == sumformula)
            {
              // save exact mass
              if(use_exact_mass)
              {
                exact_mass_precursor = spec_it->getMZ();
              }
              // remove precursor ms2 entry
              if(exclude_ms2_precursor)
              {
                transition_spectrum.erase(transition_spectrum.begin() + spec_index);
                transition_spectrum.getStringDataArrays()[0].erase(transition_spectrum.getStringDataArrays()[0].begin() + spec_index);
                transition_spectrum.getFloatDataArrays()[0].erase(transition_spectrum.getFloatDataArrays()[0].begin() + spec_index);
                break; // if last element have to break if not iterator will go out of range
              }
            }
          }
        }

        // find max and min intensity peak
        max_int = max_element(transition_spectrum.begin(),transition_spectrum.end(), intensityLess_)->getIntensity();
        min_int = min_element(transition_spectrum.begin(),transition_spectrum.end(), intensityLess_)->getIntensity();

        // no peaks or all peaks have same intensity (single peak / noise)
        if (min_int >= max_int)
        {
          continue;
        }

        vector<TargetedExperimentHelper::RetentionTime> v_cmp_rt;
        TargetedExperimentHelper::RetentionTime cmp_rt;
        cmp_rt.setRT(feature_rt);
        v_cmp_rt.push_back(cmp_rt);
        cmp.rts = v_cmp_rt;

        cmp.id = String(transition_group_counter) + "_" + description + "_" + file_counter;
        cmp.setMetaValue("CompoundName", description);
        cmp.smiles_string = "NA";

        cmp.molecular_formula = sumformula;
        cmp.setMetaValue("Adducts", adduct);

        // threshold should be at x % of the maximum intensity
        // hard minimal threshold of min_int * 1.1
        float threshold_transition = max_int * (transition_threshold / 100);
        float threshold_noise = min_int * 1.1;
        int transition_counter = 0;

        // extract current StringDataArry with annotations/explanations;
        OpenMS::DataArrays::StringDataArray explanation_array = transition_spectrum.getStringDataArrays()[0];

        // check which entry is saved in the FloatDataArry (control for "use_exact_mass")
        LOG_DEBUG << transition_spectrum.getFloatDataArrays()[0].getName() << " is not used to build the assay library." << std::endl;

        // here ms2 spectra information is used
        for (auto spec_it = transition_spectrum.begin();
                  spec_it != transition_spectrum.end();
                  ++spec_it)
        {
          ReactionMonitoringTransition rmt;
          rmt.clearMetaInfo();
          int peak_index = spec_it - transition_spectrum.begin();

          float current_int = spec_it->getIntensity();
          double current_mz = spec_it->getMZ();
          String current_explanation = explanation_array[peak_index];

          // write row for each transition
          // current int has to be higher than transition thresold and should not be smaller than threshold noise
          if (current_int > threshold_transition && current_int > threshold_noise)
          {
            float rel_int = current_int / max_int;

            rmt.setPrecursorMZ((use_exact_mass && exact_mass_precursor != 0.0)  ? exact_mass_precursor : it->first.pmass);
            rmt.setProductMZ(current_mz);
            rmt.setLibraryIntensity(rel_int);

            rmt.setCompoundRef(String(transition_group_counter) + "_" + description + "_" + file_counter);
            rmt.setNativeID(String(transition_group_counter) + "_" + String(transition_counter) + "_" + description + "_" + file_counter);

            rmt.setMetaValue("annotation", DataValue(current_explanation));

            v_rmt.push_back(rmt);
            transition_counter += 1;
          }
        }

        transition_group_counter += 1;
        PotentialTransitions pts;
        pts.precursor_mz = (use_exact_mass && exact_mass_precursor != 0.0)  ? exact_mass_precursor : it->first.pmass;
        pts.precursor_rt = it->first.rt;
        pts.precursor_int = 0;
        pts.precursor_charge = charge;
        //pts.transition_quality = ;
        pts.compound_name = description;
        pts.compound_adduct = adduct;
        pts.potential_cmp = cmp;
        pts.potential_rmts = v_rmt;
        v_pts.push_back(pts);
      }
      return v_pts;
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    // param AssayGeneratorMetabo
    StringList in = getStringList_("in");
    StringList id = getStringList_("in_id");
    String out = getStringOption_("out");
    String method = getStringOption_("method");
    bool method_consensus_spectrum = method == "consensus_spectrum" ? true : false;
    bool use_exact_mass = getFlag_("use_exact_mass");
    bool use_fragment_annotation = getFlag_("use_fragment_annotation");
    bool exclude_ms2_precursor = getFlag_("exclude_ms2_precursor");

    int min_transitions = getIntOption_("min_transitions");
    int max_transitions = getIntOption_("max_transitions");

    double precursor_rt_tol = getDoubleOption_("precursor_rt_tolerance");
    double pre_recal_win = getDoubleOption_("precursor_recalibration_window");
    String pre_recal_win_unit = getStringOption_("precursor_recalibration_window_unit");
    bool ppm_recal = pre_recal_win_unit == "ppm" ? true : false;

    double precursor_mz_distance = getDoubleOption_("precursor_mz_distance");

    double cosine_sim_threshold = getDoubleOption_("cosine_similarity_threshold");
    double transition_threshold = getDoubleOption_("transition_threshold");
    bool use_known_unknowns = getFlag_("use_known_unknowns");

    // param deisotoper
    bool use_deisotoper = getFlag_("deisotoping:use_deisotoper");
    double fragment_tolerance = getDoubleOption_("deisotoping:fragment_tolerance");
    String fragment_unit = getStringOption_("deisotoping:fragment_unit");
    bool fragment_unit_ppm = fragment_unit == "ppm" ? true : false;
    int min_charge = getIntOption_("deisotoping:min_charge");
    int max_charge = getIntOption_("deisotoping:max_charge");
    unsigned int min_isopeaks = getIntOption_("deisotoping:min_isopeaks");
    unsigned int max_isopeaks = getIntOption_("deisotoping:max_isopeaks");
    bool keep_only_deisotoped = getFlag_("deisotoping:keep_only_deisotoped");
    bool annotate_charge = getFlag_("deisotoping:annotate_charge");

    // param SiriusAdapterAlgorithm
    String executable = getStringOption_("executable");
    Param combined; 
    SiriusAdapterAlgorithm sirius_algo;
    Param preprocessing = getParam_().copy("preprocessing", false);
    Param sirius = getParam_().copy("sirius", false);
    combined.insert("", preprocessing);
    combined.insert("", sirius);
    sirius_algo.setParameters(combined);

    //-------------------------------------------------------------
    // input and check
    //-------------------------------------------------------------

    // check size of .mzML & .featureXML input
    if (in.size() != id.size())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Number of .mzML do not match to the number of .featureXML files. \n Please check and provide the corresponding files.");
    }

    vector<PotentialTransitions> v_pts;

    // iterate over all the files
    for (unsigned file_counter = 0; file_counter < in.size(); file_counter++)
    {
      // load mzML
      MzMLFile mzml;
      PeakMap spectra;
      mzml.load(in[file_counter], spectra);

      // load featurexml
      FeatureXMLFile fxml;
      FeatureMap feature_map;
      fxml.load(id[file_counter], feature_map);

      // need featureXML with Sourcefile have a look additional to ams
      StringList mzml_primary_path;
      StringList featurexml_primary_path;
      spectra.getPrimaryMSRunPath(mzml_primary_path);
      feature_map.getPrimaryMSRunPath(featurexml_primary_path);

      if (mzml_primary_path != featurexml_primary_path)
      {
        throw Exception::MissingInformation(__FILE__,
                                            __LINE__,
                                            OPENMS_PRETTY_FUNCTION,
                                            "Path of the original input file do not match in the .mzML and .featureXML files. \n Please check and provide the corresponding files.");
      }

      // determine type of spectral data (profile or centroided)
      SpectrumSettings::SpectrumType spectrum_type = spectra[0].getType();

      if (spectrum_type == SpectrumSettings::PROFILE)
      {
        if (!getFlag_("force"))
        {
          throw OpenMS::Exception::FileEmpty(__FILE__,
                                             __LINE__,
                                             __FUNCTION__,
                                             "Error: Profile data provided but centroided spectra expected. ");
        }
      }

      //-------------------------------------------------------------
      // Processing
      //-------------------------------------------------------------

      // sort spectra
      spectra.sortSpectra();

      // check if correct featureXML is given and set use_known_unkowns parameter if no id information is available
      const std::vector<DataProcessing> &processing = feature_map.getDataProcessing();
      for (auto it = processing.begin(); it != processing.end(); ++it)
      {
        if (it->getSoftware().getName() == "FeatureFinderMetabo")
        {
          // if id information is missing set use_known_unknowns to true
          if (feature_map.getProteinIdentifications().empty())
          {
            use_known_unknowns = true;
            LOG_INFO << "Due to the use of data without previous identification "
                     << "use_known_unknowns will be switched on." << std::endl;
          }
        }
      }

      // annotate and recalibrate precursor mz and intensity
      vector<double> delta_mzs;
      vector<double> mzs;
      vector<double> rts;
      PrecursorCorrection::correctToHighestIntensityMS1Peak(spectra, pre_recal_win, ppm_recal, delta_mzs, mzs, rts);

      // always use preprocessing: 
      // run masstrace filter and feature mapping
      vector<FeatureMap> v_fp; // copy FeatureMap via push_back
      KDTreeFeatureMaps fp_map_kd; // reference to *basefeature in vector<FeatureMap>
      FeatureMapping::FeatureToMs2Indices feature_mapping; // reference to *basefeature in vector<FeatureMap>
      SiriusAdapterAlgorithm::preprocessingSirius(id[file_counter],
                                                  spectra,
                                                  v_fp,
                                                  fp_map_kd,
                                                  sirius_algo,
                                                  feature_mapping);
    
      // filter known_unkowns based on description (UNKNOWN)
      std::map<const BaseFeature*, std::vector<size_t>> feature_ms2_spectra_map = feature_mapping.assignedMS2;
      std::map<const BaseFeature*, std::vector<size_t>> known_features; 
      if (!use_known_unknowns)
      {
        for (auto it = feature_ms2_spectra_map.begin(); it != feature_ms2_spectra_map.end(); ++it)
        {
          const BaseFeature *feature = it->first;
          if (!(feature->getPeptideIdentifications().empty()) &&
              !(feature->getPeptideIdentifications()[0].getHits().empty()))
              {
                String description;
                // one hit is enough for prefiltering
                description = feature->getPeptideIdentifications()[0].getHits()[0].getMetaValue("description");
                // change format of description [name] to name
                description.erase(remove_if(begin(description),
                                            end(description),
                                            [](char c) { return c == '[' || c == ']'; }), end(description));
                known_features.insert({it->first, it->second});
              }
        }
        feature_mapping.assignedMS2 = known_features;
      }

      vector< pair <SiriusMSFile::CompoundInfo, MSSpectrum> > v_cmp_spec;
      if (use_fragment_annotation)
      {
        // make temporary files
        SiriusAdapterAlgorithm::SiriusTmpStruct sirius_tmp = SiriusAdapterAlgorithm::constructSiriusTmpStruct();
        String tmp_dir = sirius_tmp.tmp_dir;
        String tmp_ms_file = sirius_tmp.tmp_ms_file;
        String tmp_out_dir = sirius_tmp.tmp_out_dir;
  
        // write msfile and store the compound information in CompoundInfo Object
        vector<SiriusMSFile::CompoundInfo> v_cmpinfo;
        bool feature_only = (sirius_algo.getFeatureOnly() == "true") ? true : false;
        bool no_mt_info = (sirius_algo.getNoMasstraceInfoIsotopePattern() == "true") ? true : false;
        int isotope_pattern_iterations = sirius_algo.getIsotopePatternIterations();
        SiriusMSFile::store(spectra,
                            tmp_ms_file,
                            feature_mapping,
                            feature_only,
                            isotope_pattern_iterations,
                            no_mt_info,
                            v_cmpinfo);

        SiriusAdapterAlgorithm::checkFeatureSpectraNumber(id[file_counter],
                                                          feature_mapping,
                                                          spectra,
                                                          sirius_algo);
  
        // calls SIRIUS and returns vector of paths to sirius folder structure
        std::vector<String> subdirs;
        String out_csifingerid;
        subdirs = SiriusAdapterAlgorithm::callSiriusQProcess(tmp_ms_file,
                                                             tmp_out_dir,
                                                             executable,
                                                             out_csifingerid,
                                                             sirius_algo);
  
        // sort vector path list
        std::sort(subdirs.begin(), subdirs.end(), extractAndCompareScanIndexLess_);
        LOG_DEBUG << subdirs.size() << " spectra were annotated using SIRIUS." << std::endl;
  
        // get Sirius FragmentAnnotion from subdirs
        vector<MSSpectrum> annotated_spectra;
        for (auto subdir : subdirs)
        {
          MSSpectrum annotated_spectrum;
          SiriusFragmentAnnotation::extractSiriusFragmentAnnotationMapping(subdir, 
                                                                           annotated_spectrum, 
                                                                           use_exact_mass);
          annotated_spectra.push_back(annotated_spectrum);
        }
        
        // clean tmp directory if debug level < 2 
        if (debug_level_ >= 2)
        {
          writeDebug_("Keeping temporary files in directory '" + tmp_dir + " and msfile at this location "+ tmp_ms_file + ". Set debug level to 1 or lower to remove them.", 2);
        }
        else
        {
          if (tmp_dir.empty() == false)
          {
            writeDebug_("Deleting temporary directory '" + tmp_dir + "'. Set debug level to 2 or higher to keep it.", 0);
            File::removeDir(tmp_dir.toQString());
          }
          if (tmp_ms_file.empty() == false)
          {
            writeDebug_("Deleting temporary msfile '" + tmp_ms_file + "'. Set debug level to 2 or higher to keep it.", 0);
            File::remove(tmp_ms_file); // remove msfile
          }
        }

        // pair compoundInfo and fragment annotation msspectrum
        // TODO: same nativeID with two descriptions, which one is the right one
        for (auto cmp : v_cmpinfo)
        {
          for (auto spec_fa : annotated_spectra)
          {
            if(std::any_of(cmp.native_ids.begin(), cmp.native_ids.end(), compare(spec_fa.getNativeID())))
            {
              v_cmp_spec.push_back(std::make_pair(cmp,spec_fa));
            }
          }
        }
      }
      else // use heuristics
      {
        if (use_deisotoper)
        {
          bool make_single_charged = false;
          for (auto& peakmap_it : spectra)
          {
            MSSpectrum& spectrum = peakmap_it;
            if (spectrum.getMSLevel() == 1) 
            {
              continue;
            }
            else 
            {
              Deisotoper::deisotopeAndSingleCharge(spectrum,
                                                  fragment_tolerance,
                                                  fragment_unit_ppm,
                                                  min_charge,
                                                  max_charge,
                                                  keep_only_deisotoped,
                                                  min_isopeaks,
                                                  max_isopeaks,
                                                  make_single_charged,
                                                  annotate_charge);
            }
          }
        }

        // remove peaks form MS2 which are at a higher mz than the precursor + 10 ppm
        for (auto& peakmap_it : spectra)
        {
          MSSpectrum& spectrum = peakmap_it;
          if (spectrum.getMSLevel() == 1) 
          {
            continue;
          }
          else 
          {
            // if peak mz higher than precursor mz set intensity to zero
            double prec_mz = spectrum.getPrecursors()[0].getMZ();
            double mass_diff = Math::ppmToMass(10.0, prec_mz);
            for (auto& spec : spectrum)
            {
              if (spec.getMZ() > prec_mz + mass_diff)
              {
                spec.setIntensity(0);
              }
            }
            spectrum.erase(remove_if(spectrum.begin(),
                                     spectrum.end(),
                                     InIntensityRange<PeakMap::PeakType>(1,
                                                                         numeric_limits<PeakMap::PeakType::IntensityType>::max(),
                                                                         true)), spectrum.end());
          }
        }
      }

      // potential transitions of one file
      vector<PotentialTransitions> tmp_pts;
      if (use_fragment_annotation)
      {
        tmp_pts = extractPotentialTransitionsFragmentAnnotation(v_cmp_spec,
                                                                transition_threshold,
                                                                use_exact_mass,
                                                                exclude_ms2_precursor,
                                                                file_counter);
      }
      else // use heuristics
      {
        tmp_pts = extractPotentialTransitions(spectra,
                                              feature_ms2_spectra_map,
                                              precursor_rt_tol,
                                              precursor_mz_distance,
                                              cosine_sim_threshold,
                                              transition_threshold,
                                              method_consensus_spectrum,
                                              exclude_ms2_precursor,
                                              file_counter);
      }
      
      // append potential transitions of one file to vector of all files
      v_pts.insert(v_pts.end(), tmp_pts.begin(),tmp_pts.end());
      
    } // end iteration over all files
    
    // sort by CompoundName
    std::sort(v_pts.begin(), v_pts.end(), compoundNameLess_);
    
    // get unique elements (CompoundName, CompoundAdduct) with highest precursor intensity
    auto uni_it = std::unique(v_pts.begin(), v_pts.end(), compoundUnique_);
    v_pts.resize(std::distance(v_pts.begin(), uni_it)); 

    // TODO: test what happens with use_known_unkowns
    // TODO: add "compoundgroup" see ProteinGroup (e.g. with different adducts - give the unique id)
    // TODO: add ID same compound with different adduct for later e.g. quantification/mapping
    
    // merge possible transitions
    vector<TargetedExperiment::Compound> v_cmp;
    vector<ReactionMonitoringTransition> v_rmt_all;
    for (auto it = v_pts.begin(); it != v_pts.end(); ++it)
    {
      v_cmp.push_back(it->potential_cmp);
      v_rmt_all.insert(v_rmt_all.end(), it->potential_rmts.begin(), it->potential_rmts.end());
    }

    // convert possible transitions to TargetedExperiment
    TargetedExperiment t_exp;
    t_exp.setCompounds(v_cmp);
    t_exp.setTransitions(v_rmt_all);

    // use MRMAssay methods for filtering
    MRMAssay assay;

    // filter: min/max transitions
    assay.detectingTransitionsCompound(t_exp, min_transitions, max_transitions);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    String extension = out.substr(out.find_last_of(".")+1);

    if (extension == "tsv")
    {
      // validate and write
      OpenMS::TransitionTSVFile::convertTargetedExperimentToTSV(out.c_str(), t_exp);
    }
    else if (extension == "traML")
    {
      // validate
      OpenMS::TransitionTSVFile::validateTargetedExperiment(t_exp);
      // write traML
      TraMLFile traml_out;
      traml_out.store(out, t_exp);
    }
    else if (extension == "pqp")
    {
      //validate 
      OpenMS::TransitionTSVFile::validateTargetedExperiment(t_exp);
      // write pqp
      TransitionPQPFile pqp_out;
      pqp_out.convertTargetedExperimentToPQP(out.c_str(), t_exp);
    }

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  TOPPAssayGeneratorMetabo tool;
  return tool.main(argc, argv);
}
/// @endcond
