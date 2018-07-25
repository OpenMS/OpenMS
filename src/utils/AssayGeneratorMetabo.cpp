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
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

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

  Internal procedure AssayGeneratorMetabo:
  1. Input mzML and featureXML
  2. Annotate precursor mz and intensity
  3. Filter feature by number of masstraces
  4. Assign precursors to specific feature
  5. Extract feature meta information (if possible)
  6. Find MS2 spectrum with highest intensity precursor for one feature
  7. Dependent on the method use the MS2 with the highest intensity precursor or a consensus spectrum
     for the transition calculation
  8. Calculate thresholds (maximum and minimum intensity for transition peak)
  9. Extract and write transitions

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SiriusAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_SiriusAdapter.html
 */

/// @cond TOPPCLASSES


/// struct to hold assay information of one row 
struct AssayRow
{
  double precursor_mz;
  double product_mz;
  float library_int;
  double normalized_rt;
  String compound_name;
  String smiles;
  String sumformula;
  String adduct;
  String transition_group_id;
  String transition_id;
  bool decoy; 
};

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

    //registerInputFileList_("in", "<file>", "", "MzML Input files used for assay library building");
    registerInputFile_("in", "<file>", "", "MzML Input file used for assay library building");
    setValidFormats_("in", ListUtils::create<String>("mzml"));

    registerInputFile_("in_id", "<file>", "", "FeatureXML Input with id information (accurate mass search)");
    setValidFormats_("in_id", ListUtils::create<String>("featurexml"));

    registerOutputFile_("out", "<file>", "", "Assay library output file");
    setValidFormats_("out", ListUtils::create<String>("tsv"));

    registerStringOption_("method", "<choice>", "highest_intensity", "",false);
    setValidStrings_("method", ListUtils::create<String>("highest_intensity,consensus_spectrum"));

    registerDoubleOption_("signal_to_noise", "<s/n ratio>", 0, "Write peaks with S/N > signal_to_noise values only", false);

    registerDoubleOption_("precursor_mz_tolerance", "<num>", 0.005, "Tolerance window for precursor selection (Feature selection in regard to the precursor)", false);
    registerStringOption_("precursor_mz_tolerance_unit", "<choice>", "Da", "Unit of the precursor_mz_tolerance", false);
    setValidStrings_("precursor_mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));

    registerDoubleOption_("precursor_mz_distance", "<num>", 0.0001, "Max m/z distance of the precursor entries of two spectra to be merged in [Da].", false);

    registerDoubleOption_("precursor_recalibration_window", "<num>", 0.1, "Tolerance window for precursor selection (Annotation of precursor mz and intensity)", false, true);
    registerStringOption_("precursor_recalibration_window_unit", "<choice>", "Da", "Unit of the precursor_mz_tolerance_annotation", false, true);
    setValidStrings_("precursor_recalibration_window_unit", ListUtils::create<String>("Da,ppm"));

    registerDoubleOption_("precursor_rt_tolerance", "<num>", 5, "Tolerance window (left and right) for precursor selection [seconds]", false);

    registerDoubleOption_("cosine_similarity_threshold", "<num>", 0.98, "Threshold for cosine similarity of MS2 spectras of same precursor used for consensus spectrum", false);

    registerIntOption_("filter_by_num_masstraces", "<num>", 2, "Features have to have at least x MassTraces", false);
    setMinInt_("filter_by_num_masstraces", 1);

    registerDoubleOption_("transition_threshold", "<num>", 10, "Further transitions need at least x% of the maximum intensity (default 10%)", false);

    registerFlag_("use_known_unknowns", "Use features without identification information", false);
  }

  // map with closest feature to index of ms2 spectra
  map<const BaseFeature*, std::vector<size_t>> mappingFeatureToMS2Index(const PeakMap & spectra, const KDTreeFeatureMaps& fp_map_kd, const double& precursor_mz_tolerance, const double& precursor_rt_tolerance, bool ppm)
  {
    map<const BaseFeature*, vector<size_t>> feature_ms2_spectra_map;

    // map precursors to closest feature and retrieve annotated metadata (if possible)
    for (size_t index = 0; index != spectra.size(); ++index)
    {
      if (spectra[index].getMSLevel() != 2) { continue; }

      // get precursor meta data (m/z, rt)
      const vector<Precursor> & pcs = spectra[index].getPrecursors();

      if (!pcs.empty())
      {
        const double mz = pcs[0].getMZ();
        const double rt = spectra[index].getRT();

        // query features in tolerance window
        vector<Size> matches;

        // get mz tolerance window
        std::pair<double,double> mz_tolerance_window = Math::getTolWindow(mz, precursor_mz_tolerance, ppm);
        fp_map_kd.queryRegion(rt - precursor_rt_tolerance, rt + precursor_rt_tolerance, mz_tolerance_window.first, mz_tolerance_window.second, matches, true);

        // no precursor matches the feature information found
        if (matches.empty()) { continue; }

        // in the case of multiple features in tolerance window, select the one closest in m/z to the precursor
        Size min_distance_feature_index(0);
        double min_distance(1e11);
        for (auto const & k_idx : matches)
        {
          const double f_mz = fp_map_kd.mz(k_idx);
          const double distance = fabs(f_mz - mz);
          if (distance < min_distance)
          {
            min_distance = distance;
            min_distance_feature_index = k_idx;
          }
        }
        const BaseFeature* min_distance_feature = fp_map_kd.feature(min_distance_feature_index);

        feature_ms2_spectra_map[min_distance_feature].push_back(index);
      }
    }
    return feature_ms2_spectra_map;
  }

  // precursor correction (highest intensity)
  Int getHighestIntensityPeakInMZRange(double test_mz, const MSSpectrum& spectrum1, double tolerance, bool ppm)
  {

    // get tolerance window and left/right iterator
    std::pair<double,double> tolerance_window = Math::getTolWindow(test_mz, tolerance, ppm);

    // Here left has to be smaller than right
    OPENMS_PRECONDITION(tolerance_window.first < tolerance_window.second, "Left has to be smaller than right");

    MSSpectrum::ConstIterator left = spectrum1.MZBegin(tolerance_window.first);
    MSSpectrum::ConstIterator right = spectrum1.MZBegin(tolerance_window.second);

    // no MS1 precursor peak in +- tolerance window found
    if (left == right)
    {
        return -1;
    }

    MSSpectrum::ConstIterator max_intensity_it = std::max_element(left, right, Peak1D::IntensityLess());

    return max_intensity_it - spectrum1.begin();
  }

  // annotate precursor intensity based on precursor spectrum and highest intensity peak in tolerance window
  void annotatePrecursorIntensity(PeakMap& spectra, double tolerance, bool ppm)
  {
    for (PeakMap::Iterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      // process only MS2 spectra
      if (s_it->getMSLevel() != 2)
      {
        continue;
      }

      MSSpectrum& spectrum = *s_it;
      vector<Precursor>& precursor = spectrum.getPrecursors();

      if (precursor.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: Invalid MS2 spectrum without precursor");
      }

      double test_mz = precursor[0].getMZ();

      // find corresponding precursor spectrum
      PeakMap::ConstIterator s_it2 = spectra.getPrecursorSpectrum(s_it);

      // no precursor spectrum found
      if (s_it2 == spectra.end())
      {
        LOG_WARN << "No MS1 spectrum was found to the specific precursor" << std::endl;
        continue;
      }

      const MSSpectrum& precursor_spectrum = *s_it2;
      Int mono_index = getHighestIntensityPeakInMZRange(test_mz, precursor_spectrum, tolerance, ppm);
      if (mono_index == -1)
      {
        LOG_WARN << "No precursor peak in MS1 spectrum found, please ensure that the tolerance is set correctly." << std::endl;
        continue;
      }
      const Peak1D& max_mono_peak = precursor_spectrum[mono_index];
      precursor[0].setMZ(max_mono_peak.getMZ());
      precursor[0].setIntensity(max_mono_peak.getIntensity());
    }
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String id = getStringOption_("in_id");
    String out = getStringOption_("out");
    String method = getStringOption_("method");
    bool method_consensus_spectrum = method == "consensus_spectrum" ? true : false;


    double sn = getDoubleOption_("signal_to_noise");

    double precursor_mz_tol = getDoubleOption_("precursor_mz_tolerance");
    String unit_prec = getStringOption_("precursor_mz_tolerance_unit");
    bool ppm_prec = unit_prec == "ppm" ? true : false;

    double precursor_rt_tol = getDoubleOption_("precursor_rt_tolerance");
    double pre_recal_win = getDoubleOption_("precursor_recalibration_window");
    String pre_recal_win_unit = getStringOption_("precursor_recalibration_window_unit");
    bool ppm_recal = pre_recal_win_unit == "ppm" ? true : false;

    double precursor_mz_distance = getDoubleOption_("precursor_mz_distance");

    double cosine_sim_threshold = getDoubleOption_("cosine_similarity_threshold");
    unsigned int num_masstrace_filter = getIntOption_("filter_by_num_masstraces");
    double transition_threshold = getDoubleOption_("transition_threshold");
    bool use_known_unknowns = getFlag_("use_known_unknowns");

    // load mzML
    MzMLFile mzml;
    PeakMap spectra;
    mzml.load(in, spectra);

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = spectra[0].getType();

    if (spectrum_type == SpectrumSettings::PROFILE)
    {
      if (!getFlag_("force"))
      {
        throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__,
                                           "Error: Profile data provided but centroided spectra expected.");
      }
    }

    // TODO: check if signal to noise filter works
    // TODO: signal to noise not working correctly (further parameter needed - win_len, min_required_elements)
    // calculate S/N values and delete data points below S/N threshold
    if (sn > 0)
    {
      SignalToNoiseEstimatorMedian<PeakMap::SpectrumType> snm;
      Param const &dc_param = getParam_().copy("algorithm:SignalToNoise:", true);
      snm.setParameters(dc_param);
      for (PeakMap::Iterator it = spectra.begin(); it != spectra.end(); ++it)
      {
        snm.init(it->begin(), it->end());
        for (PeakMap::SpectrumType::Iterator spec = it->begin(); spec != it->end(); ++spec)
        {
          if (snm.getSignalToNoise(spec) < sn)
            spec->setIntensity(0);
        }
        it->erase(remove_if(it->begin(),
                            it->end(),
                            InIntensityRange<PeakMap::PeakType>(1,
                                                                numeric_limits<PeakMap::PeakType::IntensityType>::max(),
                                                                true)), it->end());
      }
    }

    // load featurexml
    FeatureXMLFile fxml;
    FeatureMap feature_map;
    fxml.load(id, feature_map);

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
        }
      }
    }

    // annotate precursor mz and intensity
    annotatePrecursorIntensity(spectra, pre_recal_win, ppm_recal);

    // filter feature by number of masstraces
    auto map_it = remove_if(feature_map.begin(), feature_map.end(),
                            [&num_masstrace_filter](const Feature &f) -> bool
                            {
                              unsigned int n_masstraces = f.getMetaValue("num_of_masstraces");
                              return n_masstraces < num_masstrace_filter;
                            });
    feature_map.erase(map_it, feature_map.end());

    KDTreeFeatureMaps fp_map_kd;
    vector<FeatureMap> v_fp;
    v_fp.push_back(feature_map);
    fp_map_kd.addMaps(v_fp);

    // read FeatureMap in KDTree for feature-precursor assignment
    // only spectra with precursors are in the map - no need to check for presence of precursors
    map<const BaseFeature *, std::vector<size_t> > feature_ms2_spectra_map = mappingFeatureToMS2Index(spectra,
                                                                                                      fp_map_kd,
                                                                                                      precursor_mz_tol,
                                                                                                      precursor_rt_tol,
                                                                                                      ppm_prec);

    std::vector<AssayRow> assaylib;
    TargetedExperiment t_exp;
    vector<TargetedExperiment::Compound> v_cmp;
    vector<ReactionMonitoringTransition> v_rmt;
    int transition_group_counter = 0;

    for (std::map<const BaseFeature *, std::vector<size_t>>::iterator it = feature_ms2_spectra_map.begin();
         it != feature_ms2_spectra_map.end();
         ++it)
    {

      String description("UNKNOWN"), sumformula("UNKNOWN"), adduct("UNKNOWN");
      StringList v_description, v_sumformula, v_adduct;

      double feature_rt;
      const BaseFeature *min_distance_feature = it->first;
      feature_rt = min_distance_feature->getRT();

      // extract metadata from featureXML
      if (!(min_distance_feature->getPeptideIdentifications().empty()) &&
          !(min_distance_feature->getPeptideIdentifications()[0].getHits().empty()))
      {
        // Accurate Mass Search may provide multiple possible Hits
        for (unsigned int i = 0; i != min_distance_feature->getPeptideIdentifications()[0].getHits().size(); ++i)
        {
          description = min_distance_feature->getPeptideIdentifications()[0].getHits()[i].getMetaValue("description");
          sumformula = min_distance_feature->getPeptideIdentifications()[0].getHits()[i].getMetaValue("chemical_formula");
          adduct = min_distance_feature->getPeptideIdentifications()[0].getHits()[i].getMetaValue("modifications");

          // change format of description [name] to name
          description.erase(std::remove_if(std::begin(description),
                                           std::end(description),
                                           [](char c) { return c == '[' || c == ']'; }), std::end(description));

          // change format of adduct information M+H;1+ -> [M+H]1+
          String adduct_prefix = adduct.prefix(';').trim();
          String adduct_suffix = adduct.suffix(';').trim();
          adduct = "[" + adduct_prefix + "]" + adduct_suffix;

          v_description.push_back(description);
          v_sumformula.push_back(sumformula);
          v_adduct.push_back(adduct);
        }

        std::cout << "-----------" << std::endl;
        for (auto d_it = v_description.begin(); d_it != v_description.end(); ++d_it)
        {
          std::cout << "des: " << *d_it << std::endl;
        }

      }

      // TODO: does it still apply with the vector?
      // check if known unknown should be used
      if (description == "UNKNOWN" && sumformula == "UNKNOWN" && adduct == "UNKNOWN" && !use_known_unknowns)
      {
        continue;
      }

      double highest_precursor_mz = 0.0;
      float highest_precursor_int = 0.0;
      MSSpectrum highest_precursor_int_spectrum;
      MSSpectrum transition_spectrum;
      String native_id;

      // find precursor/spectrum with highest intensity precursor
      std::vector<size_t> index = it->second;
      for (std::vector<size_t>::iterator index_it = index.begin(); index_it != index.end(); ++index_it)
      {
        const MSSpectrum &spectrum = spectra[*index_it];
        const vector<Precursor> &precursor = spectrum.getPrecursors();

        // get m/z and intensity of precursor
        // only spectra with precursors are in the map, therefore no need to check for their presence
        double precursor_mz = precursor[0].getMZ();
        float precursor_int = precursor[0].getIntensity();

        native_id = spectrum.getNativeID();

        // spectrum with highest intensity precursor
        if (precursor_int > highest_precursor_int)
        {
          highest_precursor_int = precursor_int;
          highest_precursor_mz = precursor_mz;
          highest_precursor_int_spectrum = spectrum;
        }
        transition_spectrum = highest_precursor_int_spectrum;
      }

      // if only one MS2 is available and the consensus method is used - jump right to the transition list calculation
      // fallback: highest intensity precursor
      if (method_consensus_spectrum && index.size() >= 2)
      {
        // transform to binned spectra
        std::vector<BinnedSpectrum> binned;
        std::vector<MSSpectrum> similar_spectra;
        MSExperiment exp;
        const BinnedSpectrum binned_highest_int(highest_precursor_int_spectrum,
                                                BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES,
                                                false,
                                                1,
                                                BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

        // calculation of contrast angle (cosine simiarity)
        for (std::vector<size_t>::iterator index_it = index.begin(); index_it != index.end(); ++index_it)
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

      // transition calculations
      // calculate max intensity peak and threshold
      float max_int = 0.0;
      float min_int = numeric_limits<float>::max();
      for (MSSpectrum::const_iterator spec_it = transition_spectrum.begin();
           spec_it != transition_spectrum.end();
           ++spec_it)
      {
        //find the max intensity peak
        if (spec_it->getIntensity() > max_int)
        {
          max_int = spec_it->getIntensity();
        }
        if (spec_it->getIntensity() < min_int)
        {
          min_int = spec_it->getIntensity();
        }
      }

      // no peaks or all peaks have same intensity (single peak / noise)
      if (min_int >= max_int)
      {
        continue;
      }

      // threshold should be at x % of the maximum intensity
      // hard minimal threshold of min_int * 1.1
      float threshold_transition = max_int * (transition_threshold / 100);
      float threshold_noise = min_int * 1.1;

      AssayRow row;
      int transition_counter = 0;

      // TODO: Other datastructure (min/max transitions) - min -> transitions with highest intensity?
      // TODO: Put stuff in vector first? then sort by intensity? maybe map Transition intensity: everything else

      // TODO: test if output is correct and add adduct information.
      // if multiple files are used other datastructure is needed
      // feature mz and rt and save transition with compound use TargetedExperiment
      // and precursor intensity information
      // compare mz and rt (what about rt shift?)
      // and compare precursor intenstiy -> use the higher one for calculation of the
      // transition or use both ms2 consensus? (not sure if that is a good idea?


      ReactionMonitoringTransition rmt;
      rmt.clearMetaInfo();

      // here ms2 spectra information is used
      for (auto spec_it = transition_spectrum.begin(); spec_it != transition_spectrum.end(); ++spec_it)
      {
        float current_int = spec_it->getIntensity();
        double current_mz = spec_it->getMZ();

        // write row for each transition
        // current int has to be higher than transition thresold and should not be smaller than threshold noise
        if (current_int > threshold_transition && current_int > threshold_noise)
        {
          float rel_int = current_int / max_int;

          TargetedExperiment::Compound cmp;

          // TODO: use TE::Helper for the retention time of the feature;
          //TargetedExperiment::Helper

          //cmp.rts // nur in OpenSWATH genutzt
          rmt.setPrecursorMZ(highest_precursor_mz);
          rmt.setProductMZ(current_mz);
          rmt.setLibraryIntensity(rel_int);

          TargetedExperimentHelper::RetentionTime cmp_rt;
          cmp_rt.setRT(feature_rt);

          if (description == "UNKNOWN")
          {
            rmt.setCompoundRef(String(transition_group_counter) + "_" + String(transition_counter) + "_" + description);
            cmp.id = String(transition_group_counter) + "_" + String(transition_counter) + "_" + description;
            rmt.setNativeID(String(transition_group_counter) + "_" + description);
          }
          else
          {
            description = ListUtils::concatenate(v_description, ",");
            rmt.setCompoundRef(String(transition_group_counter) + "_" + String(transition_counter) + "_" + description);
            cmp.id = String(transition_group_counter) + "_" + String(transition_counter) + "_" + description;
            rmt.setNativeID(String(transition_group_counter) + "_" + String(transition_counter) + "_" + description);
          }
          // smiles to supported in AccurateMassSearch
          row.smiles = "none";
          if (sumformula == "UNKNOWN")
          {
            cmp.molecular_formula = sumformula;
          }
          else
          {
            sumformula = ListUtils::concatenate(v_sumformula, ",");
            cmp.molecular_formula = sumformula;
          }
          if (adduct == "UNKNOWN")
          {
            cmp.setMetaValue("adducts", adduct);
          }
          else
          {
            adduct = ListUtils::concatenate(v_adduct, ",");
            cmp.setMetaValue("adducts", adduct);
          }
          v_cmp.push_back(cmp);
          v_rmt.push_back(rmt);
        }
        transition_counter += 1;
      }
      transition_group_counter += 1;
    }
    // TODO: add adduct information in TargetedExperiment::Compound
    t_exp.setCompounds(v_cmp);
    t_exp.setTransitions(v_rmt);

    OpenMS::TransitionTSVFile::convertTargetedExperimentToTSV(out.c_str(), t_exp);
    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  TOPPAssayGeneratorMetabo tool;
  return tool.main(argc, argv);
}
/// @endcond
