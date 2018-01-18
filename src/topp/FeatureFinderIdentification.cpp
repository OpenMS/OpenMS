// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/SVM/SimpleSVM.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Container.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopeDistribution.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <cstdlib> // for "rand"
#include <ctime> // for "time" (seeding of random number generator)
#include <iterator> // for "inserter", "back_inserter"

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_FeatureFinderIdentification FeatureFinderIdentification

   @brief Detects features in MS1 data based on peptide identifications.

   <CENTER>
     <table>
       <tr>
         <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
         <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderIdentification \f$ \longrightarrow \f$</td>
         <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
       </tr>
       <tr>
         <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes (optional) </td>
         <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_ProteinQuantifier</td>
       </tr>
       <tr>
         <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
       </tr>
     </table>
   </CENTER>

   This tool detects quantitative features in MS1 data based on information from peptide identifications (derived from MS2 spectra).
   It uses algorithms for targeted data analysis from the OpenSWATH pipeline.

   The aim is to detect features that enable the quantification of (ideally) all peptides in the identification input.
   This is based on the following principle: When a high-confidence identification (ID) of a peptide was made based on an MS2 spectrum from a certain (precursor) position in the LC-MS map, this indicates that the particular peptide is present at that position, so a feature for it should be detectable there.

   @note It is important that only high-confidence (i.e. reliable) peptide identifications are used as input!

   Targeted data analysis on the MS1 level uses OpenSWATH algorithms and follows roughly the steps outlined below.

   <B>Use of inferred ("external") IDs</B>

   The situation becomes more complicated when several LC-MS/MS runs from related samples of a label-free experiment are considered.
   In order to quantify a larger fraction of the peptides/proteins in the samples, it is desirable to infer peptide identifications across runs.
   Ideally, all peptides identified in any of the runs should be quantified in each and every run.
   However, for feature detection of inferred ("external") IDs, the following problems arise:
   First, retention times may be shifted between the run being quantified and the run that gave rise to the ID.
   Such shifts can be corrected (see @ref TOPP_MapAlignerIdentification), but only to an extent.
   Thus, the RT location of the inferred ID may not necessarily lie within the RT range of the correct feature.
   Second, since the peptide in question was not directly identified in the run being quantified, it may not actually be present in detectable amounts in that sample, e.g. due to differential regulation of the corresponding protein.
   There is thus a risk of introducing false-positive features.

   FeatureFinderIdentification deals with these challenges by explicitly distinguishing between internal IDs (derived from the LC-MS/MS run being quantified) and external IDs (inferred from related runs).
   Features derived from internal IDs give rise to a training dataset for an SVM classifier.
   The SVM is then used to predict which feature candidates derived from external IDs are most likely to be correct.
   See steps 4 and 5 below for more details.

   <B>1. Assay generation</B>

   Feature detection is based on assays for identified peptides, each of which incorporates the retention time (RT), mass-to-charge ratio (m/z), and isotopic distribution (derived from the sequence) of a peptide.
   Peptides with different modifications are considered different peptides.
   One assay will be generated for every combination of (modified) peptide sequence, charge state, and RT region that has been identified.
   The RT regions arise by pooling all identifications of the same peptide, considering a window of size @p extract:rt_window around every RT location that gave rise to an ID, and then merging overlapping windows.

   <B>2. Ion chromatogram extraction</B>

   Ion chromatograms (XICs) are extracted from the LC-MS data (parameter @p in).
   One XIC per isotope in an assay is generated, with the corresponding m/z value and RT range (variable, depending on the RT region of the assay).

   @see @ref TOPP_OpenSwathChromatogramExtractor

   <B>3. Feature detection</B>

   Next feature candidates - typically several per assay - are detected in the XICs and scored.
   A variety of scores for different quality aspects are calculated by OpenSWATH.

   @see @ref TOPP_OpenSwathAnalyzer

   <B>4. Feature classification</B>

   Feature candidates derived from assays with "internal" IDs are classed as "negative" (candidates without matching internal IDs), "positive" (the single best candidate per assay with matching internal IDs), and "ambiguous" (other candidates with matching internal IDs).
   If "external" IDs were given as input, features based on them are initially classed as "unknown".
   Also in this case, a support vector machine (SVM) is trained on the "positive" and "negative" candidates, to distinguish between the two classes based on the different OpenSWATH quality scores (plus an RT deviation score).
   After parameter optimization by cross-validation, the resulting SVM is used to predict the probability of "unknown" feature candidates being positives.

   <B>5. Feature filtering</B>

   Feature candidates are filtered so that at most one feature per peptide and charge state remains.
   For assays with internal IDs, only candidates previously classed as "positive" are kept.
   For assays based solely on external IDs, the feature candidate with the highest SVM probability is selected and kept (possibly subject to the @p svm:min_prob threshold).

   <B>6. Elution model fitting</B>

   Elution models can be fitted to the features to improve the quantification.
   For robustness, one model is fitted to all isotopic mass traces of a feature in parallel.
   A symmetric (Gaussian) and an asymmetric (exponential-Gaussian hybrid) model type are available.
   The fitted models are checked for plausibility before they are accepted.

   Finally the results (feature maps, parameter @p out) are returned.

   @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

   <B>The command line parameters of this tool are:</B>
   @verbinclude TOPP_FeatureFinderIdentification.cli
   <B>INI file documentation of this tool:</B>
   @htmlinclude TOPP_FeatureFinderIdentification.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPFeatureFinderIdentification :
  public TOPPBase
{
public:

  // TODO
  // cppcheck-suppress uninitMemberVar
  TOPPFeatureFinderIdentification() :
    TOPPBase("FeatureFinderIdentification", "Detects features in MS1 data based on peptide identifications.")
  {
    // available scores: initialPeakQuality,total_xic,peak_apices_sum,var_xcorr_coelution,var_xcorr_coelution_weighted,var_xcorr_shape,var_xcorr_shape_weighted,var_library_corr,var_library_rmsd,var_library_sangle,var_library_rootmeansquare,var_library_manhattan,var_library_dotprod,var_intensity_score,nr_peaks,sn_ratio,var_log_sn_score,var_elution_model_fit_score,xx_lda_prelim_score,var_isotope_correlation_score,var_isotope_overlap_score,var_massdev_score,var_massdev_score_weighted,var_bseries_score,var_yseries_score,var_dotprod_score,var_manhatt_score,main_var_xx_swath_prelim_score,xx_swath_prelim_score
    // exclude some redundant/uninformative scores:
    // @TODO: intensity bias introduced by "peak_apices_sum"?
    score_metavalues_ = "peak_apices_sum,var_xcorr_coelution,var_xcorr_shape,var_library_sangle,var_intensity_score,sn_ratio,var_log_sn_score,var_elution_model_fit_score,xx_lda_prelim_score,var_isotope_correlation_score,var_isotope_overlap_score,var_massdev_score,main_var_xx_swath_prelim_score";
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file: LC-MS raw data");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("id", "<file>", "", "Input file: Peptide identifications derived directly from 'in'");
    setValidFormats_("id", ListUtils::create<String>("idXML"));
    registerInputFile_("id_ext", "<file>", "", "Input file: 'External' peptide identifications (e.g. from aligned runs)", false);
    setValidFormats_("id_ext", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file: Features");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerOutputFile_("lib_out", "<file>", "", "Output file: Assay library", false);
    setValidFormats_("lib_out", ListUtils::create<String>("traML"));
    registerOutputFile_("chrom_out", "<file>", "", "Output file: Chromatograms", false);
    setValidFormats_("chrom_out", ListUtils::create<String>("mzML"));
    registerOutputFile_("candidates_out", "<file>", "", "Output file: Feature candidates (before filtering and model fitting)", false);
    setValidFormats_("candidates_out", ListUtils::create<String>("featureXML"));
    registerInputFile_("candidates_in", "<file>", "", "Input file: Feature candidates from a previous run. If set, only feature classification and elution model fitting are carried out, if enabled. Many parameters are ignored.", false, true);
    setValidFormats_("candidates_in", ListUtils::create<String>("featureXML"));

    registerTOPPSubsection_("extract", "Parameters for ion chromatogram extraction");
    registerDoubleOption_("extract:mz_window", "<value>", 10.0, "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)", false);
    setMinFloat_("extract:mz_window", 0.0);
    registerIntOption_("extract:n_isotopes", "<number>", 2, "Number of isotopes to include in each peptide assay.", false);
    setMinInt_("extract:n_isotopes", 2);
    registerDoubleOption_("extract:isotope_pmin", "<value>", 0.0, "Minimum probability for an isotope to be included in the assay for a peptide. If set, this parameter takes precedence over 'extract:n_isotopes'.", false, true);
    setMinFloat_("extract:isotope_pmin", 0.0);
    setMaxFloat_("extract:isotope_pmin", 1.0);
    registerDoubleOption_("extract:rt_quantile", "<value>", 0.95, "Quantile of the RT deviations between aligned internal and external IDs to use for scaling the RT extraction window", false, true);
    setMinFloat_("extract:rt_quantile", 0.0);
    setMaxFloat_("extract:rt_quantile", 1.0);
    registerDoubleOption_("extract:rt_window", "<value>", 0.0, "RT window size (in sec.) for chromatogram extraction. If set, this parameter takes precedence over 'extract:rt_quantile'.", false, true);
    setMinFloat_("extract:rt_window", 0.0);

    registerTOPPSubsection_("detect", "Parameters for detecting features in extracted ion chromatograms");
    registerDoubleOption_("detect:peak_width", "<value>", 60.0, "Expected elution peak width in seconds, for smoothing (Gauss filter). Also determines the RT extration window, unless set explicitly via 'extract:rt_window'.", false);
    setMinFloat_("detect:peak_width", 0.0);
    registerDoubleOption_("detect:min_peak_width", "<value>", 0.2, "Minimum elution peak width. Absolute value in seconds if 1 or greater, else relative to 'peak_width'.", false, true);
    setMinFloat_("detect:min_peak_width", 0.0);
    registerDoubleOption_("detect:signal_to_noise", "<value>", 0.8, "Signal-to-noise threshold for OpenSWATH feature detection", false, true);
    setMinFloat_("detect:signal_to_noise", 0.1);
    registerDoubleOption_("detect:mapping_tolerance", "<value>", 0.0, "RT tolerance (plus/minus) for mapping peptide IDs to features. Absolute value in seconds if 1 or greater, else relative to the RT span of the feature.", false);
    setMinFloat_("detect:mapping_tolerance", 0.0);

    // parameters for SVM classification:
    registerTOPPSubsection_("svm", "Parameters for scoring features using a support vector machine (SVM)");
    registerIntOption_("svm:samples", "<number>", 0, "Number of observations to use for training ('0' for all)", false);
    setMinInt_("svm:samples", 0);
    registerFlag_("svm:no_selection", "By default, roughly the same number of positive and negative observations, with the same intensity distribution, are selected for training. This aims to reduce biases, but also reduces the amount of training data. Set this flag to skip this procedure and consider all available observations (subject to 'svm:samples').");
    registerOutputFile_("svm:xval_out", "<file>", "", "Output file: SVM cross-validation (parameter optimization) results", false);
    setValidFormats_("svm:xval_out", ListUtils::create<String>("csv"));
    Param svm_params;
    svm_params.insert("svm:", SimpleSVM().getParameters());
    registerFullParam_(svm_params);
    registerStringOption_("svm:predictors", "<list>", score_metavalues_, "Names of OpenSWATH scores to use as predictors for the SVM (comma-separated list)", false, true);
    registerDoubleOption_("svm:min_prob", "<value>", 0.0, "Minimum probability of correctness, as predicted by the SVM, required to retain a feature candidate", false, true);
    setMinFloat_("svm:min_prob", 0.0);
    setMaxFloat_("svm:min_prob", 1.0);

    // parameters for model fitting (via ElutionModelFitter):
    registerTOPPSubsection_("model", "Parameters for fitting elution models to features");
    StringList models = ListUtils::create<String>("symmetric,asymmetric,none");
    registerStringOption_("model:type", "<choice>", models[0], "Type of elution model to fit to features", false);
    setValidStrings_("model:type", models);
    Param emf_params;
    emf_params.insert("model:", ElutionModelFitter().getParameters());
    emf_params.remove("model:asymmetric");
    registerFullParam_(emf_params);
  }

  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTrace MassTrace;
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTraces MassTraces;

  // mapping: RT (not necessarily unique) -> pointer to peptide
  typedef multimap<double, PeptideIdentification*> RTMap;
  // mapping: charge -> internal/external: (RT -> pointer to peptide)
  typedef map<Int, pair<RTMap, RTMap> > ChargeMap;
  // mapping: sequence -> charge -> internal/external ID information
  typedef map<AASequence, ChargeMap> PeptideMap;
  // mapping: peptide ref. -> int./ext.: (RT -> pointer to peptide)
  typedef map<String, pair<RTMap, RTMap> > PeptideRefRTMap;

  // region in RT in which a peptide elutes:
  struct RTRegion
  {
    double start, end;
    ChargeMap ids; // internal/external peptide IDs (per charge) in this region
  };

  // predicate for filtering features by overall quality:
  struct FeatureFilterQuality
  {
    bool operator()(const Feature& feature)
    {
      return feature.getOverallQuality() == 0.0;
    }
  } feature_filter_quality_;

  // predicate for filtering features by assigned peptides:
  struct FeatureFilterPeptides
  {
    bool operator()(const Feature& feature)
    {
      return feature.getPeptideIdentifications().empty();
    }
  } feature_filter_peptides_;

  // comparison functor for (unassigned) peptide IDs
  struct PeptideCompare
  {
    bool operator()(const PeptideIdentification& p1,
                    const PeptideIdentification& p2)
    {
      const String& seq1 = p1.getHits()[0].getSequence().toString();
      const String& seq2 = p2.getHits()[0].getSequence().toString();
      if (seq1 == seq2)
      {
        Int charge1 = p1.getHits()[0].getCharge();
        Int charge2 = p2.getHits()[0].getCharge();
        if (charge1 == charge2)
        {
          return p1.getRT() < p2.getRT();
        }
        return charge1 < charge2;
      }
      return seq1 < seq2;
    }
  } peptide_compare_;

  // comparison functor for features
  struct FeatureCompare
  {
    bool operator()(const Feature& f1, const Feature& f2)
    {
      const String& ref1 = f1.getMetaValue("PeptideRef");
      const String& ref2 = f2.getMetaValue("PeptideRef");
      if (ref1 == ref2)
      {
        return f1.getRT() < f2.getRT();
      }
      return ref1 < ref2;
    }
  } feature_compare_;

  PeakMap ms_data_; // input LC-MS data
  PeakMap chrom_data_; // accumulated chromatograms (XICs)
  bool keep_chromatograms_; // keep chromatogram data for output?
  TargetedExperiment library_; // accumulated assays for peptides
  bool keep_library_; // keep assay data for output?
  String score_metavalues_; // names of scores to use as SVM features
  // SVM probability -> number of pos./neg. features (for FDR calculation):
  map<double, pair<Size, Size> > svm_probs_internal_;
  // SVM probabilities for "external" features (for FDR calculation):
  multiset<double> svm_probs_external_;
  Size n_internal_features_; // internal feature counter (for FDR calculation)
  Size n_external_features_; // external feature counter (for FDR calculation)
  // TransformationDescription trafo_; // RT transformation (to range 0-1)
  TransformationDescription trafo_external_; // transform. to external RT scale
  double rt_window_; // RT window width
  double mz_window_; // m/z window width
  bool mz_window_ppm_; // m/z window width is given in PPM (not Da)?
  double isotope_pmin_; // min. isotope probability for peptide assay
  Size n_isotopes_; // number of isotopes for peptide assay
  map<String, double> isotope_probs_; // isotope probabilities of transitions
  double mapping_tolerance_; // RT tolerance for mapping IDs to features
  Size n_parts_; // number of partitions for SVM cross-validation
  Size n_samples_; // number of samples for SVM training
  MRMFeatureFinderScoring feat_finder_; // OpenSWATH feature finder
  ProgressLogger prog_log_;


  /// generate transitions (isotopic traces) for a peptide ion and add them to the library:
  void generateTransitions_(const String& peptide_id, double mz, Int charge,
                            const IsotopeDistribution& iso_dist)
  {
    // go through different isotopes:
    Size counter = 0;
    for (IsotopeDistribution::ConstIterator iso_it = iso_dist.begin();
         iso_it != iso_dist.end(); ++iso_it, ++counter)
    {
      ReactionMonitoringTransition transition;
      String annotation = "i" + String(counter + 1);
      String transition_name = peptide_id + "_" + annotation;

      transition.setNativeID(transition_name);
      transition.setPrecursorMZ(mz);
      transition.setProductMZ(mz + Constants::C13C12_MASSDIFF_U * 
                              float(counter) / charge);
      transition.setLibraryIntensity(iso_it->getIntensity());
      transition.setMetaValue("annotation", annotation);
      transition.setPeptideRef(peptide_id);
      library_.addTransition(transition);
      isotope_probs_[transition_name] = iso_it->getIntensity();
    }
  }


  void addPeptideRT_(TargetedExperiment::Peptide& peptide, double rt)
  {
    TargetedExperiment::RetentionTime te_rt;
    te_rt.setRT(rt);
    te_rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::NORMALIZED;
    peptide.rts.push_back(te_rt);
  }


  /// get regions in which peptide elutes (ideally only one) by clustering RT elution times
  void getRTRegions_(ChargeMap& peptide_data, vector<RTRegion>& rt_regions)
  {
    // use RTs from all charge states here to get a more complete picture:
    vector<double> rts;
    for (ChargeMap::iterator cm_it = peptide_data.begin();
         cm_it != peptide_data.end(); ++cm_it)
    {
      // "internal" IDs:
      for (RTMap::iterator rt_it = cm_it->second.first.begin();
           rt_it != cm_it->second.first.end(); ++rt_it)
      {
        rts.push_back(rt_it->first);
      }
      // "external" IDs:
      for (RTMap::iterator rt_it = cm_it->second.second.begin();
           rt_it != cm_it->second.second.end(); ++rt_it)
      {
        rts.push_back(rt_it->first);
      }
    }
    sort(rts.begin(), rts.end());
    double rt_tolerance = rt_window_ / 2.0;

    for (vector<double>::iterator rt_it = rts.begin(); rt_it != rts.end();
         ++rt_it)
    {
      // create a new region?
      if (rt_regions.empty() || (rt_regions.back().end < *rt_it - rt_tolerance))
      {
        RTRegion region;
        region.start = *rt_it - rt_tolerance;
        // TODO
        // cppcheck-suppress uninitStructMember
        rt_regions.push_back(region);
      }
      rt_regions.back().end = *rt_it + rt_tolerance;
    }

    // sort the peptide IDs into the regions:
    for (ChargeMap::iterator cm_it = peptide_data.begin();
         cm_it != peptide_data.end(); ++cm_it)
    {
      // regions are sorted by RT, as are IDs, so just iterate linearly:
      vector<RTRegion>::iterator reg_it = rt_regions.begin();
      // "internal" IDs:
      for (RTMap::iterator rt_it = cm_it->second.first.begin();
           rt_it != cm_it->second.first.end(); ++rt_it)
      {
        while (rt_it->first > reg_it->end) ++reg_it;
        reg_it->ids[cm_it->first].first.insert(*rt_it);
      }
      reg_it = rt_regions.begin(); // reset to start
      // "external" IDs:
      for (RTMap::iterator rt_it = cm_it->second.second.begin();
           rt_it != cm_it->second.second.end(); ++rt_it)
      {
        while (rt_it->first > reg_it->end) ++reg_it;
        reg_it->ids[cm_it->first].second.insert(*rt_it);
      }
      // ID references no longer needed (now stored in the RT regions):
      cm_it->second.first.clear();
      cm_it->second.second.clear();
    }
  }

  void annotateFeaturesFinalizeAssay_(
    FeatureMap& features, map<Size, vector<PeptideIdentification*> >& feat_ids,
    RTMap& rt_internal)
  {
    set<PeptideIdentification*> assigned_ids;
    if (!feat_ids.empty())
    {
      // find the "best" feature (with the most IDs):
      Size best_index = 0;
      Size best_count = 0;
      for (map<Size, vector<PeptideIdentification*> >::iterator fi_it =
             feat_ids.begin(); fi_it != feat_ids.end(); ++fi_it)
      {
        Size current_index = fi_it->first;
        Size current_count = fi_it->second.size();
        if ((current_count > best_count) ||
            ((current_count == best_count) && // break ties by intensity
             (features[current_index].getIntensity() >
              features[best_index].getIntensity())))
        {
          best_count = current_count;
          best_index = current_index;
        }
      }
      // assign IDs:
      if (best_count > 0)
      {
        // we define the (one) feature with most matching IDs as correct:
        features[best_index].setMetaValue("feature_class", "positive");
        features[best_index].getPeptideIdentifications().resize(best_count);
        for (Size i = 0; i < best_count; ++i)
        {
          features[best_index].getPeptideIdentifications()[i] =
              *(feat_ids[best_index][i]);
        }
        assigned_ids.insert(feat_ids[best_index].begin(),
                            feat_ids[best_index].end());
      }
    }
    // store unassigned IDs from the current RT region:
    for (RTMap::const_iterator rt_it = rt_internal.begin();
         rt_it != rt_internal.end(); ++rt_it)
    {
      if (!assigned_ids.count(rt_it->second))
      {
        const PeptideIdentification& pep_id = *(rt_it->second);
        features.getUnassignedPeptideIdentifications().push_back(pep_id);
      }
    }
    // clean-up:
    feat_ids.clear();
    rt_internal.clear();
  }

  /// annotate identified features with m/z, isotope probabilities, etc.
  void annotateFeatures_(FeatureMap& features, PeptideRefRTMap& ref_rt_map)
  {
    String previous_ref, peptide_ref;
    RTMap transformed_internal;
    Size i = 0;
    map<Size, vector<PeptideIdentification*> > feat_ids;
    for (FeatureMap::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it, ++i)
    {
      feat_it->setMZ(feat_it->getMetaValue("PrecursorMZ"));
      feat_it->setCharge(feat_it->getPeptideIdentifications()[0].getHits()[0].
                         getCharge());
      ensureConvexHulls_(*feat_it);
      // remove "fake" IDs generated by OpenSWATH (they would be removed with
      // a warning when writing output, because of missing protein
      // identification with corresponding identifier):
      feat_it->getPeptideIdentifications().clear();
      // annotate subordinates with theoretical isotope intensities:
      for (vector<Feature>::iterator sub_it =
             feat_it->getSubordinates().begin(); sub_it !=
             feat_it->getSubordinates().end(); ++sub_it)
      {
        String native_id = sub_it->getMetaValue("native_id");
        sub_it->setMetaValue("isotope_probability", isotope_probs_[native_id]);
      }

      peptide_ref = feat_it->getMetaValue("PeptideRef");
      // remove region number, if present:
      Size pos_slash = peptide_ref.rfind('/');
      Size pos_colon = peptide_ref.find(':', pos_slash + 2);
      peptide_ref = peptide_ref.substr(0, pos_colon);

      if (peptide_ref != previous_ref)
      {
        if (!previous_ref.empty())
        {
          annotateFeaturesFinalizeAssay_(
            features, feat_ids, ref_rt_map[previous_ref].first);
        }
        previous_ref = peptide_ref;
      }

      RTMap& rt_internal = ref_rt_map[peptide_ref].first;
      RTMap& rt_external = ref_rt_map[peptide_ref].second;

      if (!rt_internal.empty()) // validate based on internal IDs
      {
        // map IDs to features (based on RT):
        double rt_min = features[i].getMetaValue("leftWidth");
        double rt_max = features[i].getMetaValue("rightWidth");
        if (mapping_tolerance_ > 0.0)
        {
          double abs_tol = mapping_tolerance_;
          if (abs_tol < 1.0)
          {
            abs_tol *= (rt_max - rt_min);
          }
          rt_min -= abs_tol;
          rt_max += abs_tol;
        }
        RTMap::const_iterator lower = rt_internal.lower_bound(rt_min);
        RTMap::const_iterator upper = rt_internal.upper_bound(rt_max);
        int id_count = 0;
        for (; lower != upper; ++lower)
        {
          feat_ids[i].push_back(lower->second);
          ++id_count;
        }
        // "total" only includes IDs from this RT region:
        feat_it->setMetaValue("n_total_ids", rt_internal.size());
        feat_it->setMetaValue("n_matching_ids", id_count);
        if (id_count > 0) // matching IDs -> feature may be correct
        {
          feat_it->setMetaValue("feature_class", "ambiguous");
        }
        else // no matching IDs -> feature is wrong
        {
          feat_it->setMetaValue("feature_class", "negative");
        }
      }
      else // only external IDs -> no validation possible
      {
        feat_it->setMetaValue("n_total_ids", 0);
        feat_it->setMetaValue("n_matching_ids", -1);
        feat_it->setMetaValue("feature_class", "unknown");
        // add "dummy" peptide identification:
        PeptideIdentification id = *(rt_external.begin()->second);
        id.clearMetaInfo();
        id.setMetaValue("FFId_category", "implied");
        id.setRT(feat_it->getRT());
        id.setMZ(feat_it->getMZ());
        // only one peptide hit per ID - see function "addPeptideToMap_":
        PeptideHit& hit = id.getHits()[0];
        hit.clearMetaInfo();
        hit.setScore(0.0);
        feat_it->getPeptideIdentifications().push_back(id);
      }

      // distance from feature to closest peptide ID:
      if (!trafo_external_.getDataPoints().empty())
      {
        // use external IDs if available, otherwise RT-transformed internal IDs
        // (but only compute the transform if necessary, once per assay!):
        if (rt_external.empty() && (transformed_internal.empty() ||
                                    (peptide_ref != previous_ref)))
        {
          transformed_internal.clear();
          for (RTMap::const_iterator it = rt_internal.begin();
               it != rt_internal.end(); ++it)
          {
            double transformed_rt = trafo_external_.apply(it->first);
            RTMap::value_type pair = make_pair(transformed_rt, it->second);
            transformed_internal.insert(transformed_internal.end(), pair);
          }
        }
        const RTMap& rt_ref = (rt_external.empty() ? transformed_internal :
                               rt_external);

        double rt_min = feat_it->getMetaValue("leftWidth");
        double rt_max = feat_it->getMetaValue("rightWidth");
        if (mapping_tolerance_ > 0.0)
        {
          double abs_tol = mapping_tolerance_;
          if (abs_tol < 1.0)
          {
            abs_tol *= (rt_max - rt_min);
          }
          rt_min -= abs_tol;
          rt_max += abs_tol;
        }
        RTMap::const_iterator lower = rt_ref.lower_bound(rt_min);
        RTMap::const_iterator upper = rt_ref.upper_bound(rt_max);
        if (lower != upper) // there's at least one ID within the feature
        {
          feat_it->setMetaValue("rt_delta", 0.0);
        }
        else // check closest ID
        {
          double rt_delta1 = numeric_limits<double>::infinity();
          if (lower != rt_ref.begin())
          {
            rt_delta1 = fabs((--lower)->first - rt_min);
          }
          double rt_delta2 = numeric_limits<double>::infinity();
          if (upper != rt_ref.end())
          {
            rt_delta2 = fabs(upper->first - rt_min);
          }
          feat_it->setMetaValue("rt_delta", min(rt_delta1, rt_delta2));
        }
      }
    }
    // set of features from the last assay:
    annotateFeaturesFinalizeAssay_(features, feat_ids,
                                   ref_rt_map[peptide_ref].first);
    // store unassigned peptide IDs from assays that did not generate any
    // feature candidates:
    for (PeptideRefRTMap::iterator ref_it = ref_rt_map.begin();
         ref_it != ref_rt_map.end(); ++ref_it)
    {
      RTMap& rt_internal = ref_it->second.first;
      if (!rt_internal.empty()) // not cleared by '...FinalizeAssay()'
      {
        for (RTMap::const_iterator rt_it = rt_internal.begin();
             rt_it != rt_internal.end(); ++rt_it)
        {
          const PeptideIdentification& pep_id = *(rt_it->second);
          features.getUnassignedPeptideIdentifications().push_back(pep_id);
        }
      }
    }
  }

  void ensureConvexHulls_(Feature& feature)
  {
    if (feature.getConvexHulls().empty()) // add hulls for mass traces
    {
      double rt_min = feature.getMetaValue("leftWidth");
      double rt_max = feature.getMetaValue("rightWidth");
      for (vector<Feature>::iterator sub_it = feature.getSubordinates().begin();
           sub_it != feature.getSubordinates().end(); ++sub_it)
      {
        double abs_mz_tol = mz_window_ / 2.0;
        if (mz_window_ppm_)
        {
          abs_mz_tol = sub_it->getMZ() * abs_mz_tol * 1.0e-6;
        }
        ConvexHull2D hull;
        hull.addPoint(DPosition<2>(rt_min, sub_it->getMZ() - abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_min, sub_it->getMZ() + abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_max, sub_it->getMZ() - abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_max, sub_it->getMZ() + abs_mz_tol));
        feature.getConvexHulls().push_back(hull);
      }
    }
  }


  void createAssayLibrary_(PeptideMap& peptide_map, PeptideRefRTMap& ref_rt_map)
  {
    set<String> protein_accessions;

    for (PeptideMap::iterator pm_it = peptide_map.begin();
         pm_it != peptide_map.end(); ++pm_it)
    {
      TargetedExperiment::Peptide peptide;

      const AASequence& seq = pm_it->first;
      LOG_DEBUG << "\nPeptide: " << seq.toString() << endl;
      peptide.sequence = seq.toString();
      // @NOTE: Technically, "TargetedExperiment::Peptide" stores the unmodified
      // sequence and the modifications separately. Unfortunately, creating the
      // modifications vector is complex and there is currently no convenient
      // conversion function (see "TargetedExperimentHelper::getAASequence" for
      // the reverse conversion). However, "Peptide" is later converted to
      // "OpenSwath::LightPeptide" anyway, and this is done via "AASequence"
      // (see "OpenSwathDataAccessHelper::convertTargetedPeptide"). So for our
      // purposes it works to just store the sequence including modifications in
      // "Peptide".

      // keep track of protein accessions:
      set<String> current_accessions;
      const pair<RTMap, RTMap>& pair = pm_it->second.begin()->second;
      const PeptideHit& hit = (pair.first.empty() ?
                               pair.second.begin()->second->getHits()[0] :
                               pair.first.begin()->second->getHits()[0]);
      current_accessions = hit.extractProteinAccessionsSet();
      protein_accessions.insert(current_accessions.begin(),
                                current_accessions.end());
      // missing protein accession would crash OpenSWATH algorithms:
      if (current_accessions.empty())
      {
        current_accessions.insert("not_available");
      }
      peptide.protein_refs = vector<String>(current_accessions.begin(),
                                            current_accessions.end());

      // get isotope distribution for peptide:
      Size n_isotopes = (isotope_pmin_ > 0.0) ? 10 : n_isotopes_;
      IsotopeDistribution iso_dist = 
        seq.getFormula(Residue::Full, 0).getIsotopeDistribution(new CoarseIsotopeDistribution(n_isotopes));
      if (isotope_pmin_ > 0.0)
      {
        iso_dist.trimLeft(isotope_pmin_);
        iso_dist.trimRight(isotope_pmin_);
        iso_dist.renormalize();
      }

      // get regions in which peptide elutes (ideally only one):
      vector<RTRegion> rt_regions;
      getRTRegions_(pm_it->second, rt_regions);
      LOG_DEBUG << "Found " << rt_regions.size() << " RT region(s)." << endl;

      // go through different charge states:
      for (ChargeMap::const_iterator cm_it = pm_it->second.begin();
           cm_it != pm_it->second.end(); ++cm_it)
      {
        Int charge = cm_it->first;
        double mz = seq.getMonoWeight(Residue::Full, charge) / charge;
        LOG_DEBUG << "Charge: " << charge << " (m/z: " << mz << ")" << endl;
        peptide.setChargeState(charge);
        String peptide_id = peptide.sequence + "/" + String(charge);

        // we want to detect one feature per peptide and charge state - if there
        // are multiple RT regions, group them together:
        peptide.setPeptideGroupLabel(peptide_id);
        peptide.rts.clear();
        Size counter = 0;
        // accumulate IDs over multiple regions:
        RTMap& internal_ids = ref_rt_map[peptide_id].first;
        RTMap& external_ids = ref_rt_map[peptide_id].second;
        for (vector<RTRegion>::iterator reg_it = rt_regions.begin();
             reg_it != rt_regions.end(); ++reg_it)
        {
          if (reg_it->ids.count(charge))
          {
            LOG_DEBUG << "Region " << counter + 1 << " (RT: "
                      << float(reg_it->start) << "-" << float(reg_it->end)
                      << ", size " << float(reg_it->end - reg_it->start) << ")"
                      << endl;

            peptide.id = peptide_id;
            if (rt_regions.size() > 1) peptide.id += ":" + String(++counter);

            // store beginning and end of RT region:
            peptide.rts.clear();
            addPeptideRT_(peptide, reg_it->start);
            addPeptideRT_(peptide, reg_it->end);
            library_.addPeptide(peptide);
            generateTransitions_(peptide.id, mz, charge, iso_dist);
          }
          internal_ids.insert(reg_it->ids[charge].first.begin(),
                              reg_it->ids[charge].first.end());
          external_ids.insert(reg_it->ids[charge].second.begin(),
                              reg_it->ids[charge].second.end());
        }
      }
    }

    // add proteins to library:
    for (set<String>::iterator acc_it = protein_accessions.begin(); 
         acc_it != protein_accessions.end(); ++acc_it)
    {
      TargetedExperiment::Protein protein;
      protein.id = *acc_it;
      library_.addProtein(protein);
    }
  }


  void addPeptideToMap_(PeptideIdentification& peptide, PeptideMap& peptide_map,
                        bool external = false)
  {
    if (peptide.getHits().empty()) return;
    peptide.sort();
    PeptideHit& hit = peptide.getHits()[0];
    peptide.getHits().resize(1);
    Int charge = hit.getCharge();
    double rt = peptide.getRT();
    RTMap::value_type pair = make_pair(rt, &peptide);
    if (!external)
    {
      peptide_map[hit.getSequence()][charge].first.insert(pair);
    }
    else
    {
      peptide_map[hit.getSequence()][charge].second.insert(pair);
    }
  }


  void checkNumObservations_(Size n_pos, Size n_neg, const String& note = "")
  {
    if (n_pos < n_parts_)
    {
      String msg = "Not enough positive observations for " + 
        String(n_parts_) + "-fold cross-validation" + note + ".";
      throw Exception::MissingInformation(__FILE__, __LINE__, 
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    if (n_neg < n_parts_)
    {
      String msg = "Not enough negative observations for " + 
        String(n_parts_) + "-fold cross-validation" + note + ".";
      throw Exception::MissingInformation(__FILE__, __LINE__, 
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
  }


  void getUnbiasedSample_(const multimap<double, pair<Size, bool> >& valid_obs,
                          map<Size, Int>& training_labels)
  {
    // Create an unbiased training sample:
    // - same number of pos./neg. observations (approx.),
    // - same intensity distribution of pos./neg. observations.
    // We use a sliding window over the set of observations, ordered by
    // intensity. At each step, we examine the proportion of both pos./neg.
    // observations in the window and select the middle element with according
    // probability. (We use an even window size, to cover the ideal case where
    // the two classes are balanced.)
    const Size window_size = 8;
    const Size half_win_size = window_size / 2;
    if (valid_obs.size() < half_win_size + 1)
    {
      String msg = "Not enough observations for intensity-bias filtering.";
      throw Exception::MissingInformation(__FILE__, __LINE__, 
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    srand(time(nullptr)); // seed random number generator
    Size n_obs[2] = {0, 0}; // counters for neg./pos. observations
    Size counts[2] = {0, 0}; // pos./neg. counts in current window
    // iterators to begin, middle and past-the-end of sliding window:
    multimap<double, pair<Size, bool> >::const_iterator begin, middle, end;
    begin = middle = end = valid_obs.begin();
    // initialize ("middle" is at beginning of sequence, so no full window):
    for (Size i = 0; i <= half_win_size; ++i, ++end)
    {
      ++counts[end->second.second]; // increase counter for pos./neg. obs.
    }
    // "i" is the index of one of the two middle values of the sliding window:
    // - in the left half of the sequence, "i" is left-middle,
    // - in the right half of the sequence, "i" is right-middle.
    // The counts are updated as "i" and the sliding window move to the right.
    for (Size i = 0; i < valid_obs.size(); ++i, ++middle)
    {
      // if count for either class is zero, we don't select anything:
      if ((counts[0] > 0) && (counts[1] > 0))
      {
        // probability thresholds for neg./pos. observations:
        double thresholds[2] = {counts[1] / float(counts[0]),
                                counts[0] / float(counts[1])};
        // check middle values:
        double rnd = rand() / double(RAND_MAX); // random num. in range 0-1
        if (rnd < thresholds[middle->second.second])
        {
          training_labels[middle->second.first] = Int(middle->second.second);
          ++n_obs[middle->second.second];
        }
      }
      // update sliding window and class counts;
      // when we reach the middle of the sequence, we keep the window in place
      // for one step, to change from "left-middle" to "right-middle":
      if (i != valid_obs.size() / 2)
      {
        // only move "begin" when "middle" has advanced far enough:
        if (i > half_win_size)
        {
          --counts[begin->second.second];
          ++begin;
        }
        // don't increment "end" beyond the defined range:
        if (end != valid_obs.end())
        {
          ++counts[end->second.second];
          ++end;
        }
      }
    }
    checkNumObservations_(n_obs[1], n_obs[0], " after bias filtering");
  }


  void getRandomSample_(map<Size, Int>& training_labels)
  {
    // @TODO: can this be done with less copying back and forth of data?
    // Pick a random subset of size "n_samples_" for training: Shuffle the whole
    // sequence, then select the first "n_samples_" elements.
    vector<Size> selection;
    selection.reserve(training_labels.size());
    for (map<Size, Int>::iterator it = training_labels.begin();
         it != training_labels.end(); ++it)
    {
      selection.push_back(it->first);
    }
    random_shuffle(selection.begin(), selection.end());
    // However, ensure that at least "n_parts_" pos./neg. observations are
    // included (for cross-validation) - there must be enough, otherwise
    // "checkNumObservations_" would have thrown an error. To this end, move
    // "n_parts_" pos. observations to the beginning of sequence, followed by
    // "n_parts_" neg. observations (pos. first - see reason below):
    Size n_obs[2] = {0, 0}; // counters for neg./pos. observations
    for (Int label = 1; label >= 0; --label)
    {
      for (Size i = n_obs[1]; i < selection.size(); ++i)
      {
        Size obs_index = selection[i];
        if (training_labels[obs_index] == label)
        {
          swap(selection[i], selection[n_obs[label]]);
          ++n_obs[label];
        }
        if (n_obs[label] == n_parts_) break;
      }
    }
    selection.resize(n_samples_);
    // copy the selected subset back:
    map<Size, Int> temp;
    for (vector<Size>::iterator it = selection.begin(); it != selection.end();
         ++it)
    {
      temp[*it] = training_labels[*it];
    }
    training_labels.swap(temp);
  }


  void classifyFeatures_(FeatureMap& features)
  {
    if (features.empty()) return;

    // get predictors for SVM:
    vector<String> predictor_names =
      ListUtils::create<String>(getStringOption_("svm:predictors"));
    if (features[0].metaValueExists("rt_delta")) // include RT feature
    {
      predictor_names.push_back("rt_delta");
    }
    // values for all featues per predictor (this way around to simplify scaling
    // of predictors):
    SimpleSVM::PredictorMap predictors;
    for (vector<String>::iterator pred_it = predictor_names.begin();
         pred_it != predictor_names.end(); ++pred_it)
    {
      predictors[*pred_it].reserve(features.size());
      for (FeatureMap::Iterator feat_it = features.begin(); 
           feat_it < features.end(); ++feat_it)
      {
        if (!feat_it->metaValueExists(*pred_it))
        {
          LOG_ERROR << "Meta value '" << *pred_it << "' missing for feature '"
                    << feat_it->getUniqueId() << "'" << endl;
          predictors.erase(*pred_it);
          break;
        }
        predictors[*pred_it].push_back(feat_it->getMetaValue(*pred_it));
      }
    }

    // get labels for SVM:
    map<Size, Int> training_labels;
    bool no_selection = getFlag_("svm:no_selection");
    // mapping (for bias correction): intensity -> (index, positive?)
    multimap<double, pair<Size, bool> > valid_obs;
    Size n_obs[2] = {0, 0}; // counters for neg./pos. observations
    for (Size feat_index = 0; feat_index < features.size(); ++feat_index)
    {
      String feature_class = features[feat_index].getMetaValue("feature_class");
      Int label = -1;
      if (feature_class == "positive") label = 1;
      else if (feature_class == "negative") label = 0;

      if (label != -1)
      {
        ++n_obs[label];
        if (!no_selection)
        {
          double intensity = features[feat_index].getIntensity();
          valid_obs.insert(make_pair(intensity, make_pair(feat_index,
                                                          bool(label))));
        }
        else
        {
          training_labels[feat_index] = label;
        }
      }
    }
    checkNumObservations_(n_obs[1], n_obs[0]);

    if (!no_selection) getUnbiasedSample_(valid_obs, training_labels);

    if (n_samples_ > 0) // limited number of samples for training
    {
      if (training_labels.size() < n_samples_)
      {
        LOG_WARN << "Warning: There are only " << training_labels.size()
                 << " valid observations for training." << endl;
      }
      else if (training_labels.size() > n_samples_)
      {
        getRandomSample_(training_labels);
      }
    }

    SimpleSVM svm;
    // set (only) the relevant parameters:
    Param svm_params = svm.getParameters();
    Logger::LogStream no_log; // suppress warnings about additional parameters
    svm_params.update(getParam_().copy("svm:", true), false, no_log);
    svm.setParameters(svm_params);
    svm.setup(predictors, training_labels);
    String xval_out = getStringOption_("svm:xval_out");
    if (!xval_out.empty()) svm.writeXvalResults(xval_out);
    if ((debug_level_ > 0) && String(svm_params.getValue("kernel")) == "linear")
    {
      map<String, double> feature_weights;
      svm.getFeatureWeights(feature_weights);
      LOG_DEBUG << "SVM feature weights:" << endl;
      for (map<String, double>::iterator it = feature_weights.begin();
           it != feature_weights.end(); ++it)
      {
        LOG_DEBUG << "- " << it->first << ": " << it->second << endl;
      }
    }

    vector<SimpleSVM::Prediction> predictions;
    svm.predict(predictions);
    OPENMS_POSTCONDITION(predictions.size() == features.size(), 
                         "SVM predictions for all features expected");
    for (Size i = 0; i < features.size(); ++i)
    {
      features[i].setMetaValue("predicted_class", predictions[i].label);
      double prob_positive = predictions[i].probabilities[1];
      features[i].setMetaValue("predicted_probability", prob_positive);
      // @TODO: store previous (OpenSWATH) overall quality in a meta value?
      features[i].setOverallQuality(prob_positive);
    }
  }


  void filterFeaturesFinalizeAssay_(Feature& best_feature, double best_quality,
                                    const double quality_cutoff)
  {
    const String& feature_class = best_feature.getMetaValue("feature_class");
    if (feature_class == "positive") // true positive prediction
    {
      svm_probs_internal_[best_quality].first++;
    }
    else if ((feature_class == "negative")  || // false positive prediction
             (feature_class == "ambiguous")) // let's be strict about this
    {
      svm_probs_internal_[best_quality].second++;
    }
    else if (feature_class == "unknown")
    {
      svm_probs_external_.insert(best_quality);
      if (best_quality >= quality_cutoff) 
      {
        best_feature.setOverallQuality(best_quality);
        ++n_external_features_;
      }
    }
  }


  void filterFeatures_(FeatureMap& features, bool classified)
  {
    if (features.empty()) return;

    if (classified)
    {
      // Remove features with class "negative" or "ambiguous", keep "positive".
      // For class "unknown", for every assay (meta value "PeptideRef"), keep
      // the feature with highest "predicted_probability" (= overall quality),
      // subject to the "svm:min_prob" threshold.
      // We mark features for removal by setting their overall quality to zero.
      n_internal_features_ = n_external_features_ = 0;
      FeatureMap::Iterator best_it = features.begin();
      const double quality_cutoff = getDoubleOption_("svm:min_prob");
      double best_quality = 0.0;
      String previous_ref;
      for (FeatureMap::Iterator it = features.begin(); it != features.end();
           ++it)
      {
        // features from same assay (same "PeptideRef") appear consecutively;
        // if this is a new assay, finalize the previous one:
        String peptide_ref = it->getMetaValue("PeptideRef");
        // remove region number, if present:
        Size pos_slash = peptide_ref.rfind('/');
        Size pos_colon = peptide_ref.find(':', pos_slash + 2);
        peptide_ref = peptide_ref.substr(0, pos_colon);

        if (peptide_ref != previous_ref)
        {
          if (!previous_ref.empty())
          {
            filterFeaturesFinalizeAssay_(*best_it, best_quality,
                                         quality_cutoff);
            best_quality = 0.0;
          }
          previous_ref = peptide_ref;
        }

        // update qualities:
        if ((it->getOverallQuality() > best_quality) ||
            // break ties by intensity:
            ((it->getOverallQuality() == best_quality) &&
             (it->getIntensity() > best_it->getIntensity())))
        {
          best_it = it;
          best_quality = it->getOverallQuality();
        }
        if (it->getMetaValue("feature_class") == "positive")
        {
          n_internal_features_++;
        }
        else
        {
          it->setOverallQuality(0.0); // gets overwritten for "best" candidate
        }
      }
      // set of features from the last assay:
      filterFeaturesFinalizeAssay_(*best_it, best_quality, quality_cutoff);

      features.erase(remove_if(features.begin(), features.end(),
                               feature_filter_quality_), features.end());
    }
    else
    {
      features.erase(remove_if(features.begin(), features.end(),
                               feature_filter_peptides_), features.end());
    }
  }


  void calculateFDR_(FeatureMap& features)
  {
    // cumulate the true/false positive counts, in decreasing probability order:
    Size n_false = 0, n_true = 0;
    for (map<double, pair<Size, Size> >::reverse_iterator prob_it =
           svm_probs_internal_.rbegin(); prob_it != svm_probs_internal_.rend();
         ++prob_it)
    {
      n_true += prob_it->second.first;
      n_false += prob_it->second.second;
      prob_it->second.first = n_true;
      prob_it->second.second = n_false;
    }

    // print FDR for features that made the cut-off:
    const double min_prob = getDoubleOption_("svm:min_prob");
    map<double, pair<Size, Size> >::iterator prob_it =
      svm_probs_internal_.lower_bound(min_prob);
    if (prob_it != svm_probs_internal_.end())
    {
      float fdr = float(prob_it->second.second) / (prob_it->second.first +
                                                   prob_it->second.second);
      LOG_INFO << "Estimated FDR of features detected based on 'external' IDs: "
               << fdr * 100.0 << "%" << endl;
      fdr = (fdr * n_external_features_) / (n_external_features_ + 
                                            n_internal_features_);
      LOG_INFO << "Estimated FDR of all detected features: " << fdr * 100.0
               << "%" << endl;
    }

    // calculate q-values:
    vector<double> qvalues;
    qvalues.reserve(svm_probs_internal_.size());
    double min_fdr = 1.0;
    for (prob_it = svm_probs_internal_.begin();
         prob_it != svm_probs_internal_.end(); ++prob_it)
    {
      double fdr = double(prob_it->second.second) / (prob_it->second.first +
                                                     prob_it->second.second);
      if (fdr < min_fdr) min_fdr = fdr;
      qvalues.push_back(min_fdr);
    }
    // record only probabilities where q-value changes:
    vector<double> fdr_probs, fdr_qvalues;
    vector<double>::iterator qv_it = qvalues.begin();
    double previous_qvalue = -1.0;
    for (prob_it = svm_probs_internal_.begin();
         prob_it != svm_probs_internal_.end(); ++prob_it, ++qv_it)
    {
      if (*qv_it != previous_qvalue)
      {
        fdr_probs.push_back(prob_it->first);
        fdr_qvalues.push_back(*qv_it);
        previous_qvalue = *qv_it;
      }
    }
    features.setMetaValue("FDR_probabilities", fdr_probs);
    features.setMetaValue("FDR_qvalues_raw", fdr_qvalues);

    // FDRs are estimated from "internal" features, but apply only to "external"
    // ones. "Internal" features are considered "correct" by definition.
    // We need to adjust the q-values to take this into account:
    multiset<double>::reverse_iterator ext_it = svm_probs_external_.rbegin();
    Size external_count = 0;
    for (Int i = fdr_probs.size() - 1; i >= 0; --i)
    {
      double cutoff = fdr_probs[i];
      while ((ext_it != svm_probs_external_.rend()) && (*ext_it >= cutoff))
      {
        ++external_count;
        ++ext_it;
      }
      fdr_qvalues[i] = (fdr_qvalues[i] * external_count) /
        (external_count + n_internal_features_);
    }
    features.setMetaValue("FDR_qvalues_corrected", fdr_qvalues);

    // @TODO: should we use "1 - qvalue" as overall quality for features?
    // assign q-values to features:
    for (FeatureMap::iterator feat_it = features.begin(); 
         feat_it != features.end(); ++feat_it)
    {
      if (feat_it->getMetaValue("feature_class") == "positive")
      {
        feat_it->setMetaValue("q-value", 0.0);
      }
      else
      {
        double prob = feat_it->getOverallQuality();
        // find highest FDR prob. that is less-or-equal to the feature prob.:
        vector<double>::iterator pos = upper_bound(fdr_probs.begin(),
                                                   fdr_probs.end(), prob);
        if (pos != fdr_probs.begin()) --pos;
        Size dist = distance(fdr_probs.begin(), pos);
        feat_it->setMetaValue("q-value", fdr_qvalues[dist]);
      }
    }
  }


  ExitCodes main_(int, const char**) override
  {
    FeatureMap features;
    PeptideMap peptide_map;
    Size n_internal_peps = 0, n_external_peps = 0;

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String out = getStringOption_("out");
    String candidates_out = getStringOption_("candidates_out");
    n_parts_ = getIntOption_("svm:xval");
    n_samples_ = getIntOption_("svm:samples");
    String elution_model = getStringOption_("model:type");
    bool with_external_ids = false;
    prog_log_.setLogType(log_type_);

    String candidates_in = getStringOption_("candidates_in");
    if (candidates_in.empty())
    {
      String in = getStringOption_("in");
      String id = getStringOption_("id");
      String id_ext = getStringOption_("id_ext");
      String lib_out = getStringOption_("lib_out");
      String chrom_out = getStringOption_("chrom_out");
      rt_window_ = getDoubleOption_("extract:rt_window");
      mz_window_ = getDoubleOption_("extract:mz_window");
      mz_window_ppm_ = mz_window_ >= 1;
      isotope_pmin_ = getDoubleOption_("extract:isotope_pmin");
      n_isotopes_ = getIntOption_("extract:n_isotopes");
      double peak_width = getDoubleOption_("detect:peak_width");
      double min_peak_width = getDoubleOption_("detect:min_peak_width");
      double signal_to_noise = getDoubleOption_("detect:signal_to_noise");
      mapping_tolerance_ = getDoubleOption_("detect:mapping_tolerance");

      if ((n_samples_ > 0) && (n_samples_ < 2 * n_parts_))
      {
        String msg = "Sample size of " + String(n_samples_) +
          " (parameter 'svm:samples') is not enough for " + String(n_parts_) +
          "-fold cross-validation (parameter 'svm:xval').";
        throw Exception::InvalidParameter(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION, msg);
      }

      //-------------------------------------------------------------
      // load input
      //-------------------------------------------------------------
      LOG_INFO << "Loading input data..." << endl;
      MzMLFile mzml;
      mzml.setLogType(log_type_);
      mzml.getOptions().addMSLevel(1);
      mzml.load(in, ms_data_);

      // initialize algorithm classes needed later:
      Param params = feat_finder_.getParameters();
      params.setValue("stop_report_after_feature", -1); // return all features
      params.setValue("Scores:use_rt_score", "false"); // RT may not be reliable
      if ((elution_model != "none") || (!candidates_out.empty()))
      {
        params.setValue("write_convex_hull", "true");
      }
      if (min_peak_width < 1.0) min_peak_width *= peak_width;
      params.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width",
                      peak_width);
      params.setValue("TransitionGroupPicker:min_peak_width", min_peak_width);
      // disabling the signal-to-noise threshold (setting the parameter to zero)
      // totally breaks the OpenSWATH feature detection (no features found)!
      params.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise",
                      signal_to_noise);
      params.setValue("TransitionGroupPicker:recalculate_peaks", "true");
      params.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
      params.setValue("TransitionGroupPicker:PeakPickerMRM:method",
                      "corrected");
      feat_finder_.setParameters(params);
      feat_finder_.setLogType(ProgressLogger::NONE);
      feat_finder_.setStrictFlag(false);

      // "internal" IDs:
      vector<PeptideIdentification> peptides;
      vector<ProteinIdentification> proteins;
      IdXMLFile().load(id, proteins, peptides);
      // "external" IDs:
      vector<PeptideIdentification> peptides_ext;
      vector<ProteinIdentification> proteins_ext;
      double rt_uncertainty = 0.0;
      if (!id_ext.empty())
      {
        with_external_ids = true;
        IdXMLFile().load(id_ext, proteins_ext, peptides_ext);
        // align internal and external IDs to estimate RT shifts:
        MapAlignmentAlgorithmIdentification aligner;
        aligner.setReference(peptides_ext); // go from interal to external scale
        vector<vector<PeptideIdentification> > aligner_peptides(1, peptides);
        vector<TransformationDescription> aligner_trafos;

        LOG_INFO << "Realigning internal and external IDs...";
        aligner.align(aligner_peptides, aligner_trafos);
        trafo_external_ = aligner_trafos[0];
        vector<double> aligned_diffs;
        trafo_external_.getDeviations(aligned_diffs);
        double quantile = getDoubleOption_("extract:rt_quantile");
        int index = max(0, int(quantile * aligned_diffs.size()) - 1);
        rt_uncertainty = aligned_diffs[index];
        try
        {
          aligner_trafos[0].fitModel("lowess");
          trafo_external_ = aligner_trafos[0];
        }
        catch (Exception::BaseException& e)
        {
          LOG_ERROR << "Error: Failed to align RTs of internal/external peptides. RT information will not be considered in the SVM classification. The original error message was:\n" << e.what() << endl;
        }
      }
      if (rt_window_ == 0.0)
      {
        // calculate RT window based on other parameters and alignment quality:
        double map_tol = mapping_tolerance_;
        if (map_tol < 1.0) map_tol *= (2 * peak_width); // relative tolerance
        rt_window_ = (rt_uncertainty + 2 * peak_width + map_tol) * 2;
        LOG_INFO << "RT window size calculated as " << rt_window_ << " seconds."
                 << endl;
      }

      //-------------------------------------------------------------
      // prepare peptide map
      //-------------------------------------------------------------
      LOG_INFO << "Preparing mapping of peptide data..." << endl;
      for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
           pep_it != peptides.end(); ++pep_it)
      {
        addPeptideToMap_(*pep_it, peptide_map);
        pep_it->setMetaValue("FFId_category", "internal");
      }
      n_internal_peps = peptide_map.size();
      for (vector<PeptideIdentification>::iterator pep_it =
             peptides_ext.begin(); pep_it != peptides_ext.end(); ++pep_it)
      {
        addPeptideToMap_(*pep_it, peptide_map, true);
        pep_it->setMetaValue("FFId_category", "external");
      }
      n_external_peps = peptide_map.size() - n_internal_peps;

      //-------------------------------------------------------------
      // run feature detection
      //-------------------------------------------------------------
      keep_library_ = !lib_out.empty();
      keep_chromatograms_ = !chrom_out.empty();

      LOG_INFO << "Creating assay library..." << endl;
      PeptideRefRTMap ref_rt_map;
      createAssayLibrary_(peptide_map, ref_rt_map);

      LOG_INFO << "Extracting chromatograms..." << endl;
      ChromatogramExtractor extractor;
      // extractor.setLogType(ProgressLogger::NONE);
      vector<OpenSwath::ChromatogramPtr> chrom_temp;
      vector<ChromatogramExtractor::ExtractionCoordinates> coords;
      extractor.prepare_coordinates(chrom_temp, coords, library_,
                                    numeric_limits<double>::quiet_NaN(), false);

      boost::shared_ptr<PeakMap> shared = boost::make_shared<PeakMap>(ms_data_);
      OpenSwath::SpectrumAccessPtr spec_temp =
        SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(shared);
      extractor.extractChromatograms(spec_temp, chrom_temp, coords, mz_window_,
                                     mz_window_ppm_, "tophat");
      extractor.return_chromatogram(chrom_temp, coords, library_, (*shared)[0],
                                    chrom_data_.getChromatograms(), false);

      LOG_DEBUG << "Extracted " << chrom_data_.getNrChromatograms()
                << " chromatogram(s)." << endl;

      LOG_INFO << "Detecting chromatographic peaks..." << endl;
      Log_info.remove(cout); // suppress status output from OpenSWATH
      feat_finder_.pickExperiment(chrom_data_, features, library_,
                                  TransformationDescription(), ms_data_);
      Log_info.insert(cout);
      LOG_INFO << "Found " << features.size() << " feature candidates in total."
               << endl;
      ms_data_.reset(); // not needed anymore, free up the memory

      // complete feature annotation:
      annotateFeatures_(features, ref_rt_map);

      // write auxiliary output:
      if (keep_library_)
      {
        TraMLFile().store(lib_out, library_);
      }
      if (keep_chromatograms_)
      {
        addDataProcessing_(chrom_data_,
                           getProcessingInfo_(DataProcessing::FILTERING));
        MzMLFile().store(chrom_out, chrom_data_);
        chrom_data_.clear(true);
      }

      features.setProteinIdentifications(proteins);
      // add external IDs (if any):
      features.getProteinIdentifications().insert(
        features.getProteinIdentifications().end(), proteins_ext.begin(),
        proteins_ext.end());
      features.getUnassignedPeptideIdentifications().insert(
        features.getUnassignedPeptideIdentifications().end(),
        peptides_ext.begin(), peptides_ext.end());

      features.ensureUniqueId();
      addDataProcessing_(features,
                         getProcessingInfo_(DataProcessing::QUANTITATION));
    }
    else
    {
      //-------------------------------------------------------------
      // load feature candidates
      //-------------------------------------------------------------
      LOG_INFO << "Reading feature candidates from a previous run..." << endl;
      FeatureXMLFile().load(candidates_in, features);
      LOG_INFO << "Found " << features.size() << " feature candidates in total."
               << endl;
      with_external_ids = (!features.empty() && 
                           features[0].metaValueExists("predicted_class"));

      // extract ID information for statistics:
      set<AASequence> internal_seqs;
      for (vector<PeptideIdentification>::iterator pep_it =
             features.getUnassignedPeptideIdentifications().begin(); pep_it !=
             features.getUnassignedPeptideIdentifications().end(); ++pep_it)
      {
        const AASequence& seq = pep_it->getHits()[0].getSequence();
        if (pep_it->getMetaValue("FFId_category") == "internal")
        {
          internal_seqs.insert(seq);
        }
        peptide_map[seq];
      }
      for (FeatureMap::Iterator feat_it = features.begin();
           feat_it != features.end(); ++feat_it)
      {
        if (feat_it->getPeptideIdentifications().empty()) continue;
        const PeptideIdentification& pep_id =
          feat_it->getPeptideIdentifications()[0];
        const AASequence& seq = pep_id.getHits()[0].getSequence();
        if (pep_id.getMetaValue("FFId_category") == "internal")
        {
          internal_seqs.insert(seq);
        }
        peptide_map[seq];
      }
      n_internal_peps = internal_seqs.size();
      n_external_peps = peptide_map.size() - internal_seqs.size();
    }

    // sort everything:
    sort(features.getUnassignedPeptideIdentifications().begin(),
         features.getUnassignedPeptideIdentifications().end(),
         peptide_compare_);
    sort(features.begin(), features.end(), feature_compare_);

    // don't do SVM stuff unless we have external data to apply the model to:
    if (with_external_ids) classifyFeatures_(features);

    if (!candidates_out.empty()) // store feature candidates
    {
      FeatureXMLFile().store(candidates_out, features);
    }

    filterFeatures_(features, with_external_ids);
    LOG_INFO << features.size() << " features left after filtering." << endl;

    if (!svm_probs_internal_.empty()) calculateFDR_(features);

    if (elution_model != "none")
    {
      ElutionModelFitter emf;
      Param emf_params = getParam_().copy("model:", true);
      emf_params.remove("type");
      emf_params.setValue("asymmetric",
                          (elution_model == "asymmetric") ? "true" : "false");
      emf.setParameters(emf_params);
      emf.fitElutionModels(features);
    }
    else if (!candidates_out.empty()) // hulls not needed, remove them
    {
      for (FeatureMap::Iterator feat_it = features.begin(); 
           feat_it != features.end(); ++feat_it)
      {
        for (vector<Feature>::iterator sub_it =
               feat_it->getSubordinates().begin(); sub_it != 
               feat_it->getSubordinates().end(); ++sub_it)
        {
          sub_it->getConvexHulls().clear();
        }
      }
    }

    //-------------------------------------------------------------
    // write output
    //-------------------------------------------------------------

    LOG_INFO << "Writing final results..." << endl;
    FeatureXMLFile().store(out, features);

    //-------------------------------------------------------------
    // statistics
    //-------------------------------------------------------------

    // same peptide sequence may be quantified based on internal and external
    // IDs if charge states differ!
    set<AASequence> quantified_internal, quantified_all;
    for (FeatureMap::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it)
    {
      const PeptideIdentification& pep_id =
        feat_it->getPeptideIdentifications()[0];
      const AASequence& seq = pep_id.getHits()[0].getSequence();
      if (feat_it->getIntensity() > 0.0)
      {
        quantified_all.insert(seq);
        if (pep_id.getMetaValue("FFId_category") == "internal")
        {
          quantified_internal.insert(seq);
        }
      }
    }
    Size n_quant_external = quantified_all.size() - quantified_internal.size();
    // If internal and external IDs for a peptide map to different RT regions,
    // it is possible that there is a quantification from the "external" region,
    // but not from the "internal" region (no matching feature) - therefore the
    // number of "missing" external peptides can be negative!
    Int n_missing_external = Int(n_external_peps) - n_quant_external;
    LOG_INFO << "\nSummary statistics (counting distinct peptides including "
      "PTMs):\n"
             << peptide_map.size() << " peptides identified ("
             << n_internal_peps << " internal, " << n_external_peps
             << " additional external)\n"
             << quantified_all.size() << " peptides with features ("
             << quantified_internal.size() << " internal, "
             << n_quant_external << " external)\n"
             << peptide_map.size() - quantified_all.size()
             << " peptides without features ("
             << n_internal_peps - quantified_internal.size() << " internal, "
             << n_missing_external << " external)\n" << endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
