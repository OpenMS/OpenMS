// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/SVM/SimpleSVM.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
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

   This tool detects quantitative features in MS1 data based on information from peptide identifications (derived from MS2 spectra). It uses algorithms for targeted data analysis from the OpenSWATH pipeline.

   The aim is to detect features that enable the quantification of (ideally) all peptides in the identification input. This is based on the following principle: When a high-confidence identification (ID) of a peptide was made based on an MS2 spectrum from a certain (precursor) position in the LC-MS map, this indicates that the particular peptide is present at that position, so a feature for it should be detectable there.

   @note It is important that only high-confidence (i.e. reliable) peptide identifications are used as input!

   Targeted data analysis on the MS1 level uses OpenSWATH algorithms and follows roughly the steps outlined below.

   <B>Use of inferred ("external") IDs</B>

   The situation becomes more complicated when several LC-MS/MS runs from related samples of a label-free experiment are considered. In order to quantify a larger fraction of the peptides/proteins in the samples, it is desirable to infer peptide identifications across runs. Ideally, all peptides identified in any of the runs should be quantified in each and every run. However, for feature detection of inferred ("external") IDs, the following problems arise: First, retention times may be shifted between the run being quantified and the run that gave rise to the ID. Such shifts can be corrected (see @ref TOPP_MapAlignerIdentification), but only to an extent. Thus, the RT location of the inferred ID may not necessarily lie within the RT range of the correct feature. Second, since the peptide in question was not directly identified in the run being quantified, it may not actually be present in detectable amounts in that sample, e.g. due to differential regulation of the corresponding protein. There is thus a risk of introducing false-positive features.

   FeatureFinderIdentification deals with these challenges by explicitly distinguishing between internal IDs (derived from the LC-MS/MS run being quantified) and external IDs (inferred from related runs). Features derived from internal IDs give rise to a training dataset for an SVM classifier. The SVM is then used to predict whether features derived from external IDs are correct or incorrect (and which of several candidates is the most likely to be correct). See steps 4 and 5 below for more details.

   <B>1. Assay generation</B>

   Feature detection is based on assays for identified peptides, each of which incorporates the retention time (RT), mass-to-charge ratio (m/z), and isotopic distribution (derived from the sequence) of a peptide. Peptides with different modifications are considered different peptides. One assay will be generated for every combination of (modified) peptide sequence, charge state, and RT region that has been identified. The RT regions arise by pooling all identifications of the same peptide, considering a window of size @p extract:rt_window around every RT location that gave rise to an ID, and then merging overlapping windows.

   <B>2. Ion chromatogram extraction</B>

   Ion chromatograms (XICs) are extracted from the LC-MS data (parameter @p in). One XIC per isotope in an assay is generated, with the corresponding m/z value and RT range (variable, depending on the RT region of the assay).

   @see @ref TOPP_OpenSwathChromatogramExtractor

   <B>3. Feature detection</B>

   Next feature candidates - typically several per assay - are detected in the XICs and scored. A variety of scores for different quality aspects are calculated by OpenSWATH.

   @see @ref TOPP_OpenSwathAnalyzer

   <B>4. Feature classification</B>
   
   Feature candidates derived from assays with "internal" IDs are classed as "false positives" (candidates without matching internal IDs), "true positives" (the single best candidate per assay with matching internal IDs), and "ambiguous" (other candidates with matching internal IDs). If "external" IDs were given as input, features based on them are initially classed as "unknown". Also in this case, a support vector machine (SVM) is trained on the "true positive" and "false positive" candidates, to distinguish between the two classes based on the different OpenSWATH quality scores. After parameter optimization by cross-validation, the resulting SVM is used to classify the "unknown" feature candidates into (presumed) true or false positives.

   <B>5. Feature filtering</B>

   Feature candidates are filtered so that at most one feature per assay remains. For assays with internal IDs, only candidates previously classed as "true positive" are kept. For assays based solely on external IDs, out of the feature candidates classified as "true positive" by the SVM, the one with the highest SVM probability is kept.

   <B>6. Elution model fitting</B>

   Elution models can be fitted to the features to improve the quantification. For robustness, one model is fitted to all isotopic mass traces of a feature in parallel. A symmetric (Gaussian) and an asymmetric (exponential-Gaussian hybrid) model type are available. The fitted models are checked for plausibility before they are accepted.

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
  TOPPFeatureFinderIdentification() :
    TOPPBase("FeatureFinderIdentification", "Detects features in MS1 data based on peptide identifications.")
  {
    rt_term_.setCVIdentifierRef("MS");
    rt_term_.setAccession("MS:1000896");
    rt_term_.setName("normalized retention time");
    // available scores: initialPeakQuality,total_xic,peak_apices_sum,var_xcorr_coelution,var_xcorr_coelution_weighted,var_xcorr_shape,var_xcorr_shape_weighted,var_library_corr,var_library_rmsd,var_library_sangle,var_library_rootmeansquare,var_library_manhattan,var_library_dotprod,var_intensity_score,nr_peaks,sn_ratio,var_log_sn_score,var_elution_model_fit_score,xx_lda_prelim_score,var_isotope_correlation_score,var_isotope_overlap_score,var_massdev_score,var_massdev_score_weighted,var_bseries_score,var_yseries_score,var_dotprod_score,var_manhatt_score,main_var_xx_swath_prelim_score,xx_swath_prelim_score
    // exclude some redundant/uninformative scores:
    // @TODO: intensity bias introduced by "peak_apices_sum"?
    score_metavalues_ = "initialPeakQuality,peak_apices_sum,var_xcorr_coelution,var_xcorr_shape,var_library_sangle,var_intensity_score,sn_ratio,var_log_sn_score,var_elution_model_fit_score,xx_lda_prelim_score,var_isotope_correlation_score,var_isotope_overlap_score,var_massdev_score,main_var_xx_swath_prelim_score,xx_swath_prelim_score";
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file: LC-MS raw data");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("id", "<file>", "", "Input file: peptide identifications derived directly from 'in'");
    setValidFormats_("id", ListUtils::create<String>("idXML"));
    registerInputFile_("id_ext", "<file>", "", "Input file: 'external' peptide identifications (e.g. from aligned runs)", false);
    setValidFormats_("id_ext", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file: features");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerOutputFile_("lib_out", "<file>", "", "Output file: assay library", false);
    setValidFormats_("lib_out", ListUtils::create<String>("traML"));
    registerOutputFile_("chrom_out", "<file>", "", "Output file: chromatograms", false);
    setValidFormats_("chrom_out", ListUtils::create<String>("mzML"));
    registerOutputFile_("trafo_out", "<file>", "", "Output file: RT transformation", false);
    setValidFormats_("trafo_out", ListUtils::create<String>("trafoXML"));
    registerOutputFile_("candidates_out", "<file>", "", "Output file: feature candidates (before filtering and model fitting)", false);
    setValidFormats_("candidates_out", ListUtils::create<String>("featureXML"));

    registerTOPPSubsection_("extract", "Parameters for ion chromatogram extraction");
    StringList refs = ListUtils::create<String>("adapt,score,intensity,median,all");
    registerDoubleOption_("extract:rt_window", "<value>", 60.0, "RT window size (in sec.) for chromatogram extraction.", false);
    setMinFloat_("extract:rt_window", 0.0);
    registerDoubleOption_("extract:mz_window", "<value>", 10.0, "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)", false);
    setMinFloat_("extract:mz_window", 0.0);
    registerDoubleOption_("extract:isotope_pmin", "<value>", 0.03, "Minimum probability for an isotope to be included in the assay for a peptide.", false);
    setMinFloat_("extract:isotope_pmin", 0.0);
    setMaxFloat_("extract:isotope_pmin", 1.0);

    registerTOPPSubsection_("detect", "Parameters for detecting features in extracted ion chromatograms");
    registerDoubleOption_("detect:peak_width", "<value>", 60.0, "Expected elution peak width in seconds, for smoothing (Gauss filter)", false);
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
    registerFlag_("svm:unbiased", "Reduce biases by choosing (roughly) the same number of positive and negative observations, with the same intensity distribution, for training. This reduces the number of available samples.");
    registerOutputFile_("svm:xval_out", "<file>", "", "Output file: SVM cross-validation (parameter optimization) results", false);
    setValidFormats_("svm:xval_out", ListUtils::create<String>("csv"));
    Param svm_params;
    svm_params.insert("svm:", SimpleSVM().getParameters());
    registerFullParam_(svm_params);

    // parameters for model fitting (via ElutionModelFitter):
    registerTOPPSubsection_("model", "Parameters for fitting elution models to features");
    StringList models = ListUtils::create<String>("none,symmetric,asymmetric");
    registerStringOption_("model:type", "<choice>", models[0], "Type of elution model to fit to features", false);
    setValidStrings_("model:type", models);
    Param emf_params;
    emf_params.insert("model:", ElutionModelFitter().getParameters());
    emf_params.remove("model:asymmetric");
    registerFullParam_(emf_params);
  }


  typedef MSExperiment<Peak1D> PeakMap;
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTrace MassTrace;
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTraces MassTraces;

  // mapping: RT (not necessarily unique) -> pointer to peptide
  typedef multimap<double, PeptideIdentification*> RTMap;
  // mapping: charge -> internal/external: (RT -> pointer to peptide)
  typedef map<Int, pair<RTMap, RTMap> > ChargeMap;
  // mapping: sequence -> charge -> internal/external ID information
  typedef map<AASequence, ChargeMap> PeptideMap;

  // region in RT in which a peptide elutes:
  struct RTRegion
  {
    double start, end;
    set<Int> charges; // charge states for which there are IDs in the region
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

  PeakMap ms_data_; // input LC-MS data
  PeakMap chrom_data_; // accumulated chromatograms (XICs)
  bool keep_chromatograms_; // keep chromatogram data for output?
  TargetedExperiment library_; // accumulated assays for peptides
  bool keep_library_; // keep assay data for output?
  CVTerm rt_term_; // controlled vocabulary term for reference RT
  String score_metavalues_;
  TransformationDescription trafo_; // RT transformation (to range 0-1)
  String reference_rt_; // value of "reference_rt" parameter
  double rt_window_; // RT window width
  double mz_window_; // m/z window width
  bool mz_window_ppm_; // m/z window width is given in PPM (not Da)?
  double isotope_pmin_; // min. isotope probability
  double mapping_tolerance_; // RT tolerance for mapping IDs to features
  Size n_parts_; // number of partitions for SVM cross-validation
  Size n_samples_; // number of samples for SVM training
  ChromatogramExtractor extractor_; // OpenSWATH chromatogram extractor
  MRMFeatureFinderScoring feat_finder_; // OpenSWATH feature finder
  ProgressLogger prog_log_;


  // remove duplicate entries from a vector
  void removeDuplicateProteins_(TargetedExperiment& library)
  {
    // no "TargetedExperiment::Protein::operator<" defined, need to improvise:
    vector<TargetedExperiment::Protein> proteins;
    set<String> ids;
    for (vector<TargetedExperiment::Protein>::const_iterator it = 
           library.getProteins().begin(); it != library.getProteins().end();
         ++it)
    {
      if (!ids.count(it->id)) // new protein
      {
        proteins.push_back(*it);
        ids.insert(it->id);
      }
    }
    library.setProteins(proteins);
  }


  // add transitions for a peptide ion to the library:
  void generateTransitions_(const String& peptide_id, double mz, Int charge,
                            const IsotopeDistribution& iso_dist,
                            vector<ReactionMonitoringTransition>& transitions)
  {
    transitions.resize(iso_dist.size());
    // go through different isotopes:
    Size counter = 0;
    for (IsotopeDistribution::ConstIterator iso_it = iso_dist.begin();
         iso_it != iso_dist.end(); ++iso_it, ++counter)
    {
      String annotation = "i" + String(counter + 1);
      String transition_name = peptide_id + "_" + annotation;

      transitions[counter].setNativeID(transition_name);
      transitions[counter].setPrecursorMZ(mz);
      transitions[counter].setProductMZ(mz + Constants::C13C12_MASSDIFF_U * 
                                        float(counter) / charge);
      transitions[counter].setLibraryIntensity(iso_it->second);
      transitions[counter].setMetaValue("annotation", annotation);
      transitions[counter].setPeptideRef(peptide_id);
    }
  }


  void setPeptideRT_(TargetedExperiment::Peptide& peptide, double rt)
  {
    peptide.rts.clear();
    rt_term_.setValue(trafo_.apply(rt));
    TargetedExperiment::RetentionTime te_rt;
    te_rt.addCVTerm(rt_term_);
    peptide.rts.push_back(te_rt);
  }



  void getRTRegions_(const ChargeMap& peptide_data, 
                     vector<RTRegion>& rt_regions)
  {
    vector<pair<double, Int> > rts; // RTs and charges
    // use RTs from all charge states here to get a more complete picture:
    for (ChargeMap::const_iterator cm_it = peptide_data.begin();
         cm_it != peptide_data.end(); ++cm_it)
    {
      // "internal" IDs:
      for (RTMap::const_iterator rt_it = cm_it->second.first.begin();
           rt_it != cm_it->second.first.end(); ++rt_it)
      {
        rts.push_back(make_pair(rt_it->first, cm_it->first));
      }
      // "external" IDs:
      for (RTMap::const_iterator rt_it = cm_it->second.second.begin();
           rt_it != cm_it->second.second.end(); ++rt_it)
      {
        rts.push_back(make_pair(rt_it->first, cm_it->first));
      }
    }
    sort(rts.begin(), rts.end());
    double rt_tolerance = rt_window_ / 2.0;

    for (vector<pair<double, Int> >::iterator rt_it = rts.begin();
         rt_it != rts.end(); ++rt_it)
    {
      // create a new region?
      if (rt_regions.empty() ||
          (rt_regions.back().end < rt_it->first - rt_tolerance))
      {
        RTRegion region;
        region.start = rt_it->first - rt_tolerance;
        rt_regions.push_back(region);
      }
      rt_regions.back().end = rt_it->first + rt_tolerance;
      rt_regions.back().charges.insert(rt_it->second);
    }
  }


  void annotateFeatures_(FeatureMap& features,
                         const pair<RTMap, RTMap>& rt_data)
  {
    if (!rt_data.first.empty()) // validate based on internal IDs
    {
      // map IDs to features (based on RT):
      map<Size, vector<PeptideIdentification*> > feat_ids;
      for (Size i = 0; i < features.size(); ++i)
      {
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
        RTMap::const_iterator lower = rt_data.first.lower_bound(rt_min);
        RTMap::const_iterator upper = rt_data.first.upper_bound(rt_max);
        int id_count = 0;
        for (; lower != upper; ++lower)
        {
          feat_ids[i].push_back(lower->second);
          ++id_count;
        }
        features[i].setMetaValue("n_total_ids", rt_data.first.size());
        features[i].setMetaValue("n_matching_ids", id_count);
        if (id_count > 0) // matching IDs -> feature may be correct
        {
          features[i].setMetaValue("feature_class", "ambiguous");
        }
        else // no matching IDs -> feature is wrong
        {
          features[i].setMetaValue("feature_class", "false_positive");
        }
      }

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
              ((current_count == best_count) && // break ties by feature quality
               (features[current_index].getOverallQuality() >
                features[best_index].getOverallQuality())))
          {
            best_count = current_count;
            best_index = current_index;
          }
        }
        // assign IDs:
        if (best_count > 0)
        {
          features[best_index].getPeptideIdentifications().resize(best_count);
          for (Size i = 0; i < best_count; ++i)
          {
            features[best_index].getPeptideIdentifications()[i] = 
              *(feat_ids[best_index][i]);
            // we define the (one) feature with most matching IDs as correct:
            features[best_index].setMetaValue("feature_class", "true_positive");
            assigned_ids.insert(feat_ids[best_index][i]);
          }
        }
      }
      // store unassigned IDs:
      for (RTMap::const_iterator rt_it = rt_data.first.begin();
           rt_it != rt_data.first.end(); ++rt_it)
      {
        if (!assigned_ids.count(rt_it->second))
        {
          const PeptideIdentification& pep_id = *(rt_it->second);
          features.getUnassignedPeptideIdentifications().push_back(pep_id);
        }
      }
    }
    else // only external IDs -> no validation possible
    {
      for (FeatureMap::iterator feat_it = features.begin(); 
           feat_it != features.end(); ++feat_it)
      {
        feat_it->setMetaValue("n_total_ids", 0);
        feat_it->setMetaValue("n_matching_ids", -1);
        feat_it->setMetaValue("feature_class", "unknown");
        // add "dummy" peptide identification:
        PeptideIdentification id = *(rt_data.second.begin()->second);
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


  void detectFeaturesOnePeptide_(const PeptideMap::value_type& peptide_data,
                                 FeatureMap& features)
  {
    TargetedExperiment library;
    TargetedExperiment::Peptide peptide;

    const AASequence& seq = peptide_data.first;
    LOG_DEBUG << "\nPeptide: " << seq.toString() << endl;
    peptide.sequence = seq.toString();
    // @NOTE: Technically, "TargetedExperiment::Peptide" stores the unmodified
    // sequence and the modifications separately. Unfortunately, creating the
    // modifications vector is complex and there is currently no convenient
    // conversion function (see "TargetedExperimentHelper::getAASequence" for
    // the reverse conversion). However, "Peptide" is later converted to
    // "OpenSwath::LightPeptide" anyway, and this is done via "AASequence" (see
    // "OpenSwathDataAccessHelper::convertTargetedPeptide"). So for our purposes
    // it works to just store the sequence including modifications in "Peptide".

    // keep track of protein accessions:
    set<String> accessions;
    const pair<RTMap, RTMap>& pair = peptide_data.second.begin()->second;
    const PeptideHit& hit = (pair.first.empty() ?
                             pair.second.begin()->second->getHits()[0] :
                             pair.first.begin()->second->getHits()[0]);
    accessions = hit.extractProteinAccessions();
    // missing protein accession would crash OpenSWATH algorithms:
    if (accessions.empty()) accessions.insert("not_available");
    peptide.protein_refs = vector<String>(accessions.begin(), accessions.end());
    for (set<String>::iterator acc_it = accessions.begin(); 
         acc_it != accessions.end(); ++acc_it)
    {
      TargetedExperiment::Protein protein;
      protein.id = *acc_it;
      library.addProtein(protein);
    }

    // get isotope distribution for peptide:
    IsotopeDistribution iso_dist = 
      seq.getFormula(Residue::Full, 0).getIsotopeDistribution(10);
    iso_dist.trimLeft(isotope_pmin_);
    iso_dist.trimRight(isotope_pmin_);
    iso_dist.renormalize();

    // get regions in which peptide elutes (ideally only one):
    vector<RTRegion> rt_regions;
    getRTRegions_(peptide_data.second, rt_regions);
    LOG_DEBUG << "Found " << rt_regions.size() << " RT region(s)." << endl;

    // go through different charge states:
    for (ChargeMap::const_iterator cm_it = peptide_data.second.begin(); 
         cm_it != peptide_data.second.end(); ++cm_it)
    {
      Int charge = cm_it->first;
      double mz = seq.getMonoWeight(Residue::Full, charge) / charge;
      LOG_DEBUG << "Charge: " << charge << " (m/z: " << mz << ")" << endl;
      peptide.setChargeState(charge);
      peptide.id = peptide.sequence + "/" + String(charge);

      // we want to detect one feature per peptide, charge state and RT region 
      // (provided there is an ID for that charge in the region) - so there is 
      // always only one peptide in the library!
      Size counter = 0;
      for (vector<RTRegion>::iterator reg_it = rt_regions.begin();
           reg_it != rt_regions.end(); ++reg_it)
      {
        if (reg_it->charges.count(charge))
        {
          LOG_DEBUG << "Region " << counter + 1 << " (RT: "
                    << float(reg_it->start) << "-" << float(reg_it->end)
                    << ", size " << float(reg_it->end - reg_it->start) << ")"
                    << endl;

          vector<TargetedExperiment::Peptide> lib_peps(1, peptide);
          if (rt_regions.size() > 1)
          {
            lib_peps[0].id = peptide.id + ":" + String(++counter);
          }
          // use center of region as RT of assay (for chrom. extraction):
          double assay_rt = (reg_it->start + reg_it->end) / 2.0;
          setPeptideRT_(lib_peps[0], assay_rt);
          library.setPeptides(lib_peps);
          vector<ReactionMonitoringTransition> transitions;
          generateTransitions_(lib_peps[0].id, mz, charge, iso_dist, 
                               transitions);
          library.setTransitions(transitions);

          if (keep_library_) library_ += library;

          // extract chromatograms:
          double rt_window = reg_it->end - reg_it->start;
          PeakMap chrom_data;
          extractor_.extractChromatograms(ms_data_, chrom_data, library,
                                          mz_window_, mz_window_ppm_, trafo_,
                                          rt_window, "tophat");
          LOG_DEBUG << "Extracted " << chrom_data.getNrChromatograms()
                    << " chromatogram(s)." << endl;

          if (keep_chromatograms_)
          {
            Size n_chrom = (chrom_data.getNrChromatograms() + 
                            chrom_data_.getNrChromatograms());
            chrom_data_.reserveSpaceChromatograms(n_chrom);
            for (vector<MSChromatogram<> >::const_iterator ch_it = 
                   chrom_data.getChromatograms().begin(); ch_it !=
                   chrom_data.getChromatograms().end(); ++ch_it)
            {
              chrom_data_.addChromatogram(*ch_it);
            }
          }

          // find chromatographic peaks:
          FeatureMap current_features;
          Log_info.remove(cout); // suppress status output from OpenSWATH
          feat_finder_.pickExperiment(chrom_data, current_features, library,
                                      trafo_, ms_data_);
          Log_info.insert(cout);
          LOG_DEBUG << "Found " << current_features.size() << " feature(s)."
                    << endl;

          // complete feature annotation:
          for (FeatureMap::Iterator feat_it = current_features.begin();
               feat_it != current_features.end(); ++feat_it)
          {
            feat_it->setMZ(mz);
            feat_it->setCharge(charge);
            ensureConvexHulls_(*feat_it);
            // remove "fake" IDs generated by OpenSWATH (they would be removed
            // with a warning when writing output, because of missing protein
            // identification with corresponding identifier):
            feat_it->getPeptideIdentifications().clear();
            // annotate subordinates with theoretical isotope intensities:
            for (vector<Feature>::iterator sub_it =
                   feat_it->getSubordinates().begin(); sub_it !=
                   feat_it->getSubordinates().end(); ++sub_it)
            {
              String native_id = sub_it->getMetaValue("native_id");
              Size index = native_id.suffix('i').toInt() - 1;
              sub_it->setMetaValue("isotope_probability",
                                   iso_dist.getContainer()[index].second);
            }
          }
          // which features are supported by "internal" IDs?
          annotateFeatures_(current_features, cm_it->second);

          features += current_features;
        }
      }
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
                                          __PRETTY_FUNCTION__, msg);
    }
    if (n_neg < n_parts_)
    {
      String msg = "Not enough negative observations for " + 
        String(n_parts_) + "-fold cross-validation" + note + ".";
      throw Exception::MissingInformation(__FILE__, __LINE__, 
                                          __PRETTY_FUNCTION__, msg);
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
                                          __PRETTY_FUNCTION__, msg);
    }
    srand(time(0)); // seed random number generator
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
    // get predictors for SVM:
    vector<String> predictor_names = 
      ListUtils::create<String>(score_metavalues_);
    // values for all featues per predictor (this way around to simplify scaling
    // of predictors):
    map<String, vector<double> > predictors;
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
    bool unbiased = getFlag_("svm:unbiased");
    // mapping (for bias correction): intensity -> (index, positive?)
    multimap<double, pair<Size, bool> > valid_obs;
    Size n_obs[2] = {0, 0}; // counters for neg./pos. observations
    for (Size feat_index = 0; feat_index < features.size(); ++feat_index)
    {
      String feature_class = features[feat_index].getMetaValue("feature_class");
      Int label = -1;
      if (feature_class == "true_positive") label = 1;
      else if (feature_class == "false_positive") label = 0;
      if (label != -1)
      {
        ++n_obs[label];
        if (unbiased)
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

    if (unbiased) getUnbiasedSample_(valid_obs, training_labels);

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
    Param svm_params = getParam_().copy("svm:", true);
    svm_params.remove("unbiased");
    svm_params.remove("xval_out");
    svm.setParameters(svm_params);
    svm.setup(predictors, training_labels);
    String xval_out = getStringOption_("svm:xval_out");
    if (!xval_out.empty()) svm.writeXvalResults(xval_out);
    vector<SimpleSVM::Prediction> predictions;
    svm.predict(predictions);
    OPENMS_POSTCONDITION(predictions.size() == features.size(), 
                         "SVM predictions for all features expected");
    for (Size i = 0; i < features.size(); ++i)
    {
      features[i].setMetaValue("predicted_class", predictions[i].label);
      double prob_positive = predictions[i].probabilities[1];
      features[i].setMetaValue("predicted_probability", prob_positive);
      features[i].setOverallQuality(prob_positive);
    }
  }


  void filterFeatures_(FeatureMap& features, bool classified)
  {
    if (features.empty()) return;
    
    if (classified)
    {
      // Remove features with class "false_pos." or "ambiguous", keep
      // "true_pos."; for class "unknown", for every assay (meta value
      // "PeptideRef"), keep the feature with "predicted_class" 1 and highest
      // "predicted_probability" (= overall quality). We mark features for
      // removal by setting their overall quality to zero.
      FeatureMap::Iterator best_it = features.begin();
      double best_quality = 0.0;
      String previous_ref = features[0].getMetaValue("PeptideRef");
      for (FeatureMap::Iterator it = features.begin(); it != features.end();
           ++it)
      {
        // features from same assay (same "PeptideRef") appear consecutively;
        // if this is a new assay, finalize the previous one:
        const String& peptide_ref = it->getMetaValue("PeptideRef");
        if (peptide_ref != previous_ref)
        {
          if (best_quality > 0.0) best_it->setOverallQuality(best_quality);
          best_quality = 0.0;
          previous_ref = peptide_ref;
        }

        // update qualities:
        const String& feature_class = it->getMetaValue("feature_class");
        if ((feature_class == "unknown") &&
            (Int(it->getMetaValue("predicted_class")) > 0) &&
            (it->getOverallQuality() > best_quality))
        {
          best_it = it;
          best_quality = it->getOverallQuality();
        }
        if (feature_class != "true_positive") it->setOverallQuality(0.0);
      }
      // set of features from the last assay:
      if (best_quality > 0.0)
      {
        best_it->setOverallQuality(best_quality);
      }

      features.erase(remove_if(features.begin(), features.end(),
                               feature_filter_quality_), features.end());
    }
    else
    {
      features.erase(remove_if(features.begin(), features.end(),
                               feature_filter_peptides_), features.end());
    }
  }
  

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String id = getStringOption_("id");
    String id_ext = getStringOption_("id_ext");
    String out = getStringOption_("out");
    String lib_out = getStringOption_("lib_out");
    String chrom_out = getStringOption_("chrom_out");
    String trafo_out = getStringOption_("trafo_out");
    String candidates_out = getStringOption_("candidates_out");
    rt_window_ = getDoubleOption_("extract:rt_window");
    mz_window_ = getDoubleOption_("extract:mz_window");
    mz_window_ppm_ = mz_window_ >= 1;
    isotope_pmin_ = getDoubleOption_("extract:isotope_pmin");
    double peak_width = getDoubleOption_("detect:peak_width");
    double min_peak_width = getDoubleOption_("detect:min_peak_width");
    double signal_to_noise = getDoubleOption_("detect:signal_to_noise");
    mapping_tolerance_ = getDoubleOption_("detect:mapping_tolerance");
    n_parts_ = getIntOption_("svm:xval");
    n_samples_ = getIntOption_("svm:samples");
    String elution_model = getStringOption_("model:type");
    prog_log_.setLogType(log_type_);

    if ((n_samples_ > 0) && (n_samples_ < 2 * n_parts_))
    {
      String msg = "Sample size of " + String(n_samples_) +
        " (parameter 'svm:samples') is not enough for " + String(n_parts_) +
        "-fold cross-validation (parameter 'svm:xval').";
      throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                        msg);
    }
    
    //-------------------------------------------------------------
    // load input
    //-------------------------------------------------------------
    LOG_INFO << "Loading input data..." << endl;
    MzMLFile mzml;
    mzml.setLogType(log_type_);
    mzml.getOptions().addMSLevel(1);
    mzml.load(in, ms_data_);
    if (reference_rt_ == "intensity") ms_data_.sortSpectra(true);

    // RT transformation to range 0-1:
    ms_data_.updateRanges();
    double min_rt = ms_data_.getMinRT(), max_rt = ms_data_.getMaxRT();
    TransformationDescription::DataPoints points;
    points.push_back(make_pair(min_rt, 0.0));
    points.push_back(make_pair(max_rt, 1.0));
    trafo_.setDataPoints(points);
    trafo_.fitModel("linear");
    if (!trafo_out.empty())
    {
      TransformationXMLFile().store(trafo_out, trafo_);
    }

    // initialize algorithm classes needed later:
    extractor_.setLogType(ProgressLogger::NONE);
    Param params = feat_finder_.getParameters();
    params.setValue("stop_report_after_feature", -1); // return all features
    params.setValue("Scores:use_rt_score", "false"); // RT may not be reliable
    if (elution_model != "none") params.setValue("write_convex_hull", "true");
    if (min_peak_width < 1.0) min_peak_width *= peak_width;
    params.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width",
                    peak_width);
    params.setValue("TransitionGroupPicker:min_peak_width", min_peak_width);
    // disabling the signal-to-noise threshold (setting the parameter to zero)
    // totally breaks the OpenSWATH feature detection (no features found)!
    params.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise",
                    signal_to_noise);
    params.setValue("TransitionGroupPicker:recalculate_peaks", "true");
    params.setValue("TransitionGroupPicker:compute_peak_quality", "true");
    params.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
    params.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
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
    if (!id_ext.empty()) IdXMLFile().load(id_ext, proteins_ext, peptides_ext);

    //-------------------------------------------------------------
    // prepare peptide map
    //-------------------------------------------------------------
    LOG_INFO << "Preparing mapping of peptide data..." << endl;
    PeptideMap peptide_map;
    for (vector<PeptideIdentification>::iterator pep_it = peptides.begin();
         pep_it != peptides.end(); ++pep_it)
    {
      addPeptideToMap_(*pep_it, peptide_map);
      pep_it->setMetaValue("FFId_category", "internal");
    }
    Size n_internal_ids = peptide_map.size();
    for (vector<PeptideIdentification>::iterator pep_it = peptides_ext.begin();
         pep_it != peptides_ext.end(); ++pep_it)
    {
      addPeptideToMap_(*pep_it, peptide_map, true);
      pep_it->setMetaValue("FFId_category", "external");
    }
    Size n_external_ids = peptide_map.size() - n_internal_ids;

    //-------------------------------------------------------------
    // run feature detection
    //-------------------------------------------------------------
    LOG_INFO << "Running feature detection..." << endl;
    FeatureMap features;
    keep_library_ = !lib_out.empty();
    keep_chromatograms_ = !chrom_out.empty();
    Size prog_counter = 0;
    prog_log_.startProgress(1, peptide_map.size(), "running feature detection");
    for (PeptideMap::iterator pm_it = peptide_map.begin();
         pm_it != peptide_map.end(); ++pm_it)
    {
      FeatureMap current_features;
      detectFeaturesOnePeptide_(*pm_it, current_features);
      features += current_features;
      prog_log_.setProgress(++prog_counter);
    }
    prog_log_.endProgress();
    LOG_DEBUG << "Found " << features.size() << " feature candidates in total."
              << endl;
    ms_data_.reset(); // not needed anymore, free up the memory

    // write auxiliary output:
    if (keep_library_)
    {
      removeDuplicateProteins_(library_);
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

    // don't do SVM stuff unless we have external data to apply the model to:
    if (!id_ext.empty()) classifyFeatures_(features);

    if (!candidates_out.empty()) // store feature candidates
    {
      FeatureXMLFile().store(candidates_out, features);
    }

    filterFeatures_(features, !id_ext.empty());
    LOG_INFO << features.size() << " features left after filtering." << endl;
    
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
    LOG_INFO << "\nSummary statistics (counting distinct peptides including "
      "PTMs):\n"
             << peptide_map.size() << " peptides identified ("
             << n_internal_ids << " internal, " << n_external_ids
             << " additional external)\n"
             << quantified_all.size() << " peptides with features ("
             << quantified_internal.size() << " internal, "
             << n_quant_external << " additional external)\n"
             << peptide_map.size() - quantified_all.size()
             << " peptides without features ("
             << n_internal_ids - quantified_internal.size() << " internal, "
             << n_external_ids - n_quant_external << " additional external)\n"
             << endl;
      
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
