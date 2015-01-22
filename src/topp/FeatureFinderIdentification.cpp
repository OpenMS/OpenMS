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
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
// #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
// #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>

#include <svm.h>

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
         <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
         <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_ProteinQuantifier</td>
       </tr>
       <tr>
         <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
       </tr>
     </table>
   </CENTER>

   This tool detects quantitative features in MS1 data based on information from peptide identifications (derived from MS2 spectra). It uses algorithms for targeted data analysis from the OpenSWATH pipeline.

   @note It is important that only high-confidence peptide identifications and centroided (peak-picked) LC-MS data are used as inputs!

   For every distinct peptide ion (defined by sequence and charge) in the input (parameter @p id), an assay is generated, incorporating the retention time (RT), mass-to-charge ratio (m/z), and isotopic distribution of the peptide. The parameter @p reference_rt controls how the RT of the assay is determined if the peptide has been observed multiple times. The relative intensities of the isotopes together with their m/z values are calculated from the sequence and charge.

   The assays are used to perform targeted data analysis on the MS1 level using OpenSWATH algorithms, in several steps:

   <B>1. Ion chromatogram extraction</B>

   First ion chromatograms (XICs) are extracted from the data (parameter @p in). For every assay, the RT range of the XICs is given by @p extract:rt_window (around the reference RT of the assay) and the m/z ranges by @p extract:mz_window (around the m/z values of all included isotopes). As an exception to this, if @p extract:reference_rt is @p adapt, a more complex procedure is used to find the RT range: A range of size @p rt_window around every relevant peptide ID is considered, overlapping ranges are joined, and the largest resulting range is used for the extraction. In that case, the reference RT for the assay is the median RT of the peptide IDs within the range.

   @see @ref TOPP_OpenSwathChromatogramExtractor

   <B>2. Feature detection</B>

   Next feature candidates are detected in the XICs and scored. The best candidate per assay according to the OpenSWATH scoring is turned into a feature.

   @see @ref TOPP_OpenSwathAnalyzer

   <B>3. Elution model fitting</B>

   Elution models can be fitted to every feature to improve the quantification. For robustness, one model is fitted to all isotopic mass traces of a feature in parallel. A symmetric (Gaussian) and an asymmetric (exponential-Gaussian hybrid) model type are available. The fitted models are checked for plausibility before they are accepted.

   Finally the results (feature maps, parameter @p out) are returned.

   @note This tool aims to report a feature for every distinct peptide ion given in the @p id input. Currently no attempt is made to filter out false-positives (although this may be possible in post-processing based on the OpenSWATH scores). If only high-confidence peptide IDs are used, that come from the same LC-MS/MS run that is being quantified, this should not be a problem. However, if e.g. inferred IDs from different runs (see @ref TOPP_MapAlignerIdentification) are included, false-positive features with arbitrary intensities may result for peptides that cannot be detected in the present data.

   <B>The command line parameters of this tool are:</B>
   @verbinclude TOPP_FeatureFinderIdentification.cli
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
    score_metavalues_ = "initialPeakQuality,total_xic,peak_apices_sum,var_xcorr_coelutio,var_xcorr_coelution_weighted,var_xcorr_shape,var_xcorr_shape_weighted,var_library_corr,var_library_rmsd,var_library_sangle,var_library_rootmeansquare,var_library_manhattan,var_library_dotprod,var_intensity_score,nr_peaks,sn_ratio,var_log_sn_score,var_elution_model_fit_score,xx_lda_prelim_score,var_isotope_correlation_score,var_isotope_overlap_score,var_massdev_score,var_massdev_score_weighted,var_bseries_score,var_yseries_score,var_dotprod_score,var_manhatt_score,main_var_xx_swath_prelim_score,xx_swath_prelim_score";
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
    registerOutputFile_("lib_out", "<file>", "", "Output file: assay library ('internal' IDs)", false);
    setValidFormats_("lib_out", ListUtils::create<String>("traML"));
    registerOutputFile_("chrom_out", "<file>", "", "Output file: chromatograms ('internal' IDs)", false);
    setValidFormats_("chrom_out", ListUtils::create<String>("mzML"));
    registerOutputFile_("trafo_out", "<file>", "", "Output file: RT transformation", false);
    setValidFormats_("trafo_out", ListUtils::create<String>("trafoXML"));

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
    registerDoubleOption_("detect:signal_to_noise", "<value>", 0.5, "Signal-to-noise threshold for OpenSWATH feature detection", false, true);
    setMinFloat_("detect:signal_to_noise", 0.1);
    registerDoubleOption_("detect:mapping_tolerance", "<value>", 10.0, "RT tolerance (plus/minus) for mapping peptide IDs to features. Absolute value in seconds if 1 or greater, else relative to the RT span of the feature.", false);
    setMinFloat_("detect:mapping_tolerance", 0.0);

    registerTOPPSubsection_("svm", "Parameters for scoring features using a support vector machine (SVM)");
    registerIntOption_("svm:xval", "<number>", 5, "Number of partitions for cross-validation", false);
    setMinInt_("svm:xval", 1);
    registerStringOption_("svm:kernel", "<choice>", "RBF", "SVM kernel", false);
    setValidStrings_("svm:kernel", ListUtils::create<String>("RBF,linear"));
    String values = "-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15";
    registerDoubleList_("svm:log2_C", "<values>", ListUtils::create<double>(values), "Values to try for the SVM parameter 'C' during parameter optimization. A value 'x' is used as 'C = 2^x'.", false);
    values = "-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3";
    registerDoubleList_("svm:log2_gamma", "<values>", ListUtils::create<double>(values), "Values to try for the SVM parameter 'gamma' during parameter optimization (RBF kernel only). A value 'x' is used as 'gamma = 2^x'.", false);

    registerTOPPSubsection_("model", "Parameters for fitting elution models to features");
    StringList models = ListUtils::create<String>("symmetric,asymmetric,none");
    registerStringOption_("model:type", "<choice>", models[0], "Type of elution model to fit to features", false);
    setValidStrings_("model:type", models);
    registerDoubleOption_("model:add_zeros", "<value>", 0.2, "Add zero-intensity points outside the feature range to constrain the model fit. This parameter sets the weight given to these points during model fitting; '0' to disable.", false, true);
    setMinFloat_("model:add_zeros", 0.0);
    registerFlag_("model:unweighted_fit", "Suppress weighting of mass traces according to theoretical intensities when fitting elution models", true);
    registerFlag_("model:no_imputation", "If fitting the elution model fails for a feature, set its intensity to zero instead of imputing a value from the initial intensity estimate", true);
    registerTOPPSubsection_("model:check", "Parameters for checking the validity of elution models (and rejecting them if necessary)");
    registerDoubleOption_("model:check:boundaries", "<value>", 0.5, "Time points corresponding to this fraction of the elution model height have to be within the data region used for model fitting", false, true);
    setMinFloat_("model:check:boundaries", 0.0);
    setMaxFloat_("model:check:boundaries", 1.0);
    registerDoubleOption_("model:check:width", "<value>", 10.0, "Upper limit for acceptable widths of elution models (Gaussian or EGH), expressed in terms of modified (median-based) z-scores; '0' to disable", false, true);
    setMinFloat_("model:check:width", 0.0);
    registerDoubleOption_("model:check:asymmetry", "<value>", 10.0, "Upper limit for acceptable asymmetry of elution models (EGH only), expressed in terms of modified (median-based) z-scores; '0' to disable", false, true);
    setMinFloat_("model:check:asymmetry", 0.0);
  }

  typedef MSExperiment<Peak1D> PeakMap;
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTrace MassTrace;
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTraces MassTraces;

  // mapping: RT (not necessarily unique) -> (pointer to peptide, is external?)
  typedef multimap<double, pair<PeptideIdentification*, bool> > RTMap;
  // mapping: charge -> RT -> pointer to peptide
  typedef map<Int, RTMap> ChargeMap;
  // mapping: sequence -> charge -> RT -> pointer to peptide
  typedef map<AASequence, ChargeMap> PeptideMap;

  // region in RT in which a peptide elutes:
  struct RTRegion
  {
    double start, end;
    map<Int, vector<double> > evidence; // mapping: charge -> RT of ID
  };

  // mapping: (SVM param. C, SVM param. gamma) -> classification performance
  typedef map<pair<double, double>, double> SVMPerformance;

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
  String elution_model_; // choice of elution model
  ChromatogramExtractor extractor_; // OpenSWATH chromatogram extractor
  MRMFeatureFinderScoring feat_finder_; // OpenSWATH feature finder


  // like "median", but returns the middle-right value for an even number of
  // values (no averaging of the middle two):
  double getMedoid_(vector<double>& sorted_values)
  {
    if (sorted_values.size() == 1) return sorted_values[0]; // common case
    vector<double>::iterator start = sorted_values.begin();
    if (sorted_values.size() % 2 == 0) ++start;
    return Math::median(start, sorted_values.end(), true);
  }


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

/*
  double calculateFitQuality_(const TraceFitter* fitter, 
                              const MassTraces& traces)
  {
    double mre = 0.0;
    double total_weights = 0.0;
    double rt_start = max(fitter->getLowerRTBound(), traces[0].peaks[0].first);
    double rt_end = min(fitter->getUpperRTBound(), 
                        traces[0].peaks.back().first);

    for (MassTraces::const_iterator tr_it = traces.begin();
         tr_it != traces.end(); ++tr_it)
    {
      for (vector<pair<double, const Peak1D*> >::const_iterator p_it = 
             tr_it->peaks.begin(); p_it != tr_it->peaks.end(); ++p_it)
      {
        double rt = p_it->first;
        if ((rt >= rt_start) && (rt <= rt_end))
        {
          double model_value = fitter->getValue(rt);
          double diff = fabs(model_value * tr_it->theoretical_int -
                             p_it->second->getIntensity());
          mre += diff / model_value;
          total_weights += tr_it->theoretical_int;
        }
      }
    }
    return mre / total_weights;
  }


  // fit models of elution profiles to all features:
  void fitElutionModels_(FeatureMap& features)
  {
    // assumptions:
    // - all features have subordinates (for the mass traces/transitions)
    // - all subordinates have one convex hull
    // - all convex hulls in one feature contain the same number (> 0) of points
    // - the y coordinates of the hull points store the intensities

    bool asymmetric = (elution_model_ == "asymmetric");
    double add_zeros = getDoubleOption_("model:add_zeros");
    bool weighted = !getFlag_("model:unweighted_fit");
    bool impute = !getFlag_("model:no_imputation");
    double check_boundaries = getDoubleOption_("model:check:boundaries");

    // prepare look-up of transitions by native ID:
    map<String, const ReactionMonitoringTransition*> trans_ids;
    for (vector<ReactionMonitoringTransition>::const_iterator trans_it =
           library_.getTransitions().begin(); trans_it !=
           library_.getTransitions().end(); ++trans_it)
    {
      trans_ids[trans_it->getNativeID()] = &(*trans_it);
    }

    TraceFitter* fitter;
    if (asymmetric)
    {
      fitter = new EGHTraceFitter();
    }
    else fitter = new GaussTraceFitter();
    if (weighted)
    {
      Param params = fitter->getDefaults();
      params.setValue("weighted", "true");
      fitter->setParameters(params);
    }

    // store model parameters to find outliers later:
    double width_limit = getDoubleOption_("model:check:width");
    double asym_limit = (asymmetric ?
                         getDoubleOption_("model:check:asymmetry") : 0.0);
    // store values redundantly - once aligned with the features in the map,
    // once only for successful models:
    vector<double> widths_all, widths_good, asym_all, asym_good;
    if (width_limit > 0)
    {
      widths_all.resize(features.size(),
                        numeric_limits<double>::quiet_NaN());
      widths_good.reserve(features.size());
    }
    if (asym_limit > 0)
    {
      asym_all.resize(features.size(), numeric_limits<double>::quiet_NaN());
      asym_good.reserve(features.size());
    }

    // collect peaks that constitute mass traces:
    LOG_DEBUG << "Fitting elution models to features:" << endl;
    Size index = 0;
    for (FeatureMap::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it, ++index)
    {
      LOG_DEBUG << String(feat_it->getMetaValue("PeptideRef")) << endl;
      double region_start = double(feat_it->getMetaValue("leftWidth"));
      double region_end = double(feat_it->getMetaValue("rightWidth"));

      vector<Peak1D> peaks;
      // reserve space once, to avoid copying and invalidating pointers:
      const Feature& sub = feat_it->getSubordinates()[0];
      Size points_per_hull = sub.getConvexHulls()[0].getHullPoints().size();
      peaks.reserve(feat_it->getSubordinates().size() * points_per_hull +
                    (add_zeros > 0.0)); // don't forget additional zero point
      MassTraces traces;
      traces.max_trace = 0;
      // need a mass trace for every transition, plus maybe one for add. zeros:
      traces.reserve(feat_it->getSubordinates().size() + (add_zeros > 0.0));
      for (vector<Feature>::iterator sub_it =
             feat_it->getSubordinates().begin(); sub_it !=
             feat_it->getSubordinates().end(); ++sub_it)
      {
        MassTrace trace;
        trace.peaks.reserve(points_per_hull);
        String native_id = sub_it->getMetaValue("native_id");
        trace.theoretical_int = trans_ids[native_id]->getLibraryIntensity();
        sub_it->setMetaValue("isotope_probability", trace.theoretical_int);
        const ConvexHull2D& hull = sub_it->getConvexHulls()[0];
        for (ConvexHull2D::PointArrayTypeConstIterator point_it =
               hull.getHullPoints().begin(); point_it !=
               hull.getHullPoints().end(); ++point_it)
        {
          double intensity = point_it->getY();
          if (intensity > 0.0) // only use non-zero intensities for fitting
          {
            Peak1D peak;
            peak.setMZ(sub_it->getMZ());
            peak.setIntensity(intensity);
            peaks.push_back(peak);
            trace.peaks.push_back(make_pair(point_it->getX(), &peaks.back()));
          }
        }
        trace.updateMaximum();
        if (!trace.peaks.empty()) traces.push_back(trace);
      }

      // find the trace with maximal intensity:
      Size max_trace = 0;
      double max_intensity = 0;
      for (Size i = 0; i < traces.size(); ++i)
      {
        if (traces[i].max_peak->getIntensity() > max_intensity)
        {
          max_trace = i;
          max_intensity = traces[i].max_peak->getIntensity();
        }
      }
      traces.max_trace = max_trace;
      // ??? (from "FeatureFinderAlgorithmPicked.h"):
      // traces.updateBaseline();
      // traces.baseline = 0.75 * traces.baseline;
      traces.baseline = 0.0;

      if (add_zeros > 0.0)
      {
        MassTrace trace;
        trace.peaks.reserve(2);
        trace.theoretical_int = add_zeros;
        Peak1D peak;
        peak.setMZ(feat_it->getSubordinates()[0].getMZ());
        peak.setIntensity(0.0);
        peaks.push_back(peak);
        double offset = 0.2 * (region_start - region_end);
        trace.peaks.push_back(make_pair(region_start - offset, &peaks.back()));
        trace.peaks.push_back(make_pair(region_end + offset, &peaks.back()));
        traces.push_back(trace);
      }

      // fit the model:
      bool fit_success = true;
      try
      {
        fitter->fit(traces);
      }
      catch (Exception::UnableToFit& except)
      {
        LOG_ERROR << "Error fitting model to feature '"
                  << feat_it->getUniqueId() << "': " << except.getName()
                  << " - " << except.getMessage() << endl;
        fit_success = false;
      }

      // record model parameters:
      double center = fitter->getCenter(), height = fitter->getHeight();
      feat_it->setMetaValue("model_height", height);
      feat_it->setMetaValue("model_FWHM", fitter->getFWHM());
      feat_it->setMetaValue("model_center", center);
      feat_it->setMetaValue("model_lower", fitter->getLowerRTBound());
      feat_it->setMetaValue("model_upper", fitter->getUpperRTBound());
      if (asymmetric)
      {
        EGHTraceFitter* egh =
          static_cast<EGHTraceFitter*>(fitter);
        feat_it->setMetaValue("model_EGH_tau", egh->getTau());
        feat_it->setMetaValue("model_EGH_sigma", egh->getSigma());
      }
      else
      {
        GaussTraceFitter* gauss =
          static_cast<GaussTraceFitter*>(fitter);
        feat_it->setMetaValue("model_Gauss_sigma", gauss->getSigma());
      }

      // goodness of fit:
      double mre = -1.0; // mean relative error
      if (fit_success)
      {
        mre = calculateFitQuality_(fitter, traces);
      }
      feat_it->setMetaValue("model_error", mre);

      // check model validity:
      double area = fitter->getArea();
      feat_it->setMetaValue("model_area", area);
      if ((area != area) || (area <= 0.0)) // x != x: test for NaN
      {
        feat_it->setMetaValue("model_status", "1 (invalid area)");
      }
      else if ((center <= region_start) || (center >= region_end))
      {
        feat_it->setMetaValue("model_status", "2 (center out of bounds)");
      }
      else if (fitter->getValue(region_start) > check_boundaries * height)
      {
        feat_it->setMetaValue("model_status", "3 (left side out of bounds)");
      }
      else if (fitter->getValue(region_end) > check_boundaries * height)
      {
        feat_it->setMetaValue("model_status", "4 (right side out of bounds)");
      }
      else
      {
        feat_it->setMetaValue("model_status", "0 (valid)");
        // store model parameters to find outliers later:
        if (asymmetric)
        {
          double sigma = feat_it->getMetaValue("model_EGH_sigma");
          double abs_tau = fabs(double(feat_it->
                                       getMetaValue("model_EGH_tau")));
          if (width_limit > 0)
          {
            // see implementation of "EGHTraceFitter::getArea":
            double width = sigma * 0.6266571 + abs_tau;
            widths_all[index] = width;
            widths_good.push_back(width);
          }
          if (asym_limit > 0)
          {
            double asymmetry = abs_tau / sigma;
            asym_all[index] = asymmetry;
            asym_good.push_back(asymmetry);
          }
        }
        else if (width_limit > 0)
        {
          double width = feat_it->getMetaValue("model_Gauss_sigma");
          widths_all[index] = width;
          widths_good.push_back(width);
        }
      }
    }
    delete fitter;

    // find outliers in model parameters:
    if (width_limit > 0)
    {
      double median_width = Math::median(widths_good.begin(),
                                         widths_good.end());
      vector<double> abs_diffs(widths_good.size());
      for (Size i = 0; i < widths_good.size(); ++i)
      {
        abs_diffs[i] = fabs(widths_good[i] - median_width);
      }
      // median absolute deviation (constant factor to approximate std. dev.):
      double mad_width = 1.4826 * Math::median(abs_diffs.begin(),
                                               abs_diffs.end());

      for (Size i = 0; i < features.size(); ++i)
      {
        double width = widths_all[i];
        if (width != width) continue; // NaN (failed model)
        double z_width = (width - median_width) / mad_width; // mod. z-score
        if (z_width > width_limit)
        {
          features[i].setMetaValue("model_status", "5 (width too large)");
          if (asym_limit > 0) // skip asymmetry check below
          {
            asym_all[i] = numeric_limits<double>::quiet_NaN();
          }
        }
        else if (z_width < -width_limit)
        {
          features[i].setMetaValue("model_status", "6 (width too small)");
          if (asym_limit > 0) // skip asymmetry check below
          {
            asym_all[i] = numeric_limits<double>::quiet_NaN();
          }
        }
      }
    }
    if (asym_limit > 0)
    {
      double median_asym = Math::median(asym_good.begin(), asym_good.end());
      vector<double> abs_diffs(asym_good.size());
      for (Size i = 0; i < asym_good.size(); ++i)
      {
        abs_diffs[i] = fabs(asym_good[i] - median_asym);
      }
      // median absolute deviation (constant factor to approximate std. dev.):
      double mad_asym = 1.4826 * Math::median(abs_diffs.begin(),
                                              abs_diffs.end());

      for (Size i = 0; i < features.size(); ++i)
      {
        double asym = asym_all[i];
        if (asym != asym) continue; // NaN (failed model)
        double z_asym = (asym - median_asym) / mad_asym; // mod. z-score
        if (z_asym > asym_limit)
        {
          features[i].setMetaValue("model_status", "7 (asymmetry too high)");
        }
        else if (z_asym < -asym_limit) // probably shouldn't happen in practice
        {
          features[i].setMetaValue("model_status", "8 (asymmetry too low)");
        }
      }
    }

    // impute approximate results for failed model fits:
    TransformationModel::DataPoints quant_values;
    vector<FeatureMap::Iterator> failed_models;
    Size model_successes = 0, model_failures = 0;

    for (FeatureMap::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it, ++index)
    {
      feat_it->setMetaValue("raw_intensity", feat_it->getIntensity());
      if (String(feat_it->getMetaValue("model_status"))[0] != '0')
      {
        if (impute) failed_models.push_back(feat_it);
        else feat_it->setIntensity(0.0);
        model_failures++;
      }
      else
      {
        double area = feat_it->getMetaValue("model_area");
        if (impute)
        { // apply log-transform to weight down high outliers:
          double raw_intensity = feat_it->getIntensity();
          LOG_DEBUG << "Successful model: x = " << raw_intensity << ", y = "
                    << area << "; log(x) = " << log(raw_intensity)
                    << ", log(y) = " << log(area) << endl;
          quant_values.push_back(make_pair(log(raw_intensity), log(area)));
        }
        feat_it->setIntensity(area);
        model_successes++;
      }
    }
    LOG_INFO << "Model fitting: " << model_successes << " successes, "
             << model_failures << " failures" << endl;

    if (impute) // impute results for cases where the model fit failed
    {
      TransformationModelLinear lm(quant_values, Param());
      double slope, intercept;
      lm.getParameters(slope, intercept);
      LOG_DEBUG << "LM slope: " << slope << ", intercept: " << intercept
                << endl;
      for (vector<FeatureMap::Iterator>::iterator it = failed_models.begin();
           it != failed_models.end(); ++it)
      {
        double area = exp(lm.evaluate(log((*it)->getIntensity())));
        (*it)->setIntensity(area);
      }
    }
  }
*/

  void getRTRegions_(const ChargeMap& peptide_data, 
                     vector<RTRegion>& rt_regions)
  {
    vector<pair<double, Int> > rts; // RTs and charges
    // use RTs from all charge states here to get a more complete picture:
    for (ChargeMap::const_iterator cm_it = peptide_data.begin();
         cm_it != peptide_data.end(); ++cm_it)
    {
      for (RTMap::const_iterator rt_it = cm_it->second.begin();
           rt_it != cm_it->second.end(); ++rt_it)
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
      rt_regions.back().evidence[rt_it->second].push_back(rt_it->first);
    }
  }


  void annotateFeatures_(FeatureMap& features, const RTMap& rt_data)
  {
    Size n_internal_ids = 0;
    for (RTMap::const_iterator rt_it = rt_data.begin(); rt_it != rt_data.end();
         ++rt_it)
    {
      if (!rt_it->second.second) n_internal_ids++;
      // @TODO: should we filter out the external IDs?
    }
    if (n_internal_ids > 0) // validate based on internal IDs
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
        RTMap::const_iterator lower = rt_data.lower_bound(rt_min);
        RTMap::const_iterator upper = rt_data.upper_bound(rt_max);
        int id_count = 0;
        for (; lower != upper; ++lower)
        {
          // only use "internal" IDs:
          if (!lower->second.second)
          {
            feat_ids[i].push_back(lower->second.first);
            ++id_count;
          }
        }
        features[i].setMetaValue("n_total_ids", n_internal_ids);
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
      for (RTMap::const_iterator rt_it = rt_data.begin();
           rt_it != rt_data.end(); ++rt_it)
      {
        if (!assigned_ids.count(rt_it->second.first))
        {
          const PeptideIdentification& pep_id = *(rt_it->second.first);
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

    // keep track of protein accessions:
    const PeptideHit& hit = 
      peptide_data.second.begin()->second.begin()->second.first->getHits()[0];
    set<String> accessions = hit.extractProteinAccessions();
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
      // (provided there is an ID for that charge in the region) - there is 
      // always only one peptide in the library!
      Size counter = 0;
      for (vector<RTRegion>::iterator reg_it = rt_regions.begin();
           reg_it != rt_regions.end(); ++reg_it)
      {
        if (reg_it->evidence.count(charge))
        {
          LOG_DEBUG << "Region " << counter + 1 << " (RT: "
                    << float(reg_it->start) << "-" << float(reg_it->end) << ")"
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
          feat_finder_.pickExperiment(chrom_data, current_features, library,
                                      trafo_, ms_data_);
          LOG_DEBUG << "Found " << current_features.size() << " feature(s)."
                    << endl;

          // complete feature annotation:
          for (FeatureMap::Iterator feat_it = current_features.begin();
               feat_it != current_features.end(); ++feat_it)
          {
            feat_it->setMZ(mz);
            feat_it->setCharge(charge);
            ensureConvexHulls_(*feat_it);
            // remove "fake" IDs generated by OpenSWATH:
            feat_it->getPeptideIdentifications().clear();
          }
          // which features are supported by IDs?
          const RTMap& rt_data = cm_it->second;
          annotateFeatures_(current_features, rt_data);

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
    pair<PeptideIdentification*, bool> pep_pair = make_pair(&peptide, external);
    peptide_map[hit.getSequence()][charge].insert(make_pair(rt, pep_pair));
  }


  void scaleSVMData_(vector<vector<double> >& predictor_values,
                     const vector<String>& predictors)
  {
    for (Size i = 0; i < predictors.size(); ++i)
    {
      vector<double>::iterator pv_begin = predictor_values[i].begin();
      vector<double>::iterator pv_end = predictor_values[i].end();
      double vmin = *min_element(pv_begin, pv_end);
      double vmax = *max_element(pv_begin, pv_end);
      if (vmin == vmax)
      {
        LOG_DEBUG << "Predictor '" + predictors[i] + "' is uninformative." 
                  << endl;
        predictor_values[i].clear();
        continue;
      }
      double range = vmax - vmin;
      for (vector<double>::iterator it = pv_begin; it != pv_end; ++it)
      {
        *it = (*it - vmin) / range;
      }
    }
  }


  void convertSVMData_(vector<vector<double> >& predictor_values,
                       vector<vector<struct svm_node> >& svm_nodes)
  {
    Size n_obs = predictor_values[0].size();
    svm_nodes.resize(n_obs);
    Size skipped_predictors = 0;
    for (Size pred_index = 0; pred_index < predictor_values.size(); 
         ++pred_index)
    {
      if (predictor_values[pred_index].empty()) // uninformative predictor
      {
        skipped_predictors++;
        continue;
      }
      for (Size obs_index = 0; obs_index < n_obs; ++obs_index)
      {
        double value = predictor_values[pred_index][obs_index];
        if (value > 0)
        {
          svm_node node = {pred_index - skipped_predictors, value};
          svm_nodes[obs_index].push_back(node);
        }
      }
    }
    svm_node final = {-1, 0.0};
    for (vector<vector<struct svm_node> >::iterator it = svm_nodes.begin();
         it != svm_nodes.end(); ++it)
    {
      it->push_back(final);
    }
  }


  void optimizeSVMParams_(vector<double> labels,
                          vector<vector<struct svm_node> >& svm_nodes, 
                          struct svm_parameter& svm_params, Size n_parts,
                          const vector<double>& log2_C,
                          const vector<double>& log2_gamma)
  {
    vector<Size> pos_obs; // positive observations (here: "true positives")
    vector<Size> neg_obs; // negative observations (here: "false positives")
    for (Size i = 0; i < labels.size(); ++i)
    {
      if (labels[i] == 0.0) neg_obs.push_back(i);
      else if (labels[i] == 1.0) pos_obs.push_back(i);
    }
    if (pos_obs.size() < n_parts)
    {
      String msg = "not enough positive observations for " + 
        String(n_parts) + "-fold cross-validation";
      throw Exception::MissingInformation(__FILE__, __LINE__, 
                                          __PRETTY_FUNCTION__, msg);
    }
    if (neg_obs.size() < n_parts)
    {
      String msg = ("not enough negative observations for " + 
                    String(n_parts) + "-fold cross-validation");
      throw Exception::MissingInformation(__FILE__, __LINE__, 
                                          __PRETTY_FUNCTION__, msg);
    }
    random_shuffle(pos_obs.begin(), pos_obs.end());
    random_shuffle(neg_obs.begin(), neg_obs.end());

    // Cross-validation: We consider positive and negative observations
    // separately ("stratified sampling"). We split the pos./neg. observations
    // evenly into N partitions. In turn, we use one partition for testing, all
    // the others for training SVMs with different parameters.
    vector<vector<Size> > partitions(n_parts);
    for (; n_parts > 0; --n_parts)
    {
      // integer division should do what we want here:
      Size pos_part = pos_obs.size() / n_parts;
      Size neg_part = neg_obs.size() / n_parts;
      // split partitions off the back of the obs. vector (simplifies removing):
      for (Size i = 1; i <= pos_part; ++i)
      {
        Size index = pos_obs[pos_obs.size() - i];
        partitions[n_parts - 1].push_back(index);
      }
      for (Size i = 1; i <= neg_part; ++i)
      {
        Size index = neg_obs[neg_obs.size() - i];
        partitions[n_parts - 1].push_back(index);
      }
      // remove parts that were just used:
      pos_obs.resize(pos_obs.size() - pos_part);
      neg_obs.resize(neg_obs.size() - neg_part);
    }

    for (Size test_set = 0; test_set < partitions.size(); ++test_set)
    {
      struct svm_problem svm_data;
      svm_data.l = (pos_obs.size() + neg_obs.size() - 
                    partitions[test_set].size());
      svm_data.y = new double[svm_data.l];
      svm_data.x = new svm_node*[svm_data.l];
      Size obs_index = 0;
      for (Size train_set = 0; train_set < partitions.size(); ++train_set)
      {
        if (train_set == test_set) continue;
        for (vector<Size>::iterator it = partitions[train_set].begin();
             it != partitions[train_set].end(); ++it)
        {
          svm_data.x[obs_index] = &(svm_nodes[*it][0]);
          svm_data.y[obs_index] = labels[*it];
          ++obs_index;
        }
      }

      // classification performance for different parameter pairs:
      SVMPerformance performance;
      for (vector<double>::const_iterator c_it = log2_C.begin(); 
           c_it != log2_C.end(); ++c_it)
      {
        svm_params.C = pow(2.0, *c_it);
        for (vector<double>::const_iterator g_it = log2_gamma.begin();
             g_it != log2_gamma.end(); ++g_it)
        {
          svm_params.gamma = pow(2.0, *g_it);
          struct svm_model* model = svm_train(&svm_data, &svm_params);
        }
      }

      delete[] svm_data.y;
      delete[] svm_data.x;
    }
       
  }


  void classifyFeatures_(FeatureMap& features)
  {
    vector<double> labels(features.size());
    for (Size feat_index = 0; feat_index < features.size(); ++feat_index)
    {
      String feature_class = features[feat_index].getMetaValue("feature_class");
      if (feature_class == "true_positive") labels[feat_index] = 1.0;
      else if (feature_class == "unknown") labels[feat_index] = -1.0;
      else if (feature_class == "ambiguous") labels[feat_index] = 0.5;
      // else labels[feat_index] = 0.0; // "false_positive", zero is default
    }

    vector<String> predictors = ListUtils::create<String>(score_metavalues_);
    // values for all featues per predictor (this way around to simplify scaling
    // of predictors):
    vector<vector<double> > predictor_values(predictors.size());
    for (Size pred_index = 0; pred_index < predictors.size(); ++pred_index)
    {
      const String& metavalue = predictors[pred_index];
      predictor_values[pred_index].resize(features.size());
      for (Size feat_index = 0; feat_index < features.size(); ++feat_index)
      {
        predictor_values[pred_index][feat_index] = 
          double(features[feat_index].getMetaValue(metavalue));
      }
    }

    scaleSVMData_(predictor_values, predictors);

    vector<vector<struct svm_node> > svm_nodes;
    convertSVMData_(predictor_values, svm_nodes);

    struct svm_parameter svm_params;
    svm_params.svm_type = C_SVC;
    String svm_kernel = getStringOption_("svm:kernel");
    svm_params.kernel_type = (svm_kernel == "RBF") ? RBF : LINEAR;
    svm_params.eps = 0.001;
    svm_params.cache_size = 100.0;
    svm_params.nr_weight = 0; // use weighting of unbalanced classes?
    svm_params.shrinking = 0; // use shrinking heuristics?
    svm_params.probability = 1;

    Size n_parts = getIntOption_("svm:xval");
    vector<double> log2_C = getDoubleList_("svm:log2_C");
    vector<double> log2_gamma = getDoubleList_("svm:log2_gamma");
    optimizeSVMParams_(labels, svm_nodes, svm_params, n_parts, log2_C, 
                       log2_gamma);
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
    rt_window_ = getDoubleOption_("extract:rt_window");
    mz_window_ = getDoubleOption_("extract:mz_window");
    mz_window_ppm_ = mz_window_ >= 1;
    isotope_pmin_ = getDoubleOption_("extract:isotope_pmin");
    double peak_width = getDoubleOption_("detect:peak_width");
    double min_peak_width = getDoubleOption_("detect:min_peak_width");
    double signal_to_noise = getDoubleOption_("detect:signal_to_noise");
    mapping_tolerance_ = getDoubleOption_("detect:mapping_tolerance");
    elution_model_ = getStringOption_("model:type");

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
    extractor_.setLogType(log_type_);
    Param params = feat_finder_.getParameters();
    params.setValue("stop_report_after_feature", -1); // return all features
    params.setValue("Scores:use_rt_score", "false"); // RT may not be reliable
    if (elution_model_ != "none") params.setValue("write_convex_hull", "true");
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
    feat_finder_.setLogType(log_type_);
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
    }
    for (vector<PeptideIdentification>::iterator pep_it = peptides_ext.begin();
         pep_it != peptides_ext.end(); ++pep_it)
    {
      addPeptideToMap_(*pep_it, peptide_map, true);
    }

    //-------------------------------------------------------------
    // run feature detection
    //-------------------------------------------------------------
    LOG_INFO << "Running feature detection..." << endl;
    FeatureMap features;
    keep_library_ = !lib_out.empty();
    keep_chromatograms_ = !chrom_out.empty();
    for (PeptideMap::iterator pm_it = peptide_map.begin();
         pm_it != peptide_map.end(); ++pm_it)
    {
      FeatureMap current_features;
      detectFeaturesOnePeptide_(*pm_it, current_features);
      features += current_features;
    }
    LOG_DEBUG << "Found " << features.size() << " features in total." << endl;
    ms_data_.reset(); // not needed anymore, free up the memory
    features.setProteinIdentifications(proteins);

    // @FIXME
    // if (elution_model_ != "none")
    // {
    //   fitElutionModels_(features);
    // }

    //-------------------------------------------------------------
    // write output
    //-------------------------------------------------------------
    if (keep_library_)
    {
      removeDuplicateProteins_(library_);
      TraMLFile().store(lib_out, library_);
      library_.clear(true);
    }
    if (keep_chromatograms_)
    {
      addDataProcessing_(chrom_data_,
                         getProcessingInfo_(DataProcessing::FILTERING));
      MzMLFile().store(chrom_out, chrom_data_);
      chrom_data_.clear(true);
    }

    LOG_INFO << "Writing results..." << endl;
    features.ensureUniqueId();
    addDataProcessing_(features,
                       getProcessingInfo_(DataProcessing::QUANTITATION));
    FeatureXMLFile().store(out, features);

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
