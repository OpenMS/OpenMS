// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h>

#include <OpenMS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/FEATUREFINDER/GaussTraceFitter.h>
#include <OpenMS/FEATUREFINDER/TraceFitter.h>

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>

#include <OpenMS/FORMAT/FileHandler.h>

#include <OpenMS/MATH/MathFunctions.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/PROCESSING/FEATURE/FeatureOverlapFilter.h>

#include <vector>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <random>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace OpenMS
{
  FeatureFinderAlgorithmMetaboIdent::FeatureFinderAlgorithmMetaboIdent() :
    DefaultParamHandler("FeatureFinderAlgorithmMetaboIdent")
  {
    defaults_.setValue("candidates_out", "", "Optional output file: Feature candidates (before filtering and model fitting).", {"output file"});

    defaults_.setValue("extract:mz_window", 10.0, "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)");
    defaults_.setMinFloat("extract:mz_window", 0.0);

    defaults_.setValue(
      "extract:rt_window", 
      0.0, 
      "RT window size (in sec.) for chromatogram extraction. If set, this parameter takes precedence over 'extract:rt_quantile'.",
      vector<string>{"advanced"});
    defaults_.setMinFloat("extract:rt_window", 0.0);

    defaults_.setValue("extract:n_isotopes", 2, "Number of isotopes to include in each peptide assay.");
    defaults_.setMinInt("extract:n_isotopes", 2);
    defaults_.setValue(
      "extract:isotope_pmin",
      0.0, 
      "Minimum probability for an isotope to be included in the assay for a peptide. If set, this parameter takes precedence over 'extract:n_isotopes'.",
      vector<string>{"advanced"});
    defaults_.setMinFloat("extract:isotope_pmin", 0.0);
    defaults_.setMaxFloat("extract:isotope_pmin", 1.0);

    defaults_.setSectionDescription("extract", "Parameters for ion chromatogram extraction");

    defaults_.setValue("detect:peak_width", 60.0, "Expected elution peak width in seconds, for smoothing (Gauss filter). Also determines the RT extration window, unless set explicitly via 'extract:rt_window'.");
    defaults_.setMinFloat("detect:peak_width", 0.0);
    defaults_.setValue(
      "detect:min_peak_width", 
      0.2, 
      "Minimum elution peak width. Absolute value in seconds if 1 or greater, else relative to 'peak_width'.",
      vector<string>{"advanced"});
    defaults_.setMinFloat("detect:min_peak_width", 0.0);

    defaults_.setValue(
      "detect:signal_to_noise", 
      0.8, 
      "Signal-to-noise threshold for OpenSWATH feature detection",
      vector<string>{"advanced"});
    defaults_.setMinFloat("detect:signal_to_noise", 0.1);

    defaults_.setSectionDescription("detect", "Parameters for detecting features in extracted ion chromatograms");  

    // parameters for model fitting (via ElutionModelFitter):
    defaults_.setValue("model:type", "symmetric", "Type of elution model to fit to features");
    defaults_.setValidStrings("model:type", {"symmetric", "asymmetric", "none"});
    defaults_.insert("model:", ElutionModelFitter().getParameters()); // copy parameters
    defaults_.remove("model:asymmetric");

    defaults_.setSectionDescription("model", "Parameters for fitting elution models to features");

    defaults_.setValue("EMGScoring:max_iteration", 100, "Maximum number of iterations for EMG fitting.");
    defaults_.setMinInt("EMGScoring:max_iteration", 1);
    defaults_.setValue("EMGScoring:init_mom", "false", "Alternative initial parameters for fitting through method of moments.");
    defaults_.setValidStrings("EMGScoring:init_mom", {"true","false"});

    defaults_.setSectionDescription("EMGScoring", "Parameters for fitting exp. mod. Gaussians to mass traces.");

    defaults_.setValue("debug", 0, "Debug level for feature detection.", vector<string>{"advanced"});
    defaults_.setMinInt("debug", 0);

    defaultsToParam_();
  }

 void FeatureFinderAlgorithmMetaboIdent::updateMembers_()
  {
    peak_width_ = param_.getValue("detect:peak_width");
    min_peak_width_ = param_.getValue("detect:min_peak_width");
    signal_to_noise_ = param_.getValue("detect:signal_to_noise");

    rt_window_ = param_.getValue("extract:rt_window");
    if (rt_window_ == 0.0)
    {
      // calculate RT window based on other parameters:
      rt_window_ = 4 * peak_width_;
      OPENMS_LOG_INFO << "RT window size calculated as " << rt_window_
                      << " seconds." << endl;
    }

    mz_window_ = param_.getValue("extract:mz_window");
    mz_window_ppm_ = mz_window_ >= 1;

    isotope_pmin_ = param_.getValue("extract:isotope_pmin");

    // extract up to 10 isotopes if minimum probability is larger than 0
    n_isotopes_ = ((isotope_pmin_ > 0.0) ?
                   10 : (int)param_.getValue("extract:n_isotopes"));

    iso_gen_.setMaxIsotope(n_isotopes_);

    elution_model_ = (string)param_.getValue("model:type");

    // debug
    debug_level_ = param_.getValue("debug");
    candidates_out_ = (string)param_.getValue("candidates_out");
  }

  void FeatureFinderAlgorithmMetaboIdent::run(const vector<FeatureFinderAlgorithmMetaboIdent::FeatureFinderMetaboIdentCompound>& metaboIdentTable, 
    FeatureMap& features, 
    const String& spectra_file)
  {
    // if proper mzML is annotated in MS data use this as reference. Otherwise, overwrite with spectra_file information.
    features.setPrimaryMSRunPath({spectra_file}, ms_data_); 
  
    if (ms_data_.empty())
    {
      OPENMS_LOG_WARN << "Warning: No MS1 scans in:"<< spectra_file << endl;      
      return;
    }

    for (const auto& c : metaboIdentTable)
    {
      addTargetToLibrary_(c.getName(), c.getFormula(), c.getMass(), c.getCharges(), c.getRTs(), c.getRTRanges(),
                      c.getIsotopeDistribution());
    }

    // initialize algorithm classes needed later:
    Param params = feat_finder_.getParameters();

    params.setValue("stop_report_after_feature", -1); // return all features
    params.setValue("EMGScoring:max_iteration", param_.getValue("EMGScoring:max_iteration")); // propagate setting to sub algorithms
    params.setValue("EMGScoring:init_mom", param_.getValue("EMGScoring:init_mom")); // propagate setting to sub algorithms
    params.setValue("Scores:use_rt_score", "false"); // RT may not be reliable
    params.setValue("Scores:use_ionseries_scores", "false"); // since FFID only uses MS1 spectra, this is useless
    params.setValue("Scores:use_ms2_isotope_scores", "false"); // since FFID only uses MS1 spectra, this is useless
    params.setValue("Scores:use_ms1_correlation", "false"); // this would be redundant to the "MS2" correlation and since
    // precursor transition = first product transition, additionally biased
    params.setValue("Scores:use_ms1_mi", "false"); // same as above. On MS1 level we basically only care about the "MS1 fullscan" scores
    //TODO for MS1 level scoring there is an additional parameter add_up_spectra with which we can add up spectra
    // around the apex, to complete isotopic envelopes (and therefore make this score more robust).

    params.setValue("write_convex_hull", "true"); // some parts of FFMId expect convex hulls

    if ((elution_model_ != "none") || (!candidates_out_.empty()))
    {
      params.setValue("Scores:use_elution_model_score", "false"); // TODO: test if this works for requantificiation
    }
    else // no elution model
    {
      params.setValue("Scores:use_elution_model_score", "true");  
    }      
    
    if (min_peak_width_ < 1.0)
    {
      min_peak_width_ *= peak_width_;
    }
    params.setValue("TransitionGroupPicker:PeakPickerChromatogram:gauss_width",
                    peak_width_);
    params.setValue("TransitionGroupPicker:min_peak_width", min_peak_width_);
    // disabling the signal-to-noise threshold (setting the parameter to zero)
    // totally breaks the OpenSWATH feature detection (no features found)!
    params.setValue("TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise",
                    signal_to_noise_);
    
    params.setValue("TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages", "false");     
    params.setValue("TransitionGroupPicker:recalculate_peaks", "true");
    params.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", -1.0);
    params.setValue("TransitionGroupPicker:PeakPickerChromatogram:method",
                    "corrected");
    feat_finder_.setParameters(params);
    feat_finder_.setLogType(ProgressLogger::NONE);
    feat_finder_.setStrictFlag(false);

    //-------------------------------------------------------------
    // run feature detection
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Extracting chromatograms..." << endl;
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

    OPENMS_LOG_DEBUG << "Extracted " << chrom_data_.getNrChromatograms()
              << " chromatogram(s)." << endl;

    OPENMS_LOG_INFO << "Detecting chromatographic peaks..." << endl;
    OpenMS_Log_info.remove(cout); // suppress status output from OpenSWATH
    feat_finder_.pickExperiment(chrom_data_, features, library_,
                                TransformationDescription(), ms_data_);
    OpenMS_Log_info.insert(cout);
    OPENMS_LOG_INFO << "Found " << features.size()
                    << " feature candidates in total." << endl;
    ms_data_.reset(); // not needed anymore, free up the memory

    // complete feature annotation:
    annotateFeatures_(features);

    // write auxiliary output:
    // features.setProteinIdentifications(proteins);
    features.ensureUniqueId();
    
    if (!candidates_out_.empty()) // store feature candidates
    {
      sort(features.begin(), features.end(), feature_compare_);
      FileHandler().storeFeatures(candidates_out_, features);
    }



    selectFeaturesFromCandidates_(features);
    OPENMS_LOG_INFO << features.size()
             << " features left after selection of best candidates." << endl;

    constexpr bool CHECK_TRACES_FOR_OVERLAP = true;

    // criterium used to select the best feature amongs overlapping ones (lower = better)
    auto FeatureComparator = [](const Feature& left, const Feature& right)
      {
        double left_rt_delta = std::abs(double(left.getMetaValue("rt_deviation")));
        double right_rt_delta = std::abs(double(right.getMetaValue("rt_deviation")));
        size_t left_intensity = left.getIntensity();
        size_t right_intensity = right.getIntensity();
        return std::tie(left_rt_delta, right_intensity) < std::tie(right_rt_delta, left_intensity); // Note: left and right intensity are swapped because here higher is better
      };

    // callback used to transfer information from an identical overlapping feature with different annotation to the representative on
    auto FeatureOverlapCallback = [](Feature& cluster_representative, Feature& overlap)
      {
        size_t best_intensity = cluster_representative.getIntensity();
        size_t overlap_intensity = overlap.getIntensity();

        if (overlap_intensity != best_intensity) return true; // early out: features are different

        // this part will nearly never be called (e.g., only completely identicial features)
        // so it is ok to perform some slow operations like querying meta values 
        double best_rt_delta = std::abs(double(cluster_representative.getMetaValue("rt_deviation")));
        double overlap_rt_delta = std::abs(double(overlap.getMetaValue("rt_deviation")));

        if (overlap_rt_delta == best_rt_delta)
        {
          double best_RT = cluster_representative.getRT();
          double overlap_RT = overlap.getRT();
          double best_MZ = cluster_representative.getMZ();
          double overlap_MZ = overlap.getMZ();

          // are the features the same? (@TODO: use "Math::approximatelyEqual"?)
          if ((overlap_MZ == best_MZ) && (overlap_RT == best_RT))
          {
            // update annotations:
            // @TODO: also adjust "formula" and "expected_rt"?
            String label = cluster_representative.getMetaValue("label");            
            label += "/" + String(overlap.getMetaValue("label"));
            cluster_representative.setMetaValue("label", label);
            StringList alt_refs;
            if (cluster_representative.metaValueExists("alt_PeptideRef"))
            {
              alt_refs = cluster_representative.getMetaValue("alt_PeptideRef");
            }
            alt_refs.push_back(overlap.getMetaValue("PeptideRef"));
            cluster_representative.setMetaValue("alt_PeptideRef", alt_refs);
          }
        }

        // annotate which features were removed because of overlap with the representative feature
        String ref = String(overlap.getMetaValue("PeptideRef")) + " (RT " +
          String(float(overlap.getRT())) + ")";

        StringList overlap_refs = cluster_representative.getMetaValue("overlap_removed", StringList{});
        overlap_refs.push_back(std::move(ref));
        cluster_representative.setMetaValue("overlap_removed", std::move(overlap_refs)); // TODO: implement setMetaValue that takes DataValue as r-value reference &&

        return true;
      };

    FeatureOverlapFilter::filter(features, FeatureComparator, FeatureOverlapCallback, CHECK_TRACES_FOR_OVERLAP);
    std::stable_sort(features.begin(), features.end(), feature_compare_); // sort by ref and rt

    if (features.empty())
    {
      OPENMS_LOG_INFO << "No features left after filtering." << endl;
      return;
    }    

    n_shared_ = addTargetAnnotations_(features);

    if (elution_model_ != "none")
    {
      ElutionModelFitter emf;
      Param emf_params = param_.copy("model:", true);
      emf_params.remove("type");
      emf_params.setValue("asymmetric",
                          (elution_model_ == "asymmetric") ? "true" : "false");
      emf.setParameters(emf_params);
      emf.fitElutionModels(features);
    }
    else if (!candidates_out_.empty()) // hulls not needed, remove them
    {
      for (Feature& feat : features)
      {
        for (Feature& sub : feat.getSubordinates())
        {
          sub.getConvexHulls().clear();
        }
      }
    }

    extractTransformations_(features);
  }

  /// Calculate mass-to-charge ratio from mass and charge
  double FeatureFinderAlgorithmMetaboIdent::calculateMZ_(double mass, Int charge) const
  {
    return (mass + charge * Constants::PROTON_MASS_U) / fabs(charge);
  }

  /// Add a target (from the input file) to the assay library
  void FeatureFinderAlgorithmMetaboIdent::addTargetToLibrary_(const String& name, const String& formula,
                           double mass, const vector<Int>& charges,
                           const vector<double>& rts,
                           vector<double> rt_ranges,
                           const vector<double>& iso_distrib)
  {
    if ((mass <= 0) && formula.empty())
    {
      OPENMS_LOG_ERROR << "Error: No mass or sum formula given for target '"
                       << name << "' - skipping this target." << endl;
      return;
    }
    if (rts.empty())
    {
      OPENMS_LOG_ERROR << "Error: No retention time (RT) given for target '"
                       << name << "' - skipping this target." << endl;
      return;
    }
    // @TODO: detect entries with same RT and m/z ("collisions")
    TargetedExperiment::Compound target;
    target.setMetaValue("name", name);
    target.molecular_formula = formula;
    EmpiricalFormula emp_formula(formula);
    bool mass_given = (mass > 0);
    if (!mass_given)
    {
      mass = emp_formula.getMonoWeight();
    }
    target.theoretical_mass = mass;
    String target_id = name + "_m" + String(float(mass));

    // get isotope distribution for target:
    IsotopeDistribution iso_dist;
    Size n_isotopes = n_isotopes_;
    if (iso_distrib.empty() || (iso_distrib[0] == 0))
    {
      if (formula.empty())
      {
        OPENMS_LOG_ERROR << "Error: No sum formula given for target '" << name
                         << "'; cannot calculate isotope distribution"
                         << " - using estimation method for peptides." << endl;
        iso_dist = iso_gen_.estimateFromPeptideWeight(mass);
      }
      else
      {
        iso_dist = emp_formula.getIsotopeDistribution(iso_gen_);
      }
    }
    else
    {
      n_isotopes = min(n_isotopes, iso_distrib.size());
      IsotopeDistribution::ContainerType probs;
      probs.reserve(n_isotopes);
      for (Size i = 0; i < n_isotopes; ++i)
      {
        probs.push_back(Peak1D(i, iso_distrib[i]));
      }
      iso_dist.set(probs);
    }
    if (isotope_pmin_ > 0.0)
    {
      iso_dist.trimLeft(isotope_pmin_);
      iso_dist.trimRight(isotope_pmin_);
    }
    iso_dist.renormalize();

    // go through different charge states:
    for (vector<Int>::const_iterator z_it = charges.begin();
         z_it != charges.end(); ++z_it)
    {
      if (*z_it == 0)
      {
        OPENMS_LOG_ERROR << "Error: Invalid charge 0 for target '" << name
                         << "' - skipping this charge." << endl;
        continue;
      }
      target.setChargeState(*z_it);
      double mz = 0.0;
      if (!mass_given) // calculate m/z from formula
      {
        emp_formula.setCharge(*z_it);
        // "EmpiricalFormula::getMonoWeight()" already includes charges:
        mz = abs(emp_formula.getMonoWeight() / *z_it);
      }
      else
      {
        mz = calculateMZ_(mass, *z_it);
      }

      // recycle to one range entry per RT:
      if (rt_ranges.empty())
      {
        rt_ranges.resize(rts.size(), 0.0);
      }
      else if (rt_ranges.size() == 1)
      {
        rt_ranges.resize(rts.size(), rt_ranges[0]);
      }

      for (Size i = 0; i < rts.size(); ++i)
      {
        target.id = target_id + "_z" + String(*z_it) + "_rt" +
          String(float(rts[i]));
        target.setMetaValue("expected_rt", rts[i]);
        target_rts_[target.id] = rts[i];

        double rt_tol = rt_ranges[i] / 2.0;
        if (rt_tol == 0)
        {
          rt_tol = rt_window_ / 2.0;
        }
        // store beginning and end of RT region:
        target.rts.clear();
        addTargetRT_(target, rts[i] - rt_tol);
        addTargetRT_(target, rts[i] + rt_tol);
        library_.addCompound(target);
        generateTransitions_(target.id, mz, *z_it, iso_dist);
      }
    }
  }

  /// Generate transitions for a target ion and add them to the library
  void FeatureFinderAlgorithmMetaboIdent::generateTransitions_(const String& target_id, double mz, Int charge,
                            const IsotopeDistribution& iso_dist)
  {
    // go through different isotopes:
    Size counter = 0;
    for (const Peak1D& iso : iso_dist)
    {
      ReactionMonitoringTransition transition;
      String annotation = "i" + String(counter);
      String transition_name = target_id + "_" + annotation;

      transition.setNativeID(transition_name);
      transition.setPrecursorMZ(mz);
      // @TODO: use accurate masses from the isotope distribution here?
      transition.setProductMZ(mz + abs(Constants::C13C12_MASSDIFF_U *
                                       float(counter) / charge));
      transition.setLibraryIntensity(iso.getIntensity());
      // transition.setMetaValue("annotation", annotation); // ???
      transition.setCompoundRef(target_id);
      library_.addTransition(transition);
      isotope_probs_[transition_name] = iso.getIntensity();
      
      ++counter;
    }
  }

  /// Helper function to add retention time to a target
  void FeatureFinderAlgorithmMetaboIdent::addTargetRT_(TargetedExperiment::Compound& target, double rt)
  {
    TargetedExperiment::RetentionTime te_rt;
    te_rt.retention_time_unit =
      TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
    te_rt.retention_time_type =
      TargetedExperimentHelper::RetentionTime::RTType::LOCAL;
    te_rt.setRT(rt);
    target.rts.push_back(te_rt);
  }


  /// Add relevant annotations/meta values to features
  void FeatureFinderAlgorithmMetaboIdent::annotateFeatures_(FeatureMap& features)
  {
    for (Feature& feat : features)
    {
      feat.setMZ(feat.getMetaValue("PrecursorMZ"));
      String ref = feat.getMetaValue("PeptideRef");
      const TargetedExperiment::Compound& compound =
        library_.getCompoundByRef(ref);
      feat.setCharge(compound.getChargeState());
      ensureConvexHulls_(feat);
      feat.getPeptideIdentifications().clear();
      feat.setMetaValue("label", compound.getMetaValue("name"));
      feat.setMetaValue("sum_formula", compound.molecular_formula);
      feat.setMetaValue("expected_rt",
                            compound.getMetaValue("expected_rt"));
      // annotate subordinates with theoretical isotope intensities:
      for (Feature& sub : feat.getSubordinates())
      {
        String native_id = sub.getMetaValue("native_id");
        sub.setMetaValue("isotope_probability", isotope_probs_[native_id]);
        sub.removeMetaValue("FeatureLevel"); // value "MS2" is misleading
      }
      // annotate num_mass_traces, required for SIRIUS
      feat.setMetaValue(Constants::UserParam::NUM_OF_MASSTRACES, feat.getSubordinates().size());
    }
    features.getProteinIdentifications().clear();
  }

  /// Create hulls for mass traces of a feature, if not already present
  void FeatureFinderAlgorithmMetaboIdent::ensureConvexHulls_(Feature& feature) const
  {
    if (feature.getConvexHulls().empty())
    {
      double rt_min = feature.getMetaValue("leftWidth");
      double rt_max = feature.getMetaValue("rightWidth");
      for (Feature& sub : feature.getSubordinates())
      {
        double abs_mz_tol = mz_window_ / 2.0;
        if (mz_window_ppm_)
        {
          abs_mz_tol = sub.getMZ() * abs_mz_tol * 1.0e-6;
        }
        ConvexHull2D hull;
        hull.addPoint(DPosition<2>(rt_min, sub.getMZ() - abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_min, sub.getMZ() + abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_max, sub.getMZ() - abs_mz_tol));
        hull.addPoint(DPosition<2>(rt_max, sub.getMZ() + abs_mz_tol));
        feature.getConvexHulls().push_back(hull);
      }
    }
  }

  /// Select the best feature for an assay from a set of candidates
  void FeatureFinderAlgorithmMetaboIdent::selectFeaturesFromCandidates_(FeatureMap& features)
  {
    String previous_ref;
    double best_rt_dist = numeric_limits<double>::infinity();
    FeatureMap::Iterator best_it = features.begin();
    for (FeatureMap::Iterator it = features.begin(); it != features.end();
         ++it)
    {
      // features from same assay (same "PeptideRef") appear consecutively:
      String ref = it->getMetaValue("PeptideRef");
      if (ref != previous_ref) // new assay
      {
        previous_ref = ref;
        best_rt_dist = rt_window_;
        best_it = it;
      }
      double target_rt = target_rts_[ref];
      double rt_min = it->getMetaValue("leftWidth");
      double rt_max = it->getMetaValue("rightWidth");
      double rt_dist = numeric_limits<double>::infinity();
      if ((rt_min <= target_rt) && (rt_max >= target_rt))
      {
        if (best_rt_dist <= 0.0)
        {
          OPENMS_LOG_WARN
            << "Warning: overlapping feature candidates for assay '" << ref
            << "'" << endl;
        }
        rt_dist = 0.0;
      }
      else if (best_rt_dist > 0.0)
      {
        rt_dist = (rt_min > target_rt) ? (rt_min - target_rt) : (target_rt -
                                                                 rt_max);
      }
      if ((rt_dist < best_rt_dist) ||
          ((rt_dist == best_rt_dist) && (it->getIntensity() >
                                         best_it->getIntensity())))
      {
        // new best candidate for this assay:
        best_rt_dist = rt_dist;
        // mark no-longer-best candidate for removal:
        if (best_it != it) best_it->setMetaValue("FFMetId_remove", "");
        best_it = it;
        best_it->setMetaValue("rt_deviation", target_rt - best_it->getRT());
      }
      else // this candidate is worse than a previous one
      {
        it->setMetaValue("FFMetId_remove", ""); // mark for removal
      }
    }
    features.erase(remove_if(features.begin(), features.end(),
                             feature_filter_), features.end());
  }

  /// Create a string of identifying information for a compound
  String FeatureFinderAlgorithmMetaboIdent::prettyPrintCompound(const TargetedExperiment::Compound& compound)
  {
    return (String(compound.getMetaValue("name")) + " (m=" +
            String(float(compound.theoretical_mass)) + ", z=" +
            String(compound.getChargeState()) + ", rt=" +
            String(float(double(compound.getMetaValue("expected_rt")))) + ")");
  }

  /// Add "peptide" identifications with information about targets to features
  Size FeatureFinderAlgorithmMetaboIdent::addTargetAnnotations_(FeatureMap& features)
  {
    Size n_shared = 0;
    set<String> found_refs;
    for (FeatureMap::Iterator it = features.begin(); it != features.end(); ++it)
    {
      found_refs.insert(it->getMetaValue("PeptideRef"));
      if (it->metaValueExists("alt_PeptideRef"))
      {
        n_shared++;
        StringList alt_refs = it->getMetaValue("alt_PeptideRef");
        found_refs.insert(alt_refs.begin(), alt_refs.end());
      }
    }
    // targets without features:
    size_t n_missing = library_.getCompounds().size() - found_refs.size();
    features.getUnassignedPeptideIdentifications().reserve(n_missing);
    for (vector<TargetedExperiment::Compound>::const_iterator it =
           library_.getCompounds().begin(); it != library_.getCompounds().end();
         ++it)
    {
      if (!found_refs.count(it->id))
      {
        PeptideIdentification peptide;
        peptide.setIdentifier("id");
        peptide.setMetaValue("label", it->getMetaValue("name"));
        peptide.setMetaValue("PeptideRef", it->id);
        peptide.setRT(it->getMetaValue("expected_rt"));
        peptide.setMZ(calculateMZ_(it->theoretical_mass, it->getChargeState()));
        features.getUnassignedPeptideIdentifications().push_back(peptide);
      }
      if (features.getUnassignedPeptideIdentifications().size() >= n_missing)
      {
        break; // found all
      }
    }
    if (n_missing)
    {
      features.getProteinIdentifications().resize(1);
      features.getProteinIdentifications()[0].setIdentifier("id");
    }
    return n_shared; // for summary statistics
  }

  void FeatureFinderAlgorithmMetaboIdent::extractTransformations_(const FeatureMap& features)
  {
    TransformationDescription::DataPoints points;
    for (const auto& f : features)
    {
      TransformationDescription::DataPoint point;
      point.first = f.getMetaValue("expected_rt");
      point.second = f.getRT();
      point.note = f.getMetaValue("PeptideRef");
      points.push_back(point);
    }
    trafo_.setDataPoints(points);
  }

  void FeatureFinderAlgorithmMetaboIdent::setMSData(const PeakMap& m)
  { 
    ms_data_ = m; 
    
    vector<MSSpectrum>& specs = ms_data_.getSpectra();

    // keep only MS1
    specs.erase(
      std::remove_if(specs.begin(), specs.end(),
        [](const MSSpectrum & s) { return s.getMSLevel() != 1; }),
      specs.end());
  }

  void FeatureFinderAlgorithmMetaboIdent::setMSData(PeakMap&& m)
  { 
    ms_data_ = std::move(m); 
    
    vector<MSSpectrum>& specs = ms_data_.getSpectra();

    // keep only MS1
    specs.erase(
      std::remove_if(specs.begin(), specs.end(),
        [](const MSSpectrum & s) { return s.getMSLevel() != 1; }),
      specs.end());
  }

}

