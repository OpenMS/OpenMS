// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Timo Sachsenberg, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmMetaboIdent.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/TraceFitter.h>

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

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
    defaults_.setValue("candidates_out", "", "Optional output file with feature candidates.", vector<string>{"output file"});

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
    String spectra_file)
  {
    // if proper mzML is annotated in MS data use this as reference. Otherwise, overwrite with spectra_file information.
    features.setPrimaryMSRunPath({spectra_file}, ms_data_); 

    for (const auto& c : metaboIdentTable)
    {
      addTargetToLibrary_(c.name, c.formula, c.mass, c.charges, c.rts, c.rt_ranges,
                      c.iso_distrib);
    }

    // initialize algorithm classes needed later:
    Param params = feat_finder_.getParameters();
    params.setValue("stop_report_after_feature", -1); // return all features
    params.setValue("Scores:use_rt_score", "false"); // RT may not be reliable
    params.setValue("write_convex_hull", "true");
    if (min_peak_width_ < 1.0)
    {
      min_peak_width_ *= peak_width_;
    }
    params.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width",
                    peak_width_);
    params.setValue("TransitionGroupPicker:min_peak_width", min_peak_width_);
    // disabling the signal-to-noise threshold (setting the parameter to zero)
    // totally breaks the OpenSWATH feature detection (no features found)!
    params.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise",
                    signal_to_noise_);
    
    params.setValue("TransitionGroupPicker:PeakPickerMRM:write_sn_log_messages", "false");     
    params.setValue("TransitionGroupPicker:recalculate_peaks", "true");
    params.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
    params.setValue("TransitionGroupPicker:PeakPickerMRM:method",
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
    
    // sort features:
    sort(features.begin(), features.end(), feature_compare_);

    if (!candidates_out_.empty()) // store feature candidates
    {
      FileHandler().storeFeatures(candidates_out_, features);
    }

    selectFeaturesFromCandidates_(features);
    OPENMS_LOG_INFO << features.size()
             << " features left after selection of best candidates." << endl;

    // get bounding boxes for all mass traces in all features:
    FeatureBoundsMap feature_bounds;
    getFeatureBounds_(features, feature_bounds);
    // find and resolve overlaps:
    vector<FeatureGroup> overlap_groups;
    findOverlappingFeatures_(features, feature_bounds, overlap_groups);
    if (overlap_groups.size() == features.size())
    {
      OPENMS_LOG_INFO << "No overlaps between features found." << endl;
    }
    else
    {
      Size n_overlap_groups = 0, n_overlap_features = 0;
      for (FeatureGroup& group : overlap_groups)
      {
        if (group.size() > 1)
        {
          n_overlap_groups++;
          n_overlap_features += group.size();
          resolveOverlappingFeatures_(group, feature_bounds);
        }
      }
      features.erase(remove_if(features.begin(), features.end(),
                             feature_filter_), features.end());
      OPENMS_LOG_INFO << features.size()
               << " features left after resolving overlaps (involving "
               << n_overlap_features << " features in " << n_overlap_groups
               << " groups)." << endl;
      if (features.empty())
      {
        OPENMS_LOG_INFO << "No features left after filtering." << endl;
      }    
    }

    if (features.empty()) return;

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

  /// Check if two sets of mass trace boundaries overlap
  bool FeatureFinderAlgorithmMetaboIdent::hasOverlappingBounds_(const vector<MassTraceBounds>& mtb1,
                             const vector<MassTraceBounds>& mtb2) const
  {
    for (const MassTraceBounds& mt1 : mtb1)
    {
      for (const MassTraceBounds& mt2 : mtb2)
      {
        if (!((mt1.rt_max < mt2.rt_min) ||
              (mt1.rt_min > mt2.rt_max) ||
              (mt1.mz_max < mt2.mz_min) ||
              (mt1.mz_min > mt2.mz_max)))
        {
          return true;
        }
      }
    }
    return false;
  }

  /// Check if a feature overlaps with a group of other features
  bool FeatureFinderAlgorithmMetaboIdent::hasOverlappingFeature_(const Feature& feature, const FeatureGroup& group,
                              const FeatureBoundsMap& feature_bounds) const
  {
    FeatureBoundsMap::const_iterator fbm_it1 =
      feature_bounds.find(feature.getUniqueId());
    // check overlaps with other features:
    for (FeatureGroup::const_iterator group_it = group.begin();
         group_it != group.end(); ++group_it)
    {
      FeatureBoundsMap::const_iterator fbm_it2 =
        feature_bounds.find((*group_it)->getUniqueId());
      // two features overlap if any of their mass traces overlap:
      if (hasOverlappingBounds_(fbm_it1->second, fbm_it2->second))
      {
        return true;
      }
    }
    return false;
  }


  /// Get bounding boxes for all mass traces in all features of a feature map
  void FeatureFinderAlgorithmMetaboIdent::getFeatureBounds_(const FeatureMap& features,
                         FeatureBoundsMap& feature_bounds)
  {
    for (const  Feature& feat : features)
    {
      for (Size i = 0; i < feat.getSubordinates().size(); ++i)
      {
        MassTraceBounds mtb;
        mtb.sub_index = i;
        const ConvexHull2D::PointArrayType& points =
          feat.getConvexHulls()[i].getHullPoints();
        mtb.mz_min = points.front().getY();
        mtb.mz_max = points.back().getY();
        const Feature& sub = feat.getSubordinates()[i];
        // convex hulls should be written out by "MRMFeatureFinderScoring" (see
        // parameter "write_convex_hull"):
        if (sub.getConvexHulls().empty())
        {
          String error = "convex hulls for mass traces missing";
          throw Exception::MissingInformation(__FILE__, __LINE__,
                                              OPENMS_PRETTY_FUNCTION, error);
        }
        const ConvexHull2D& hull = sub.getConvexHulls()[0];
        // find beginning of mass trace (non-zero intensity):
        if (hull.getHullPoints().empty())
        {
          continue;
        }
        double rt_min = hull.getHullPoints().back().getX();
        for (ConvexHull2D::PointArrayType::const_iterator p_it =
               hull.getHullPoints().begin(); p_it != hull.getHullPoints().end();
             ++p_it)
        {
          if (p_it->getY() > 0)
          {
            rt_min = p_it->getX();
            break;
          }
        }
        // find end of mass trace (non-zero intensity):
        double rt_max = hull.getHullPoints().front().getX();
        for (ConvexHull2D::PointArrayType::const_reverse_iterator p_it =
               hull.getHullPoints().rbegin(); p_it !=
               hull.getHullPoints().rend(); ++p_it)
        {
          if (p_it->getX() < rt_min)
          {
            break;
          }
          if (p_it->getY() > 0)
          {
            rt_max = p_it->getX();
            break;
          }
        }
        if (rt_min > rt_max)
        {
          continue; // no peak -> skip
        }
        mtb.rt_min = rt_min;
        mtb.rt_max = rt_max;
        feature_bounds[feat.getUniqueId()].push_back(mtb);
      }
    }
  }

  /// Partition features of a feature map into groups of overlapping features
  void FeatureFinderAlgorithmMetaboIdent::findOverlappingFeatures_(FeatureMap& features,
                                const FeatureBoundsMap& feature_bounds,
                                vector<FeatureGroup>& overlap_groups)
  {
    for (Feature& feat : features)
    {
      // @TODO: make this more efficient?
      vector<FeatureGroup> current_overlaps;
      vector<FeatureGroup> no_overlaps;
      for (const FeatureGroup& group : overlap_groups)
      {
        if (hasOverlappingFeature_(feat, group, feature_bounds))
        {
          current_overlaps.push_back(group);
        }
        else
        {
          no_overlaps.push_back(group);
        }
      }
      if (current_overlaps.empty()) // make new group for current feature
      {
        FeatureGroup new_group(1, &(feat));
        no_overlaps.push_back(new_group);
      }
      else // merge all groups that overlap the current feature, then add it
      {
        FeatureGroup& merged = current_overlaps.front();
        for (vector<FeatureGroup>::const_iterator group_it =
               ++current_overlaps.begin(); group_it != current_overlaps.end();
             ++group_it)
        {
          merged.insert(merged.end(), group_it->begin(), group_it->end());
        }
        merged.push_back(&feat);
        no_overlaps.push_back(merged);
      }
      overlap_groups.swap(no_overlaps);
    }
  }

  /// Resolve overlapping features by picking the best and removing all others
  void FeatureFinderAlgorithmMetaboIdent::resolveOverlappingFeatures_(FeatureGroup& group,
                                   const FeatureBoundsMap& feature_bounds)
  {
    if (debug_level_ > 0)
    {
      String msg = "Overlapping features: ";
      for (FeatureGroup::const_iterator it = group.begin(); it != group.end();
           ++it)
      {
        if (it != group.begin())
        {
          msg += ", ";
        }
        msg += String((*it)->getMetaValue("PeptideRef")) + " (RT " +
          String(float((*it)->getRT())) + ")";
      }
      OPENMS_LOG_DEBUG << msg << endl;
    }

    Feature* best_feature = 0;
    while (!group.empty())
    {
      double best_rt_delta = numeric_limits<double>::infinity();
      // best feature is the one with min. RT deviation to target:
      for (FeatureGroup::const_iterator it = group.begin(); it != group.end();
           ++it)
      {
        double rt_delta = abs(double((*it)->getMetaValue("rt_deviation")));
        if ((rt_delta < best_rt_delta) ||
            ((rt_delta == best_rt_delta) && ((*it)->getIntensity() >
                                             best_feature->getIntensity())))
        {
          best_rt_delta = rt_delta;
          best_feature = *it;
        }
        else if ((rt_delta == best_rt_delta) && ((*it)->getIntensity() ==
                                                 best_feature->getIntensity()))
        {
          // are the features the same? (@TODO: use "Math::approximatelyEqual"?)
          if (((*it)->getRT() == best_feature->getRT()) &&
              ((*it)->getMZ() == best_feature->getMZ()))
          {
            // update annotations:
            // @TODO: also adjust "formula" and "expected_rt"?
            String label = best_feature->getMetaValue("label");
            label += "/" + String((*it)->getMetaValue("label"));
            best_feature->setMetaValue("label", label);
            StringList alt_refs;
            if (best_feature->metaValueExists("alt_PeptideRef"))
            {
              alt_refs = best_feature->getMetaValue("alt_PeptideRef");
            }
            alt_refs.push_back((*it)->getMetaValue("PeptideRef"));
            best_feature->setMetaValue("alt_PeptideRef", alt_refs);
          }
          else
          {
            OPENMS_LOG_WARN
              << "Warning: cannot decide between equally good feature candidates; picking the first one of "
              << best_feature->getMetaValue("PeptideRef") << " (RT "
              << float(best_feature->getRT()) << ") and "
              << (*it)->getMetaValue("PeptideRef") << " (RT "
              << float((*it)->getRT()) << ")." << endl;
          }
        }
      }
      // we have found a "best" feature, now remove other features that overlap:
      FeatureGroup no_overlaps;
      FeatureBoundsMap::const_iterator fbm_it1 =
        feature_bounds.find(best_feature->getUniqueId());
      for (FeatureGroup::const_iterator it = group.begin(); it != group.end();
           ++it)
      {
        if (*it == best_feature)
        {
          continue;
        }
        FeatureBoundsMap::const_iterator fbm_it2 =
          feature_bounds.find((*it)->getUniqueId());
        if (hasOverlappingBounds_(fbm_it1->second, fbm_it2->second))
        {
          // keep a record of the feature that is getting removed:
          String ref = String((*it)->getMetaValue("PeptideRef")) + " (RT " +
            String(float((*it)->getRT())) + ")";
          StringList overlap_refs;
          if (best_feature->metaValueExists("overlap_removed"))
          {
            overlap_refs = best_feature->getMetaValue("overlap_removed");
          }
          overlap_refs.push_back(ref);
          best_feature->setMetaValue("overlap_removed", overlap_refs);
          (*it)->setMetaValue("FFMetId_remove", ""); // mark for removal
        }
        else
        {
          no_overlaps.push_back(*it);
        }
      }
      group.swap(no_overlaps);
    }
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

