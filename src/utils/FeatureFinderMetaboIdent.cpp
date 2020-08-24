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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ElutionModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussTraceFitter.h>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page UTILS_FeatureFinderMetaboIdent FeatureFinderMetaboIdent

   @brief Detects features in MS1 data corresponding to small molecule identifications.

   <CENTER>
     <table>
       <tr>
         <td ALIGN="center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
         <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ FeatureFinderMetaboIdent \f$ \longrightarrow \f$</td>
         <td ALIGN="center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
       </tr>
       <tr>
         <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes (optional) </td>
         <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_TextExporter</td>
       </tr>
     </table>
   </CENTER>

   This tool detects quantitative features in MS1 data for a list of targets, typically small molecule/metabolite identifications.
   It uses algorithms for targeted data analysis from the OpenSWATH pipeline.

   @note This tool is still experimental!

   @see @ref TOPP_FeatureFinderIdentification - targeted feature detection based on peptide identifications.

   <B>Input format</B>

   The targets to quantify have to be specified in a tab-separated text file that is passed via the @p id parameter.
   This file has to start with the following header line, defining its columns:
   <pre>
   <TT>CompoundName    SumFormula    Mass    Charge    RetentionTime    RetentionTimeRange    IsoDistribution</TT>
   </pre>

   Every subsequent line defines a target.
   (Except lines starting with "#", which are considered as comments and skipped.)
   The following requirements apply:
   - @p CompoundName: unique name for the target compound
   - @p SumFormula: chemical sum formula (see @ref OpenMS::EmpiricalFormula), optional
   - @p Mass: neutral mass; if zero calculated from @p Formula
   - @p Charge: charge state, or comma-separated list of multiple charges
   - @p RetentionTime: retention time (RT), or comma-separated list of multiple RTs
   - @p RetentionTimeRange: RT window around @p RetentionTime for chromatogram extraction, either one value or one per @p RT entry; if zero parameter @p extract:rt_window is used
   - @p IsoDistribution: comma-separated list of relative abundances of isotopologues (see @ref OpenMS::IsotopeDistribution); if zero calculated from @p Formula

   In the simplest case, only @p CompoundName, @p SumFormula, @p Charge and @p RetentionTime need to be given, all other values may be zero.
   Every combination of compound (mass), RT and charge defines one target for feature detection.

   <B>Output format</B>

   The main output (parameter @p out) is a featureXML file containing the detected features, with annotations in meta data entries.
   This file can be visualized in TOPPView - perhaps most usefully as a layer on top of the LC-MS data that gave rise to it.
   Compound annotations of features (@p Name entries from the @p id input) can be shown by clicking the "Show feature annotation" button in the tool bar and selecting "Label meta data".
   Positions of targets for which no feature was detected can be shown by clicking the "Show unassigned peptide identifications" button and selecting "Show label meta data".

   To export the data from the featureXML file to a tabular text file (CSV), use @ref TOPP_TextExporter with the options @p no_ids and <TT>feature:add_metavalues 0</TT> (to include all meta data annotations).
   In the result, the information from the @p CompoundName, @p SumFormula, @p Charge and @p RetentionTime columns from the input will be in the @p label, @p sum_formula, @p charge and @p expected_rt columns, respectively.

   <B>The command line parameters of this tool are:</B>
   @verbinclude UTILS_FeatureFinderMetaboIdent.cli
   <B>INI file documentation of this tool:</B>
   @htmlinclude UTILS_FeatureFinderMetaboIdent.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPFeatureFinderMetaboIdent :
  public TOPPBase
{
public:
  TOPPFeatureFinderMetaboIdent() :
    TOPPBase("FeatureFinderMetaboIdent", "Detects features in MS1 data based on metabolite identifications.", false),
    keep_chromatograms_(false), keep_library_(false), rt_window_(0.0),
    mz_window_(0.0), mz_window_ppm_(false), isotope_pmin_(0.0), n_isotopes_(0)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file: LC-MS raw data");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("id", "<file>", "", "Input file: Metabolite identifications");
    setValidFormats_("id", ListUtils::create<String>("tsv"));
    registerOutputFile_("out", "<file>", "", "Output file: Features");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    registerOutputFile_("lib_out", "<file>", "", "Output file: Assay library", false);
    setValidFormats_("lib_out", ListUtils::create<String>("traML"));
    registerOutputFile_("chrom_out", "<file>", "", "Output file: Chromatograms", false);
    setValidFormats_("chrom_out", ListUtils::create<String>("mzML"));
    registerOutputFile_("candidates_out", "<file>", "", "Output file: Feature candidates (before filtering and model fitting)", false);
    registerOutputFile_("trafo_out", "<file>", "", "Output file: Retention times (expected vs. observed)", false);
    setValidFormats_("trafo_out", ListUtils::create<String>("trafoXML"));
    setValidFormats_("candidates_out", ListUtils::create<String>("featureXML"));

    registerTOPPSubsection_("extract", "Parameters for ion chromatogram extraction");
    registerDoubleOption_("extract:mz_window", "<value>", 10.0, "m/z window size for chromatogram extraction (unit: ppm if 1 or greater, else Da/Th)", false);
    setMinFloat_("extract:mz_window", 0.0);
    registerIntOption_("extract:n_isotopes", "<number>", 2, "Number of isotopes to include in each assay.", false);
    setMinInt_("extract:n_isotopes", 2);
    registerDoubleOption_("extract:isotope_pmin", "<value>", 0.0, "Minimum probability for an isotope to be included in the assay for a compound. If set, this parameter takes precedence over 'extract:n_isotopes'.", false, true);
    setMinFloat_("extract:isotope_pmin", 0.0);
    setMaxFloat_("extract:isotope_pmin", 1.0);
    registerDoubleOption_("extract:rt_window", "<value>", 0.0, "RT window size (in sec.) for chromatogram extraction. If zero, calculated based on 'detect:peak_width'.", false);
    setMinFloat_("extract:rt_window", 0.0);

    registerTOPPSubsection_("detect", "Parameters for detecting features in extracted ion chromatograms");
    registerDoubleOption_("detect:peak_width", "<value>", 5.0, "Expected elution peak width in seconds, for smoothing (Gauss filter). Also determines the RT extration window, unless set explicitly via 'extract:rt_window'.", false);
    setMinFloat_("detect:peak_width", 0.0);
    registerDoubleOption_("detect:min_peak_width", "<value>", 0.2, "Minimum elution peak width. Absolute value in seconds if 1 or greater, else relative to 'peak_width'.", false, true);
    setMinFloat_("detect:min_peak_width", 0.0);
    registerDoubleOption_("detect:signal_to_noise", "<value>", 0.8, "Signal-to-noise threshold for OpenSWATH feature detection", false, true);
    setMinFloat_("detect:signal_to_noise", 0.1);

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

  typedef vector<Feature*> FeatureGroup; ///< group of (overlapping) features

  /// Boundaries for a mass trace in a feature
  struct MassTraceBounds
  {
    Size sub_index;
    double rt_min, rt_max, mz_min, mz_max;
  };

  /// Boundaries for all mass traces per feature
  typedef map<UInt64, vector<MassTraceBounds> > FeatureBoundsMap;

  /// Predicate for filtering features by overall quality
  struct FeatureFilterQuality
  {
    bool operator()(const Feature& feature)
    {
      return feature.metaValueExists("FFMetId_remove");
    }
  } feature_filter_;

  /// Comparison functor for features
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

  PeakMap ms_data_; ///< input LC-MS data
  PeakMap chrom_data_; ///< accumulated chromatograms (XICs)
  bool keep_chromatograms_; ///< keep chromatogram data for output?
  TargetedExperiment library_; ///< accumulated assays for targets
  bool keep_library_; ///< keep assay data for output?
  double rt_window_; ///< RT window width
  double mz_window_; ///< m/z window width
  bool mz_window_ppm_; ///< m/z window width is given in PPM (not Da)?
  double isotope_pmin_; ///< min. isotope probability for MS1 assay
  Size n_isotopes_; ///< number of isotopes for MS1 assay
  CoarseIsotopePatternGenerator iso_gen_; ///< isotope pattern generator
  map<String, double> isotope_probs_; ///< isotope probabilities of transitions
  map<String, double> target_rts_; ///< RTs of targets (assays)
  MRMFeatureFinderScoring feat_finder_; ///< OpenSWATH feature finder
  ProgressLogger prog_log_; ///< progress logger

  /// Read input file with information about targets
  void readTargets_(const String& in_path)
  {
    const string header =
      "CompoundName\tSumFormula\tMass\tCharge\tRetentionTime\tRetentionTimeRange\tIsoDistribution";
    ifstream source(in_path.c_str());
    if (!source.is_open())
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, in_path);
    }
    string line;
    getline(source, line);
    if (!String(line).hasPrefix(header))
    {
      String msg = "expected header line starting with: '" + header + "'";
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  line, msg);
    }
    Size line_count = 1;
    set<String> names;
    while (getline(source, line))
    {
      line_count++;
      if (line[0] == '#') continue; // skip comments
      vector<String> parts = ListUtils::create<String>(line, '\t'); // split
      if (parts.size() < 7)
      {
        OPENMS_LOG_ERROR
          << "Error: Expected 7 tab-separated fields, found only "
          << parts.size() << " in line " << line_count
          << " - skipping this line." << endl;
        continue;
      }
      String name = parts[0];
      if (name.empty())
      {
        OPENMS_LOG_ERROR << "Error: Empty name field in input line "
                         << line_count << " - skipping this line." << endl;
        continue;
      }
      if (!names.insert(name).second) // @TODO: is this check necessary?
      {
        OPENMS_LOG_ERROR << "Error: Duplicate name '" << name
                         << "' in input line " << line_count
                         << " - skipping this line." << endl;
        continue;
      }
      String formula = parts[1];
      double mass = parts[2].toDouble();
      vector<Int> charges = ListUtils::create<Int>(parts[3]);
      vector<double> rts = ListUtils::create<double>(parts[4]);
      vector<double> rt_ranges = ListUtils::create<double>(parts[5]);
      vector<double> iso_distrib = ListUtils::create<double>(parts[6]);
      addTargetToLibrary_(name, formula, mass, charges, rts, rt_ranges,
                          iso_distrib);
    }
  }

  /// Calculate mass-to-charge ratio from mass and charge
  double calculateMZ_(double mass, Int charge)
  {
    return (mass + charge * Constants::PROTON_MASS_U) / fabs(charge);
  }

  /// Add a target (from the input file) to the assay library
  void addTargetToLibrary_(const String& name, const String& formula,
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
    if (!mass_given) mass = emp_formula.getMonoWeight();
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
        if (rt_tol == 0) rt_tol = rt_window_ / 2.0;

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
  void generateTransitions_(const String& target_id, double mz, Int charge,
                            const IsotopeDistribution& iso_dist)
  {
    // go through different isotopes:
    Size counter = 0;
    for (IsotopeDistribution::ConstIterator iso_it = iso_dist.begin();
         iso_it != iso_dist.end(); ++iso_it, ++counter)
    {
      ReactionMonitoringTransition transition;
      String annotation = "i" + String(counter);
      String transition_name = target_id + "_" + annotation;

      transition.setNativeID(transition_name);
      transition.setPrecursorMZ(mz);
      // @TODO: use accurate masses from the isotope distribution here?
      transition.setProductMZ(mz + abs(Constants::C13C12_MASSDIFF_U *
                                       float(counter) / charge));
      transition.setLibraryIntensity(iso_it->getIntensity());
      // transition.setMetaValue("annotation", annotation); // ???
      transition.setCompoundRef(target_id);
      library_.addTransition(transition);
      isotope_probs_[transition_name] = iso_it->getIntensity();
    }
  }

  /// Helper function to add retention time to a target
  void addTargetRT_(TargetedExperiment::Compound& target, double rt)
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
  bool hasOverlappingBounds_(const vector<MassTraceBounds>& mtb1,
                             const vector<MassTraceBounds>& mtb2)
  {
    for (vector<MassTraceBounds>::const_iterator mtb1_it = mtb1.begin();
         mtb1_it != mtb1.end(); ++mtb1_it)
    {
      for (vector<MassTraceBounds>::const_iterator mtb2_it = mtb2.begin();
           mtb2_it != mtb2.end(); ++mtb2_it)
      {
        if (!((mtb1_it->rt_max < mtb2_it->rt_min) ||
              (mtb1_it->rt_min > mtb2_it->rt_max) ||
              (mtb1_it->mz_max < mtb2_it->mz_min) ||
              (mtb1_it->mz_min > mtb2_it->mz_max)))
        {
          return true;
        }
      }
    }
    return false;
  }

  /// Check if a feature overlaps with a group of other features
  bool hasOverlappingFeature_(const Feature& feature, const FeatureGroup& group,
                              const FeatureBoundsMap& feature_bounds)
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
  void getFeatureBounds_(const FeatureMap& features,
                         FeatureBoundsMap& feature_bounds)
  {
    for (FeatureMap::ConstIterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it)
    {
      for (Size i = 0; i < feat_it->getSubordinates().size(); ++i)
      {
        MassTraceBounds mtb;
        mtb.sub_index = i;
        const ConvexHull2D::PointArrayType& points =
          feat_it->getConvexHulls()[i].getHullPoints();
        mtb.mz_min = points.front().getY();
        mtb.mz_max = points.back().getY();
        const Feature& sub = feat_it->getSubordinates()[i];
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
        if (hull.getHullPoints().empty()) continue;
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
          if (p_it->getX() < rt_min) break;
          if (p_it->getY() > 0)
          {
            rt_max = p_it->getX();
            break;
          }
        }
        if (rt_min > rt_max) continue; // no peak -> skip
        mtb.rt_min = rt_min;
        mtb.rt_max = rt_max;
        feature_bounds[feat_it->getUniqueId()].push_back(mtb);
      }
    }
  }

  /// Partition features of a feature map into groups of overlapping features
  void findOverlappingFeatures_(FeatureMap& features,
                                const FeatureBoundsMap& feature_bounds,
                                vector<FeatureGroup>& overlap_groups)
  {
    for (FeatureMap::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it)
    {
      // @TODO: make this more efficient?
      vector<FeatureGroup> current_overlaps;
      vector<FeatureGroup> no_overlaps;
      for (vector<FeatureGroup>::const_iterator group_it =
             overlap_groups.begin(); group_it != overlap_groups.end();
           ++group_it)
      {
        if (hasOverlappingFeature_(*feat_it, *group_it, feature_bounds))
        {
          current_overlaps.push_back(*group_it);
        }
        else
        {
          no_overlaps.push_back(*group_it);
        }
      }
      if (current_overlaps.empty()) // make new group for current feature
      {
        FeatureGroup new_group(1, &(*feat_it));
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
        merged.push_back(&(*feat_it));
        no_overlaps.push_back(merged);
      }
      overlap_groups.swap(no_overlaps);
    }
  }

  /// Resolve overlapping features by picking the best and removing all others
  void resolveOverlappingFeatures_(FeatureGroup& group,
                                   const FeatureBoundsMap& feature_bounds)
  {
    if (debug_level_ > 0)
    {
      String msg = "Overlapping features: ";
      for (FeatureGroup::const_iterator it = group.begin(); it != group.end();
           ++it)
      {
        if (it != group.begin()) msg += ", ";
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
        if (*it == best_feature) continue;
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
  void annotateFeatures_(FeatureMap& features)
  {
    for (FeatureMap::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it)
    {
      feat_it->setMZ(feat_it->getMetaValue("PrecursorMZ"));
      String ref = feat_it->getMetaValue("PeptideRef");
      const TargetedExperiment::Compound& compound =
        library_.getCompoundByRef(ref);
      feat_it->setCharge(compound.getChargeState());
      ensureConvexHulls_(*feat_it);
      feat_it->getPeptideIdentifications().clear();
      feat_it->setMetaValue("label", compound.getMetaValue("name"));
      feat_it->setMetaValue("sum_formula", compound.molecular_formula);
      feat_it->setMetaValue("expected_rt",
                            compound.getMetaValue("expected_rt"));
      // annotate subordinates with theoretical isotope intensities:
      for (vector<Feature>::iterator sub_it =
             feat_it->getSubordinates().begin(); sub_it !=
             feat_it->getSubordinates().end(); ++sub_it)
      {
        String native_id = sub_it->getMetaValue("native_id");
        sub_it->setMetaValue("isotope_probability", isotope_probs_[native_id]);
        sub_it->removeMetaValue("FeatureLevel"); // value "MS2" is misleading
      }
    }
    features.getProteinIdentifications().clear();
  }

  /// Create hulls for mass traces of a feature, if not already present
  void ensureConvexHulls_(Feature& feature)
  {
    if (feature.getConvexHulls().empty())
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

  /// Select the best feature for an assay from a set of candidates
  void selectFeaturesFromCandidates_(FeatureMap& features)
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
  String prettyPrintCompound_(const TargetedExperiment::Compound& compound)
  {
    return (String(compound.getMetaValue("name")) + " (m=" +
            String(float(compound.theoretical_mass)) + ", z=" +
            String(compound.getChargeState()) + ", rt=" +
            String(float(double(compound.getMetaValue("expected_rt")))) + ")");
  }

  /// Add "peptide" identifications with information about targets to features
  Size addTargetAnnotations_(FeatureMap& features)
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
    Size n_missing = library_.getCompounds().size() - found_refs.size();
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

  /// Main function
  ExitCodes main_(int, const char**) override
  {
    FeatureMap features;

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String id = getStringOption_("id");
    String out = getStringOption_("out");
    String candidates_out = getStringOption_("candidates_out");
    String lib_out = getStringOption_("lib_out");
    String chrom_out = getStringOption_("chrom_out");
    String trafo_out = getStringOption_("trafo_out");
    rt_window_ = getDoubleOption_("extract:rt_window");
    mz_window_ = getDoubleOption_("extract:mz_window");
    mz_window_ppm_ = mz_window_ >= 1;
    isotope_pmin_ = getDoubleOption_("extract:isotope_pmin");
    n_isotopes_ = ((isotope_pmin_ > 0.0) ?
                   10 : getIntOption_("extract:n_isotopes"));
    iso_gen_.setMaxIsotope(n_isotopes_);
    double peak_width = getDoubleOption_("detect:peak_width");
    double min_peak_width = getDoubleOption_("detect:min_peak_width");
    double signal_to_noise = getDoubleOption_("detect:signal_to_noise");
    String elution_model = getStringOption_("model:type");
    prog_log_.setLogType(log_type_);

    if (rt_window_ == 0.0)
    {
      // calculate RT window based on other parameters:
      rt_window_ = 4 * peak_width;
      OPENMS_LOG_INFO << "RT window size calculated as " << rt_window_
                      << " seconds." << endl;
    }

    //-------------------------------------------------------------
    // load input
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Loading targets and creating assay library..." << endl;
    readTargets_(id);

    OPENMS_LOG_INFO << "Loading input LC-MS data..." << endl;
    MzMLFile mzml;
    mzml.setLogType(log_type_);
    mzml.getOptions().addMSLevel(1);
    mzml.load(in, ms_data_);

    // initialize algorithm classes needed later:
    Param params = feat_finder_.getParameters();
    params.setValue("stop_report_after_feature", -1); // return all features
    params.setValue("Scores:use_rt_score", "false"); // RT may not be reliable
    params.setValue("write_convex_hull", "true");
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

    //-------------------------------------------------------------
    // run feature detection
    //-------------------------------------------------------------
    keep_library_ = !lib_out.empty();
    keep_chromatograms_ = !chrom_out.empty();

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

    // features.setProteinIdentifications(proteins);
    features.ensureUniqueId();
    addDataProcessing_(features,
                       getProcessingInfo_(DataProcessing::QUANTITATION));

    // sort features:
    sort(features.begin(), features.end(), feature_compare_);

    if (!candidates_out.empty()) // store feature candidates
    {
      FeatureXMLFile().store(candidates_out, features);
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
      for (vector<FeatureGroup>::iterator group_it = overlap_groups.begin();
           group_it != overlap_groups.end(); ++group_it)
      {
        if (group_it->size() > 1)
        {
          n_overlap_groups++;
          n_overlap_features += group_it->size();
          resolveOverlappingFeatures_(*group_it, feature_bounds);
        }
      }
      features.erase(remove_if(features.begin(), features.end(),
                             feature_filter_), features.end());
      OPENMS_LOG_INFO << features.size()
               << " features left after resolving overlaps (involving "
               << n_overlap_features << " features in " << n_overlap_groups
               << " groups)." << endl;
    }
    Size n_shared = addTargetAnnotations_(features);

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

    OPENMS_LOG_INFO << "Writing final results..." << endl;
    FeatureXMLFile().store(out, features);

    if (!trafo_out.empty()) // write expected vs. observed retention times
    {
      TransformationDescription trafo;
      TransformationDescription::DataPoints points;
      for (FeatureMap::ConstIterator it = features.begin();
           it != features.end(); ++it)
      {
        TransformationDescription::DataPoint point;
        point.first = it->getMetaValue("expected_rt");
        point.second = it->getRT();
        point.note = it->getMetaValue("PeptideRef");
        points.push_back(point);
      }
      trafo.setDataPoints(points);
      TransformationXMLFile().store(trafo_out, trafo);
    }

    //-------------------------------------------------------------
    // statistics
    //-------------------------------------------------------------

    Size n_missing = features.getUnassignedPeptideIdentifications().size();
    OPENMS_LOG_INFO << "\nSummary statistics:\n"
             << library_.getCompounds().size() << " targets specified\n"
             << features.size() << " features found\n"
             << n_shared << " features with multiple target annotations\n"
             << n_missing << " targets without features";
    const Size n_examples = 5;
    if (n_missing)
    {
      OPENMS_LOG_INFO << ":";
      for (Size i = 0;
           ((i < features.getUnassignedPeptideIdentifications().size()) &&
            (i < n_examples)); ++i)
      {
        const PeptideIdentification& id =
          features.getUnassignedPeptideIdentifications()[i];
        const TargetedExperiment::Compound& compound =
          library_.getCompoundByRef(id.getMetaValue("PeptideRef"));
        OPENMS_LOG_INFO << "\n- " << prettyPrintCompound_(compound);
      }
      if (n_missing > n_examples)
      {
        OPENMS_LOG_INFO << "\n- ... (" << n_missing - n_examples << " more)";
      }
    }
    OPENMS_LOG_INFO << "\n" << endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderMetaboIdent tool;
  return tool.main(argc, argv);
}

/// @endcond
