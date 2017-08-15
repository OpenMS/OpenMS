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
   @page TOPP_FeatureFinderMetaboIdent FeatureFinderMetaboIdent

   @brief Detects features in MS1 data based on metabolite identifications.

   <CENTER>
     <table>
       <tr>
         <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
         <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ FeatureFinderMetaboIdent \f$ \longrightarrow \f$</td>
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

   This tool detects quantitative features in MS1 data based on information from metabolite identifications.
   It uses algorithms for targeted data analysis from the OpenSWATH pipeline.

   The aim is to detect features that enable the quantification of (ideally) all entries in the identification input.

   <B>The command line parameters of this tool are:</B>
   @verbinclude TOPP_FeatureFinderMetaboIdent.cli
   <B>INI file documentation of this tool:</B>
   @htmlinclude TOPP_FeatureFinderMetaboIdent.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPFeatureFinderMetaboIdent :
  public TOPPBase
{
public:
  TOPPFeatureFinderMetaboIdent() :
    TOPPBase("FeatureFinderMetaboIdent", "Detects features in MS1 data based on metabolite identifications.")
  {
    rt_term_.setCVIdentifierRef("MS");
    rt_term_.setAccession("MS:1000896");
    rt_term_.setName("normalized retention time");
  }

protected:

  void registerOptionsAndFlags_()
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

  // predicate for filtering features by overall quality:
  struct FeatureFilterQuality
  {
    bool operator()(const Feature& feature)
    {
      return feature.metaValueExists("FFMetId_remove");
    }
  } feature_filter_;

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
  TargetedExperiment library_; // accumulated assays for targets
  bool keep_library_; // keep assay data for output?
  CVTerm rt_term_; // controlled vocabulary term for reference RT
  double rt_window_; // RT window width
  double mz_window_; // m/z window width
  bool mz_window_ppm_; // m/z window width is given in PPM (not Da)?
  double isotope_pmin_; // min. isotope probability for MS1 assay
  Size n_isotopes_; // number of isotopes for MS1 assay
  map<String, double> isotope_probs_; // isotope probabilities of transitions
  map<String, double> target_rts_; // RTs of targets (assays)
  MRMFeatureFinderScoring feat_finder_; // OpenSWATH feature finder
  ProgressLogger prog_log_;


  // read input with information about targets:
  void readTargets_(const String& in_path)
  {
    const string header =
      "Name\tFormula\tMass\tCharge\tRT\tRT_range\tIso_distrib";
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
    String name, formula, temp_charge, temp_rt, temp_range, temp_iso;
    double mass;
    Size line_count = 1;
    set<String> names;
    while (getline(source, line))
    {
      line_count++;
      if (line[0] == '#') continue; // skip comments
      stringstream ss(line);
      ss >> name;
      if (name.empty())
      {
        LOG_ERROR << "Error: Empty name field in input line " << line_count
                  << " - skipping this line." << endl;
        continue;
      }
      if (!names.insert(name).second) // @TODO: is this check necessary?
      {
        LOG_ERROR << "Error: Duplicate name '" << name << "' in input line "
                  << line_count << " - skipping this line." << endl;
        continue;
      }
      ss >> formula >> mass >> temp_charge >> temp_rt >> temp_range >> temp_iso;
      vector<Int> charges = ListUtils::create<Int>(temp_charge);
      vector<double> rts = ListUtils::create<double>(temp_rt);
      vector<double> rt_ranges = ListUtils::create<double>(temp_range);
      vector<double> iso_distrib = ListUtils::create<double>(temp_iso);
      addTargetToLibrary_(name, formula, mass, charges, rts, rt_ranges,
                          iso_distrib);
    }
  }


  void addTargetToLibrary_(const String& name, const String& formula,
                           double mass, const vector<Int>& charges,
                           const vector<double>& rts,
                           vector<double> rt_ranges,
                           const vector<double>& iso_distrib)
  {
    if ((mass <= 0) && formula.empty())
    {
      LOG_ERROR << "Error: No mass or sum formula given for target '" << name
                << "' - skipping this target." << endl;
      return;
    }
    if (rts.empty())
    {
      LOG_ERROR << "Error: No retention time (RT) given for target '" << name
                << "' - skipping this target." << endl;
      return;
    }
    // @TODO: detect entries with same RT and m/z ("collisions")
    TargetedExperiment::Compound target;
    target.molecular_formula = formula;
    EmpiricalFormula emp_formula(formula);
    bool mass_given = (mass > 0);
    if (!mass_given) mass = emp_formula.getMonoWeight();
    target.theoretical_mass = mass;
    String target_id = name + "_m" + String(float(mass));

    // get isotope distribution for target:
    IsotopeDistribution iso_dist;
    Size n_isotopes = (isotope_pmin_ > 0.0) ? 10 : n_isotopes_;
    if (iso_distrib.empty() || (iso_distrib[0] == 0))
    {
      if (formula.empty())
      {
        LOG_ERROR << "Error: No sum formula given for target '" << name
                  << "'; cannot calculate isotope distribution"
                  << " - using estimation method for peptides." << endl;
        iso_dist.estimateFromPeptideWeight(mass);
      }
      else
      {
        iso_dist = emp_formula.getIsotopeDistribution(n_isotopes);
      }
    }
    else
    {
      n_isotopes = min(n_isotopes, iso_distrib.size());
      IsotopeDistribution::ContainerType probs;
      probs.reserve(n_isotopes);
      for (Size i = 0; i < n_isotopes; ++i)
      {
        probs.push_back(make_pair(i, iso_distrib[i]));
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
        LOG_ERROR << "Error: Invalid charge 0 for target '" << name
                  << "' - skipping this charge." << endl;
        continue;
      }
      target.setChargeState(*z_it);
      double mz = 0.0;
      if (!mass_given) // calculate m/z from formula
      {
        emp_formula.setCharge(*z_it);
        mz = emp_formula.getMonoWeight() / *z_it;
      }
      else
      {
        mz = (mass + *z_it * Constants::PROTON_MASS_U) / *z_it;
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


  // generate transitions for a target ion and add them to the library:
  void generateTransitions_(const String& target_id, double mz, Int charge,
                            const IsotopeDistribution& iso_dist)
  {
    // go through different isotopes:
    Size counter = 0;
    for (IsotopeDistribution::ConstIterator iso_it = iso_dist.begin();
         iso_it != iso_dist.end(); ++iso_it, ++counter)
    {
      ReactionMonitoringTransition transition;
      String annotation = "i" + String(counter + 1);
      String transition_name = target_id + "_" + annotation;

      transition.setNativeID(transition_name);
      transition.setPrecursorMZ(mz);
      transition.setProductMZ(mz + Constants::C13C12_MASSDIFF_U *
                              float(counter) / charge);
      transition.setLibraryIntensity(iso_it->second);
      // transition.setMetaValue("annotation", annotation); // ???
      transition.setCompoundRef(target_id);
      library_.addTransition(transition);
      isotope_probs_[transition_name] = iso_it->second;
    }
  }


  void addTargetRT_(TargetedExperiment::Compound& target, double rt)
  {
    rt_term_.setValue(rt);
    TargetedExperiment::RetentionTime te_rt;
    te_rt.addCVTerm(rt_term_);
    target.rts.push_back(te_rt);
  }


  void annotateFeatures_(FeatureMap& features)
  {
    for (FeatureMap::Iterator feat_it = features.begin();
         feat_it != features.end(); ++feat_it)
    {
      feat_it->setMZ(feat_it->getMetaValue("PrecursorMZ"));
      String ref = feat_it->getMetaValue("PeptideRef");
      Int z = library_.getCompoundByRef(ref).getChargeState();
      feat_it->setCharge(z);
      ensureConvexHulls_(*feat_it);
      feat_it->getPeptideIdentifications().clear();
      // annotate subordinates with theoretical isotope intensities:
      for (vector<Feature>::iterator sub_it =
             feat_it->getSubordinates().begin(); sub_it !=
             feat_it->getSubordinates().end(); ++sub_it)
      {
        String native_id = sub_it->getMetaValue("native_id");
        sub_it->setMetaValue("isotope_probability", isotope_probs_[native_id]);
      }
    }
    features.getProteinIdentifications().clear();
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


  void filterFeatures_(FeatureMap& features)
  {
    String previous_ref;
    double best_rt_dist = rt_window_;
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
      double rt_dist = rt_window_;
      if ((rt_min <= target_rt) && (rt_max >= target_rt))
      {
        if (best_rt_dist <= 0.0)
        {
          String msg = "overlapping feature candidates for assay '" + ref + "'";
          throw Exception::Precondition(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION, msg);
        }
        rt_dist = 0.0;
      }
      else if (best_rt_dist > 0.0)
      {
        rt_dist = (rt_min > target_rt) ? (rt_min - target_rt) : (target_rt -
                                                                 rt_max);
      }
      if (rt_dist < best_rt_dist) // new best candidate for this assay
      {
        best_rt_dist = rt_dist;
        if (best_it != it) best_it->setMetaValue("FFMetId_remove", "");
        best_it = it;
      }
      else // this candidate is worse than a previous one
      {
        it->setMetaValue("FFMetId_remove", "");
      }
    }
    features.erase(remove_if(features.begin(), features.end(),
                             feature_filter_), features.end());
  }


  ExitCodes main_(int, const char**)
  {
    FeatureMap features;

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String out = getStringOption_("out");
    String candidates_out = getStringOption_("candidates_out");
    String elution_model = getStringOption_("model:type");
    prog_log_.setLogType(log_type_);

    String candidates_in = getStringOption_("candidates_in");
    if (candidates_in.empty())
    {
      String in = getStringOption_("in");
      String id = getStringOption_("id");
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

      if (rt_window_ == 0.0)
      {
        // calculate RT window based on other parameters:
        rt_window_ = 4 * peak_width;
        LOG_INFO << "RT window size calculated as " << rt_window_ << " seconds."
                 << endl;
      }

      // read target IDs and create assay library:
      LOG_INFO << "Creating assay library..." << endl;
      readTargets_(id);

      //-------------------------------------------------------------
      // run feature detection
      //-------------------------------------------------------------
      keep_library_ = !lib_out.empty();
      keep_chromatograms_ = !chrom_out.empty();

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
    }

    // sort features:
    sort(features.begin(), features.end(), feature_compare_);

    if (!candidates_out.empty()) // store feature candidates
    {
      FeatureXMLFile().store(candidates_out, features);
    }

    filterFeatures_(features);
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

    LOG_INFO << "\nSummary statistics:\n"
             << library_.getCompounds().size() << " targets specified\n"
             << features.size() << " features found\n"
             << library_.getCompounds().size() - features.size()
             << " targets without features\n" << endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPFeatureFinderMetaboIdent tool;
  return tool.main(argc, argv);
}

/// @endcond
