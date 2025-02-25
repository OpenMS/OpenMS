// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Tom Lukas Lankenau $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>

#include <OpenMS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/FEATUREFINDER/GaussTraceFitter.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <QtCore/QDir>

#include <boost/math/special_functions/fpclassify.hpp> // isnan

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS
{
  FeatureFinderAlgorithmPicked::FeatureFinderAlgorithmPicked() :
    DefaultParamHandler("FeatureFinderAlgorithmPicked")
  {
    //debugging
    defaults_.setValue("write_debug", "false", "When debug mode is activated, several files with intermediate results are written to the folder 'debug' (do not use in parallel mode).");
    defaults_.setValidStrings("write_debug", {"true","false"});
    //intensity
    defaults_.setValue("intensity:bins", 10, "Number of bins per dimension (RT and m/z). The higher this value, the more local the intensity significance score is.\nThis parameter should be decreased, if the algorithm is used on small regions of a map.");
    defaults_.setMinInt("intensity:bins", 1);
    defaults_.setSectionDescription("intensity", "Settings for the calculation of a score indicating if a peak's intensity is significant in the local environment (between 0 and 1)");
    //mass trace search parameters
    defaults_.setValue("mass_trace:mz_tolerance", 0.03, "Tolerated m/z deviation of peaks belonging to the same mass trace.\nIt should be larger than the m/z resolution of the instrument.\nThis value must be smaller than that 1/charge_high!");
    defaults_.setMinFloat("mass_trace:mz_tolerance", 0.0);
    defaults_.setValue("mass_trace:min_spectra", 10, "Number of spectra that have to show a similar peak mass in a mass trace.");
    defaults_.setMinInt("mass_trace:min_spectra", 1);
    defaults_.setValue("mass_trace:max_missing", 1, "Number of consecutive spectra where a high mass deviation or missing peak is acceptable.\nThis parameter should be well below 'min_spectra'!");
    defaults_.setMinInt("mass_trace:max_missing", 0);
    defaults_.setValue("mass_trace:slope_bound", 0.1, "The maximum slope of mass trace intensities when extending from the highest peak.\nThis parameter is important to separate overlapping elution peaks.\nIt should be increased if feature elution profiles fluctuate a lot.");
    defaults_.setMinFloat("mass_trace:slope_bound", 0.0);
    defaults_.setSectionDescription("mass_trace", "Settings for the calculation of a score indicating if a peak is part of a mass trace (between 0 and 1).");
    //Isotopic pattern search parameters
    defaults_.setValue("isotopic_pattern:charge_low", 1, "Lowest charge to search for.");
    defaults_.setMinInt("isotopic_pattern:charge_low", 1);
    defaults_.setValue("isotopic_pattern:charge_high", 4, "Highest charge to search for.");
    defaults_.setMinInt("isotopic_pattern:charge_high", 1);
    defaults_.setValue("isotopic_pattern:mz_tolerance", 0.03, "Tolerated m/z deviation from the theoretical isotopic pattern.\nIt should be larger than the m/z resolution of the instrument.\nThis value must be smaller than that 1/charge_high!");
    defaults_.setMinFloat("isotopic_pattern:mz_tolerance", 0.0);
    defaults_.setValue("isotopic_pattern:intensity_percentage", 10.0, "Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity must be present.", {"advanced"});
    defaults_.setMinFloat("isotopic_pattern:intensity_percentage", 0.0);
    defaults_.setMaxFloat("isotopic_pattern:intensity_percentage", 100.0);
    defaults_.setValue("isotopic_pattern:intensity_percentage_optional", 0.1, "Isotopic peaks that contribute more than this percentage to the overall isotope pattern intensity can be missing.", {"advanced"});
    defaults_.setMinFloat("isotopic_pattern:intensity_percentage_optional", 0.0);
    defaults_.setMaxFloat("isotopic_pattern:intensity_percentage_optional", 100.0);
    defaults_.setValue("isotopic_pattern:optional_fit_improvement", 2.0, "Minimal percental improvement of isotope fit to allow leaving out an optional peak.", {"advanced"});
    defaults_.setMinFloat("isotopic_pattern:optional_fit_improvement", 0.0);
    defaults_.setMaxFloat("isotopic_pattern:optional_fit_improvement", 100.0);
    defaults_.setValue("isotopic_pattern:mass_window_width", 25.0, "Window width in Dalton for precalculation of estimated isotope distributions.", {"advanced"});
    defaults_.setMinFloat("isotopic_pattern:mass_window_width", 1.0);
    defaults_.setMaxFloat("isotopic_pattern:mass_window_width", 200.0);
    defaults_.setValue("isotopic_pattern:abundance_12C", 98.93, "Rel. abundance of the light carbon. Modify if labeled.", {"advanced"});
    defaults_.setMinFloat("isotopic_pattern:abundance_12C", 0.0);
    defaults_.setMaxFloat("isotopic_pattern:abundance_12C", 100.0);
    defaults_.setValue("isotopic_pattern:abundance_14N", 99.632, "Rel. abundance of the light nitrogen. Modify if labeled.", {"advanced"});
    defaults_.setMinFloat("isotopic_pattern:abundance_14N", 0.0);
    defaults_.setMaxFloat("isotopic_pattern:abundance_14N", 100.0);

    defaults_.setSectionDescription("isotopic_pattern", "Settings for the calculation of a score indicating if a peak is part of a isotopic pattern (between 0 and 1).");
    //Seed settings
    defaults_.setValue("seed:min_score", 0.8, "Minimum seed score a peak has to reach to be used as seed.\nThe seed score is the geometric mean of intensity score, mass trace score and isotope pattern score.\nIf your features show a large deviation from the averagene isotope distribution or from an gaussian elution profile, lower this score.");
    defaults_.setMinFloat("seed:min_score", 0.0);
    defaults_.setMaxFloat("seed:min_score", 1.0);
    defaults_.setSectionDescription("seed", "Settings that determine which peaks are considered a seed");
    //Fitting settings
    defaults_.setValue("fit:max_iterations", 500, "Maximum number of iterations of the fit.", {"advanced"});
    defaults_.setMinInt("fit:max_iterations", 1);
    defaults_.setSectionDescription("fit", "Settings for the model fitting");
    //Feature settings
    defaults_.setValue("feature:min_score", 0.7, "Feature score threshold for a feature to be reported.\nThe feature score is the geometric mean of the average relative deviation and the correlation between the model and the observed peaks.");
    defaults_.setMinFloat("feature:min_score", 0.0);
    defaults_.setMaxFloat("feature:min_score", 1.0);
    defaults_.setValue("feature:min_isotope_fit", 0.8, "Minimum isotope fit of the feature before model fitting.", {"advanced"});
    defaults_.setMinFloat("feature:min_isotope_fit", 0.0);
    defaults_.setMaxFloat("feature:min_isotope_fit", 1.0);
    defaults_.setValue("feature:min_trace_score", 0.5, "Trace score threshold.\nTraces below this threshold are removed after the model fitting.\nThis parameter is important for features that overlap in m/z dimension.", {"advanced"});
    defaults_.setMinFloat("feature:min_trace_score", 0.0);
    defaults_.setMaxFloat("feature:min_trace_score", 1.0);
    defaults_.setValue("feature:min_rt_span", 0.333, "Minimum RT span in relation to extended area that has to remain after model fitting.", {"advanced"});
    defaults_.setMinFloat("feature:min_rt_span", 0.0);
    defaults_.setMaxFloat("feature:min_rt_span", 1.0);
    defaults_.setValue("feature:max_rt_span", 2.5, "Maximum RT span in relation to extended area that the model is allowed to have.", {"advanced"});
    defaults_.setMinFloat("feature:max_rt_span", 0.5);
    defaults_.setValue("feature:rt_shape", "symmetric", "Choose model used for RT profile fitting. If set to symmetric a gauss shape is used, in case of asymmetric an EGH shape is used.", {"advanced"});
    defaults_.setValidStrings("feature:rt_shape", {"symmetric","asymmetric"});
    defaults_.setValue("feature:max_intersection", 0.35, "Maximum allowed intersection of features.", {"advanced"});
    defaults_.setMinFloat("feature:max_intersection", 0.0);
    defaults_.setMaxFloat("feature:max_intersection", 1.0);
    defaults_.setValue("feature:reported_mz", "monoisotopic", "The mass type that is reported for features.\n'maximum' returns the m/z value of the highest mass trace.\n'average' returns the intensity-weighted average m/z value of all contained peaks.\n'monoisotopic' returns the monoisotopic m/z value derived from the fitted isotope model.");
    defaults_.setValidStrings("feature:reported_mz", {"maximum","average","monoisotopic"});
    defaults_.setSectionDescription("feature", "Settings for the features (intensity, quality assessment, ...)");
    //user-specified seed settings
    defaults_.setValue("user-seed:rt_tolerance", 5.0, "Allowed RT deviation of seeds from the user-specified seed position.");
    defaults_.setMinFloat("user-seed:rt_tolerance", 0.0);
    defaults_.setValue("user-seed:mz_tolerance", 1.1, "Allowed m/z deviation of seeds from the user-specified seed position.");
    defaults_.setMinFloat("user-seed:mz_tolerance", 0.0);
    defaults_.setValue("user-seed:min_score", 0.5, "Overwrites 'seed:min_score' for user-specified seeds. The cutoff is typically a bit lower in this case.");
    defaults_.setMinFloat("user-seed:min_score", 0.0);
    defaults_.setMaxFloat("user-seed:min_score", 1.0);
    defaults_.setSectionDescription("user-seed", "Settings for user-specified seeds.");
    //advanced/debugging settings
    defaults_.setValue("advanced:pseudo_rt_shift", 500.0, "Pseudo RT shift used when .", {"advanced"});
    defaults_.setMinFloat("advanced:pseudo_rt_shift", 1.0);
    this->defaultsToParam_();
  }

  void FeatureFinderAlgorithmPicked::setSeeds(const FeatureMap& seeds)
  {
    seeds_ = seeds;
  }

  void FeatureFinderAlgorithmPicked::setData(const MSExperiment& map, FeatureMap& features)
  {
    map_ = map;
    features_ = &features;
  }

  void FeatureFinderAlgorithmPicked::run()
  {
    //-------------------------------------------------------------------------
    // General initialization
    //---------------------------------------------------------------------------

    //quality estimation
    double min_feature_score = param_.getValue("feature:min_score");
    // charges to look at
    SignedSize charge_low = (Int)param_.getValue("isotopic_pattern:charge_low");
    SignedSize charge_high = (Int)param_.getValue("isotopic_pattern:charge_high");
    //fitting settings
    UInt max_iterations = param_.getValue("fit:max_iterations");

    Size max_isotopes = 20;

    //check if non-natural isotopic abundances are set. If so modify
    double abundance_12C = param_.getValue("isotopic_pattern:abundance_12C");
    double abundance_14N = param_.getValue("isotopic_pattern:abundance_14N");

    const Element* carbon_const = ElementDB::getInstance()->getElement("Carbon");
    Element* carbon = const_cast<Element*>(carbon_const);

    if (param_.getValue("isotopic_pattern:abundance_12C") != defaults_.getValue("isotopic_pattern:abundance_12C"))
    {
      max_isotopes += 1000; // Why?
      IsotopeDistribution isotopes;
      isotopes.insert(12, abundance_12C / 100.0);
      isotopes.insert(13, 1.0 - (abundance_12C / 100.0));
      carbon->setIsotopeDistribution(isotopes);
    }

    const Element* nitrogen_const = ElementDB::getInstance()->getElement("Nitrogen");
    Element* nitrogen = const_cast<Element*>(nitrogen_const);

    if (param_.getValue("isotopic_pattern:abundance_14N") != defaults_.getValue("isotopic_pattern:abundance_14N"))
    {
      max_isotopes += 1000; // Why?
      IsotopeDistribution isotopes;
      isotopes.insert(14, abundance_14N / 100.0);
      isotopes.insert(15, 1.0 - (abundance_14N / 100.0));
      nitrogen->setIsotopeDistribution(isotopes);
    }

    // initialize trace fitter parameters here to avoid
    // https://github.com/OpenMS/OpenMS/issues/147
    Param trace_fitter_params;
    trace_fitter_params.setValue("max_iteration", max_iterations);

    //flag for user-specified seed mode
    bool user_seeds = (!seeds_.empty());
    if (user_seeds)
    {
      seeds_.sortByMZ();
    }
    double user_rt_tol = param_.getValue("user-seed:rt_tolerance");
    double user_mz_tol = param_.getValue("user-seed:mz_tolerance");
    double user_seed_score = param_.getValue("user-seed:min_score");

    //reserve space for calculated scores
    UInt charge_count = charge_high - charge_low + 1;
    for (auto& s : map_)
    {
      Size scan_size = s.size();
      s.getFloatDataArrays().resize(3 + 2 * charge_count);
      s.getFloatDataArrays()[0].setName("trace_score");
      s.getFloatDataArrays()[0].assign(scan_size, 0.0);
      s.getFloatDataArrays()[1].setName("intensity_score");
      s.getFloatDataArrays()[1].assign(scan_size, 0.0);
      s.getFloatDataArrays()[2].setName("local_max");
      s.getFloatDataArrays()[2].assign(scan_size, 0.0);
      //create isotope pattern score arrays
      UInt charge = charge_low;
      for (Size i = 3; i < 3 + charge_count; ++i)
      {
        s.getFloatDataArrays()[i].setName(String("pattern_score_") + charge);
        s.getFloatDataArrays()[i].assign(scan_size, 0.0);
        ++charge;
      }
      //create overall score arrays
      charge = charge_low;
      for (Size i = 3 + charge_count; i < 3 + 2 * charge_count; ++i)
      {
        s.getFloatDataArrays()[i].setName(String("overall_score_") + charge);
        s.getFloatDataArrays()[i].assign(scan_size, 0.0);
        ++charge;
      }
    }

    debug_ = param_.getValue("write_debug").toBool();
    //clean up / create folders for debug information
    if (debug_)
    {
      QDir dir(".");
      dir.mkpath("debug/features");
      log_.open("debug/log.txt");
    }

    //---------------------------------------------------------------------------
    // Step 1:
    // Precalculate intensity scores for peaks
    //---------------------------------------------------------------------------
    if (debug_) log_ << "Precalculating intensity thresholds ..." << std::endl;
    //new scope to make local variables disappear
    {
      startProgress(0, intensity_bins_ * intensity_bins_, "Precalculating intensity scores");
      double rt_start = map_.getMinRT();
      double mz_start = map_.getMinMZ();
      intensity_rt_step_ = (map_.getMaxRT() - rt_start) / (double)intensity_bins_;
      intensity_mz_step_ = (map_.getMaxMZ() - mz_start) / (double)intensity_bins_;
      intensity_thresholds_.resize(intensity_bins_);
      for (Size rt = 0; rt < intensity_bins_; ++rt)
      {
        intensity_thresholds_[rt].resize(intensity_bins_);
        double min_rt = rt_start + rt * intensity_rt_step_;
        double max_rt = rt_start + (rt + 1) * intensity_rt_step_;
        std::vector<double> tmp;
        for (Size mz = 0; mz < intensity_bins_; ++mz)
        {
          setProgress(rt * intensity_bins_ + mz);
          double min_mz = mz_start + mz * intensity_mz_step_;
          double max_mz = mz_start + (mz + 1) * intensity_mz_step_;
          //std::cout << "rt range: " << min_rt << " - " << max_rt << std::endl;
          //std::cout << "mz range: " << min_mz << " - " << max_mz << std::endl;
          tmp.clear();
          for (MapType::ConstAreaIterator it = map_.areaBeginConst(min_rt, max_rt, min_mz, max_mz); it != map_.areaEndConst(); ++it)
          {
            tmp.push_back(it->getIntensity());
          }
          //init vector
          intensity_thresholds_[rt][mz].assign(21, 0.0);
          //store quantiles (20)
          if (!tmp.empty())
          {
            std::sort(tmp.begin(), tmp.end());
            for (Size i = 0; i < 21; ++i)
            {
              Size index = (Size)std::floor(0.05 * i * (tmp.size() - 1));
              intensity_thresholds_[rt][mz][i] = tmp[index];
            }
          }
        }
      }

      //store intensity score in PeakInfo
      for (Size s = 0; s < map_.size(); ++s)
      {
        for (Size p = 0; p < map_[s].size(); ++p)
        {
          map_[s].getFloatDataArrays()[1][p] = intensityScore_(s, p);
        }
      }
      endProgress();
    }

    //---------------------------------------------------------------------------
    // Step 2:
    // Precalculate mass trace scores and local trace maximum for each peak
    //---------------------------------------------------------------------------
    //new scope to make local variables disappear
    {
      Size end_iteration = map_.size() - std::min((Size)min_spectra_, map_.size());
      startProgress(min_spectra_, end_iteration, "Precalculating mass trace scores");
      // skip first and last scans since we cannot extend the mass traces there
      for (Size s = min_spectra_; s < end_iteration; ++s)
      {
        setProgress(s);
        SpectrumType& spectrum = map_[s];
        //iterate over all peaks of the scan
        for (Size p = 0; p < spectrum.size(); ++p)
        {
          double trace_score = 0.0;

          double pos = spectrum[p].getMZ();
          float intensity = spectrum[p].getIntensity();

          //if(debug_) log_ << std::endl << "Peak: " << pos << std::endl;
          bool is_max_peak = true; // checking the maximum intensity peaks -> use them later as feature seeds.
          for (Size i = 1; i <= min_spectra_; ++i)
          {
            SpectrumType& next_spectrum = map_[s + i];
            if (!next_spectrum.empty()) // There are peaks in the spectrum
            {
              Size spec_index = next_spectrum.findNearest(pos);
              double position_score = positionScore_(pos, next_spectrum[spec_index].getMZ(), trace_tolerance_);
              if (position_score > 0 && next_spectrum[spec_index].getIntensity() > intensity) is_max_peak = false;
              trace_score += position_score;
            }
          }
          for (Size i = 1; i <= min_spectra_; ++i)
          {
            SpectrumType& next_spectrum = map_[s - i];
            if (!next_spectrum.empty()) // There are peaks in the spectrum
            {
              Size spec_index = next_spectrum.findNearest(pos);
              double position_score = positionScore_(pos, next_spectrum[spec_index].getMZ(), trace_tolerance_);
              if (position_score > 0 && next_spectrum[spec_index].getIntensity() > intensity)
              {
                is_max_peak = false;
              }
              trace_score += position_score;
            }
          }
          //Calculate a consensus score out of the scores calculated before
          trace_score /= 2 * min_spectra_;

          //store final score for later use
          spectrum.getFloatDataArrays()[0][p] = trace_score;
          spectrum.getFloatDataArrays()[2][p] = is_max_peak;
        }
      }
      endProgress();
    }

    //---------------------------------------------------------------------------
    // Step 2.5:
    // Precalculate isotope distributions for interesting mass ranges
    //---------------------------------------------------------------------------
    //new scope to make local variables disappear
    {
      double max_mass = map_.getMaxMZ() * charge_high;
      Size num_isotopes = std::ceil(max_mass / mass_window_width_) + 1;
      startProgress(0, num_isotopes, "Precalculating isotope distributions");

      //reserve enough space
      isotope_distributions_.resize(num_isotopes);

      //calculate distribution if necessary
      for (Size index = 0; index < num_isotopes; ++index)
      {
        //if(debug_) log_ << "Calculating iso dist for mass: " << 0.5*mass_window_width_ + index * mass_window_width_ << std::endl;
        CoarseIsotopePatternGenerator solver(max_isotopes);
        auto d = solver.estimateFromPeptideWeight(0.5 * mass_window_width_ + index * mass_window_width_);
        //trim left and right. And store the number of isotopes on the left, to reconstruct the monoisotopic peak
        Size size_before = d.size();
        d.trimLeft(intensity_percentage_optional_);
        isotope_distributions_[index].trimmed_left = size_before - d.size();
        d.trimRight(intensity_percentage_optional_);

        for (auto& peak : d)
        {
          isotope_distributions_[index].intensity.push_back(peak.getIntensity());
          //if(debug_) log_ << " - " << it->second << std::endl;
        }

        //determine the number of optional peaks at the beginning/end
        Size begin = 0;
        Size end = 0;
        bool is_begin = true;
        bool is_end = false;

        for (double i : isotope_distributions_[index].intensity)
        {
          if (i < intensity_percentage_)
          {
            if (!is_end && !is_begin)
            {
              is_end = true;
            }
            if (is_begin)
            {
              ++begin;
            }
            else if (is_end)
            {
              ++end;
            }
          }
          else if (is_begin)
          {
            is_begin = false;
          }
        }
        isotope_distributions_[index].optional_begin = begin;
        isotope_distributions_[index].optional_end = end;
        //scale the distribution to a maximum of 1
        double max = 0.0;
        for (double i : isotope_distributions_[index].intensity)
        {
          if (i > max)
          {
            max = i;
          }
        }
        isotope_distributions_[index].max = max;
        for (double& i : isotope_distributions_[index].intensity)
        {
          i /= max;
        }

        //if(debug_) log_ << " - optional begin/end:" << begin << " / " << end << std::endl;
      }

      endProgress();
    }

    //-------------------------------------------------------------------------
    // Step 3:
    // Charge loop (create seeds and features for each charge separately)
    //-------------------------------------------------------------------------
    Int plot_nr_global = -1;   // counter for the number of plots (debug info)
    Int feature_nr_global = 0; // counter for the number of features (debug info)
    for (SignedSize c = charge_low; c <= charge_high; ++c)
    {
      UInt meta_index_isotope = 3 + c - charge_low;
      UInt meta_index_overall = 3 + charge_count + c - charge_low;

      Size feature_candidates = 0;
      std::vector<Seed> seeds;

      //-----------------------------------------------------------
      // Step 3.1: Precalculate IsotopePattern score
      //-----------------------------------------------------------
      startProgress(0, map_.size(), String("Calculating isotope pattern scores for charge ") + String(c));
      for (Size s = 0; s < map_.size(); ++s)
      {
        setProgress(s);
        const SpectrumType& spectrum = map_[s];
        for (Size p = 0; p < spectrum.size(); ++p)
        {
          double mz = spectrum[p].getMZ();

          //get isotope distribution for this mass
          const TheoreticalIsotopePattern& isotopes = getIsotopeDistribution_(mz * c);
          //determine highest peak in isotope distribution
          Size max_isotope = std::max_element(isotopes.intensity.begin(), isotopes.intensity.end()) - isotopes.intensity.begin();
          //Look up expected isotopic peaks (in the current spectrum or adjacent spectra)
          Size peak_index = spectrum.findNearest(mz - ((double)(isotopes.size() + 1) / c));
          IsotopePattern pattern(isotopes.size());

          for (Size i = 0; i < isotopes.size(); ++i)
          {
            double isotope_pos = mz + ((double)i - max_isotope) / c;
            findIsotope_(isotope_pos, s, pattern, i, peak_index);
          }

          double pattern_score = isotopeScore_(isotopes, pattern, true);

          //update pattern scores of all contained peaks (if necessary)
          if (pattern_score > 0.0)
          {
            for (Size i = 0; i < pattern.peak.size(); ++i)
            {
              if (pattern.peak[i] >= 0 && pattern_score > map_[pattern.spectrum[i]].getFloatDataArrays()[meta_index_isotope][pattern.peak[i]])
              {
                map_[pattern.spectrum[i]].getFloatDataArrays()[meta_index_isotope][pattern.peak[i]] = pattern_score;
              }
            }
          }
        }
      }
      endProgress();
      //-----------------------------------------------------------
      // Step 3.2:
      // Find seeds for this charge
      //-----------------------------------------------------------
      Size end_of_iteration = map_.size() - std::min((Size)min_spectra_, map_.size());
      startProgress(min_spectra_, end_of_iteration, String("Finding seeds for charge ") + String(c));

      double min_seed_score = param_.getValue("seed:min_score");
      //do nothing for the first few and last few spectra as the scans required to search for traces are missing
      for (Size s = min_spectra_; s < end_of_iteration; ++s)
      {
        setProgress(s);

        //iterate over peaks
        for (Size p = 0; p < map_[s].size(); ++p)
        {
          FloatDataArrays& meta = map_[s].getFloatDataArrays();
          double overall_score = std::pow(meta[0][p] * meta[1][p] * meta[meta_index_isotope][p], 1.0f / 3.0f);
          meta[meta_index_overall][p] = overall_score;

          //add seed to vector if certain conditions are fulfilled
          if (meta[2][p] != 0.0) // local maximum of mass trace is prerequisite for all features
          {
            //automatic seeds: overall score greater than the min seed score
            if (!user_seeds && overall_score >= min_seed_score)
            {
              Seed seed;
              seed.spectrum = s;
              seed.peak = p;
              seed.intensity = map_[s][p].getIntensity();
              seeds.push_back(seed);
            }
            //user-specified seeds: overall score greater than USER min seed score
            else if (user_seeds && overall_score >= user_seed_score)
            {
              //only consider seeds, if they are near a user-specified seed
              Feature tmp;
              tmp.setMZ(map_[s][p].getMZ() - user_mz_tol);
              for (FeatureMap::const_iterator it = std::lower_bound(seeds_.begin(), seeds_.end(), tmp, Feature::MZLess()); it < seeds_.end(); ++it)
              {
                if (it->getMZ() > map_[s][p].getMZ() + user_mz_tol)
                {
                  break;
                }
                if (fabs(it->getMZ() - map_[s][p].getMZ()) < user_mz_tol && fabs(it->getRT() - map_[s].getRT()) < user_rt_tol)
                {
                  Seed seed;
                  seed.spectrum = s;
                  seed.peak = p;
                  seed.intensity = map_[s][p].getIntensity();
                  seeds.push_back(seed);
                  break;
                }
              }
            }
          }
        }
      }
      //sort seeds according to intensity
      std::sort(seeds.rbegin(), seeds.rend());
      //create and store seeds map and selected peak map
      if (debug_)
      {
        //seeds
        FeatureMap seed_map;
        seed_map.reserve(seeds.size());
        for (auto& seed : seeds)
        {
          Size spectrum = seed.spectrum;
          Size peak = seed.peak;
          const FloatDataArrays& meta = map_[spectrum].getFloatDataArrays();
          Feature tmp;
          tmp.setIntensity(seed.intensity);
          tmp.setOverallQuality(meta[meta_index_overall][peak]);
          tmp.setRT(map_[spectrum].getRT());
          tmp.setMZ(map_[spectrum][peak].getMZ());
          tmp.setMetaValue("intensity_score", meta[1][peak]);
          tmp.setMetaValue("pattern_score", meta[meta_index_isotope][peak]);
          tmp.setMetaValue("trace_score", meta[0][peak]);
          seed_map.push_back(tmp);
        }
        FileHandler().storeFeatures(String("debug/seeds_") + String(c) + ".featureXML", seed_map);
      }

      endProgress();
      std::cout << "Found " << seeds.size() << " seeds for charge " << c << "." << std::endl;

      //------------------------------------------------------------------
      // Step 3.3:
      // Extension of seeds
      //------------------------------------------------------------------

      // We do not want to store features whose seeds lie within other
      // features with higher intensity. We thus store this information in
      // the map seeds_in_features which contains for each seed i a vector
      // of other seeds that are contained in the corresponding feature i.
      //
      // The features are stored in an temporary feature map until it is
      // decided whether they are contained within a seed of higher
      // intensity.
      std::map<Size, std::vector<Size>> seeds_in_features;
      typedef std::map<Size, Feature> FeatureMapType;
      FeatureMapType tmp_feature_map;
      int gl_progress = 0;
      startProgress(0, seeds.size(), String("Extending seeds for charge ") + String(c));

#pragma omp parallel for
      for (SignedSize i = 0; i < (SignedSize)seeds.size(); ++i)
      {
        //------------------------------------------------------------------
        // Step 3.3.1:
        // Extend all mass traces
        //------------------------------------------------------------------

        const SpectrumType& spectrum = map_[seeds[i].spectrum];
        const PeakType& peak = spectrum[seeds[i].peak];

        IF_MASTERTHREAD
        {
          setProgress(gl_progress++);

          if (debug_)
          {
            log_ << std::endl << "Seed " << i << ":" << std::endl;
            //If the intensity is zero this seed is already uses in another feature
            log_ << " - Int: " << peak.getIntensity() << std::endl;
            log_ << " - RT: " << spectrum.getRT() << std::endl;
            log_ << " - MZ: " << peak.getMZ() << std::endl;
          }
        }

        //----------------------------------------------------------------
        //Find best fitting isotope pattern for this charge (using averagine)
        IsotopePattern best_pattern(0);
        double isotope_fit_quality = findBestIsotopeFit_(seeds[i], c, best_pattern);

        if (isotope_fit_quality < min_isotope_fit_)
        {
          abort_(seeds[i], "Could not find good enough isotope pattern containing the seed");
          continue;
        }
        //extend the convex hull in RT dimension (starting from the trace peaks)
        MassTraces traces;
        traces.reserve(best_pattern.peak.size());
        extendMassTraces_(best_pattern, traces, meta_index_overall);

        //check if the traces are still valid
        double seed_mz = map_[seeds[i].spectrum][seeds[i].peak].getMZ();

        if (!traces.isValid(seed_mz, trace_tolerance_))
        {
          abort_(seeds[i], "Could not extend seed");
          continue;
        }

        //------------------------------------------------------------------
        // Step 3.3.2:
        // Gauss/EGH fit (first fit to find the feature boundaries)
        //------------------------------------------------------------------
        Int plot_nr = -1;


#pragma omp critical(FeatureFinderAlgorithmPicked_PLOTNR)
        {
          plot_nr = ++plot_nr_global;
        }

        //------------------------------------------------------------------

        //TODO try fit with baseline term once more
        //baseline estimate
        traces.updateBaseline();
        traces.baseline = 0.75 * traces.baseline;

        traces[traces.max_trace].updateMaximum();

        //choose fitter
        double egh_tau = 0.0;
        std::shared_ptr<TraceFitter> fitter = chooseTraceFitter_(egh_tau);

        fitter->setParameters(trace_fitter_params);
        fitter->fit(traces);

#if 0
        TraceFitter<PeakType>* alt_fitter = new GaussTraceFitter<PeakType>();
        Param alt_p;
        alt_p.setValue("max_iteration", max_iterations);

        alt_fitter->setParameters(alt_p);
        alt_fitter->fit(traces);

        OPENMS_LOG_DEBUG << "EGH:   " << fitter->getCenter() << " " << fitter->getHeight() << std::endl;
        OPENMS_LOG_DEBUG << "GAUSS: " << alt_fitter->getCenter() << " " << alt_fitter->getHeight() << std::endl;
#endif
        // what should come out
        // left "sigma"
        // right "sigma"
        // x0 .. "center" position of RT fit
        // height .. "height" of RT fit

        //------------------------------------------------------------------

        //------------------------------------------------------------------
        // Step 3.3.3:
        // Crop feature according to RT fit (2.5*sigma) and remove badly fitting traces
        //------------------------------------------------------------------
        MassTraces new_traces;
        cropFeature_(fitter, traces, new_traces);

        //------------------------------------------------------------------
        // Step 3.3.4:
        // Check if feature is ok
        //------------------------------------------------------------------
        String error_msg = "";

        double fit_score = 0.0;
        double correlation = 0.0;
        double final_score = 0.0;

        int number_of_datapoints = 0;

        bool feature_ok = checkFeatureQuality_(fitter, new_traces, seed_mz, min_feature_score, error_msg, fit_score, correlation, final_score);

        {
          //write debug output of feature
          if (debug_)
          {
#pragma omp critical(FeatureFinderAlgorithmPicked_DEBUG)
            writeFeatureDebugInfo_(fitter, traces, new_traces, feature_ok, error_msg, final_score, plot_nr, peak);
          }
        }


        //validity output
        if (!feature_ok)
        {
          abort_(seeds[i], error_msg);
          continue;
        }
        traces = new_traces;

        //------------------------------------------------------------------
        // Step 3.3.5:
        // Feature creation
        //------------------------------------------------------------------
        Feature f;
        //set label
        f.setMetaValue(3, plot_nr);
        f.setCharge(c);
        f.setOverallQuality(final_score);
        f.setMetaValue("score_fit", fit_score);
        f.setMetaValue("score_correlation", correlation);
        f.setRT(fitter->getCenter());
        f.setWidth(fitter->getFWHM());

        // metavalue num_of_datapoints
        for (size_t t = 0; t < traces.size(); ++t)
        {
          number_of_datapoints += traces[t].peaks.size();
        }
        f.setMetaValue(Constants::UserParam::NUM_OF_DATAPOINTS, number_of_datapoints);

        //Extract some of the model parameters.
        if (egh_tau != 0.0)
        {
          egh_tau = (std::dynamic_pointer_cast<EGHTraceFitter>(fitter))->getTau();
          f.setMetaValue("EGH_tau", egh_tau);
          f.setMetaValue("EGH_height", (std::dynamic_pointer_cast<EGHTraceFitter>(fitter))->getHeight());
          f.setMetaValue("EGH_sigma", (std::dynamic_pointer_cast<EGHTraceFitter>(fitter))->getSigma());
        }

        //Calculate the mass of the feature: maximum, average, monoisotopic
        if (reported_mz_ == "maximum")
        {
          f.setMZ(traces[traces.getTheoreticalmaxPosition()].getAvgMZ());
        }
        else if (reported_mz_ == "average")
        {
          double total_intensity = 0.0;
          double average_mz = 0.0;
          for (Size t = 0; t < traces.size(); ++t)
          {
            for (auto& p : traces[t].peaks)
            {
              average_mz += p.second->getMZ() * p.second->getIntensity();
              total_intensity += p.second->getIntensity();
            }
          }
          average_mz /= total_intensity;
          f.setMZ(average_mz);
        }
        else if (reported_mz_ == "monoisotopic")
        {
          double mono_mz = traces[traces.getTheoreticalmaxPosition()].getAvgMZ();
          mono_mz -= (Constants::PROTON_MASS_U / c) * (traces.getTheoreticalmaxPosition() + best_pattern.theoretical_pattern.trimmed_left);
          f.setMZ(mono_mz);
        }

        // Calculate intensity based on model only
        // - the model does not include the baseline, so we ignore it here
        // - as we scaled the isotope distribution to
        f.setIntensity(fitter->getArea() / getIsotopeDistribution_(f.getMZ()).max);

        //add convex hulls of mass traces
        for (Size j = 0; j < traces.size(); ++j)
        {
          f.getConvexHulls().push_back(traces[j].getConvexhull());
        }

#pragma omp critical(FeatureFinderAlgorithmPicked_TMPFEATUREMAP)
        {
          tmp_feature_map[i] = f;
        }

        //----------------------------------------------------------------
        //Remember all seeds that lie inside the convex hull of the new feature
        DBoundingBox<2> bb = f.getConvexHull().getBoundingBox();
        for (Size j = i + 1; j < seeds.size(); ++j)
        {
          double rt = map_[seeds[j].spectrum].getRT();
          double mz = map_[seeds[j].spectrum][seeds[j].peak].getMZ();
          if (bb.encloses(rt, mz) && f.encloses(rt, mz))
          {
#pragma omp critical(FeatureFinderAlgorithmPicked_SEEDSINFEATURES)
            {
              seeds_in_features[i].push_back(j);
            }
          }
        }
      } //end of OPENMP over seeds

      // Here we have to evaluate which seeds are already contained in
      // features of seeds with higher intensities. Only if the seed is not
      // used in any feature with higher intensity, we can add it to the
      // features_ list.
      std::vector<Size> seeds_contained;
      for (auto& f : tmp_feature_map)
      {
        Size seed_nr = f.first;
        bool is_used = false;
        for (Size i : seeds_contained)
        {
          if (seed_nr == i)
          {
            is_used = true;
            break;
          }
        }
        if (!is_used)
        {
          ++feature_candidates;

          //re-set label
          f.second.setMetaValue(3, feature_nr_global);
          ++feature_nr_global;
          features_->push_back(f.second);

          std::vector<Size> curr_seed = seeds_in_features[seed_nr];
          for (Size k : curr_seed)
          {
            seeds_contained.push_back(k);
          }
        }
      }

      IF_MASTERTHREAD endProgress();
      std::cout << "Found " << feature_candidates << " feature candidates for charge " << c << "." << std::endl;
    }
    // END OPENMP

    //------------------------------------------------------------------
    //Step 4:
    //Resolve contradicting and overlapping features
    //------------------------------------------------------------------
    startProgress(0, features_->size() * features_->size(), "Resolving overlapping features");
    if (debug_) log_ << "Resolving intersecting features (" << features_->size() << " candidates)" << std::endl;
    //sort features according to m/z in order to speed up the resolution
    features_->sortByMZ();
    //precalculate BBs and maximum mz span
    std::vector<DBoundingBox<2> > bbs(features_->size());
    double max_mz_span = 0.0;

    for (Size i = 0; i < features_->size(); ++i)
    {
      bbs[i] = (*features_)[i].getConvexHull().getBoundingBox();
      if (bbs[i].height() > max_mz_span)
      {
        max_mz_span = bbs[i].height();
      }
    }

    Size removed(0);
    //intersect
    for (Size i = 0; i < features_->size(); ++i)
    {
      Feature& f1((*features_)[i]);
      for (Size j = i + 1; j < features_->size(); ++j)
      {
        setProgress(i * features_->size() + j);
        Feature& f2((*features_)[j]);
        //features that are more than 2 times the maximum m/z span apart do not overlap => abort
        if (f2.getMZ() - f1.getMZ() > 2.0 * max_mz_span)
        {
          break;
        }
        //do nothing if one of the features is already removed
        if (f1.getIntensity() == 0.0 || f2.getIntensity() == 0.0)
        {
          continue;
        }
        //do nothing if the overall convex hulls do not overlap
        if (!bbs[i].intersects(bbs[j]))
        {
          continue;
        }
        //act depending on the intersection
        double intersection = intersection_(f1, f2);

        if (intersection >= max_feature_intersection_)
        {
          ++removed;

          if (debug_)
          {
            log_ << " - Intersection (" << (i + 1) << "/" << (j + 1) << "): " << intersection << std::endl;
          }
          if (f1.getCharge() == f2.getCharge())
          {
            if (f1.getIntensity() * f1.getOverallQuality() > f2.getIntensity() * f2.getOverallQuality())
            {
              if (debug_)
              {
                log_ << "   - same charge -> removing duplicate " << (j + 1) << std::endl;
              }
              f1.getSubordinates().push_back(f2);
              f2.setIntensity(0.0);
            }
            else
            {
              if (debug_)
              {
                log_ << "   - same charge -> removing duplicate " << (i + 1) << std::endl;
              }
              f2.getSubordinates().push_back(f1);
              f1.setIntensity(0.0);
            }
          }
          else if (f2.getCharge() % f1.getCharge() == 0)
          {
            if (debug_)
            {
              log_ << "   - different charge (one is the multiple of the other) -> removing lower charge " << (i + 1) << std::endl;
            }
            f2.getSubordinates().push_back(f1);
            f1.setIntensity(0.0);
          }
          else if (f1.getCharge() % f2.getCharge() == 0)
          {
            if (debug_)
            {
              log_ << "   - different charge (one is the multiple of the other) -> removing lower charge " << (i + 1) << std::endl;
            }
            f1.getSubordinates().push_back(f2);
            f2.setIntensity(0.0);
          }
          else
          {
            if (f1.getOverallQuality() > f2.getOverallQuality())
            {
              if (debug_)
              {
                log_ << "   - different charge -> removing lower score " << (j + 1) << std::endl;
              }
              f1.getSubordinates().push_back(f2);
              f2.setIntensity(0.0);
            }
            else
            {
              if (debug_)
              {
                log_ << "   - different charge -> removing lower score " << (i + 1) << std::endl;
              }
              f2.getSubordinates().push_back(f1);
              f1.setIntensity(0.0);
            }
          }
        }
      }
    }
    OPENMS_LOG_INFO << "Removed " << removed << " overlapping features." << std::endl;
    // finally remove features with intensity 0
    FeatureMap tmp;
    tmp.reserve(features_->size());
    for (Size i = 0; i < features_->size(); ++i)
    {
      if (features_->operator[](i).getIntensity() != 0.0)
      {
        tmp.push_back(features_->operator[](i));
      }
    }
    tmp.swapFeaturesOnly(*features_);
    // sort features by intensity
    features_->sortByIntensity(true);
    endProgress();

    // Abort reasons
    OPENMS_LOG_INFO << '\n';
    OPENMS_LOG_INFO << "Info: reasons for not finalizing a feature during its construction:\n";
    for (const auto& reason : aborts_)
    {
      OPENMS_LOG_INFO << " - " << reason.first << ": " << reason.second << " times\n";
    }

    OPENMS_LOG_INFO << "\n" << features_->size() << " features found." << std::endl;

    if (debug_)
    {
      //store map of abort reasons for failed seeds
      FeatureMap abort_map;
      abort_map.reserve(abort_reasons_.size());
      Size counter = 0;
      for (std::map<Seed, String>::iterator it2 = abort_reasons_.begin(); it2 != abort_reasons_.end(); ++it2, ++counter)
      {
        Feature f;
        f.setRT(map_[it2->first.spectrum].getRT());
        f.setMZ(map_[it2->first.spectrum][it2->first.peak].getMZ());
        f.setIntensity(map_[it2->first.spectrum][it2->first.peak].getIntensity());
        f.setMetaValue("label", it2->second);
        f.setUniqueId(counter); // ID = index
        abort_map.push_back(f);
      }
      abort_map.setUniqueId();
      FileHandler().storeFeatures("debug/abort_reasons.featureXML", abort_map);

      //store input map with calculated scores (without overall score)
      for (auto& s : map_)
      {
        s.getFloatDataArrays().erase(s.getFloatDataArrays().begin() + 2);
      }
      FileHandler().storeExperiment("debug/input.mzML", map_, {FileTypes::MZML});
    }

  }

  void FeatureFinderAlgorithmPicked::run(PeakMap& input_map, FeatureMap& features, const Param& param, const FeatureMap& seeds)
  {
    // Nothing to do if there is no data
    if (input_map.empty())
    {
      features.clear(true);
      return;
    }

    // check input
    {
      // We need updated ranges => check number of peaks
      if (input_map.getSize() == 0)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "FeatureFinder needs updated ranges on input map. Aborting.");
      }

      // We need MS1 data only => check levels
      if (input_map.getMSLevels().size() != 1 || input_map.getMSLevels()[0] != 1)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "FeatureFinder can only operate on MS level 1 data. Please do not use MS/MS data. Aborting.");
      }

      //Check if the peaks are sorted according to m/z
      if (!input_map.isSorted(true))
      {
        OPENMS_LOG_WARN << "Input map is not sorted by RT and m/z! This is done now, before applying the algorithm!" << std::endl;
        input_map.sortSpectra(true);
        input_map.sortChromatograms(true);
      }

      for (Size s = 0; s < input_map.size(); ++s)
      {
        if (input_map[s].empty())
        {
          continue;
        }
        if (input_map[s][0].getMZ() < 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "FeatureFinder can only operate on spectra that contain peaks with positive m/z values. Filter the data accordingly beforehand! Aborting.");
        }
      }
    }

    // do the work
    setParameters(param);
    setData(input_map, features);
    setSeeds(seeds);
    run();

    //report RT apex spectrum index and native ID for each feature
    for (Size i = 0; i < features.size(); ++i)
    {
      //index
      Size spectrum_index = input_map.RTBegin(features[i].getRT()) - input_map.begin();
      features[i].setMetaValue("spectrum_index", spectrum_index);
      //native id
      if (spectrum_index < input_map.size())
      {
        String native_id = input_map[spectrum_index].getNativeID();
        features[i].setMetaValue("spectrum_native_id", native_id);
      }
      else
      {
        /// @todo that happens sometimes using IsotopeWaveletFeatureFinder (Rene, Marc, Andreas, Clemens)
        std::cerr << "FeatureFinderAlgorithm_impl, line=" << __LINE__ << "; FixMe this cannot be, but happens" << std::endl;
      }
    }
  }

  void FeatureFinderAlgorithmPicked::updateMembers_()
  {
    pattern_tolerance_ = param_.getValue("mass_trace:mz_tolerance");
    trace_tolerance_ = param_.getValue("isotopic_pattern:mz_tolerance");
    min_spectra_ = (UInt) std::floor((double)param_.getValue("mass_trace:min_spectra") * 0.5);
    max_missing_trace_peaks_ = param_.getValue("mass_trace:max_missing");
    slope_bound_ = param_.getValue("mass_trace:slope_bound");
    intensity_percentage_ = (double)param_.getValue("isotopic_pattern:intensity_percentage") / 100.0;
    intensity_percentage_optional_ = (double)param_.getValue("isotopic_pattern:intensity_percentage_optional") / 100.0;
    optional_fit_improvement_ = (double)param_.getValue("isotopic_pattern:optional_fit_improvement") / 100.0;
    mass_window_width_ = param_.getValue("isotopic_pattern:mass_window_width");
    intensity_bins_ =  param_.getValue("intensity:bins");
    min_isotope_fit_ = param_.getValue("feature:min_isotope_fit");
    min_trace_score_ = param_.getValue("feature:min_trace_score");
    min_rt_span_ = param_.getValue("feature:min_rt_span");
    max_rt_span_ = param_.getValue("feature:max_rt_span");
    max_feature_intersection_ = param_.getValue("feature:max_intersection");
    reported_mz_ = param_.getValue("feature:reported_mz").toString();
  }

  /// Writes the abort reason to the log file and counts occurrences for each reason
  void FeatureFinderAlgorithmPicked::abort_(const Seed& seed, const String& reason)
  {
    if (debug_)
    {
      log_ << "Abort: " << reason << std::endl;
    }
    aborts_[reason]++;
    if (debug_)
    {
      abort_reasons_[seed] = reason;
    }
  }

  double FeatureFinderAlgorithmPicked::intersection_(const Feature& f1, const Feature& f2) const
  {
    //calculate the RT range sum of feature 1
    double s1 = 0.0;
    const std::vector<ConvexHull2D>& hulls1 = f1.getConvexHulls();
    for (const auto& i : hulls1)
    {
      s1 += i.getBoundingBox().width();
    }

    //calculate the RT range sum of feature 2
    double s2 = 0.0;
    const std::vector<ConvexHull2D>& hulls2 = f2.getConvexHulls();
    for (const auto& j : hulls2)
    {
      s2 += j.getBoundingBox().width();
    }

    //calculate overlap
    double overlap = 0.0;
    for (const auto& i : hulls1)
    {
      DBoundingBox<2> bb1 = i.getBoundingBox();
      for (const auto& j : hulls2)
      {
        DBoundingBox<2> bb2 = j.getBoundingBox();
        if (bb1.intersects(bb2))
        {
          if (bb1.minPosition()[0] <= bb2.minPosition()[0] &&
              bb1.maxPosition()[0] >= bb2.maxPosition()[0]) //bb1 contains bb2
          {
            overlap += bb2.width();
          }
          else if (bb2.minPosition()[0] <= bb1.minPosition()[0] &&
                   bb2.maxPosition()[0] >= bb1.maxPosition()[0]) //bb2 contains bb1
          {
            overlap += bb1.width();
          }
          else if (bb1.minPosition()[0] <= bb2.minPosition()[0] &&
                   bb1.maxPosition()[0] <= bb2.maxPosition()[0]) //the end of bb1 overlaps with bb2
          {
            overlap += bb1.maxPosition()[0] - bb2.minPosition()[0];
          }
          else if (bb2.minPosition()[0] <= bb1.minPosition()[0] &&
                   bb2.maxPosition()[0] <= bb1.maxPosition()[0]) //the end of bb2 overlaps with bb1
          {
            overlap += bb2.maxPosition()[0] - bb1.minPosition()[0];
          }
        }
      }
    }

    return overlap / std::min(s1, s2);
  }

  const FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern& FeatureFinderAlgorithmPicked::getIsotopeDistribution_(double mass) const
  {
    //calculate index in the vector
    Size index = (Size) std::floor(mass / mass_window_width_);

    if (index >= isotope_distributions_.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "IsotopeDistribution not precalculated. Maximum allowed index is " + String(isotope_distributions_.size()), String(index));
    }

    //Return distribution
    return isotope_distributions_[index];
  }

  double FeatureFinderAlgorithmPicked::findBestIsotopeFit_(const Seed& center, UInt charge, IsotopePattern& best_pattern) const
  {
    if (debug_)
    {
      log_ << "Testing isotope patterns for charge " << charge << ": " << std::endl;
    }
    const SpectrumType& spectrum = map_[center.spectrum];
    const TheoreticalIsotopePattern& isotopes = getIsotopeDistribution_(spectrum[center.peak].getMZ() * charge);
    if (debug_)
    {
      log_ << " - Seed: " << center.peak << " (mz:" << spectrum[center.peak].getMZ() << ")" << std::endl;
    }
    //Find m/z boundaries of search space (linear search as this is local and we have the center already)
    double mass_window = (double)(isotopes.size() + 1) / (double)charge;
    if (debug_)
    {
      log_ << " - Mass window: " << mass_window << std::endl;
    }
    Size end = center.peak;
    while (end < spectrum.size() &&
           spectrum[end].getMZ() < spectrum[center.peak].getMZ() + mass_window)
    {
      ++end;
    }
    --end;

    //search begin
    SignedSize begin = center.peak;
    while (begin >= 0 &&
           spectrum[begin].getMZ() > spectrum[center.peak].getMZ() - mass_window)
    {
      --begin;
    }
    ++begin;
    if (debug_)
    {
      log_ << " - Begin: " << begin << " (mz:" << spectrum[begin].getMZ() << ")" << std::endl;
    }
    if (debug_)
    {
      log_ << " - End: " << end << " (mz:" << spectrum[end].getMZ() << ")" << std::endl;
    }
    //fit isotope distribution to peaks
    double max_score = 0.0;
    for (Size start = begin; start <= end; ++start)
    {
      //find isotope peaks for the current start peak
      Size peak_index = start;
      IsotopePattern pattern(isotopes.size());
      if (debug_)
      {
        log_ << " - Fitting at " << start << " (mz:" << spectrum[start].getMZ() << ")" << std::endl;
      }
      for (Size iso = 0; iso < isotopes.size(); ++iso)
      {
        double pos = spectrum[start].getMZ() + iso / (double)charge;
        findIsotope_(pos, center.spectrum, pattern, iso, peak_index);
      }

      //check if the seed is contained, otherwise abort
      bool seed_contained = false;
      for (Size iso = 0; iso < pattern.peak.size(); ++iso)
      {
        if (pattern.peak[iso] == (Int)center.peak && pattern.spectrum[iso] == center.spectrum)
        {
          seed_contained = true;
          break;
        }
      }
      if (!seed_contained)
      {
        if (debug_)
        {
          log_ << "   - aborting: seed is not contained!" << std::endl;
        }
        continue;
      }

      double score = isotopeScore_(isotopes, pattern, false);

      //check if the seed is still contained, otherwise abort
      seed_contained = false;
      for (Size iso = 0; iso < pattern.peak.size(); ++iso)
      {
        if (pattern.peak[iso]     == (Int)center.peak &&
            pattern.spectrum[iso] == center.spectrum)
        {
          seed_contained = true;
          break;
        }
      }
      if (!seed_contained)
      {
        if (debug_)
        {
          log_ << "   - aborting: seed was removed during isotope fit!" << std::endl;
        }
        continue;
      }

      if (debug_)
      {
        log_ << "   - final score: " << score << std::endl;
      }
      if (score > max_score)
      {
        max_score = score;
        best_pattern = pattern;
      }
    }
    if (debug_)
    {
      log_ << " - best score              : " << max_score << std::endl;
    }
    best_pattern.theoretical_pattern = isotopes;
    return max_score;
  }

  void FeatureFinderAlgorithmPicked::extendMassTraces_(const IsotopePattern& pattern, MassTraces& traces, Size meta_index_overall) const
  {
    //find index of the trace with the maximum intensity
    double max_int =  0.0;
    Size max_trace_index = 0;
    for (Size p = 0; p < pattern.peak.size(); ++p)
    {
      if (pattern.peak[p] < 0)
      {
        continue; //skip missing and removed traces
      }
      if (map_[pattern.spectrum[p]][pattern.peak[p]].getIntensity() > max_int)
      {
        max_int = map_[pattern.spectrum[p]][pattern.peak[p]].getIntensity();
        max_trace_index = p;
      }
    }

    //extend the maximum intensity trace to determine the boundaries in RT dimension
    Size start_index = pattern.spectrum[max_trace_index];
    const PeakType* start_peak = &(map_[pattern.spectrum[max_trace_index]][pattern.peak[max_trace_index]]);
    double start_mz = start_peak->getMZ();
    double start_rt = map_[start_index].getRT();
    if (debug_)
    {
      log_ << " - Trace " << max_trace_index << " (maximum intensity)" << std::endl;
    }
    if (debug_)
    {
      log_ << "   - extending from: " << map_[start_index].getRT() << " / " << start_mz << " (int: " << start_peak->getIntensity() << ")" << std::endl;
    }
    //initialize the trace and extend
    MassTrace max_trace;
    max_trace.peaks.emplace_back(start_rt, start_peak);
    extendMassTrace_(max_trace, start_index, start_mz, false, meta_index_overall);
    extendMassTrace_(max_trace, start_index, start_mz, true, meta_index_overall);

    double rt_max = max_trace.peaks.back().first;
    double rt_min = max_trace.peaks.begin()->first;
    if (debug_)
    {
      log_ << "   - rt bounds: " << rt_min << "-" << rt_max << std::endl;
    }
    //Abort if too few peak were found
    if (!max_trace.isValid() || max_trace.peaks.size() < 2 * min_spectra_ - max_missing_trace_peaks_)
    {
      if (debug_)
      {
        log_ << "   - could not extend trace with maximum intensity => abort" << std::endl;
      }
      return;
    }
    for (Size p = 0; p < pattern.peak.size(); ++p)
    {
      if (debug_)
      {
        log_ << " - Trace " << p << std::endl;
      }
      if (p == max_trace_index)
      {
        if (debug_)
        {
          log_ << "   - previously extended maximum trace" << std::endl;
        }
        traces.push_back(std::move(max_trace));
        traces.back().theoretical_int = pattern.theoretical_pattern.intensity[p];
        traces.max_trace = traces.size() - 1;
        continue;
      }
      Seed starting_peak;
      starting_peak.spectrum = pattern.spectrum[p];
      starting_peak.peak = pattern.peak[p];
      if (pattern.peak[p] == -2)
      {
        if (debug_)
        {
          log_ << "   - removed during isotope fit" << std::endl;
        }
        continue;
      }
      else if (pattern.peak[p] == -1)
      {
        if (debug_)
        {
          log_ << "   - missing" << std::endl;
        }
        continue;
      }
      starting_peak.intensity = map_[starting_peak.spectrum][starting_peak.peak].getIntensity();
      if (debug_) log_ << "   - trace seed: " << map_[starting_peak.spectrum].getRT() << " / " << map_[starting_peak.spectrum][starting_peak.peak].getMZ() << " (int: " << map_[starting_peak.spectrum][starting_peak.peak].getIntensity() << ")" << std::endl;

      //search for nearby maximum of the mass trace as the extension assumes that it starts at the maximum
      Size begin = std::max((Size)0, starting_peak.spectrum - min_spectra_);
      Size end = std::min(starting_peak.spectrum + min_spectra_, (Size)map_.size());
      double mz = map_[starting_peak.spectrum][starting_peak.peak].getMZ();
      double inte = map_[starting_peak.spectrum][starting_peak.peak].getIntensity();
      for (Size spectrum_index = begin; spectrum_index < end; ++spectrum_index)
      {
        //find better seeds (no-empty scan/low mz diff/higher intensity)
        SignedSize peak_index = -1;
        if (!map_[spectrum_index].empty())
        {
          peak_index = map_[spectrum_index].findNearest(map_[starting_peak.spectrum][starting_peak.peak].getMZ());
        }

        if (peak_index < 0 ||
            map_[spectrum_index][peak_index].getIntensity() <= inte ||
            std::fabs(mz - map_[spectrum_index][peak_index].getMZ()) >= pattern_tolerance_
            )
        {
          continue;
        }

        starting_peak.spectrum = spectrum_index;
        starting_peak.peak = peak_index;
        inte = map_[spectrum_index][peak_index].getIntensity();
      }
      if (debug_)
      {
        log_ << "   - extending from: " << map_[starting_peak.spectrum].getRT() << " / " << map_[starting_peak.spectrum][starting_peak.peak].getMZ() << " (int: " << map_[starting_peak.spectrum][starting_peak.peak].getIntensity() << ")" << std::endl;
      }
      //------------------------------------------------------------------
      //Extend seed to a mass trace
      MassTrace trace;
      const PeakType* seed = &(map_[starting_peak.spectrum][starting_peak.peak]);
      //initialize trace with seed data and extend
      trace.peaks.emplace_back(map_[starting_peak.spectrum].getRT(), seed);
      extendMassTrace_(trace, starting_peak.spectrum, seed->getMZ(), false, meta_index_overall, rt_min, rt_max);
      extendMassTrace_(trace, starting_peak.spectrum, seed->getMZ(), true, meta_index_overall, rt_min, rt_max);

      //check if enough peaks were found
      if (!trace.isValid())
      {
        if (debug_)
        {
          log_ << "   - could not extend trace " << std::endl;
        }
        //Missing traces in the middle of a pattern are not acceptable => fix this
        if (p < traces.max_trace)
        {
          traces.clear(); //remove earlier traces
          continue;
        }
        else if (p > traces.max_trace)
        {
          break; //no more traces are possible
        }
      }
      traces.push_back(trace);
      traces.back().theoretical_int = pattern.theoretical_pattern.intensity[p];
    }
  }

  void FeatureFinderAlgorithmPicked::extendMassTrace_(MassTrace& trace, SignedSize spectrum_index, double mz, bool increase_rt, Size meta_index_overall, double min_rt, double max_rt) const
  {
    //Reverse peaks if we run the method for the second time (to keep them in chronological order)
    if (increase_rt)
    {
      ++spectrum_index;
      std::reverse(trace.peaks.begin(), trace.peaks.end());
    }
    else
    {
      --spectrum_index;
    }

    //check if boundaries are set
    bool boundaries = false;
    if (max_rt != min_rt)
    {
      boundaries = true;
    }

    //Relax slope threshold if there is a hard boundary for the extension
    double current_slope_bound = (1.0 + (double)boundaries) * slope_bound_;

    Size delta_count = min_spectra_;
    std::vector<double> deltas(delta_count - 1, 0);

    double last_observed_intensity = trace.peaks.back().second->getIntensity();

    UInt missing_peaks = 0;
    Size peaks_before_extension = trace.peaks.size();
    String abort_reason = "";

    while ((!increase_rt && spectrum_index >= 0) || (increase_rt && spectrum_index < (SignedSize)map_.size()))
    {
      if (boundaries &&
          ((!increase_rt && map_[spectrum_index].getRT() < min_rt) ||
           (increase_rt && map_[spectrum_index].getRT() > max_rt))
          )
      {
        abort_reason = "Hit upper/lower boundary";
        break;
      }

      SignedSize peak_index = -1;

      if (!map_[spectrum_index].empty())
      {
        peak_index = map_[spectrum_index].findNearest(mz);
      }

      // check if the peak is "missing"
      if (
        peak_index < 0 // no peak found
         || map_[spectrum_index].getFloatDataArrays()[meta_index_overall][peak_index] < 0.01 // overall score is to low
         || positionScore_(mz, map_[spectrum_index][peak_index].getMZ(), trace_tolerance_) == 0.0 // deviation of mz is too big
        )
      {
        ++missing_peaks;

        if (missing_peaks > max_missing_trace_peaks_)
        {
          abort_reason = "too many peaks missing";
          break;
        }
      }
      else
      {
        missing_peaks = 0;

        //add found peak to trace
        trace.peaks.emplace_back(map_[spectrum_index].getRT(), &(map_[spectrum_index][peak_index]));

        //update deltas and intensities
        deltas.push_back((map_[spectrum_index][peak_index].getIntensity() - last_observed_intensity) / last_observed_intensity);
        last_observed_intensity = map_[spectrum_index][peak_index].getIntensity();

        //Abort if the average delta is too big (as intensity increases then)
        double average_delta = std::accumulate(deltas.end() - delta_count, deltas.end(), 0.0) / (double)delta_count;
        if (average_delta > current_slope_bound)
        {
          abort_reason = String("Average delta above threshold: ") + average_delta + "/" + current_slope_bound;

          //remove last peaks as we extended too far
          Size remove = std::min((Size)(trace.peaks.size() - peaks_before_extension), delta_count - 1);
          trace.peaks.erase(trace.peaks.end() - remove, trace.peaks.end());
          break;
        }
      }

      //increase/decrease scan index
      if (increase_rt)
      {
        ++spectrum_index;
      }
      else
      {
        --spectrum_index;
      }
    }
    if (debug_)
    {
      log_ << "   - Added " << (trace.peaks.size() - peaks_before_extension) << " peaks (abort: " << abort_reason << ")" << std::endl;
    }
  }

  Size FeatureFinderAlgorithmPicked::nearest_(double pos, const MSSpectrum& spec, Size start) const
  {
    Size index = start;
    double distance = std::fabs(pos - spec[index].getMZ());
    ++index;
    while (index < spec.size())
    {
      double new_distance = std::fabs(pos - spec[index].getMZ());
      if (new_distance < distance)
      {
        distance = new_distance;
        ++index;
      }
      else
      {
        break;
      }
    }
    return --index;
  }

  void FeatureFinderAlgorithmPicked::findIsotope_(double pos, Size spectrum_index, IsotopePattern& pattern, Size pattern_index, Size& peak_index) const
  {
    if (debug_)
    {
      log_ << "   - Isotope " << pattern_index << ": ";
    }
    double intensity = 0.0;
    double pos_score = 0.0;
    UInt matches = 0;

    //search in the center spectrum
    const SpectrumType& spectrum = map_[spectrum_index];
    peak_index = nearest_(pos, spectrum, peak_index);
    double this_mz_score = positionScore_(pos, spectrum[peak_index].getMZ(), pattern_tolerance_);
    pattern.theoretical_mz[pattern_index] = pos;

    if (this_mz_score != 0.0)
    {
      if (debug_)
      {
        log_ << String::number(spectrum[peak_index].getIntensity(), 1) << " ";
      }
      pattern.peak[pattern_index] = peak_index;
      pattern.spectrum[pattern_index] = spectrum_index;
      intensity += spectrum[peak_index].getIntensity();
      pos_score += this_mz_score;
      ++matches;
    }

    //previous spectrum
    if (spectrum_index != 0 && !map_[spectrum_index - 1].empty())
    {
      const SpectrumType& spectrum_before = map_[spectrum_index - 1];
      Size index_before = spectrum_before.findNearest(pos);
      double mz_score = positionScore_(pos, spectrum_before[index_before].getMZ(), pattern_tolerance_);
      if (mz_score != 0.0)
      {
        if (debug_) log_ << String::number(spectrum_before[index_before].getIntensity(), 1) << "b ";
        intensity += spectrum_before[index_before].getIntensity();
        pos_score += mz_score;
        ++matches;

        if (pattern.peak[pattern_index] == -1)
        {
          pattern.peak[pattern_index] = index_before;
          pattern.spectrum[pattern_index] = spectrum_index - 1;
        }
      }
    }

    //next spectrum
    if (spectrum_index != map_.size() - 1 && !map_[spectrum_index + 1].empty())
    {
      const SpectrumType& spectrum_after = map_[spectrum_index + 1];
      Size index_after = spectrum_after.findNearest(pos);
      double mz_score = positionScore_(pos, spectrum_after[index_after].getMZ(), pattern_tolerance_);
      if (mz_score != 0.0)
      {
        if (debug_) log_ << String::number(spectrum_after[index_after].getIntensity(), 1) << "a ";
        intensity += spectrum_after[index_after].getIntensity();
        pos_score += mz_score;
        ++matches;

        if (pattern.peak[pattern_index] == -1)
        {
          pattern.peak[pattern_index] = index_after;
          pattern.spectrum[pattern_index] = spectrum_index + 1;
        }
      }
    }

    //no isotope found
    if (matches == 0)
    {
      if (debug_)
      {
        log_ << " missing" << std::endl;
      }
      pattern.peak[pattern_index] = -1;
      pattern.mz_score[pattern_index] = 0.0;
      pattern.intensity[pattern_index] = 0.0;
    }
    else
    {
      if (debug_)
      {
        log_ << "=> " << intensity / matches << std::endl;
      }
      pattern.mz_score[pattern_index] = pos_score / matches;
      pattern.intensity[pattern_index] = intensity / matches;
    }
  }

  double FeatureFinderAlgorithmPicked::positionScore_(double pos1, double pos2, double allowed_deviation) const
  {
    double diff = fabs(pos1 - pos2);
    if (diff <= 0.5 * allowed_deviation)
    {
      return 0.1 * (0.5 * allowed_deviation - diff) / (0.5 * allowed_deviation) + 0.9;
    }
    else if (diff <= allowed_deviation)
    {
      return 0.9 * (allowed_deviation - diff) / (0.5 * allowed_deviation);
    }
    return 0.0;
  }

  /// Calculates a score between 0 and 1 for the correlation between theoretical and found isotope pattern
  double FeatureFinderAlgorithmPicked::isotopeScore_(const TheoreticalIsotopePattern& isotopes, IsotopePattern& pattern, bool consider_mz_distances) const
  {
    if (debug_) log_ << "   - fitting " << pattern.intensity.size() << " peaks" << std::endl;
    //Abort if a core peak is missing
    for (Size iso = 0 + isotopes.optional_begin; iso < pattern.peak.size() - isotopes.optional_end; ++iso)
    {
      if (pattern.peak[iso] == -1)
      {
        if (debug_)
        {
          log_ << "   - aborting: core peak is missing" << std::endl;
        }
        return 0.0;
      }
    }
    //Find best isotope fit
    // - try to leave out optional isotope peaks to improve the fit
    // - do not allow gaps inside the pattern
    double best_int_score = 0.01; //Not 0 as this would result in problems when checking for the percental improvement
    Size best_begin = 0;
    for (Size i = isotopes.optional_begin; i > 0; --i)
    {
      if (pattern.peak[i - 1] == -1)
      {
        best_begin = i;
        break;
      }
    }
    Size best_end = 0;
    for (Size i = isotopes.optional_end; i > 0; --i)
    {
      if (pattern.peak[pattern.peak.size() - i] == -1)
      {
        best_end = i;
        break;
      }
    }
    if (debug_)
    {
      log_ << "   - best_begin/end: " << best_begin << "/" << best_end << std::endl;
    }
    for (Size b = best_begin; b <= isotopes.optional_begin; ++b)
    {
      for (Size e = best_end; e <= isotopes.optional_end; ++e)
      {
        //Make sure we have more than 2 peaks (unless in the first loop iteration, there we allow two points)
        if (isotopes.size() - b - e > 2 || (b == best_begin &&
                                            e == best_end &&
                                            isotopes.size() - b - e > 1))
        {
          double int_score = Math::pearsonCorrelationCoefficient(isotopes.intensity.begin() + b, isotopes.intensity.end() - e, pattern.intensity.begin() + b, pattern.intensity.end() - e);
          if (std::isnan(int_score))
          {
            int_score = 0.0;
          }
          if (isotopes.size() - b - e == 2 && int_score > min_isotope_fit_)
          {
            int_score = min_isotope_fit_; //special case for the first loop iteration (otherwise the score is 1)
          }
          if (debug_)
          {
            log_ << "   - fit (" << b << "/" << e << "): " << int_score;
          }
          if (int_score / best_int_score >= 1.0 + optional_fit_improvement_)
          {
            if (debug_)
            {
              log_ << " - new best fit ";
            }
            best_int_score = int_score;
            best_begin = b;
            best_end = e;
          }
          if (debug_)
          {
            log_ << std::endl;
          }
        }
      }
    }

    //if the best fit is empty, abort
    if (pattern.mz_score.size() - best_begin - best_end == 0)
    {
      return 0.0;
    }

    //remove left out peaks from the beginning
    for (Size i = 0; i < best_begin; ++i)
    {
      pattern.peak[i] = -2;
      pattern.intensity[i] = 0.0;
      pattern.mz_score[i] = 0.0;
    }
    //remove left out peaks from the end
    for (Size i = 0; i < best_end; ++i)
    {
      pattern.peak[isotopes.size() - 1 - i] = -2;
      pattern.intensity[isotopes.size() - 1 - i] = 0.0;
      pattern.mz_score[isotopes.size() - 1 - i] = 0.0;
    }
    //calculate m/z score (if required)
    if (consider_mz_distances)
    {
      best_int_score *= std::accumulate(pattern.mz_score.begin() + best_begin, pattern.mz_score.end() - best_end, 0.0) / (pattern.mz_score.size() - best_begin - best_end);
    }

    //return final score
    OPENMS_POSTCONDITION(best_int_score >= 0.0, (String("Internal error: Isotope score (") + best_int_score + ") should be >=0.0").c_str())
    OPENMS_POSTCONDITION(best_int_score <= 1.0, (String("Internal error: Isotope score (") + best_int_score + ") should be <=1.0").c_str())
    return best_int_score;
  }

  double FeatureFinderAlgorithmPicked::intensityScore_(Size spectrum, Size peak) const
  {
    // calculate (half) bin numbers
    double intensity  = map_[spectrum][peak].getIntensity();
    double rt = map_[spectrum].getRT();
    double mz = map_[spectrum][peak].getMZ();
    double rt_min = map_.getMinRT();
    double mz_min = map_.getMinMZ();
    UInt rt_bin = std::min(2 * intensity_bins_ - 1, (UInt) std::floor((rt - rt_min) / intensity_rt_step_ * 2.0));
    UInt mz_bin = std::min(2 * intensity_bins_ - 1, (UInt) std::floor((mz - mz_min) / intensity_mz_step_ * 2.0));
    // determine mz bins
    UInt ml, mh;
    if (mz_bin == 0 || mz_bin == 2 * intensity_bins_ - 1)
    {
      ml = mz_bin / 2;
      mh = mz_bin / 2;
    }
    else if (Math::isOdd(mz_bin))
    {
      ml = mz_bin / 2;
      mh = mz_bin / 2 + 1;
    }
    else
    {
      ml = mz_bin / 2 - 1;
      mh = mz_bin / 2;
    }
    // determine rt bins
    UInt rl, rh;
    if (rt_bin == 0 || rt_bin == 2 * intensity_bins_ - 1)
    {
      rl = rt_bin / 2;
      rh = rt_bin / 2;
    }
    else if (Math::isOdd(rt_bin))
    {
      rl = rt_bin / 2;
      rh = rt_bin / 2 + 1;
    }
    else
    {
      rl = rt_bin / 2 - 1;
      rh = rt_bin / 2;
    }
    // calculate distances to surrounding bin centers (normalized to [0,1])
    double drl = std::fabs(rt_min + (0.5 + rl) * intensity_rt_step_ - rt) / intensity_rt_step_;
    double drh = std::fabs(rt_min + (0.5 + rh) * intensity_rt_step_ - rt) / intensity_rt_step_;
    double dml = std::fabs(mz_min + (0.5 + ml) * intensity_mz_step_ - mz) / intensity_mz_step_;
    double dmh = std::fabs(mz_min + (0.5 + mh) * intensity_mz_step_ - mz) / intensity_mz_step_;
    // Calculate weights for the intensity scores based on the distances to the
    // bin center(the nearer to better)
    double d1 = std::sqrt(std::pow(1.0 - drl, 2) + std::pow(1.0 - dml, 2));
    double d2 = std::sqrt(std::pow(1.0 - drh, 2) + std::pow(1.0 - dml, 2));
    double d3 = std::sqrt(std::pow(1.0 - drl, 2) + std::pow(1.0 - dmh, 2));
    double d4 = std::sqrt(std::pow(1.0 - drh, 2) + std::pow(1.0 - dmh, 2));
    double d_sum = d1 + d2 + d3 + d4;
    // Final score .. intensityScore in the surrounding bins, weighted by the distance of the
    // bin center to the peak
    double final = intensityScore_(rl, ml, intensity) * (d1 / d_sum)
                   + intensityScore_(rh, ml, intensity) * (d2 / d_sum)
                   + intensityScore_(rl, mh, intensity) * (d3 / d_sum)
                   + intensityScore_(rh, mh, intensity) * (d4 / d_sum);

    OPENMS_POSTCONDITION(final >= 0.0, (String("Internal error: Intensity score (") + final + ") should be >=0.0").c_str())
    OPENMS_POSTCONDITION(final <= 1.0001, (String("Internal error: Intensity score (") + final + ") should be <=1.0").c_str())
    return final;
  }

  std::unique_ptr<TraceFitter> FeatureFinderAlgorithmPicked::chooseTraceFitter_(double& tau)
  {
    // choose fitter
    if (param_.getValue("feature:rt_shape") == "asymmetric")
    {
      OPENMS_LOG_DEBUG << "use asymmetric rt peak shape" << std::endl;
      tau = -1.0;
      return std::make_unique<EGHTraceFitter>();
    }
    else // if (param_.getValue("feature:rt_shape") == "symmetric")
    {
      OPENMS_LOG_DEBUG << "use symmetric rt peak shape" << std::endl;
      return std::make_unique<GaussTraceFitter>();
    }
  }

  double FeatureFinderAlgorithmPicked::intensityScore_(Size rt_bin, Size mz_bin, double intensity) const
  {
    // interpolate score value according to quantiles(20)
    const std::vector<double>& quantiles20 = intensity_thresholds_[rt_bin][mz_bin];
    // get iterator pointing to quantile that is >= intensity
    std::vector<double>::const_iterator it = std::lower_bound(quantiles20.begin(), quantiles20.end(), intensity);
    // bigger than the biggest value => return 1.0
    if (it == quantiles20.end())
    {
      return 1.0;
    }
    // interpolate inside the bin
    double bin_score = 0.0;
    if (it == quantiles20.begin())
    {
      bin_score = 0.05 * intensity / *it;
    }
    else
    {
      // (intensity - vigintile_low) / (vigintile_high - vigintile_low)
      bin_score = 0.05 * (intensity - *(it - 1)) / (*it - *(it - 1));
    }

    double final = bin_score +
                   0.05 * ((it - quantiles20.begin()) - 1.0); // determine position of lower bound in the vector

    //fix numerical problems
    if (final < 0.0)
    {
      final = 0.0;
    }
    if (final > 1.0)
    {
      final = 1.0;
    }
    // final = 1/20 * [ index(vigintile_low) + (intensity-vigintile_low) / (vigintile_high - vigintile_low) ]
    return final;
  }

  void FeatureFinderAlgorithmPicked::cropFeature_(const std::shared_ptr<TraceFitter>& fitter,
                                                  const MassTraces& traces,
                                                  MassTraces& new_traces)
  {
    double low_bound = fitter->getLowerRTBound();
    double high_bound = fitter->getUpperRTBound();

    if (debug_)
    {
      log_ << "    => RT bounds: " << low_bound << " - " << high_bound << std::endl;
    }
    std::vector<double> v_theo, v_real;
    for (Size t = 0; t < traces.size(); ++t)
    {
      const MassTrace& trace = traces[t];
      if (debug_)
      {
        log_ << "   - Trace " << t << ": (" << trace.theoretical_int << ")" << std::endl;
      }
      MassTrace new_trace;
      //compute average relative deviation and correlation
      double deviation = 0.0;
      v_theo.clear();
      v_real.clear();
      for (Size k = 0; k < trace.peaks.size(); ++k)
      {
        //consider peaks when inside RT bounds only
        if (trace.peaks[k].first >= low_bound && trace.peaks[k].first <= high_bound)
        {
          new_trace.peaks.push_back(trace.peaks[k]);

          double theo = traces.baseline + fitter->computeTheoretical(trace, k);
          v_theo.push_back(theo);
          double real = trace.peaks[k].second->getIntensity();
          v_real.push_back(real);
          deviation += std::fabs(real - theo) / theo;
        }
      }
      double fit_score = 0.0;
      double correlation = 0.0;
      double final_score = 0.0;
      if (!new_trace.peaks.empty())
      {
        fit_score = deviation / new_trace.peaks.size();
        correlation = std::max(0.0, Math::pearsonCorrelationCoefficient(v_theo.begin(), v_theo.end(), v_real.begin(), v_real.end()));
        final_score = std::sqrt(correlation * std::max(0.0, 1.0 - fit_score));
      }
      if (debug_)
      {
        log_ << "     - peaks: " << new_trace.peaks.size() << " / " << trace.peaks.size() << " - relative deviation: " << fit_score << " - correlation: " << correlation << " - final score: " << correlation << std::endl;
      }
      //remove badly fitting traces
      if (!new_trace.isValid() || final_score < min_trace_score_)
      {
        if (t < traces.max_trace)
        {
          new_traces = MassTraces();
          if (debug_)
          {
            log_ << "     - removed this and previous traces due to bad fit" << std::endl;
          }
          new_traces.clear(); //remove earlier traces
          continue;
        }
        else if (t == traces.max_trace)
        {
          new_traces = MassTraces();
          if (debug_)
          {
            log_ << "     - aborting (max trace was removed)" << std::endl;
          }
          break;
        }
        else if (t > traces.max_trace)
        {
          if (debug_)
          {
            log_ << "     - removed due to bad fit => omitting the rest" << std::endl;
          }
          break; //no more traces are possible
        }
      }
      //add new trace
      else
      {
        new_trace.theoretical_int = trace.theoretical_int;
        new_traces.push_back(new_trace);
        if (t == traces.max_trace)
        {
          new_traces.max_trace = new_traces.size() - 1;
        }
      }
    }
    new_traces.baseline = traces.baseline;
  }

  bool FeatureFinderAlgorithmPicked::checkFeatureQuality_(const std::shared_ptr<TraceFitter>& fitter,
                                                          MassTraces& feature_traces,
                                                          const double& seed_mz, const double& min_feature_score,
                                                          String& error_msg, double& fit_score, double& correlation, double& final_score)
  {
    //check if the sigma fit was ok (if it is larger than 'max_rt_span')
    // 5.0 * sigma > max_rt_span_ * region_rt_span
    if (fitter->checkMaximalRTSpan(max_rt_span_))
    {
      error_msg = "Invalid fit: Fitted model is bigger than 'max_rt_span'";
      return false;
    }

    //check if the feature is valid
    if (!feature_traces.isValid(seed_mz, trace_tolerance_))
    {
      error_msg = "Invalid feature after fit - too few traces or peaks left";
      return false;
    }

    //check if x0 is inside feature bounds
    {
      std::pair<double, double> rt_bounds = feature_traces.getRTBounds();
      if (fitter->getCenter() < rt_bounds.first || fitter->getCenter() > rt_bounds.second)
      {
        error_msg = "Invalid fit: Center outside of feature bounds";
        return false;
      }
    }

    //check if the remaining traces fill out at least 'min_rt_span' of the RT span
    {
      std::pair<double, double> rt_bounds = feature_traces.getRTBounds();
      if (fitter->checkMinimalRTSpan(rt_bounds, min_rt_span_))
      {
        error_msg = "Invalid fit: Less than 'min_rt_span' left after fit";
        return false;
      }
    }

    //check if feature quality is high enough (average relative deviation and correlation of the whole feature)
    {
      std::vector<double> v_theo, v_real;
      double deviation = 0.0;
      for (Size t = 0; t < feature_traces.size(); ++t)
      {
        MassTrace& trace = feature_traces[t];
        for (Size k = 0; k < trace.peaks.size(); ++k)
        {
          // was double theo = new_traces.baseline + trace.theoretical_int *  height * exp(-0.5 * pow(trace.peaks[k].first - x0, 2) / pow(sigma, 2) );
          double theo = feature_traces.baseline + fitter->computeTheoretical(trace, k);
          v_theo.push_back(theo);
          double real = trace.peaks[k].second->getIntensity();
          v_real.push_back(real);
          deviation += std::fabs(real - theo) / theo;
        }
      }
      fit_score = std::max(0.0, 1.0 - (deviation / feature_traces.getPeakCount()));
      correlation = std::max(0.0, Math::pearsonCorrelationCoefficient(v_theo.begin(), v_theo.end(), v_real.begin(), v_real.end()));
      final_score = std::sqrt(correlation * fit_score);

      //quality output
      if (debug_)
      {
        log_ << "Quality estimation:" << std::endl;
        log_ << " - relative deviation: " << fit_score << std::endl;
        log_ << " - correlation: " << correlation << std::endl;
        log_ << " => final score: " << final_score << std::endl;
      }

      if (final_score < min_feature_score)
      {
        error_msg = "Feature quality too low after fit";
        return false;
      }


    }

    return true;
  }

  void FeatureFinderAlgorithmPicked::writeFeatureDebugInfo_(const std::shared_ptr<TraceFitter>& fitter,
                                                            const MassTraces& traces,
                                                            const MassTraces& new_traces,
                                                            bool feature_ok, const String& error_msg, const double final_score, const Int plot_nr, const PeakType& peak,
                                                            const String& path)
  {

    double pseudo_rt_shift = param_.getValue("debug:pseudo_rt_shift");
    String script;
    {
      TextFile tf;
      //gnuplot script
      script = String("plot \"") + path + plot_nr + ".dta\" title 'before fit (RT: " +  String::number(fitter->getCenter(), 2) + " m/z: " +  String::number(peak.getMZ(), 4) + ")' with points 1";
      //feature before fit
      for (Size k = 0; k < traces.size(); ++k)
      {
        for (Size j = 0; j < traces[k].peaks.size(); ++j)
        {
          tf.addLine(String(pseudo_rt_shift * k + traces[k].peaks[j].first) + "\t" + traces[k].peaks[j].second->getIntensity());
        }
      }
      tf.store(path + plot_nr + ".dta");
    }

    {
      //fitted feature
      if (new_traces.getPeakCount() != 0)
      {
        TextFile tf_new_trace;
        for (Size k = 0; k < new_traces.size(); ++k)
        {
          for (Size j = 0; j < new_traces[k].peaks.size(); ++j)
          {
            tf_new_trace.addLine(String(pseudo_rt_shift * k + new_traces[k].peaks[j].first) + "\t" + new_traces[k].peaks[j].second->getIntensity());
          }
        }

        tf_new_trace.store(path + plot_nr + "_cropped.dta");
        script = script + ", \"" + path + plot_nr + "_cropped.dta\" title 'feature ";

        if (!feature_ok)
        {
          script = script + " - " + error_msg;
        }
        else
        {
          script = script + (features_->size() + 1) + " (score: " +  String::number(final_score, 3) + ")";
        }
        script = script + "' with points 3";
      }
    }

    {
      //fitted functions
      TextFile tf_fitted_func;
      for (Size k = 0; k < traces.size(); ++k)
      {
        char fun = 'f';
        fun += (char)k;
        tf_fitted_func.addLine(fitter->getGnuplotFormula(traces[k], fun, traces.baseline, pseudo_rt_shift * k));
        //tf.push_back(String(fun)+"(x)= " + traces.baseline + " + " + fitter->getGnuplotFormula(traces[k], pseudo_rt_shift * k));
        script =  script + ", " + fun + "(x) title 'Trace " + k + " (m/z: " + String::number(traces[k].getAvgMZ(), 4) + ")'";
      }

      //output
      tf_fitted_func.addLine("set xlabel \"pseudo RT (mass traces side-by-side)\"");
      tf_fitted_func.addLine("set ylabel \"intensity\"");
      tf_fitted_func.addLine("set samples 1000");
      tf_fitted_func.addLine(script);
      tf_fitted_func.addLine("pause -1");
      tf_fitted_func.store(path + plot_nr + ".plot");
    }
  }

}
