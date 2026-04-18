// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h> // integrateWindow
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAPrescoring.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/MATH/MathFunctions.h> // getPPM
#include <OpenMS/MATH/StatisticFunctions.h>

#include <numeric>
#include <algorithm>
#include <functional>
#include <iostream>
#include <unordered_map>

#include <cmath> // for isnan
#include <utility>

const double C13C12_MASSDIFF_U = 1.0033548;

namespace OpenMS
{
  namespace
  {
    void hashCombine(std::size_t& seed, std::size_t value)
    {
      seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    struct IntegratedDIAWindow
    {
      double mz = -1.0;
      double im = -1.0;
      double intensity = 0.0;
      bool found = false;
    };

    struct DIAWindowKey
    {
      double mz_min = 0.0;
      double mz_max = 0.0;

      bool operator==(const DIAWindowKey& rhs) const
      {
        return mz_min == rhs.mz_min &&
               mz_max == rhs.mz_max;
      }
    };

    struct DIAWindowKeyHash
    {
      std::size_t operator()(const DIAWindowKey& key) const
      {
        std::size_t seed = std::hash<double>{}(key.mz_min);
        hashCombine(seed, std::hash<double>{}(key.mz_max));
        return seed;
      }
    };

    struct DIAIntegrationCache
    {
      std::vector<const OpenSwath::Spectrum*> spectra_key;
      std::unordered_map<DIAWindowKey, IntegratedDIAWindow, DIAWindowKeyHash> windows;
      bool valid = false;
      bool im_empty = true;
      double im_min = 0.0;
      double im_max = 0.0;
      bool centroided = false;

      bool matches(const SpectrumSequence& spectrum, const RangeMobility& im_range, bool centroided_in) const
      {
        if (!valid || centroided != centroided_in || im_empty != im_range.isEmpty() || spectra_key.size() != spectrum.size())
        {
          return false;
        }
        if (!im_empty && (im_min != im_range.getMin() || im_max != im_range.getMax()))
        {
          return false;
        }
        for (Size i = 0; i < spectrum.size(); ++i)
        {
          if (spectra_key[i] != spectrum[i].get())
          {
            return false;
          }
        }
        return true;
      }

      void begin(const SpectrumSequence& spectrum, const RangeMobility& im_range, bool centroided_in, bool clear_windows)
      {
        if (matches(spectrum, im_range, centroided_in))
        {
          if (clear_windows)
          {
            windows.clear();
          }
          return;
        }

        spectra_key.clear();
        spectra_key.reserve(spectrum.size());
        for (const OpenSwath::SpectrumPtr& spec : spectrum)
        {
          spectra_key.push_back(spec.get());
        }
        im_empty = im_range.isEmpty();
        if (!im_empty)
        {
          im_min = im_range.getMin();
          im_max = im_range.getMax();
        }
        else
        {
          im_min = 0.0;
          im_max = 0.0;
        }
        centroided = centroided_in;
        windows.clear();
        valid = true;
      }

      const IntegratedDIAWindow& integrate(const SpectrumSequence& spectrum,
                                           const RangeMZ& mz_range,
                                           const RangeMobility& im_range,
                                           bool centroided_in)
      {
        begin(spectrum, im_range, centroided_in, false);
        const DIAWindowKey key{mz_range.getMin(), mz_range.getMax()};
        auto cached = windows.find(key);
        if (cached != windows.end())
        {
          return cached->second;
        }

        IntegratedDIAWindow result;
        result.found = DIAHelpers::integrateWindow(spectrum, result.mz, result.im, result.intensity, mz_range, im_range, centroided_in);
        return windows.emplace(key, result).first->second;
      }
    };

    thread_local DIAIntegrationCache dia_integration_cache;

    struct DIAScoringPool
    {
      std::vector<double> isotopes_int;
      std::vector<double> exp_isotopes_int;
      std::vector<double> prescore_int_exp;
      std::vector<double> bseries;
      std::vector<double> yseries;
      std::string byseries_sequence;
      int byseries_charge = 0;
      bool byseries_valid = false;

      void reset()
      {
        isotopes_int.clear();
        exp_isotopes_int.clear();
        prescore_int_exp.clear();
      }
    };
    thread_local DIAScoringPool dias_pool;

    void scorePreparedCached(const SpectrumSequence& spec,
                             const DiaPrescore::PreparedSpectrum& prepared,
                             double dia_extract_window,
                             const RangeMobility& im_range,
                             double& dotprod,
                             double& manhattan)
    {
      std::vector<double>& int_exp = dias_pool.prescore_int_exp;
      int_exp.clear();
      if (prepared.mz_theor.empty())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No windows supplied!");
      }
      if (!spec.empty())
      {
        int_exp.reserve(prepared.mz_theor.size());
        for (double mz_theor : prepared.mz_theor)
        {
          RangeMZ mz_range(mz_theor);
          mz_range.minSpanIfSingular(dia_extract_window);
          const IntegratedDIAWindow& window = dia_integration_cache.integrate(spec, mz_range, im_range, false);
          int_exp.push_back(window.found ? window.intensity : 0.0);
        }
      }

      std::transform(int_exp.cbegin(), int_exp.cend(), int_exp.begin(), [](double val) { return std::sqrt(val); });

      double int_exp_total = std::accumulate(int_exp.cbegin(), int_exp.cend(), 0.0);
      OpenSwath::normalize(int_exp, int_exp_total, int_exp);

      manhattan = OpenSwath::manhattanDist(int_exp.cbegin(), int_exp.cend(), prepared.int_theor.cbegin());

      double int_exp_euclid_norm = OpenSwath::norm(int_exp.cbegin(), int_exp.cend());
      OpenSwath::normalize(int_exp, int_exp_euclid_norm, int_exp);

      dotprod = OpenSwath::dotProd(int_exp.cbegin(), int_exp.cend(), prepared.int_theor_neg.cbegin());
      dotprod = (((dotprod - prepared.neg_val) * 2.) / (prepared.pos_val - prepared.neg_val)) - 1.;
    }

    struct DiaPrescoreThreadCache
    {
      std::vector<double> product_mz;
      std::vector<double> library_intensity;
      std::vector<int> fragment_charge;
      double dia_extract_window = 0.0;
      int nr_isotopes = 0;
      int nr_charges = 0;
      bool valid = false;
      DiaPrescore::PreparedSpectrum prepared;
    };
    thread_local DiaPrescoreThreadCache dia_prescore_cache;

    bool diaPrescoreCacheMatches(const DiaPrescoreThreadCache& cache,
                                 const std::vector<OpenSwath::LightTransition>& transitions,
                                 double dia_extract_window,
                                 int nr_isotopes,
                                 int nr_charges)
    {
      if (!cache.valid ||
          cache.dia_extract_window != dia_extract_window ||
          cache.nr_isotopes != nr_isotopes ||
          cache.nr_charges != nr_charges ||
          cache.product_mz.size() != transitions.size())
      {
        return false;
      }

      for (Size i = 0; i < transitions.size(); ++i)
      {
        if (cache.product_mz[i] != transitions[i].getProductMZ() ||
            cache.library_intensity[i] != transitions[i].getLibraryIntensity() ||
            cache.fragment_charge[i] != transitions[i].fragment_charge)
        {
          return false;
        }
      }
      return true;
    }

    void updateDiaPrescoreCacheKey(DiaPrescoreThreadCache& cache,
                                   const std::vector<OpenSwath::LightTransition>& transitions,
                                   double dia_extract_window,
                                   int nr_isotopes,
                                   int nr_charges)
    {
      cache.valid = false;
      cache.product_mz.clear();
      cache.library_intensity.clear();
      cache.fragment_charge.clear();
      cache.product_mz.reserve(transitions.size());
      cache.library_intensity.reserve(transitions.size());
      cache.fragment_charge.reserve(transitions.size());
      for (const OpenSwath::LightTransition& transition : transitions)
      {
        cache.product_mz.push_back(transition.getProductMZ());
        cache.library_intensity.push_back(transition.getLibraryIntensity());
        cache.fragment_charge.push_back(transition.fragment_charge);
      }
      cache.dia_extract_window = dia_extract_window;
      cache.nr_isotopes = nr_isotopes;
      cache.nr_charges = nr_charges;
      cache.valid = true;
    }

    struct AveragineIsotopeDistributionKey
    {
      double neutral_mass;
      int nr_isotopes;

      bool operator==(const AveragineIsotopeDistributionKey& rhs) const
      {
        return neutral_mass == rhs.neutral_mass &&
               nr_isotopes == rhs.nr_isotopes;
      }
    };

    struct AveragineIsotopeDistributionKeyHash
    {
      std::size_t operator()(const AveragineIsotopeDistributionKey& key) const
      {
        std::size_t seed = std::hash<double>{}(key.neutral_mass);
        hashCombine(seed, std::hash<int>{}(key.nr_isotopes));
        return seed;
      }
    };

    const IsotopeDistribution& getCachedAveragineIsotopeDistribution(double neutral_mass, int nr_isotopes)
    {
      static thread_local std::unordered_map<AveragineIsotopeDistributionKey, IsotopeDistribution, AveragineIsotopeDistributionKeyHash> isotope_cache;

      const AveragineIsotopeDistributionKey key{neutral_mass, nr_isotopes};
      const auto cached = isotope_cache.find(key);
      if (cached != isotope_cache.end())
      {
        return cached->second;
      }

      CoarseIsotopePatternGenerator solver(nr_isotopes);
      return isotope_cache.emplace(key, solver.estimateFromPeptideWeight(neutral_mass)).first->second;
    }
  }

  DIAScoring::DIAScoring() :
    DefaultParamHandler("DIAScoring")
  {

    defaults_.setValue("dia_extraction_window", 0.05, "DIA extraction window in Th or ppm.");
    defaults_.setMinFloat("dia_extraction_window", 0.0);
    defaults_.setValue("dia_extraction_unit", "Th", "DIA extraction window unit");
    defaults_.setValidStrings("dia_extraction_unit", {"Th","ppm"});
    defaults_.setValue("dia_centroided", "false", "Use centroided DIA data.");
    defaults_.setValidStrings("dia_centroided", {"true","false"});
    defaults_.setValue("dia_byseries_intensity_min", 300.0, "DIA b/y series minimum intensity to consider.");
    defaults_.setMinFloat("dia_byseries_intensity_min", 0.0);
    defaults_.setValue("dia_byseries_ppm_diff", 10.0, "DIA b/y series minimal difference in ppm to consider.");
    defaults_.setMinFloat("dia_byseries_ppm_diff", 0.0);

    defaults_.setValue("dia_nr_isotopes", 4, "DIA number of isotopes to consider.");
    defaults_.setMinInt("dia_nr_isotopes", 0);
    defaults_.setValue("dia_nr_charges", 4, "DIA number of charges to consider.");
    defaults_.setMinInt("dia_nr_charges", 0);

    defaults_.setValue("peak_before_mono_max_ppm_diff", 20.0, "DIA maximal difference in ppm to count a peak at lower m/z when searching for evidence that a peak might not be monoisotopic.");
    defaults_.setMinFloat("peak_before_mono_max_ppm_diff", 0.0);

    // write defaults into Param object param_
    defaultsToParam_();

    // for void getBYSeries
    {
      generator = new TheoreticalSpectrumGenerator();
      Param p;
      p.setValue("add_metainfo", "true",
          "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
      generator->setParameters(p);
  }

    // for simulateSpectrumFromAASequence
    //  Param p;
    //  p.setValue("add_metainfo", "false",
    //      "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    //  p.setValue("add_precursor_peaks", "true", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    //  generator->setParameters(p);
  }

  DIAScoring::~DIAScoring()
  {
    delete generator;
  }

  void DIAScoring::updateMembers_()
  {
    dia_extract_window_ = (double)param_.getValue("dia_extraction_window");
    dia_extraction_ppm_ = param_.getValue("dia_extraction_unit") == "ppm";
    dia_centroided_ = param_.getValue("dia_centroided").toBool();
    dia_byseries_intensity_min_ = (double)param_.getValue("dia_byseries_intensity_min");
    dia_byseries_ppm_diff_ = (double)param_.getValue("dia_byseries_ppm_diff");

    dia_nr_isotopes_ = (int)param_.getValue("dia_nr_isotopes");
    dia_nr_charges_ = (int)param_.getValue("dia_nr_charges");
    peak_before_mono_max_ppm_diff_ = (double)param_.getValue("peak_before_mono_max_ppm_diff");
  }

  ///////////////////////////////////////////////////////////////////////////
  // DIA / SWATH scoring

  void DIAScoring::dia_isotope_scores(const std::vector<TransitionType>& transitions, std::vector<SpectrumPtrType>& spectrum,
                                      OpenSwath::IMRMFeature* mrmfeature, const RangeMobility& im_range, double& isotope_corr, double& isotope_overlap) const
  {
    isotope_corr = 0;
    isotope_overlap = 0;
    diaIsotopeScoresSub_(transitions, spectrum, mrmfeature, im_range, isotope_corr, isotope_overlap);
  }

  void DIAScoring::dia_massdiff_score(const std::vector<TransitionType>& transitions,
                                      const SpectrumSequence& spectrum,
                                      const std::vector<double>& normalized_library_intensity,
                                      const RangeMobility& im_range,
                                      double& ppm_score,
                                      double& ppm_score_weighted,
                                      std::vector<double>& diff_ppm) const
  {
    // Calculate the difference of the theoretical mass and the actually measured mass
    ppm_score = 0;
    ppm_score_weighted = 0;
    diff_ppm.clear();
    dia_integration_cache.begin(spectrum, im_range, dia_centroided_, true);
    size_t n_observed = 0;
    for (std::size_t k = 0; k < transitions.size(); k++)
    {
      const TransitionType& transition = transitions[k];

      RangeMZ mz_range = DIAHelpers::createMZRangePPM(transition.getProductMZ(), dia_extract_window_, dia_extraction_ppm_);
      const IntegratedDIAWindow& window = dia_integration_cache.integrate(spectrum, mz_range, im_range, dia_centroided_);
      // Continue if no signal was found - we therefore don't make a statement
      // about the mass difference if no signal is present.
      if (!window.found)
      {
        diff_ppm.push_back(-1); // if no signal is found than we set the ppm to -1
        continue;
      }

      ++n_observed;
      double ppm = Math::getPPM(window.mz, transition.getProductMZ());
      diff_ppm.push_back(ppm);
      ppm_score += std::fabs(ppm);
      ppm_score_weighted += std::fabs(ppm) * normalized_library_intensity[k];
#ifdef MRMSCORING_TESTING
      std::cout << " weighted int of the peak is " << window.mz << " diff is in ppm " << ppm
                << " thus append " << ppm * ppm << " or weighted " << ppm * normalized_library_intensity[k] << '\n';
#endif
    }

    // Only average over transitions where signal was actually found
    ppm_score = (n_observed > 0) ? ppm_score / n_observed : 0.0;
  }

  bool DIAScoring::dia_ms1_massdiff_score(double precursor_mz, const SpectrumSequence& spectrum,
                                          const RangeMobility& im_range, double& ppm_score) const
  {
    ppm_score = -1;
    {
      // Calculate the difference of the theoretical mass and the actually measured mass
      RangeMZ mz_range = DIAHelpers::createMZRangePPM(precursor_mz, dia_extract_window_, dia_extraction_ppm_);
      dia_integration_cache.begin(spectrum, im_range, dia_centroided_, true);
      const IntegratedDIAWindow& window = dia_integration_cache.integrate(spectrum, mz_range, im_range, dia_centroided_);

      // Catch if no signal was found and replace it with the most extreme
      // value. Otherwise, calculate the difference in ppm.
      if (!window.found)
      {
        ppm_score = Math::getPPMAbs(precursor_mz + mz_range.getSpan(), precursor_mz);
        return false;
      }
      else
      {
        ppm_score = Math::getPPMAbs(window.mz, precursor_mz);
        return true;
      }
    }
  }

  /// Precursor isotope scores
  void DIAScoring::dia_ms1_isotope_scores(double precursor_mz, const std::vector<SpectrumPtrType>& spectrum,
                                          RangeMobility& im_range, double& isotope_corr, double& isotope_overlap, const EmpiricalFormula& sum_formula) const
  {
    dia_integration_cache.begin(spectrum, im_range, dia_centroided_, true);
    // although precursor_mz can be received from the empirical formula (if non-empty), the actual precursor could be
    // slightly different. And also for compounds, usually the neutral sum_formula without adducts is given.
    // Therefore calculate the isotopes based on the formula but place them at precursor_mz
    dias_pool.exp_isotopes_int.clear();
    getIsotopeIntysFromExpSpec_(precursor_mz, spectrum, sum_formula.getCharge(), im_range, dias_pool.exp_isotopes_int);
    double max_ratio = 0;
    int nr_occurrences = 0;

    // calculate the scores:
    // isotope correlation (forward) and the isotope overlap (backward) scores
    isotope_corr = scoreIsotopePattern_(dias_pool.exp_isotopes_int, sum_formula);
    largePeaksBeforeFirstIsotope_(spectrum, precursor_mz, dias_pool.exp_isotopes_int[0], nr_occurrences, max_ratio, im_range);
    isotope_overlap = max_ratio;
  }

  void DIAScoring::getIsotopeIntysFromExpSpec_(double precursor_mz, const SpectrumSequence& spectrum, int charge_state, const RangeMobility& im_range,
                            std::vector<double>& isotopes_int) const
  {
    double abs_charge = std::fabs(static_cast<double>(charge_state));
    for (int iso = 0; iso <= dia_nr_isotopes_; ++iso)
    {
      RangeMZ mz_range = DIAHelpers::createMZRangePPM(precursor_mz + iso * C13C12_MASSDIFF_U / abs_charge, dia_extract_window_, dia_extraction_ppm_);
      const IntegratedDIAWindow& window = dia_integration_cache.integrate(spectrum, mz_range, im_range, dia_centroided_);
      isotopes_int.push_back(window.intensity);
    }
  }

  void DIAScoring::dia_ms1_isotope_scores_averagine(double precursor_mz, const SpectrumSequence& spectrum, int charge_state, RangeMobility& im_range,
                                                    double& isotope_corr, double& isotope_overlap) const
  {
    dia_integration_cache.begin(spectrum, im_range, dia_centroided_, true);
    dias_pool.exp_isotopes_int.clear();
    getIsotopeIntysFromExpSpec_(precursor_mz, spectrum, charge_state, im_range, dias_pool.exp_isotopes_int);
    // NOTE: this is a rough estimate of the neutral mz value since we would not know the charge carrier for negative ions
    const IsotopeDistribution& isotope_dist = getCachedAveragineIsotopeDistribution(
      std::fabs(precursor_mz * charge_state),
      static_cast<int>(dia_nr_isotopes_) + 1);

    double max_ratio;
    int nr_occurrences;
    // calculate the scores:
    // isotope correlation (forward) and the isotope overlap (backward) scores
    isotope_corr = scoreIsotopePattern_(dias_pool.exp_isotopes_int, isotope_dist);
    largePeaksBeforeFirstIsotope_(spectrum, precursor_mz, dias_pool.exp_isotopes_int[0], nr_occurrences, max_ratio, im_range);
    isotope_overlap = max_ratio;
  }

  void DIAScoring::dia_by_ion_score(const SpectrumSequence& spectrum,
                                    AASequence& sequence, int charge, const RangeMobility& im_range, double& bseries_score,
                                    double& yseries_score) const
  {
    bseries_score = 0;
    yseries_score = 0;
    OPENMS_PRECONDITION(charge > 0, "Charge is a positive integer"); // for peptides, charge should be positive

    const std::string sequence_key = sequence.toString();
    if (!dias_pool.byseries_valid ||
        dias_pool.byseries_charge != charge ||
        dias_pool.byseries_sequence != sequence_key)
    {
      dias_pool.bseries.clear();
      dias_pool.yseries.clear();
      OpenMS::DIAHelpers::getBYSeries(sequence, dias_pool.bseries, dias_pool.yseries, generator, charge);
      dias_pool.byseries_sequence = sequence_key;
      dias_pool.byseries_charge = charge;
      dias_pool.byseries_valid = true;
    }
    for (const auto& b_ion_mz : dias_pool.bseries)
    {
      RangeMZ mz_range = DIAHelpers::createMZRangePPM(b_ion_mz, dia_extract_window_, dia_extraction_ppm_);

      const IntegratedDIAWindow& window = dia_integration_cache.integrate(spectrum, mz_range, im_range, dia_centroided_);
      double ppmdiff = Math::getPPMAbs(window.mz, b_ion_mz);
      if (window.found && ppmdiff < dia_byseries_ppm_diff_ && window.intensity > dia_byseries_intensity_min_)
      {
        bseries_score++;
      }
    }
    for (const auto& y_ion_mz : dias_pool.yseries)
    {
      RangeMZ mz_range = DIAHelpers::createMZRangePPM(y_ion_mz, dia_extract_window_, dia_extraction_ppm_);

      const IntegratedDIAWindow& window = dia_integration_cache.integrate(spectrum, mz_range, im_range, dia_centroided_);
      double ppmdiff = Math::getPPMAbs(window.mz, y_ion_mz);
      if (window.found && ppmdiff < dia_byseries_ppm_diff_ && window.intensity > dia_byseries_intensity_min_)
      {
        yseries_score++;
      }
    }
  }

  void DIAScoring::score_with_isotopes(SpectrumSequence& spectrum, const std::vector<TransitionType>& transitions, const RangeMobility& im_range, double& dotprod, double& manhattan) const
  {
    DiaPrescoreThreadCache& cache = dia_prescore_cache;
    if (!diaPrescoreCacheMatches(cache, transitions, dia_extract_window_, dia_nr_isotopes_, dia_nr_charges_))
    {
      OpenMS::DiaPrescore dp(dia_extract_window_, dia_nr_isotopes_, dia_nr_charges_);
      dp.prepare(transitions, cache.prepared);
      updateDiaPrescoreCacheKey(cache, transitions, dia_extract_window_, dia_nr_isotopes_, dia_nr_charges_);
    }
    scorePreparedCached(spectrum, cache.prepared, dia_extract_window_, im_range, dotprod, manhattan);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Private methods

  void DIAScoring::diaIsotopeScoresSub_(const std::vector<TransitionType>& transitions, const SpectrumSequence& spectrum,
                                        OpenSwath::IMRMFeature* mrmfeature,
                                        const RangeMobility& im_range,
                                        double& isotope_corr,
                                        double& isotope_overlap) const
  {
    dias_pool.isotopes_int.clear();
    std::vector<double>& isotopes_int = dias_pool.isotopes_int;
    double max_ratio;
    int nr_occurences;
    const double feature_intensity = mrmfeature->getIntensity();
    const auto* openms_feature = dynamic_cast<const OpenMS::MRMFeatureOpenMS*>(mrmfeature);
    const bool use_openms_index = openms_feature != nullptr && openms_feature->hasCachedFeatureIds();
    for (Size k = 0; k < transitions.size(); k++)
    {
      isotopes_int.clear();
      const TransitionType& transition = transitions[k];
      double transition_intensity = 0.0;
      if (use_openms_index)
      {
        transition_intensity = openms_feature->getFeatureIntensity(k);
      }
      else if (openms_feature != nullptr)
      {
        transition_intensity = openms_feature->getFeatureIntensity(transition.transition_name, k);
      }
      else
      {
        transition_intensity = mrmfeature->getFeature(transition.transition_name)->getIntensity();
      }
      const double rel_intensity = transition_intensity / feature_intensity;

      // If no charge is given, we assume it to be 1
      int putative_fragment_charge = 1;
      if (transition.fragment_charge != 0)
      {
        putative_fragment_charge = transition.fragment_charge;
      }

      // collect the potential isotopes of this peak
      double abs_charge = std::fabs(static_cast<double>(putative_fragment_charge));
      for (int iso = 0; iso <= dia_nr_isotopes_; ++iso)
      {
        RangeMZ mz_range = DIAHelpers::createMZRangePPM(transition.getProductMZ() + iso * C13C12_MASSDIFF_U / abs_charge, dia_extract_window_, dia_extraction_ppm_);
        const IntegratedDIAWindow& window = dia_integration_cache.integrate(spectrum, mz_range, im_range, dia_centroided_);
        isotopes_int.push_back(window.intensity);
      }

      // calculate the scores:
      // isotope correlation (forward) and the isotope overlap (backward) scores
      double score = scoreIsotopePattern_(isotopes_int, transition.getProductMZ(), putative_fragment_charge);
      isotope_corr += score * rel_intensity;
      largePeaksBeforeFirstIsotope_(spectrum, transition.getProductMZ(), isotopes_int[0], nr_occurences, max_ratio, im_range);
      isotope_overlap += nr_occurences * rel_intensity;
    }
  }

  void DIAScoring::largePeaksBeforeFirstIsotope_(const SpectrumSequence& spectrum, double mono_mz, double mono_int, int& nr_occurences, double& max_ratio, const RangeMobility& im_range) const
  {
    nr_occurences = 0;
    max_ratio = 0.0;

    for (int ch = 1; ch <= dia_nr_charges_; ++ch)
    {
      double center = mono_mz - C13C12_MASSDIFF_U / (double) ch;
      RangeMZ mz_range = DIAHelpers::createMZRangePPM(center, dia_extract_window_, dia_extraction_ppm_);

      const IntegratedDIAWindow& window = dia_integration_cache.integrate(spectrum, mz_range, im_range, dia_centroided_);
      // Continue if no signal was found - we therefore don't make a statement
      // about the mass difference if no signal is present.
      if (!window.found)
      {
        continue;
      }

      // Compute ratio between the (presumed) monoisotopic peak intensity and the now found peak
      double ratio;
      if (mono_int != 0)
      {
        ratio = window.intensity / mono_int;
      }
      else
      {
        ratio = 0;
      }
      if (ratio > max_ratio) {max_ratio = ratio;}

      double ddiff_ppm = std::fabs(window.mz - center) * 1e6 / center;

      // FEATURE we should fit a theoretical distribution to see whether we really are a secondary peak
      if (ratio > 1 && ddiff_ppm < peak_before_mono_max_ppm_diff_)
      {
        //isotope_overlap += 1.0 * rel_intensity;

        nr_occurences += 1; // we count how often this happens...

#ifdef MRMSCORING_TESTING
        std::cout << " _ overlap diff ppm  " << ddiff_ppm << " and inten ratio " << ratio << " with " << mono_int << '\n';
#endif
      }
    }
  }

  double DIAScoring::scoreIsotopePattern_(const std::vector<double>& isotopes_int,
                                          double product_mz,
                                          int putative_fragment_charge) const
  {
    OPENMS_PRECONDITION(putative_fragment_charge != 0, "Charge needs to be set to != 0"); // charge can be positive and negative

    // create the theoretical distribution from the peptide weight
    // NOTE: this is a rough estimate of the neutral mz value since we would not know the charge carrier for negative ions
    const IsotopeDistribution& isotope_dist = getCachedAveragineIsotopeDistribution(
      std::fabs(product_mz * putative_fragment_charge),
      static_cast<int>(dia_nr_isotopes_) + 1);

    return scoreIsotopePattern_(isotopes_int, isotope_dist);
  } //end of dia_isotope_corr_sub

  double DIAScoring::scoreIsotopePattern_(const std::vector<double>& isotopes_int,
                                          const EmpiricalFormula& empf) const
  {
    return scoreIsotopePattern_(isotopes_int,
                                empf.getIsotopeDistribution(CoarseIsotopePatternGenerator(dia_nr_isotopes_ + 1)));
  }

  double DIAScoring::scoreIsotopePattern_(const std::vector<double>& isotopes_int,
                                          const IsotopeDistribution& isotope_dist) const
  {
    typedef OpenMS::FeatureFinderAlgorithmPickedHelperStructs::TheoreticalIsotopePattern TheoreticalIsotopePattern;

    TheoreticalIsotopePattern isotopes;

    for (IsotopeDistribution::ConstIterator it = isotope_dist.begin(); it != isotope_dist.end(); ++it)
    {
      isotopes.intensity.push_back(it->getIntensity());
    }
    isotopes.optional_begin = 0;
    isotopes.optional_end = dia_nr_isotopes_;

    // scale the distribution to a maximum of 1
    double max = 0.0;
    for (Size i = 0; i < isotopes.intensity.size(); ++i)
    {
      if (isotopes.intensity[i] > max)
      {
        max = isotopes.intensity[i];
      }
    }
    isotopes.max = max;
    if (max == 0.) max = 1.;
    for (Size i = 0; i < isotopes.intensity.size(); ++i)
    {
      isotopes.intensity[i] /= max;
    }
    isotopes.trimmed_left = 0;

    // score the pattern against a theoretical one
    OPENMS_POSTCONDITION(isotopes_int.size() == isotopes.intensity.size(), "Vectors for pearson correlation do not have the same size.");
    double int_score = OpenMS::Math::pearsonCorrelationCoefficient(isotopes_int.begin(), isotopes_int.end(), isotopes.intensity.begin(), isotopes.intensity.end());
    if (std::isnan(int_score))
    {
      int_score = 0;
    }
    return int_score;
  } //end of dia_isotope_corr_sub
}
