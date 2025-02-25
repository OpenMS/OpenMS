// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/FEATUREFINDER/MultiplexClustering.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/ML/CLUSTERING/GridBasedCluster.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>
#include <OpenMS/PROCESSING/MISC/SplinePackage.h>
#include <OpenMS/PROCESSING/MISC/SplineInterpolatedPeaks.h>

#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/ML/REGRESSION/LinearRegressionWithoutIntercept.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <vector>
#include <numeric>
#include <fstream>
#include <iostream>
#include <ostream>
#include <algorithm>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/classification.hpp>

//#define DEBUG

using namespace std;

namespace OpenMS
{
  FeatureFinderMultiplexAlgorithm::FeatureFinderMultiplexAlgorithm() :
    DefaultParamHandler("FeatureFinderMultiplexAlgorithm")
  {
    // parameter section: algorithm
    defaults_.setValue("algorithm:labels", "[][Lys8,Arg10]", "Labels used for labelling the samples. If the sample is unlabelled (i.e. you want to detect only single peptide features) please leave this parameter empty. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL");
    defaults_.setValue("algorithm:charge", "1:4", "Range of charge states in the sample, i.e. min charge : max charge.");
    defaults_.setValue("algorithm:isotopes_per_peptide", "3:6", "Range of isotopes per peptide in the sample. For example 3:6, if isotopic peptide patterns in the sample consist of either three, four, five or six isotopic peaks. ", {"advanced"});
    defaults_.setValue("algorithm:rt_typical", 40.0, "Typical retention time [s] over which a characteristic peptide elutes. (This is not an upper bound. Peptides that elute for longer will be reported.)");
    defaults_.setMinFloat("algorithm:rt_typical", 0.0);
    defaults_.setValue("algorithm:rt_band", 0.0, "The algorithm searches for characteristic isotopic peak patterns, spectrum by spectrum. For some low-intensity peptides, an important peak might be missing in one spectrum but be present in one of the neighbouring ones. The algorithm takes a bundle of neighbouring spectra with width rt_band into account. For example with rt_band = 0, all characteristic isotopic peaks have to be present in one and the same spectrum. As rt_band increases, the sensitivity of the algorithm but also the likelihood of false detections increases.");
    defaults_.setMinFloat("algorithm:rt_band", 0.0);
    defaults_.setValue("algorithm:rt_min", 2.0, "Lower bound for the retention time [s]. (Any peptides seen for a shorter time period are not reported.)");
    defaults_.setMinFloat("algorithm:rt_min", 0.0);
    defaults_.setValue("algorithm:mz_tolerance", 6.0, "m/z tolerance for search of peak patterns.");
    defaults_.setMinFloat("algorithm:mz_tolerance", 0.0);
    defaults_.setValue("algorithm:mz_unit", "ppm", "Unit of the 'mz_tolerance' parameter.");
    defaults_.setValidStrings("algorithm:mz_unit", {"Da","ppm"});
    defaults_.setValue("algorithm:intensity_cutoff", 1000.0, "Lower bound for the intensity of isotopic peaks.");
    defaults_.setMinFloat("algorithm:intensity_cutoff", 0.0);
    defaults_.setValue("algorithm:peptide_similarity", 0.5, "Two peptides in a multiplet are expected to have the same isotopic pattern. This parameter is a lower bound on their similarity.");
    defaults_.setMinFloat("algorithm:peptide_similarity", -1.0);
    defaults_.setMaxFloat("algorithm:peptide_similarity", 1.0);
    defaults_.setValue("algorithm:averagine_similarity", 0.4, "The isotopic pattern of a peptide should resemble the averagine model at this m/z position. This parameter is a lower bound on similarity between measured isotopic pattern and the averagine model.");
    defaults_.setMinFloat("algorithm:averagine_similarity", -1.0);
    defaults_.setMaxFloat("algorithm:averagine_similarity", 1.0);
    defaults_.setValue("algorithm:averagine_similarity_scaling", 0.95, "Let x denote this scaling factor, and p the averagine similarity parameter. For the detection of single peptides, the averagine parameter p is replaced by p' = p + x(1-p), i.e. x = 0 -> p' = p and x = 1 -> p' = 1. (For knock_out = true, peptide doublets and singlets are detected simultaneously. For singlets, the peptide similarity filter is irreleavant. In order to compensate for this 'missing filter', the averagine parameter p is replaced by the more restrictive p' when searching for singlets.)", {"advanced"});
    defaults_.setMinFloat("algorithm:averagine_similarity_scaling", 0.0);
    defaults_.setMaxFloat("algorithm:averagine_similarity_scaling", 1.0);
    defaults_.setValue("algorithm:missed_cleavages", 0, "Maximum number of missed cleavages due to incomplete digestion. (Only relevant if enzymatic cutting site coincides with labelling site. For example, Arg/Lys in the case of trypsin digestion and SILAC labelling.)");
    defaults_.setMinInt("algorithm:missed_cleavages", 0);
    defaults_.setValue("algorithm:spectrum_type", "automatic", "Type of MS1 spectra in input mzML file. 'automatic' determines the spectrum type directly from the input mzML file.", {"advanced"});
    defaults_.setValidStrings("algorithm:spectrum_type", {"profile","centroid","automatic"});
    defaults_.setValue("algorithm:averagine_type","peptide","The type of averagine to use, currently RNA, DNA or peptide", {"advanced"});
    defaults_.setValidStrings("algorithm:averagine_type", {"peptide","RNA","DNA"});
    defaults_.setValue("algorithm:knock_out", "false", "Is it likely that knock-outs are present? (Supported for doublex, triplex and quadruplex experiments only.)", {"advanced"});
    defaults_.setValidStrings("algorithm:knock_out", {"true","false"});

    defaults_.setSectionDescription("algorithm", "algorithmic parameters");

    // parameter section: labels
    defaults_.setValue("labels:Arg6", 12.08, "description", {"advanced"});

    defaults_.setSectionDescription("labels", "mass shifts for all possible labels");

    MultiplexDeltaMassesGenerator generator;
    const Param& p = generator.getParameters();
    for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
    {
      String label_name = "labels:";
      label_name += it->name;

      defaults_.setValue(label_name, it->value, it->description, {"advanced"});
      defaults_.setMinFloat(label_name, 0.0);

      label_mass_shift_.insert(make_pair(it->name, it->value));
    }

    // parameter section: algorithm, get selected charge range
    String charge_string = defaults_.getValue("algorithm:charge").toString();
    charge_min_ = charge_string.prefix(':').toInt();
    charge_max_ = charge_string.suffix(':').toInt();
    if (charge_min_ > charge_max_)
    {
      swap(charge_min_, charge_max_);
    }

    // parameter section: algorithm, get isotopes per peptide range
    String isotopes_per_peptide_string = defaults_.getValue("algorithm:isotopes_per_peptide").toString();
    isotopes_per_peptide_min_ = isotopes_per_peptide_string.prefix(':').toInt();
    isotopes_per_peptide_max_ = isotopes_per_peptide_string.suffix(':').toInt();
    if (isotopes_per_peptide_min_ > isotopes_per_peptide_max_)
    {
      swap(isotopes_per_peptide_min_, isotopes_per_peptide_max_);
    }

    centroided_ = false;
  }

  /**
   * @brief order of charge states
   *
   * 2+ 3+ 4+ 1+ 5+ 6+ ...
   *
   * Order charge states by the likelihood of their occurrence, i.e. we search for the most likely charge states first.
   */
  static size_t orderCharge(int charge)
  {
    if ((1 < charge) && (charge < 5))
    {
      return (charge - 1);
    }
    else if (charge == 1)
    {
      return 4;
    }
    else
    {
      return charge;
    }
  }

  /**
   * @brief comparator of peak patterns
   *
   * The comperator determines in which order the peak patterns are searched for.
   * First we check the number of mass shifts (triplets before doublets before singlets).
   * Then we check the first mass shift (for example 6 Da before 12 Da i.e. misscleavage).
   * Finally we check for charges (2+ before 1+, most likely first).
   *
   * @param pattern1    first peak pattern
   * @param pattern2    second peak pattern
   *
   * @return true if pattern1 should be searched before pattern2
   */
  static bool lessPattern(const MultiplexIsotopicPeakPattern& pattern1, const MultiplexIsotopicPeakPattern& pattern2)
  {
    if (pattern1.getMassShiftCount() == pattern2.getMassShiftCount())
    {
      // The first mass shift is by definition always zero.
      if ((pattern1.getMassShiftCount() > 1) && (pattern2.getMassShiftCount() > 1))
      {
        if (pattern1.getMassShiftAt(1) == pattern2.getMassShiftAt(1))
        {
          // 2+ before 3+ before 4+ before 1+ before 5+ before 6+ etc.
          return orderCharge(pattern1.getCharge()) < orderCharge(pattern2.getCharge());
        }
        else
        {
          return pattern1.getMassShiftAt(1) < pattern2.getMassShiftAt(1);
        }
      }
      else
      {
        // 2+ before 3+ before 4+ before 1+ before 5+ before 6+ etc.
        return orderCharge(pattern1.getCharge()) < orderCharge(pattern2.getCharge());
      }
    }
    else
    {
      // triplets before doublets before singlets
      return pattern1.getMassShiftCount() > pattern2.getMassShiftCount();
    }
  }

  std::vector<MultiplexIsotopicPeakPattern> FeatureFinderMultiplexAlgorithm::generatePeakPatterns_(int charge_min, int charge_max, int peaks_per_peptide_max, const std::vector<MultiplexDeltaMasses>& mass_pattern_list)
  {
    std::vector<MultiplexIsotopicPeakPattern> list;

    // iterate over all charge states
    for (int c = charge_max; c >= charge_min; --c)
    {
      // iterate over all mass shifts
      for (unsigned i = 0; i < mass_pattern_list.size(); ++i)
      {
        MultiplexIsotopicPeakPattern pattern(c, peaks_per_peptide_max, mass_pattern_list[i], i);
        list.push_back(pattern);
      }
    }

    sort(list.begin(), list.end(), lessPattern);

    return list;
  }

  void FeatureFinderMultiplexAlgorithm::correctPeptideIntensities_(const MultiplexIsotopicPeakPattern& pattern, std::map<size_t, SplinePackage>& spline_chromatograms, const std::vector<double>& rt_peptide, std::vector<double>& intensity_peptide) const
  {
    // determine ratios through linear regression
    // (In most labelled mass spectrometry experiments, the fold change i.e. ratio and not the individual peptide intensities
    // are of primary interest. For that reason, we determine the ratios from interpolated chromatogram data points directly,
    // and then correct the current ones.)

    std::vector<double> ratios;    // light/light, medium/light, heavy/light etc.
    ratios.push_back(1.0);
    // loop over peptides
    for (size_t peptide = 1; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      std::vector<double> intensities1;
      std::vector<double> intensities2;

      // loop over isotopes
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {

        // find splines for the mass traces of the lightest and other peptide
        size_t idx_1 = isotope;
        size_t idx_2 = peptide * isotopes_per_peptide_max_ + isotope;
        if ((spline_chromatograms.find(idx_1) == spline_chromatograms.end()) || (spline_chromatograms.find(idx_2) == spline_chromatograms.end()))
        {
           continue;
        }

        std::map<size_t, SplinePackage>::iterator it1 = spline_chromatograms.find(idx_1);
        std::map<size_t, SplinePackage>::iterator it2 = spline_chromatograms.find(idx_2);

        double rt_min = std::min(it1->second.getPosMin(), it2->second.getPosMin());
        double rt_max = std::max(it1->second.getPosMax(), it2->second.getPosMax());
        double rt_step_width = 0.7 * std::min(it1->second.getPosStepWidth(), it2->second.getPosStepWidth());

        for (double rt = rt_min; rt < rt_max; rt += rt_step_width)
        {
          double intensity1 = it1->second.eval(rt);
          double intensity2 = it2->second.eval(rt + rt_peptide[peptide] - rt_peptide[0]);    // Take RT shifts between peptide into account to find corresponding intensities.

          // Use only if we land within both chromatograms i.e. non-zero intensities.
          if ((intensity1 > 0) && (intensity2 > 0))
          {
            intensities1.push_back(intensity1);
            intensities2.push_back(intensity2);
          }
        }

      }

      // We require at least five data points for a reliable linear interpolation.
      if (intensities1.size() > 5)
      {
        OpenMS::Math::LinearRegressionWithoutIntercept linreg;
        linreg.addData(intensities1, intensities2);

        ratios.push_back(linreg.getSlope());
      }
      else
      {
        ratios.push_back(-1.0);
      }

    }

    // correction for doublets
    if ((pattern.getMassShiftCount() == 2) && (ratios[1] > 0))
    {
      double intensity0 = (intensity_peptide[0] + ratios[1] * intensity_peptide[1]) / (1 + ratios[1] * ratios[1]);
      double intensity1 = ratios[1] * intensity0;

      intensity_peptide[0] = intensity0;
      intensity_peptide[1] = intensity1;
    }
    // correction for triplets or higher multiplets
    else if ((pattern.getMassShiftCount() > 2))
    {
      for (size_t peptide = 1; peptide < pattern.getMassShiftCount(); ++peptide)
      {
        if (ratios[peptide] > 0)
        {
          intensity_peptide[peptide] = ratios[peptide] * intensity_peptide[0];
        }
      }
    }

  }

  std::vector<double> FeatureFinderMultiplexAlgorithm::determinePeptideIntensitiesCentroided_(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteCentroided >& satellites)
  {
    // determine peptide intensities and RT shift between the peptides
    // i.e. first determine the RT centre of mass for each peptide
    std::vector<double> rt_peptide;
    std::vector<double> intensity_peptide;

    std::map<size_t, SplinePackage> spline_chromatograms;

    PeakIntegrator pi;
    Param param = pi.getDefaults();
    param.setValue("integration_type","trapezoid");    // intensity_sum, simpson, trapezoid (Note that 'simpson' may lead to negative area-under-the-curve.)
    pi.setParameters(param);

    // loop over peptides
    for (size_t peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // coordinates of the peptide feature
      // RT is the intensity-average of all satellites peaks of all (!) mass traces
      double rt(0);
      double rt_min(0);
      double rt_max(0);
      double intensity_sum_simple(0);    // for intensity-averaged rt
      double intensity_sum(0);

      // loop over isotopes i.e. mass traces of the peptide
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {
        // find satellites for this isotope i.e. mass trace
        size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
        std::pair<std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator, std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator> satellites_isotope;
        satellites_isotope = satellites.equal_range(idx);

        MSChromatogram chromatogram;
        // loop over satellites for this isotope i.e. mass trace
        for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it = satellites_isotope.first; satellite_it != satellites_isotope.second; ++satellite_it)
        {
          // find indices of the peak
          size_t rt_idx = (satellite_it->second).getRTidx();
          size_t mz_idx = (satellite_it->second).getMZidx();

          // find peak itself
          MSExperiment::ConstIterator it_rt = exp_centroid_.begin();
          std::advance(it_rt, rt_idx);
          MSSpectrum::ConstIterator it_mz = it_rt->begin();
          std::advance(it_mz, mz_idx);

          double rt_temp = it_rt->getRT();
          double intensity_temp = it_mz->getIntensity();

          if ((peptide + isotope == 0) || (rt_temp < rt_min))
          {
            rt_min = rt_temp;
          }

          if ((peptide + isotope == 0) || (rt_temp > rt_max))
          {
            rt_max = rt_temp;
          }

          rt += rt_temp * intensity_temp;
          intensity_sum_simple += intensity_temp;

          chromatogram.push_back(ChromatogramPeak(rt_temp, intensity_temp));
        }

        chromatogram.sortByPosition();

        // construct spline interpolations for later use
        // Reliable spline interpolation only for 5 or more data points in chromatogram.
        if (chromatogram.size() > 5)
        {
          std::vector<double> rt;
          std::vector<double> intensity;
          for (const auto &it : chromatogram)
          {
            rt.push_back(it.getRT());
            intensity.push_back(it.getIntensity());
          }
          spline_chromatograms.insert(std::make_pair(idx, SplinePackage(rt, intensity)));
        }

        if (chromatogram.size() > 1)
        {
          double rt_start = chromatogram.begin()->getPos();
          double rt_end = chromatogram.back().getPos();

          PeakIntegrator::PeakArea pa = pi.integratePeak(chromatogram, rt_start, rt_end);
          intensity_sum += pa.area;
        }
        else if (chromatogram.size() == 1)
        {
          intensity_sum += chromatogram.begin()->getIntensity();
        }
      }

      rt /= intensity_sum_simple;
      rt_peptide.push_back(rt);
      if (intensity_sum == 0 || (rt_max - rt_min < static_cast<double>(param_.getValue("algorithm:rt_min"))))
      {
        intensity_sum = -1.0;
      }
      intensity_peptide.push_back(intensity_sum);
    }

    // If any of the peptide intensities could not be determined (i.e. -1) then there is no need for further corrections.
    if (std::find(intensity_peptide.begin(), intensity_peptide.end(), -1.0) != intensity_peptide.end())
    {
      return intensity_peptide;
    }

    // If the pattern searched for peptide singlets, then there are no further corrections possible.
    if (pattern.getMassShiftCount() < 2)
    {
      return intensity_peptide;
    }

    correctPeptideIntensities_(pattern, spline_chromatograms, rt_peptide, intensity_peptide);

    return intensity_peptide;
  }

  std::vector<double> FeatureFinderMultiplexAlgorithm::determinePeptideIntensitiesProfile_(const MultiplexIsotopicPeakPattern& pattern, const std::multimap<size_t, MultiplexSatelliteProfile >& satellites)
  {
    // determine peptide intensities and RT shift between the peptides
    // i.e. first determine the RT centre of mass for each peptide
    std::vector<double> rt_peptide;
    std::vector<double> intensity_peptide;

    std::map<size_t, SplinePackage> spline_chromatograms;

    PeakIntegrator pi;
    Param param = pi.getDefaults();
    param.setValue("integration_type","trapezoid");    // intensity_sum, simpson, trapezoid (Note that 'simpson' may lead to negative area-under-the-curve.)
    pi.setParameters(param);

    // loop over peptides
    for (size_t peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
    {
      // coordinates of the peptide feature
      // RT is the intensity-average of all satellites peaks of all (!) mass traces
      double rt(0);
      double rt_min(0);
      double rt_max(0);
      double intensity_sum_simple(0);    // for intensity-averaged rt
      double intensity_sum(0);

      // loop over isotopes i.e. mass traces of the peptide
      for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
      {
        // find satellites for this isotope i.e. mass trace
        size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
        std::pair<std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator, std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator> satellites_isotope;
        satellites_isotope = satellites.equal_range(idx);

        MSChromatogram chromatogram;
        // loop over satellites for this isotope i.e. mass trace
        for (std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator satellite_it = satellites_isotope.first; satellite_it != satellites_isotope.second; ++satellite_it)
        {
          double rt_temp = (satellite_it->second).getRT();
          double intensity_temp = (satellite_it->second).getIntensity();

          if ((peptide + isotope == 0) || (rt_temp < rt_min))
          {
            rt_min = rt_temp;
          }

          if ((peptide + isotope == 0) || (rt_temp > rt_max))
          {
            rt_max = rt_temp;
          }

          rt += rt_temp * intensity_temp;
          intensity_sum_simple += intensity_temp;

          chromatogram.push_back(ChromatogramPeak(rt_temp, intensity_temp));
        }

        makePeakPositionUnique(chromatogram, IntensityAveragingMethod::MEDIAN);

        // construct spline interpolations for later use
        // Reliable spline interpolation only for 5 or more data points in chromatogram.
        if (chromatogram.size() > 5)
        {
          std::vector<double> rt;
          std::vector<double> intensity;
          for (const auto &it : chromatogram)
          {
            rt.push_back(it.getRT());
            intensity.push_back(it.getIntensity());
          }
          spline_chromatograms.insert(std::make_pair(idx, SplinePackage(rt, intensity)));
        }

        if (chromatogram.size() > 1)
        {
          // Positions are already sorted in makePeakPositionUnique(), i.e. sortByPosition() not necessary.
          double rt_start = chromatogram.begin()->getPos();
          double rt_end = chromatogram.back().getPos();

          PeakIntegrator::PeakArea pa = pi.integratePeak(chromatogram, rt_start, rt_end);
          intensity_sum += pa.area;
        }
        else if (chromatogram.size() == 1)
        {
          intensity_sum += chromatogram.begin()->getIntensity();
        }

      }

      rt /= intensity_sum_simple;
      rt_peptide.push_back(rt);
      if (intensity_sum == 0 || (rt_max - rt_min < static_cast<double>(param_.getValue("algorithm:rt_min"))))
      {
        intensity_sum = -1.0;
      }
      intensity_peptide.push_back(intensity_sum);
    }

    // If any of the peptide intensities could not be determined (i.e. -1) then there is no need for further corrections.
    if (std::find(intensity_peptide.begin(), intensity_peptide.end(), -1.0) != intensity_peptide.end())
    {
      return intensity_peptide;
    }

    // If the pattern searched for peptide singlets, then there are no further corrections possible.
    if (pattern.getMassShiftCount() < 2)
    {
      return intensity_peptide;
    }

    correctPeptideIntensities_(pattern, spline_chromatograms, rt_peptide, intensity_peptide);

    return intensity_peptide;
  }

  void FeatureFinderMultiplexAlgorithm::generateMapsCentroided_(const std::vector<MultiplexIsotopicPeakPattern>& patterns, const std::vector<MultiplexFilteredMSExperiment>& filter_results, std::vector<std::map<int, GridBasedCluster> >& cluster_results)
  {
    // loop over peak patterns
    for (unsigned pattern = 0; pattern < patterns.size(); ++pattern)
    {
      // loop over clusters
      for (std::map<int, GridBasedCluster>::const_iterator cluster_it = cluster_results[pattern].begin(); cluster_it != cluster_results[pattern].end(); ++cluster_it)
      {
        GridBasedCluster cluster = cluster_it->second;
        std::vector<int> points = cluster.getPoints();

        // Construct a satellite set for the complete peptide multiplet
        // Make sure there are no duplicates, i.e. the same satellite from different filtered peaks.
        std::multimap<size_t, MultiplexSatelliteCentroided > satellites;
        // loop over points in cluster
        for (std::vector<int>::const_iterator point_it = points.begin(); point_it != points.end(); ++point_it)
        {
          MultiplexFilteredPeak peak = filter_results[pattern].getPeak(*point_it);
          // loop over satellites of the peak
          for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it = peak.getSatellites().begin(); satellite_it != peak.getSatellites().end(); ++satellite_it)
          {
            // check if this satellite (i.e. these indices) are already in the set
            bool satellite_in_set = false;
            for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it_2 = satellites.begin(); satellite_it_2 != satellites.end(); ++satellite_it_2)
            {
              if ((satellite_it_2->second.getRTidx() == satellite_it->second.getRTidx()) && (satellite_it_2->second.getMZidx() == satellite_it->second.getMZidx()))
              {
                satellite_in_set = true;
                break;
              }
            }
            if (satellite_in_set)
            {
              break;
            }

            satellites.insert(std::make_pair(satellite_it->first, MultiplexSatelliteCentroided(satellite_it->second.getRTidx(), satellite_it->second.getMZidx())));
          }
        }

        // determine peptide intensities
        std::vector<double> peptide_intensities = determinePeptideIntensitiesCentroided_(patterns[pattern], satellites);

        // If no reliable peptide intensity can be determined, we do not report the peptide multiplet.
        if (std::find(peptide_intensities.begin(), peptide_intensities.end(), -1.0) != peptide_intensities.end())
        {
          continue;
        }

        std::vector<Feature> features;
        ConsensusFeature consensus;
        bool abort = false;

        // construct the feature and consensus maps
        // loop over peptides
        for (size_t peptide = 0; (peptide < patterns[pattern].getMassShiftCount() && !abort); ++peptide)
        {
          // coordinates of the peptide feature
          // RT is the intensity-average of all satellites peaks of the mono-isotopic mass trace
          // m/z is the intensity-average of all satellites peaks of the mono-isotopic mass trace
          Feature feature;
          double rt(0);
          double mz(0);
          double intensity_sum(0);

          // loop over isotopes i.e. mass traces of the peptide
          for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
          {
            // find satellites for this isotope i.e. mass trace
            size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
            std::pair<std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator, std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator> satellites_isotope;
            satellites_isotope = satellites.equal_range(idx);

            DBoundingBox<2> mass_trace;

            // loop over satellites for this isotope i.e. mass trace
            for (std::multimap<size_t, MultiplexSatelliteCentroided >::const_iterator satellite_it = satellites_isotope.first; satellite_it != satellites_isotope.second; ++satellite_it)
            {
              // find indices of the peak
              size_t rt_idx = (satellite_it->second).getRTidx();
              size_t mz_idx = (satellite_it->second).getMZidx();

              // find peak itself
              MSExperiment::ConstIterator it_rt = exp_centroid_.begin();
              std::advance(it_rt, rt_idx);
              MSSpectrum::ConstIterator it_mz = it_rt->begin();
              std::advance(it_mz, mz_idx);

              if (isotope == 0)
              {
                rt += it_rt->getRT() * it_mz->getIntensity();
                mz += it_mz->getMZ() * it_mz->getIntensity();
                intensity_sum += it_mz->getIntensity();
              }

              mass_trace.enlarge(it_rt->getRT(), it_mz->getMZ());
            }

            if ((mass_trace.width() == 0) || (mass_trace.height() == 0))
            {
              // The mass trace contains only a single point. Add a small margin around
              // the point, otherwise the mass trace is considered empty and not drawn.
              // TODO: Remove the magic number for the margin.
              mass_trace.enlarge(mass_trace.minX() - 0.01, mass_trace.minY() - 0.01);
              mass_trace.enlarge(mass_trace.maxX() + 0.01, mass_trace.maxY() + 0.01);
            }

            if (!(mass_trace.isEmpty()))
            {
              ConvexHull2D hull;
              hull.addPoint(DPosition<2>(mass_trace.minX(), mass_trace.minY()));
              hull.addPoint(DPosition<2>(mass_trace.minX(), mass_trace.maxY()));
              hull.addPoint(DPosition<2>(mass_trace.maxX(), mass_trace.minY()));
              hull.addPoint(DPosition<2>(mass_trace.maxX(), mass_trace.maxY()));
              feature.getConvexHulls().push_back(hull);
            }
          }

          if (intensity_sum <= 0)
          {
            continue;
          }
          rt /= intensity_sum;
          mz /= intensity_sum;

          feature.setRT(rt);
          feature.setMZ(mz);
          feature.setIntensity(peptide_intensities[peptide]);
          feature.setCharge(patterns[pattern].getCharge());
          feature.setOverallQuality(1.0);

          // Check that the feature eluted long enough.
          // DBoundingBox<2> box = feature.getConvexHull().getBoundingBox();    // convex hull of the entire peptide feature
          DBoundingBox<2> box = feature.getConvexHulls()[0].getBoundingBox();    // convex hull of the mono-isotopic mass trace
          if (box.maxX() - box.minX() < static_cast<double>(param_.getValue("algorithm:rt_min")))
          {
            abort = true;
            break;
          }

          features.push_back(feature);

          if (peptide == 0)
          {
            // The first/lightest peptide acts as anchor of the peptide multiplet consensus.
            // All peptide feature handles are connected to this point.
            consensus.setRT(rt);
            consensus.setMZ(mz);
            consensus.setIntensity(peptide_intensities[peptide]);
            consensus.setCharge(patterns[pattern].getCharge());
            consensus.setQuality(1.0);
          }

          FeatureHandle feature_handle;
          feature_handle.setRT(rt);
          feature_handle.setMZ(mz);
          feature_handle.setIntensity(peptide_intensities[peptide]);
          feature_handle.setCharge(patterns[pattern].getCharge());
          feature_handle.setMapIndex(peptide);
          //feature_handle.setUniqueId(&UniqueIdInterface::setUniqueId);    // TODO: Do we need to set unique ID?
          consensus.insert(feature_handle);
          consensus_map_.getColumnHeaders()[peptide].size++;
        }

        if (!abort)
        {
          consensus_map_.push_back(consensus);
          for (Feature& feat : features)
          {
            feature_map_.push_back(feat);
          }
        }

      }

    }

  }

  void FeatureFinderMultiplexAlgorithm::generateMapsProfile_(const std::vector<MultiplexIsotopicPeakPattern>& patterns, const std::vector<MultiplexFilteredMSExperiment>& filter_results, const std::vector<std::map<int, GridBasedCluster> >& cluster_results)
  {
    // progress logger
    unsigned progress = 0;
    startProgress(0, patterns.size(), "constructing maps");

    // loop over peak patterns
    for (unsigned pattern = 0; pattern < patterns.size(); ++pattern)
    {
      setProgress(++progress);

      // loop over clusters
      for (std::map<int, GridBasedCluster>::const_iterator cluster_it = cluster_results[pattern].begin(); cluster_it != cluster_results[pattern].end(); ++cluster_it)
      {
        GridBasedCluster cluster = cluster_it->second;
        std::vector<int> points = cluster.getPoints();

        // Construct a satellite set for the complete peptide multiplet
        // Make sure there are no duplicates, i.e. the same satellite from different filtered peaks.
        std::multimap<size_t, MultiplexSatelliteProfile > satellites;
        // loop over points in cluster
        for (std::vector<int>::const_iterator point_it = points.begin(); point_it != points.end(); ++point_it)
        {
          MultiplexFilteredPeak peak = filter_results[pattern].getPeak(*point_it);
          // loop over satellites of the peak
          for (std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator satellite_it = peak.getSatellitesProfile().begin(); satellite_it != peak.getSatellitesProfile().end(); ++satellite_it)
          {
            satellites.insert(std::make_pair(satellite_it->first, MultiplexSatelliteProfile(satellite_it->second.getRT(), satellite_it->second.getMZ(), satellite_it->second.getIntensity())));
          }
        }

        // determine peptide intensities
        std::vector<double> peptide_intensities = determinePeptideIntensitiesProfile_(patterns[pattern], satellites);

        // If no reliable peptide intensity can be determined for one of the peptides, we do not report the peptide multiplet.
        if (std::find(peptide_intensities.begin(), peptide_intensities.end(), -1.0) != peptide_intensities.end())
        {
          continue;
        }

        std::vector<Feature> features;
        ConsensusFeature consensus;
        bool abort = false;

        // construct the feature and consensus maps
        // loop over peptides
        for (size_t peptide = 0; (peptide < patterns[pattern].getMassShiftCount() && !abort); ++peptide)
        {
          // coordinates of the peptide feature
          // RT is the intensity-average of all satellites peaks of the mono-isotopic mass trace
          // m/z is the intensity-average of all satellites peaks of the mono-isotopic mass trace
          Feature feature;
          double rt(0);
          double mz(0);
          double intensity_sum(0);

          // loop over isotopes i.e. mass traces of the peptide
          for (size_t isotope = 0; isotope < isotopes_per_peptide_max_; ++isotope)
          {
            // find satellites for this isotope i.e. mass trace
            size_t idx = peptide * isotopes_per_peptide_max_ + isotope;
            std::pair<std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator, std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator> satellites_isotope;
            satellites_isotope = satellites.equal_range(idx);

            DBoundingBox<2> mass_trace;

            // loop over satellites for this isotope i.e. mass trace
            for (std::multimap<size_t, MultiplexSatelliteProfile >::const_iterator satellite_it = satellites_isotope.first; satellite_it != satellites_isotope.second; ++satellite_it)
            {
              if (isotope == 0)
              {
                // Satellites of zero intensity makes sense (borders of peaks), but mess up feature/consensus construction.
                double intensity_temp = (satellite_it->second).getIntensity() + 0.0001;

                rt += (satellite_it->second).getRT() * intensity_temp;
                mz += (satellite_it->second).getMZ() * intensity_temp;
                intensity_sum += intensity_temp;
              }

              mass_trace.enlarge((satellite_it->second).getRT(), (satellite_it->second).getMZ());
            }

            if ((mass_trace.width() == 0) || (mass_trace.height() == 0))
            {
              // The mass trace contains only a single point. Add a small margin around
              // the point, otherwise the mass trace is considered empty and not drawn.
              // TODO: Remove the magic number for the margin.
              mass_trace.enlarge(mass_trace.minX() - 0.01, mass_trace.minY() - 0.01);
              mass_trace.enlarge(mass_trace.maxX() + 0.01, mass_trace.maxY() + 0.01);
            }

            if (!(mass_trace.isEmpty()))
            {
              ConvexHull2D hull;
              hull.addPoint(DPosition<2>(mass_trace.minX(), mass_trace.minY()));
              hull.addPoint(DPosition<2>(mass_trace.minX(), mass_trace.maxY()));
              hull.addPoint(DPosition<2>(mass_trace.maxX(), mass_trace.minY()));
              hull.addPoint(DPosition<2>(mass_trace.maxX(), mass_trace.maxY()));
              feature.getConvexHulls().push_back(hull);
            }
          }

          rt /= intensity_sum;
          mz /= intensity_sum;

          feature.setRT(rt);
          feature.setMZ(mz);
          feature.setIntensity(peptide_intensities[peptide]);
          feature.setCharge(patterns[pattern].getCharge());
          feature.setOverallQuality(1.0);

          // Check that the feature eluted long enough.
          // DBoundingBox<2> box = feature.getConvexHull().getBoundingBox();    // convex hull of the entire peptide feature
          DBoundingBox<2> box = feature.getConvexHulls()[0].getBoundingBox();    // convex hull of the mono-isotopic mass trace
          if (box.maxX() - box.minX() < static_cast<double>(param_.getValue("algorithm:rt_min")))
          {
            abort = true;
            break;
          }

          features.push_back(feature);

          if (peptide == 0)
          {
            // The first/lightest peptide acts as anchor of the peptide multiplet consensus.
            // All peptide feature handles are connected to this point.
            consensus.setRT(rt);
            consensus.setMZ(mz);
            consensus.setIntensity(peptide_intensities[peptide]);
            consensus.setCharge(patterns[pattern].getCharge());
            consensus.setQuality(1.0);
          }

          FeatureHandle feature_handle;
          feature_handle.setRT(rt);
          feature_handle.setMZ(mz);
          feature_handle.setIntensity(peptide_intensities[peptide]);
          feature_handle.setCharge(patterns[pattern].getCharge());
          feature_handle.setMapIndex(peptide);
          //feature_handle.setUniqueId(&UniqueIdInterface::setUniqueId);    // TODO: Do we need to set unique ID?
          consensus.insert(feature_handle);
          consensus_map_.getColumnHeaders()[peptide].size++;
        }

        if (!abort)
        {
          consensus_map_.push_back(consensus);
          for (Feature& feat : features)
          {
            feature_map_.push_back(feat);
          }
        }

      }

    }

    endProgress();
  }

  void FeatureFinderMultiplexAlgorithm::run(MSExperiment& exp, bool progress)
  {
    // parameter section: algorithm, get selected charge range
    String charge_string = param_.getValue("algorithm:charge").toString();
    charge_min_ = charge_string.prefix(':').toInt();
    charge_max_ = charge_string.suffix(':').toInt();
    if (charge_min_ > charge_max_)
    {
      swap(charge_min_, charge_max_);
    }

    // parameter section: algorithm, get isotopes per peptide range
    String isotopes_per_peptide_string = param_.getValue("algorithm:isotopes_per_peptide").toString();
    isotopes_per_peptide_min_ = isotopes_per_peptide_string.prefix(':').toInt();
    isotopes_per_peptide_max_ = isotopes_per_peptide_string.suffix(':').toInt();
    if (isotopes_per_peptide_min_ > isotopes_per_peptide_max_)
    {
      swap(isotopes_per_peptide_min_, isotopes_per_peptide_max_);
    }

    // parameter section: labels, get all mass shifts
    label_mass_shift_["Arg6"] = param_.getValue("labels:Arg6");
    label_mass_shift_["Arg10"] = param_.getValue("labels:Arg10");
    label_mass_shift_["Lys4"] = param_.getValue("labels:Lys4");
    label_mass_shift_["Lys6"] = param_.getValue("labels:Lys6");
    label_mass_shift_["Lys8"] = param_.getValue("labels:Lys8");
    label_mass_shift_["Leu3"] = param_.getValue("labels:Leu3");
    label_mass_shift_["Dimethyl0"] = param_.getValue("labels:Dimethyl0");
    label_mass_shift_["Dimethyl4"] = param_.getValue("labels:Dimethyl4");
    label_mass_shift_["Dimethyl6"] = param_.getValue("labels:Dimethyl6");
    label_mass_shift_["Dimethyl8"] = param_.getValue("labels:Dimethyl8");
    label_mass_shift_["ICPL0"] = param_.getValue("labels:ICPL0");
    label_mass_shift_["ICPL4"] = param_.getValue("labels:ICPL4");
    label_mass_shift_["ICPL6"] = param_.getValue("labels:ICPL6");
    label_mass_shift_["ICPL10"] = param_.getValue("labels:ICPL10");

    progress_ = progress;

    // check for empty experimental data
    if (exp.getSpectra().empty())
    {
      throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS1 spectra in input file.");
    }

    //TODO allow skipping?

    // update m/z and RT ranges
    exp.updateRanges();

    // sort according to RT and MZ
    exp.sortSpectra();

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type;

    if (param_.getValue("algorithm:spectrum_type") == "automatic")
    {
      spectrum_type = exp[0].getType(true);
      // The following means that UNKNOWN will be handled as profile.
      centroided_ = (spectrum_type == SpectrumSettings::CENTROID);
    }
    else if (param_.getValue("algorithm:spectrum_type") == "centroid")
    {
      centroided_ = true;
    }
    else  // "profile"
    {
      centroided_ = false;
    }

    // store experiment in member variables
    if (centroided_)
    {
      exp.swap(exp_centroid_);
      // exp_profile_ will never be used.
    }
    else
    {
      exp.swap(exp_profile_);
      // exp_centroid_ will be constructed later on.
    }

    /**
     * pick peaks (if input data are in profile mode)
     */
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s; // peak boundaries for spectra
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c; // peak boundaries for chromatograms

    if (!centroided_)
    {
      PeakPickerHiRes picker;
      Param param = picker.getParameters();
      picker.setLogType(getLogType());
      param.setValue("ms_levels", ListUtils::create<Int>("1"));
      param.setValue("signal_to_noise", 0.0); // signal-to-noise estimation switched off
      picker.setParameters(param);

      picker.pickExperiment(exp_profile_, exp_centroid_, boundaries_exp_s, boundaries_exp_c);
    }

    /**
     * generate peak patterns for subsequent filtering step
     */
    MultiplexDeltaMassesGenerator generator = MultiplexDeltaMassesGenerator(param_.getValue("algorithm:labels").toString(), param_.getValue("algorithm:missed_cleavages"), label_mass_shift_);
    if (param_.getValue("algorithm:knock_out") == "true")
    {
      generator.generateKnockoutDeltaMasses();
    }

    #ifdef DEBUG
    generator.printSamplesLabelsList(std::cout);
    generator.printDeltaMassesList(std::cout);
    #endif

    std::vector<MultiplexDeltaMasses> masses = generator.getDeltaMassesList();
    std::vector<MultiplexIsotopicPeakPattern> patterns = generatePeakPatterns_(charge_min_, charge_max_, isotopes_per_peptide_max_, masses);

    // Switch off averagine_similarity_scaling if we search for single peptide features only.
    // (This scaling parameter is only relevant if we search for multiplets and (!) singlets.)
    double averagine_similarity_scaling;
    std::vector<std::vector<String> > list = generator.getSamplesLabelsList();
    if (list.size() == 1 && list[0].size() == 1 && list[0][0] == "no_label")
    {
      // search for singlets only
      averagine_similarity_scaling = 0.0;
    }
    else
    {
      // search for multiplets (and optionally singlets, if knock_out switched on)
      averagine_similarity_scaling = param_.getValue("algorithm:averagine_similarity_scaling");
    }

    if (centroided_)
    {
      // centroided data

      /**
       * filter for peak patterns
       */
      MultiplexFilteringCentroided filtering(exp_centroid_, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, param_.getValue("algorithm:intensity_cutoff"), param_.getValue("algorithm:rt_band"), param_.getValue("algorithm:mz_tolerance"), (param_.getValue("algorithm:mz_unit") == "ppm"), param_.getValue("algorithm:peptide_similarity"), param_.getValue("algorithm:averagine_similarity"), averagine_similarity_scaling, param_.getValue("algorithm:averagine_type").toString());
      filtering.setLogType(getLogType());
      std::vector<MultiplexFilteredMSExperiment> filter_results = filtering.filter();

      /**
       * cluster filter results
       */
      MultiplexClustering clustering(exp_centroid_, param_.getValue("algorithm:mz_tolerance"), (param_.getValue("algorithm:mz_unit") == "ppm"), param_.getValue("algorithm:rt_typical"));
      clustering.setLogType(getLogType());
      std::vector<std::map<int, GridBasedCluster> > cluster_results = clustering.cluster(filter_results);

      /**
       * construct feature and consensus maps i.e. the final results
       */
      filtering.getCentroidedExperiment().swap(exp_centroid_);
      generateMapsCentroided_(patterns, filter_results, cluster_results);
    }
    else
    {
      // profile data

      /**
       * filter for peak patterns
       */
      MultiplexFilteringProfile filtering(exp_profile_, exp_centroid_, boundaries_exp_s, patterns, isotopes_per_peptide_min_, isotopes_per_peptide_max_, param_.getValue("algorithm:intensity_cutoff"), param_.getValue("algorithm:rt_band"), param_.getValue("algorithm:mz_tolerance"), (param_.getValue("algorithm:mz_unit") == "ppm"), param_.getValue("algorithm:peptide_similarity"), param_.getValue("algorithm:averagine_similarity"), averagine_similarity_scaling, param_.getValue("algorithm:averagine_type").toString());
      filtering.setLogType(getLogType());
      std::vector<MultiplexFilteredMSExperiment> filter_results = filtering.filter();
      exp_blacklist_ = filtering.getBlacklist();

      /**
       * cluster filter results
       */
      MultiplexClustering clustering(exp_profile_, exp_centroid_, boundaries_exp_s, param_.getValue("algorithm:rt_typical"));
      clustering.setLogType(getLogType());
      std::vector<std::map<int, GridBasedCluster> > cluster_results = clustering.cluster(filter_results);

      /**
       * construct feature and consensus maps i.e. the final results
       */
      filtering.getCentroidedExperiment().swap(exp_centroid_);
      filtering.getPeakBoundaries().swap(boundaries_exp_s);
      generateMapsProfile_(patterns, filter_results, cluster_results);
    }

    // finalize consensus map

    //TODO only if sample labels are not empty
    consensus_map_.setExperimentType("labeled_MS1");
    consensus_map_.sortByPosition();
    consensus_map_.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    Size i{0};
    for (auto & ch : consensus_map_.getColumnHeaders())
    {
      ch.second.setMetaValue("channel_id", i);
      ++i;
    }

    // construct sample_labels
    std::vector<std::vector<String> > samples_labels;
    std::vector<String> temp_samples;

    String labels(param_.getValue("algorithm:labels").toString());
    boost::replace_all(labels, "[]", "no_label");
    boost::replace_all(labels, "()", "no_label");
    boost::replace_all(labels, "{}", "no_label");
    boost::split(temp_samples, labels, boost::is_any_of("[](){}")); // any bracket allowed to separate samples

    for (unsigned i = 0; i < temp_samples.size(); ++i)
    {
      if (!temp_samples[i].empty())
      {
        if (temp_samples[i]=="no_label")
        {
          vector<String> temp_labels;
          temp_labels.emplace_back("no_label");
          samples_labels.push_back(temp_labels);
        }
        else
        {
          vector<String> temp_labels;
          boost::split(temp_labels, temp_samples[i], boost::is_any_of(",;: ")); // various separators allowed to separate labels
          samples_labels.push_back(temp_labels);
        }
      }
    }

    if (samples_labels.empty())
    {
      vector<String> temp_labels;
      temp_labels.emplace_back("no_label");
      samples_labels.push_back(temp_labels);
    }

    // annotate maps
    for (unsigned i = 0; i < samples_labels.size(); ++i)
    {
      ConsensusMap::ColumnHeader& desc = consensus_map_.getColumnHeaders()[i];

      if (param_.getValue("algorithm:knock_out") == "true")
      {
        // With knock-outs present, the correct labels can only be determined during ID mapping.
        // For now, we simply store a unique identifier.
        std::stringstream stream;
        stream << "label " << i;
        desc.label = stream.str();
      }
      else
      {
        String label_string;
        for (unsigned j = 0; j < samples_labels[i].size(); ++j)
        {
          label_string.append(samples_labels[i][j]);
        }
        desc.label = label_string;
      }
    }

    // finalize feature map
    feature_map_.sortByPosition();
    feature_map_.applyMemberFunction(&UniqueIdInterface::setUniqueId);
  }

  FeatureMap& FeatureFinderMultiplexAlgorithm::getFeatureMap()
  {
    return feature_map_;
  }

  ConsensusMap& FeatureFinderMultiplexAlgorithm::getConsensusMap()
  {
    return consensus_map_;
  }

  MSExperiment& FeatureFinderMultiplexAlgorithm::getBlacklist()
  {
    return exp_blacklist_;
  }
}
