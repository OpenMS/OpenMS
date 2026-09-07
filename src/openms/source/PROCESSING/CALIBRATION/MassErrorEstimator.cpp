// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Athenea Valls Bleher $
// $Authors: Athenea Valls Bleher $
// --------------------------------------------------------------------------


#include <OpenMS/PROCESSING/CALIBRATION/MassErrorEstimator.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <set>

namespace OpenMS 
{

std::vector<MarkerPeak> MassErrorEstimator::findPolysiloxaneCandidates (const MSExperiment& experiment)
{
    std::vector<MarkerPeak> candidates;
    double mz_window = Math::ppmToMass(ppm_tolerance_, reference_mz_); 
    
    for (const auto& spec : experiment.getSpectra())
    {
        if (spec.getMSLevel() != 1) continue;

        // estimate local noise level from nearby peaks
        std::vector<double> intensities; 
        for(const auto& peak: spec)
        {
            if(peak.getMZ() < 444.0 || peak.getMZ() > 446.0) continue;
            intensities.push_back(peak.getIntensity());
        }
        if (intensities.empty()) continue;

        std::sort(intensities.begin(),intensities.end());

        // use median intensity for robust noise estimation
        double median = 0.0;
        Size n = intensities.size();
        if (n % 2 == 0)
        {
            median = (intensities[n / 2 - 1] + intensities[n / 2]) / 2.0;
        }
        else
        {
            median = intensities[n / 2];
        }

        for (const auto& peak : spec)
        {
            if (peak.getMZ() < reference_mz_ - mz_window || peak.getMZ()> reference_mz_ + mz_window) continue;
            if (peak.getIntensity() < median * intensity_factor_)continue;
            candidates.push_back({peak.getMZ(), spec.getRT(), peak.getIntensity()});
        }
    }
    return candidates;
}

std::vector<MarkerPeak> MassErrorEstimator::RTdistribution (const std::vector<MarkerPeak>& candidates)
{
    if (candidates.empty()) return {};
    std::vector<MarkerPeak> sorted = candidates; 

    std::sort(sorted.begin(), sorted.end(), 
            [](const MarkerPeak& a, const MarkerPeak& b)
            {
                return a.rt < b.rt;
            });

    // select the most intense peak per RT bin to avoid oversampling
    std::vector<MarkerPeak> selected;
    double min_rt = sorted.front().rt;
    double max_rt = sorted.back().rt;
    double bin_width = (max_rt - min_rt) / rt_bins_; 
    for (Size i = 0; i < Size(rt_bins_); ++i)
    {
        double bin_start = min_rt + i * bin_width;
        double bin_end   = min_rt + (i + 1) * bin_width;
        bool found = false;
        MarkerPeak best;
        double best_intensity = -1.0;
        for (const auto& p : sorted)
        {
            if (p.rt < bin_start || p.rt >= bin_end) continue;
            if (!found || p.intensity > best_intensity)
            {
                best = p;
                best_intensity = p.intensity;
                found = true;
            }
        }
        if (found)
        {
            selected.push_back(best);
        }
    }
    return selected;
}

double MassErrorEstimator::calculateMean (const std::vector<MarkerPeak>& selected, double ref)
{
    if (selected.empty())
    {
        return 0.0;
    }

    std::vector<double> ppm_errors;
    double sum= 0.0;
    for (const auto& c : selected)
    {
        double ppm_error = Math::getPPM(c.mz, ref);
        sum += ppm_error;
        ppm_errors.push_back(ppm_error);
    }
    double mean = sum / ppm_errors.size();
    return mean; 
}

double MassErrorEstimator::calculateSD (const std::vector<MarkerPeak>& selected, double ref)
{
    if (selected.size() < 2)
    {
        return 0.0;
    }
    std::vector<double> ppm_errors;
    double sum = 0.0;
    double mean = calculateMean(selected, ref);
    for (const auto& c: selected)
    {
        double ppm_error=Math::getPPM(c.mz, ref);
        ppm_errors.push_back(ppm_error);
    }
    for (double x : ppm_errors)
    {
        sum += (x - mean) * (x - mean);
    }
    double sd = std::sqrt(sum / (ppm_errors.size()-1));
    return sd;

}

void MassErrorEstimator::estimate (const MSExperiment& experiment)
{
    std::vector<MarkerPeak> polysiloxane_candidates = MassErrorEstimator::findPolysiloxaneCandidates(experiment);

    if (polysiloxane_candidates.empty())
    {
        std::cout << "No polysiloxane contaminant found --> skipping\n";
    }
    else
    {
        std::vector<MarkerPeak> dist= MassErrorEstimator::RTdistribution(polysiloxane_candidates);

        std::set<double> unique_mz;
        for (const auto& c : dist)
        {
            unique_mz.insert(c.mz);
        }
        // insufficient unique m/z values may indicate centroid quantization
        // or lockmass correction to the contaminant peak
        if (unique_mz.size() < min_unique_mz_)
        {
            std::cerr << "Error: too few unique m/z values detected --> " 
                    << "mass error estimation based on polysiloxane contaminant unreliable\n";
            return;
        }

        double mean = calculateMean(dist,reference_mz_); 

        double sd = calculateSD (dist,reference_mz_);

        std::cout<< "mean (Error Bias) = " << mean << " ppm" << "\n"
                 << "standard deviation (Error Scatter) = " << sd << " ppm" << "\n";

        double window_beginning = mean - 3 * sd;
        double window_end = mean + 3 * sd;
        std::cout<< "The recommended mass tolerance window is: ["<<window_beginning<<" , "<<window_end<<"] and was calculated with factor 3 "<<  "\n";

    }
}
}