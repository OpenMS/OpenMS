// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FILTERING/DATAREDUCTION/PeakPattern.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FilterResult.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MultiplexFiltering.h>
#include <OpenMS/MATH/MISC/Spline2d.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

	MultiplexFiltering::MultiplexFiltering(MSExperiment<Peak1D> exp_profile, MSExperiment<RichPeak1D> exp_picked, std::vector<PeakPattern> patterns, int peaks_per_peptide_max, double mz_tolerance, bool mz_tolerance_unit)
    : exp_profile_(exp_profile), exp_picked_(exp_picked), patterns_(patterns), peaks_per_peptide_max_(peaks_per_peptide_max), mz_tolerance_(mz_tolerance), mz_tolerance_unit_(mz_tolerance_unit)
	{		
        if (exp_profile_.size() != exp_picked_.size())
        {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Profile and centroided data do not contain same number of spectra.");
        }
	}
    
    std::vector<FilterResult> MultiplexFiltering::filter()
    {
        // CONSTANTS
        
        //cout << "number of raw spectra = " << exp_profile_.size() << "\n";        
        //cout << "number of picked spectra = " << exp_picked_.size() << "\n";
        
                
        // list of filter results for each peak pattern
        std::vector<FilterResult> filter_results;
        
        // loop over patterns
        for (unsigned pattern = 0; pattern < patterns_.size(); ++pattern)
        {
            cout << "peak pattern " << pattern << "\n";
            
            // data structure storing peaks which pass all filters
            FilterResult * result = new FilterResult();
            
            // m/z position passing all filters (is rejected by a particular filter)
            std::vector<std::vector<double> > debug_filtered;
            std::vector<std::vector<double> > debug_rejected;
            
            // loop over (picked) spectra
            MSExperiment<Peak1D>::Iterator it_rt2 = exp_profile_.begin();
            for (MSExperiment<RichPeak1D>::Iterator it_rt = exp_picked_.begin(); it_rt != exp_picked_.end(); ++it_rt)
            {
                double rt = it_rt->getRT();
                double rt2 = it_rt2->getRT();
                //cout << "    RT = " << rt << "    RT2 = " << rt2 << "\n";
                
                // spline fit profile data
                SplineSpectrum * spline = new SplineSpectrum(*it_rt2);
                SplineSpectrum::Navigator nav = (*spline).getNavigator();
                
                // vectors of peak details 
                std::vector<double> peak_position;
                std::vector<double> peak_min;
                std::vector<double> peak_max;
                std::vector<double> peak_intensity;
                for (MSSpectrum<RichPeak1D>::Iterator it_mz = it_rt->begin(); it_mz != it_rt->end(); ++it_mz)
                {
                    peak_position.push_back(it_mz->getMZ());
                    peak_min.push_back(it_mz->getMetaValue("mzMin"));
                    peak_max.push_back(it_mz->getMetaValue("mzMax"));
                    peak_intensity.push_back(it_mz->getIntensity());
                }
                
                // iterate over data points in spectrum (mz)
                for (unsigned peak = 0; peak < peak_position.size(); ++peak)
                {
                    
                    /**
                     * Filter (1): m/z position and blacklist filter
                     * Are there non-black peaks with the expected relative m/z shifts?
                     */
                    std::vector<double> mz_shifts_actual;    // actual m/z shifts (differ slightly from expected m/z shifts)
                    std::vector<int> mz_shifts_actual_indices;
                    int peaks_found_in_all_peptides = positionsAndBlacklistFilter(patterns_[pattern], peak_position, peak, mz_shifts_actual, mz_shifts_actual_indices);
                   
                    /**
                     * Filter (2): blunt intensity filter
                     * Are the mono-isotopic peak intensities of all peptides above the cutoff?
                     */
                     
                    /**
                     * Filter (3): non-local intensity filter
                     * Are the spline interpolated intensities at m/z above the threshold?
                     */
                     
                    /**
                     * Filter (4): zeroth peak filter
                     * There should not be a significant peak to the left of the mono-isotopic
                     * (i.e. first) peak.
                     */
                     
                    /**
                     * Filter (5): peptide similarity filter
                     * How similar are the isotope patterns of the peptides?
                     */
                     
                    /**
                     * Filter (6): averagine similarity filter
                     * Does each individual isotope pattern resemble a peptide?
                     */
                     
                }                
             
                ++it_rt2;
            }
            
            // add results of this pattern to list
            filter_results.push_back(*result);
        }
        
        return filter_results;
    }
    
    int MultiplexFiltering::positionsAndBlacklistFilter(PeakPattern pattern, std::vector<double> peak_position, int peak, std::vector<double> & mz_shifts_actual, std::vector<int> & mz_shifts_actual_indices)
    {
        // Try to find peaks at the expected m/z positions
        // loop over expected m/z shifts of a peak pattern
        for (unsigned mzPosition = 0; mzPosition < pattern.getMzShiftCount(); ++mzPosition)
        {
            double scaling = 1;
            if (mzPosition % (peaks_per_peptide_max_ + 1) == 0)
            {
                // Let us be more lenient when looking for zeroths peaks
                // i.e. allow for an increased deviation between expected m/z position and the actual one
                scaling = 2;
            }
            
            int index = getPeakIndex(peak_position, peak, peak_position[peak] + pattern.getMzShiftAt(mzPosition), scaling);
            
            if (index != -1)
            {
                mz_shifts_actual.push_back(peak_position[index] - peak_position[peak]);
                mz_shifts_actual_indices.push_back(index);
            }
            else
            {
                mz_shifts_actual.push_back(-1000);
                mz_shifts_actual_indices.push_back(-1);
            }

        }
        
        // remove peaks which run into the next peptide
        // i.e. the isotopic peak of one peptide lies to the right of the mono-isotopic peak of the next one
        for (unsigned peptide = 0; peptide < pattern.getMassShiftCount() - 1; ++peptide)
        {
            double mzShiftNextPeptide = mz_shifts_actual[(peptide + 1) * (peaks_per_peptide_max_ + 1) + 1];    // m/z shift of the mono-isotopic peak of the following peptide
            if (mzShiftNextPeptide > 0)
            {
                for (int isotope = 0; isotope < peaks_per_peptide_max_; ++isotope)
                {
                    int mzPosition = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1;    // index in m/z shift list
                    if (mz_shifts_actual[mzPosition] >= mzShiftNextPeptide)
                    {
                        mz_shifts_actual[mzPosition] = -1000;
                        mz_shifts_actual_indices[mzPosition] = -1;
                    }
                }
            }
        }
        
        // remove blacklisted peaks
        // loop over isotopes in peptides
        for (int isotope = 0; isotope < peaks_per_peptide_max_; ++isotope)
        {
            // loop over peptides
            for (int peptide = 0; peptide < (int) pattern.getMassShiftCount(); ++peptide)
            {
                int mzPosition = peptide * (peaks_per_peptide_max_ + 1) + isotope +1;    // index in m/z shift list
                int peakIndex = mz_shifts_actual_indices[mzPosition];    // index of the peak in the spectrum
                if (peakIndex != -1)
                {
                    bool black = true;
                    bool blackException = true;
                }
            }
        }
        
        // count how many isotopic peaks seen simultaneously in all of the peptides
        // and (optionally) remove peaks following missing ones
        
        return 3;
    }
    
    int MultiplexFiltering::getPeakIndex(std::vector<double> peak_position, int start, double mz, double scaling)
    {
        std::vector<int> validIndex;    // indices of valid peaks that lie within the ppm range of the expected peak
        std::vector<double> validDeviation;    // ppm deviations between expected and (valid) actual peaks
        
        if (peak_position[start] < mz)
        {                                                                                                   
            for (unsigned i = start; i < peak_position.size(); ++i)
            {
                double mzMin;
                double mzMax;
                if (mz_tolerance_unit_)
                {
                    mzMin = (1 - scaling * mz_tolerance_ / 1000000) * peak_position[i];
                    mzMax = (1 + scaling * mz_tolerance_ / 1000000) * peak_position[i];
                }
                else
                {
                    mzMin = peak_position[i] - scaling * mz_tolerance_;
                    mzMax = peak_position[i] + scaling * mz_tolerance_;
                }
                
                if (mz >= mzMin && mz <= mzMax)
                {
                    validIndex.push_back(i);
                    validDeviation.push_back(std::abs(mz - peak_position[i]) / mz * 1000000);
                }
                if (mz < peak_position[i])
                {
                    break;
                }
            }
        }
        else
        {
            for (int i = start; i >= 0; --i)
            {
                double mzMin;
                double mzMax;
                if (mz_tolerance_unit_)
                {
                    mzMin = (1 - scaling * mz_tolerance_ / 1000000) * peak_position[i];
                    mzMax = (1 + scaling * mz_tolerance_ / 1000000) * peak_position[i];
                }
                else
                {
                    mzMin = peak_position[i] - scaling * mz_tolerance_;
                    mzMax = peak_position[i] + scaling * mz_tolerance_;
                }
                
                if (mz >= mzMin && mz <= mzMax)
                {
                    validIndex.push_back(i);
                    validDeviation.push_back(std::abs(mz - peak_position[i]) / mz * 1000000);
                }
                if (mz > peak_position[i])
                {
                    break;
                }
            }
        }
        
        if (validIndex.size() == 0)
        {
            return -1;
        }
        else
        {
            // find best index
            int bestIndex = -1;
            double bestDeviation = validDeviation[0];
            for (unsigned i = 1; i < validIndex.size(); ++i)
            {
                if (validDeviation[i] < bestDeviation)
                {
                    bestIndex = validIndex[i];
                    bestDeviation = validDeviation[i];
                }
            }
            
            return bestIndex;
        }
    }
        
}
