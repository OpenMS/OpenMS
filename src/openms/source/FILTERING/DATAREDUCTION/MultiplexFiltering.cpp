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
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/FILTERING/DATAREDUCTION/PeakPattern.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FilterResult.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MultiplexFiltering.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

	MultiplexFiltering::MultiplexFiltering(MSExperiment<Peak1D> exp_profile, MSExperiment<Peak1D> exp_picked, vector<vector<PeakPickerHiRes::PeakBoundary> > boundaries, std::vector<PeakPattern> patterns, int peaks_per_peptide_min, int peaks_per_peptide_max, bool missing_peaks, double intensity_cutoff, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, bool debug)
    : exp_profile_(exp_profile), exp_picked_(exp_picked), boundaries_(boundaries), patterns_(patterns), peaks_per_peptide_min_(peaks_per_peptide_min), peaks_per_peptide_max_(peaks_per_peptide_max), missing_peaks_(missing_peaks), intensity_cutoff_(intensity_cutoff), mz_tolerance_(mz_tolerance), mz_tolerance_unit_(mz_tolerance_unit), peptide_similarity_(peptide_similarity), averagine_similarity_(averagine_similarity), debug_(debug)
	{		
        if (exp_profile_.size() != exp_picked_.size())
        {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Profile and centroided data do not contain same number of spectra.");
        }
        
        if (exp_picked_.size() != boundaries_.size())
        {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Centroided data and the corresponding list of peak boundaries do not contain same number of spectra.");
        }
        
        // fill peak registry and initialise blacklist
        MSExperiment<Peak1D>::Iterator it_rt;
        vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries;
        for (it_rt = exp_picked_.begin(), it_rt_boundaries = boundaries_.begin();
            it_rt < exp_picked_.end() && it_rt_boundaries < boundaries_.end();
            ++it_rt, ++it_rt_boundaries)
        {
            int index = it_rt - exp_picked_.begin();
            
            vector<PeakReference> registry_spec;
            vector<BlackListEntry> blacklist_spec;
            for (MSSpectrum<Peak1D>::Iterator it_mz = it_rt->begin(); it_mz < it_rt->end(); ++it_mz)
            {
                // peak registry
                PeakReference reference;
                if (index > 0)
                {
                    reference.index_in_last_spectrum = getPeakIndex(index - 1, it_mz->getMZ(), 1.0);
                }
                else
                {
                    reference.index_in_last_spectrum = -1;
                }
                if (index + 1 < exp_picked_.size())
                {
                    reference.index_in_next_spectrum = getPeakIndex(index + 1, it_mz->getMZ(), 1.0);
                }
                else
                {
                    reference.index_in_next_spectrum = -1;
                }
                registry_spec.push_back(reference);

                // blacklist
                BlackListEntry entry;
                entry.black = false;
                entry.black_exception_mass_shift_index = 0;
                entry.black_exception_charge = 0;
                entry.black_exception_mz_position = 0;
                blacklist_spec.push_back(entry);
            }
            registry_.push_back(registry_spec);
            blacklist_.push_back(blacklist_spec);
        }
       
	}
    
    vector<FilterResult> MultiplexFiltering::filter()
    {
        // list of filter results for each peak pattern
        vector<FilterResult> filter_results;
        
        // loop over patterns
        for (unsigned pattern = 0; pattern < patterns_.size(); ++pattern)
        {
            cout << "peak pattern " << pattern << "\n";
            
            // data structure storing peaks which pass all filters
            FilterResult result;
            
            // m/z position passing all filters (is rejected by a particular filter)
            vector<vector<double> > debug_filtered;
            vector<vector<double> > debug_rejected;
            
            // iterate over spectra
            MSExperiment<Peak1D>::Iterator it_rt_profile;
            MSExperiment<Peak1D>::Iterator it_rt_picked;
            vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries;
            for (it_rt_profile = exp_profile_.begin(), it_rt_picked = exp_picked_.begin(), it_rt_boundaries = boundaries_.begin();
                it_rt_profile < exp_profile_.end() && it_rt_picked < exp_picked_.end() && it_rt_boundaries < boundaries_.end();
                ++it_rt_profile, ++it_rt_picked, ++it_rt_boundaries)
            {
                if ((*it_rt_picked).size() != (*it_rt_boundaries).size())
                {
                    throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"Number of peaks and number of peak boundaries differ.");
                }

                int spectrum = it_rt_profile - exp_profile_.begin();    // index of the spectrum in exp_profile_, exp_picked_ and boundaries_
                double rt_profile = it_rt_profile->getRT();
                double rt_picked = it_rt_picked->getRT();
                //cout << "    RT (profile) = " << rt_profile << "    RT (picked) = " << rt_picked << "\n";
                
                // spline fit profile data
                SplineSpectrum spline(*it_rt_profile);
                SplineSpectrum::Navigator nav = spline.getNavigator();
                
                // vectors of peak details 
                vector<double> peak_position;
                vector<double> peak_min;
                vector<double> peak_max;
                vector<double> peak_intensity;
                MSSpectrum<Peak1D>::Iterator it_mz;
                vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundary;
                for (it_mz = it_rt_picked->begin(), it_mz_boundary = it_rt_boundaries->begin();
                    it_mz < it_rt_picked->end(), it_mz_boundary < it_rt_boundaries->end();
                    ++it_mz, ++it_mz_boundary)
                {
                    peak_position.push_back(it_mz->getMZ());
                    peak_min.push_back((*it_mz_boundary).mz_min);
                    peak_max.push_back((*it_mz_boundary).mz_max);
                    peak_intensity.push_back(it_mz->getIntensity());
                    //cout << "m/z = " << it_mz->getMZ() << "  [" << (*it_mz_boundary).mz_min << ", " << (*it_mz_boundary).mz_max << "]\n";
                }
                
                // iterate over peaks in spectrum (mz)
                for (unsigned peak = 0; peak < peak_position.size(); ++peak)
                {
                    
                    /**
                     * Filter (1): m/z position and blacklist filter
                     * Are there non-black peaks with the expected relative m/z shifts?
                     */
                    std::vector<double> mz_shifts_actual;    // actual m/z shifts (differ slightly from expected m/z shifts)
                    std::vector<int> mz_shifts_actual_indices;    // peak indices in the spectrum corresponding to the actual m/z shifts
                    int peaks_found_in_all_peptides = positionsAndBlacklistFilter(patterns_[pattern], spectrum, peak_position, peak, mz_shifts_actual, mz_shifts_actual_indices);
                    if (peaks_found_in_all_peptides < peaks_per_peptide_min_)
                    {
                        if (debug_)
                        {
                            vector<double> rt_mz_flag;
                            rt_mz_flag.push_back(rt_picked);
                            rt_mz_flag.push_back(peak_position[peak]);
                            rt_mz_flag.push_back(1);    // filter 1 failed
                            debug_rejected.push_back(rt_mz_flag);
                        }
                        //continue;
                    }
                   
                    /**
                     * Filter (2): blunt intensity filter
                     * Are the mono-isotopic peak intensities of all peptides above the cutoff?
                     */
                    bool bluntVeto = monoIsotopicPeakIntensityFilter(patterns_[pattern], spectrum, mz_shifts_actual_indices);
                    if (bluntVeto)
                    {
                        if (debug_)
                        {
                            vector<double> rt_mz_flag;
                            rt_mz_flag.push_back(rt_picked);
                            rt_mz_flag.push_back(peak_position[peak]);
                            rt_mz_flag.push_back(2);    // filter 2 failed
                            debug_rejected.push_back(rt_mz_flag);
                        }
                        //continue;
                    }
                     
                    // Arrangement of peaks looks promising. Now scan through the spline fitted data.
                    vector<double> raw_entries;    // raw data points of this peak that will pass the remaining filters
                    bool blacklisted = false;    // Has this peak already been blacklisted?
                    for (double mz = peak_min[peak]; mz < peak_max[peak]; mz = nav.getNextMz(mz))
                    {
                        //cout << "m/z = " << mz << "\n";
                        
                        /**
                         * Filter (3): non-local intensity filter
                         * Are the spline interpolated intensities at m/z above the threshold?
                         */
                        vector<double> intensities_actual;    // spline interpolated intensities @ m/z + actual m/z shift
                        int peaks_found_in_all_peptides_spline = nonLocalIntensityFilter(patterns_[pattern], spectrum, mz_shifts_actual, mz_shifts_actual_indices, nav, intensities_actual, peaks_found_in_all_peptides, mz);
                        if (peaks_found_in_all_peptides_spline < peaks_per_peptide_min_)
                        {
                            if (debug_)
                            {
                                vector<double> rt_mz_flag;
                                rt_mz_flag.push_back(rt_picked);
                                rt_mz_flag.push_back(mz);
                                rt_mz_flag.push_back(3);    // filter 3 failed
                                debug_rejected.push_back(rt_mz_flag);
                            }
                            //continue;
                        }
                         
                        /**
                         * Filter (4): zeroth peak filter
                         * There should not be a significant peak to the left of the mono-isotopic
                         * (i.e. first) peak.
                         */
                        bool zero_peak_veto = zerothPeakVetoFilter(patterns_[pattern], intensities_actual);
                        if (zero_peak_veto)
                        {
                            if (debug_)
                            {
                                vector<double> rt_mz_flag;
                                rt_mz_flag.push_back(rt_picked);
                                rt_mz_flag.push_back(mz);
                                rt_mz_flag.push_back(4);    // filter 4 failed
                                debug_rejected.push_back(rt_mz_flag);
                            }
                            //continue;
                        }
                         
                        /**
                         * Filter (5): peptide similarity filter
                         * How similar are the isotope patterns of the peptides?
                         */
                        vector<double> isotope_pattern_1;
                        vector<double> isotope_pattern_2;
                        bool peptide_similarity_veto = peptideSimilarityFilter(patterns_[pattern], intensities_actual, peaks_found_in_all_peptides_spline, isotope_pattern_1, isotope_pattern_2);
                        if (peptide_similarity_veto)
                        {
                            if (debug_)
                            {
                                vector<double> rt_mz_flag;
                                rt_mz_flag.push_back(rt_picked);
                                rt_mz_flag.push_back(mz);
                                rt_mz_flag.push_back(5);    // filter 5 failed
                                debug_rejected.push_back(rt_mz_flag);
                            }
                            //continue;
                        }
                         
                        /**
                         * Filter (6): averagine similarity filter
                         * Does each individual isotope pattern resemble a peptide?
                         */
                        bool averagine_similarity_veto = averagineSimilarityFilter(patterns_[pattern], intensities_actual, peaks_found_in_all_peptides_spline, mz);
                        if (averagine_similarity_veto)
                        {
                            if (debug_)
                            {
                                vector<double> rt_mz_flag;
                                rt_mz_flag.push_back(rt_picked);
                                rt_mz_flag.push_back(mz);
                                rt_mz_flag.push_back(6);    // filter 6 failed
                                debug_rejected.push_back(rt_mz_flag);
                            }
                            //continue;
                        }
                    }
                     
                }                
             
            }
            
            // add results of this pattern to list
            filter_results.push_back(result);
        }
        
        return filter_results;
    }
    
    int MultiplexFiltering::positionsAndBlacklistFilter(PeakPattern pattern, int spectrum, std::vector<double> peak_position, int peak, std::vector<double> & mz_shifts_actual, std::vector<int> & mz_shifts_actual_indices)
    {
        // Try to find peaks at the expected m/z positions
        // loop over expected m/z shifts of a peak pattern
        for (unsigned mz_position = 0; mz_position < pattern.getMzShiftCount(); ++mz_position)
        {
            double scaling = 1;
            if (mz_position % (peaks_per_peptide_max_ + 1) == 0)
            {
                // Let us be more lenient when looking for zeroths peaks
                // i.e. allow for an increased deviation between expected m/z position and the actual one
                scaling = 2;
            }
            
            int index = getPeakIndex(peak_position, peak, peak_position[peak] + pattern.getMzShiftAt(mz_position), scaling);
            
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
                    int mz_position = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1;    // index in m/z shift list
                    if (mz_shifts_actual[mz_position] >= mzShiftNextPeptide)
                    {
                        mz_shifts_actual[mz_position] = -1000;
                        mz_shifts_actual_indices[mz_position] = -1;
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
                int mz_position = peptide * (peaks_per_peptide_max_ + 1) + isotope + 1;    // index in m/z shift list
                int peak_index = mz_shifts_actual_indices[mz_position];    // index of the peak in the spectrum
                if (peak_index != -1)
                {
                    bool black = blacklist_[spectrum][peak_index].black;
                    bool black_exception = blacklist_[spectrum][peak_index].black_exception_mass_shift_index == pattern.getMassShiftIndex() &&
                        blacklist_[spectrum][peak_index].black_exception_charge == pattern.getCharge() &&
                        blacklist_[spectrum][peak_index].black_exception_mz_position == mz_position;
                    if (black && !black_exception)
                    {
                        mz_shifts_actual[mz_position] = -1000;
                        mz_shifts_actual_indices[mz_position] = -1;
                    }
                }
            }
        }
        
        // count how many isotopic peaks seen simultaneously in all of the peptides
        // and (optionally) remove peaks following missing ones
        int peaks_found_in_all_peptides = peaks_per_peptide_max_;
        for (int peptide = 0; peptide < (int) pattern.getMassShiftCount(); ++peptide)
        {
            bool missing_peak_seen = false;
            for (int isotope = 0; isotope < peaks_per_peptide_max_; ++isotope)
            {
                int mz_position = peptide * (peaks_per_peptide_max_ +1) + isotope + 1;    // index in m/z shift list
                int index = mz_shifts_actual_indices[mz_position];    // peak index in spectrum
                if (index == -1)
                {
                    missing_peak_seen = true;
                    peaks_found_in_all_peptides = min(peaks_found_in_all_peptides, isotope);
                }
                // if missing peaks are not allowed and we have already encountered one,
                // delete all higher isotopic peaks
                if (missing_peaks_ && missing_peak_seen)
                {
                    mz_shifts_actual[mz_position] = -1000;
                    mz_shifts_actual_indices[mz_position] = -1;
                }
            }
        }
        
        return peaks_found_in_all_peptides;
    }
    
    bool MultiplexFiltering::monoIsotopicPeakIntensityFilter(PeakPattern pattern, int spectrum_index, std::vector<int> & mz_shifts_actual_indices)
    {
        MSExperiment<Peak1D>::Iterator it_rt = exp_picked_.begin() + spectrum_index;
        for (int peptide = 0; peptide < (int) pattern.getMassShiftCount(); ++peptide)
        {
            int peak_index = mz_shifts_actual_indices[peptide * (peaks_per_peptide_max_ + 1) +1];
            MSSpectrum<Peak1D>::Iterator it_mz = it_rt->begin() + peak_index;
            if (it_mz->getIntensity() < intensity_cutoff_)
            {
                return true;
            }
        }
        return false;
    }
    
    int MultiplexFiltering::nonLocalIntensityFilter(PeakPattern pattern, int spectrum_index, std::vector<double> & mz_shifts_actual, std::vector<int> & mz_shifts_actual_indices, SplineSpectrum::Navigator nav, std::vector<double> & intensities_actual, int peaks_found_in_all_peptides, double mz)
    {
        // calculate intensities
        for (int i = 0; i < (int) mz_shifts_actual_indices.size(); ++i)
        {
            if (mz_shifts_actual_indices[i] != -1)
            {
                intensities_actual.push_back(nav.eval(mz + mz_shifts_actual[i]));
            }
            else
            {
                intensities_actual.push_back(-1000.0);
            }
        }
        
        // peaks_found_in_all_peptides is the number of peaks in this region. At this particular m/z
        // some of the intensities might be below the cutoff. => return isotope
        for (int isotope = 0; isotope < peaks_found_in_all_peptides; ++isotope)
        {
            bool seen_in_all_peptides = true;
            for (int peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
            {
                if (intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1] < intensity_cutoff_)
                {
                    seen_in_all_peptides = false;
                }
            }
            if (!seen_in_all_peptides)
            {
                return isotope;
            }
        }
        
        return peaks_found_in_all_peptides;
    }
    
    bool MultiplexFiltering::zerothPeakVetoFilter(PeakPattern pattern, vector<double> & intensities_actual)
    {
        for (int peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
        {
            // scaling factor for the zeroth peak intensity
            // (The zeroth peak is problematic if its intensity exceeds zero_scaling * intensity of mono-isotopic peak.)
            double zero_scaling = 0.7;
            if (intensities_actual[peptide * (peaks_per_peptide_max_ + 1)] > zero_scaling * intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + 1])
            {
                return true;
            }
        }
        
        return false;
    }
    
    bool MultiplexFiltering::peptideSimilarityFilter(PeakPattern pattern, vector<double> & intensities_actual, int peaks_found_in_all_peptides_spline, vector<double> isotope_pattern_1, vector<double> isotope_pattern_2)
    {
        for (int peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
        {
            for (int isotope = 0; isotope < peaks_found_in_all_peptides_spline; ++isotope)
            {
                isotope_pattern_1.push_back(intensities_actual[isotope + 1]);
                isotope_pattern_2.push_back(intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1]);
            }
            if (getPatternSimilarity(isotope_pattern_1, isotope_pattern_2) < peptide_similarity_)
            {
                return true;
            }
        }
        
        return false;
    }
    
    bool MultiplexFiltering::averagineSimilarityFilter(PeakPattern pattern, vector<double> & intensities_actual, int peaks_found_in_all_peptides_spline, double mz)
    {
        for (int peptide = 0; peptide < pattern.getMassShiftCount(); ++peptide)
        {
            vector<double> isotope_pattern;
            for (int isotope = 0; isotope < peaks_found_in_all_peptides_spline; ++isotope)
            {
                isotope_pattern.push_back(intensities_actual[peptide * (peaks_per_peptide_max_ + 1) + isotope + 1]);
            }
            if (getAveragineSimilarity(isotope_pattern, mz * pattern.getCharge()) < averagine_similarity_)
            {
                return true;
            }
        }
        
        return false;
    }
    
    int MultiplexFiltering::getPeakIndex(int spectrum_index, double mz, double scaling)
    {
        MSExperiment<Peak1D>::Iterator it_rt = exp_picked_.begin() + spectrum_index;
        vector<vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries = boundaries_.begin() + spectrum_index;
        
        MSSpectrum<Peak1D>::Iterator it_mz;
        vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundaries;
        for (it_mz = it_rt->begin(), it_mz_boundaries = it_rt_boundaries->begin();
            it_mz < it_rt->end(), it_mz_boundaries < it_rt_boundaries->end();
            ++it_mz, ++it_mz_boundaries)
        {
            if (mz >= scaling * (*it_mz_boundaries).mz_min + (1 - scaling) * it_mz->getMZ() &&
                mz <= scaling * (*it_mz_boundaries).mz_max + (1 - scaling) * it_mz->getMZ())
            {
                return it_mz - it_rt->begin();
            }
            if (mz < scaling * (*it_mz_boundaries).mz_min + (1 - scaling) * it_mz->getMZ())
            {
                return -1;
            }
        }
    
        return -1;
    }
    
    int MultiplexFiltering::getPeakIndex(std::vector<double> peak_position, int start, double mz, double scaling)
    {
        vector<int> valid_index;    // indices of valid peaks that lie within the ppm range of the expected peak
        vector<double> valid_deviation;    // ppm deviations between expected and (valid) actual peaks
        
        if (peak_position[start] < mz)
        {                                                                                                   
            for (unsigned i = start; i < peak_position.size(); ++i)
            {
                double mz_min;
                double mz_max;
                if (mz_tolerance_unit_)
                {
                    mz_min = (1 - scaling * mz_tolerance_ / 1000000) * peak_position[i];
                    mz_max = (1 + scaling * mz_tolerance_ / 1000000) * peak_position[i];
                }
                else
                {
                    mz_min = peak_position[i] - scaling * mz_tolerance_;
                    mz_max = peak_position[i] + scaling * mz_tolerance_;
                }
                
                if (mz >= mz_min && mz <= mz_max)
                {
                    valid_index.push_back(i);
                    valid_deviation.push_back(std::abs(mz - peak_position[i]) / mz * 1000000);
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
                double mz_min;
                double mz_max;
                if (mz_tolerance_unit_)
                {
                    mz_min = (1 - scaling * mz_tolerance_ / 1000000) * peak_position[i];
                    mz_max = (1 + scaling * mz_tolerance_ / 1000000) * peak_position[i];
                }
                else
                {
                    mz_min = peak_position[i] - scaling * mz_tolerance_;
                    mz_max = peak_position[i] + scaling * mz_tolerance_;
                }
                
                if (mz >= mz_min && mz <= mz_max)
                {
                    valid_index.push_back(i);
                    valid_deviation.push_back(std::abs(mz - peak_position[i]) / mz * 1000000);
                }
                if (mz > peak_position[i])
                {
                    break;
                }
            }
        }
        
        if (valid_index.size() == 0)
        {
            return -1;
        }
        else
        {
            // find best index
            int best_index = -1;
            double best_deviation = valid_deviation[0];
            for (unsigned i = 1; i < valid_index.size(); ++i)
            {
                if (valid_deviation[i] < best_deviation)
                {
                    best_index = valid_index[i];
                    best_deviation = valid_deviation[i];
                }
            }
            
            return best_index;
        }
    }
    
    double MultiplexFiltering::getPatternSimilarity(vector<double> pattern1, vector<double> pattern2)
    {
        return OpenMS::Math::pearsonCorrelationCoefficient(pattern1.begin(), pattern1.end(), pattern2.begin(), pattern2.end());
    }
    
    double MultiplexFiltering::getAveragineSimilarity(vector<double> pattern, double m)
    {
        return 1.0;
    }

}
