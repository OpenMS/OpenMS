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

#ifndef OPENMS_FILTERING_DATAREDUCTION_MULTIPLEXFILTERING_H
#define OPENMS_FILTERING_DATAREDUCTION_MULTIPLEXFILTERING_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/FILTERING/DATAREDUCTION/PeakPattern.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FilterResult.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
    /**
     * @brief filters peak and profile data for peak patterns
     * 
     * The algorithm searches for patterns of multiple peptides in the data.
     * The peptides appear as characteristic patterns of isotopic peaks in
     * MS1 spectra. We first search the centroided data, and optionally in
     * a second step the spline interpolated profile data. For each
     * peak pattern the algorithm generates a filter result.
     * 
     * @see PeakPattern
     * @see FilterResult
     */
    class OPENMS_DLLAPI MultiplexFiltering
    {        
        public:
        /// structure for peak position in neighbouring spectra
        struct PeakReference
        {
            int index_in_last_spectrum;
            int index_in_next_spectrum;
        };

        /// structure for peak blacklisting
        struct BlackListEntry
        {
            bool black;
            int black_exception_mass_shift_index;
            int black_exception_charge;
            int black_exception_mz_position;
        };

        /**
         * @brief constructor
         * 
         * @param exp_profile    experimental data in profile mode
         * @param exp_picked    experimental data in centroid mode
         * @param boundaries    peak boundaries for exp_picked
         * @param patterns    patterns of isotopic peaks to be searched for
         * @param peaks_per_peptide_min    minimum number of isotopic peaks in peptides
         * @param peaks_per_peptide_max    maximum number of isotopic peaks in peptides
         * @param missing_peaks    flag for missing peaks
         * @param mz_tolerance    error margin in m/z for matching expected patterns to experimental data
         * @param mz_tolerance_unit    unit for mz_tolerance, ppm (true), Da (false)
         * @param debug    debug mode
         */
        MultiplexFiltering(MSExperiment<Peak1D> exp_profile, MSExperiment<Peak1D> exp_picked, std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries, std::vector<PeakPattern> patterns, int peaks_per_peptide_min, int peaks_per_peptide_max, bool missing_peaks, double mz_tolerance, bool mz_tolerance_unit, bool debug);
        
        /**
         * @brief filter for patterns
         * (generates a filter result for each pattern)
         */
        std::vector<FilterResult> filter();
        
        private:
        /**
         * @brief position and blacklist filter
         * 
         * @param pattern
         * @param peak_position
         * @param peak
         * @param mz_shifts_actual
         * @param mz_shifts_actual_indices
         */
        int positionsAndBlacklistFilter(PeakPattern pattern, int spectrum, std::vector<double> peak_position, int peak, std::vector<double> & mz_shifts_actual, std::vector<int> & mz_shifts_actual_indices);
        
        /**
         * @brief mono-isotopic peak intensity filter
         * 
         * @param ...
         */
        bool monoIsotopicPeakIntensityFilter(PeakPattern pattern);
        
        /**
         * @brief returns the index of a peak at m/z
         * (for initialisation of peak registry)
         * 
         * @param spectrum_index    index of the spectrum in exp_picked_ and boundaries_
         * @param mz
         * @param scaling
         */
        int getPeakIndex(int spectrum_index, double mz, double scaling);
        
        /**
         * @brief returns the index of a peak at m/z
         * (finds not only a valid peak, i.e. within certain m/z deviation, but the best of the valid peaks)
         * 
         * @param peak_position
         * @param start
         * @param mz
         * @param scaling
         */
        int getPeakIndex(std::vector<double> peak_position, int start, double mz, double scaling);
        
        /**
         * @brief profile and centroided experimental data
         */
        MSExperiment<Peak1D> exp_profile_;
        MSExperiment<Peak1D> exp_picked_;
        std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_;
        
         /**
         * @brief auxiliary structs for navigation and blacklisting 
         */
        std::vector<std::vector<PeakReference> > registry_;
        std::vector<std::vector<BlackListEntry> > blacklist_;
       
        /**
         * @brief list of peak patterns
         */
        std::vector<PeakPattern> patterns_;

        /**
         * @brief minimum number of isotopic peaks per peptide
         */
        int peaks_per_peptide_min_;

        /**
         * @brief maximum number of isotopic peaks per peptide
         */
        int peaks_per_peptide_max_;

        /**
         * @brief flag for missing peaks
         */
        bool missing_peaks_;

        /**
         * @brief m/z shift tolerance
         */
        double mz_tolerance_;

        /**
         * @brief unit for m/z shift tolerance (ppm - true, Da - false)
         */
        bool mz_tolerance_unit_;
        
        /**
         * @brief debug mode
         */
        bool debug_;
        
   };
  
}

#endif /* MULTIPLEXFILTERING_H_ */
