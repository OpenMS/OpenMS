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

#ifndef OPENMS_FILTERING_DATAREDUCTION_MULTIPLEXCLUSTERING_H
#define OPENMS_FILTERING_DATAREDUCTION_MULTIPLEXCLUSTERING_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
//#include <OpenMS/FILTERING/DATAREDUCTION/PeakPattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResult.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/COMPARISON/CLUSTERING/MultiplexCluster.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
    /**
     * @brief clusters results from multiplex filtering
     * 
     * The multiplex filtering algorithm identified regions in the picked and
     * profile data that correspond to peptide features. This clustering algorithm
     * takes these filter results as input and groups data points that belong to
     * the same peptide features. It makes use of the general purpose hierarchical
     * clustering implementation LocalClustering. 
     * 
     * @see MultiplexFiltering
     * @see LocalClustering
     */
    class OPENMS_DLLAPI MultiplexClustering
    {        
        public:
        /**
         * @brief constructor
         * 
         * @param exp_profile    experimental data in profile mode
         * @param exp_picked    experimental data in centroid mode
         * @param boundaries    peak boundaries for exp_picked
         * @param rt_typical    elution time of a characteristic peptide in the sample
         * @param rt_minimum    shortest elution time i.e. all peptides appearing for a shorter time are being ignored
         * @param debug    debug mode
         * 
         * @throw Exception::IllegalArgument if centroided data and the corresponding list of peak boundaries do not contain same number of spectra
         */
        MultiplexClustering(MSExperiment<Peak1D> exp_profile, MSExperiment<Peak1D> exp_picked, std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries, double rt_typical, double rt_minimum, bool debug);
        
        /**
         * @brief cluster filter results
         * Data points are grouped into clusteres. Each cluster contains data about one peptide multiplet.
         * 
         * @param filter_results    data points relevant for peptide multiplets i.e. output from multiplex filtering
         * 
         * @return cluster results (cluster ID, details about cluster including list of filter result IDs belonging to the cluster)
         */
        std::vector<std::map<int,MultiplexCluster> > cluster(std::vector<MultiplexFilterResult> filter_results);
        
        /**
         * @brief rough estimation of the peak width at m/z
         * 
         * Based on the peaks of the dataset (peak position & width), the typical peak width is estimated for arbitrary m/z.
         * The peak width is assumed to be a linear funtion of m/z.
         * 
         * TO DO: Use Lowess instead of Linear Regression.
         */
        class OPENMS_DLLAPI PeakWidthEstimator
        {
            public:
            /**
            * @brief constructor
            * 
            * @param exp_picked    experimental data in centroid mode
            * @param boundaries    peak boundaries for exp_picked
            * 
            * @throw Exception::IllegalArgument if centroided data and the corresponding list of peak boundaries do not contain same number of spectra
            */
            PeakWidthEstimator(MSExperiment<Peak1D> exp_picked, std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries);
        
            /**
            * @brief returns the estimated peak width at m/z
            * 
            * @throw Exception::InvalidValue if the peak width estimation returns a negative value.
            */
            double getPeakWidth(double mz);
        
            private:        
            /// hide default constructor
            PeakWidthEstimator();
        
            /**
            * @brief m/z range of peak width interpolation
            */
            double mz_min_;
            double mz_max_;
            
            /**
             * @brief slope and intercept for peak width interpolation
             */
            double slope_;
            double intercept_;
        
        };

        private:
        /**
         * @brief structure for debug output
         * 
         * Position of a filter result data point
         * and ID of the cluster it belongs to.
         */
        struct DebugPoint
        {
            double rt;
            double mz;
            double cluster;
        };
        
        /**
         * @brief write debug output
         * 
         * Debug data written to file.
         * 
         * @param points    data points for debug output
         * @param pattern    pattern ID
         */
        void writeDebug(std::vector<DebugPoint> points, int pattern) const;

        /**
         * @brief returns a colour for ID
         * 
         * @param c    integer ID
         * @return string for one of 50 colours
         */
        String getColour(int c) const;

        /**
         * @brief grid spacing for clustering
         */
        std::vector<double> grid_spacing_mz_;
        std::vector<double> grid_spacing_rt_;
        
        /**
         * @brief scaling in y-direction for clustering
         */
        double rt_scaling_;
        
        /**
         * @brief typical retention time
         */
        double rt_typical_;
        
        /**
         * @brief minimum retention time
         */
        double rt_minimum_;
        
        /**
         * @brief debug mode
         */
        bool debug_;
   };
  
}

#endif /* MULTIPLEXCLUSTERING_H_ */
