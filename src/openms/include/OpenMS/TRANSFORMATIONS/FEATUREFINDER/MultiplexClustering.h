// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/MISC/BSpline2d.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/COMPARISON/CLUSTERING/GridBasedCluster.h>

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
    class OPENMS_DLLAPI MultiplexClustering :
        public ProgressLogger
    {        
        public:
        /**
         * @brief cluster centre, cluster bounding box, grid index
         */
        typedef GridBasedCluster::Point Point;    ///< DPosition<2>

        /**
         * @brief constructor
         * 
         * @param exp_profile    experimental data in profile mode
         * @param exp_picked    experimental data in centroid mode
         * @param boundaries    peak boundaries for exp_picked
         * @param rt_typical    elution time of a characteristic peptide in the sample
         * @param rt_minimum    shortest elution time i.e. all peptides appearing for a shorter time are being ignored
         * 
         * @throw Exception::IllegalArgument if centroided data and the corresponding list of peak boundaries do not contain same number of spectra
         */
        MultiplexClustering(const MSExperiment& exp_profile, const MSExperiment& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, double rt_typical, double rt_minimum);
        
        /**
         * @brief constructor
         * 
         * @param exp    experimental data in centroid mode
         * @param mz_tolerance    margin in m/z with which the centres of the same peak in different spectra my shift (or 'jitter')
         * @param mz_tolerance_unit    unit for mz_tolerance, ppm (true), Da (false)
         * @param rt_typical    elution time of a characteristic peptide in the sample
         * @param rt_minimum    shortest elution time i.e. all peptides appearing for a shorter time are being ignored
         * 
         * @throw Exception::IllegalArgument if centroided data and the corresponding list of peak boundaries do not contain same number of spectra
         */
        MultiplexClustering(const MSExperiment& exp, double mz_tolerance, bool mz_tolerance_unit, double rt_typical, double rt_minimum);
        
        /**
         * @brief cluster filter results
         * Data points are grouped into clusters. Each cluster contains data about one peptide multiplet.
         * 
         * @param filter_results    data points relevant for peptide multiplets i.e. output from multiplex filtering
         * 
         * @return cluster results (cluster ID, details about cluster including list of filter result IDs belonging to the cluster)
         */
        std::vector<std::map<int,GridBasedCluster> > cluster(const std::vector<MultiplexFilteredMSExperiment>& filter_results);
        
        /**
         * @brief scaled Euclidean distance for clustering
         */
        class OPENMS_DLLAPI MultiplexDistance
        {
            public:
            /**
            * @brief constructor
            * 
            * @param rt_scaling    scaling of RT coordinates before calculating Euclidean distance 
            */
            MultiplexDistance(double rt_scaling);
       
            /**
             * @brief constructor
             */
            MultiplexDistance();
            
            /**
             * @brief returns Euclidean distance
             * 
             * @param p1    first point in the (m/z,RT) plane
             * @param p2    second point in the (m/z,RT) plane
             * @return distance
             */
            double operator()(Point p1, Point p2);
     
            private:
            double rt_scaling_;
            
        };

        private:
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
        
   };
  
}

