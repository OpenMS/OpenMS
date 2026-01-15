// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/MISC/BSpline2d.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/FEATUREFINDER/MultiplexFiltering.h>
#include <OpenMS/ML/CLUSTERING/GridBasedCluster.h>

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
         * 
         * @throw Exception::IllegalArgument if centroided data and the corresponding list of peak boundaries do not contain same number of spectra
         */
        MultiplexClustering(const MSExperiment& exp_profile, const MSExperiment& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, double rt_typical);
        
        /**
         * @brief constructor
         * 
         * @param exp    experimental data in centroid mode
         * @param mz_tolerance    margin in m/z with which the centres of the same peak in different spectra my shift (or 'jitter')
         * @param mz_tolerance_unit    unit for mz_tolerance, ppm (true), Da (false)
         * @param rt_typical    elution time of a characteristic peptide in the sample
         * 
         * @throw Exception::IllegalArgument if centroided data and the corresponding list of peak boundaries do not contain same number of spectra
         */
        MultiplexClustering(const MSExperiment& exp, double mz_tolerance, bool mz_tolerance_unit, double rt_typical);
        
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
            double operator()(const Point& p1, const Point& p2) const;
     
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
        //unused
        //double rt_minimum_;
        
   };
  
}

