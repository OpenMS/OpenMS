// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredPeak.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
    /**
     * @brief data structure storing all peaks (and optionally their raw data points)
     * of an experiment corresponding to one specific peak pattern
     * 
     * @see MultiplexPeakPattern
     */
    class OPENMS_DLLAPI MultiplexFilteredMSExperiment
    {
        public:
        /**
         * @brief constructor
         */
        MultiplexFilteredMSExperiment();
        
        /**
         * @brief adds a single peak to the results
         */
        void addPeak(const MultiplexFilteredPeak& peak);
        
        /**
         * @brief returns a single peak from the results
         */
        MultiplexFilteredPeak getPeak(size_t i) const;
               
        /**
         * @brief returns m/z of a single peak
         */
        double getMZ(size_t i) const;
        
        /**
         * @brief returns m/z positions of all peaks
         */
        std::vector<double> getMZ() const;
        
        /**
         * @brief returns RT of a single peak
         */
        double getRT(size_t i) const;
        
        /**
         * @brief returns RT of all peaks
         */
        std::vector<double> getRT() const;
        
        /**
         * @brief returns number of peaks in the result
         */
        size_t size() const;
        
        /**
         * @brief write all peaks to a consensusXML file
         * 
         * @param exp_picked   original (i.e. not white) centroided experimental data
         * @param debug_out    file name of the debug output
         */
        void writeDebugOutput(const MSExperiment& exp_picked, const String& debug_out) const;
        
        private:
        /**
         * @brief peaks which passed the peak pattern filter
         */
        std::vector<MultiplexFilteredPeak> result_;

   };
  
}

