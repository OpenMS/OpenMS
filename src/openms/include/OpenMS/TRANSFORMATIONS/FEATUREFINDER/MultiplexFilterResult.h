// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTERRESULT_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTERRESULT_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultRaw.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
    /**
     * @brief data structure storing all peaks (and optionally their raw data points)
     * corresponding to one specific peak pattern
     * 
     * @see MultiplexPeakPattern
     */
    class OPENMS_DLLAPI MultiplexFilterResult
    {
        public:
        /**
         * @brief constructor
         */
        MultiplexFilterResult();
        
        /**
         * @brief adds a single peak to the results
         */
        void addFilterResultPeak(double mz, double rt, std::vector<double> mzShifts, std::vector<double> intensities, std::vector<MultiplexFilterResultRaw> result);
        
        /**
         * @brief returns a single peak from the results
         */
        MultiplexFilterResultPeak getFilterResultPeak(int i) const;
        
        /**
         * @brief returns a single raw data point from peak i in the result
         */
        MultiplexFilterResultRaw getFilterResultRaw(int i, int j) const;
        
        /**
         * @brief returns m/z of a single peak
         */
        double getMZ(int i) const;
        
        /**
         * @brief returns m/z positions of all peaks
         */
        std::vector<double> getMZ() const;
        
        /**
         * @brief returns RT of a single peak
         */
        double getRT(int i) const;
        
        /**
         * @brief returns RT of all peaks
         */
        std::vector<double> getRT() const;
        
        /**
         * @brief returns number of peaks in the result
         */
        int size() const;
        
        private:
        /**
         * @brief peaks which passed the peak pattern filter
         */
        std::vector<MultiplexFilterResultPeak> result_;

   };
  
}

#endif /* MULTIPLEXFILTERRESULT_H_ */
