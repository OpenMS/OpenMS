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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTERRESULTPEAK_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_MULTIPLEXFILTERRESULTPEAK_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexFilterResultRaw.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
  /**
   * @brief data structure storing a single peak that passed all filters
   * 
   * Each peak filter result corresponds to a successful search for a particular
   * peak pattern in the centroided data. The actual m/z shifts seen in the filter
   * result might differ from the theoretical shifts listed in the peak pattern.
   * 
   * @see MultiplexPeakPattern
   */
  class OPENMS_DLLAPI MultiplexFilterResultPeak
  {
    public:
    /**
     * @brief constructor
     */
    MultiplexFilterResultPeak(double mz, double rt, std::vector<double> mz_shifts,
                              std::vector<double> intensities, std::vector<MultiplexFilterResultRaw> rawDataPoints);

     /**
     * @brief returns m/z of the peak
     */
     double getMZ() const;
     
     /**
     * @brief returns RT of the peak
     */
     double getRT() const;
     
    /**
     * @brief returns m/z shifts
     */
    std::vector<double> getMZShifts() const;

    /**
     * @brief returns intensities
     */
    std::vector<double> getIntensities() const;
    
    /**
     * @brief returns the number of raw data points belonging to the peak
     */
     int size() const;
     
     /**
     * @brief returns a single raw data point belonging to the peak
     */
     MultiplexFilterResultRaw getFilterResultRaw(int i) const;
     
     private:
    /**
     * @brief position of the peak
     */
    double mz_;
    double rt_;

    /**
     * @brief m/z shifts at which peaks corresponding to a pattern were found
     */
    std::vector<double> mz_shifts_;

    /**
     * @brief peak intensities at mz_ + mz_shifts_
     */
    std::vector<double> intensities_;

    /**
     * @brief (optional) raw data points corresponding to the peak
     */
    std::vector<MultiplexFilterResultRaw> raw_data_points_;
 
  };
  
}

#endif /* MULTIPLEXFILTERRESULTPEAK_H */
