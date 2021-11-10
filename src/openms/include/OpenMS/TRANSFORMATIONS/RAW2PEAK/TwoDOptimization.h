// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

//#define DEBUG_2D
#undef DEBUG_2D

#ifdef DEBUG_2D
#include <iostream>
#include <fstream>
#endif

#include <vector>
#include <utility>
#include <cmath>
#include <set>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/PeakIndex.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>


#ifndef OPENMS_SYSTEM_STOPWATCH_H
#endif

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>

namespace OpenMS
{

  /**
    @brief This class provides the two-dimensional optimization of the picked peak parameters.

    Given the picked peaks, this class optimizes the peak parameters of each isotope pattern using
    a non-linear optimization. The peaks of adjacent scans are adjusted to achieve that a peak occurring in
    several scans has always the same m/z position. For the optimization the Levenberg-Marquardt algorithm
    provided from the Eigen is used. The optimized parameters are the m/z values,
    the left and right width, which shall be equal for a peak in all scans,
    and the peaks' heights.

    @todo Works only with defined types due to pointers to the data in the optimization namespace! Change that or remove templates (Alexandra)

    @htmlinclude OpenMS_TwoDOptimization.parameters
  */
  class OPENMS_DLLAPI TwoDOptimization :
    public DefaultParamHandler
  {
public:

    typedef MSExperiment::const_iterator InputSpectrumIterator;

    /// Constructor
    TwoDOptimization();

    /// Copy constructor
    TwoDOptimization(const TwoDOptimization& opt);

    /// Destructor
    ~TwoDOptimization() override{}

    /// Assignment operator
    TwoDOptimization& operator=(const TwoDOptimization& opt);


    ///Non-mutable access to the matching epsilon
    inline double getMZTolerance() const {return tolerance_mz_; }
    ///Mutable access to the matching epsilon
    inline void setMZTolerance(double tolerance_mz)
    {
      tolerance_mz_ = tolerance_mz;
      param_.setValue("2d:tolerance_mz", tolerance_mz);
    }

    ///Non-mutable access to the maximal peak distance in a cluster
    inline double getMaxPeakDistance() const {return max_peak_distance_; }
    ///Mutable access to the maximal peak distance in a cluster
    inline void setMaxPeakDistance(double max_peak_distance)
    {
      max_peak_distance_ = max_peak_distance;
      param_.setValue("2d:max_peak_distance", max_peak_distance);
    }

    ///Non-mutable access to the maximal number of iterations
    inline UInt getMaxIterations() const {return max_iteration_; }
    ///Mutable access to the  maximal number of iterations
    inline void setMaxIterations(UInt max_iteration)
    {
      max_iteration_ = max_iteration;
      param_.setValue("iterations", max_iteration);
    }

    ///Non-mutable access to the minimal number of adjacent scans
    inline const OptimizationFunctions::PenaltyFactorsIntensity& getPenalties() const {return penalties_; }
    ///Mutable access to the minimal number of adjacent scans
    inline void setPenalties(OptimizationFunctions::PenaltyFactorsIntensity& penalties)
    {
      penalties_ = penalties;
      param_.setValue("penalties:position", penalties.pos);
      param_.setValue("penalties:height", penalties.height);
      param_.setValue("penalties:left_width", penalties.lWidth);
      param_.setValue("penalties:right_width", penalties.rWidth);
    }

    /**
        @brief Find two dimensional peak clusters and optimize their peak parameters

        @note For the peak spectra, the following meta data arrays (see MSSpectrum) have to be present and have to be named just as listed here:
        - intensity (index:1)
        - leftWidth (index:3)
        - rightWidth (index:4)
        - peakShape (index:5)

        @param first begin of the raw data spectra iterator range
        @param last end of the raw data spectra iterator range
        @param ms_exp peak map corresponding to the raw data in the range from @p first to @p last
        @param real2D flag if the optimization should be two dimensional or on each scan separately
        @exception Exception::IllegalArgument is thrown if required meta information from peak picking is missing (area, shape, left width, right width) or if the input data is invalid in some other way

    */
    void optimize(InputSpectrumIterator first,
                  InputSpectrumIterator last,
                  PeakMap& ms_exp, bool real2D = true);


protected:
    /// Helper struct (contains the size of an area and a raw data container)
    struct Data
    {
      std::vector<std::pair<SignedSize, SignedSize> > signal2D;
      std::multimap<double, IsotopeCluster>::iterator iso_map_iter;
      Size total_nr_peaks;
      std::map<Int, std::vector<PeakIndex> > matching_peaks;
      PeakMap picked_peaks;
      PeakMap::ConstIterator raw_data_first;
      OptimizationFunctions::PenaltyFactorsIntensity penalties;
      std::vector<double> positions;
      std::vector<double> signal;
    };

    class OPENMS_DLLAPI TwoDOptFunctor
    {
    public:
      int inputs() const { return m_inputs; }
      int values() const { return m_values; }

      TwoDOptFunctor(unsigned dimensions, unsigned num_data_points, const TwoDOptimization::Data* data)
      : m_inputs(dimensions), m_values(num_data_points), m_data(data) {}

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec);
      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &J);

    private:
      const int m_inputs, m_values;
      const TwoDOptimization::Data* m_data;
    };


    /// stores the retention time of each isotopic cluster
    std::multimap<double, IsotopeCluster> iso_map_;

    /// Pointer to the current region
    std::multimap<double, IsotopeCluster>::const_iterator curr_region_;

    /// upper bound for distance between two peaks belonging to the same region
    double max_peak_distance_;

    /// threshold for the difference in the peak position of two matching peaks
    double tolerance_mz_;

    /// Indices of peaks in the adjacent scans matching peaks in the scan with no. ref_scan
    //  std::map<Int, std::vector<PeakMap::SpectrumType::Iterator > > matching_peaks_;
    std::map<Int, std::vector<PeakIndex> > matching_peaks_;


    /// Convergence Parameter: Maximal number of iterations
    UInt max_iteration_;

    /// Optimization considering all scans of a cluster or optimization of each scan separately
    bool real_2D_;


    /// Penalty factors for some parameters in the optimization
    OptimizationFunctions::PenaltyFactorsIntensity penalties_;


    /**
         @name Auxiliary Functions for the search of matching regions
    */
    //@{
    std::vector<double>::iterator searchInScan_(std::vector<double>::iterator scan_begin,
                                                    std::vector<double>::iterator scan_end,
                                                    double current_mz);

    /** Performs 2D optimization of all regions */
    void optimizeRegions_(InputSpectrumIterator& first,
                          InputSpectrumIterator& last,
                          PeakMap& ms_exp);

    /** Performs an optimization of all regions by calling OptimizePick */
    void optimizeRegionsScanwise_(InputSpectrumIterator& first,
                                  InputSpectrumIterator& last,
                                  PeakMap& ms_exp);


    /// Get the indices of the first and last raw data point of this region
    void getRegionEndpoints_(PeakMap& exp,
                             InputSpectrumIterator& first,
                             InputSpectrumIterator& last,
                             Size iso_map_idx,
                             double noise_level,
                             TwoDOptimization::Data& d);

    /// Identify matching peak in a peak cluster
    void findMatchingPeaks_(std::multimap<double, IsotopeCluster>::iterator& it,
                            PeakMap& ms_exp);

    //@}

    /// update members method from DefaultParamHandler to update the members
    void updateMembers_() override;
  };

}

