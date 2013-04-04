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
// $Authors: Steffen Sass, Holger Plattfaut, Bastian Blank $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_SILACFILTERING_H
#define OPENMS_FILTERING_DATAREDUCTION_SILACFILTERING_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <list>
#include <map>
#include <vector>

namespace OpenMS
{
  class SILACFilter;

  /**
    @brief Filtering for SILAC data.

    This filtering can be used to extract SILAC features from an MS experiment.
    Several SILACFilters can be added to the filtering to search for specific SILAC patterns.

    @see SILACFilter
  */
  class OPENMS_DLLAPI SILACFiltering :
    public ProgressLogger
  {
public:
    typedef std::vector<SILACFilter> Filters;

    /**
     * @brief holds all filters used in the filtering
     */
    Filters filters_;

    /**
     * @brief Wrapper class for spectrum interpolation
     */
    class OPENMS_DLLAPI SpectrumInterpolation
    {
private:
      gsl_interp_accel * current_;
      gsl_spline * spline_;

public:
      SpectrumInterpolation(const MSSpectrum<> &, const SILACFiltering &);
      ~SpectrumInterpolation();

      DoubleReal operator()(DoubleReal mz) const
      {
        return gsl_spline_eval(spline_, mz, current_);
      }

    };

private:
    /**
     * @brief minimal intensity of SILAC features
     */
    DoubleReal intensity_cutoff_;

    /**
     * @brief raw data
     */
    MSExperiment<Peak1D> & exp_;

    /**
     * @brief picked data
     */
    MSExperiment<Peak1D> picked_exp_;

    /**
     * @brief picked data seeds
     */
    MSExperiment<Peak1D> picked_exp_seeds_;

    /**
     * Filename base for debugging output
     */
    const String debug_filebase_;

    /**
     * @brief pick data seeds
     */
    void pickSeeds_();

    /**
     * @brief apply filtering to picked data seeds
     */
    void filterSeeds_();

public:

    /// peak-width equation
    const PeakWidthEstimator::Result peak_width;

    /**
     * @brief detailed constructor
     * @param exp raw data
     * @param intensity_cutoff minimal intensity of SILAC features
     */
    SILACFiltering(MSExperiment<Peak1D> & exp, const PeakWidthEstimator::Result &, const DoubleReal intensity_cutoff, const String debug_filebase_ = "");

    /**
     * @brief adds a new filter to the filtering
     * @param filter filter to add
     */
    void addFilter(SILACFilter & filter);

    /**
     * @brief starts the filtering based on the added filters
     */
    void filterDataPoints();

    /**
     * @brief structure for blacklist
     * @param range m/z and RT interval to be blacklisted
     * @param charge charge of the generating filter
     * @param mass_separations mass separations of the generating filter
     * @param relative_peak_position m/z position of the blacklisted area relative to the mono-isotopic peak of the unlabelled peptide
     */
    struct BlacklistEntry
    {
      DRange<2> range;
      Int charge;
      std::vector<DoubleReal> mass_separations;
      DoubleReal relative_peak_position;
    };

    /**
     * @brief holds the range that is blacklisted for other filters and the filter that generated the blacklist entry
     */
    std::multimap<DoubleReal, BlacklistEntry> blacklist;
  };
}

#endif /* SILACFILTERING_H_ */
