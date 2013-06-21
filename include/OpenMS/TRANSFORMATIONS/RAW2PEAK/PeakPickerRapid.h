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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERRAPID_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERRAPID_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>

#include <boost/dynamic_bitset.hpp>

#include <map>
#include <algorithm>
#include <limits>

#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#undef DEBUG_DECONV
namespace OpenMS
{
/**
   @brief This class implements a fast peak-picking algorithm best suited for high resolution MS data (FT-ICR-MS, Orbitrap). In high resolution data, the signals of ions with similar mass-to-charge ratios (m/z) exhibit little or no overlapping and therefore allow for a clear separation. Furthermore, ion signals tend to show well-defined peak shapes with narrow peak width.

   This peak-picking algorithm detects ion signals in raw data and reconstructs the corresponding peak shape by cubic spline interpolation. Signal detection depends on the signal-to-noise ratio which is adjustable by the user (see parameter signal_to_noise). A picked peak's m/z and intensity value is given by the maximum of the underlying peak spline.

   So far, this peak picker was mainly tested on high resolution data. With appropriate preprocessing steps (e.g. noise reduction and baseline subtraction), it might be also applied to low resolution data.

   @htmlinclude OpenMS_PeakPickerRapid.parameters

   @note The peaks must be sorted according to ascending m/z!

   @ingroup PeakPicking
  */


class OPENMS_DLLAPI CmpPeakByIntensity
{
public:
    template <typename PeakType>
    bool operator()(PeakType x, PeakType y) const
    {
        return x.getIntensity() > y.getIntensity();
    }
};

class OPENMS_DLLAPI PeakPickerRapid
        : public DefaultParamHandler,
        public ProgressLogger
{
public:
    /// Constructor
    PeakPickerRapid();

    /// Destructor
    virtual ~PeakPickerRapid();

    template <typename PeakType>
    bool computeTPG(const PeakType& p1, const PeakType& p2, const PeakType& p3, DoubleReal& mu, DoubleReal& sigma, DoubleReal& area, DoubleReal& height) const
    {            
        const DoubleReal x1(p1.getMZ());
        const DoubleReal y1(std::log(p1.getIntensity()));
        const DoubleReal x2(p2.getMZ());
        const DoubleReal y2(std::log(p2.getIntensity()));
        const DoubleReal x3(p3.getMZ());
        const DoubleReal y3(std::log(p3.getIntensity()));

        DoubleReal D = (x1-x2)*(x1-x3)*(x2-x3);
        DoubleReal alpha = (x3*(y2-y1) + x2*(y1-y3) + x1*(y3-y2)) / D;
        DoubleReal beta = (x3*x3*(y1-y2) + x2*x2*(y3-y1) + x1*x1*(y2-y3)) / D;
        DoubleReal gamma = (y1*x2*x3*(x2-x3) + y2*x3*x1*(x3-x1) + y3*x1*x2*(x1-x2)) / D;

        mu = -beta/(2.0*alpha);
        DoubleReal c_square = -1.0 / alpha;
        DoubleReal sigma_square = c_square / 2.0;
        height = std::exp(gamma + mu * mu / c_square);
        area = height / std::sqrt(2.0 * M_PI * sigma_square);
        sigma = std::sqrt(sigma_square);

        return (area != std::numeric_limits<DoubleReal>::infinity());
    }

    /**
    @brief Applies the peak-picking algorithm to a single spectrum (MSSpectrum). The resulting picked peaks are written to the output spectrum.
    */
    template <typename PeakType>
    void pick(const MSSpectrum<PeakType>& cinput, MSSpectrum<PeakType>& output) 
    {
        MSSpectrum<PeakType> input = cinput;
        threshold_mower_.filterPeakSpectrum(input);
        input.sortByPosition();

        // copy meta data of the input spectrum
        output.clear(true);
        output.SpectrumSettings::operator=(input);
        output.MetaInfoInterface::operator=(input);
        output.setRT(input.getRT());
        output.setMSLevel(input.getMSLevel());
        output.setName(input.getName());
        output.setType(SpectrumSettings::PEAKS);

        bool intensity_type_area = param_.getValue("intensity_type") == "peakarea" ? true : false;

        // find local maxima in raw data
        for (Size i = 2; i < input.size() - 2; ++i)
        {
            DoubleReal central_peak_mz = input[i].getMZ(), central_peak_int = input[i].getIntensity();

            DoubleReal l1_neighbor_mz = input[i-1].getMZ(), l1_neighbor_int = input[i-1].getIntensity();
            DoubleReal r1_neighbor_mz = input[i+1].getMZ(), r1_neighbor_int = input[i+1].getIntensity();

            DoubleReal l2_neighbor_mz = input[i-2].getMZ(), l2_neighbor_int = input[i-2].getIntensity();
            DoubleReal r2_neighbor_mz = input[i+2].getMZ(), r2_neighbor_int = input[i+2].getIntensity();

            // MZ spacing sanity checks
            DoubleReal l1_to_central = std::fabs(central_peak_mz - l1_neighbor_mz);
            DoubleReal l2_to_l1 = std::fabs(l1_neighbor_mz - l2_neighbor_mz);

            DoubleReal central_to_r1 = std::fabs(r1_neighbor_mz - central_peak_mz);
            DoubleReal r1_to_r2 = std::fabs(r2_neighbor_mz - r1_neighbor_mz);

            DoubleReal min_spacing = (l1_to_central < central_to_r1) ? l1_to_central : central_to_r1;
                
            // look for peak cores meeting MZ and intensity/SNT criteria
            if (central_peak_int > 1.0 && l1_neighbor_int > 1.0 && l2_neighbor_int > 1.0 && r1_neighbor_int > 1.0 && r2_neighbor_int > 1.0
                    && l1_to_central < 1.5*min_spacing
                    && l2_to_l1 < 1.5*min_spacing
                    && (l2_neighbor_int < l1_neighbor_int && l1_neighbor_int < central_peak_int)
                    && central_to_r1 < 1.5*min_spacing
                    && r1_to_r2 < 1.5*min_spacing
                    && (r2_neighbor_int < r1_neighbor_int && r1_neighbor_int < central_peak_int)
                    )
            {
                // potential triple
                DoubleReal mu(0.0), sigma(0.0), area(0.0), height(0.0);

                bool compOK = computeTPG(input[i - 1], input[i], input[i + 1], mu, sigma, area, height);

                // save picked pick into output spectrum
                if (compOK)
                {
                  PeakType peak;
                  peak.setMZ(mu);

                  DoubleReal output_intensity = intensity_type_area ? area : height;

                  peak.setIntensity(output_intensity);
                  output.push_back(peak);
                }

                // jump over raw data points that have been considered already
                i = i + 1;
            }
        }

        return;
    }

    /**
    @brief Applies the peak-picking algorithm to a map (MSExperiment). 
    
    This method picks peaks for each scan in the map consecutively. The
    resulting picked peaks are written to the output map.
    */
    template <typename PeakType>
    void pickExperiment(MSExperiment<PeakType>& input, MSExperiment<PeakType>& output)
    {
        // make sure that output is clear
        output.clear(true);

        // copy experimental settings
        static_cast<ExperimentalSettings&>(output) = input;

        // resize output with respect to input
        output.resize(input.size());

        bool ms1_only = param_.getValue("ms1_only").toBool();
        Size progress = 0;

        startProgress(0,input.size(),"picking peaks");
        for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
        {
            if (ms1_only && (input[scan_idx].getMSLevel() != 1))
            {
                output[scan_idx] = input[scan_idx];
            }
            else
            {
                pick(input[scan_idx], output[scan_idx]);
            }
            setProgress(++progress);
        }
        endProgress();

        return;
    }

protected:
    //docu in base class
    void updateMembers_();

    ThresholdMower threshold_mower_;
}; // end PeakPickerRapid

}// namespace OpenMS

#endif
