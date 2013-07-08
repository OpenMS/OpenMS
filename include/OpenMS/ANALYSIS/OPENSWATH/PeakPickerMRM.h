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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_PEAKPICKERMRM_H
#define OPENMS_ANALYSIS_OPENSWATH_PEAKPICKERMRM_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>

#ifdef WITH_CRAWDAD
#include <CrawdadWrapper.h>
#endif

namespace OpenMS
{

  /**
  @brief The PeakPickerMRM finds peaks a single chromatogram.

  @htmlinclude OpenMS_PeakPickerMRM.parameters

  It uses the PeakPickerHiRes internally to find interesting seed candidates.
  These candidates are then expanded and a right/left border of the peak is
  searched.
  Additionally, overlapping peaks can be removed.

  */

  class OPENMS_DLLAPI PeakPickerMRM :
    public DefaultParamHandler
  {

public:

    // this is the type in which we store the chromatograms for this analysis
    typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram; 

    //@{
    /// Constructor
    PeakPickerMRM();

    /// Destructor
    ~PeakPickerMRM() {}
    //@}

#ifdef WITH_CRAWDAD
    void pickChromatogramCrowdad(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom)
    {
      std::cout << " using crawdad " << std::endl;

      std::vector<double> time;
      std::vector<double> intensity;
      for (Size i = 0; i < chromatogram.size(); i++)
      {
        time.push_back(chromatogram[i].getRT());
        intensity.push_back(chromatogram[i].getIntensity());
      }

      CrawdadWrapper crawdad_pp;
      crawdad_pp.SetChromatogram(time, intensity);
      std::vector<crawpeaks::SlimCrawPeak> result = crawdad_pp.CalcPeaks();

      picked_chrom.getFloatDataArrays().clear();
      picked_chrom.getFloatDataArrays().resize(3);
      picked_chrom.getFloatDataArrays()[0].setName("IntegratedIntensity");
      picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
      picked_chrom.getFloatDataArrays()[2].setName("rightWidth");
      for(std::vector<crawpeaks::SlimCrawPeak>::iterator it = result.begin(); it != result.end(); it++)
      {

        ChromatogramPeak p;
        p.setRT( chromatogram[it->peak_rt_idx].getRT() );
        p.setIntensity( it->peak_area ); //chromatogram[it->peak_rt_idx].getIntensity() );

        picked_chrom.getFloatDataArrays()[0].push_back( it->peak_area );
        picked_chrom.getFloatDataArrays()[1].push_back( chromatogram[it->start_rt_idx].getRT() );
        picked_chrom.getFloatDataArrays()[2].push_back( chromatogram[it->stop_rt_idx].getRT() );
        /*
        int peak_rt_idx, start_rt_idx, stop_rt_idx, max_rt_idx;
        int mz_idx;
        int len;
        float fwhm;
        bool fwhm_calculated_ok;
        float bg_area;
        float raw_area; // total area under the curve, including background
        float peak_area;  
        float bgslope;
        ///cutoff level for extending past the peak

        ///maximum height, calculated above background
        float peak_height;
        float raw_height;

        */

        LOG_DEBUG << "Found peak at " << p.getRT() << " and "  << chromatogram[it->peak_rt_idx].getIntensity()
        << " with borders " << chromatogram[it->start_rt_idx].getRT() << " " << chromatogram[it->stop_rt_idx].getRT()  <<  " (" << chromatogram[it->start_rt_idx].getRT() - chromatogram[it->stop_rt_idx].getRT() << ") " 
        << it->peak_area << " weighted RT " << /* weighted_mz << */ std::endl;

        picked_chrom.push_back(p);
        
      }

    }
#else
    void pickChromatogramCrowdad(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
        "PeakPickerMRM was not compiled with crawdad, please choose a different algorithm!");
    }
#endif

    /**
      @brief Finds peaks in a single chromatogram and annotates left/right borders

      It uses a modified algorithm of the PeakPickerHiRes

      This function will return a smoothed chromatogram and a picked chromatogram
    */
    void pickChromatogram(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom);

    void pickAndSmoothChromatogram(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom);

protected:

    /**
      @brief Helper function to find the closest peak in a chromatogram to "target_rt" 

      The search will start from the index current_peak, so the function is
      assuming the closest peak is to the right of current_peak.

      It will return the index of the closest peak in the chromatogram.
    */
    Size findClosestPeak_(const RichPeakChromatogram& chromatogram, double target_rt, Size current_peak = 0);

    /**
      @brief Helper function to remove overlapping peaks in a single Chromatogram

    */
    void removeOverlappingPeaks_(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom);


    /// Synchronize members with param class
    void updateMembers_();

    /// Assignment operator is protected for algorithm
    PeakPickerMRM& operator=(const PeakPickerMRM& rhs);

    // Members
    UInt sgolay_frame_length_;
    UInt sgolay_polynomial_order_;
    DoubleReal gauss_width_;
    bool use_gauss_;
    bool remove_overlapping_;

    DoubleReal peak_width_;
    DoubleReal signal_to_noise_;

    DoubleReal sn_win_len_;
    UInt sn_bin_count_;
    String method_;
  };
}

#endif
