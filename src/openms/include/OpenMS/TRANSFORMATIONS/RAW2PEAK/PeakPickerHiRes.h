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
// $Maintainer: Erhan Kenar, Lars Nilse $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERHIRES_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERHIRES_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>


#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#undef DEBUG_DECONV
namespace OpenMS
{
  /**
    @brief This class implements a fast peak-picking algorithm best suited for
    high resolution MS data (FT-ICR-MS, Orbitrap). In high resolution data, the
    signals of ions with similar mass-to-charge ratios (m/z) exhibit little or
    no overlapping and therefore allow for a clear separation. Furthermore, ion
    signals tend to show well-defined peak shapes with narrow peak width.

    This peak-picking algorithm detects ion signals in raw data and
    reconstructs the corresponding peak shape by cubic spline interpolation.
    Signal detection depends on the signal-to-noise ratio which is adjustable
    by the user (see parameter signal_to_noise). A picked peak's m/z and
    intensity value is given by the maximum of the underlying peak spline.

    So far, this peak picker was mainly tested on high resolution data. With
    appropriate preprocessing steps (e.g. noise reduction and baseline
    subtraction), it might be also applied to low resolution data.

    @htmlinclude OpenMS_PeakPickerHiRes.parameters

    @note The peaks must be sorted according to ascending m/z!

    @ingroup PeakPicking
  */
  class OPENMS_DLLAPI PeakPickerHiRes :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Constructor
    PeakPickerHiRes();

    /// Destructor
    virtual ~PeakPickerHiRes();
    
    /// structure for peak boundaries
    struct PeakBoundary
    {
        double mzMin;
        double mzMax;
    };

    /**
     * @brief Applies the peak-picking algorithm to a single spectrum
     * (MSSpectrum). The resulting picked peaks are written to the output
     * spectrum.
     * 
     * @param input  input spectrum in profile mode
     * @param output  output spectrum with picked peaks
     */
    template <typename PeakType>
    void pick(const MSSpectrum<PeakType> & input, MSSpectrum<PeakType> & output) const
    {
        std::vector<PeakBoundary> boundaries;
        pick(input, output, boundaries, false);
    }

    /**
     * @brief Applies the peak-picking algorithm to a single spectrum
     * (MSSpectrum). The resulting picked peaks are written to the output
     * spectrum. Peak boundaries are written to a separate structure.
     * 
     * @param input  input spectrum in profile mode
     * @param output  output spectrum with picked peaks
     * @param boundaries  boundaries of the picked peaks
     * @param detectShoulds  optional flag for the detection of peak shoulders (extrema in first derivative)
     */
    template <typename PeakType>
    void pick(const MSSpectrum<PeakType> & input, MSSpectrum<PeakType> & output, std::vector<PeakBoundary> & boundaries, bool detect_shoulders) const
    {
      // copy meta data of the input spectrum
      output.clear(true);
      output.SpectrumSettings::operator=(input);
      output.MetaInfoInterface::operator=(input);
      output.setRT(input.getRT());
      output.setMSLevel(input.getMSLevel());
      output.setName(input.getName());
      output.setType(SpectrumSettings::PEAKS);

      // don't pick a spectrum with less than 5 data points
      if (input.size() < 5) return;

      // signal-to-noise estimation
      SignalToNoiseEstimatorMedian<MSSpectrum<PeakType> > snt;
      snt.setParameters(param_.copy("SignalToNoise:", true));

      if (signal_to_noise_ > 0.0)
      {
        snt.init(input);
      }

      // find local maxima in raw data
      for (Size i = 2; i < input.size() - 2; ++i)
      {
        double central_peak_mz = input[i].getMZ(), central_peak_int = input[i].getIntensity();
        double left_neighbor_mz = input[i - 1].getMZ(), left_neighbor_int = input[i - 1].getIntensity();
        double right_neighbor_mz = input[i + 1].getMZ(), right_neighbor_int = input[i + 1].getIntensity();

        //do not interpolate when the left or right support is a zero-data-point
        if(std::fabs(left_neighbor_int) < std::numeric_limits<double>::epsilon() )
          continue;
        if(std::fabs(right_neighbor_int) < std::numeric_limits<double>::epsilon() )
          continue;

        // MZ spacing sanity checks
        double left_to_central = std::fabs(central_peak_mz - left_neighbor_mz);
        double central_to_right = std::fabs(right_neighbor_mz - central_peak_mz);
        double min_spacing = (left_to_central < central_to_right) ? left_to_central : central_to_right;

        double act_snt = 0.0, act_snt_l1 = 0.0, act_snt_r1 = 0.0;

        if (signal_to_noise_ > 0.0)
        {
          act_snt = snt.getSignalToNoise(input[i]);
          act_snt_l1 = snt.getSignalToNoise(input[i - 1]);
          act_snt_r1 = snt.getSignalToNoise(input[i + 1]);
        }

        // look for peak cores meeting MZ and intensity/SNT criteria
        if (act_snt >= signal_to_noise_
           && left_to_central < spacing_difference_ * min_spacing
           && central_peak_int > left_neighbor_int
           && act_snt_l1 >= signal_to_noise_
           && central_to_right < spacing_difference_ * min_spacing
           && central_peak_int > right_neighbor_int
           && act_snt_r1 >= signal_to_noise_)
        {
          // special case: if a peak core is surrounded by more intense
          // satellite peaks (indicates oscillation rather than
          // real peaks) -> remove

          double act_snt_l2 = 0.0, act_snt_r2 = 0.0;

          if (signal_to_noise_ > 0.0)
          {
            act_snt_l2 = snt.getSignalToNoise(input[i - 2]);
            act_snt_r2 = snt.getSignalToNoise(input[i + 2]);
          }

          //checking signal-to-noise?
          if ((i > 1
              && std::fabs(left_neighbor_mz - input[i - 2].getMZ()) < spacing_difference_ * min_spacing
              && left_neighbor_int < input[i - 2].getIntensity()
              && act_snt_l2 >= signal_to_noise_)
             &&
              ((i + 2) < input.size()
              && std::fabs(input[i + 2].getMZ() - right_neighbor_mz) < spacing_difference_ * min_spacing
              && right_neighbor_int < input[i + 2].getIntensity()
              && act_snt_r2 >= signal_to_noise_)
              )
          {
            ++i;
            continue;
          }


          std::map<double, double> peak_raw_data;

          peak_raw_data[central_peak_mz] = central_peak_int;
          peak_raw_data[left_neighbor_mz] = left_neighbor_int;
          peak_raw_data[right_neighbor_mz] = right_neighbor_int;

          // peak core found, now extend it
          // to the left
          Size k = 2;

          Size missing_left(0);
          Size left_boundary(i-1);    // index of the left boundary for the spline interpolation

          while ( k <= i//prevent underflow
                && (i - k + 1) > 0
                && (missing_left < 2)
                && input[i - k].getIntensity() <= peak_raw_data.begin()->second)
          {

            double act_snt_lk = 0.0;

            if (signal_to_noise_ > 0.0)
            {
              act_snt_lk = snt.getSignalToNoise(input[i - k]);
            }


            if (act_snt_lk >= signal_to_noise_ && std::fabs(input[i - k].getMZ() - peak_raw_data.begin()->first) < spacing_difference_ * min_spacing)
            {
              peak_raw_data[input[i - k].getMZ()] = input[i - k].getIntensity();
            }
            else
            {
              peak_raw_data[input[i - k].getMZ()] = input[i - k].getIntensity();
              ++missing_left;
            }

            left_boundary = i - k;
            ++k;

          }

          // to the right
          k = 2;
          
          Size missing_right(0);
          Size right_boundary(i+1);    // index of the left boundary for the spline interpolation
          
          while ((i + k) < input.size()
                && (missing_right < 2)
                && input[i + k].getIntensity() <= peak_raw_data.rbegin()->second)
          {

            double act_snt_rk = 0.0;

            if (signal_to_noise_ > 0.0)
            {
              act_snt_rk = snt.getSignalToNoise(input[i + k]);
            }

            if (act_snt_rk >= signal_to_noise_ && std::fabs(input[i + k].getMZ() - peak_raw_data.rbegin()->first) < spacing_difference_ * min_spacing)
            {
              peak_raw_data[input[i + k].getMZ()] = input[i + k].getIntensity();
            }
            else
            {
              peak_raw_data[input[i + k].getMZ()] = input[i + k].getIntensity();
              ++missing_right;
            }

            right_boundary = i + k;
            ++k;
          }

          //skip if the minimal number of 3 points for fitting is not reached
          if(peak_raw_data.size() < 4)
            continue;

          CubicSpline2d peak_spline (peak_raw_data);

          // calculate maximum by evaluating the spline's 1st derivative
          // (bisection method)
          double max_peak_mz = central_peak_mz;
          double max_peak_int = central_peak_int;
          double threshold = 0.000001;
          double lefthand = left_neighbor_mz;
          double righthand = right_neighbor_mz;

          bool lefthand_sign = 1;
          double eps = std::numeric_limits<double>::epsilon();

          // bisection
          do
          {
            double mid = (lefthand + righthand) / 2;
            double midpoint_deriv_val = peak_spline.derivatives(mid, 1);

            // if deriv nearly zero then maximum already found
            if (!(std::fabs(midpoint_deriv_val) > eps))
            {
              break;
            }

            bool midpoint_sign = (midpoint_deriv_val < 0.0) ? 0 : 1;

            if (lefthand_sign ^ midpoint_sign)
            {
              righthand = mid;
            }
            else
            {
              lefthand = mid;
            }
          }
          while (std::fabs(lefthand - righthand) > threshold);

          // sanity check?
          max_peak_mz = (lefthand + righthand) / 2;
          max_peak_int = peak_spline.eval( max_peak_mz );

          if (detect_shoulders)
          {
              std::vector<double> peak_raw_data_mz;
              std::vector<double> peak_raw_data_int;
              for (std::map< double, double, std::less< double > >::const_iterator iter = peak_raw_data.begin(); iter != peak_raw_data.end(); ++iter)
              {
                peak_raw_data_mz.push_back(iter->first);
                peak_raw_data_int.push_back(iter->second);
              }
              
              int j(0);
              while (peak_raw_data_mz[j] < max_peak_mz)
              {
                  ++j;
              }
              
              // Is there a shoulder, i.e. a maximum in the first derivative, to the right?
              int rs(0);
              for (int l = j+1; l < (int) peak_raw_data_mz.size()-1 -1; ++l)    // Do not trust shoulders close to the borders, hence -1.
              {
                  double left = peak_spline.derivatives(peak_raw_data_mz[l-1], 1);
                  double centre = peak_spline.derivatives(peak_raw_data_mz[l], 1);
                  double right = peak_spline.derivatives(peak_raw_data_mz[l+1], 1);
                  if (left < centre && centre > right)
                  {
                      rs = l;
                      break;    // We consider only shoulders right next to the peak maximum. Multiple shoulders on one side are ignored.
                  }
              }

              // Is there a shoulder, i.e. a minimum in the first derivative, to the left?
              int ls(0);
              for (int l = j-1; l > 0; --l)    // Do not trust shoulders close to the borders, hence >0.
              {
                  double left = peak_spline.derivatives(peak_raw_data_mz[l-1], 1);
                  double centre = peak_spline.derivatives(peak_raw_data_mz[l], 1);
                  double right = peak_spline.derivatives(peak_raw_data_mz[l+1], 1);
                  if (left > centre && centre < right)
                  {
                      ls = l;
                      break;    // We consider only shoulders right next to the peak maximum. Multiple shoulders on one side are ignored.
                  }
              }

              PeakType peak;
              PeakBoundary peak_boundary;

              // add central peak to results
              peak.setMZ(max_peak_mz);
              peak.setIntensity(max_peak_int);
              output.push_back(peak);
              // Peak boundaries are where the (possible) shoulders start.
              if (ls==0)
              {
                  peak_boundary.mzMin = input[left_boundary].getMZ();
              }
              else
              {
                  peak_boundary.mzMin = (peak_raw_data_mz[ls]+peak_raw_data_mz[ls+1])/2;
              }
              if (rs==0)
              {
                  peak_boundary.mzMax = input[right_boundary].getMZ();
              }
              else
              {
                  peak_boundary.mzMax = (peak_raw_data_mz[rs]+peak_raw_data_mz[rs-1])/2;
              }
              boundaries.push_back(peak_boundary);
              
              // add left peak to results
              if (ls!=0)
              {
                  // m/z centre of the shoulder is the intensity-weighted average of the m/z
                  double mz(0);
                  double summed_intensities(0);
                  for (int s=0; s<=ls; ++s)
                  {
                      mz += peak_raw_data_mz[s]*peak_raw_data_int[s];
                      summed_intensities += peak_raw_data_int[s];
                  }
                  //peak.setMZ(mz/summed_intensities);
                  // intensity of the shoulder is the spline interpolation at the m/z centre
                  //peak.setIntensity(std::max(0.0,peak_spline.eval(mz/summed_intensities)));
                  peak.setMZ(peak_raw_data_mz[ls]);
                  peak.setIntensity(peak_raw_data_int[ls]);
                  output.push_back(peak);
                  
                  peak_boundary.mzMin = input[left_boundary].getMZ();
                  peak_boundary.mzMax = (peak_raw_data_mz[ls]+peak_raw_data_mz[ls+1])/2;
                  boundaries.push_back(peak_boundary);
              }

              // add right peak to results
              if (rs!=0)
              {
                  // m/z centre of the shoulder is the intensity-weighted average of the m/z
                  double mz(0);
                  double summed_intensities(0);
                  for (int s = rs; s < (int) peak_raw_data_mz.size(); ++s)
                  {
                      mz += peak_raw_data_mz[s]*peak_raw_data_int[s];
                      summed_intensities += peak_raw_data_int[s];
                  }
                  //peak.setMZ(mz/summed_intensities);
                  // intensity of the shoulder is the spline interpolation at the m/z centre
                  //peak.setIntensity(std::max(0.0,peak_spline.eval(mz/summed_intensities)));
                  peak.setMZ(peak_raw_data_mz[rs]);
                  peak.setIntensity(peak_raw_data_int[rs]);
                  output.push_back(peak);
                  
                  peak_boundary.mzMin = (peak_raw_data_mz[rs]+peak_raw_data_mz[rs-1])/2;
                  peak_boundary.mzMax = input[right_boundary].getMZ();
                  boundaries.push_back(peak_boundary);
              }

          }
          else
          {
              // save picked pick into output spectrum
              PeakType peak;
              PeakBoundary peak_boundary;
              
              peak.setMZ(max_peak_mz);
              peak.setIntensity(max_peak_int);
              output.push_back(peak);
              peak_boundary.mzMin = input[left_boundary].getMZ();
              peak_boundary.mzMax = input[right_boundary].getMZ();
              boundaries.push_back(peak_boundary);
         }

          // jump over raw data points that have been considered already
          i = i + k - 1;
        }
      }

      return;
    }

     /**
     * @brief Applies the peak-picking algorithm to a single chromatogram
     * (MSChromatogram). The resulting picked peaks are written to the output chromatogram.
     * 
     * @param input  input chromatogram in profile mode
     * @param output  output chromatogram with picked peaks
     */
    template <typename PeakType>
    void pick(const MSChromatogram<PeakType> & input, MSChromatogram<PeakType> & output) const
    {
        std::vector<PeakBoundary> boundaries;
        pick(input, output, boundaries, false);
    }
    
    /**
     * @brief Applies the peak-picking algorithm to a single chromatogram
     * (MSChromatogram). The resulting picked peaks are written to the output chromatogram.
     * 
     * @param input  input chromatogram in profile mode
     * @param output  output chromatogram with picked peaks
     * @param boundaries  boundaries of the picked peaks
     * @param detectShoulds  optional flag for the detection of peak shoulders (extrema in first derivative)
     */
    template <typename PeakType>
    void pick(const MSChromatogram<PeakType> & input, MSChromatogram<PeakType> & output, std::vector<PeakBoundary> & boundaries, bool shoulders) const
    {
      // copy meta data of the input chromatogram
      output.clear(true);
      output.ChromatogramSettings::operator=(input);
      output.MetaInfoInterface::operator=(input);
      output.setName(input.getName());

      MSSpectrum<PeakType> input_spectrum;
      MSSpectrum<PeakType> output_spectrum;
      for (typename MSChromatogram<PeakType>::const_iterator it = input.begin(); it != input.end(); ++it)
      {
        input_spectrum.push_back(*it);
      }
      pick(input_spectrum, output_spectrum, boundaries, shoulders);
      for (typename MSSpectrum<PeakType>::const_iterator it = output_spectrum.begin(); it != output_spectrum.end(); ++it)
      {
        output.push_back(*it);
      }

    }

    /**
     * @brief Applies the peak-picking algorithm to a map (MSExperiment). This
     * method picks peaks for each scan in the map consecutively. The resulting
     * picked peaks are written to the output map.
     * 
     * @param input  input map in profile mode
     * @param output  output map with picked peaks
     */
    template <typename PeakType, typename ChromatogramPeakT>
    void pickExperiment(const MSExperiment<PeakType, ChromatogramPeakT> & input, MSExperiment<PeakType, ChromatogramPeakT> & output) const
    {
        std::vector<std::vector<PeakBoundary> > boundaries_spec;
        std::vector<std::vector<PeakBoundary> > boundaries_chrom;
        pickExperiment(input, output, boundaries_spec, boundaries_chrom);
    }

    /**
     * @brief Applies the peak-picking algorithm to a map (MSExperiment). This
     * method picks peaks for each scan in the map consecutively. The resulting
     * picked peaks are written to the output map.
     * 
     * @param input  input map in profile mode
     * @param output  output map with picked peaks
     * @param boundaries_spec  boundaries of the picked peaks in spectra
     * @param boundaries_chrom  boundaries of the picked peaks in chromatograms
     * @param detectShoulds  optional flag for the detection of peak shoulders (extrema in first derivative)
     */
    template <typename PeakType, typename ChromatogramPeakT>
    void pickExperiment(const MSExperiment<PeakType, ChromatogramPeakT> & input, MSExperiment<PeakType, ChromatogramPeakT> & output, std::vector<std::vector<PeakBoundary> > & boundaries_spec, std::vector<std::vector<PeakBoundary> > & boundaries_chrom) const
    {
      // make sure that output is clear
      output.clear(true);

      // copy experimental settings
      static_cast<ExperimentalSettings &>(output) = input;

      // resize output with respect to input
      output.resize(input.size());

      bool ms1_only = param_.getValue("ms1_only").toBool();
      bool detect_shoulders = param_.getValue("detect_shoulders").toBool();
      Size progress = 0;

      startProgress(0, input.size() + input.getChromatograms().size(), "picking peaks");
      for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
      {
        if (ms1_only && (input[scan_idx].getMSLevel() != 1))
        {
          output[scan_idx] = input[scan_idx];
        }
        else
        {
          std::vector<PeakBoundary> boundaries_s;    // peak boundaries of a single spectrum
          pick(input[scan_idx], output[scan_idx], boundaries_s, detect_shoulders);
          boundaries_spec.push_back(boundaries_s);
        }
        setProgress(++progress);
      }
      for (Size i = 0; i < input.getChromatograms().size(); ++i)
      {
        MSChromatogram<ChromatogramPeakT> chromatogram;
        std::vector<PeakBoundary> boundaries_c;    // peak boundaries of a single chromatogram
        pick(input.getChromatograms()[i], chromatogram, boundaries_c, detect_shoulders);
        output.addChromatogram(chromatogram);
        boundaries_chrom.push_back(boundaries_c);
        setProgress(++progress);
      }

      endProgress();

      return;
    }

    /**
      @brief Applies the peak-picking algorithm to a map (MSExperiment). This
      method picks peaks for each scan in the map consecutively. The resulting
      picked peaks are written to the output map.

      Currently we have to give up const-correctness but we know that everything on disc is constant
    */
    template <typename PeakType, typename ChromatogramPeakT>
    void pickExperiment(/* const */ OnDiscMSExperiment<PeakType, ChromatogramPeakT> & input, MSExperiment<PeakType, ChromatogramPeakT> & output) const
    {
      // make sure that output is clear
      output.clear(true);

      // copy experimental settings
      static_cast<ExperimentalSettings &>(output) = *input.getExperimentalSettings();

      // resize output with respect to input
      output.resize(input.size());

      bool ms1_only = param_.getValue("ms1_only").toBool();
      Size progress = 0;

      startProgress(0, input.size() + input.getNrChromatograms(), "picking peaks");
      for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
      {
        if (ms1_only && (input[scan_idx].getMSLevel() != 1))
        {
          output[scan_idx] = input[scan_idx];
        }
        else
        {
          MSSpectrum<PeakType> s = input[scan_idx];
          s.sortByPosition();
          pick(s, output[scan_idx]);
        }
        setProgress(++progress);
      }
      for (Size i = 0; i < input.getNrChromatograms(); ++i)
      {
        MSChromatogram<ChromatogramPeakT> chromatogram;
        pick(input.getChromatogram(i), chromatogram);
        output.addChromatogram(chromatogram);
        setProgress(++progress);
      }

      endProgress();

      return;
    }

protected:
    // signal-to-noise parameter
    double signal_to_noise_;

    // maximal spacing difference
    double spacing_difference_;

    // docu in base class
    void updateMembers_();

  }; // end PeakPickerHiRes

} // namespace OpenMS

#endif
