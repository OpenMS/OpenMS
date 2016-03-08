// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Chris Bielow $
// $Authors: Alexandra Zerck, Chris Bielow $
// --------------------------------------------------------------------------


#ifndef OPENMS_FILTERING_CALIBRATION_INTERNALCALIBRATION_H
#define OPENMS_FILTERING_CALIBRATION_INTERNALCALIBRATION_H
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
// TODO remove dependency from file reader here!
#include <OpenMS/FORMAT/TransformationXMLFile.h>

namespace OpenMS
{

  /**
    @brief A simple calibration method using linear interpolation of given reference masses.

    This class implements a simple calibration method: given a list of reference masses,
    the relative errors of the peaks in the data are approximated by linear interpolation and
    subtracted from the data.

    @htmlinclude OpenMS_InternalCalibration.parameters

    @ingroup SignalProcessing
  */
  class OPENMS_DLLAPI InternalCalibration :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    InternalCalibration();

    /// Destructor
    ~InternalCalibration(){}

    /**
      @brief Calibrate a peak map using given reference masses with a separate calibration function for each spectrum.

      The calibration function is calculated for each spectrum
      separately. If not enough reference masses are found for a spectrum it is left uncalibrated.
      For the matching of the reference masses and the peaks, the parameter mz_tolerance is used to
      calculate a window around the reference masses. If more than one peak is found within this window the
      closest peak is taken.

      @param exp The peak map to calibrate
      @param ref_masses The reference m/z values
    */
    template <typename InputPeakType>
    void calibrateMapSpectrumwise(MSExperiment<InputPeakType>& exp, std::vector<double>& ref_masses)
    {
#ifdef DEBUG_CALIBRATION
      std::cout.precision(writtenDigits<double>(0.0));
#endif
      if (exp.empty())
      {
        std::cout << "Input is empty." << std::endl;
        return;
      }

      if (exp[0].getType() != SpectrumSettings::PEAKS)
      {
        std::cout << "Attention: this function is assuming peak data." << std::endl;
      }

      Size num_ref_peaks = ref_masses.size();
      bool use_ppm = (param_.getValue("mz_tolerance_unit") == "ppm") ? true : false;
      double mz_tol = param_.getValue("mz_tolerance");
      startProgress(0, exp.size(), "calibrate spectra");
      // for each spectrum
      for (Size spec = 0; spec < exp.size(); ++spec)
      {
        // calibrate only MS1 spectra
        if (exp[spec].getMSLevel() != 1)
        {
          continue;
        }


        std::vector<double> corr_masses, found_ref_masses;
        UInt corr_peaks = 0;
        for (Size peak = 0; peak <  exp[spec].size(); ++peak)
        {
          for (Size ref_peak = 0; ref_peak < num_ref_peaks; ++ref_peak)
          {
            if (!use_ppm && fabs(exp[spec][peak].getMZ() - ref_masses[ref_peak]) <  mz_tol)
            {
              found_ref_masses.push_back(ref_masses[ref_peak]);
              corr_masses.push_back(exp[spec][peak].getMZ());
              ++corr_peaks;
              break;
            }
            else if (use_ppm && fabs(exp[spec][peak].getMZ() - ref_masses[ref_peak]) / ref_masses[ref_peak] * 1e6 <  mz_tol)
            {
              found_ref_masses.push_back(ref_masses[ref_peak]);
              corr_masses.push_back(exp[spec][peak].getMZ());
              ++corr_peaks;
              break;
            }
          }
        }
        if (corr_peaks < 2)
        {
          std::cout << "spec: " << spec
            << " less than 2 reference masses were detected within a reasonable error range\n";
          std::cout << "This spectrum cannot be calibrated!\n";
          continue;
        }

        makeLinearRegression_(corr_masses, found_ref_masses);

        // now calibrate the whole spectrum
        applyTransformation_<InputPeakType>(exp[spec]);
        setProgress(spec);
      }  // for(Size spec=0;spec <  exp.size(); ++spec)
      endProgress();
    }

    /**
     @brief Calibrate a peak map using given reference masses with one calibration function for the whole map.

     The calibration function is calculated for the whole map.
     For the matching of the reference masses and the peaks the parameter mz_tolerance is used to
     calculate a window around the reference masses. If more than one peak is found within this window the
     closest peak is taken.

     @param exp The peak map to calibrate
     @param ref_masses The reference m/z values
     @param trafo_file_name Output file to store the transformation function
    */
    template <typename InputPeakType>
    void calibrateMapGlobally(MSExperiment<InputPeakType>& exp, 
      const std::vector<double>& ref_masses,
      const String& trafo_file_name = "")
    {
      if (exp.empty())
      {
        std::cout << "Input is empty." << std::endl;
        return;
      }

      if (exp[0].getType() != SpectrumSettings::PEAKS)
      {
        std::cout << "Attention: this function is assuming peak data." << std::endl;
      }

      Size num_ref_peaks = ref_masses.size();
      bool use_ppm = (param_.getValue("mz_tolerance_unit") == "ppm") ? true : false;
      double mz_tol = param_.getValue("mz_tolerance");
      startProgress(0, exp.size(), "finding calibration masses");
      std::vector<double> corr_masses, found_ref_masses;
      UInt corr_peaks = 0;
      // for each spectrum
      for (Size spec = 0; spec < exp.size(); ++spec)
      {
        // obtain calibration points only from MS1 spectra
        if (exp[spec].getMSLevel() != 1)
        {
          continue;
        }
        for (Size peak = 0; peak < exp[spec].size(); ++peak)
        {
          for (Size ref_peak = 0; ref_peak < num_ref_peaks; ++ref_peak)
          {
            if (!use_ppm &&  fabs(exp[spec][peak].getMZ() - ref_masses[ref_peak]) <  mz_tol)
            {
              found_ref_masses.push_back(ref_masses[ref_peak]);
              corr_masses.push_back(exp[spec][peak].getMZ());
              ++corr_peaks;
              break;
            }
            else if (use_ppm &&  fabs(exp[spec][peak].getMZ() - ref_masses[ref_peak]) / ref_masses[ref_peak] * 1e6 <  mz_tol)
            {
              found_ref_masses.push_back(ref_masses[ref_peak]);
              corr_masses.push_back(exp[spec][peak].getMZ());
              ++corr_peaks;
              break;
            }
          }
        }
        setProgress(spec);
      }
      endProgress();

      if (corr_peaks < 2)
      {
        std::cout << "Less than 2 reference masses were detected within a reasonable error range\n";
        std::cout << "This spectrum cannot be calibrated!\n";
        return;
      }

      // calculate the (linear) calibration function
      makeLinearRegression_(corr_masses, found_ref_masses);

      applyTransformation_<InputPeakType>(exp);

      if (trafo_file_name != "")
      {
        TransformationXMLFile().store(trafo_file_name, trafo_);
      }
    }

    /**
     @brief Calibrate a peak map using given reference ids with one calibration function for the whole map.

     Calibrate a map using given peptide identifications. The calibration function is calculated for the whole map.
     The m/z-values of the reference identifications are calculated through the given sequence and charge of the peptide.
     For the matching of the reference masses and the peaks the parameter mz_tolerance is used to
     calculate a window around the reference masses. If more than one peak is found within this window the
     closest peak is taken.

     @param exp The peak map to calibrate
     @param ref_ids The reference peptide identifications
     @param trafo_file_name Output file to store the transformation function
    */
    template <typename InputPeakType>
    void calibrateMapGlobally(MSExperiment<InputPeakType>& exp,
      const std::vector<PeptideIdentification>& ref_ids,
      const String& trafo_file_name = "")
    {
      bool use_ppm = param_.getValue("mz_tolerance_unit") == "ppm" ? true : false;
      double mz_tolerance = param_.getValue("mz_tolerance");
      if (exp.empty())
      {
        std::cout << "Input is empty." << std::endl;
        return;
      }

      if (exp[0].getType() != SpectrumSettings::PEAKS)
      {
        std::cout << "Attention: this function is assuming peak data." << std::endl;
      }
      // check if the ids contain meta information about the peak positions
      checkReferenceIds_(ref_ids);

      std::vector<double> theoretical_masses, observed_masses;
      for (Size p_id = 0; p_id < ref_ids.size(); ++p_id)
      {
        for (Size p_h = 0; p_h < ref_ids[p_id].getHits().size(); ++p_h)
        {
          Int charge = ref_ids[p_id].getHits()[p_h].getCharge();
          double theo_mass = ref_ids[p_id].getHits()[p_h].getSequence().getMonoWeight(Residue::Full, charge) / (double)charge;
          // first find corresponding ms1-spectrum
          typename MSExperiment<InputPeakType>::ConstIterator rt_iter = exp.RTBegin(ref_ids[p_id].getRT());
          while (rt_iter != exp.begin() && rt_iter->getMSLevel() != 1)
          {
            --rt_iter;
          }
          // now find closest peak
          typename MSSpectrum<InputPeakType>::ConstIterator mz_iter = rt_iter->MZBegin(ref_ids[p_id].getMZ());
          //std::cout << mz_iter->getMZ() <<" "<<(double)ref_ids[p_id].getMZ()<<"\t";
          double dist = ref_ids[p_id].getMZ() - mz_iter->getMZ();
          //std::cout << dist << "\t";
          if ((mz_iter + 1) != rt_iter->end()
            && fabs((mz_iter + 1)->getMZ() - ref_ids[p_id].getMZ()) < fabs(dist)
            && mz_iter != rt_iter->begin()
            && fabs((mz_iter - 1)->getMZ() - ref_ids[p_id].getMZ()) < fabs((mz_iter + 1)->getMZ() - ref_ids[p_id].getMZ()))  // if mz_iter +1 has smaller dist than mz_iter and mz_iter-1
          {
            if ((use_ppm &&
              fabs((mz_iter + 1)->getMZ() - ref_ids[p_id].getMZ()) / ref_ids[p_id].getMZ() * 1e06 < mz_tolerance) ||
              (!use_ppm && fabs((mz_iter + 1)->getMZ() - ref_ids[p_id].getMZ()) < mz_tolerance))
            {
              //std::cout <<(mz_iter +1)->getMZ() - ref_ids[p_id].getMZ()<<"\t";
              observed_masses.push_back((mz_iter + 1)->getMZ());
              theoretical_masses.push_back(theo_mass);
              //std::cout << (mz_iter +1)->getMZ() << " ~ "<<theo_mass << " charge: "<<ref_ids[p_id].getHits()[p_h].getCharge()
              //<< "\tplus 1"<< std::endl;
            }
          }
          else if (mz_iter != rt_iter->begin()
            && fabs((mz_iter - 1)->getMZ() - ref_ids[p_id].getMZ()) < fabs(dist))                        // if mz_iter-1 has smaller dist than mz_iter
          {
            if ((use_ppm &&
              fabs((mz_iter - 1)->getMZ() - ref_ids[p_id].getMZ()) / ref_ids[p_id].getMZ() * 1e06 < mz_tolerance) ||
              (!use_ppm && fabs((mz_iter - 1)->getMZ() - ref_ids[p_id].getMZ()) < mz_tolerance))
            {
              //std::cout <<(mz_iter -1)->getMZ() - ref_ids[p_id].getMZ()<<"\t";
              observed_masses.push_back((mz_iter - 1)->getMZ());
              theoretical_masses.push_back(theo_mass);
              //std::cout << (mz_iter -1)->getMZ() << " ~ "<<theo_mass << " charge: "<<ref_ids[p_id].getHits()[p_h].getCharge()
              //<< "\tminus 1"<< std::endl;
            }
          }
          else
          {
            if ((use_ppm &&
              fabs((mz_iter)->getMZ() - ref_ids[p_id].getMZ()) / ref_ids[p_id].getMZ() * 1e06 < mz_tolerance) ||
              (!use_ppm && fabs((mz_iter)->getMZ() - ref_ids[p_id].getMZ()) < mz_tolerance))
            {

              observed_masses.push_back(mz_iter->getMZ());
              theoretical_masses.push_back(theo_mass);
              //                                      std::cout <<"\t"<< mz_iter->getMZ() << " ~ "<<theo_mass<< " charge: "<<ref_ids[p_id].getHits()[p_h].getCharge()
              //                                                          << "\tat mz_iter"<< std::endl;
            }
          }
        }
      }

      makeLinearRegression_(observed_masses, theoretical_masses);

      applyTransformation_<InputPeakType>(exp);

      if (trafo_file_name != "")
      {
        TransformationXMLFile().store(trafo_file_name, trafo_);
      }
    }

    /**
     @brief Calibrate an annotated feature map with one calibration function for the whole map.

     Calibrate an annotated (!) feature map using the features' identifications. The calibration function is calculated for the whole map.
     The m/z-values of the reference identifications are calculated through the given sequence and charge of the peptide.
    
     @param feature_map The feature map to calibrate (annotated with peptide ids)
     @param trafo_file_name Output file to store the transformation function
    */
    void calibrateMapGlobally(FeatureMap& feature_map, const String& trafo_file_name = "");

    /**
     @brief Calibrate a feature map using given reference ids with one calibration function for the whole map.

     Calibrate a feature map using given peptide identifications. The calibration function is calculated for the whole map.
     Even if the features are already annotated with peptide ids these annotations are ignored for the calibration, only the reference ids are used.
     The m/z-values of the reference identifications are calculated through the given sequence and charge of the peptide.
     The reference ids are mapped onto the FeatureMap using IDMapper with the mz_tolerance and rt_tolerance parameters.

     @param feature_map The feature map to calibrate
     @param ref_ids the reference peptide identifications
     @param trafo_file_name file where the transformation function of the calibration is stored
    */
    void calibrateMapGlobally(FeatureMap& feature_map, std::vector<PeptideIdentification>& ref_ids, const String& trafo_file_name = "");


protected:

    /// the actual calibration function
    void makeLinearRegression_(const std::vector<double>& observed_masses, const std::vector<double>& theoretical_masses);

    /// check if reference ids contain RT and MZ information as meta values
    void checkReferenceIds_(const std::vector<PeptideIdentification>& pep_ids);

    /// check if reference ids contain RT and MZ information as meta values
    void checkReferenceIds_(const FeatureMap& feature_map);

    /// apply transformation to all features (including subordinates and convex hulls)
    void applyTransformation_(FeatureMap& feature_map);


    template <typename InputPeakType>
    void applyTransformation_(typename MSExperiment<InputPeakType>::SpectrumType& spec)
    {
      // calibrate only MS1 spectra
      if (spec.getMSLevel() != 1)
      {
        return;
      }

      for (unsigned int peak = 0; peak <  spec.size(); ++peak)
      {
#ifdef DEBUG_CALIBRATION
        std::cout << spec[peak].getMZ() << "\t";
#endif
        double mz = spec[peak].getMZ();
        mz = trafo_.apply(mz);
        spec[peak].setMZ(mz);

#ifdef DEBUG_CALIBRATION
        std::cout << spec[peak].getMZ() << std::endl;
#endif

      }
    }

    template <typename InputPeakType>
    void applyTransformation_(MSExperiment<InputPeakType>& exp)
    {
      startProgress(0, exp.size(), "applying calibration to data");
      // apply the calibration function to each peak
      for (Size spec = 0; spec < exp.size(); ++spec)
      {
        applyTransformation_<InputPeakType>(exp[spec]);
        setProgress(spec);
      }  // for(Size spec=0;spec <  exp.size(); ++spec)
      endProgress();

    }

    /// here the transformation is stored
    TransformationDescription trafo_;
  }; // class InternalCalibration
  
} // namespace OpenMS

#endif // OPENMS_FILTERING_CALIBRATION_INTERNALCALIBRATION_H
