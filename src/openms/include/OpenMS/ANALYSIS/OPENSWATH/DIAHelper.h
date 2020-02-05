// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Witold Wolski, Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>

#include <vector>

namespace OpenMS
{
  class TheoreticalSpectrumGenerator;
  namespace DIAHelpers
  {

    /**
      @brief Helper functions for the DIA scoring of OpenSWATH
    */
    ///@{

    /**
      @brief Integrate intensity in a spectrum from start to end

      This function will integrate the intensity in a spectrum between mz_start
      and mz_end, returning the total intensity and an intensity-weighted m/z
      value.

      @note If there is no signal, mz will be set to -1 and intensity to 0
      @return Returns true if a signal was found (and false if no signal was found)

    */
    OPENMS_DLLAPI bool integrateWindow(const OpenSwath::SpectrumPtr spectrum, double mz_start,
                                       double mz_end, double& mz, double& intensity, bool centroided = false);

    /**
      @brief Integrate intensities in a spectrum from start to end
    */
    OPENMS_DLLAPI void integrateWindows(const OpenSwath::SpectrumPtr spectrum, //!< [in] Spectrum
                                        const std::vector<double>& windows_center, //!< [in] center location
                                        double width,
                                        std::vector<double>& integrated_windows_intensity,
                                        std::vector<double>& integrated_windows_mz,
                                        bool remove_zero = false);

    /**
      @brief Integrate intensity in an ion mobility spectrum from start to end

      This function will integrate the intensity in a spectrum between mz_start
      and mz_end, returning the total intensity and an intensity-weighted drift
      time value.

      @note If there is no signal, mz will be set to -1 and intensity to 0
    */
    OPENMS_DLLAPI void integrateDriftSpectrum(OpenSwath::SpectrumPtr spectrum,
                                              double mz_start,
                                              double mz_end,
                                              double & im,
                                              double & intensity,
                                              double drift_start,
                                              double drift_end);

    /**
      @brief Adjust left/right window based on window and whether its ppm or not
    */
    OPENMS_DLLAPI void adjustExtractionWindow(double& right, double& left, const double& mz_extract_window, const bool& mz_extraction_ppm);

    /// compute the b and y series masses for a given AASequence
    OPENMS_DLLAPI void getBYSeries(const AASequence& a,
                                   std::vector<double>& bseries,
                                   std::vector<double>& yseries,
                                   TheoreticalSpectrumGenerator const * g,
                                   UInt charge = 1u);

    /// for SWATH -- get the theoretical b and y series masses for a sequence
    OPENMS_DLLAPI void getTheorMasses(const AASequence& a,
                                      std::vector<double>& masses,
                                      TheoreticalSpectrumGenerator const * g,
                                      UInt charge = 1u);

    /// get averagine distribution given mass
    OPENMS_DLLAPI void getAveragineIsotopeDistribution(const double product_mz,
                                         std::vector<std::pair<double, double> >& isotopesSpec, 
                                         const double charge = 1.,
                                         const int nr_isotopes = 4,
                                         const double mannmass = 1.00048);

    /// simulate spectrum from AASequence
    OPENMS_DLLAPI void simulateSpectrumFromAASequence(const AASequence& aa,
                                        std::vector<double>& firstIsotopeMasses, //[out]
                                        std::vector<std::pair<double, double> >& isotopeMasses, //[out]
                                        TheoreticalSpectrumGenerator const * g,
                                        double charge = 1.);

    /// modify masses by charge
    OPENMS_DLLAPI void modifyMassesByCharge(const std::vector<std::pair<double, double> >& masses, //![in]
                              std::vector<std::pair<double, double> >& modmass, //!< [out]
                              double charge = 1.);

    /// add negative pre-isotope weights to spectrum
    OPENMS_DLLAPI void addPreisotopeWeights(const std::vector<double>& firstIsotopeMasses,
                              std::vector<std::pair<double, double> >& isotopeSpec, // output
                              UInt nrpeaks = 2, //nr of pre-isotope peaks
                              double preIsotopePeaksWeight = -0.5, // weight of pre-isotope peaks
                              double mannmass = 1.000482, //
                              double charge = 1.);

    /// given an experimental spectrum add isotope pattern.
    OPENMS_DLLAPI void addIsotopes2Spec(const std::vector<std::pair<double, double> >& spec,
                          std::vector<std::pair<double, double> >& isotopeMasses, //[out]
                          double charge = 1.);

    /// sorts vector of pairs by first
    OPENMS_DLLAPI void sortByFirst(std::vector<std::pair<double, double> >& tmp);
    /// extract first from vector of pairs
    OPENMS_DLLAPI void extractFirst(const std::vector<std::pair<double, double> >& peaks, std::vector<double>& mass);
    /// extract second from vector of pairs
    OPENMS_DLLAPI void extractSecond(const std::vector<std::pair<double, double> >& peaks, std::vector<double>& mass);
    
    ///}@
  }
}

