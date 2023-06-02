//--------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <iomanip>

namespace OpenMS
{
  class PeakGroup;

  /**
       @brief A class representing a deconvolved spectrum.
       DeconvolvedSpectrum consists of PeakGroups representing masses.
       For MSn n>1, a PeakGroup representing the precursor mass is also added in this class. Properly assigning a precursor mass
       from the original precursor peak and its deconvolution result is very important in top down proteomics. This assignment is
       performed here for conventional datasets. But for FLASHIda acquired datasets, the assignment is already done by FLASHIda.
       So this class simply use the results from FLASHIda log file for assignment. The parsing of FLASHIda log file is done
       in FLASHDeconv tool class.
  */
  class OPENMS_DLLAPI DeconvolvedSpectrum
  {
  public:
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// default constructor
    DeconvolvedSpectrum() = default;

    /**
       @brief Constructor for DeconvolvedSpectrum. Takes the spectrum and scan number calculated from outside
       @param scan_number scan number of the spectrum
  */
    explicit DeconvolvedSpectrum(int scan_number);

    /// default destructor
    ~DeconvolvedSpectrum() = default;

    /// copy constructor
    DeconvolvedSpectrum(const DeconvolvedSpectrum&) = default;

    /// move constructor
    DeconvolvedSpectrum(DeconvolvedSpectrum&& other) noexcept = default;

    /// assignment operator
    DeconvolvedSpectrum& operator=(const DeconvolvedSpectrum& deconvolved_spectrum) = default;


    /// Convert DeconvolvedSpectrum to MSSpectrum (e.g., used to store in mzML format).
    /// @param to_charge the charge of each peak in mzml output.
    /// @param min_ms_level the minimum MS level. If the original spec had an MS level lower than @p min_ms_level the precursor information of the returned spectrum is set to this value.
    /// @param tol the ppm tolerance
    /// @param retain_undeconvolved if set, undeconvolved peaks in the original peaks are output (assuming their abs charge == 1 and m/zs are adjusted with the to_charge parameter)
    MSSpectrum toSpectrum(int to_charge, uint min_ms_level, double tol = 10.0, bool retain_undeconvolved = false);

    /// original spectrum getter
    const MSSpectrum& getOriginalSpectrum() const;

    /// get precursor peak group for MSn (n>1) spectrum. It returns an empty peak group if no peak group is registered (by registerPrecursor)
    PeakGroup& getPrecursorPeakGroup();

    /// precursor charge getter (set in registerPrecursor)
    int getPrecursorCharge() const;

    /// get precursor peak
    const Precursor& getPrecursor() const;

    /// get possible max mass of the deconvolved masses - for MS1, max mass specified by user
    /// for MSn, min value between max mass specified by the user and precursor mass
    /// @param max_mass the max mass specified by the user
    double getCurrentMaxMass(double max_mass) const;

    /// get possible min mass of the deconvolved masses - for MS1, min mass specified by user
    /// for MSn, 50.0
    /// @param min_mass the min mass specified by the user
    double getCurrentMinMass(double min_mass) const;

    /// get possible max charge of the deconvolved masses - for MS1, max charge specified by user
    /// for MSn, min value between max charge specified by the user and precursor charge
    /// @param max_abs_charge the max absolute value of the charge specified by the user
    int getCurrentMaxAbsCharge(int max_abs_charge) const;

    /// get scan number of the original spectrum
    int getScanNumber() const;

    /// get precursor scan number - only if it is registered. Otherwise return 0
    int getPrecursorScanNumber() const;

    /// get activation method
    const Precursor::ActivationMethod& getActivationMethod() const;

    /// set precursor for MSn for n>1
    void setPrecursor(const Precursor& precursor);

    /// set precursor peak intensity
    void setPrecursorIntensity(float i);

    /// set precursor scan number
    void setPrecursorScanNumber(int scan_number);

    /// set activation method
    void setActivationMethod(const Precursor::ActivationMethod& method);

    /// set precursor peakGroup
    void setPrecursorPeakGroup(const PeakGroup& pg);

    /// original spectrum setter
    void setOriginalSpectrum(const MSSpectrum& spec);

    /// set peak groups in this spectrum
    void setPeakGroups(std::vector<PeakGroup>& x);

    /// iterators and vector operators for std::vector<PeakGroup> peak_groups_ in this spectrum
    std::vector<PeakGroup>::const_iterator begin() const noexcept;
    std::vector<PeakGroup>::const_iterator end() const noexcept;

    std::vector<PeakGroup>::iterator begin() noexcept;
    std::vector<PeakGroup>::iterator end() noexcept;

    const PeakGroup& operator[](Size i) const;
    PeakGroup& operator[](Size i);
    void push_back(const PeakGroup& pg);
    Size size() const noexcept;
    void clear();
    void reserve(Size n);
    bool empty() const;

    /// sort by deconvolved monoisotopic masses
    void sort();
    /// sort by Qscore of peakGroups
    void sortByQscore();

  private:
    /// peak groups (deconvolved masses)
    std::vector<PeakGroup> peak_groups_;
    /// the original raw spectrum (not deconvolved)
    MSSpectrum spec_;
    /// precursor peakGroup (or mass)
    PeakGroup precursor_peak_group_;
    /// precursor raw peak (not deconvolved one)
    Precursor precursor_peak_;
    /// activation method for file output
    Precursor::ActivationMethod activation_method_;
    /// scan number and precursor scan number
    int scan_number_ = 0, precursor_scan_number_ = 0;
  };
} // namespace OpenMS