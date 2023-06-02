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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/KERNEL/Peak1D.h>


namespace OpenMS
{
  class PeakGroup;

  /**
@brief   Qvalue : contains functions to calculate Qvalues from deconvolution quality score
@ingroup Topdown
*/

  class OPENMS_DLLAPI Qvalue
  {
  public:
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// Calculate and perform a batch update of peak group qvalues using Qscores of target and dummy peak groups in deconvolved spectra, when FDR report is necessary.
    /// @param deconvolved_spectra target deconvolved spectra
    /// @param deconvolved_decoy_spectra decoy deconvolved spectra
    void static updatePeakGroupQvalues(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<DeconvolvedSpectrum>& deconvolved_decoy_spectra);

  private:
    /// get a bin number given qvalue. qvalue is calculated per bin (bin from 0 to 1).
    static uint getBinNumber(float qscore, uint total_bin_number);
    /// get the qvalue corresponding to a bin number
    static float getBinValue(uint bin_number, uint total_bin_number);
    /// get the Qscore distribution
    static std::vector<float> getDistribution(const std::vector<float>& qscores, uint bin_number);
    /// get the weights of different dummy types.
    static std::vector<float> getDistributionWeights(const std::vector<float>& mixed_dist, const std::vector<std::vector<float>>& comp_dists, uint num_iterations = 100);
  };
} // namespace OpenMS