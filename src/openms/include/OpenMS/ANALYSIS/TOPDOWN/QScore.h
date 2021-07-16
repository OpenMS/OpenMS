//--------------------------------------------------------------------------
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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>

namespace OpenMS
{
  class PeakGroup;

  /**
@brief   QScore : quality score for precursors. This class is being updated. For now, simply it calculate the QScore using a fixed weight vector. But afterwards, the training
   part for the QScore should be added in here.
@ingroup Topdown
*/

  class OPENMS_DLLAPI QScore
  {
  public:
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// get QScore for a peak group of specific abs_charge
    static double getQScore(const PeakGroup *pg, const int abs_charge);

    /// function to generate attribute tsv file for weka interface (for now)
    static void writeAttTsv(const int scan_number,
                            const String &acc,
                            const int proID,
                            const double rt,
                            const int pscan,
                            const double pmass,
                            const double pmz,
                            const double fintensity,
                            PeakGroup &pg,
                            const int fr,
                            const int lr,
                            const int charge,
                            const double precursor_intensity,
                            const std::vector<double> ptm_mass,
        //const std::vector<int> ptm_start,
        //const std::vector<int> ptm_end,
                            const bool is_identified,
                            const double e_value,
                            const double q_value,
                            const FLASHDeconvHelperStructs::PrecalculatedAveragine &avg,
                            std::fstream &f,
                            bool write_detail = false);

    /// write header for attirbute tsv file
    static void writeAttHeader(std::fstream &f, bool write_detail = false);

  private:
    /// convert a peak group to a feature vector for QScore calculation
    static std::vector<double> toFeatureVector_(const PeakGroup *pg, const int abs_charge);
  };
}
