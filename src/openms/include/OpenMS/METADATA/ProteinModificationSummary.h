// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#pragma once

#include <map>

#include <OpenMS/CHEMISTRY/ResidueModification.h>

namespace OpenMS
{
  /**
   * @brief Maps a protein position to all observed modifications and associated statistics 
   * 
   * For example, to store that position 10 maps to Oxidation (M) which was observed in 123 PSMs.
   * 
   */
  struct OPENMS_DLLAPI ProteinModificationSummary
  {
    /// basic modification statistic
    struct OPENMS_DLLAPI Statistics
    {
      bool operator==(const Statistics& rhs) const;
      size_t count = -1;  ///< total number of PSMs supporting the modification at this position
      double frequency = -1.0; ///< PSMs with modification / total number of PSMs
      double FLR = -1.0; ///< false localization rate
      double probability = -1.0; ///< (localization) probability
    };

    /// comparison operator
    bool operator==(const ProteinModificationSummary& rhs) const;

    using ModificationsToStatistics = std::map<ResidueModification, Statistics>;
    using AALevelModificationSummary = std::map<size_t, ModificationsToStatistics>;

    /// position -> modification -> statistic (counts, etc.)
    AALevelModificationSummary AALevelSummary;
  };
}

