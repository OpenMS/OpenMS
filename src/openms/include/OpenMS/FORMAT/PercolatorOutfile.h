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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>

#include <vector>

namespace OpenMS
{

  /**
     @brief Class for reading Percolator tab-delimited output files.

     For PSM-level output, the file extension should be ".psms".
  */
  class OPENMS_DLLAPI PercolatorOutfile
  {

  public:

    /// Types of Percolator scores
    enum ScoreType { QVALUE, POSTERRPROB, SCORE, SIZE_OF_SCORETYPE };

    /// Names of Percolator scores (to match ScoreType)
    static const std::string score_type_names[SIZE_OF_SCORETYPE];

    /// Return a score type given its name
    static enum ScoreType getScoreType(String score_type_name);

    /// Constructor
    PercolatorOutfile();

    /// Loads a Percolator output file
    void load(const String& filename, ProteinIdentification& proteins,
              std::vector<PeptideIdentification>& peptides,
              SpectrumMetaDataLookup& lookup,
              enum ScoreType output_score = QVALUE);

  private:
    /// Converts the peptide string to an 'AASequence' instance
    void getPeptideSequence_(String peptide, AASequence& seq) const;

    /// Resolve cases where N-terminal modifications may be misassigned to the first residue (for X! Tandem results)
    void resolveMisassignedNTermMods_(String& peptide) const;
  };

} // namespace OpenMS

