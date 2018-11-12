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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>

#include <unordered_set>

namespace OpenMS
{
  class OPENMS_DLLAPI IDMergerAlgorithm:
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    explicit IDMergerAlgorithm (const String& runIdentifier = "merged");

    /// Insert (=move and clear) a run with its peptide IDs into the internal merged data structures,
    /// based on the initial mapping from fileorigins to new run
    void insertRun(std::vector<ProteinIdentification>& prots,
        std::vector<PeptideIdentification>& peps);

    /// Return the merged results and reset/clear all internal data
    void returnResultsAndClear(ProteinIdentification& prots,
                   std::vector<PeptideIdentification>& peps);

  private:
    String getNewIdentifier_() const;
    void copySearchParams_(ProteinIdentification& from, ProteinIdentification& to);

    /// Checks consistency of search engines and settings across runs before merging.
    /// @param protRuns The runs to check (first = implicit reference)
    /// @param experiment_type allow some mismatches in case of other experiment types (e.g. SILAC)
    /// @return all same? TODO: a merged RunDescription about what to put in the new runs (e.g. for SILAC)
    /// @throws BaseException for disagreeing settings
    bool checkOldRunConsistency_(
        const std::vector<ProteinIdentification>& protRuns,
        const String& experiment_type) const;

    /// Same as above, if you want to use a specific reference
    /// @param protRuns The runs to check (first = reference)
    /// @param ref A possibly external protein run reference
    /// @param experiment_type allow some mismatches in case of other experiment types (e.g. SILAC)
    /// @return all same? TODO: a merged RunDescription about what to put in the new runs (e.g. for SILAC)
    /// @throws BaseException for disagreeing settings
    bool checkOldRunConsistency_(
        const std::vector<ProteinIdentification>& protRuns,
        const ProteinIdentification& ref,
        const String& experiment_type) const;


    void movePepIDsAndRefProteinsToResult_(
        std::vector<PeptideIdentification>& pepIDs,
        std::vector<ProteinIdentification>& oldProtRuns
    );

    ProteinIdentification protResult;
    std::vector<PeptideIdentification> pepResult;
    std::unordered_set<std::string> proteinsCollected;
    bool filled = false;
    std::map<String, Size> fileOriginToIdx;
    String id;
  };
} // namespace OpenMS
