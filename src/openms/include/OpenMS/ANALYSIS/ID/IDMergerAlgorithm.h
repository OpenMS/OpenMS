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

  //TODO add params for checking consistency
  //TODO add another subclass that does score-aware merging?
  class OPENMS_DLLAPI IDMergerAlgorithm:
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    explicit IDMergerAlgorithm (const String& runIdentifier = "merged");

    /// Insert (=move and clear) a run with its peptide IDs into the internal merged data structures,
    /// based on the initial mapping from fileorigins to new run
    void insertRuns(std::vector<ProteinIdentification>&& prots,
                    std::vector<PeptideIdentification>&& peps);
    void insertRuns(const std::vector<ProteinIdentification>& prots,
                    const std::vector<PeptideIdentification>& peps);
    //TODO add methods to just insert prots or just peps. Especially makes sense if you do re-indexing anyway,
    // then you do not need the proteins. But then we need origin information. Either externally in form of a
    // String or StringList (like the one from ProteinID.getPrimaryMSRunPath). Or by having the file annotated
    // at the PeptideID (with getBasename maybe?)
    // Current solution would be to clear the ProteinIdentification if you do not need the proteins and add all the
    // necessary information about origin(s) to this ProteinIdentification.

    /// Return the merged results and reset/clear all internal data
    void returnResultsAndClear(ProteinIdentification& prots,
                   std::vector<PeptideIdentification>& peps);

  private:
    String getNewIdentifier_() const;
    static void copySearchParams_(const ProteinIdentification& from, ProteinIdentification& to);

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


    /*void movePepIDsAndRefProteinsToResult_(
        std::vector<PeptideIdentification>&& pepIDs,
        std::vector<ProteinIdentification>&& oldProtRuns
    );*/

    void insertProteinIDs_(
        std::vector<ProteinIdentification>&& oldProtRuns
    );

    void updateAndMovePepIDs_(
        std::vector<PeptideIdentification>&& pepIDs,
        const std::map<String, Size>& runIDToRunIdx,
        const std::vector<StringList>& originFiles,
        bool annotate_origin
    );

    void movePepIDsAndRefProteinsToResultFaster_(
        std::vector<PeptideIdentification>&& pepIDs,
        std::vector<ProteinIdentification>&& oldProtRuns
    );

    ProteinIdentification protResult;
    std::vector<PeptideIdentification> pepResult;

    static size_t accessionHash(const ProteinHit& p){
      return std::hash<String>()(p.getAccession());
    }
    static bool accessionEqual(const ProteinHit& p1, const ProteinHit& p2){
      return p1.getAccession() == p2.getAccession();
    }
    using hash_type = std::size_t (*)(const ProteinHit&);
    using equal_type = bool (*)(const ProteinHit&, const ProteinHit&);
    std::unordered_set<ProteinHit, hash_type, equal_type> proteinsCollectedHits;

    bool filled = false;
    std::map<String, Size> fileOriginToIdx;
    String id;
  };
} // namespace OpenMS
