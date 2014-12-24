// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_QUANTITATION_PROTEININFERENCE_H
#define OPENMS_ANALYSIS_QUANTITATION_PROTEININFERENCE_H

#include <vector>

#include <OpenMS/KERNEL/Peak2D.h>

namespace OpenMS
{

  class PeptideHit;
  class ConsensusMap;

  /**
    @brief [experimental class] given a peptide quantitation, infer corresponding protein quantities

    Infers protein ratios from peptide ratios (currently using unique peptides only).
    Use the IDMapper class to add protein and peptide information to a
    quantitative ConsensusMap prior to this step.
  */
  class OPENMS_DLLAPI ProteinInference
  {

public:

    typedef Peak2D::IntensityType IntensityType;

    /// Constructor
    ProteinInference();

    /// copy constructor
    ProteinInference(const ProteinInference & cp);

    /// assignment operator
    ProteinInference & operator=(const ProteinInference & rhs);

    /**
        @brief given a peptide quantitation, infer corresponding protein quantities

        Infers protein ratios from peptide ratios (currently using unique peptides only).
        Use the IDMapper class to add protein and peptide information to a
        quantitative ConsensusMap prior to this step.

        @param consensus_map Peptide quantitation with ProteinIdentifications attached, where
                     Protein quantitation will be attached
        @param reference_map Index of (iTRAQ) reference channel within the consensus map

        @throws Exception::MissingInformation if Protein/PeptideIdentifications are missing
    */
    void infer(ConsensusMap & consensus_map, const UInt reference_map);


protected:

    void infer_(ConsensusMap & consensus_map,
                const size_t protein_idenfication_index,
                const UInt reference_map);

    bool sortByUnique_(std::vector<PeptideHit> & peptide_hits_local, const bool is_higher_score_better);

  };   // !class

} // !namespace

#endif // OPENMS_ANALYSIS_QUANTITATION_PROTEININFERENCE_H
