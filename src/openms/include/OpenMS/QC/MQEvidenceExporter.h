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
// $Maintainer: Chris Bielow$
// $Authors: Valentin Noske, Vincent Musch$
// --------------------------------------------------------------------------

#pragma once

#include <fstream>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>

class OPENMS_DLLAPI MQEvidence{
private:
    std::fstream file_;
    int id_;
    std::map<OpenMS::String,OpenMS::UInt64> protein_id_;

    /**
  @brief ASHORTDESCRIPTION_HERE

    Writes the header of MQEvidence.txt in file (Names of columns)
*/
    void exportHeader();

    /**
      @brief ASHORTDESCRIPTION_HERE

        Gets Protein and Returns corresponding number.

      LONGER_DESCRIPTION HERE

        In MQEvidence.txt every protein gets a random but distinct number. The first protein
        a peptide is mapped to get the number one and so on. By this means, it can be seen very easily
        which peptides are mapped to the same protein.

      @param String that is a description of a protein

      @return Returns distinct number for every Protein that is part of the file.
    */
    OpenMS::UInt64 proteinGroupID(const OpenMS::String &protein);

    /**
      @brief ASHORTDESCRIPTION_HERE

        Creates map that has information which FeatureID is mapped to which ConsensusFeature

      @param cmap ConsensusMap for purpos

      @return Returns map, the index is a FeatureID, the value is the index of the ConsensusFeature
      in the vector of ConsensusMap
    */
    std::map<OpenMS::UInt64, OpenMS::Size> mapFeatureIDtoConsensusID(const OpenMS::ConsensusMap & cmap);

    /**
      @brief ASHORTDESCRIPTION_HERE

        Gives a pointer best peptideHit referring to the score of a vector PeptideIdentifications
        and the index of PeptideIdentification the best PeptideHit belongs to

      @param vector of PeptideIdentifications

      @return Returns pair of the index of the PeptideIdentification and the pointer to the best hit
    */
    std::pair<OpenMS::Size , const OpenMS::PeptideHit *> getBestPeptideHit(const std::vector<OpenMS::PeptideIdentification>& pep_ids);

    /**
      @brief ASHORTDESCRIPTION_HERE

        Set the given Pointer  to the best hit of the PeptideIdentifications of the Feature

      LONGER_DESCRIPTION HERE

        Set the given Pointer to the best hit of the PeptideIdentifications of the Feature.
        If there are no PeptideIdentifications, the best hit of the Feature cannot be found in corresponding ConsensusFeature
        or the best hit does not have a sequence, the functions returns false to show that something went wrong.
        getBestPeptideHit is used.

      @param cmap, c_feature_number, UIDs and mp_f for comparing Feature and ConsensusFeature, Feature to get information about the PeptideIdenfications.
      best for saving the information about the best PeptideHit.

      @return Returns true if the functions succeeds.
    */
    bool exportFeaturePepID(
            const OpenMS::ConsensusMap &cmap,
            const OpenMS::Int64 c_feature_number,
            const std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>> &UIDs,
            const OpenMS::ProteinIdentification::Mapping &mp_f,
            const OpenMS::Feature & f,
            const OpenMS::PeptideHit* & best);

    /**
      @brief ASHORTDESCRIPTION_HERE

        Set the given Pointer to the best hit of the PeptideIdentifications of the ConsensusFeature

      LONGER_DESCRIPTION HERE

      Set the given Pointer to the best hit of the PeptideIdentifications of the ConsensusFeature
      If there are no PeptideIdentifications, the best hit of the Feature cannot be found in corresponding ConsensusFeature
      or the best hit does not have a sequence, the functions returns false to show that something went wrong.
      getBestPeptideHit is used.

      @param cmap and c_feature_number are used to identify the Consensusfeature and extract the information about the PeptideIdentifications
      best for saving the information about the best PeptideHit.

      @return Returns true if the functions succeeds.
    */
    bool exportConsensusFeaturePepID(
            const OpenMS::ConsensusMap &cmap,
            const OpenMS::Int64 c_feature_number,
            const OpenMS::PeptideHit* & best);

    /**
      @brief Export one Feature as a row in MQEvidence.txt

        Export one Feature as a row in MQEvidence.txt.
        exportFeaturePepID and exportConsensusFeaturePepID are used.
        If there are problems (missing data or confusing data), no output will be generated.

      @param cmap ConsensusMap for purpose XYZ...

    */
    void exportRowFromFeature(
            const OpenMS::Feature &f,
            const OpenMS::ConsensusMap &cmap,
            const OpenMS::Int64 c_feature_number,
            const OpenMS::String &raw_file,
            const std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>> &UIDs,
            const OpenMS::ProteinIdentification::Mapping &mp_f);

public:
        /**
      @brief ASHORTDESCRIPTION_HERE

         Creates file MQEvidence.txt in given path, opens f_stream and calls exportHeader

      @param String that is the path where file has to be stored

    */
    explicit MQEvidence(const OpenMS::String & p);

    /**
      @brief ASHORTDESCRIPTION_HERE.

        Closes f_stream
    */
    ~MQEvidence();

        /**
      @brief ASHORTDESCRIPTION_HERE.

         Checks if file is writable

      @return Returns true if file is writable
    */
    bool isValid();

    /**
      @brief ASHORTDESCRIPTION_HERE

      LONGER_DESCRIPTION HERE

      @param cmap ConsensusMap for purpose XYZ...

      @return Returns true if ...
    */
    void exportFeatureMapTotxt(
            const OpenMS::FeatureMap & feature_map,
            const OpenMS::ConsensusMap& cmap);
};
