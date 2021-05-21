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

class OPENMS_DLLAPI MQEvidence
{
private:
    std::fstream file_; // Stream where the data is added to create MQEvidence file
    int id_; // number of rows in file to give each row a specific id
    std::map<OpenMS::String, OpenMS::UInt64> protein_id_; // map that maps each Feature to the index of the associated ConsensusFeature in the ConsensusMap

    /**
  @brief Writes the header of MQEvidence.txt in file (Names of columns)
*/
    void exportHeader_();

    /**
      @brief Gets Protein and Returns corresponding number.

        In MQEvidence.txt every protein gets a random but distinct number. The first protein
        a peptide is mapped to get the number one and so on. By this means, it can be seen very easily
        which peptides are mapped to the same protein.

      @param String that is a description of a protein

      @return Returns distinct number for every Protein that is part of the file.
    */
    OpenMS::UInt64 proteinGroupID_(const OpenMS::String &protein);

    /**
      @brief Creates map that has information which FeatureUID is mapped to which ConsensusFeature

      @param cmap ConsensusMap for purpos

      @return Returns map, the index is a FeatureID, the value is the index of the ConsensusFeature
      in the vector of ConsensusMap
    */
    std::map<OpenMS::UInt64, OpenMS::Size> makeFeatureUIDtoConsensusMapIndex_(const OpenMS::ConsensusMap &cmap);

    /**
      @brief Checks if Feature has valid PeptideIdentifications

        If there are no PeptideIdentifications or the best hit of the Feature cannot be found in corresponding ConsensusFeature,
        the functions returns false to show that something went wrong.

      @param cmap, c_feature_number, UIDs and mp_f for comparing Feature and ConsensusFeature, f is used to extract PeptideIdentifications

      @return Returns true if the PeptideIdentifications are valid
    */
    bool hasValidPepID_(
            const OpenMS::Feature &f,
            const OpenMS::Int64 c_feature_number,
            const std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>> &UIDs,
            const OpenMS::ProteinIdentification::Mapping &mp_f);

    /**
      @brief Checks if ConsensusFeature has valid PeptideIdentifications

      If there are no PeptideIdentifications,
      the functions returns false to show that something went wrong.

      @param cf is used to extract PeptideIdentifications

      @return Returns true if the PeptideIdentifications are valid
    */
    bool hasPeptideIdentifications_(const OpenMS::ConsensusFeature &cf);

    /**
      @brief Export one Feature as a row in MQEvidence.txt

        Export one Feature as a row in MQEvidence.txt.
        hasValidPepID and HasPeptideIdentifications are used.
        If there are problems (missing data or confusing data), no output will be generated.

      @param f, cmap and c_number_feature to extract data out of Feature or ConsensusFeature,
      raw_file is specifying the raw_file the feature belongs to,
      c_feature_number, mp_f and UIDs are used in hasValidPepID

    */
    void exportRowFromFeature_(
            const OpenMS::Feature &f,
            const OpenMS::ConsensusMap &cmap,
            const OpenMS::Int64 c_feature_number,
            const OpenMS::String &raw_file,
            const std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>> &UIDs,
            const OpenMS::ProteinIdentification::Mapping &mp_f);

public:
    /**
  @brief Creates MQEvidence object and file MQEvidence.txt in given path, opens f_stream and calls exportHeader

  @param String that is the path where file has to be stored

*/
    explicit MQEvidence(const OpenMS::String &p);

    /**
      @brief Closes f_stream and destruct MQEvidence object
    */
    ~MQEvidence();

    /**
  @brief Checks if file is writable

  @return Returns true if file is writable
*/
    bool isValid();

    /**
      @brief ASHORTDESCRIPTION_HERE

      To export the Features of the FeatureMap the exportFeatureToRow function is used.
      If there are problems (missing data or confusing data) with one feature, it will be skipped.

      @param feature_map and cmap are used to extract data.
    */
    void exportFeatureMap(
            const OpenMS::FeatureMap &feature_map,
            const OpenMS::ConsensusMap &cmap);
};
