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
    /**
@brief Builds an MaxQuant Evidence.txt

  This class is closely related to QualityControl, it creates an evidence file similar
  to an MaxQuant evidence.txt. But not all columns of a MaxQuant file get exported.
  By the construction of an object, the column names of the evidence values are added to the file.
  For the construction a valid path is needed (check out constructor) where the evidence.txt can be stored.
  To fill the output file with data from the MS/MS run use the exportFeatureMap function,
  it needs a FeatureMap and the matching ConsensusMap as an input.
  To check if the created file is writable use the function isValid.

    @ingroup Metadata
*/
{
private:
    std::fstream file_; // Stream where the data is added to create MQEvidence file
    int id_; // number of rows in file to give each row a specific id
    std::map<OpenMS::String, OpenMS::Size> protein_id_; // map that maps each Feature to the index of the associated ConsensusFeature in the ConsensusMap

    /**
  @brief Writes the header of MQEvidence.txt in file (Names of columns)
*/
    void exportHeader_();

    /**
      @brief returns the MaxQuant unique evidence number of a protein accession

        In MQEvidence.txt every protein gets a random but distinct number. The first protein
        a peptide is mapped to get the number one and so on. By this means, it can be seen very easily
        which peptides are mapped to the same protein.

      @param String that is a description of a protein

      @return Returns distinct number for every Protein that is part of the file.
    */
    OpenMS::Size proteinGroupID_(const OpenMS::String &protein);

    /**
      @brief Creates map that has the information which FeatureUID is mapped to which ConsensusFeature in ConsensusMap

      @param cmap ConsensusMap that includes ConsensusFeatures

      @return Returns map, the index is a FeatureID, the value is the index of the ConsensusFeature
      in the vector of ConsensusMap
    */
    std::map<OpenMS::Size, OpenMS::Size> makeFeatureUIDtoConsensusMapIndex_(const OpenMS::ConsensusMap &cmap);

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
  @brief Creates MQEvidence object and file MQEvidence.txt in given path

        If the path for the constructor is empty (path not valid), no file is created.
        If the creation of the fstream object is successful a constant header is added to the file

  @param String that is the path where file has to be stored

*/
    explicit MQEvidence(const OpenMS::String &path);

    /**
      @brief Closes f_stream and destruct MQEvidence object
    */
    ~MQEvidence();

    /**
        @brief Checks if file is writable

                can be called anytime after MQEvidence is constructed

        @return Returns true if file is writable
    */
    bool isValid();

    /**
      @brief Exports a FeatureMap to the evidence file

      export relevant evidence data from the featureMap and the ConsensusMAp to the evidence file.
      Features without a PeptideIdentification or with a PeptideIdentification who doesnt match with
      the PeptideIdentification of the ConsensusFeature from the ConsensusMap get skipped.
      If the file not Valid the tool throws an exception.

      @param feature_map and the cmap, which contains the fmap, are used to extract data.
    */
    void exportFeatureMap(
            const OpenMS::FeatureMap &feature_map,
            const OpenMS::ConsensusMap &cmap);
};