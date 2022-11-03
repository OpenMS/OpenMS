// --------------------------------------------------------------------------
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
// $Maintainer: Chris Bielow$
// $Authors: Valentin Noske, Vincent Musch$
// --------------------------------------------------------------------------

#pragma once

#include <fstream>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/QC/MQExporterHelper.h>

#include <map>

class OPENMS_DLLAPI MQEvidence
/**
@brief Builds a MaxQuant Evidence.txt

This class is closely related to QualityControl, it creates an evidence.txt similar
to a MaxQuant evidence.txt. But not all columns of a MaxQuant file get exported.
By the construction of an object, the column names of the evidence values are added to the evidence.txt.
For the construction a valid path is needed (check out constructor) where the evidence.txt can be stored.
To fill the output evidence.txt with data from the MS/MS run use the exportFeatureMap function,
it needs a FeatureMap and the matching ConsensusMap as an input.
To check if the created evidence.txt is writable use the function isValid.

@ingroup Metadata
*/
{
private:
  std::fstream file_; ///< Stream where the data is added to create evidence.txt
  OpenMS::Size id_ = 0; ///< number of rows in evidence.txt to give each row a specific id
  OpenMS::String filename_; ///< path and name of the evidence.txt

  /**
    @brief Writes the header of evidence.txt (Names of columns)
  */
  void exportHeader_();

  /**
    @brief Export one Feature as a row in MQEvidence.txt

      If the feature has no PepID's or the corresponding CF has no PepIDs,
      no row will be exported

    @param f Feature to extract evidence data
    @param cmap ConsensusMap to extract evidence data if Feature has no valid PeptideIdentifications
    @param c_feature_number Index of corresponding ConsensusFeature in ConsensusMap
    @param raw_file is specifying the raw_file the feature belongs to
    @param UIDs UIDs of all PeptideIdentifications of the ConsensusMap
    @param mp_f Mapping between the FeatureMap and ProteinIdentifications for the UID
           from PeptideIdenfitication::buildUIDfromAllPepIds
    @param exp MS Experiment holds evidence data to extract
    @param prot_map Mapping a protein_accession to its description(proteinname, genename...)
  */
  void exportRowFromFeature_(
    const OpenMS::Feature& f,
    const OpenMS::ConsensusMap& cmap,
    const OpenMS::Size c_feature_number,
    const OpenMS::String& raw_file,
    const std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>>& UIDs,
    const OpenMS::ProteinIdentification::Mapping& mp_f,
    const OpenMS::MSExperiment& exp= {},
    const std::map<OpenMS::String,OpenMS::String>& prot_map = {});

public:
  /**
    @brief Creates MQEvidence object and evidence.txt in given path

          If the path for the constructor is empty (path not valid), no evidence.txt is created.
          If the creation of the fstream object is successful a constant header is added to the evidence.txt
          If the path does not exist, it will be created

    @throw Exception::FileNotWritable if evidence.txt could not be created

    @param path that is the path where evidence.txt has to be stored

  */
  explicit MQEvidence(const OpenMS::String& path);

  /**
    @brief Closes f_stream
  */
  ~MQEvidence();

  /**
    @brief Exports a FeatureMap to the evidence.txt

      Exports one row per feature from the FeatureMap to the evidence.txt file.

    @throw Exception::FileNotWritable if evidence.txt is not writable
    @throw Exception::MissingInformation if Feature_map has no corresponding ConsensusFeature

    @param feature_map which contains Features to extract evidence data
    @param cmap ConsensusMap to extract evidence data if Feature has no valid PeptideIdentifications
    @param exp MS Experiment holds evidence data to extract
    @param prot_map Mapping a protein_accession to its description(proteinname, genename...)
  */
  void exportFeatureMap(
    const OpenMS::FeatureMap& feature_map,
    const OpenMS::ConsensusMap& cmap,
    const OpenMS::MSExperiment& exp= {},
    const std::map<OpenMS::String,OpenMS::String>& prot_map = {});
};