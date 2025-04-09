// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Virginia Rossow, Lenny Kovac, Hendrik Beschorner$
// --------------------------------------------------------------------------

#pragma once

#include <fstream>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MathFunctions.h>

class OPENMS_DLLAPI MQExporterHelper
/**
@brief Helper class for common functions and NON trivial values needed for exporting MaxQuant outputs

@ingroup Metadata
*/
{
public:

  struct MQCommonOutputs
  {
    std::stringstream modifications;
    char acetyl;
    std::stringstream oxidation;
    std::stringstream gene_names;
    std::stringstream protein_names;
    std::stringstream msms_mz;
    std::stringstream mass_error_ppm;
    std::stringstream mass_error_da;
    std::stringstream uncalibrated_mass_error_ppm;
    std::stringstream uncalibrated_mass_error_da;
    std::stringstream uncalibrated_calibrated_mz_ppm;
    std::stringstream uncalibrated_calibrated_mz_mda;
    std::stringstream base_peak_fraction;

    // common columns in msms and evividence exporter
    //file_ << "Sequence" << "\t"; maybe, relativ trivial
    //file_ << "Length" << "\t";
    //file_ << "Modifications" << "\t"; implementieren
    // file_ << "Modified sequence" << "\t"; implementieren
    //file_ << "Acetyl (Protein N-term)" << "\t"; implementieren
    //file_ << "Oxidation (M)" << "\t"; implementieren
    //file_ << "Missed cleavages" << "\t"; trivial
    //file_ << "Proteins" << "\t"; trivial
    //file_ << "Gene Names" << "\t"; // in progress, aber implementieren
    //file_ << "Protein Names" << "\t"; // in progress, aber implementieren
    //file_ << "Type" << "\t"; TODO different type
    //file_ << "Raw file" << "\t"; trivial
    //file_ << "MS/MS m/z" << "\t"; implementieren TODO is m/z in MSMS MS/MS m/z
    //file_ << "Charge" << "\t"; trivial
    //file_ << "m/z" << "\t"; trivial TODO
    //file_ << "Mass" << "\t"; trivial
    //file_ << "Mass Error [ppm]" << "\t"; vielleicht, beim einen halt noch calibrated dabei
    //file_ << "Mass Error [Da]" << "\t"; vielleicht, beim einen halt noch calibrated dabei
    //file_ << "Retention time" << "\t"; trivial
    //file_ << "Fraction of total spectrum" << "\t"; trivial
    //file_ << "Base peak fraction" << "\t"; trvial
    //file_ << "PEP" << "\t"; trivial
    //file_ << "MS/MS Scan Number" << "\t"; trivial
    //file_ << "Score" << "\t"; trivial
    //file_ << "Delta score" << "\t"; trivial
    //file_ << "Reverse" << "\t";
    //file_ << "id" << "\t"; ?
    //file_ << "Protein group IDs" << "\n"; trivial

    explicit MQCommonOutputs(
      const OpenMS::Feature& f,
      const OpenMS::ConsensusMap& cmap,
      const OpenMS::Size c_feature_number,
      const std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>>& UIDs,
      const OpenMS::ProteinIdentification::Mapping& mp_f,
      const OpenMS::MSExperiment& exp,
      const std::map<OpenMS::String,OpenMS::String>& prot_mapper);
  };

  /**
  @brief Extract a gene name from a protein description by looking for the substring 'GN='

  If no such substring exists, an empty string is returned.
*/
static OpenMS::String extractGeneName(const OpenMS::String& prot_description);

  /**
  @brief Returns a unique ID (number) for each distinct protein accession, or creates a new ID by augmenting the given database.

      Obtains a unique, consecutive number for each distinct protein, which can
      be used as a protein ID in the MaxQuant output files (in lack of a proper
      proteingroup ID which maps to proteinGroups.txt)

   @param database A map from accession to ID (which can be augmented by this function)
   @param protein_accession The protein accession which needs translation to an ID

   @return The ID for the @p protein_accession

 */
 static OpenMS::Size proteinGroupID_(std::map<OpenMS::String, OpenMS::Size>& database,
                                     const OpenMS::String& protein_accession);

  /**
    @brief Creates map that has the information which FeatureUID is mapped to which ConsensusFeature in ConsensusMap

    @throw Exception::Precondition if FeatureHandle exists twice in ConsensusMap

    @param cmap ConsensusMap that includes ConsensusFeatures

    @return Returns map, the index is a FeatureID, the value is the index of the ConsensusFeature
    in the vector of ConsensusMap
  */
  static std::map<OpenMS::Size, OpenMS::Size> makeFeatureUIDtoConsensusMapIndex_(const OpenMS::ConsensusMap& cmap);

  /**
    @brief Checks if Feature has valid PeptideIdentifications

      If there are no PeptideIdentifications or the best hit of the Feature cannot be found in corresponding ConsensusFeature,
      the functions returns false to show that something went wrong.

    @param f Feature to extract PeptideIdentifications
    @param c_feature_number Index of corresponding ConsensusFeature in ConsensusMap
    @param UIDs UIDs of all PeptideIdentifications of the ConsensusMap
    @param mp_f Mapping between the FeatureMap and ProteinIdentifications for the UID

    @return Returns true if the PeptideIdentifications exist and are valid
  */
  static bool hasValidPepID_(
    const OpenMS::Feature& f,
    const OpenMS::Size c_feature_number,
    const std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>>& UIDs,
    const OpenMS::ProteinIdentification::Mapping& mp_f);

  /**
    @brief Checks if ConsensusFeature has valid PeptideIdentifications

    If there are no PeptideIdentifications,
    the functions returns false to show that something went wrong.

    @param cf is used to extract PeptideIdentifications

    @return Returns true if the ConsensusFeature has any PepIDs; otherwise false
  */
  static bool hasPeptideIdentifications_(const OpenMS::ConsensusFeature& cf);


  /**
      @brief Checks if file is writable
             (i.e. the path in the ctor was not empty and could be created)

      @return Returns true if file is writable
  */
  static bool isValid(const std::string& filename_);
};