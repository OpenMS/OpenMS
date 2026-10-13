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
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <optional>

class MQExporterHelper
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

    /// The columns of feature @p f with its best hit @p best_hit (see bestIdentification())
    explicit MQCommonOutputs(
      const OpenMS::Feature& f,
      const OpenMS::PeptideHit& best_hit,
      const OpenMS::MSExperiment& exp,
      const std::map<std::string,std::string>& prot_mapper);
  };

  /// An identification that a feature links, as the MaxQuant exporters read it
  struct LinkedIdentification
  {
    std::string file;                            ///< path of the file the identification was made from (its source)
    OpenMS::IdentificationData::Observation observation; ///< spectrum reference, RT, m/z and meta values
    std::optional<OpenMS::PeptideHit> best;      ///< the match with the best primary score; none without scored matches
  };

  /**
    @brief The identifications that @p feature links in @p data, best first

    Sorted by the primary score of their best match, in the direction of the score (of equal scores, in the order of
    BaseFeature::getLinkedIdentifications()); identifications without a scored match go last.
  */
  static std::vector<LinkedIdentification> linkedIdentifications(const OpenMS::BaseFeature& feature, const OpenMS::IdentificationData& data);

  /**
    @brief The identification to export for a feature and how it was found (@p type), or nullptr

    The best identification of the feature (@p f_ids, see linkedIdentifications()) if it has a hit and the consensus
    feature links an identification of the same spectrum (@p cf_ids), as "MULTI-MSMS"; else the best identification of
    the consensus feature if it has a hit, as "MULTI-MATCH". The result has a hit (LinkedIdentification::best).
  */
  static const LinkedIdentification* bestIdentification(const std::vector<LinkedIdentification>& f_ids, const std::vector<LinkedIdentification>& cf_ids, std::string& type);

  /**
  @brief Extract a gene name from a protein description by looking for the substring 'GN='

  If no such substring exists, an empty string is returned.
*/
static std::string extractGeneName(const std::string& prot_description);

  /**
  @brief Returns a unique ID (number) for each distinct protein accession, or creates a new ID by augmenting the given database.

      Obtains a unique, consecutive number for each distinct protein, which can
      be used as a protein ID in the MaxQuant output files (in lack of a proper
      proteingroup ID which maps to proteinGroups.txt)

   @param[in,out] database A map from accession to ID (which can be augmented by this function)
   @param[in] protein_accession The protein accession which needs translation to an ID

   @return The ID for the @p protein_accession

 */
 static OpenMS::Size proteinGroupID_(std::map<std::string, OpenMS::Size>& database,
                                     const std::string& protein_accession);

  /**
    @brief Creates map that has the information which FeatureUID is mapped to which ConsensusFeature in ConsensusMap

    @throw Exception::Precondition if FeatureHandle exists twice in ConsensusMap

    @param[in] cmap ConsensusMap that includes ConsensusFeatures

    @return Returns map, the index is a FeatureID, the value is the index of the ConsensusFeature
    in the vector of ConsensusMap
  */
  static std::map<OpenMS::Size, OpenMS::Size> makeFeatureUIDtoConsensusMapIndex_(const OpenMS::ConsensusMap& cmap);

  /**
      @brief Checks if file is writable
             (i.e. the path in the ctor was not empty and could be created)

      @return Returns true if file is writable
  */
  static bool isValid(const std::string& filename_);
};