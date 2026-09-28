// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------

#pragma once

#include <vector>
#include <cmath>
#include <string>
#include <map>
// #include <algorithm>
#include <limits>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

namespace OpenMS
{
    /**
        @brief Percolator feature set and integration helper

        This class contains functions to handle (compute, aggregate, integrate)
        Percolator features. This includes the calculation or extraction of
        Percolator features depending on the search engine(s) for later use with
        PercolatorAdapter. It also includes handling the reintegration of the
        percolator result into the set of Identifications.
    */

    class OPENMS_DLLAPI PercolatorFeatureSetHelper
    {

    public:
        /**
          @brief mergeMULTISEProteinIds
          @param[in] all_protein_ids ProteinIdentification vector to be merged into
          @param[in] new_protein_ids ProteinIdentification vector to merge

          Merges a ProteinIdentification run into another: takes over the search parameters of the
          first run, sets the search engine to "multiple" once runs from differing search engines are
          combined, and unions the ProteinHits by accession.
         */
        static void mergeMULTISEProteinIds(std::vector<ProteinIdentification>& all_protein_ids, std::vector<ProteinIdentification>& new_protein_ids);


        /**
          @brief addMSGFFeatures
          @param[in] peptide_ids PeptideIdentification vector to create Percolator features in
          @param[in] feature_set register of added features
         
          Creates and adds MSGF+ specific Percolator features and registers them in feature_set. MSGF+ should be run with the addFeatures flag enabled.
         */
        static void addMSGFFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set);

        /**
          @brief addXTANDEMFeatures
          @param[in] peptide_ids PeptideIdentification vector to create Percolator features in
          @param[in] feature_set register of added features
         
          Creates and adds X!Tandem specific Percolator features and registers them in feature_set
         */
        static void addXTANDEMFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set);

        /**
          @brief addCOMETFeatures
          @param[in] peptide_ids PeptideIdentification vector to create Percolator features in
          @param[in] feature_set register of added features
         
          Creates and adds Comet specific Percolator features and registers them in feature_set
         */
        static void addCOMETFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set);

        /**
          @brief addANDESFeatures
          @param[in] peptide_ids PeptideIdentification vector to create Percolator features in
          @param[in] feature_set register of added features

          Collects and registers the andes search engine specific Percolator features.
          andes annotates every PSM with its own numeric features as MetaValues named
          with the "andes:" prefix (e.g. "andes:RankScore", "andes:RawScore", ...). Rather
          than hard-coding the (currently 47) feature names here, all numeric "andes:"-prefixed
          MetaValues found on the hits are collected, keeping the andes side the single source
          of truth for its feature set.
         */
        static void addANDESFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set);

        /**
          @brief addMASCOTFeatures
          @param[in] peptide_ids PeptideIdentification vector to create Percolator features in
          @param[in] feature_set register of added features
         
          Creates and adds Mascot specific Percolator features and registers them in feature_set
         */
        static void addMASCOTFeatures(PeptideIdentificationList& peptide_ids, StringList& feature_set);

        /**
          @brief checkExtraFeatures
          @param[in] psms the vector of PeptideHit to be checked
          @param[in] extra_features the list of requested extra features
         
          checks and removes requested extra Percolator features that are actually unavailable (to compute)
         */
        static void checkExtraFeatures(const std::vector<PeptideHit> &psms, StringList& extra_features);

        /**
         * @brief addMSFraggerFeatures
         * @param[in] extra_features register of added features
         *
         * Registers the MSFragger specific Percolator features in extra_features.
         */
        static void addMSFRAGGERFeatures(StringList& extra_features);
        

    protected:
        /// Rescales the fragment features to penalize features calculated by few ions, adapted from MSGFtoPercolator
        static double rescaleFragmentFeature_(double featureValue, int NumMatchedMainIons);

        /// helper function for assigning the frequently occurring feature delta score
        static void assignDeltaScore_(std::vector<PeptideHit>& hits, const std::string& score_ref, const std::string& output_ref);

        /// For accession dependent sorting of ProteinHits
        struct lq_ProteinHit
        {
          inline bool operator() (const ProteinHit& h1, const ProteinHit& h2)
          {
            return (h1.getAccession() < h2.getAccession());
          }
        };

    };

} //namespace OpenMS


