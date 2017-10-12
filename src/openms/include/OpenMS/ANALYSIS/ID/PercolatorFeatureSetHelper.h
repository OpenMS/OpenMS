// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_PERCOLATORFEATURESETHELPER_H
#define OPENMS_ANALYSIS_ID_PERCOLATORFEATURESETHELPER_H

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>
#include <limits>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <boost/lexical_cast.hpp>

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
         * @brief concatMULTISEPeptideIds
         * @param all_peptide_ids PeptideIdentification vector to append to
         * @param new_peptide_ids PeptideIdentification vector to be appended
         * @param search_engine search engine to depend on for feature creation
         *
         * Appends a vector of PeptideIdentification to another and prepares Percolator features in MetaInfo (With the respective key "CONCAT:" + search_engine).
         */
        static void concatMULTISEPeptideIds(std::vector<PeptideIdentification>& all_peptide_ids, std::vector<PeptideIdentification>& new_peptide_ids, String search_engine);

        /**
         * @brief mergeMULTISEPeptideIds
         * @param all_peptide_ids PeptideIdentification vector to be merged into
         * @param new_peptide_ids PeptideIdentification vector to merge
         * @param search_engine search engine to create features from their scores
         *
         * Merges a vector of PeptideIdentification into another and prepares the merged MetaInfo and scores for collection in addMULTISEFeatures for feature registration.
         */
        static void mergeMULTISEPeptideIds(std::vector<PeptideIdentification>& all_peptide_ids, std::vector<PeptideIdentification>& new_peptide_ids, String search_engine);

        /**
         * @brief mergeMULTISEProteinIds
         * @param all_protein_ids ProteinIdentification vector to be merged into
         * @param new_protein_ids ProteinIdentification vector to merge
         *
         * Concatenates SearchParameter of multiple search engine runs and merges PeptideEvidences, collects used search engines in MetaInfo for collection in addMULTISEFeatures for feature registration.
         */
        static void mergeMULTISEProteinIds(std::vector<ProteinIdentification>& all_protein_ids, std::vector<ProteinIdentification>& new_protein_ids);
        

        /**
         * @brief addMSGFFeatures
         * @param peptide_ids PeptideIdentification vector to create Percolator features in
         * @param feature_set register of added features
         *
         * Creates and adds MSGF+ specific Percolator features and registers them in feature_set. MSGF+ should be run with the addFeatures flag enabled.
         */
        static void addMSGFFeatures(std::vector<PeptideIdentification>& peptide_ids, StringList& feature_set);

        /**
         * @brief addXTANDEMFeatures
         * @param peptide_ids PeptideIdentification vector to create Percolator features in
         * @param feature_set register of added features
         *
         * Creates and adds X!Tandem specific Percolator features and registers them in feature_set
         */
        static void addXTANDEMFeatures(std::vector<PeptideIdentification>& peptide_ids, StringList& feature_set);

        /**
         * @brief addCOMETFeatures
         * @param peptide_ids PeptideIdentification vector to create Percolator features in
         * @param feature_set register of added features
         *
         * Creates and adds Comet specific Percolator features and registers them in feature_set
         */
        static void addCOMETFeatures(std::vector<PeptideIdentification>& peptide_ids, StringList& feature_set);

        /**
         * @brief addMASCOTFeatures
         * @param peptide_ids PeptideIdentification vector to create Percolator features in
         * @param feature_set register of added features
         *
         * Creates and adds Mascot specific Percolator features and registers them in feature_set
         */
        static void addMASCOTFeatures(std::vector<PeptideIdentification>& peptide_ids, StringList& feature_set);

        /**
         * @brief addMULTISEFeatures
         * @param peptide_ids PeptideIdentification vector to create Percolator features in
         * @param search_engines_used the list of search engines to be considered
         * @param feature_set register of added features
         * @param complete_only will only add features for PeptideIdentifications where all given search engines identified something
         * @param limits_imputation
         *
         * Adds multiple search engine specific Percolator features and registers them in feature_set
         */
        static void addMULTISEFeatures(std::vector<PeptideIdentification>& peptide_ids, StringList& search_engines_used, StringList& feature_set, bool complete_only = true, bool limits_imputation = false);

        /**
         * @brief addCONCATSEFeatures
         * @param peptide_id_list PeptideIdentification vector to create Percolator features in
         * @param search_engines_used the list of search engines to be considered
         * @param feature_set register of added features
         *
         * Adds multiple search engine specific Percolator features and registers them in feature_set
        */
        static void addCONCATSEFeatures(std::vector<PeptideIdentification>& peptide_id_list, StringList& search_engines_used, StringList& feature_set);

        /**
         * @brief checkExtraFeatures
         * @param psms the vector of PeptideHit to be checked
         * @param extra_features the list of requested extra features
         *
         * checks and removes requested extra Percolator features that are actually unavailable (to compute)
         */
        static void checkExtraFeatures(const std::vector<PeptideHit> &psms, StringList& extra_features);
        

    protected:
        /// Rescales the fragment features to penalize features calculated by few ions, adapted from MSGFtoPercolator
        static double rescaleFragmentFeature_(double featureValue, int NumMatchedMainIons);

        /// helper function for assigning the frequently occurring feature delta score
        static void assignDeltaScore_(std::vector<PeptideHit>& hits, String score_ref, String output_ref);

        /// gets the scan identifier to merge by
        static String getScanMergeKey_(std::vector<PeptideIdentification>::iterator it, std::vector<PeptideIdentification>::iterator start);

        /// For accession dependent sorting of ProteinHits
        struct lq_ProteinHit
        {
          inline bool operator() (const ProteinHit& h1, const ProteinHit& h2)
          {
            return (h1.getAccession() < h2.getAccession());
          }
        };

        /// For accession dependent sorting of PeptideEvidences
        struct lq_PeptideEvidence
        {
          inline bool operator() (const PeptideEvidence& h1, const PeptideEvidence& h2)
          {
            return (h1.getProteinAccession() < h2.getProteinAccession());
          }
        };

    };

} //namespace OpenMS

#endif //OPENMS_ANALYSIS_ID_PERCOLATORFEATURESETHELPER_H

