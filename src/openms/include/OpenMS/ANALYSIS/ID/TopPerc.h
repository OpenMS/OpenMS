// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#ifndef OPENMS_ANALYSIS_ID_TOPPERC_H
#define OPENMS_ANALYSIS_ID_TOPPERC_H

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <map>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>

namespace OpenMS
{
    class OPENMS_DLLAPI TopPerc
    {

    public:
        static void concatMULTISEPeptideIds(std::vector<PeptideIdentification>& all_peptide_ids, std::vector<PeptideIdentification>& new_peptide_ids, String search_engine);
        static void mergeMULTISEPeptideIds(std::vector<PeptideIdentification>& all_peptide_ids, std::vector<PeptideIdentification>& new_peptide_ids);
        static void mergeMULTISEProteinIds(std::vector<ProteinIdentification>& all_protein_ids, std::vector<ProteinIdentification>& new_protein_ids);
        
        static void addMSGFFeatures(std::vector<PeptideIdentification>& peptide_ids, StringList& feature_set);
        static void addXTANDEMFeatures(std::vector<PeptideIdentification>& peptide_ids, StringList& feature_set);
        static void addCOMETFeatures(std::vector<PeptideIdentification>& peptide_ids, StringList& feature_set);
        static void addMASCOTFeatures(std::vector<PeptideIdentification>& peptide_ids, StringList& feature_set);
        static void addMULTISEFeatures(std::vector<PeptideIdentification>& peptide_ids, StringList& search_engines_used, StringList& feature_set);
        static void addCONCATSEFeatures(std::vector<PeptideIdentification>& peptide_id_list, StringList& search_engines_used, StringList& feature_set);
        
        struct lq_ProteinHit
        {
          inline bool operator() (const ProteinHit& h1, const ProteinHit& h2)
          {
            return (h1.getAccession() < h2.getAccession());
          }
        };

        struct lq_PeptideEvidence
        {
          inline bool operator() (const PeptideEvidence& h1, const PeptideEvidence& h2)
          {
            return (h1.getProteinAccession() < h2.getProteinAccession());
          }
        };
        

    protected:
        TopPerc();
        virtual ~TopPerc();
        
        static double rescaleFragmentFeature_(double featureValue, int NumMatchedMainIons);
        static void assignDeltaScore_(std::vector<PeptideHit>& hits, String score_ref, String output_ref);
        static bool hasMHCEnd_(String peptide);
        static String getScanMergeKey_(std::vector<PeptideIdentification>::iterator it, std::vector<PeptideIdentification>::iterator start);

    };

} //namespace OpenMS

#endif //OPENMS_ANALYSIS_ID_TOPPERC_H

