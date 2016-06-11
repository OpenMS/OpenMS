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
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_TOPPERC_H
#define OPENMS_ANALYSIS_ID_TOPPERC_H

#include <vector>
#include <iostream>
#include <cmath>
#include <string>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>

namespace OpenMS
{
    class OPENMS_DLLAPI TopPerc
    {
    public:
      struct PercolatorResult
      {
        String PSMId;
        double score;
        double qvalue;
        double posterior_error_prob;
        String peptide;
        char preAA;
        char postAA;
        StringList proteinIds;

        PercolatorResult(const String& pid, const double s, const double q, const String& p, const char pre, const char pos, const StringList& pl):
            PSMId (pid),
            score (s),
            qvalue (q),
            peptide (p),
            preAA (pre),
            postAA (pos),
            proteinIds (pl)
        {
        }

        PercolatorResult(StringList& row):
        proteinIds()
        {
          // peptide sequence
          StringList pep;
          row[4].split(".", pep);
          //TODO test pep size 3
          peptide = pep[1];
          preAA = pep[0]=="-"?'[':pep[0].c_str()[0];  // const char PeptideEvidence::N_TERMINAL_AA = '[';
          postAA = pep[2]=="-"?']':pep[2].c_str()[0]; // const char PeptideEvidence::C_TERMINAL_AA = ']';
          // SVM-score
          score = row[1].toDouble();
          // q-Value
          qvalue = row[2].toDouble();
          // PEP
          posterior_error_prob = row[3].toDouble();
          // scannr. as written in preparePIN
          PSMId = row[0];
          proteinIds = std::vector<String>(row.begin()+5,row.end());
        }

        bool operator!=(const TopPerc::PercolatorResult& rhs) const
        {
          if (PSMId != rhs.PSMId || score != rhs.score || qvalue != rhs.qvalue ||
              posterior_error_prob != rhs.posterior_error_prob || peptide != rhs.peptide ||
              proteinIds != rhs.proteinIds)
            return true;
          return false;
        }

        bool operator==(const TopPerc::PercolatorResult& rhs) const
        {
          return !(operator !=(rhs));
        }
    };

    public:
        static bool isEnz(const char& n, const char& c, std::string& enz);
        static void prepareCUSTOMpin(std::vector<PeptideIdentification>& peptide_ids, TextFile& txt, std::vector<String>& user_param_features, char out_sep='\t');
        static void prepareMSGFpin(std::vector<PeptideIdentification>& peptide_ids, std::string& enz, TextFile& txt, int min_charge, int max_charge, bool addMHC = false, char out_sep='\t');
        static void prepareXTANDEMpin(std::vector<PeptideIdentification>& peptide_ids, std::string& enz, TextFile& txt, int min_charge, int max_charge, char out_sep='\t');
        static void prepareCOMETpin(std::vector<PeptideIdentification>& peptide_ids, std::string& enz, TextFile& txt, int min_charge, int max_charge, char out_sep='\t');
        static void prepareMASCOTpin(std::vector<PeptideIdentification>& peptide_ids, std::string& enz, TextFile& txt, int min_charge, int max_charge, char out_sep='\t');
        static void prepareMULTIpin(std::vector<PeptideIdentification>& peptide_ids, ProteinIdentification& protein_id, std::string& enz, TextFile& txt, int min_charge, int max_charge, char out_sep='\t');
        static void prepareCONCATpin(std::vector<std::vector<PeptideIdentification> >& peptide_id_list, std::vector<std::vector<ProteinIdentification> >& protein_id_list, std::string& enz, TextFile& txt, int min_charge, int max_charge, char out_sep='\t');
        static size_t countEnzymatic(String peptide, std::string& enz);
        static double rescaleFragmentFeature(double featureValue, int NumMatchedMainIons);
        static String getScanIdentifier(std::vector<PeptideIdentification>::iterator it, std::vector<PeptideIdentification>::iterator start);
        static void assignDeltaScore(std::vector<PeptideHit>& hits, String score_ref);
        static void mergeMULTIids(std::vector<std::vector<ProteinIdentification> >& protein_ids_list, std::vector<std::vector<PeptideIdentification> >& peptide_ids_list, bool skip_checks=false);

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

    private:
        TopPerc();
        virtual ~TopPerc();



    };

} //namespace OpenMS

#endif //OPENMS_ANALYSIS_ID_TOPPERC_H

