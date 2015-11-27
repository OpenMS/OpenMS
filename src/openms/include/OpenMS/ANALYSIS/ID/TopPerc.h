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
        static bool isEnz(const char& n, const char& c, std::string& enz);
        static void prepareCUSTOMpin(std::vector<PeptideIdentification>& peptide_ids, TextFile& txt, std::vector<String>& user_param_features, char out_sep='\t');
        static void prepareMSGFpin(std::vector<PeptideIdentification>& peptide_ids, std::string& enz, TextFile& txt, int minCharge, int maxCharge, bool addMHC = false, char out_sep='\t');
        static void prepareXTANDEMpin(std::vector<PeptideIdentification>& peptide_ids, std::string& enz, TextFile& txt, int minCharge, int maxCharge, char out_sep='\t');
        static void prepareCOMETpin(std::vector<PeptideIdentification>& peptide_ids, std::string& enz, TextFile& txt, int minCharge, int maxCharge, char out_sep='\t');
        static void prepareMASCOTpin(std::vector<PeptideIdentification>& peptide_ids, std::string& enz, TextFile& txt, int minCharge, int maxCharge, char out_sep='\t');
        static size_t countEnzymatic(String peptide, std::string& enz);
        static double rescaleFragmentFeature(double featureValue, int NumMatchedMainIons);
        static String getScanIdentifier(std::vector<PeptideIdentification>::iterator it, std::vector<PeptideIdentification>::iterator start);
    private:
        TopPerc();
        virtual ~TopPerc();
    };

} //namespace OpenMS

#endif //OPENMS_ANALYSIS_ID_TOPPERC_H

