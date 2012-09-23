// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Alexandra Zerck $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_ANALYSIS_TARGETED_PSPROTEININFERENCE_H
#define OPENMS_ANALYSIS_TARGETED_PSPROTEININFERENCE_H

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>

namespace OpenMS
{

  /**
       @brief This class implements protein inference for the precursor ion selection strategies.





  */
  class OPENMS_DLLAPI PSProteinInference
  {
public:

    PSProteinInference();

    virtual ~PSProteinInference();


    Size findMinimalProteinList(const std::vector<PeptideIdentification> & peptide_ids);

    void calculateProteinProbabilities(const std::vector<PeptideIdentification> & ids);

//     DoubleReal getProteinProbability(const String& acc,const std::vector<String>& accessions, const std::vector<DoubleReal>& probabilities);

    DoubleReal getProteinProbability(const String & acc);

    bool isProteinInMinimalList(const String & acc);
    Int getNumberOfProtIds(DoubleReal protein_id_threshold);
    Int getNumberOfProtIdsPeptideRule(Int min_peptides, std::map<String, std::set<String> > & prot_id_counter);

    void setSolver(LPWrapper::SOLVER solver)
    {
      solver_ = solver;
    }

    LPWrapper::SOLVER getSolver()
    {
      return solver_;
    }

private:
    std::vector<String> minimal_protein_list_accessions_;
    std::vector<String> accessions_;
    std::vector<DoubleReal> probabilities_;
    LPWrapper::SOLVER solver_;
  };

}



#endif // #ifndef OPENMS_ANALYSIS_TARGETED_PSPROTEININFERENCE_H
