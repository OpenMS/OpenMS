// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Alexandra Zerck $
// --------------------------------------------------------------------------
//

#pragma once

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

//     double getProteinProbability(const String& acc,const std::vector<String>& accessions, const std::vector<double>& probabilities);

    double getProteinProbability(const String & acc);

    bool isProteinInMinimalList(const String & acc);
    Int getNumberOfProtIds(double protein_id_threshold);
    Int getNumberOfProtIdsPeptideRule(Int min_peptides, std::map<String, std::set<String> > & prot_id_counter);

    LPWrapper::SOLVER getSolver()
    {
      return solver_;
    }

private:
    std::vector<String> minimal_protein_list_accessions_;
    std::vector<String> accessions_;
    std::vector<double> probabilities_;
    LPWrapper::SOLVER solver_;
  };

}



