// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

    
    Size findMinimalProteinList(const std::vector<PeptideIdentification>& peptide_ids);

    void calculateProteinProbabilities(const std::vector<PeptideIdentification>& ids);
    
//     DoubleReal getProteinProbability(const String& acc,const std::vector<String>& accessions, const std::vector<DoubleReal>& probabilities);

    DoubleReal getProteinProbability(const String& acc);

    bool isProteinInMinimalList(const String& acc);
    Int getNumberOfProtIds(DoubleReal protein_id_threshold);
    Int getNumberOfProtIdsPeptideRule(Int min_peptides,	std::map<String,std::set<String> >& prot_id_counter);
    
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
