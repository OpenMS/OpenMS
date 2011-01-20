// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_FALSEDISCOVERYRATE_H
#define OPENMS_ANALYSIS_ID_FALSEDISCOVERYRATE_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Calculates an FDR from identifications
    
		Either two runs of forward and decoy database identification or
		one run containing both (with marks) can be used to annotate
		each of the peptide hits with a FDR.
		
    Also q-values can be reported instead of p-values.
    q-values are basically only adjusted p-values, also ranging from 0 to 1, with lower values being preferable.
    When looking at the list of hits ordered by q-values, then a hit with q-value of @em x means that there is an
    @x*100 percent chance that all hits with a q-value <= @em x are a false positive hit.

		@todo implement combined searches properly (Andreas)
		@improvement implement charge state separated fdr/q-values (Andreas)

		@htmlinclude OpenMS_FalseDiscoveryRate.parameters

		@ingroup Analysis_ID
  */
  class OPENMS_DLLAPI FalseDiscoveryRate
  	: public DefaultParamHandler
  {
  	public: 
	  	///Default constructor
	  	FalseDiscoveryRate();
  		
  		/**
  			@brief Calculates the FDR of two runs, a forward run and a decoy run on peptide level

				@param fwd_ids forward peptide identifications
				@param rev_ids reverse peptide identifications
  		*/
  		void apply(std::vector<PeptideIdentification>& fwd_ids, std::vector<PeptideIdentification>& rev_ids);

			/**
				@brief Calculates the FDR of one run from a concatenated sequence db search

				@param ids peptide identifications, containing target and decoy hits
			*/
			void apply(std::vector<PeptideIdentification>& id);

			/**
			 	@brief Calculates the FDR of two runs, a forward run and decoy run on protein level

				@param fwd_ids forward protein identifications
				@param rev_ids reverse protein identifications
			*/
			void apply(std::vector<ProteinIdentification>& fwd_ids, std::vector<ProteinIdentification>& rev_ids);

			/** 
				@brief Calculate the FDR of one run from a concatenated sequence db search
				
				
				@param ids protein identifications, containing target and decoy hits
			*/
			void apply(std::vector<ProteinIdentification>& ids);

  	private:
  		///Not implemented
  		FalseDiscoveryRate(const FalseDiscoveryRate&);
  		
			///Not implemented
			FalseDiscoveryRate& operator = (const FalseDiscoveryRate&);

			/// calculates the fdr stored into fdrs, given two vectors of scores
			void calculateFDRs_(Map<DoubleReal, DoubleReal>& score_to_fdr, std::vector<DoubleReal>& target_scores, std::vector<DoubleReal>& decoy_scores, bool q_value, bool higher_score_better);

  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_FALSEDISCOVERYRATE_H
