// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_FALSEIDENTIFICATIONRATE_H
#define OPENMS_ANALYSIS_ID_FALSEIDENTIFICATIONRATE_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Calculates an FDR from identifications
    
		Either two runs of forward and decoy database identification or
		one run containing both (with marks) can be used to annotate
		each of the peptide hits with a FDR.
		
		@ref FalseDiscoveryRate_Parameters are explained on a separate page.

		@ingroup Analysis_ID
  */
  class FalseDiscoveryRate
  	: public DefaultParamHandler
  {
  	public: 
	  	///Default constructor
	  	FalseDiscoveryRate();
  		
  		/**
  			@brief Calculates the FDR of two runs, a forward run and a decoy run
  			
  			
  		*/
  		void apply(std::vector<PeptideIdentification>& fwd_ids, std::vector<PeptideIdentification>& rev_ids);

			/**
				@brief Calculates the FDR of the run containing both, decoy and forward hits

				@argument marker gives the marker used in the protein title to inidicate a decoy protein
			*/
			void apply(std::vector<PeptideIdentification>& ids);
  		
  	private:
  		///Not implemented
  		FalseDiscoveryRate(const FalseDiscoveryRate&);
  		
			///Not implemented
			FalseDiscoveryRate& operator = (const FalseDiscoveryRate&);
			
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_FALSEDISCOVERYRATE_H
