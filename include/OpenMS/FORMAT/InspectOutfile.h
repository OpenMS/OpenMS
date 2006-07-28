// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: InspectOutfile.h,v 1.0 2006/07/12 15:58:59 martinlangwisch Exp $
// $Author: martinlangwisch $
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_INSPECTOUTFILE_H
#define OPENMS_FORMAT_INSPECTOUTFILE_H

#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/FORMAT/InspectFile.h>

namespace OpenMS
{
  /**
    @brief Representation of an Inspect outfile
    
    This class serves to read in an Inspect outfile. The information can be 
    retrieved via the >> operator. 
  
  	@ingroup FileIO
  */
  class InspectOutfile: 
		InspectFile
  {
    public:
      /// Constructor
			InspectOutfile();
			
			/// load the results of an InsPecT search
      void load(const std::string& result_filename, std::vector< Identification >&	identifications, ProteinIdentification&	protein_identification, std::vector< float >& 	precursor_retention_times, std::vector< float >& precursor_mz_values, const double& p_value_threshold, const std::string& database_filename, const std::string& database_path, const std::string& sequence_filename, std::string index_filename = "") throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument);
      
		protected:			
			/// get the accession and accession type of a protein
			void getACAndACType(String line, std::string& accession, std::string& accession_type) throw (Exception::ParseError);
			
			/// given a vector of peptide hits, either insert the new peptide hit or update its ProteinHits, returns whether an update took place
			bool updatePeptideHits(PeptideHit& peptide_hit, std::vector< PeptideHit >& peptide_hits);
		
   };
	
} //namespace OpenMS

#endif // OPENMS_FORMAT_INSPECTOUTFILE_H
