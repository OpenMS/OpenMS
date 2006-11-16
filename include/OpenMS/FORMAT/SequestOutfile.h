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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_SEQUESTOUTFILE_H
#define OPENMS_FORMAT_SEQUESTOUTFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

namespace OpenMS
{
	/**
		@brief Representation of an Sequest outfile.
		
		@todo move ensurePathChar (rename to 'ensurePathEnding') to String and add test(Martin)
		@todo check whether every getline may throw an exception (Martin)
		
		@ingroup FileIO
	*/
	class SequestOutfile
	{
		public:
			/// Constructor
			SequestOutfile();
			
			 /**
				@brief loads data from a Mascot outfile
				
				@param filename the file to be loaded
				@param identifications the identifications
				@param protein_identifications the protein identifications
				@param precursor_mz_values the mz values of the precursors corresponding to the identifications
				@param p_value_threshold the significance level (for the peptide hit scores)
				
				This class serves to read in a Sequest outfile. The information can be
				retrieved via the load function.
				
				@ingroup FileIO
				*/
			void
			load(
				const string& result_filename,
				vector< Identification >&	identifications,
				ProteinIdentification&	protein_identification,
// 				vector< Real >& precursor_retention_times,
				vector< Real >& precursor_mz_values,
				const Real& p_value_threshold,
				const string& database = "",
				const string& snd_database = "")
			throw (
				Exception::FileNotFound,
				Exception::ParseError);
			
		//protected:
			/// retrieve columns from a Sequest outfile line
			bool
			getColumns(
				const String& line,
				vector< String >& substrings,
				UnsignedInt number_of_columns,
				UnsignedInt reference_column);
			
			/// retrieve sequences from a FASTA database
			void
			getSequences(
				const String& database_filename,
				const map< String, UnsignedInt >& ac_position_map,
				vector< String >& sequences,
				vector< pair< String, UnsignedInt > >& found,
				map< String, UnsignedInt >& not_found)
			throw (
				Exception::FileNotFound);
				
			/// retrieve the accession type and accession number from a protein description line
			/// (e.g. from FASTA line: >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus], get ac:AAD44166.1 ac type: GenBank)
			void
			getACAndACType(
				String line,
				string& accession,
				string& accession_type);
			
			/// either insert the new peptide hit or update it's protein indices
			bool
			updatePeptideHits(
				PeptideHit& peptide_hit,
				vector< PeptideHit >& peptide_hits);
   };
	
} //namespace OpenMS

#endif // OPENMS_FORMAT_SEQUESTOUTFILE_H
