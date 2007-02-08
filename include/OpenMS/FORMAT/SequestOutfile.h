// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

namespace OpenMS
{
	/**
    @brief Representation of a Sequest output file
    
    This class serves to read in a Sequest outfile. The information can be
    retrieved via the load function.
  	
  	@ingroup FileIO
  */
	class SequestOutfile
	{
		public:
			/// Constructor
			SequestOutfile();
			
			 /**
				@brief loads data from a Inspect outfile
				
				@param result_filename the file to be loaded
				@param identifications the identifications
				@param protein_identification the protein identifications
				@param p_value_threshold the significance level (for the peptide hit scores)
			  @param pvalues a list with the pvalues of the peptides (pvalues computed with peptide prophet)
			  @param database the database used for the search
				
				This class serves to read in a Sequest outfile. The information can be
				retrieved via the load function.
			*/
			void load(const std::string& result_filename, std::vector< IdentificationData >& identifications, ProteinIdentification& protein_identification, const Real& p_value_threshold, std::vector< Real >& pvalues, const std::string& database = "") throw (Exception::FileNotFound, Exception::ParseError);

			
			void finishSummaryHtml(const std::string& summary_filename) throw (Exception::UnableToCreateFile);

			/// write a
			void out2SummaryHtml(std::string out_filename, const std::string& summary_filename, const std::string& database_filename, bool& append) throw(Exception::FileNotFound, Exception::ParseError, Exception::UnableToCreateFile);

			std::map< String, std::vector< Real > > getPeptidePValues(const std::string& out_dir, const std::string& prob_filename) throw (Exception::FileNotFound);
			
			/// retrieve columns from a Sequest outfile line
			bool getColumns(const String& line, std::vector< String >& substrings, UnsignedInt number_of_columns, UnsignedInt reference_column);
			
			/// retrieve sequences from a FASTA database
			void getSequences(const String& database_filename, const std::map< String, UnsignedInt >& ac_position_map, std::vector< String >& sequences, std::vector< std::pair< String, UnsignedInt > >& found, std::map< String, UnsignedInt >& not_found) throw (Exception::FileNotFound);
				
			/// retrieve the accession type and accession number from a protein description line
			/// (e.g. from FASTA line: >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus], get ac:AAD44166.1 ac type: GenBank)
			void getACAndACType(String line, std::string& accession, std::string& accession_type);
			
			/// either insert the new peptide hit or update it's protein indices
			bool updatePeptideHits(PeptideHit& peptide_hit, std::vector< PeptideHit >& peptide_hits);

			void readOutHeader(const std::string& result_filename, DateTime& datetime, Real& precursor_mz_value, SignedInt& charge, UnsignedInt& precursor_mass_type, UnsignedInt& ion_mass_type, SignedInt& number_column, SignedInt& rank_sp_column, SignedInt& id_column, SignedInt& mh_column, SignedInt& delta_cn_column, SignedInt& xcorr_column, SignedInt& sp_column, SignedInt& sf_column, SignedInt& ions_column, SignedInt& reference_column, SignedInt& peptide_column, SignedInt& score_column, UnsignedInt& number_of_columns, UnsignedInt& displayed_peptides) throw(Exception::FileNotFound, Exception::ParseError);
   };
	
} //namespace OpenMS

#endif // OPENMS_FORMAT_SEQUESTOUTFILE_H
