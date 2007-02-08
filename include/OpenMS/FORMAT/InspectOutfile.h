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

#ifndef OPENMS_FORMAT_INSPECTOUTFILE_H
#define OPENMS_FORMAT_INSPECTOUTFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace OpenMS
{
	/**
		@brief Representation of an Inspect outfile
		
		This class serves to read in an Inspect outfile and write an AnalysisXML file
		
		@ingroup FileIO
	*/
	class InspectOutfile
	{
		public:
			/// Constructor
			InspectOutfile();
			
			/// load the results of an Inspect search
			std::vector< UnsignedInt > load(const std::string& result_filename, std::vector< IdentificationData >& identifications, ProteinIdentification& protein_identification, Real p_value_threshold) throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument);
			
			std::vector< UnsignedInt > getWantedRecords(const std::string& result_filename, Real p_value_threshold) throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument);

			/// generates a trie database from another one, using the wanted records only
			void compressTrieDB(const std::string& database_filename, const std::string& index_filename, std::vector< UnsignedInt >& wanted_records, const std::string& snd_database_filename, const std::string& snd_index_filename, bool append = false) throw (Exception::FileNotFound, Exception::ParseError, Exception::UnableToCreateFile);

			/// generates a trie database from a given one (the type of database is determined by getLabels)
			void generateTrieDB(const std::string& source_database_filename, const std::string& database_filename, const std::string& index_filename, bool append = false, const std::string species = "") throw (Exception::FileNotFound, Exception::ParseError, Exception::UnableToCreateFile);
			

			/// retrieve the accession type and accession number from a protein description line
			/// (e.g. from FASTA line: >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus], get ac:AAD44166.1 ac type: GenBank)
			void getACAndACType(String line, std::string& accession, std::string& accession_type);
			
			/// either insert the new peptide hit or update it's protein indices
			bool updatePeptideHits(PeptideHit& peptide_hit, std::vector< PeptideHit >& peptide_hits);

			/// retvrieve the precursor retention time and mz value
			void getPrecursorRTandMZ(const std::vector< std::pair< String, std::vector< UnsignedInt > > >& files_and_scan_numbers, std::vector< IdentificationData >& ids) throw(Exception::ParseError);

			/// retrieve the labes of a given database (at the moment FASTA and Swissprot)
			void getLabels(const std::string& source_database_filename, std::string& ac_label, std::string& sequence_start_label, std::string& sequence_end_label, std::string& comment_label, std::string& species_label) throw (Exception::FileNotFound, Exception::ParseError);

			/// retrieve sequences from a trie database
			std::vector< UnsignedInt > getSequences(const std::string& database_filename, const std::map< UnsignedInt, UnsignedInt >& wanted_records, std::vector< String >& sequences) throw (Exception::FileNotFound);

			template< typename PeakT > void getExperiment(MSExperiment< PeakT >& exp, String& type, const String& in_filename) throw(Exception::ParseError);

		protected:
			/// a record in the index file that belongs to a trie database consists of three parts
			/// 1) the protein's position in the original database
			/// 2) the proteins's position in the trie database
			/// 3) the name of the protein (the line with the accession identifier)
			static const UnsignedInt db_pos_length_; ///< length of 1)
			static const UnsignedInt trie_db_pos_length_; ///< length of 2)
			static const UnsignedInt protein_name_length_; ///< length of 3)
			static const UnsignedInt record_length_; ///< length of the whole record
			static const char trie_delimiter_; ///< the sequences in the trie database are delimited by this character
			static const std::string score_type_;///< type of score
	};
	
} //namespace OpenMS

#endif // OPENMS_FORMAT_INSPECTOUTFILE_H
