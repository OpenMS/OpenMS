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
//#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/DATASTRUCTURES/String.h>
//#include <OpenMS/FORMAT/MzDataFile.h>
//#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Identification.h>
//#include <OpenMS/METADATA/PeptideHit.h>
//#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
//#include <OpenMS/FORMAT/FileHandler.h>

//#include <fstream>
//#include <iostream>
//#include <map>
//#include <set>
//#include <vector>

namespace OpenMS
{
	/**
		@brief Representation of an Inspect outfile
		
		This class serves to read in an Inspect outfile and write an IdXML file
		
  	@todo Handle Modifications (Andreas)
		
		@ingroup FileIO
	*/
	class InspectOutfile
	{
		public:
			/// Constructor
			InspectOutfile();
			
			/// load the results of an Inspect search
			std::vector<UInt> load(const String& result_filename, std::vector<PeptideIdentification>& peptide_identifications, Identification& protein_identification, Real p_value_threshold) throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument);
			
			std::vector< UInt > getWantedRecords(const String& result_filename, Real p_value_threshold) throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument);

			/// generates a trie database from another one, using the wanted records only
			void compressTrieDB(const String& database_filename, const String& index_filename, std::vector<UInt>& wanted_records, const String& snd_database_filename, const String& snd_index_filename, bool append = false) throw (Exception::FileNotFound, Exception::ParseError, Exception::UnableToCreateFile);

			/// generates a trie database from a given one (the type of database is determined by getLabels)
			void generateTrieDB(const String& source_database_filename, const String& database_filename, const String& index_filename, bool append = false, const String species = "") throw (Exception::FileNotFound, Exception::ParseError, Exception::UnableToCreateFile);
			

			/// retrieve the accession type and accession number from a protein description line
			/// (e.g. from FASTA line: >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus], get ac:AAD44166.1 ac type: GenBank)
			void getACAndACType(String line, String& accession, String& accession_type);
			
			/// either insert the new peptide hit or update it's protein indices
			//bool updatePeptideHits(PeptideHit& peptide_hit, std::vector<PeptideHit>& peptide_hits);

			/// retvrieve the precursor retention time and mz value
			void getPrecursorRTandMZ(const std::vector<std::pair<String, std::vector<UInt> > >& files_and_scan_numbers, std::vector<PeptideIdentification>& ids) throw(Exception::ParseError);

			/// retrieve the labes of a given database (at the moment FASTA and Swissprot)
			void getLabels(const String& source_database_filename, String& ac_label, String& sequence_start_label, String& sequence_end_label, String& comment_label, String& species_label) throw (Exception::FileNotFound, Exception::ParseError);

			/// retrieve sequences from a trie database
			std::vector<UInt> getSequences(const String& database_filename, const std::map<UInt, UInt>& wanted_records, std::vector<String>& sequences) throw (Exception::FileNotFound);

			template< typename PeakT > void getExperiment(MSExperiment< PeakT >& exp, String& type, const String& in_filename) throw(Exception::ParseError);

		protected:
			/// a record in the index file that belongs to a trie database consists of three parts
			/// 1) the protein's position in the original database
			/// 2) the proteins's position in the trie database
			/// 3) the name of the protein (the line with the accession identifier)
			static const UInt db_pos_length_; ///< length of 1)
			static const UInt trie_db_pos_length_; ///< length of 2)
			static const UInt protein_name_length_; ///< length of 3)
			static const UInt record_length_; ///< length of the whole record
			static const char trie_delimiter_; ///< the sequences in the trie database are delimited by this character
			static const String score_type_;///< type of score
	};
	
} //namespace OpenMS

#endif // OPENMS_FORMAT_INSPECTOUTFILE_H
