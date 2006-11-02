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

#ifndef OPENMS_FORMAT_INSPECTFILE_H
#define OPENMS_FORMAT_INSPECTFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <fstream>
#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>

namespace OpenMS
{
	/**
		@brief 
		
		@todo add docu
		@todo write test(Martin)
	*/
	class InspectFile
	{
		public:
			/// @todo move to String (Martin)
			void remove(std::string& s, const std::string& unwanted_character);
			
			/// erase all whitespaces in a string
			/// @todo move to String (Martin)
			void removeWhitespaces(std::string& s);
			
			/// used to generate a trie database from another kind of database
			void generateTrieDB(const std::string& source_filename, const std::string& source_path, const std::string& database_path, std::vector< unsigned int > wanted_records = std::vector< unsigned int >(), std::string database_filename = "", std::string index_filename = "", bool append = false, std::string species = "None") throw (Exception::FileNotFound, Exception::ParseError);
			
			/// compress a trie database to contain the wanted records only
			void compressTrieDB(const std::string database_filename, std::string index_filename, std::string database_path, std::vector< unsigned int > wanted_records, std::string second_database_filename = "", std::string second_index_filename = "", std::string second_database_path = "", bool append = false) throw (Exception::FileNotFound, Exception::ParseError);
			
			/// get the sequence, accession and accession type for some proteins from a database
			void getSequenceAndACandACType(const std::string& database_filename, std::vector< unsigned int > wanted_records, std::vector< std::vector< String > >& protein_info, const std::string& ac_label, const std::string& sequence_start_label, const std::string& sequence_end_label, const std::string& comment_label, std::string species_label , std::string species = "None") throw (Exception::FileNotFound, Exception::ParseError);
			
			/**
				@brief retrieve the labels used in a file
				
				id line
				label for the start of the sequence
				label for the end of the sequence
				label for the line defining the species
			*/
			void getLabels(const std::string& source_filename, std::string& ac_label, std::string& sequence_start_label, std::string& sequence_end_label, std::string& comment_label, std::string& species_label) throw (Exception::FileNotFound, Exception::ParseError);
			
			/// retrieve sequences from a trie database
			void getSequences(std::string database_path, const std::string& database_filename, std::string index_filename, const std::vector< unsigned int > wanted_records, std::vector<std::string>& sequences) throw (Exception::FileNotFound, Exception::ParseError);
			
			/// if there's no character marking the name as a directory it is set
			void ensurePathChar(std::string& path, char path_char = '/');
			
			bool setTempDatabaseFilename(const std::string& temp_database_filename);
			const std::string& getTempDatabaseFilename() const;
			
			bool setTempIndexFilename(const std::string& temp_index_filename);
			const std::string& getTempIndexFilename() const;
			
			bool setSecondTempDatabaseFilename(const std::string& temp_second_database_filename);
			const std::string& getSecondTempDatabaseFilename() const;
			
			bool setSecondTempIndexFilename(const std::string& temp_second_index_filename);
			const std::string& getSecondTempIndexFilename() const;
			
			void setUseTempFiles(const bool use_temp_files);
			const bool getUseTempFiles() const;
		
		protected:
			/**
				@brief can convert any file format fullfilling the following conditions to the ".trie" database format (and creating a corresponding index file)
				
				for each record in the file there has to be
				- a unique* label for the id line (the label and the id have to be seperated by at least one whitespace)
				- a unique* label for the start of the sequence (a line preceeding the sequence lines!)
				- a unique* label for the end of the sequence (a line succeeding the sequence lines!)
				- there may be a unique label for the line defining the species
				the lines labeled correspondingly have to start with the label (whitespaces before are allowed)
				* unique means there must not be more than one line with this label, there may be one line (and thus one label) for several things, though; e.g. in a FASTA file a ">" marks the id, the start and the end of a sequence:
				> very_nice_protein_indeed
				MYVERYNICEPRQTEIN
				> some_other_protein
				...
			*/
			void generateTrieDB_(const std::string& source_filename, const std::string& source_path, const std::string& database_path, const std::string& ac_label, const std::string& sequence_start_label, const std::string& sequence_end_label, const std::string& comment_label, std::string species_label = "", std::string species = "None", std::vector< unsigned int > wanted_records = std::vector< unsigned int >(), std::string database_filename = "", std::string index_filename = "", bool append = false) throw (Exception::FileNotFound, Exception::ParseError);
			
			static const unsigned int index_peptide_name_length_; ///< the length of a peptide name in the index file
			static const unsigned int index_db_record_length_; ///< the length of the originial database position in the index file
			static const unsigned int index_trie_record_length_; ///< the length of the trie database position in the index file
			static const unsigned int index_record_length_; ///< the length of a record (the sum of the three parameters above)
			static const char trie_delimiter_; ///< the character that delimits protein sequences in a trie database
			static const std::string score_type_;
			
			static std::string temp_database_filename_;
			static std::string temp_index_filename_;
			static std::string temp_second_database_filename_;
			static std::string temp_second_index_filename_;
			static bool use_temp_files_;
		
	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_INSPECTFILE_H
