/// -*- Mode: C++; tab-width: 2; -*-
/// vi: set ts=2:
///
/// --------------------------------------------------------------------------
///                   OpenMS Mass Spectrometry Framework
/// --------------------------------------------------------------------------
///  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
///
///  This library is free software; you can redistribute it and/or
///  modify it under the terms of the GNU Lesser General Public
///  License as published by the Free Software Foundation; either
///  version 2.1 of the License, or (at your option) any later version.
///
///  This library is distributed in the hope that it will be useful,
///  but WITHOUT ANY WARRANTY; without even the implied warranty of
///  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
///  Lesser General Public License for more details.
///
///  You should have received a copy of the GNU Lesser General Public
///  License along with this library; if not, write to the Free Software
///  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
///
/// --------------------------------------------------------------------------
/// $Id: InspectFile.h,v 1.0 2006/07/12 15:58:59 martinlangwisch Exp $
/// $Author: martinlangwisch $
/// $Maintainer: Martin Langwisch $
/// --------------------------------------------------------------------------

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
	class InspectFile
	{
		public:
			
			/// erase all given characters in a string ### in string-klasse einbauen?
			void erase(std::string& s, const std::string& unwanted_character);
			
			/// erase all whitespaces in a string ### in string-klasse einbauen?
			void erase_whitespaces(std::string& s);
			
			/// used to generate a trie database from another kind of database
			void compressor(const std::string& source_filename, const std::string& source_path, const std::string& database_path, std::vector< unsigned int > wanted_records = std::vector< unsigned int >(), std::string database_filename = "", std::string index_filename = "", bool append = false, std::string species = "None") throw (Exception::FileNotFound, Exception::ParseError);
			
			void compressTrieDB(const std::string database_filename, std::string index_filename, std::string database_path, std::vector< unsigned int > wanted_records, std::string second_database_filename = "", std::string second_index_filename = "", std::string second_database_path = "", bool append = false) throw (Exception::FileNotFound, Exception::ParseError);
			
			void getSeparators(const std::string& source_filename, std::string& ac_label, std::string& sequence_start_label, std::string& sequence_end_label, std::string& comment_label, std::string& species_label) throw (Exception::FileNotFound, Exception::ParseError);
			
			void getSequences(std::string database_path, const std::string& database_filename, std::string index_filename, const std::vector< unsigned int > wanted_records, std::vector<std::string>& sequences) throw (Exception::FileNotFound, Exception::ParseError);
			
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
			
			void ensurePathChar(std::string& path, char path_char = '/');
		
		protected:
			/* can convert any file format fullfilling the following conditions to the ".trie" database format (and creating a corresponding index file)
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
			void compressor_(const std::string& source_filename, const std::string& source_path, const std::string& database_path, const std::string& ac_label, const std::string& sequence_start_label, const std::string& sequence_end_label, const std::string& comment_label, std::string species_label = "", std::string species = "None", std::vector< unsigned int > wanted_records = std::vector< unsigned int >(), std::string database_filename = "", std::string index_filename = "", bool append = false) throw (Exception::FileNotFound, Exception::ParseError);
			
			static const unsigned int index_peptide_name_length_;
			static const unsigned int index_db_record_length_;
			static const unsigned int index_trie_record_length_;
			static const unsigned int index_record_length_;
			static const char trie_delimiter_;
			static const std::string score_type_;
			
			static std::string temp_database_filename_;
			static std::string temp_index_filename_;
			static std::string temp_second_database_filename_;
			static std::string temp_second_index_filename_;
			static bool use_temp_files_;
		
	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_INSPECTFILE_H
