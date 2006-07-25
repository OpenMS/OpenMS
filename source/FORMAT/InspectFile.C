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
// $Id: InspectFile.C,v 1.0 2006/07/25 13:46:15 martinlangwisch Exp $
// $Author: martinlangwisch $
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/InspectFile.h>

namespace OpenMS
{
	void InspectFile::remove(std::string& s, const std::string& unwanted_characters)
	{
		std::string::iterator start = s.begin();
		std::string::iterator end = s.begin();
		while ( (end != s.end()) && (start != s.end()) )
		{
			// find the next unwanted character
			while ( (start != s.end()) && (unwanted_characters.find(*start, 0) == std::string::npos ) )
			{
				++start;
			}
			end = start;
		
			// find the next symbol that is not an unwanted character
			while ( (end != s.end()) && (unwanted_characters.find(*end, 0) != std::string::npos ) )
			{
				++end;
			}
			// remove the unwanted characters
			start = s.erase(start, end);
			// find the next unwanted character
		}
	}
	
	void InspectFile::removeWhitespaces(std::string& s)
	{
		std::string::iterator start = s.begin();
		std::string::iterator end = s.begin();
		while ( (end != s.end()) && (start != s.end()) )
		{
			// find the next symbol that is not a whitespace or an asterisk
			while ( (end != s.end()) && (*end < 33) )
			{
				++end;
			}
			// remove the whitespaces or asterisks
			start = s.erase(start, end);
			// find the next whitespace
			while ( (start != s.end()) && (*start > 32) )
			{
				++start;
			}
			end = start;
		}
	}
	
	void InspectFile::compressTrieDB(const std::string database_filename, std::string index_filename, std::string database_path, std::vector< unsigned int > wanted_records, std::string second_database_filename_, std::string second_index_filename, std::string second_database_path, bool append) throw (Exception::FileNotFound, Exception::ParseError)
	{
		ensurePathChar(database_path);
		std::string path_and_file = database_path + database_filename;
		std::ifstream database_file( path_and_file.c_str(), std::ios::in | std::ios::binary );
		if ( !database_file ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
		// if no index filename was given, assume that it is the same as database_filename but with ending ".index"
		if ( index_filename.empty() ) index_filename = database_filename.substr(0, database_filename.length()-4) + "index";
		
		path_and_file = database_path+index_filename;
		std::ifstream index_file(path_and_file.c_str(), std::ios::in | std::ios::binary);
		if ( !index_file ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, index_filename);
		
		String second_database_filename(second_database_filename_);
		
		// if no name for the second database is given, the name of the database with ending "snd.trie" is used
		if ( second_database_filename.empty() )
		{
			if ( getUseTempFiles() ) second_database_filename = getSecondTempDatabaseFilename();
			else
			{
				second_database_filename = database_filename;
				if ( second_database_filename.hasSuffix(".trie") ) second_database_filename.insert(second_database_filename.length()-1-4, ".snd");
				else second_database_filename.append(".snd.trie");
			}
		}
		// if no name for the index is given, the name of the second database with ending ".index" is used
		if ( second_index_filename.empty() )
		{
			if ( getUseTempFiles() ) second_index_filename = getSecondTempIndexFilename();
			else
			{
				second_index_filename = second_database_filename;
				second_index_filename.replace(second_database_filename.length()-4, 4, "index");
			}
		}
		if ( second_database_path.empty() )
		{
			second_database_path = database_path;
		}
		else ensurePathChar(second_database_path);
		
		if ( (database_filename == second_database_filename) && (database_path == second_database_path) && (index_filename == second_index_filename) )
		{
			database_file.close();
			index_file.close();
			return;
		}
		
		char openmode[3] = "wb";
		openmode[2] = 0;
		if ( append ) openmode[1] = 'a';
		path_and_file = second_database_path+second_database_filename;
		FILE* second_database_file = fopen( path_and_file.c_str(), openmode);
		if ( second_database_file == NULL ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, second_database_filename);
		
		// if in append mode, check whether the file is empty and if so, change to write mode
		if ( append )
		{
			fseek(second_database_file, 0, SEEK_END);
			if ( ftell(second_database_file) )	fseek(second_database_file, 0, SEEK_SET);
			else
			{
				fclose(second_database_file);
				append = false;
				openmode[0] = 'w';
				second_database_file = fopen(path_and_file.c_str(), openmode);
				if ( second_database_file == NULL ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, second_database_filename);
			}
		}
		
		path_and_file = second_database_path+second_index_filename;
		FILE* second_index_file = fopen( path_and_file.c_str(), openmode);
		if ( second_index_file == NULL ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, second_index_filename);
		
		// determine the length of the index file
		index_file.seekg(0, std::ios::end);
		unsigned int file_length = index_file.tellg();
		index_file.seekg(0, std::ios::beg);
		
		// if all records are selected, just copy the files
		if ( wanted_records.empty() )
		{
			database_file.seekg(0, std::ios::end);
			unsigned int db_file_length = database_file.tellg();
			database_file.seekg(0, std::ios::beg);
			
			void* buffer = malloc( std::max(db_file_length, file_length) );
			database_file.readsome((char*) buffer, db_file_length );
			database_file.close();
			
			fwrite(buffer, db_file_length, 1, second_database_file);
			fclose(second_database_file);
			
			index_file.readsome((char*) buffer, file_length);
			index_file.close();
			
			fwrite(buffer, file_length, 1, second_index_file);
			fclose(second_index_file);
			
			free(buffer);
			
			return;
		}
		
		// write the protein sequences to the new database
		char* protein_name = new char[index_peptide_name_length_+1];
		protein_name[index_peptide_name_length_] = 0;
		char* file_pos_from_index_file = new char[index_db_record_length_];
		std::stringbuf sequence;
		unsigned int index_pos, second_database_file_pos, database_pos;
		
		for ( unsigned int i = 0; i < wanted_records.size(); ++i)
		{
			// get the according record in the index file (a record has length 8+4+80 and the database position is stored after 8 bytes in a record and is 4 bytes long)
			index_pos = wanted_records[i] * (index_record_length_);
			index_file.seekg(index_pos);
			// if the file is too short
			if ( file_length < (index_pos+index_trie_record_length_+index_db_record_length_+index_peptide_name_length_) )
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "index file is too short!", index_filename);
			}
			
			// read the index of the original file and write it to the new database
			index_file.readsome(file_pos_from_index_file, index_db_record_length_);
			fwrite(file_pos_from_index_file, index_db_record_length_, 1, second_index_file);
			
			// all but the first sequence are preceded by an asterisk
			if ( i || append ) fputc('*', second_database_file);
			
			// write the actual position in the new database
			second_database_file_pos = ftell(second_database_file);
			fwrite(&second_database_file_pos, index_trie_record_length_, 1, second_index_file);
			
			// read the sequence;
			index_file.readsome((char*) &database_pos, index_trie_record_length_);
			database_file.seekg(database_pos);
			database_file.get(sequence, trie_delimiter_);
			fputs(sequence.str().c_str(), second_database_file);
			sequence.str("");
			
			// read the protein name and write it to the new database
			index_file.readsome(protein_name, 80);
			fwrite(protein_name, 80, 1, second_index_file);
		}
		
		delete(protein_name);
		delete(file_pos_from_index_file);
		
		index_file.close();
		database_file.close();
		fclose(second_database_file);
		fclose(second_index_file);
	}
	
	void InspectFile::generateTrieDB_(const std::string& source_filename, const std::string& source_path, const std::string& database_path, const std::string& ac_label, const std::string& sequence_start_label, const std::string& sequence_end_label, const std::string& comment_label, std::string species_label, std::string species, std::vector< unsigned int > wanted_records, std::string database_filename, std::string index_filename, bool append) throw (Exception::FileNotFound, Exception::ParseError)
	{
		// if no database name is given, the name of the source file plus ending ".trie" is used
		if ( database_filename.empty() )
		{
			if ( getUseTempFiles() ) database_filename = getTempDatabaseFilename();
			else database_filename = source_filename+".trie";
		}
		// if no index name is given, the name of the source file plus ending ".index" is used
		if ( index_filename.empty() )
		{
			if ( getUseTempFiles() ) index_filename = getTempIndexFilename();
			else index_filename = source_filename+".index";
		}
		std::string path_and_file = source_path;
		ensurePathChar(path_and_file);
		path_and_file.append(source_filename);
		std::ifstream source_file(path_and_file.c_str(), std::ifstream::in | std::ifstream::binary );
		if ( !source_file ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, source_filename);
		
		char openmode[3] = "wb";
		openmode[2] = 0;
		if ( append ) openmode[0] = 'a';
		
		path_and_file = database_path;
		ensurePathChar(path_and_file);
		path_and_file.append(database_filename);
		FILE* database_file = fopen(path_and_file.c_str(), openmode);
		if ( database_file == NULL ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
		
		// if in append mode, check whether the file is empty and if so set append to false
		if ( append )
		{
			//fseek(database_file, 0, SEEK_END);
			if ( !ftell(database_file) )	append = false;
			//if ( ftell(database_file) )	fseek(database_file, 0, SEEK_SET);
			/*else
			{
				fclose(database_file);
				append = false;
				openmode[0] = 'w';
				database_file = fopen(path_and_file.c_str(), openmode);
				if ( database_file == NULL ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
			}*/
		}
		
		path_and_file = database_path;
		ensurePathChar(path_and_file);
		path_and_file.append(index_filename);
		FILE* index_file = fopen(path_and_file.c_str(), openmode);
		if ( index_file == NULL ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, index_filename);
		unsigned char ac_flag = 1;
		unsigned char species_flag = 2;
		unsigned char sequence_flag = 4;
		unsigned char record_flags = (species == "None") ? species_flag : 0; // what records have already been set
		unsigned int written_records = false; // whether records have already been written
		unsigned int seen_records = 0;
		unsigned int pos; // the position in a line
		unsigned long long int source_file_pos; // the position of the current line in the source file (has to be determined before reading in the corresponding line)
		unsigned long long int source_file_pos_buffer = source_file.tellg();
		unsigned int database_file_pos; // the position of the start of the current line in the database file
		String line; // the current line read
		String sequence; // the sequence of a record
		char* ac = new char[index_peptide_name_length_]; //  the name of the record
		
		std::vector< unsigned int >::const_iterator i = wanted_records.begin();
		
		while ( getline(source_file, line) && ( wanted_records.empty() || (written_records < wanted_records.size())) )
		{
			if ( !line.empty() ) line.resize(line.length()-1);
			line.trim();
			
			// empty and comment lines are skipped
			if ( line.empty() || line.hasPrefix(comment_label) ) continue;
			
			// read the sequence
			if ( record_flags==(ac_flag|species_flag|sequence_flag) )
			{
				if ( line.hasPrefix(sequence_end_label) )
				{
					// if the sequence is not empty, the record has the correct form
					if ( !sequence.empty() )
					{
						// write the index and the database (for the wanted entries only)
						if ( wanted_records.empty() || (*i == seen_records) )
						{
							if ( !wanted_records.empty() ) ++i;
							// all but the first record in the database are preceded by an asterisk (if in append mode an asterisk has to be put at any time)
							if ( written_records || append ) fputc('*', database_file);
							
							database_file_pos = ftell(database_file);
							// a record in the index file consists of eight bytes for the record position in the source file, four bytes for the record position in the database file and 80 bytes for the record name (so one record is 92 bytes long, there are no separators in the index file!
							fwrite(&source_file_pos, sizeof(unsigned long long int), 1, index_file);
							fwrite(&database_file_pos, sizeof(unsigned int), 1, index_file);
							fwrite(ac, 1, index_peptide_name_length_, index_file);
							fputs(sequence.c_str(), database_file);
							++written_records;
						}
						
						// stop when the all wanted records are found
						if ( !wanted_records.empty() && (i == wanted_records.end()) ) break;
					}
					// if the sequence is empty
					else throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "emtpy sequence!" , source_filename);
					
					sequence.clear();
					// set back the record flags for a new record
					record_flags = (species == "None") ? species_flag : 0;
					++seen_records;
				}
				else
				{
					// erase all whitespaces from the sequence
					removeWhitespaces(line);
					remove(line, "*");
					// save this part of the sequence
					sequence.append(line);
				}
			}
			if ( !(record_flags&sequence_flag) )
			{
				if ( line.hasPrefix(ac_label) )
				{
					// find the beginning of the accession
					pos = ac_label.length();
					// discard the whitespaces after the label
					while ( (line.length() > pos) && (line[pos] < 33) )
					{
						++pos;
					}
					if ( pos != line.length() )
					{
						// clear the ac-string
						memset(ac, 0, index_peptide_name_length_);
						// read at most 80 characters from the record name
						memcpy(ac, line.substr(pos, std::min((size_t) index_peptide_name_length_, line.length()-pos)).c_str(), std::min((size_t) index_peptide_name_length_, line.length()-pos));
						source_file_pos = source_file_pos_buffer;
					}
					else
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "record has no accession!" , source_filename);
					}
					record_flags |= ac_flag; // set the ac flag
				}
				// if a species line is found and a species (other than "None") is given, check whether this record is from the wanted species ### None auch global setzten
				// (if no species is given, the species label ("") is always found, but the compare will always give zero as the species then will be "None")
				if ( line.hasPrefix(species_label) && (record_flags==ac_flag) )
				{
					pos = species_label.length();
					if ( line.find(species, pos) != std::string::npos ) record_flags |= species_flag;
					
					// if it's not from the wanted species, skip the record
					if ( !(record_flags&species_flag) ) record_flags = 0;
				}
				// if the beginning of the sequence is found
				if ( line.hasPrefix(sequence_start_label) && (record_flags&ac_flag) && (record_flags&species_flag) ) record_flags |= sequence_flag;
			}
			source_file_pos_buffer = source_file.tellg();
		} // source file read
		
		// if the last record has no sequence end label, the sequence has to be appended nevertheless
		if ( record_flags==(ac_flag|species_flag|sequence_flag) )
		{
			// if the sequence is not empty, the record has the correct form
			if ( !sequence.empty() )
			{
				// write the index and the database
				if ( wanted_records.empty() || (*i == seen_records) )
				{
					++i;
					// all but the first record in the database are preceded by an asterisk
					if ( written_records ) fputc('*', database_file);
					
					database_file_pos = ftell(database_file);
					// a record in the index file consists of eight bytes for the record position in the source file, four bytes for the record position in the database file and 80 bytes for the record name (so one record is 92 bytes long, there are no separators in the index file!
					fwrite(&source_file_pos, sizeof(unsigned long long int), 1, index_file);
					fwrite(&database_file_pos, sizeof(unsigned int), 1, index_file);
					fwrite(ac, 1, index_peptide_name_length_, index_file);
					remove(sequence, "*");
					fputs(sequence.c_str(), database_file);
				}
			}
			// if the sequence is empty
			else
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "emtpy sequence!" , source_filename);
			}
		}
		
		delete(ac);
		
		// close the filestreams
		source_file.close();
		source_file.clear();
		fclose(database_file);
		fclose(index_file);
	}
	
	void InspectFile::generateTrieDB(const std::string& source_filename, const std::string& source_path, const std::string& database_path, std::vector< unsigned int > wanted_records, std::string database_filename, std::string index_filename, bool append, std::string species) throw (Exception::FileNotFound, Exception::ParseError)
	{
		std::string ac_label, sequence_start_label, sequence_end_label, comment_label, species_label;
		
		std::string path_and_file = source_path;
		ensurePathChar(path_and_file);
		path_and_file.append(source_filename);
		getLabels(path_and_file, ac_label, sequence_start_label, sequence_end_label, comment_label, species_label);
		
		if ( sequence_start_label == std::string(1, trie_delimiter_) ) compressTrieDB(source_filename, "", source_path, wanted_records, database_filename, index_filename, database_path, append);
		else generateTrieDB_(source_filename, source_path, database_path, ac_label, sequence_start_label, sequence_end_label, comment_label, species_label, species, wanted_records, database_filename, index_filename, append);
	}
	
	void InspectFile::getLabels(const std::string& source_filename, std::string& ac_label, std::string& sequence_start_label, std::string& sequence_end_label, std::string& comment_label, std::string& species_label) throw (Exception::FileNotFound, Exception::ParseError)
	{
		std::ifstream source_file(source_filename.c_str());
		if ( !source_file ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, source_filename);
		
		String line;
		while ( getline(source_file, line) && (sequence_start_label.empty()) )
		{
			if ( !line.empty() ) line.resize(line.length()-1);
			if ( line.trim().empty() ) continue;
			else if ( line.hasPrefix("#") ) continue;
			
			else if ( line.hasPrefix(">") )
			{
				ac_label = ">";
				sequence_start_label = ">";
				sequence_end_label = ">";
				comment_label = ";";
			}
			else if ( line.hasPrefix("SQ") )
			{
				ac_label = "AC";
				sequence_start_label = "SQ";
				sequence_end_label = "//";
				comment_label = "CC";
				species_label = "OS";
			}
			// a trie database file consists but of one line
			else if ( source_file.eof() )
			{
				ac_label = "*";
				sequence_start_label = "*";
				sequence_end_label = "*";
				comment_label = "#";
			}
		}
		source_file.close();
		source_file.clear();
		
		// if no known start seperator is found
		if ( sequence_start_label.empty() )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "database has unknown file format (neither trie nor FASTA nor swissprot)" , source_filename);
		}
	}

	void InspectFile::getSequences(std::string database_path, const std::string& database_filename_, std::string index_filename, const std::vector< unsigned int > wanted_records, std::vector<std::string>& sequences) throw (Exception::FileNotFound, Exception::ParseError)
	{
		String database_filename(database_filename_);
		ensurePathChar(database_path);
		std::string path_and_file = database_path+database_filename;
		std::ifstream database_file( path_and_file.c_str(), std::ios::in | std::ios::binary );
		if ( !database_file ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
		// if no index file is given and a trie database is used, check whether threre's a corresponding index file
		if ( index_filename.empty() && database_filename.hasSuffix(".trie") )
		{
			index_filename = database_filename.substr(0, database_filename.length()-4);
			index_filename.append("index");
		}
		path_and_file = database_path+index_filename;
		std::ifstream index_file(path_and_file.c_str(), std::ios::in | std::ios::binary);
		if ( !index_file ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, index_filename);
		
		// go through the trie database and extract the sequences according to the record numbers in the record vector
		// determine the length of the index file
		index_file.seekg (0, std::ios::end);
		unsigned int file_length = index_file.tellg();
		index_file.seekg (0, std::ios::beg);
		
		unsigned int index_pos;
		unsigned int database_pos;
		char buffer[11];
		std::stringbuf sequence;
		
		// with index file
		if ( index_file )
		{
			for ( unsigned int i = 0; i < wanted_records.size(); ++i)
			{
				// get the according record in the index file (a record has length 8+4+80 and the database position is stored after 8 bytes in a record and is 4 bytes long)
				index_pos = wanted_records[i] * (index_record_length_) + index_db_record_length_;
				
				// if the file is too short
				if ( file_length < (index_pos+index_trie_record_length_) )
				{
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "index file is too short!", index_filename);
				}
				index_file.seekg(index_pos);
				index_file.readsome(buffer, index_trie_record_length_);
				memcpy(&database_pos, buffer, index_trie_record_length_);
				
				// read the sequence;
				database_file.seekg(database_pos);
				database_file.get(sequence, trie_delimiter_);
				
				// save it in the corresponding protein hit
				sequences.push_back(sequence.str());
				// read out the delimiter
				database_file.get();
				// clear the streambuffer
				sequence.str("");
			}
			
			index_file.close();
		}
		// without index file
		else
		{
			unsigned int i = 0;
			for ( std::vector< unsigned int >::const_iterator it = wanted_records.begin(); it != wanted_records.end(); ++it)
			{
				// skip the unwanted proteins
				while ( i < *it )
				{
					database_file.get(sequence, trie_delimiter_);
					database_file.get();
					sequence.str("");
					++i;
				}
				// read the sequence;
				database_file.get(sequence, trie_delimiter_);
				// save it in the corresponding protein hit
				sequences.push_back(sequence.str());
				// read out the delimiter
				database_file.get();
				// clear the streambuffer
				sequence.str("");
				++i;
			}
		}
		
		database_file.close();
	}
	
	bool InspectFile::setTempDatabaseFilename(const std::string& temp_database_filename)
	{
		if ( temp_database_filename_.empty() )
		{
			temp_database_filename_ = temp_database_filename;
			return true;
		}
		return false;
	}
	const std::string& InspectFile::getTempDatabaseFilename() const
	{
		return temp_database_filename_;
	}
	
	bool InspectFile::setTempIndexFilename(const std::string& temp_index_filename)
	{
		if ( temp_index_filename_.empty() )
		{
			temp_index_filename_ = temp_index_filename;
			return true;
		}
		return false;
	}
	const std::string& InspectFile::getTempIndexFilename() const
	{
		return temp_index_filename_;
	}
	
	bool InspectFile::setSecondTempDatabaseFilename(const std::string& temp_second_database_filename)
	{
		if ( temp_second_database_filename_.empty() )
		{
			temp_second_database_filename_ = temp_second_database_filename;
			return true;
		}
		return false;
	}
	const std::string& InspectFile::getSecondTempDatabaseFilename() const
	{
		return temp_second_database_filename_;
	}
	
	bool InspectFile::setSecondTempIndexFilename(const std::string& temp_second_index_filename)
	{
		if ( temp_second_index_filename_.empty() )
		{
			temp_second_index_filename_ = temp_second_index_filename;
			return true;
		}
		return false;
	}
	const std::string& InspectFile::getSecondTempIndexFilename() const
	{
		return temp_second_index_filename_;
	}
	
	void InspectFile::setUseTempFiles(const bool use_temp_files)
	{
		use_temp_files_ = use_temp_files;
	}
	const bool InspectFile::getUseTempFiles() const
	{
		return use_temp_files_;
	}
	
	void InspectFile::ensurePathChar(std::string& path, char path_char)
	{
		if ( !path.empty() && (std::string("/\\").find(path[path.length()-1], 0) == std::string::npos) ) path.append(1, path_char);
	}
	
	std::string InspectFile::temp_database_filename_ = "";
	std::string InspectFile::temp_index_filename_ = "";
	std::string InspectFile::temp_second_database_filename_ = "";
	std::string InspectFile::temp_second_index_filename_ = "";
	bool InspectFile::use_temp_files_ = 0;
	
	const unsigned int InspectFile::index_peptide_name_length_ = 80;
	const unsigned int InspectFile::index_db_record_length_ = 8;
	const unsigned int InspectFile::index_trie_record_length_ = 4;
	const unsigned int InspectFile::index_record_length_ = index_peptide_name_length_ + index_db_record_length_ + index_trie_record_length_;
	const char InspectFile::trie_delimiter_ = '*';
	const std::string InspectFile::score_type_ = "MQScore";

} // namespace OpenMS
