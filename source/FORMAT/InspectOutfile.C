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

#include <OpenMS/FORMAT/InspectOutfile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>

using namespace std;

namespace OpenMS 
{
	InspectOutfile::InspectOutfile()
	{}
	
	vector< UnsignedInt >
	InspectOutfile::load(
		const string& result_filename,
		vector< IdentificationData >& identifications,
		ProteinIdentification& protein_identification,
		Real p_value_threshold)
	throw (
		Exception::FileNotFound,
		Exception::ParseError,
		Exception::IllegalArgument)
	{
		vector< ProteinHit > protein_hits;

		// check whether the p_value is correct
		if ( (p_value_threshold < 0) || (p_value_threshold > 1) )
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "p_value_threshold");
		}

		// get the header
		enum columns
		{
			spectrum_file_column,
			scan_column,
			peptide_column,
			protein_column,
			charge_column,
			MQ_score_column,
			cut_score_column,
			intense_by_column,
			by_present_column,
			number_of_tryptic_termini_column,
			p_value_column,
			delta_score_column,
			delta_score_other_column,
			record_number_column,
			DB_file_pos_column,
			spec_file_pos_column
		};
		UnsignedInt number_of_columns = 16;
		String line;
		vector<String> substrings;

		//	 record number, position in protein_hits
		map< UnsignedInt, UnsignedInt > rn_position_map;
		Identification* query = NULL;
		PeptideHit peptide_hit;
		vector< PeptideHit >::iterator pep_hit_i;
		DateTime datetime;
		datetime.now();
		ProteinHit protein_hit;
		vector< pair< String, String > >::iterator prot_hit_i1, prot_hit_i2;
		string accession, accession_type, spectrum_file;
		UnsignedInt record_number, scan_number, start, end;
		UnsignedInt rank = 0;
		UnsignedInt peptide_hits;
		UnsignedInt line_number = 0; // used to report in which line an error occured
		UnsignedInt scans = 0;
		vector< UnsignedInt > corrupted_lines;
		// to get the precursor retention time and mz values later, save the filename and the numbers of the scans
		vector< pair< String, vector< UnsignedInt > > > files_and_scan_numbers;
		vector< UnsignedInt >* scan_numbers = NULL;
		
		ifstream result_file(result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}

		while ( getline(result_file, line) )
		{
			++line_number;
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.split('\t', substrings);

			// check whether the line has enough columns
			if ( substrings.size() != number_of_columns )
			{
				if ( line_number == 1 )
				{
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename + " doesn't seem to be an inspect output file!", result_filename);
				}
				corrupted_lines.push_back(line_number);
				continue;
			}

			// the is a header which is skipped
			if ( substrings[0] == "#SpectrumFile" ) continue;

			// take only those peptides whose p-value is less or equal the given threshold
			if ( substrings[p_value_column].toFloat() > p_value_threshold ) continue;

			// get accession number and type
			getACAndACType(substrings[protein_column], accession, accession_type);

			record_number = substrings[record_number_column].toInt();
			
			// if a new protein is found, get the rank and insert it
			if ( rn_position_map.find(record_number) == rn_position_map.end() )
			{
				// map the record number to the size in protein hits
				rn_position_map[record_number] = protein_hits.size();

				protein_hit.clear();
				protein_hit.setRank(rn_position_map.size());
				protein_hit.setAccession(accession);
				protein_hit.setAccessionType(accession_type);
				//protein_hit.setScore(0.0);
				//protein_hit.setScoreType(score_type_);
				protein_hits.push_back(protein_hit);
			}
			
			// if a new query is found, insert it into the vector
			// the first time, the condition is always fullfilled because spectrum_file is ""
			if ( (substrings[spectrum_file_column] != spectrum_file) || ((UnsignedInt) substrings[scan_column].toInt() != scan_number) )
			{
				identifications.push_back(IdentificationData());
				query = &(identifications.back().id);
				
				query->setPeptideSignificanceThreshold(p_value_threshold);
				query->setDateTime(datetime);
				rank = 0;
				
				if ( substrings[spectrum_file_column] != spectrum_file )
				{
					files_and_scan_numbers.push_back(make_pair(substrings[spectrum_file_column], vector< UnsignedInt >()));
					scan_numbers = &(files_and_scan_numbers.back().second);
				}
				
				spectrum_file = substrings[spectrum_file_column];
				scan_number = substrings[scan_column].toInt();
				
				scan_numbers->push_back(scan_number);
				++scans;
			}
			
			// get the peptide infos from the new peptide and insert it
			peptide_hit.clear();
			peptide_hit.setCharge(substrings[charge_column].toInt());
			peptide_hit.setScore(substrings[MQ_score_column].toFloat());
			peptide_hit.setScoreType("Inspect");
			start = substrings[peptide_column].find('.')+1;
			end = substrings[peptide_column].find_last_of('.');
			peptide_hit.setSequence(substrings[peptide_column].substr(start, end-start));
			peptide_hit.setRank(++rank);
			peptide_hit.addProteinIndex(datetime, accession);
			
			peptide_hits = query->getPeptideHits().size();
			updatePeptideHits(peptide_hit, query->getPeptideHits());
			rank -= ( peptide_hits == query->getPeptideHits().size() );
		}
		// result file read
		result_file.close();
		result_file.clear();
		
//		// search the sequence of the proteins
//		if ( !protein_hits.empty() )
//		{
//			vector< String > sequences;
//			getSequences(database_filename, rn_position_map, sequences);
//
//			// set the retrieved sequences
//			vector< String >::const_iterator s_i = sequences.begin();
//			for ( map< UnsignedInt, UnsignedInt >::const_iterator rn_i = rn_position_map.begin(); rn_i != rn_position_map.end(); ++rn_i, ++s_i )
//			{
//				protein_hits[rn_i->second].setSequence(*s_i);
//			}
//			sequences.clear();
//		}
		
		// get the precursor retention times and mz values
		getPrecursorRTandMZ(files_and_scan_numbers, identifications);

		// if there's but one query the protein hits are inserted there instead of a ProteinIdentification object
		if ( identifications.size() == 1 )
		{
			query->setProteinHits(protein_hits);
			query->setDateTime(datetime);
		}
		
		protein_identification.setProteinHits(protein_hits);
		protein_identification.setDateTime(datetime);
		
		protein_hits.clear();
		peptide_hit.clear();
		protein_hit.clear();
		
		rn_position_map.clear();
		return corrupted_lines;
  }

	vector< UnsignedInt >
	InspectOutfile::getSequences(
		const string& database_filename,
		const map< UnsignedInt, UnsignedInt >& wanted_records, // < record number, number of protein in a vector >
		vector< String >& sequences)
	throw (
		Exception::FileNotFound)
	{
		ifstream database(database_filename.c_str());
		if ( !database )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
		}

		vector< UnsignedInt > not_found;
		UnsignedInt seen_records = 0;
		stringbuf sequence;
		database.seekg(0, ios::end);
		streampos sp = database.tellg();
		database.seekg(0, ios::beg);

		for ( map< UnsignedInt, UnsignedInt >::const_iterator wr_i = wanted_records.begin(); wr_i !=  wanted_records.end(); ++wr_i )
		{
			for ( ; seen_records < wr_i->first; ++seen_records )
			{
				database.ignore(sp, trie_delimiter_);
			}
			database.get(sequence, trie_delimiter_);
			sequences.push_back(sequence.str());
			if ( sequences.back().empty() ) not_found.push_back(wr_i->first);
			sequence.str("");
		}

		// close the filestreams
		database.close();
		database.clear();
		
		return not_found;
	}
	
	void
	InspectOutfile::getACAndACType(
		String line,
		string& accession,
		string& accession_type)
	{
		accession.clear();
		accession_type.clear();
		pair<string, string> p;
		// if it's a FASTA line
		if ( line.hasPrefix(">") ) line.erase(0,1);
		if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
		line.trim();
		
		// if it's a swissprot accession
		if ( line.hasPrefix("tr") || line.hasPrefix("sp") )
		{
			accession = line.substr(3, line.find('|', 3)-3);
			accession_type = "SwissProt";
		}
		else if ( line.hasPrefix("gi") )
		{
			size_t snd = line.find('|', 3);
			size_t third(0);
			if ( snd != string::npos )
			{
				third = line.find('|', ++snd) + 1;
				
				accession = line.substr(third, line.find('|', third)-third);
				accession_type = line.substr(snd, third-1-snd);
			}
			if ( accession_type == "gb" ) accession_type = "GenBank";
			else if ( accession_type == "emb" ) accession_type = "EMBL";
			else if ( accession_type == "dbj" ) accession_type = "DDBJ";
			else if ( accession_type == "ref" ) accession_type = "NCBI";
			else if ( (accession_type == "sp") || (accession_type == "tr") ) accession_type = "SwissProt";
			else if ( accession_type == "gnl" )
			{
				accession_type = accession;
				snd = line.find('|', third);
				third = line.find('|', ++snd);
				if ( third != string::npos ) accession = line.substr(snd, third-snd);
				else
				{
					third = line.find(' ', snd);
					if ( third != string::npos ) accession = line.substr(snd, third-snd);
					else accession = line.substr(snd);
				}
			}
			else
			{
				accession_type = "gi";
				if ( snd != string::npos ) accession = line.substr(3, snd-4);
				else
				{
					if ( snd == string::npos ) snd = line.find(' ', 3);
					if ( snd != string::npos ) accession = line.substr(3, snd-3);
					else accession = line.substr(3);
				}
			}
		}
		else if ( line.hasPrefix("ref") )
		{
			accession = line.substr(4, line.find('|', 4) - 4);
			accession_type = "NCBI";
		}
// 		// if it's a swissprot line
// 		else if ( line.hasPrefix("AC") )
// 		{
// 			line.erase(0,2);
// 			accession = line.trim();
// 			accession_type = "SwissProt";
// 		}
		else if ( line.hasPrefix("gnl") )
		{
			line.erase(0,3);
			accession_type = line.substr(0, line.find('|', 0));
			accession = line.substr(accession_type.length()+1);
		}
		else if ( line.hasPrefix("lcl") )
		{
			line.erase(0,4);
			accession_type = "lcl";
			accession = line;
		}
		else
		{
			size_t pos1 = line.find('(', 0);
			size_t pos2;
			if ( pos1 != string::npos )
			{
				pos2 = line.find(')', ++pos1);
				if ( pos2 != string::npos )
				{
					accession = line.substr(pos1, pos2 - pos1);
					if ( (accession.size() == 6) && (String("OPQ").find(accession[0], 0) != string::npos) ) accession_type = "SwissProt";
					else accession.clear();
				}
			}
			if ( accession.empty() )
			{
				pos1 = line.find('|');
				accession = line.substr(0, pos1);
				if ( (accession.size() == 6) && (String("OPQ").find(accession[0], 0) != string::npos) ) accession_type = "SwissProt";
				else
				{
					pos1 = line.find(' ');
					accession = line.substr(0, pos1);
					if ( (accession.size() == 6) && (String("OPQ").find(accession[0], 0) != string::npos) ) accession_type = "SwissProt";
					else
					{
						accession = line.substr(0, 6);
						if ( String("OPQ").find(accession[0], 0) != string::npos ) accession_type = "SwissProt";
						else accession.clear();
					}
				}
			}
		}
		if ( accession.empty() )
		{
			accession = line.trim();
			accession_type = "unknown";
		}
	}
	
	bool
	InspectOutfile::updatePeptideHits(
		PeptideHit& peptide_hit,
		vector< PeptideHit >& peptide_hits)
	{
		// a peptide hit may only be inserted if it's score type matches the one of the existing hits
		if ( (peptide_hits.empty()) || (peptide_hits[0].getScoreType() == peptide_hit.getScoreType()) )
		{
			// search for the peptide hit
			vector< PeptideHit >::iterator pep_hit_i;
			for ( pep_hit_i = peptide_hits.begin(); pep_hit_i != peptide_hits.end(); ++pep_hit_i)
			{
				if ( (pep_hit_i->getSequence() == peptide_hit.getSequence()) && (pep_hit_i->getScore() == peptide_hit.getScore()) ) break;
			}
			// if a new peptide is found, insert it
			if ( pep_hit_i == peptide_hits.end() )
			{
				peptide_hits.push_back(peptide_hit);
				return true;
			}
			// if the peptide has already been inserted, insert additional protein hits
			else
			{
				vector< pair< String, String > >::iterator prot_hit_i1, prot_hit_i2;
				// remove all protein hits from the peptide that are already in the list
				for ( prot_hit_i1 = peptide_hit.getProteinIndices().begin(); prot_hit_i1 != peptide_hit.getProteinIndices().end(); )
				{
					prot_hit_i2 = find(pep_hit_i->getProteinIndices().begin(), pep_hit_i->getProteinIndices().end(), *prot_hit_i1);
					if ( prot_hit_i2 != pep_hit_i->getProteinIndices().end() ) peptide_hit.getProteinIndices().erase(prot_hit_i1);
					else ++prot_hit_i1;
				}
				// add the additional protein hits
				for ( prot_hit_i2 = peptide_hit.getProteinIndices().begin(); prot_hit_i2 != peptide_hit.getProteinIndices().end(); ++prot_hit_i2 )
				{
					pep_hit_i->addProteinIndex(*prot_hit_i2);
				}
				return true;
			}
		}
		return false;
	}
	
	void
	InspectOutfile::getPrecursorRTandMZ(
		const vector< pair< String, vector< UnsignedInt > > >& files_and_scan_numbers,
		vector< IdentificationData >& ids)
	throw(
		Exception::ParseError)
	{
		MSExperiment<> experiment;
		String type;
		
		UnsignedInt pos = 0;
		for ( vector< pair< String, vector< UnsignedInt > > >::const_iterator fs_i = files_and_scan_numbers.begin(); fs_i != files_and_scan_numbers.end(); ++fs_i )
		{
			getExperiment(experiment, type, fs_i->first); // may throw an exception if the filetype could not be determined
			
			if ( experiment.size() < fs_i->second.back() )
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Not enought scans in file! (" + String(experiment.size()) + " available, should be " + String(fs_i->second.back()) + ")", fs_i->first);
			}
			
			for ( vector< UnsignedInt >::const_iterator scan_i = fs_i->second.begin(); scan_i != fs_i->second.end(); ++scan_i )
			{
				ids[pos].mz = experiment[*scan_i - 1].getPrecursorPeak().getPosition()[0];
				ids[pos++].rt = experiment[*scan_i - 1].getRetentionTime();
			}
		}
	}
	
	void
	InspectOutfile::compressTrieDB(
		const string& database_filename,
		const string& index_filename,
		vector< UnsignedInt >& wanted_records,
		const string& snd_database_filename,
		const string& snd_index_filename,
		bool append)
	throw (
		Exception::FileNotFound,
		Exception::ParseError,
		Exception::UnableToCreateFile)
	{
		if ( (database_filename == snd_database_filename) || (index_filename == snd_index_filename) )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Same filename can not be used for original and second database!", index_filename);
		}
		
		ifstream database( database_filename.c_str());
		if ( !database )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
		}
		
		ifstream index(index_filename.c_str(), ios::in | ios::binary);
		if ( !index )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, index_filename);
		}
		
		// determine the length of the index file
		index.seekg(0, ios::end);
		streampos index_length = index.tellg();
		index.seekg(0, ios::beg);
		bool empty_records = wanted_records.empty();
		if ( wanted_records.empty() )
		{
			for ( UnsignedInt i = 0; i < index_length / record_length_; ++i ) wanted_records.push_back(i);
		}
		
		// take the wanted records, copy their sequences to the new db and write the index file accordingly
		ofstream snd_database;
		if ( append ) snd_database.open(snd_database_filename.c_str(), std::ios::out | std::ios::app);
		else snd_database.open(snd_database_filename.c_str(), std::ios::out | std::ios::trunc);
		
		ofstream snd_index;
		if ( append ) snd_index.open(snd_index_filename.c_str(), std::ios::out | std::ios::binary | std::ios::app);
		else snd_index.open(snd_index_filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
		
		char* index_record = new char[record_length_]; // to copy one record from the index file
		UnsignedInt database_pos, snd_database_pos; // their sizes HAVE TO BE 4 bytes
		stringbuf sequence;
		streampos index_pos;
		
		for ( vector< UnsignedInt >::const_iterator wr_i = wanted_records.begin(); wr_i != wanted_records.end(); ++wr_i )
		{
			// get the according record in the index file
			if ( index_length < (*wr_i + 1) * record_length_ ) // if the file is too short
			{
				delete(index_record);
				database.close();
				index.close();
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "index file is too short!", index_filename);
			}
			index.seekg((*wr_i) * record_length_);
			index.read(index_record, record_length_);
			
			// all but the first sequence are preceded by an asterisk
			if ( append ) snd_database.put(trie_delimiter_);
			append = true;
			
			// go to the beginning of the sequence
			memcpy(&database_pos, index_record + db_pos_length_, trie_db_pos_length_);
			database.seekg(database_pos);
			
			// store the corresponding index for the second database
			snd_database_pos = snd_database.tellp(); // get the position in the second database
			memcpy(index_record + db_pos_length_, &snd_database_pos, trie_db_pos_length_); // and copy to its place in the index record
			snd_index.write((char*) index_record, record_length_); // because only the trie-db position changed, not the position in the original database, nor the protein name
			
			// store the sequence
			database.get(sequence, trie_delimiter_);
			snd_database << sequence.str();
			sequence.str("");
		}
		
		if ( empty_records ) wanted_records.clear();
		delete(index_record);
		index.close();
		database.close();
		snd_database.close();
		snd_index.close();
	}
	
	void
	InspectOutfile::generateTrieDB(
		const string& source_database_filename,
		const string& database_filename,
		const string& index_filename,
		bool append,
		const string species)
	throw (
		Exception::FileNotFound,
		Exception::ParseError,
		Exception::UnableToCreateFile)
	{
		ifstream source_database(source_database_filename.c_str());
		if ( !source_database )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, source_database_filename);
		}
		
		// get the labels
		string ac_label, sequence_start_label, sequence_end_label, comment_label, species_label;
		getLabels(source_database_filename, ac_label, sequence_start_label, sequence_end_label, comment_label, species_label);

		ofstream database;
		if ( append ) database.open(database_filename.c_str(), ios::app | ios::out );
		else database.open(database_filename.c_str());
		if ( !database )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
		}
		ofstream index;
		if ( append ) index.open(index_filename.c_str(), ios::app | ios::out | ios::binary );
		else index.open(index_filename.c_str(), ios::out | ios::binary );
		if ( !index )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, index_filename);
		}
		
		// using flags to mark what has already been read
		// the flags
		unsigned char ac_flag = 1;
		unsigned char species_flag = !species.empty()*2; // if no species is given, take all proteins
		unsigned char sequence_flag = 4;
		// the value
		unsigned char record_flags = 0;
		
		string::size_type pos; // the position in a line
		unsigned long long source_database_pos = source_database.tellg(); // the start of a protein in the source database
		unsigned long long source_database_pos_buffer = 0; // because you don't know whether a new protein starts unless the line is read, the actual position is buffered before any new getline
		UnsignedInt database_pos;
		String line, sequence, protein_name;
		char* record = new char[record_length_]; // a record in the index file
		char* protein_name_pos = record + db_pos_length_ + trie_db_pos_length_;
		
		while ( getline(source_database, line) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.trim();
			
			// empty and comment lines are skipped
			if ( line.empty() || line.hasPrefix(comment_label) )
			{
				source_database_pos_buffer = source_database.tellg();
				continue;
			}
			
			// read the sequence if the accession and the species have been read already
			if ( record_flags == (ac_flag | species_flag | sequence_flag) )
			{
				if ( !line.hasPrefix(sequence_end_label) ) // if it is still the same protein, append the sequence
				{
					line.trim(); // erase all whitespaces from the sequence
					line.remove(trie_delimiter_);
					// save this part of the sequence
					sequence.append(line);
				}
				else // if a new protein is found, write down the old one
				{
					// if the sequence is not empty, the record has the correct form
					if ( !sequence.empty() )
					{
						// all but the first record in the database are preceded by an asterisk (if in append mode an asterisk has to be put at any time)
						if ( append ) database.put('*');
						database_pos = database.tellp();
						
						// write the record
						memcpy(record, &source_database_pos, db_pos_length_); // source database position
						memcpy(record + db_pos_length_, &database_pos, trie_db_pos_length_); // database position
						index.write(record, record_length_);
						// protein name / accession has already been written
						database << sequence;
						source_database_pos = source_database_pos_buffer; // the position of the start of the new protein
						append = true;
					}
					sequence.clear();
					
					// set back the record flags for a new record
					record_flags = 0;
				}
			}
			
			// if not reading the sequence
			if ( !(record_flags & sequence_flag) )
			{
				if ( line.hasPrefix(ac_label) )
				{
					pos = ac_label.length(); // find the beginning of the accession

					while ( (line.length() > pos) && (line[pos] < 33) ) ++pos; // discard the whitespaces after the label
					if ( pos != line.length() ) // if no accession is found, skip this protein
					{
						memset(protein_name_pos, 0, protein_name_length_); // clear the protein name
						// read at most protein_name_length_ characters from the record name and write them to the record
						protein_name = line.substr(pos, protein_name_length_);
						protein_name.substitute('>', '}');
						memcpy(protein_name_pos, protein_name.c_str(), protein_name.length());

						record_flags |= ac_flag; // set the ac flag
					}
					else record_flags = 0;
				}
				// if a species line is found and an accession has already been found, check whether this record is from the wanted species, if not, skip it
				if ( species_flag && line.hasPrefix(species_label) && (record_flags == ac_flag) )
				{
					pos = species_label.length();
					if ( line.find(species, pos) != string::npos ) record_flags |= species_flag;
					else record_flags = 0;
				}
				// if the beginning of the sequence is found and accession and correct species have been found
				if ( line.hasPrefix(sequence_start_label) && ((record_flags & (ac_flag | species_flag)) == (ac_flag | species_flag)) ) record_flags |= sequence_flag;
			}
			source_database_pos_buffer = source_database.tellg();
		}
		// source file read
		source_database.close();
		source_database.clear();
		
		// if the last record has no sequence end label, the sequence has to be appended nevertheless (e.g. FASTA)
		if ( record_flags == (ac_flag | species_flag | sequence_flag) && !sequence.empty() )
		{
			// all but the first record in the database are preceded by an asterisk (if in append mode an asterisk has to be put at any time)
			if ( append ) database.put('*');
			database_pos = database.tellp();
			
			// write the record
			memcpy(record, &source_database_pos, db_pos_length_); // source database position
			memcpy(record + db_pos_length_, &database_pos, trie_db_pos_length_); // database position
			index.write(record, record_length_);
			// protein name / accession has already been written
			database << sequence;
			append = true;
		}
		
		delete(record);
		
		// close the filestreams
		database.close();
		database.clear();
		index.close();
		index.clear();
	}
	
	void
	InspectOutfile::getLabels(
		const string& source_database_filename,
		string& ac_label,
		string& sequence_start_label,
		string& sequence_end_label,
		string& comment_label,
		string& species_label)
	throw (
		Exception::FileNotFound,
		Exception::ParseError)
	{
		ifstream source_database(source_database_filename.c_str());
		if ( !source_database ) throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, source_database_filename);
		
		String line;
		while ( getline(source_database, line) && (sequence_start_label.empty()) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			if ( line.trim().empty() ) continue;
			
			else if ( line.hasPrefix(">") )
			{
				ac_label = ">";
				sequence_start_label = ">";
				sequence_end_label = ">";
				comment_label = ";";
				species_label = ">";
			}
			else if ( line.hasPrefix("SQ") )
			{
				ac_label = "AC";
				sequence_start_label = "SQ";
				sequence_end_label = "//";
				comment_label = "CC";
				species_label = "OS";
			}
		}
		source_database.close();
		source_database.clear();
		
		// if no known start seperator is found
		if ( sequence_start_label.empty() )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "database has unknown file format (neither trie nor FASTA nor swissprot)" , source_database_filename);
		}
	}

	vector< UnsignedInt >
	InspectOutfile::getWantedRecords(
		const string& result_filename,
		Real p_value_threshold)
	throw (
		Exception::FileNotFound,
		Exception::ParseError,
		Exception::IllegalArgument)
	{
		// check whether the p_value is correct
		if ( (p_value_threshold < 0) || (p_value_threshold > 1) )
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "p_value_threshold");
		}
		
		// get the header
		enum columns
		{
			spectrum_file_column,
			scan_column,
			peptide_column,
			protein_column,
			charge_column,
			MQ_score_column,
			cut_score_column,
			intense_by_column,
			by_present_column,
			number_of_tryptic_termini_column,
			p_value_column,
			delta_score_column,
			delta_score_other_column,
			record_number_column,
			DB_file_pos_column,
			spec_file_pos_column
		};
		UnsignedInt number_of_columns = 16;
		String line;
		vector<String> substrings;
		
		set< UnsignedInt > wanted_records_set;
		vector< UnsignedInt > wanted_records;
		
		ifstream result_file(result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}
		
		while ( getline(result_file, line) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.split('\t', substrings);

			// if the version Inspect.20060620.zip is used, there is a header which is skipped
			if ( substrings[0] == "#SpectrumFile" ) continue;

			// check whether the line has enough columns
			if ( substrings.size() != number_of_columns ) continue;
			
			// take only those peptides whose p-value is less or equal the given threshold
			if ( substrings[p_value_column].toFloat() > p_value_threshold ) continue;
			
			wanted_records_set.insert(substrings[record_number_column].toInt());
		}
		
		result_file.close();
		result_file.clear();
		
		for ( set< UnsignedInt >::const_iterator rn_i = wanted_records_set.begin(); rn_i != wanted_records_set.end(); ++rn_i )
		{
			wanted_records.push_back(*rn_i);
		}
		
		return wanted_records;
	}

	template< typename PeakT >
	void
	InspectOutfile::getExperiment(
		MSExperiment< PeakT >& exp,
		String& type,
		const String& in_filename)
	throw(
		Exception::ParseError)
	{
		type.clear();
		exp.reset();
		//input file type
		FileHandler fh;
		FileHandler::Type in_type = fh.getTypeByContent(in_filename);
		if (in_type==FileHandler::UNKNOWN)
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Could not determine type of the file. Aborting!" , in_filename);
		}
		type = fh.typeToName(in_type);
		fh.loadExperiment(in_filename, exp, in_type);
	}
	
	const UnsignedInt InspectOutfile::db_pos_length_ = 8;
	const UnsignedInt InspectOutfile::trie_db_pos_length_ = 4;
	const UnsignedInt InspectOutfile::protein_name_length_ = 80;
	const UnsignedInt InspectOutfile::record_length_ = db_pos_length_ + trie_db_pos_length_ + protein_name_length_;
	const char InspectOutfile::trie_delimiter_ = '*';
	const string InspectOutfile::score_type_ = "Inspect";
	
} //namespace OpenMS
