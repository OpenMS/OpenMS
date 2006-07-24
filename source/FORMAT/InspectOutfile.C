/// -*- Mode: C++; tab-width: 2; -*-
/// vi: set ts=2:
//
/// --------------------------------------------------------------------------
///                   OpenMS Mass Spectrometry Framework
/// --------------------------------------------------------------------------
///  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
///  This library is free software; you can redistribute it and/or
///  modify it under the terms of the GNU Lesser General Public
///  License as published by the Free Software Foundation; either
///  version 2.1 of the License, or (at your option) any later version.
//
///  This library is distributed in the hope that it will be useful,
///  but WITHOUT ANY WARRANTY; without even the implied warranty of
///  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
///  Lesser General Public License for more details.
//
///  You should have received a copy of the GNU Lesser General Public
///  License along with this library; if not, write to the Free Software
///  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
/// --------------------------------------------------------------------------
/// $Id: InspectOutfile.C,v 1.0 2006/07/12 15:58:59 martinlangwisch Exp $
/// $Author: martinlangwisch $
/// $Maintainer: Martin Langwisch $
/// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/InspectOutfile.h>

namespace OpenMS 
{
	InspectOutfile::InspectOutfile()
    : Outfile()
	{}
	
	InspectOutfile::InspectOutfile(const InspectOutfile& source)
    : Outfile(source)
	{}
	
  InspectOutfile::InspectOutfile(const std::string& result_filename, const std::string& database_filename_, const std::string& database_path, std::string index_filename_)
  	throw (Exception::FileNotFound, Exception::ParseError)
    : Outfile()
  {
		/// (0) preparations
		/// open the result and database file
		String database_filename(database_filename_);
		String index_filename(index_filename_);
		std::ifstream result_file(result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}
		std::string path_and_file = database_path;
		ensurePathChar(path_and_file);
		path_and_file.append(database_filename);
		std::ifstream database_file( path_and_file.c_str(), std::ios::in | std::ios::binary );
		if ( !database_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
		}
		
		std::string start_seperator, buffer;
		getSeparators(path_and_file, buffer, start_seperator, buffer, buffer, buffer);
		database_file.close();
		database_file.clear();
		
		/// get the header
		String line;
		std::vector<String> substrings;
		
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
			unused_column,
			p_value_column,
			delta_score_column,
			delta_score_other_column,
			record_number_column,
			DB_file_pos_column,
			spec_file_pos_column
		};
		unsigned int number_of_columns = 16;
		
		/// map the protein hits according to their record number in the result file
		//					record number			position in protein_hits_
		std::map< unsigned int, unsigned int > rn_position_map;
		Identification* query;
		PeptideHit peptide_hit;
		std::vector< PeptideHit >::iterator pep_hit_i;
		DateTime datetime;
		datetime.now();
		ProteinHit protein_hit;
		std::vector< std::pair< String, String > >::iterator prot_hit_i1, prot_hit_i2;
		std::string accession, accession_type, spectrum_file;
		unsigned int record_number, scan_number;
		unsigned int rank = 0;
		unsigned int max_record_number = 0;
		unsigned int line_number = 0; /// used to report in which line an error occured
		
		while ( getline(result_file, line) )
		{
			if ( !line.empty() ) line.resize(line.length()-1);
			++line_number;
			line.split('\t', substrings);
			
			/// check whether the line has enough columns
			if (substrings.size() < number_of_columns )
			{
				char buffer[10];
				sprintf(buffer, "%i", line_number);
				std::string error_message = "wrong number of columns in row ";
				error_message.append(buffer);
				error_message.append("! (");
				sprintf(buffer, "%i", substrings.size());
				error_message.append(buffer);
				error_message.append(" present, should be ");
				sprintf(buffer, "%i", number_of_columns);
				error_message.append(buffer);
				error_message.append(")");
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, error_message.c_str() , result_filename);
			}
			
			/// (1.0) if a new query is found, insert it into the vector and start a new one
			if ( (substrings[spectrum_file_column] != spectrum_file) || ((unsigned int) atoi(substrings[scan_column].c_str()) != scan_number) )
			{
				queries_.push_back(Identification());
				query = &queries_.back();
				spectrum_file = substrings[spectrum_file_column];
				scan_number = atoi(substrings[scan_column].c_str());
				
				query->setCharge(atoi(substrings[charge_column].c_str()));
				query->setDateTime(datetime);
				query->setPeptideSignificanceThreshold(0.0);
				query->setProteinSignificanceThreshold(0.0);
				precursor_retention_times_.push_back(0.0);
				precursor_mz_values_.push_back(0.0);
				rank = 0;
			}
			
			record_number = atoi(substrings[record_number_column].c_str());
			/// get accession number and type
			get_ac_and_ac_type(substrings[protein_column], result_filename, accession, accession_type);
			
			/// (1.1)  if a new protein is found, get all the information and insert it
			if ( rn_position_map.find(record_number) == rn_position_map.end() )
			{
				max_record_number = std::max(max_record_number, record_number);
				
				protein_hit.clear();
				protein_hit.setAccession(accession);
				protein_hit.setAccessionType(accession_type);
/// ### score einfach zusammenrechnen?
				//protein_hit.setScore(0.0);
				//protein_hit.setScoreType("MQScore");
				
				rn_position_map.insert(std::make_pair(record_number, protein_hits_.size()));
				protein_hit.setRank(rn_position_map.size());
				protein_hits_.push_back(protein_hit);
			}
			
			/// (1.2) get the peptide infos from the new peptide and insert it
			peptide_hit.clear();
			peptide_hit.setScore(atof(substrings[MQ_score_column].c_str()));
			peptide_hit.setScoreType("MQScore");
			peptide_hit.setSequence(substrings[peptide_column]);
			peptide_hit.setRank(++rank);
			peptide_hit.addProteinIndex(datetime, accession);
			
			rank -= updatePeptideHits(peptide_hit, query->getPeptideHits());
			updatePeptideHits(peptide_hit, peptide_hits_);
		} /// result file read
		result_file.close();
		
		/// get the sequences of the protein
		/// make a vector of the records to get the sequences of the proteins
		std::vector< unsigned int > record_vector;
		for ( std::map< unsigned int, unsigned int >::iterator i = rn_position_map.begin(); i != rn_position_map.end(); ++i)
		{
			record_vector.push_back(i->first);
		}
		/// it it's no trie database generate one from the database and use this one
		if ( start_seperator != std::string(1, trie_delimiter_) )
		{
			std::string old_database_filename = database_filename;
			database_filename = getTempDatabaseFilename();
			index_filename = getTempIndexFilename();
			
			compressor(old_database_filename, database_path, database_path, record_vector, database_filename, index_filename);
			
			/// change the vector so it contains the positions of the wanted records in the new database (which is just their number, but this way the record vector can be used for 'old', as well as just generated trie databases)
			for ( unsigned int i = 0; i < record_vector.size(); ++i)
			{
				record_vector[i] = i;
			}
		}
		
		/// retrieve the sequences
		std::vector<std::string> sequences;
		getSequences(database_path, database_filename, index_filename, record_vector, sequences);
		
		unsigned int i = 0;
		for ( std::map< unsigned int, unsigned int >::iterator protein_i = rn_position_map.begin(); protein_i != rn_position_map.end(); ++protein_i)
		{
			protein_hits_[protein_i->second].setSequence(sequences[i]);
			++i;
		}
		sequences.clear();
		
		/// if there's but one query the protein hits are inserted there instead of a ProteinIdentification object
		if ( queries_.empty() )
		{
			query->setProteinHits(protein_hits_);
			query->setDateTime(datetime);
			query->setPeptideSignificanceThreshold(0.0);
			query->setProteinSignificanceThreshold(0.0);
		}
		
		protein_ids_.setProteinHits(protein_hits_);
		protein_ids_.setDateTime(datetime);
		
		peptide_hit.clear();
		protein_hit.clear();
		
		rn_position_map.clear();
		record_vector.clear();
		
		curr_peptide_hit_ = peptide_hits_.begin();
		curr_protein_hit_ = protein_hits_.begin();
		curr_query_ = queries_.begin();

		/// now that all information are read
 		ok_ = true;
  }

} //namespace OpenMS
