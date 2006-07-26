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
// $Id: InspectOutfile.C,v 1.0 2006/07/25 13:46:15 martinlangwisch Exp $
// $Author: martinlangwisch $
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/InspectOutfile.h>

namespace OpenMS 
{
  InspectOutfile::InspectOutfile(const std::string& result_filename, const std::string& database_filename_, const std::string& database_path, const double& p_value_threshold, std::string index_filename_)
  	throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument)
    : Outfile()
  {
		// (0) preparations

		// check whether the p_value is correct
		if ( (p_value_threshold < 0) || (p_value_threshold > 1) )
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "p_value_threshold");
		}

		// open the result and database file
		String database_filename(database_filename_);
		String index_filename(index_filename_);
		std::ifstream result_file( result_filename.c_str());
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
		getLabels(path_and_file, buffer, start_seperator, buffer, buffer, buffer);
		database_file.close();
		database_file.clear();
		
		// get the header
		String line;
		//unsigned int number_of_columns;
		std::vector<String> substrings;
		/*if ( getline(result_file, line) )
		{
			if ( !line.empty() ) line.resize(line.length()-1);
			line.split('\t', substrings);
			number_of_columns = substrings.size();
		}
		else // Inspect search was not successfull
		{
			return;
		}
		
		// get the number of some columns
		int spectrum_file_column = -1;
		int scan_column = -1;
		int peptide_column = -1;
		int protein_column = -1;
		int charge_column = -1;
		int MQScore_column = -1;
		int record_number_column = -1;
		
		for ( std::vector< String >::const_iterator iter = substrings.begin(); iter != substrings.end(); ++iter)
		{
			if ( (*iter) == "#SpectrumFile" )
			{
				spectrum_file_column = (iter - substrings.begin());
			}
			else if ( (*iter) == "Scan#" )
			{
				scan_column = (iter - substrings.begin());
			}
			else if ( (*iter) == "Annotation" )
			{
				peptide_column = (iter - substrings.begin());
			}
			else if ( (*iter) == "Protein" )
			{
				protein_column = (iter - substrings.begin());
			}
			else if ( (*iter) == "Charge" )
			{
				charge_column = (iter - substrings.begin());
			}
			else if ( (*iter) == "MQScore" )
			{
				MQScore_column = (iter - substrings.begin());
			}
			else if ( (*iter) == "RecordNumber" )
			{
				record_number_column = (iter - substrings.begin());
			}
		}
		
		// check whether the columns are available in the header
		if ( (spectrum_file_column == -1) || (scan_column == -1) || (peptide_column == -1) || (protein_column == -1) || (charge_column == -1) || (MQScore_column == -1) ||  (record_number_column == -1))
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "at least one of the columns '#SpectrumFile', 'Scan#', 'Annotation', 'Protein', 'Charge', 'MQScore' or  'RecordNumber' is missing!" , result_filename);
		}*/
		
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
		
		// map the protein hits according to their record number in the result file
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
		unsigned int record_number, scan_number, start, end;
		unsigned int rank = 0;
		unsigned int max_record_number = 0;
		unsigned int line_number = 0; // used to report in which line an error occured
		
		while ( getline(result_file, line) )
		{
			if ( !line.empty() ) line.resize(line.length()-1);
			++line_number;
			line.split('\t', substrings);
			
			// check whether the line has enough columns
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
			
			// if the version Inspect.20060620.zip is used, there is a header
			if ( substrings[0] == "#SpectrumFile" ) continue;
			
			// take only those peptides whose p-value is less or equal the given threshold
			if ( atof(substrings[p_value_column].c_str()) > p_value_threshold ) continue;			

			// (1.0) if a new query is found, insert it into the vector and start a new one
			if ( (substrings[spectrum_file_column] != spectrum_file) || ((unsigned int) atoi(substrings[scan_column].c_str()) != scan_number) )
			{
				queries_.push_back(Identification());
				query = &queries_.back();
				spectrum_file = substrings[spectrum_file_column];
				scan_number = atoi(substrings[scan_column].c_str());
				
				query->setCharge(atoi(substrings[charge_column].c_str()));
				query->setPeptideSignificanceThreshold(p_value_threshold);
				rank = 0;
			}
			
			record_number = atoi(substrings[record_number_column].c_str());
			// get accession number and type
			getACAndACType(substrings[protein_column], accession, accession_type);
			
			// (1.1)  if a new protein is found, get all the information and insert it
			if ( rn_position_map.find(record_number) == rn_position_map.end() )
			{
				max_record_number = std::max(max_record_number, record_number);
				
				protein_hit.clear();
				protein_hit.setAccession(accession);
				protein_hit.setAccessionType(accession_type);
// ### score einfach zusammenrechnen?
				//protein_hit.setScore(0.0);
				//protein_hit.setScoreType(score_type_);
				
				rn_position_map.insert(std::make_pair(record_number, protein_hits_.size()));
				protein_hit.setRank(rn_position_map.size());
				protein_hits_.push_back(protein_hit);
			}
			
			// (1.2) get the peptide infos from the new peptide and insert it
			peptide_hit.clear();
			peptide_hit.setScore(atof(substrings[MQ_score_column].c_str()));
			peptide_hit.setScoreType(score_type_);
			start = substrings[peptide_column].find('.')+1;
			end = substrings[peptide_column].find_last_of('.');
			peptide_hit.setSequence(substrings[peptide_column].substr(start, end-start));
			peptide_hit.setRank(++rank);
			peptide_hit.addProteinIndex(datetime, accession);
			
			rank -= updatePeptideHits(peptide_hit, query->getPeptideHits());
			updatePeptideHits(peptide_hit, peptide_hits_);
		} // result file read
		
		// get the sequences of the protein
		// make a vector of the records to get the sequences of the proteins
		std::vector< unsigned int > record_vector;
		for ( std::map< unsigned int, unsigned int >::iterator i = rn_position_map.begin(); i != rn_position_map.end(); ++i)
		{
			record_vector.push_back(i->first);
		}
		// it it's no trie database generate one from the database and use this one
		if ( start_seperator != std::string(1, trie_delimiter_) )
		{
			database_file.close();
			database_file.clear();
			std::string old_database_filename = database_filename;
			database_filename = getTempDatabaseFilename();
			index_filename = getTempIndexFilename();
			
			generateTrieDB(old_database_filename, database_path, database_path, record_vector, database_filename, index_filename);
			
			// change the vector so it contains the positions of the wanted records in the new database (which is just their number, but this way the record vector can be used for 'old', as well as just generated trie databases)
			for ( unsigned int i = 0; i < record_vector.size(); ++i)
			{
				record_vector[i] = i;
			}
		}
		
		std::vector<std::string> sequences;
		getSequences(database_path, database_filename, index_filename, record_vector, sequences);
		
		unsigned int i = 0;
		for ( std::map< unsigned int, unsigned int >::iterator protein_i = rn_position_map.begin(); protein_i != rn_position_map.end(); ++protein_i)
		{
			protein_hits_[protein_i->second].setSequence(sequences[i]);
			++i;
		}
		sequences.clear();
		
		// if there's but one query the protein hits are inserted there instead of a ProteinIdentification object
		if ( queries_.empty() )
		{
			query->setProteinHits(protein_hits_);
			query->setDateTime(datetime);
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

		// now that all information are read
 		ok_ = true;
  }

} //namespace OpenMS
