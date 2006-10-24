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

#include <OpenMS/FORMAT/InspectOutfile.h>

namespace OpenMS 
{
	InspectOutfile::InspectOutfile()
	{}
	
  void InspectOutfile::load(const std::string& result_filename, std::vector< Identification >&	identifications, ProteinIdentification&	protein_identification, std::vector< float >& 	precursor_retention_times, std::vector< float >& precursor_mz_values, const double& p_value_threshold, const double& score_value_threshold, const std::string& database_filename, const std::string& database_path, const std::string& sequence_filename, std::string index_filename) throw (Exception::FileNotFound, Exception::ParseError, Exception::IllegalArgument)
  {
		// (0) preparations
		std::vector< PeptideHit > peptide_hits;
		std::vector< ProteinHit > protein_hits;

		// check whether the p_value is correct
		if ( (p_value_threshold < 0) || (p_value_threshold > 1) )
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "p_value_threshold");
		}

		// open the result and database file(s)
		std::ifstream result_file(result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}
		
		// get the header
		String line;
		//unsigned int number_of_columns;
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
		bool from_fasta = (!sequence_filename.empty());
		bool missing_column;
		if ( (!database_filename.empty()) && from_fasta )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "trie AND FASTA database given (only one of both is allowed).", result_filename);
		}
		
		// map the protein hits according to their record number in the result file
		//					record number			position in protein_hits
		std::map< std::pair< bool, int >, int > rn_position_map;
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
		unsigned int line_number = 0; // used to report in which line an error occured
		double p_value;
		
		// workaround for a bug in inspect
		// if there is at least one line with a missing protein column, the record numbers are one too high
		bool false_record_number = false;
		
		while ( getline(result_file, line) && !false_record_number )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.split('\t', substrings);
			
			// check whether the line has enough columns (a line from a fasta-db does not include the protein name)
			if ( substrings.size() < number_of_columns - 1 ) continue;
			false_record_number = ( line_number && (substrings.size() != line_number) );
			line_number = substrings.size();
		}
		result_file.close();
		result_file.clear();
		result_file.open(result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}
		line_number = 0;
		
		// get all the proteins and make a vector of those who come from a FASTA file and those who come from a trie database
		std::vector< unsigned int > FASTA_proteins, trie_proteins;
		
		while ( getline(result_file, line) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			++line_number;
			line.split('\t', substrings);
			
			// check whether the line has enough columns (a line from a fasta-db does not include the protein name)
			missing_column = ( substrings.size() == number_of_columns-1 );
			if ( substrings.size() < number_of_columns - missing_column ) continue;
			
			// if the version Inspect.20060620.zip is used, there is a header
			if ( substrings[0] == "#SpectrumFile" ) continue;
			
			// take only those peptides whose p-value is less or equal the given threshold (if no p-value is found, take the protein if it's MQ score is above the given threshold
			// workaround for a bug in inspect. sometimes, when the protein column is missing, there's an extra zeros-column before the p-value column
			if ( missing_column && (substrings[p_value_column - missing_column].length() >= 5) && (substrings[p_value_column - missing_column].substr(0,5) == "0.000") )
			{
				p_value = atof(substrings[p_value_column - missing_column].c_str());
				
				if ( (p_value >= 0) && (p_value <= 1) && (p_value > p_value_threshold) ) continue;
			}
			else if ( ((substrings[p_value_column - missing_column] == "nan") && (atof(substrings[MQ_score_column - missing_column].c_str()) < score_value_threshold)) || (atof(substrings[p_value_column - missing_column].c_str()) > p_value_threshold) ) continue;
			
			// if there's a missing column, the record number is one too high
			record_number = atoi(substrings[record_number_column - missing_column].c_str()) - false_record_number;
			
			// (1.1)  if a new protein is found, get the rank and insert it
			if ( rn_position_map.find(std::make_pair(from_fasta, record_number)) == rn_position_map.end() )
			{
				protein_hit.clear();
				//protein_hit.setScore(0.0);
				//protein_hit.setScoreType(score_type_);
				rn_position_map.insert(std::make_pair(std::make_pair(from_fasta, record_number), protein_hits.size()));
				
				// get accession number and type
				if ( !from_fasta )
				{
					getACAndACType(substrings[protein_column], accession, accession_type);
					protein_hit.setAccession(accession);
					protein_hit.setAccessionType(accession_type);
				}
				
				protein_hit.setRank(rn_position_map.size());
				protein_hits.push_back(protein_hit);
			}
		} // result file read
		result_file.close();
		result_file.clear();
		
		for ( std::map< std::pair< bool, int >, int >::const_iterator i = rn_position_map.begin(); i != rn_position_map.end(); ++i )
		{
			if ( i->first.first ) FASTA_proteins.push_back(i->first.second);
			else trie_proteins.push_back(i->first.second);
		}
		
		// search the sequence, accession and accession type of the proteins from a FASTA file
		if ( !FASTA_proteins.empty() )
		{
			std::vector< std::vector< String > > protein_info;
			std::string ac_label, sequence_start_label, sequence_end_label, comment_label, species_label;
			
			getLabels(sequence_filename, ac_label, sequence_start_label, sequence_end_label, comment_label, species_label);
			
			getSequenceAndACandACType(sequence_filename, FASTA_proteins, protein_info, ac_label, sequence_start_label, sequence_end_label, comment_label, species_label);
			
			// set the retrieved protein info
			std::vector< std::vector< String > >::const_iterator p_i = protein_info.begin();
			for ( std::vector< unsigned int >::const_iterator i = FASTA_proteins.begin(); i != FASTA_proteins.end(); ++i, ++p_i )
			{
				protein_hits[rn_position_map[std::make_pair(true, *i)]].setSequence((*p_i)[0]);
				getACAndACType((*p_i)[1], accession, accession_type);
				protein_hits[rn_position_map[std::make_pair(true, *i)]].setAccession(accession);
				protein_hits[rn_position_map[std::make_pair(true, *i)]].setAccessionType(accession_type);
			}
			protein_info.clear();
			FASTA_proteins.clear();
		}
		
		result_file.open( result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}
		
		// to get the precursor retention time and mz values
		//                      filename     scan numbers
		std::vector< std::pair< String, std::vector< unsigned int > > > files_and_scan_numbers;
		std::vector< unsigned int >* scan_numbers;
		while ( getline(result_file, line) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			++line_number;
			line.split('\t', substrings);
			
			// check whether the line has enough columns (a line from a fasta-db does not include the protein name)
			missing_column = ( substrings.size() == number_of_columns-1 );
			if ( substrings.size() < number_of_columns - missing_column ) continue;
			
			// if the version Inspect.20060620.zip is used, there is a header
			if ( substrings[0] == "#SpectrumFile" ) continue;
			
			// workaround for a bug in inspect. sometimes, when the protein column is missing, there's an extra zeros-column before the p-value column
			if ( missing_column && (substrings[p_value_column - missing_column].length() >= 5) && (substrings[p_value_column - missing_column].substr(0,5) == "0.000") )
			{
				p_value = atof(substrings[p_value_column - missing_column].c_str());
				
				if ( (p_value >= 0) && (p_value <= 1) && (p_value > p_value_threshold) ) continue;
			}
			else if ( ((substrings[p_value_column - missing_column] == "nan") && (atof(substrings[MQ_score_column - missing_column].c_str()) < score_value_threshold)) || (atof(substrings[p_value_column - missing_column].c_str()) > p_value_threshold) ) continue;
			
			// (1.0) if a new query is found, insert it into the vector
			if ( (substrings[spectrum_file_column] != spectrum_file) || ((unsigned int) atoi(substrings[scan_column].c_str()) != scan_number) )
			{
				identifications.push_back(Identification());
				query = &identifications.back();
				
				query->setCharge(atoi(substrings[charge_column - missing_column].c_str()));
				query->setPeptideSignificanceThreshold(p_value_threshold);
				query->setDateTime(datetime);
				rank = 0;
				
				if ( substrings[spectrum_file_column] != spectrum_file )
				{
					files_and_scan_numbers.push_back(std::make_pair(substrings[spectrum_file_column], std::vector< unsigned int >()));
					scan_numbers = &(files_and_scan_numbers.back().second);
				}
				
				spectrum_file = substrings[spectrum_file_column];
				scan_number = atoi(substrings[scan_column].c_str());
				
				scan_numbers->push_back(scan_number);
				
				precursor_retention_times.push_back(0);
				precursor_mz_values.push_back(0);
			}
			
			record_number = atoi(substrings[record_number_column - missing_column].c_str()) - false_record_number;
			
			// (1.2) get the peptide infos from the new peptide and insert it
			peptide_hit.clear();
			peptide_hit.setScore(atof(substrings[MQ_score_column - missing_column].c_str()));
			peptide_hit.setScoreType(score_type_);
			start = substrings[peptide_column].find('.')+1;
			end = substrings[peptide_column].find_last_of('.');
			peptide_hit.setSequence(substrings[peptide_column].substr(start, end-start));
			peptide_hit.setRank(++rank);
			peptide_hit.addProteinIndex(datetime, protein_hits[rn_position_map[std::make_pair(from_fasta, record_number)]].getAccession());
			
			rank -= updatePeptideHits(peptide_hit, query->getPeptideHits());
			updatePeptideHits(peptide_hit, peptide_hits);
		} // result file read
		result_file.close();
		
		// get the sequences of the trie proteins
		if ( !trie_proteins.empty() )
		{
			std::vector< std::string > sequences;
			getSequences(database_path, database_filename, index_filename, trie_proteins, sequences);
			
			std::vector< std::string >::const_iterator p_i = sequences.begin();
			for ( std::vector< unsigned int >::const_iterator i = trie_proteins.begin(); i != trie_proteins.end(); ++i, ++p_i)
			{
				protein_hits[rn_position_map[std::make_pair(true, *i)]].setSequence(*p_i);
			}
			sequences.clear();
			trie_proteins.clear();
		}
		
		// get the precursor retention times and mz values
		getPrecursorRTandMZ(files_and_scan_numbers, precursor_retention_times, precursor_mz_values);
		
		// if there's but one query the protein hits are inserted there instead of a ProteinIdentification object
		if ( identifications.size() == 1 )
		{
			query->setProteinHits(protein_hits);
			query->setDateTime(datetime);
		}
		
		protein_identification.setProteinHits(protein_hits);
		protein_identification.setDateTime(datetime);
		
		peptide_hits.clear();
		protein_hits.clear();
		peptide_hit.clear();
		protein_hit.clear();
		
		rn_position_map.clear();
  }

	void InspectOutfile::getACAndACType(String line, std::string& accession, std::string& accession_type) throw (Exception::ParseError)
	{
		std::pair<std::string, std::string> p;
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
			unsigned int snd = line.find('|', 3);
			unsigned int third = line.find('|', ++snd);
			++third;
			accession = line.substr(third, line.find('|', third)-third);
			accession_type = line.substr(snd, third-1-snd);
			if ( accession_type == "gb" ) accession_type = "GenBank";
			else if ( accession_type == "emb" ) accession_type = "EMBL";
			else if ( accession_type == "dbj" ) accession_type = "DDBJ";
			else if ( accession_type == "ref" ) accession_type = "NCBI";
			else if ( (accession_type == "sp") || (accession_type == "tr") ) accession_type = "SwissProt";
		}
		else if ( line.hasPrefix("ref") )
		{
			accession = line.substr(4, line.find('|', 4) - 4);
			accession_type = "NCBI";
		}
		// if it's a swissprot line
		else if ( line.hasPrefix("AC") )
		{
			line.erase(0,2);
			accession = line.trim();
			accession_type = "SwissProt";
		}
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
			std::string::size_type pos1 = line.find('(', 0);
			std::string::size_type pos2;
			if ( pos1 != std::string::npos )
			{
				pos2 = line.find(')', ++pos1);
				if ( pos2 != std::string::npos )
				{
					accession = line.substr(pos1, pos2 - pos1);
					if ( (accession.size() == 6) && (String("OPQ").find(accession[0], 0) != std::string::npos) ) accession_type = "SwissProt";
					else accession.clear();
				}
			}
		}
		if ( accession.empty() )
		{
			accession = line.trim();
			accession_type = "unknown";
		}
	}
	
	bool InspectOutfile::updatePeptideHits(PeptideHit& peptide_hit, std::vector< PeptideHit >& peptide_hits)
	{
		std::vector< PeptideHit >::iterator pep_hit_i;
		for ( pep_hit_i = peptide_hits.begin(); pep_hit_i != peptide_hits.end(); ++pep_hit_i)
		{
			if ( (pep_hit_i->getSequence() == peptide_hit.getSequence()) && (pep_hit_i->getScore() == peptide_hit.getScore()) ) break;
		}
		
		// a peptide hit may only be inserted if it's score type matches the one of the existing hits
		if ( (peptide_hits.empty()) || (peptide_hits[0].getScoreType() == peptide_hit.getScoreType()) )
		{
			// if a new peptide is found, insert it
			if ( pep_hit_i == peptide_hits.end() )
			{
				peptide_hits.push_back(peptide_hit);
				return false;
			}
			// if the peptide has already been inserted, insert additional protein hits
			else
			{
				std::vector< std::pair< String, String > >::iterator prot_hit_i1, prot_hit_i2;
				// remove all protein hits from the peptide that are already in the list
				for ( prot_hit_i1 = peptide_hit.getProteinIndices().begin(); prot_hit_i1 != peptide_hit.getProteinIndices().end(); )
				{
					prot_hit_i2 = std::find(pep_hit_i->getProteinIndices().begin(), pep_hit_i->getProteinIndices().end(), *prot_hit_i1);
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
	
	void InspectOutfile::getPrecursorRTandMZ(const std::vector< std::pair< String, std::vector< unsigned int > > >& files_and_scan_numbers, std::vector< float >&  precursor_retention_times, std::vector< float >& precursor_mz_values)
	{
		std::ifstream ifs;
		String line;
		bool found_scan;
		String search_for = "<scan num=\"";
		unsigned int prefix_length = search_for.length();
		String rt = "retentionTime=\"PT";
		unsigned int rt_length = rt.length();
		String mz = "basePeakMz=\"";
		unsigned int mz_length = mz.length();
		unsigned int pos;
		std::vector< float >::iterator rt_i = precursor_retention_times.begin();
		std::vector< float >::iterator mz_i = precursor_mz_values.begin();
		
		for ( std::vector< std::pair< String, std::vector< unsigned int > > >::const_iterator file_and_scan_numbers = files_and_scan_numbers.begin(); file_and_scan_numbers != files_and_scan_numbers.end(); ++file_and_scan_numbers )
		{
			ifs.open(file_and_scan_numbers->first.c_str());
			if ( ifs )
			{
				for ( std::vector< unsigned int >::const_iterator scan_number = file_and_scan_numbers->second.begin(); scan_number != file_and_scan_numbers->second.end(); ++scan_number, ++rt_i, ++mz_i )
				{
					found_scan = false;
					while ( getline(ifs, line) && !found_scan )
					{
						if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
						line.trim();
						if ( (line.hasPrefix(search_for)) && ((*scan_number) == (unsigned int) atoi(line.substr(prefix_length, line.find('\"', prefix_length) - prefix_length).c_str())) )
						{
							pos = line.find(rt, prefix_length);
							(*rt_i) = atof(line.substr(pos + rt_length , line.find('\"', pos + rt_length) - pos - rt_length).c_str());
							pos = line.find(mz, prefix_length);
							(*mz_i) = atof(line.substr(pos + mz_length , line.find('\"', pos + mz_length) - pos - mz_length).c_str());
							found_scan = true;
						}
					}
				}
				ifs.close();
			}
			else
			{
				rt_i += file_and_scan_numbers->second.size();
				mz_i += file_and_scan_numbers->second.size();
			}
			ifs.clear();
		}
	}

} //namespace OpenMS
