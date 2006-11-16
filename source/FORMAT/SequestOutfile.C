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

#include <OpenMS/FORMAT/SequestOutfile.h>

using namespace std;

namespace OpenMS 
{
	
	SequestOutfile::SequestOutfile()
	{
	}

  void SequestOutfile::load(const string& result_filename, vector< Identification >&	identifications, ProteinIdentification&	protein_identification, vector< float >& precursor_retention_times, vector< float >& precursor_mz_values, const double& p_value_threshold, const string& database, const string& snd_database) throw (Exception::FileNotFound, Exception::ParseError)
  {
  	// generally used variables
		String line, buffer;
		vector< String > substrings;
		
		// map the protein hits according to their accession number in the result file
		map< String, unsigned int > ac_position_map;
		
		// get the protein hits that have already been found in another out-file
		vector< ProteinHit > protein_hits = protein_identification.getProteinHits();
		// and insert them into the map
		for ( vector< ProteinHit >::const_iterator phit_i = protein_hits.begin(); phit_i != protein_hits.end(); ++phit_i )
		{
			ac_position_map.insert(make_pair(phit_i->getAccession(), ac_position_map.size()));
		}
		
		// (0) preparations
		// open the result
		ifstream result_file( result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}
		
		// get the date and time
		DateTime datetime, datetime_empty;
		datetime.clear();
		datetime_empty.clear();
		// search for the first line with a date: mm/dd/yyyy
		while ( (getline(result_file, line)) && (datetime == datetime_empty) )
		{
			if ( !line.empty() ) line.resize(line.length()-1);
			line.trim();
			line.split(',', substrings);
			
			if ( (!substrings.empty()) && (substrings[0].size() == 10) && (substrings[0][2] == '/') && (substrings[0][5] == '/') )
			{
				datetime.setDate(substrings[0]);
				buffer = substrings[1];
				buffer.trim();
				buffer.split(' ', substrings);
				// 12:00 = 12:00 PM; 24:00 = 12:00 AM
				int hour = atoi(substrings[0].substr(0,2).c_str());
				if ( (hour == 12) && (substrings[1] == "AM") ) substrings[0].string::replace(0, 2, "00");
				else if ( (hour != 12) && (substrings[1] == "PM") )
				{
					hour += 12;
					substrings[0].string::replace(0, 2, String(hour));
					string x = "abc";
					x.string::replace(0,3, "def");
				}
				substrings[0].append(":00");
				datetime.setTime(substrings[0]);
			}
		}
		if ( datetime == datetime_empty )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "no time found!" , result_filename);
		}
		
		// get the precursor mass
		if ( !line.empty() ) line.resize(line.length()-1);
		line.trim();
		
		Identification* query;
		if ( line.hasPrefix("(M+H)+ mass = ") )
		{
			line.erase(0, String("(M+H)+ mass = ").length());
			line.split(' ', substrings);
			precursor_mz_values.push_back(atof(substrings[0].c_str()));
			
			identifications.push_back(Identification());
			query = &identifications.back();
			query->setCharge(atoi(substrings[3].substr(1,2).c_str()));
		}
		else
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "no precursor mass found!" , result_filename);
		}
		
		vector< String > databases;
		
		databases.push_back(database);
		if ( !snd_database.empty() ) databases.push_back(snd_database);

		while ( getline(result_file, line) ) // skip all lines until the one with 'display top'
		{
			if ( !line.empty() ) line.resize(line.length()-1);
			line.trim();
			if ( line.hasPrefix("display top") ) break;
		}
		// get the number of peptides displayed
		line.split(',', substrings);
		if ( !substrings[0].hasPrefix("display top ") )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "No number of displayed peptides found!" , result_filename);
		}
		unsigned int displayed_peptides = String("display top ").length();
		displayed_peptides = atoi(substrings[0].substr(displayed_peptides, substrings[0].find('/', displayed_peptides)).c_str());
		
		// skip the next three lines
		while ( getline(result_file, line) )
		{
			if ( !line.empty() ) line.resize(line.length()-1);
			line.trim();
			if ( line.hasPrefix("#") ) break;
		}
		if ( !line.hasPrefix("#") )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Not enough lines in file (no table header found)!" , result_filename);
		}
		
		// retrieve the number of lines and get the needed column numbers
		int number_column = -1;
		int rank_Sp_column = -1;
		int id_column = -1;
		int MH_column = -1;
		int delt_Cn_column = -1;
		int XCorr_column = -1;
		int Sp_column = -1;
		int Sf_column = -1;
		//int P_column = -1;
		int ions_column = -1;
		int reference_column = -1;
		int peptide_column = -1;
		
		// get the single columns (seperated by "  ")
		split(line, "  ", substrings);
		
		// remove the empty strings
		for ( vector< String >::iterator i = substrings.begin(); i != substrings.end(); )
		{
			i->trim();
			if ( i->empty() ) i = substrings.erase(i);
			else ++i;
		}
		unsigned int number_of_columns = substrings.size();
		
		// get the numbers of the columns
		for ( vector< String >::const_iterator iter = substrings.begin(); iter != substrings.end(); ++iter)
		{
			if ( !iter->compare("#") )
			{
				number_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Rank/Sp") )
			{
				rank_Sp_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Id#") )
			{
				id_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("(M+H)+") )
			{
				MH_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("deltCn") )
			{
				delt_Cn_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("XCorr") )
			{
				XCorr_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Sp") )
			{
				Sp_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Sf") )
			{
				Sf_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Ions") )
			{
				ions_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Reference") )
			{
				reference_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Peptide") )
			{
				peptide_column = (iter - substrings.begin());
			}
		}
		// check whether the columns are available in the table header
		if ( (number_column == -1) || (rank_Sp_column == -1) ||  (id_column == -1) || (MH_column == -1) || (delt_Cn_column == -1) || (XCorr_column == -1) || (Sp_column == -1) || (ions_column == -1) || (reference_column == -1) || (peptide_column == -1) )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "at least one of the columns '#', 'Rank/Sp', 'Id#', '(M+H)+', 'deltCn', 'XCorr', 'Sp', 'Ions', 'Reference' or 'Peptide' is missing!" , result_filename);
		}
		// skip the next line
		if ( !getline(result_file, line) )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "not enough lines in file!" , result_filename);
		}
		if ( !line.empty() ) line.resize(line.length()-1);
		
		PeptideHit peptide_hit;
		ProteinHit protein_hit;
		string accession, accession_type;
		unsigned int line_number = 0; // used to report in which line an error occured
		
		int proteins_per_peptide;
		int score_column = (Sf_column == -1) ? Sp_column : Sf_column;
		String score_type = (Sf_column == -1) ? "SEQUEST prelim." : "SEQUEST";
		
		for ( unsigned int i = 0; i < displayed_peptides; ++i )
		{
			if ( !getline(result_file, line) )
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "not enough lines in file!" , result_filename);
			}
			if ( !line.empty() ) line.resize(line.length()-1);
			
			if ( !getColumns(line, substrings, number_of_columns, reference_column) ) break; // if less peptides were found than may be displayed, break
			++line_number;
			
			// check whether the line has enough columns
			if (substrings.size() < number_of_columns )
			{
				stringstream error_message;
				error_message << "wrong number of columns in row " << line_number << " " << line << " " << "! (" 
											<< substrings.size() << " present, should be " << number_of_columns << ")";
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, error_message.str().c_str() , result_filename);
			}
			
			// get the peptide information and insert it
			peptide_hit.clear();
			peptide_hit.setScore(atof(substrings[score_column].c_str()));
			peptide_hit.setScoreType(score_type);
			
			peptide_hit.setSequence(substrings[peptide_column].substr(2, substrings[peptide_column].length()-4));
			peptide_hit.setRank(atoi(substrings[rank_Sp_column].substr(0, substrings[rank_Sp_column].find('/')).c_str()));
			
			// get the protein information
			
			// check whether there are multiple proteins that belong to this peptide
			if ( substrings[reference_column].find_last_of('+') != string::npos )
			{
				// save the number of multiple proteins
				proteins_per_peptide = atoi(substrings[reference_column].substr(substrings[reference_column].find_last_of('+')).c_str());
				// and remove this number from the string
				substrings[reference_column].resize(substrings[reference_column].find_last_of('+'));
			}
			else proteins_per_peptide = 0;
			
			protein_hit.clear();
			getACAndACType(substrings[reference_column], accession, accession_type);
			protein_hit.setAccession(accession);
			protein_hit.setAccessionType(accession_type);
			protein_hit.setRank(ac_position_map.size());
			// ### score einfach zusammenrechnen?
			//protein_hit.setScore(0.0);
			//protein_hit.setScoreType(score_type);
			
			if ( ac_position_map.insert(make_pair(accession, protein_hits.size())).second ) protein_hits.push_back(protein_hit);
			
			peptide_hit.addProteinIndex(datetime, accession);
			
			for ( int i = 0; i < proteins_per_peptide; ++i )
			{
				getline(result_file, line);
				if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
				line.trim();
				split(line, "  ", substrings);
				
				protein_hit.clear();
				getACAndACType(substrings[1], accession, accession_type);
				protein_hit.setAccession(accession);
				protein_hit.setAccessionType(accession_type);
				protein_hit.setRank(ac_position_map.size());
				// ### score einfach zusammenrechnen?
				//protein_hit.setScore(0.0);
				//protein_hit.setScoreType(score_type);
				
				if ( ac_position_map.insert(make_pair(accession, protein_hits.size())).second ) protein_hits.push_back(protein_hit);
				
				peptide_hit.addProteinIndex(datetime, accession);
			}
			
			updatePeptideHits(peptide_hit, query->getPeptideHits());
		} // result file read
		result_file.close();
		
		// get the sequences of the protein
		vector< String > sequences;
		map< String, unsigned int > not_found;
		vector< pair< String, unsigned int > > found;
		for ( vector< String >::const_iterator db_i = databases.begin(); db_i != databases.end(); ++db_i )
		{
			getSequences(*db_i, ac_position_map, sequences, found, not_found);
			
			if ( db_i != databases.begin() ) ac_position_map = not_found;
		}
		
		vector< String >::const_iterator seq_i = sequences.begin();
		for ( vector< pair< String, unsigned int > >::const_iterator protein_i = found.begin(); protein_i != found.end(); ++protein_i, ++seq_i) protein_hits[protein_i->second].setSequence(*seq_i);
		
		protein_identification.setProteinHits(protein_hits);
		protein_identification.setDateTime(datetime);
		protein_hits.clear();
		
		sequences.clear();
		found.clear();
		not_found.clear();
		
		// if there's but one query the protein hits are inserted there instead of a ProteinIdentification object
		if ( identifications.size() == 1 )
		{
			query->setProteinHits(protein_hits);
			query->setDateTime(datetime);
			query->setPeptideSignificanceThreshold(p_value_threshold);
		}
		else // if several out-files were used, the first query owns a vector with protein hits
		{
			identifications.front().setProteinHits(vector< ProteinHit >());
		}
		
		peptide_hit.clear();
		protein_hit.clear();
		
		ac_position_map.clear();
  }
	
	bool SequestOutfile::split(const String& s, const String& splitter, vector<String>& substrings)
	{
		substrings.clear();
		
		string::size_type start = 0;
		string::size_type end = s.find(splitter, start);
		
		while ( end != string::npos )
		{
			substrings.push_back(s.substr(start, end-start));
			start = end + splitter.length();
			end = s.find(splitter, start);
		}
		substrings.push_back(s.substr(start));
		
		return substrings.empty();
	}
	
	// get the columns from a line
	bool SequestOutfile::getColumns(const String& line, vector< String >& substrings, unsigned int number_of_columns, unsigned int reference_column)
	{
		String buffer;

		if ( line.empty() ) return false;
		
		line.split(' ', substrings);

		// remove any empty strings
		for ( vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); )
		{
			s_i->trim();
			if ( s_i->empty() ) substrings.erase(s_i);
			else ++s_i;
		}

		for ( vector< String >::iterator s_i = substrings.begin(); s_i != substrings.end(); )
		{
		// if there are three columns, the middle one being a '/', they are merged
			if ( s_i+1 != substrings.end() )
			{
				if ( ((*(s_i+1)) == "/") && (s_i+2 != substrings.end()) )
				{
					s_i->append(*(s_i+1));
					s_i->append(*(s_i+2));
					substrings.erase(s_i+2);
					substrings.erase(s_i+1);
				}
				// if there are two columns, and the first ends with, or the second starts with a '/', they are merged
				else if ( (*(s_i+1))[0] == '/' )
				{
					s_i->append(*(s_i+1));
					substrings.erase(s_i+1);
				}
				else if ( (*(s_i))[s_i->length()-1] == '/' )
				{
					s_i->append(*(s_i+1));
					substrings.erase(s_i+1);
				}
			// if there are two columns and the second is a number preceeded by a '+', they are merged
				else if ( (*(s_i+1))[0] == '+' )
				{
					bool is_digit = true;
					for ( unsigned int i = 1; i < (s_i+1)->length(); ++i ) is_digit &= isdigit((*(s_i+1))[i]);
					if ( is_digit && ((s_i+1)->length()-1) )
					{
						s_i->append(*(s_i+1));
						substrings.erase(s_i+1);
					}
				}
				else ++s_i;
			}
			else ++s_i;
		}

		// if there are more columns than should be, there were spaces in the protein column
		for ( vector< String >::iterator s_i = substrings.begin()+reference_column; substrings.size() > number_of_columns; )
		{
			s_i->append(*(s_i+1));
			substrings.erase(s_i+1);
		}
		
		return true;
	}

	// retrieve the sequences
	void SequestOutfile::getSequences(const String& database_filename, map< String, unsigned int > ac_position_map, vector< String >& sequences, vector< pair< String, unsigned int > >& found, map< String, unsigned int >& not_found) throw (Exception::FileNotFound)
	{
		sequences.clear();
		ifstream ifs(database_filename.c_str());
		if ( !ifs )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
		}
		
		String line, accession, accession_type, sequence;
		
		not_found = ac_position_map;
		map< String, unsigned int >::iterator i = not_found.begin();
		while ( getline(ifs, line) && !not_found.empty() )
		{
			if ( line.hasPrefix(">") )
			{
				if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
				getACAndACType(line, accession, accession_type);
				break;
			}
		}
		while ( getline(ifs, line) && !not_found.empty() )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			i = not_found.find(accession);
			if ( i != not_found.end() )
			{
				// get the sequence
				sequence = line;
				while ( getline(ifs, line) )
				{
					if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
					if ( !line.hasPrefix(">") )	sequence.append(line);
					else break;
				}
				sequences.push_back(sequence);
				found.push_back(*i);
				not_found.erase(i);
				getACAndACType(line, accession, accession_type);
			}
			else // search for the next entry
			{
				while ( getline(ifs, line) )
				{
					if ( line.hasPrefix(">") )
					{
						if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
						getACAndACType(line, accession, accession_type);
						break;
					}
				}
			}
		}
		
		ifs.close();
		ifs.clear();
	}
	
	void SequestOutfile::getACAndACType(String line, string& accession, string& accession_type) throw (Exception::ParseError)
	{
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
			string::size_type snd = line.find('|', 3);
			string::size_type third = 0;
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
			string::size_type pos1 = line.find('(', 0);
			string::size_type pos2;
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
		}
		if ( accession.empty() )
		{
			accession = line.trim();
			accession_type = "unknown";
		}
	}
	
	bool SequestOutfile::updatePeptideHits(PeptideHit& peptide_hit, vector< PeptideHit >& peptide_hits)
	{
		vector< PeptideHit >::iterator pep_hit_i;
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
	
} //namespace OpenMS
