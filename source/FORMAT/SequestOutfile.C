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
/// $Id: SequestOutfile.C,v 1.0 2006/07/12 15:58:59 martinlangwisch Exp $
/// $Author: martinlangwisch $
/// $Maintainer: Martin Langwisch $
/// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SequestOutfile.h>

namespace OpenMS 
{
	SequestOutfile::SequestOutfile()
    : Outfile()
	{}
	
	SequestOutfile::SequestOutfile(const SequestOutfile& source)
    : Outfile(source)
	{}
	
  SequestOutfile::SequestOutfile(const std::string& result_filename)
  	throw (Exception::FileNotFound, Exception::ParseError)
		: Outfile()
  {
  	/// generally used variables
		String line, buffer;
		std::vector< String > substrings;
		
		/// (0) preparations
		/// open the result
		std::ifstream result_file( result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}
		
		/// get the date and time
		DateTime datetime, datetime_empty;
		datetime.clear();
		datetime_empty.clear();
		// ### gucken, ob immer was gelesen wird, sonst exception schmeißen
		/// search for the first line with a date: mm/dd/yyyy
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
				/// 12:00 = 12:00 PM; 24:00 = 12:00 AM
				int hour = atoi(substrings[0].substr(0,2).c_str());
				if ( (hour == 12) && (substrings[1] == "AM") ) substrings[0].replace(0, 2, "00");
				else if ( (hour != 12) && (substrings[1] == "PM") )
				{
					hour += 12;
					substrings[0].replace(0, 2, String(hour));
				}
				substrings[0].append(":00");
				datetime.setTime(substrings[0]);
			}
		}
		if ( datetime == datetime_empty )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "no time found!" , result_filename);
		}
		datetime.get(buffer);
		
		/// get the precursor mass
		if ( !line.empty() ) line.resize(line.length()-1);
		line.trim();
		
		Identification* query;
		if ( line.hasPrefix("(M+H)+ mass = ") )
		{
			line.erase(0, String("(M+H)+ mass = ").length());
			line.split(' ', substrings);
			precursor_mz_values_.push_back(atof(substrings[0].c_str()));
			
			queries_.push_back(Identification());
			query = &queries_.back();
			query->setCharge(atoi(substrings[3].substr(1,2).c_str()));
		}
		else
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "no precursor mass found!" , result_filename);
		}
		
		/*enum columns
		{
			number_column,
			rank_Sp_column,
			id_column,
			MH_column,
			delt_Cn_column,
			XCorr_column,
			Sp_column,
			Sf_column,
			P_column,
			ions_column,
			reference_column,
			peptide_column
		};
		unsigned int number_of_columns = 12;*/
		
		/// skip the next four lines
		for ( unsigned int i = 0; i < 4; ++i)
		{
			if ( !getline(result_file, line) )
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "not enough lines in file!" , result_filename);
			}
		}
		
		/// get the number of peptides displayed
		if ( !line.empty() ) line.resize(line.length()-1);
		line.split(',', substrings);
		if ( !substrings[0].hasPrefix("display top ") )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "no number of displayed peptides found!" , result_filename);
		}
		unsigned int displayed_peptides = String("display top ").length();
		displayed_peptides = atoi(substrings[0].substr(displayed_peptides, substrings[0].find('/', displayed_peptides)));
		
		/// skip the next three lines
		for ( unsigned int i = 0; i < 3; ++i)
		{
			if ( !getline(result_file, line) )
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "not enough lines in file!" , result_filename);
			}
		}
		
		/// retrieve the number of lines and get the needed column numbers //### noch die anderen, nicht benötigten spalten rausnehmen
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
		
		/// get the single columns (seperated by "  ")
		if ( !line.empty() ) line.resize(line.length()-1);
		split(line, "  ", substrings);
		
		for ( std::vector< String >::iterator i = substrings.begin(); i != substrings.end(); )
		{
			i->trim();
			if ( i->empty() ) i = substrings.erase(i);
			else ++i;
		}
		unsigned int number_of_columns = substrings.size();
		
		/// get the numbers of the columns
		for ( std::vector< String >::const_iterator iter = substrings.begin(); iter != substrings.end(); ++iter)
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
		/// check whether the columns are available in the table header
		if ( (number_column == -1) || (rank_Sp_column == -1) ||  (id_column == -1) || (MH_column == -1) || (delt_Cn_column == -1) || (XCorr_column == -1) || (Sp_column == -1) || (ions_column == -1) || (reference_column == -1) || (peptide_column == -1) )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "at least one of the columns '#', 'Rank/Sp', 'Id#', '(M+H)+', 'deltCn', 'XCorr', 'Sp', 'Ions', 'Reference' or 'Peptide' is missing!" , result_filename);
		}
		
		/// map the protein hits according to their accession number in the result file
		std::map< String, unsigned int > acdt_position_map;
		PeptideHit peptide_hit;
		ProteinHit protein_hit;
		std::string accession, accession_type;
		unsigned int id_number;
		unsigned int line_number = 0; /// used to report in which line an error occured
		
		int proteins_per_peptide = 1;
		
		for ( unsigned int i = 0; i < displayed_peptides; ++i )
		{
			getColumns(result_file, substrings, results_filename);
			++line_number;
			
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
			
			/// get the peptide information and insert it
			peptide_hit.clear();
			if ( Sf_column != -1 )
			{
				peptide_hit.setScore(atof(substrings[Sf_column].c_str()));
				peptide_hit.setScoreType("Sf");
			}
			else
			{
				peptide_hit.setScore(atof(substrings[Sp_column].c_str()));
				peptide_hit.setScoreType("Sp");
			}
			peptide_hit.setSequence(substrings[peptide_column]);
			peptide_hit.setRank(atoi(substrings[rank_Sp_column].substr(0, substrings[rank_Sp_column].find('/')).c_str()));
			
			/// get the protein information
			protein_hit.clear();
			get_ac_and_ac_type(substrings[reference_column], result_filename, accession, accession_type);
			protein_hit.setAccession(accession);
			protein_hit.setAccessionType(accession_type);
			protein_hit.setRank(acdt_position_map.size());
			// ### score einfach zusammenrechnen?
			//protein_hit.setScore(0.0);
			//protein_hit.setScoreType("MQScore");
			
			datetime.get(buffer);
			if ( acdt_position_map.insert(std::mape_pair(String(accession+"\t"+buffer), protein_hits_.size())->second ) protein_hits_.push_back(protein_hit);
			
			peptide_hit.addProteinIndex(datetime, accession);
			
			/// check whether there are multiple proteins that belong to this peptide
			proteins_per_peptide = atoi(substrings[reference_column].substr(substrings[reference_colum].find_last_of(' ')));
			for ( unsigned int i = 0; i < proteins_per_peptide; ++i )
			{
				getline(result_file, line);
				line.resize(line.length()-1);
				line.trim();
				split(line, "  ", substrings);
				buffer = substrings[1];
				
				protein_hit.clear();
				get_ac_and_ac_type(buffer, result_filename, accession, accession_type);
				protein_hit.setAccession(accession);
				protein_hit.setAccessionType(accession_type);
				protein_hit.setRank(acdt_position_map.size());
				// ### score einfach zusammenrechnen?
				//protein_hit.setScore(0.0);
				//protein_hit.setScoreType("MQScore");
				
				datetime.get(buffer);
				if ( acdt_position_map.insert(std::mape_pair(String(accession+"\t"+buffer), protein_hits_.size())->second ) protein_hits_.push_back(protein_hit);
				
				peptide_hit.addProteinIndex(datetime, accession);
			}
			
			updatePeptideHits(peptide_hit, peptide_hits_);
		} /// result file read
		
		result_file.close();
		
		/// get the sequences of the protein
		std::vector< String > sequences;
		std::map< String, unsigned int > found, not_found;
		for ( std::vector< String >::const_iterator db_i = databases.begin(); db_i != databases.end(); ++db_i )
		{
			getSequences_(database_path, database_filename, acdt_position_map, sequences, found, not_found);
			
			if ( db_i != databases.begin() )
		}
		
		unsigned int i = 0;
		for ( std::map< String, unsigned int >::iterator protein_i = found.begin(); protein_i != found.end(); ++protein_i)
		{
			protein_hits_[protein_i->second].setSequence(sequences[i]);
			++i;
		}
		sequences.clear();
		found.clear();
		
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
 		ok_ = true;*/
  }
	
	bool SequestOutfile::split_(const String& s, const String& splitter, std::vector<String>& substrings)
	{
		substrings.clear();
		
		unsigned int start = 0;
		unsigned int end = s.find(splitter, start);
		
		while ( end != std::string::npos )
		{
			substrings.push_back(s.substr(start, end-start));
			start = end + splitter.length();
			end = s.find(splitter, start);
		}
		substrings.push_back(s.substr(start));
		
		return substrings.empty();
	}
	
	/// get the columns from a line
	void SequestOutfile::getColumns_(std::ifstream& result_file, std::vector< String >& substrings, String& result_filename)
	{
		if ( !getline(result_file, line) )
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "not enough lines in file!" , result_filename);
		}
		if ( !line.empty() ) line.resize(line.length()-1);
		
		line.trim();
		
		split_(line, "  ", substrings);
std::cout << line << "#blubb: " << substrings.size() << std::endl;
		
		std::vector< String > substrings2;
		char* cp = new char[line.length()];
		float f1, f2;
		bool multiple_columns;
		
		/// remove the empty substrings or merge substrings if one ends with '/'
		for ( std::vector< String >::iterator i = substrings.begin(); i != substrings.end(); )
		{
			i->trim();
			if ( i->empty() ) i = substrings.erase(i);
			else
			{
				/// if two numbers are seperated by one whitespace only, they are assumed to belong to seperate columns
				i->split(' ', substrings2);
				multiple_columns = false;
				std::vector< String >::iterator i2 = i;
				buffer.clear();
std::cout << *i << "\t" << substrings2.size() << std::endl;
				for ( std::vector< String >::const_iterator j = substrings2.begin(); ( (j != (substrings2.end())) && (j != (substrings2.end()-1)) ); ++j)
				{
					if( sscanf(String((*j)+" "+(*(j+1))).c_str(), "%f %f%s", &f1, &f2, cp) == 2 )
					{
						if ( !buffer.empty() ) i2 = substrings.insert(++i2, buffer);
						if ( (*i2) != *j ) i2 = substrings.insert(++i2, *j);
						i2 = substrings.insert(++i2, *(j+1));
						buffer.clear();
						multiple_columns = true;
					}
					else
					{
						buffer.append(*j);
					}
				}
std::cout << "272: " << multiple_columns << std::endl;
				if ( multiple_columns ) substrings.erase(i);
				
				/// merge substrings if one ends with / and both are numbers
				if ( ((*i)[i->length()-1] == '/') && (i != substrings.end()-1) && (sscanf(String(i->substr(0, i->length()-1)+" "+(*(i+1))).c_str(), "%f %f%s", &f1, &f2, cp)  == 2) )
				{
					(i+1)->trim();
					i->append(*(i+1));
					substrings.erase(i+1);
				}
				
				++i;
			}
		}
		delete(cp);
			
for ( std::vector< String >::iterator i = substrings.begin(); i != substrings.end(); ++i)
{
	std::cout << *i << "$";
}
std::cout << std::endl;
	}

	/// retrieve the sequences
	void SequestOutfile::getSequences_(const String& database_path, const String& database_filename, std::map< String, unsigned int > acdt_position_map, std::vector< String >& sequences, std::map< String, unsigned int >& found, std::map< String, unsigned int >& not_found)
	{
		FASTAType data;
		String path_and_file = database_path;
		ensurePathChar(path_and_file);
		path_and_file.append(database_filename);
		FASTAFile().load(path_and_file, data);
		sequences.clear();
		std::vector< std::pair< String, String > >::iterator fa_i;
		
		for ( std::set< String >::const_iterator i = acdt_position_map.begin(); i != acdt_position_map.end(); ++i )
		{
			fa_i = find(data.begin(), data.end(), i->substr(0, i->find('\t')));
			if ( fa_i != data.end() )
			{
				sequences.push_back(fa_i->second);
				found.insert(*i);
			}
			else not_found.insert(*i);
		}
	}
	
} //namespace OpenMS
