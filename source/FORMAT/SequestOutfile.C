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

#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;

namespace OpenMS 
{
	SequestOutfile::SequestOutfile():
		out2summary_number(0)
	{}
	
  void
	SequestOutfile::load(
		const string& result_filename,
		vector< IdentificationData >&	identifications,
		ProteinIdentification&	protein_identification,
		const Real& p_value_threshold,
		vector< Real >& pvalues,
		const string& database)
	throw(
		Exception::FileNotFound,
		Exception::ParseError)
  {
		// if no p_values were computed take all peptides
		bool no_pvalues = pvalues.empty();
		if ( no_pvalues ) pvalues.push_back(0.0); // to make sure pvalues.end() is never reached
		
  	// generally used variables
		String line, buffer;
		vector< String > substrings;
		
		// map the protein hits according to their accession number in the result file
		map< String, UnsignedInt > ac_position_map;
		
		// get the protein hits that have already been found in another out-file
		vector< ProteinHit > protein_hits = protein_identification.getProteinHits();
		// and insert them SignedInto the map
		for ( vector< ProteinHit >::const_iterator phit_i = protein_hits.begin(); phit_i != protein_hits.end(); ++phit_i )
		{
			ac_position_map.insert(make_pair(phit_i->getAccession(), ac_position_map.size()));
		}
		
		String
			accession,
			accession_type,
			score_type;
		
		DateTime datetime;
		Real precursor_mz_value;
		UnsignedInt
			precursor_mass_type,
			ion_mass_type,
			number_of_columns,
			displayed_peptides,
			proteins_per_peptide;
			
		SignedInt
			charge,
			number_column,
			rank_sp_column,
			id_column,
			mh_column,
			delta_cn_column,
			xcorr_column,
			sp_column,
			sf_column,
			ions_column,
			reference_column,
			peptide_column,
			score_column;
		
		readOutHeader(result_filename, datetime, precursor_mz_value, charge, precursor_mass_type, ion_mass_type, number_column, rank_sp_column, id_column, mh_column, delta_cn_column, xcorr_column, sp_column, sf_column, ions_column, reference_column, peptide_column, score_column, number_of_columns, displayed_peptides);
		
		// open the result
		ifstream result_file(result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}
		
		UnsignedInt line_number = 1; // used to report in which line an error occured
		while ( getline(result_file, line) ) // skip all lines until the one with '---'
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.trim();
			++line_number;
			if ( line.hasPrefix("---") ) break;
		}
		
		Identification* query;
		identifications.push_back(IdentificationData());
		identifications.back().mz = precursor_mz_value;
		
		query = &(identifications.back().id);
		
// 		vector< String > databases;
// 		databases.push_back(database);
		
		PeptideHit peptide_hit;
		ProteinHit protein_hit;
		
		score_type = (sf_column == -1) ? "SEQUEST prelim." : "SEQUEST";

		if ( no_pvalues ) pvalues.insert(pvalues.end(), displayed_peptides, 0.0);
		vector< Real >::const_iterator p_value = pvalues.begin();
		
		for ( UnsignedInt i = 0; i < displayed_peptides; ++i, ++p_value )
		{
			++line_number;
			// if less peptides were found than may be displayed, break
			if ( !getline(result_file, line) ) break;
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.trim();
			if ( line.empty() ) continue; // skip empty lines
			
			getColumns(line, substrings, number_of_columns, reference_column);
			
			// check whether the line has enough columns
			if (substrings.size() < number_of_columns )
			{
				stringstream error_message;
				error_message << "Wrong number of columns in line " << line_number << "! (" << substrings.size() << " present, should be " << number_of_columns << ")";
				result_file.close();
				result_file.clear();
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, error_message.str().c_str() , result_filename);
			}
			
			// check whether there are multiple proteins that belong to this peptide
			if ( substrings[reference_column].find_last_of('+') != string::npos )
			{
				// save the number of multiple proteins
				proteins_per_peptide = atoi(substrings[reference_column].substr(substrings[reference_column].find_last_of('+')).c_str());
				// and remove this number from the string
				substrings[reference_column].resize(substrings[reference_column].find_last_of('+'));
			}
			else proteins_per_peptide = 0;
			// get the peptide information and insert it
			if ( p_value != pvalues.end() && *p_value <= p_value_threshold )
			{
				peptide_hit.clear();
				peptide_hit.setScore(atof(substrings[score_column].c_str()));
				peptide_hit.setScoreType(score_type);
				peptide_hit.setCharge(charge);
				
				peptide_hit.setSequence(substrings[peptide_column].substr(2, substrings[peptide_column].length()-4));
				peptide_hit.setRank(atoi(substrings[rank_sp_column].substr(0, substrings[rank_sp_column].find('/')).c_str()));
				
				// get the protein information
				protein_hit.clear();
				getACAndACType(substrings[reference_column], accession, accession_type);
				protein_hit.setAccession(accession);
				protein_hit.setAccessionType(accession_type);
				protein_hit.setRank(ac_position_map.size());
				// score einfach zusammenrechnen?
	//			protein_hit.setScore(0.0);
	//			protein_hit.setScoreType(score_type);
				
				if ( ac_position_map.insert(make_pair(accession, protein_hits.size())).second ) protein_hits.push_back(protein_hit);
				
				peptide_hit.addProteinIndex(datetime, accession);
				
				for ( UnsignedInt i = 0; i < proteins_per_peptide; ++i )
				{
					getline(result_file, line);
					if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
					line.trim();
					// all these lines look like '0  accession', e.g. '0  gi|1584947|prf||2123446B gamma sar'
					if ( !line.hasPrefix("0  ") ) // if the line doesn't look like that
					{
						stringstream error_message;
						error_message << "Line " << line_number << " doesn't look like a line with additional found proteins! (Should look like this: 0  gi|1584947|prf||2123446B gamma sar)";
						result_file.close();
						result_file.clear();
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, error_message.str().c_str() , result_filename);
					}
					line.erase(0,3);
					
					protein_hit.clear();
					getACAndACType(line, accession, accession_type);
					protein_hit.setAccession(accession);
					protein_hit.setAccessionType(accession_type);
					protein_hit.setRank(ac_position_map.size());
					// score einfach zusammenrechnen?
	//				protein_hit.setScore(0.0);
	//				protein_hit.setScoreType(score_type);
					
					if ( ac_position_map.insert(make_pair(accession, protein_hits.size())).second ) protein_hits.push_back(protein_hit);
					
					peptide_hit.addProteinIndex(datetime, accession);
				}
				
				updatePeptideHits(peptide_hit, query->getPeptideHits());
			}
			else
			{
				for ( UnsignedInt i = 0; i < proteins_per_peptide; ++i ) getline(result_file, line);
			}
		}
		result_file.close();
		if ( no_pvalues ) pvalues.clear();
		
		// get the sequences of the protein
//		vector< String > sequences;
//		map< String, UnsignedInt > not_found;
//		vector< pair< String, UnsignedInt > > found;
//		for ( vector< String >::const_iterator db_i = databases.begin(); db_i != databases.end(); ++db_i )
//		{
//			getSequences(*db_i, ac_position_map, sequences, found, not_found);
//			
//			if ( db_i != databases.begin() ) ac_position_map = not_found;
//		}
//
//		vector< String >::const_iterator seq_i = sequences.begin();
//		for ( vector< pair< String, UnsignedInt > >::const_iterator protein_i = found.begin(); protein_i != found.end(); ++protein_i, ++seq_i)
//		{
//			protein_hits[protein_i->second].setSequence(*seq_i);
//		}
//		
//		sequences.clear();
//		found.clear();
//		not_found.clear();
		
		if ( protein_hits.empty() ) identifications.pop_back();
		else
		{
			protein_identification.setProteinHits(protein_hits);
			protein_identification.setDateTime(datetime);
			protein_hits.clear();
		}
		
		// if there's but one query the protein hits are inserted there instead of a ProteinIdentification object
		if ( !identifications.empty() )
		{
			if ( identifications.size() == 1 )
			{
				query->setProteinHits(protein_hits);
				query->setDateTime(datetime);
				query->setPeptideSignificanceThreshold(p_value_threshold);
			}
			else identifications.front().id.setProteinHits(vector< ProteinHit >()); // if several out-files were used, the first query owns a vector with protein hits
		}
		
		peptide_hit.clear();
		protein_hit.clear();
		
		ac_position_map.clear();
  }
	
	// get the columns from a line
	bool
	SequestOutfile::getColumns(
		const String& line,
		vector< String >& substrings,
		UnsignedInt number_of_columns,
		UnsignedInt reference_column)
	{
		String buffer;

		if ( line.empty() ) return false;
		
		line.split(' ', substrings);

		// remove any empty strings
		substrings.erase(remove(substrings.begin(),substrings.end(),""),substrings.end());

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
					for ( UnsignedInt i = 1; i < (s_i+1)->length(); ++i ) is_digit &= (bool) isdigit((*(s_i+1))[i]);
					if ( is_digit && ((s_i+1)->length()-1) )
					{
						s_i->append(*(s_i+1));
						substrings.erase(s_i+1);
					}
					else ++s_i;
				}
				else ++s_i;
			}
			else ++s_i;
		}

		// if there are more columns than should be, there were spaces in the protein column
		for ( vector< String >::iterator s_i = substrings.begin()+reference_column; substrings.size() > number_of_columns; )
		{
			s_i->append(" ");
			s_i->append(*(s_i+1));
			substrings.erase(s_i+1);
		}
		
		return true;
	}

	// retrieve the sequences
	void
	SequestOutfile::getSequences(
		const String& database_filename,
		const map< String, UnsignedInt >& ac_position_map,
		vector< String >& sequences,
		vector< pair< String, UnsignedInt > >& found,
		map< String, UnsignedInt >& not_found)
	throw (
		Exception::FileNotFound)
	{
		ifstream database_file(database_filename.c_str());
		if ( !database_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, database_filename);
		}
		
		String line, accession, accession_type, sequence;
		not_found = ac_position_map;
		map< String, UnsignedInt >::iterator nf_i = not_found.end();
		while ( getline(database_file, line) && !not_found.empty() )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.trim();
			
			// empty and comment lines are skipped
			if ( line.empty() || line.hasPrefix(";") ) continue;

			// the sequence belonging to the predecessing protein ('>') is stored, so when a new protein ('>') is found, save the sequence of the old protein
			if ( line.hasPrefix(">") )
			{
				getACAndACType(line, accession, accession_type);
				if ( nf_i != not_found.end() )
				{
					sequences.push_back(sequence);
					found.push_back(*nf_i);
					not_found.erase(nf_i);
				}
				nf_i = not_found.find(accession); // for the first protein in the database, there's no predecessing protein
				sequence.clear();
			}
			else if ( nf_i != not_found.end() ) sequence.append(line);
		}
		if ( nf_i != not_found.end() )
		{
			sequences.push_back(sequence);
			found.push_back(*nf_i);
			not_found.erase(nf_i);
		}
		
		database_file.close();
		database_file.clear();
	}
	
	void
	SequestOutfile::getACAndACType(
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
	SequestOutfile::updatePeptideHits(
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
	SequestOutfile::readOutHeader(
		const string& result_filename,
		DateTime& datetime,
		Real& precursor_mz_value,
		SignedInt& charge,
		UnsignedInt& precursor_mass_type,
		UnsignedInt& ion_mass_type,
		SignedInt& number_column,
		SignedInt& rank_sp_column,
		SignedInt& id_column,
		SignedInt& mh_column,
		SignedInt& delta_cn_column,
		SignedInt& xcorr_column,
		SignedInt& sp_column,
		SignedInt& sf_column,
//		SignedInt& P_column,
		SignedInt& ions_column,
		SignedInt& reference_column,
		SignedInt& peptide_column,
		SignedInt& score_column,
		UnsignedInt& number_of_columns,
		UnsignedInt& displayed_peptides)
	throw(
		Exception::FileNotFound,
		Exception::ParseError)
	{
		// open the result
		ifstream result_file( result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}

		String line, buffer;
		vector< String > substrings;
		
		// get the date and time
		DateTime datetime_empty;
		datetime.clear();
		datetime_empty.clear();
		// search for the first line with a date: mm/dd/yyyy
		while ( (getline(result_file, line)) && (datetime == datetime_empty) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.trim();
			line.split(',', substrings);
			
			if ( (!substrings.empty()) && (substrings[0].size() == 10) && (substrings[0][2] == '/') && (substrings[0][5] == '/') )
			{
				datetime.setDate(substrings[0]);
				buffer = substrings[1];
				buffer.trim();
				buffer.split(' ', substrings);
				// 12:00 = 12:00 PM; 24:00 = 12:00 AM
				SignedInt hour = atoi(substrings[0].substr(0,2).c_str());
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
			result_file.close();
			result_file.clear();
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "no time found!" , result_filename);
		}
		
		// get the precursor mass
		if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
		line.trim();

		if ( line.hasPrefix("(M+H)+ mass = ") )
		{
			line.erase(0, String("(M+H)+ mass = ").length());
			line.split(' ', substrings);
			precursor_mz_value = substrings[0].toFloat();
			charge = substrings[3].substr(1,2).toInt();
			line = *(--substrings.end());
			line.split('/', substrings);
			substrings[0].toUpper();
			substrings[1].toUpper();
			precursor_mass_type = ( substrings[0] == "MONO" ) ? 1 : 0;
			ion_mass_type = ( substrings[1] == "MONO" ) ? 1 : 0;
		}
		else
		{
			result_file.close();
			result_file.clear();
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "no precursor mass found!" , result_filename);
		}
		
		while ( getline(result_file, line) ) // skip all lines until the one with 'display top'
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.trim();
			if ( line.hasPrefix("display top") ) break;
		}
		// get the number of peptides displayed
		line.split(',', substrings);
		if ( !substrings[0].hasPrefix("display top ") )
		{
			result_file.close();
			result_file.clear();
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "No number of displayed peptides found!" , result_filename);
		}
		displayed_peptides = String("display top ").length();
		displayed_peptides = atoi(substrings[0].substr(displayed_peptides, substrings[0].find('/', displayed_peptides)).c_str());
		
		// skip the next three lines
		while ( getline(result_file, line) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.trim();
			if ( line.hasPrefix("#") ) break;
		}
		if ( !line.hasPrefix("#") )
		{
			result_file.close();
			result_file.clear();
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Not enough lines in file (no table header found)!" , result_filename);
		}
		
		// retrieve the number of lines and get the needed column numbers
		number_column = -1;
		rank_sp_column = -1;
		id_column = -1;
		mh_column = -1;
		delta_cn_column = -1;
		xcorr_column = -1;
		sp_column = -1;
		sf_column = -1;
		//P_column = -1;
		ions_column = -1;
		reference_column = -1;
		peptide_column = -1;
		
		// get the single columns
		line.split(' ', substrings);
		
		// remove the empty strings
		for ( vector< String >::iterator i = substrings.begin(); i != substrings.end(); )
		{
			i->trim();
			if ( i->empty() ) i = substrings.erase(i);
			else ++i;
		}
		number_of_columns = substrings.size();
		
		// get the numbers of the columns
		for ( vector< String >::const_iterator iter = substrings.begin(); iter != substrings.end(); ++iter)
		{
			if ( !iter->compare("#") )
			{
				number_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Rank/Sp") )
			{
				rank_sp_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Id#") )
			{
				id_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("(M+H)+") )
			{
				mh_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("deltCn") )
			{
				delta_cn_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("XCorr") )
			{
				xcorr_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Sp") )
			{
				sp_column = (iter - substrings.begin());
			}
			else if ( !iter->compare("Sf") )
			{
				sf_column = (iter - substrings.begin());
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
		if ( (number_column == -1) || (rank_sp_column == -1) ||  (id_column == -1) || (mh_column == -1) || (delta_cn_column == -1) || (xcorr_column == -1) || (sp_column == -1) || (ions_column == -1) || (reference_column == -1) || (peptide_column == -1) )
		{
			result_file.close();
			result_file.clear();
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "at least one of the columns '#', 'Rank/Sp', 'Id#', '(M+H)+', 'deltCn', 'XCorr', 'Sp', 'Ions', 'Reference' or 'Peptide' is missing!" , result_filename);
		}
		
		score_column = (sf_column == -1) ? sp_column : sf_column;
		
		result_file.close();
		result_file.clear();
	}

	void
	SequestOutfile::out2SummaryHtml(
		const string& out_filename,
		const string& summary_filename,
		const string& database_filename)
	throw(
		Exception::FileNotFound,
		Exception::ParseError,
		Exception::UnableToCreateFile)
	{
		ofstream summary;
		// write the fileheader, if not in append mode (one fileheader only)
		if ( !out2summary_number )
		{
			++out2summary_number;
			summary.open(summary_filename.c_str());
			if ( !summary )
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, summary_filename);
			}
			summary << "<HTML>" << endl << "<HEAD><TITLE>HTML-SUMMARY</TITLE></HEAD>" << endl << "<BODY BGCOLOR=\"#FFFFFF\">" << endl << "<PRE><FONT COLOR=\"green\">HTML-SUMMARY</FONT>" << endl  << endl << "<FONT COLOR=\"green\">   #    File         MH+                  XCorr    dCn      Sp    RSp    Ions   Ref             Sequence</FONT>" << endl << "<FONT COLOR=\"green\">  ----  ---------  -------------------    ------  -----   ------  ---   ------  -----------  ----------</FONT>" << endl;
		}
		else summary.open(summary_filename.c_str(), ios::out | ios::app);
		if ( !summary )
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, summary_filename);
		}
		
		DateTime datetime;
		Real
			precursor_mz_value(0),
			last_delta_cn(0);
			UnsignedInt
			precursor_mass_type(0),
			ion_mass_type(0),
			number_of_columns(0),
			displayed_peptides(0);
		
		SignedInt
			charge(0),
			number_column(0),
			rank_sp_column(0),
			id_column(0),
			mh_column(0),
			delta_cn_column(0),
			xcorr_column(0),
			sp_column(0),
			sf_column(0),
			ions_column(0),
			reference_column(0),
			peptide_column(0),
			score_column(0);
			
		readOutHeader(out_filename, datetime, precursor_mz_value, charge, precursor_mass_type, ion_mass_type, number_column, rank_sp_column, id_column, mh_column, delta_cn_column, xcorr_column, sp_column, sf_column, ions_column, reference_column, peptide_column, score_column, number_of_columns, displayed_peptides);
		
		// open the result file
		ifstream out_file(out_filename.c_str());
		if ( !out_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, out_filename);
		}

		String out_filename_base = File::basename(out_filename);
		out_filename_base.erase(out_filename_base.length() - 4);

		String line, line_buffer, peptide;
		stringstream line_ss;
		vector< String > substrings;
		UnsignedInt line_number = 0;
		while ( getline(out_file, line) ) // skip all lines until the one with '---'
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.trim();
			++line_number;
			if ( line.hasPrefix("---") ) break;
		}
		
		UnsignedInt proteins_per_peptide;
		UnsignedInt i = 0;
		for ( ; i < displayed_peptides; )
		{
			if ( !getline(out_file, line) )
			{
// 				out_file.close();
// 				out_file.clear();
// 				summary.close();
// 				summary.clear();
// 				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "not enough lines in file!" , out_filename);
				break;
			}
			++line_number;
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.trim();
			if ( line.empty() ) continue;
			
			getColumns(line, substrings, number_of_columns, reference_column);
			
			// check whether the line has enough columns
			if (substrings.size() < number_of_columns )
			{
				stringstream error_message;
				error_message << "Wrong number of columns in line " << line_number << "! (" << substrings.size() << " present, should be " << number_of_columns << ")";
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, error_message.str().c_str() , out_filename);
			}
			
			// check whether there are multiple proteins that belong to this peptide
			if ( substrings[reference_column].find_last_of('+') != string::npos )
			{
				// save the number of multiple proteins
				proteins_per_peptide = atoi(substrings[reference_column].substr(substrings[reference_column].find_last_of('+')).c_str());
				// and remove this number from the string
				substrings[reference_column].resize(substrings[reference_column].find_last_of('+'));
			}
			else proteins_per_peptide = 0;
			if ( i ) line_buffer.replace(line_buffer.find("###"), 3, String(substrings[delta_cn_column].toFloat() - last_delta_cn)); // the delta_cn for the the preceeding value is computed using this value
			++i;
			summary << line_buffer;
			last_delta_cn = substrings[delta_cn_column].toFloat();
			
			peptide = substrings[peptide_column].substr(2, substrings[peptide_column].length() - 4);
			
			substrings[reference_column].remove('>');
			line_ss.str("");
			line_ss << "     " << (out2summary_number + i) << "  <A HREF=\"\">" << out_filename_base << "</A>  " << substrings[mh_column] << " (" << precursor_mz_value - substrings[mh_column].toFloat() << ")  " << substrings[xcorr_column] << "  " << "###" << " " << substrings[sp_column] << " " << substrings[rank_sp_column].substr(substrings[rank_sp_column].find('/') + 1) << " <A HREF=\"\">" << substrings[ions_column] << "</A>   <A HREF=\"&amp;Db=" << database_filename << "&amp;\">" <<  substrings[reference_column] << "</A>   " << substrings[peptide_column][0] << ".<A HREF=\"http://www.ncbi.nlm.nih.gov\">" << peptide << "</A>." << substrings[peptide_column][substrings[peptide_column].length() - 1] << endl;
			
			line_buffer = line_ss.str();
			
			for ( UnsignedInt prot = 0; prot < proteins_per_peptide; ++prot )
			{
				getline(out_file, line);
			}
		}
		if ( i == 1 )
		{
			line_buffer.replace(line_buffer.find("###"), 3, "1"); // the delta_cn for the the preceeding value is computed using this value
			summary << line_buffer;
		}
		out2summary_number += displayed_peptides;
		
		out_file.close();
		out_file.clear();
		summary.close();
		summary.clear();
	}
	
	void
	SequestOutfile::finishSummaryHtml(
		const string& summary_filename
	)
	throw(
		Exception::UnableToCreateFile)
	{
		ofstream summary(summary_filename.c_str(), ios::out | ios::app);
		if ( !summary )
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, summary_filename);
		}
		summary << "</BODY></HTML>" << endl;
		summary.close();
		summary.clear();
	}
	
	map< String, vector< Real > >
	SequestOutfile::getPeptidePValues(
		const string& out_dir,
		const string& prob_filename)
	throw(
		Exception::FileNotFound)
	{
		ifstream prob_file(prob_filename.c_str());
		if ( !prob_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, prob_filename);
		}
unsigned int n = 0;
		// the probability file has one line for each peptide
		// each line consists of filename\tprobability\tstuff/negonly
		String line, filename;
		vector< String > substrings;
		map< String, vector< Real > > filenames_and_pvalues;
		vector< Real >* pvalues(0);
		while ( getline(prob_file, line) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.split('\t', substrings);
			substrings[0].insert(0, out_dir);
			substrings[0].append(".out");
			if ( filename != substrings[0] ) // if a new filename is found, insert a vector and set the pointer accordingly
			{
				filename = substrings[0];
				filenames_and_pvalues[filename] = vector< Real >();
				pvalues = &filenames_and_pvalues[filename];
			}
			if ( substrings[2] == "negonly" ) pvalues->push_back(1.0);
			else
			{
				pvalues->push_back(1.0 - substrings[1].toFloat());
				if ( 1.0 - substrings[1].toFloat() <= 0.05 ) ++n;
			}
		}
		
		return filenames_and_pvalues;
	}
} //namespace OpenMS
