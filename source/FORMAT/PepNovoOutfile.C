// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/FORMAT/PepNovoOutfile.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <fstream>
#include <limits>

using namespace std;

namespace OpenMS
{
	PepNovoOutfile::PepNovoOutfile() {}

	PepNovoOutfile::PepNovoOutfile(const PepNovoOutfile&) {}

	PepNovoOutfile::~PepNovoOutfile() {}

	PepNovoOutfile& PepNovoOutfile::operator=(const PepNovoOutfile&)
	{
		return *this;
	}

	bool PepNovoOutfile::operator==(const PepNovoOutfile&) const
	{
		return true;
	}

	void
	PepNovoOutfile::load(
		const string& result_filename,
		vector< PeptideIdentification >&	peptide_identifications,
		ProteinIdentification& protein_identification,
		const Real& score_threshold,
		const map< String, Real >& rt_and_index)
	{
		// generally used variables
		vector< String > substrings;
		map< String, Int > columns;
		PeptideHit peptide_hit;

		String
			line,
			score_type,
			version,
			identifier,
			filename,
			sequence,
			sequence_with_mods;

		DateTime datetime = DateTime::now(); // there's no date given from PepNovo
		protein_identification.setDateTime(datetime);

		peptide_identifications.clear();
		PeptideIdentification peptide_identification;
		protein_identification = ProteinIdentification();

		// open the result
		ifstream result_file(result_filename.c_str());
		if ( !result_file )
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, result_filename);
		}

		UInt line_number(0); // used to report in which line an error occured

		score_type = "PepNovo3";
		version = "v3.00";

		if ( protein_identification.getSearchEngineVersion().empty() )
		{
			protein_identification.setSearchEngine("PepNovo");
			protein_identification.setSearchEngineVersion(version);
		}
		identifier = protein_identification.getSearchEngine() + "_" + datetime.getDate();
		protein_identification.setIdentifier(identifier);

		while ( getline(result_file, line) )
		{
			if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
			line.trim();
			++line_number;
			if ( line.hasPrefix(">> ") ) // >> 1 /home/shared/pepnovo/4611_raw_ms2_picked.mzXML.1001.2.dta
			{
				if ( !peptide_identifications.empty() || !peptide_identification.getHits().empty() ) peptide_identifications.push_back(peptide_identification);
				peptide_identification = PeptideIdentification();
				
				String index = File::basename(line.substr(line.find(' ', strlen(">> ")) + 1));
				if ( rt_and_index.find(index) != rt_and_index.end() ) peptide_identification.setMetaValue("RT",  rt_and_index.find(index)->second);
				else peptide_identification.setMetaValue("RT", 0);
				
				peptide_identification.setSignificanceThreshold(score_threshold);
				peptide_identification.setScoreType(score_type);
				peptide_identification.setIdentifier(identifier);
			}
			else if ( line.hasPrefix("#Index") ) // #Index  Prob    Score   N-mass  C-Mass  [M+H]   Charge  Sequence
			{
				if ( columns.empty() ) // map the column names to their column number
				{
					line.split('\t', substrings);
					for ( vector< String >::const_iterator s_i = substrings.begin(); s_i != substrings.end(); ++s_i )
					{
						if ( (*s_i) == "#Index" ) columns["Index"] = s_i - substrings.begin();
						else if ( (*s_i) == "RnkScr" ) columns["RnkScr"] = s_i - substrings.begin();
						else if ( (*s_i) == "PnvScr" ) columns["PnvScr"] = s_i - substrings.begin();
						else if ( (*s_i) == "N-Gap" ) columns["N-Gap"] = s_i - substrings.begin();
						else if ( (*s_i) == "C-Gap" ) columns["C-Gap"] = s_i - substrings.begin();
						else if ( (*s_i) == "[M+H]" ) columns["[M+H]"] = s_i - substrings.begin();
						else if ( (*s_i) == "Charge" ) columns["Charge"] = s_i - substrings.begin();
						else if ( (*s_i) == "Sequence" ) columns["Sequence"] = s_i - substrings.begin();
					}

					if ( columns.size() != 8 )
					{
						result_file.close();
						result_file.clear();
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Not enough columns in file in line " + String(line_number) + String(" (should be 8)!"), result_filename);
					}
				}

				while ( getline(result_file, line) )
				{
					++line_number;
					if ( !line.empty() && (line[line.length()-1] < 33) ) line.resize(line.length()-1);
					line.trim();

					if ( line.empty() ) break;

					line.split('\t', substrings);
					if ( !substrings.empty() )
					{
						if ( substrings.size() != 8 )
						{
							result_file.close();
							result_file.clear();
							throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Not enough columns in file in line " + String(line_number) + String(" (should be 8)!"), result_filename);
						}
						if ( 1 - substrings[columns["Prob"]].toFloat() <= score_threshold + std::numeric_limits<Real>::epsilon() )
						{
							peptide_hit = PeptideHit();
							peptide_hit.setCharge(substrings[columns["Charge"]].toInt());
							peptide_hit.setRank(substrings[columns["Index"]].toInt() + 1);
							peptide_hit.setScore(substrings[columns["PnvScr"]].toFloat());
							peptide_hit.setMetaValue("RnkScr", substrings[columns["RnkScr"]].toFloat());
							peptide_hit.setMetaValue("N-Gap", substrings[columns["N-Gap"]].toFloat());
							peptide_hit.setMetaValue("C-Gap", substrings[columns["C-Gap"]].toFloat());
							//remove modifications (small characters and anything that's not in the alphabet)
							sequence_with_mods = substrings[columns["Sequence"]];
							sequence.clear();
							for ( String::ConstIterator c_i = sequence_with_mods.begin(); c_i != sequence_with_mods.end(); ++c_i )
							{
								if ( (bool) isalpha(*c_i) && (bool) isupper(*c_i) ) sequence.append(1, *c_i);
							}
							peptide_hit.setSequence(sequence);
							peptide_identification.insertHit(peptide_hit);
						}
					}
				}
				peptide_identification.setMetaValue("MZ", substrings[columns["[M+H]"]].toFloat());
			}
		}
		if ( !peptide_identifications.empty() || !peptide_identification.getHits().empty() ) peptide_identifications.push_back(peptide_identification);
		
		result_file.close();
		result_file.clear();
	}

	void
	PepNovoOutfile::getSearchEngineAndVersion(
		const String& pepnovo_output_without_parameters_filename,
		ProteinIdentification& protein_identification)
	{
		ifstream pepnovo_output_without_parameters(pepnovo_output_without_parameters_filename.c_str());
		if (!pepnovo_output_without_parameters)
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, pepnovo_output_without_parameters_filename);
		}

		// searching for something like this: PepNovo v1.03
		String line;
		vector< String > substrings;
		while (getline(pepnovo_output_without_parameters, line))
		{
			if (!line.empty() && (line[line.length()-1] < 33)) line.resize(line.length() - 1);
			line.trim();
			if ( line.empty() ) continue;
			if ( !line.hasPrefix("PepNovo") ) continue;
			else
			{
				line.split(' ', substrings);
				if ( substrings.empty() ) substrings.push_back(line);
				protein_identification.setSearchEngine(substrings[0]);
				if ( substrings.size() >= 2 ) protein_identification.setSearchEngineVersion(substrings[1]);
				return;
			}
		}
	}

} //namespace OpenMS
