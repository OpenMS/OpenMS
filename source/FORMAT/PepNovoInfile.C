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

#include <OpenMS/FORMAT/PepNovoInfile.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/PTMXMLFile.h>

#include <algorithm>
#include <set>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

namespace OpenMS
{

	PepNovoInfile::PepNovoInfile()
	{}
	
	PepNovoInfile::PepNovoInfile(const PepNovoInfile& pepnovo_infile)
	{
		PTMname_residues_mass_type_ = pepnovo_infile.getModifications();
	}

	PepNovoInfile::~PepNovoInfile()
	{
		PTMname_residues_mass_type_.clear();
	}
	
	PepNovoInfile& PepNovoInfile::operator=(const PepNovoInfile& pepnovo_infile)
	{
		if ( this != &pepnovo_infile ) PTMname_residues_mass_type_ = pepnovo_infile.getModifications();
		return *this;
	}
	
	bool PepNovoInfile::operator==(const PepNovoInfile& pepnovo_infile) const
	{
		if ( this != &pepnovo_infile )
		{
			return ( PTMname_residues_mass_type_ == pepnovo_infile.getModifications() );
		}
		return true;
	}

	String PepNovoInfile::store(const String& filename)
	{
		ofstream ofs(filename.c_str());
		if ( !ofs )
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
		stringstream file_content;
		String line, abbreviation_string;
		set< String > abbreviations;
		
		stringstream fixed_ptms, optional_ptms, cterminal_ptms, nterminal_ptms;
		stringstream* sstream_p = NULL;
		String residues, abbreviation;
		String rounded_mass;
		Int counter(-1);
		// first write the fixed ptms
		for ( map< String, vector< String > >::iterator mods_i = PTMname_residues_mass_type_.begin(); mods_i != PTMname_residues_mass_type_.end(); ++mods_i )
		{
			residues = mods_i->second[0];
			rounded_mass = mods_i->second[1].toFloat()< 0 ? String((Int) (mods_i->second[1].toFloat() - .5)) : "+" + String((Int) (mods_i->second[1].toFloat() + .5));
			if ( residues == "CTERM" )
			{
				sstream_p = &cterminal_ptms;
				residues = "C_TERM";
				mods_i->second[2] = "OPTIONAL";
			}
			else if ( residues == "NTERM" )
			{
				sstream_p = &nterminal_ptms;
				residues = "N_TERM";
				mods_i->second[2] = "OPTIONAL";
			}
			else if ( mods_i->second[2] == "FIX" )
			{
				sstream_p = &fixed_ptms;
				mods_i->second[2] = "FIXED";
			}
			else if ( mods_i->second[2] == "OPT" )
			{
				sstream_p = &optional_ptms;
				mods_i->second[2] = "OPTIONAL";
			}
			
			if ( residues == "C_TERM" || residues == "N_TERM" )
			{
				if ( residues == "C_TERM" ) abbreviation = "§" + rounded_mass;
				if ( residues == "N_TERM" ) abbreviation = "^" + rounded_mass;
				counter = 0;
				// make sure no abbreviation appears more than once
				while ( abbreviations.find(abbreviation + "_" + String(counter)) != abbreviations.end() ) ++counter;
				if ( counter ) abbreviation.append("_" + String(counter));
				abbreviations.insert(abbreviation);
				(*sstream_p) << residues << "   " << mods_i->second[1] << "   " << mods_i->second[2] << "   " + residues + "   " << abbreviation << "   " << mods_i->first << "\n";
			}
			else
			{
				for ( String::const_iterator residue_i = residues.begin(); residue_i != residues.end(); ++residue_i )
				{
					counter = 0;
					abbreviation = *residue_i + rounded_mass;
					// make sure no abbreviation appears more than once
					while ( abbreviations.find(abbreviation + "_" + String(counter)) != abbreviations.end() ) ++counter;
					if ( counter ) abbreviation.append("_" + String(counter));
					abbreviations.insert(abbreviation);
					(*sstream_p) << *residue_i << "   " << mods_i->second[1] << "   " << mods_i->second[2] << "   ALL   " << abbreviation << "   " << mods_i->first << "\n";
				}
			}
		}

		file_content << "# Fixed PTMs" << "\n" << "#AA  offset      type    locations  symbol  PTM name" << "\n" << fixed_ptms.str() << "\n" <<  "# Optional PTMs" << "\n" << optional_ptms.str() << "\n" << "# Terminal PTMs" << "\n" << cterminal_ptms.str() << nterminal_ptms.str();

		for ( set< String >::const_iterator abbreviation_i = abbreviations.begin(); abbreviation_i != abbreviations.end(); ++abbreviation_i )
		{
			abbreviation_string.append(*abbreviation_i + ":");
		}
		abbreviation_string.erase(abbreviation_string.length() - 1); // remove the last ':'
		
		ofs << file_content.str();
		ofs.close();
		ofs.clear();
		
		return abbreviation_string;
	}

	void PepNovoInfile::handlePTMs(const String& modification_line, const String& modifications_filename, const bool monoisotopic)
	{
		PTMname_residues_mass_type_.clear();
		// to store the information about modifications from the ptm xml file
		map< String, pair< String, String > > ptm_informations;
		if ( !modification_line.empty() ) // if modifications are used look whether whether composition and residues (and type and name) is given, the name (type) is used (then the modifications file is needed) or only the mass and residues (and type and name) is given
		{
			vector< String > modifications, mod_parts;
			modification_line.split(':', modifications); // get the single modifications
			
			// to get masses from a formula
			EmpiricalFormula add_formula, substract_formula;
			
			String types = "OPT#FIX#";
			String name, residues, mass, type;
			
			// 0 - mass; 1 - composition; 2 - ptm name
			Int mass_or_composition_or_name(-1);
			
			for ( vector< String >::const_iterator mod_i = modifications.begin(); mod_i != modifications.end(); ++mod_i )
			{
				if ( mod_i->empty() ) continue;
				// clear the formulae
				add_formula = substract_formula = name = residues = mass = type = "";
				
				// get the single parts of the modification string
				mod_i->split(',', mod_parts);
				mass_or_composition_or_name = -1;
				
				// check whether the first part is a mass, composition or name
				
				// check whether it is a mass
				try
				{
					mass = mod_parts.front();
					// to check whether the first part is a mass, it is converted into a float and then back into a string and compared to the given string
					// remove + signs because they don't appear in a float
					if ( mass.hasPrefix("+") ) mass.erase(0, 1);
					if ( mass.hasSuffix("+") ) mass.erase(mass.length() - 1, 1);
					if ( mass.hasSuffix("-") ) // a - sign at the end will not be converted
					{
						mass.erase(mass.length() - 1, 1);
						mass.insert(0, "-");
					}
					// if it is a mass
					if ( String(mass.toFloat()) == mass ) mass_or_composition_or_name = 0;
				}
				catch ( Exception::ConversionError c_e ){ mass_or_composition_or_name = -1; }
				
				// check whether it is a name (look it up in the corresponding file)
				if ( mass_or_composition_or_name == -1 )
				{
					if ( ptm_informations.empty() ) // if the ptm xml file has not been read yet, read it
					{
						if ( !File::exists(modifications_filename) )
						{
							throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, modifications_filename);
						}
						if ( !File::readable(modifications_filename) )
						{
							throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, modifications_filename);
						}
						
						// getting all available modifications from a file
						PTMXMLFile().load(modifications_filename, ptm_informations);
					}
					// if the modification cannot be found
					if ( ptm_informations.find(mod_parts.front()) != ptm_informations.end() )
					{
						mass = ptm_informations[mod_parts.front()].first; // composition
						residues = ptm_informations[mod_parts.front()].second; // residues
						name = mod_parts.front(); // name
						
						mass_or_composition_or_name = 2;
					}
				}
				
				// check whether it's an empirical formula / if a composition was given, get the mass
				if ( mass_or_composition_or_name == -1 ) mass = mod_parts.front();
				if ( mass_or_composition_or_name == -1 || mass_or_composition_or_name == 2 )
				{
					// check whether there is a positive and a negative formula
					String::size_type pos = mass.find("-");
					try
					{
						if ( pos != String::npos )
						{
							add_formula = mass.substr(0, pos);
							substract_formula = mass.substr(++pos);
						}
						else
						{
							add_formula = mass;
						}
						// sum up the masses
						if ( monoisotopic ) mass = String(add_formula.getMonoWeight() - substract_formula.getMonoWeight());
						else mass = String(add_formula.getAverageWeight() - substract_formula.getAverageWeight());
						if ( mass_or_composition_or_name == -1 ) mass_or_composition_or_name = 1;
					}
					catch ( Exception::ParseError pe )
					{
						PTMname_residues_mass_type_.clear();
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, *mod_i, "There's something wrong with this modification. Aborting!");
					}
				}
				
				// now get the residues
				mod_parts.erase(mod_parts.begin());
				if ( mass_or_composition_or_name < 2 )
				{
					if ( mod_parts.empty() )
					{
						PTMname_residues_mass_type_.clear();
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, *mod_i, "No residues for modification given. Aborting!");
					}
					
					// get the residues
					residues = mod_parts.front();
					residues.substitute('*', 'X');
					residues.toUpper();
					mod_parts.erase(mod_parts.begin());
				}
				
				// get the type
				if ( mod_parts.empty() ) type = "OPT";
				else
				{
					type = mod_parts.front();
					type.toUpper();
					if ( types.find(type) != String::npos ) mod_parts.erase(mod_parts.begin());
					else type = "OPT";
				}
				
				if ( mod_parts.size() > 1 )
				{
					PTMname_residues_mass_type_.clear();
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, *mod_i, "There's something wrong with the type of this modification. Aborting!");
				}
				
				// get the name
				if ( mass_or_composition_or_name < 2 )
				{
					if ( mod_parts.empty() ) name = "PTM_" + String(PTMname_residues_mass_type_.size());
					else name = mod_parts.front();
				}
				
				// insert the modification
				if ( PTMname_residues_mass_type_.find(name) == PTMname_residues_mass_type_.end() )
				{
					PTMname_residues_mass_type_[name] = vector< String >(3);
					PTMname_residues_mass_type_[name][0] = residues;
					// mass must not have more than 5 digits after the . (otherwise the test may fail)
					PTMname_residues_mass_type_[name][1] = mass.substr(0, mass.find(".") + 6);
					PTMname_residues_mass_type_[name][2] = type;
				}
				else
				{
					PTMname_residues_mass_type_.clear();
					throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, *mod_i, "There's already a modification with this name. Aborting!");
				}
			}
		}
	}

	const map< String, vector< String > >& PepNovoInfile::getModifications() const {return PTMname_residues_mass_type_;}
}
