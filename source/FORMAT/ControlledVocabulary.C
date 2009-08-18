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
// $Maintainer: Marc Sturm, Andreas Bertsch  $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ControlledVocabulary.h>

#include <fstream>
#include <iostream>
#include <utility>
#include <algorithm>

using namespace std;

namespace OpenMS 
{
	
	ControlledVocabulary::ControlledVocabulary()
		: terms_(),
			name_("")
	{
		
	}

	ControlledVocabulary::~ControlledVocabulary()
	{
		
	}
	
	void ControlledVocabulary::loadFromOBO(const String& name, const String& filename)
	{
		bool in_term = false;
		name_ = name;
		
		ifstream is(filename.c_str());
    if (!is)
    {
     	throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
		
		String line, line_wo_spaces;
		CVTerm term;
	    
    //parse file
    while(getline(is,line,'\n'))
    {
			line.trim();
			line_wo_spaces = line;
			line_wo_spaces.removeWhitespaces();
			
			//do nothing for empty lines
			if (line=="") continue;
			
			//********************************************************************************
			//stanza line
			if (line_wo_spaces[0]=='[')
			{
				//[term] stanza
				if (line_wo_spaces.toLower()=="[term]") //new term
				{
					in_term = true;
					if (term.id!="") //store last term
					{
						terms_[term.id] = term;
					}
	
					//clear temporary term members
					term = CVTerm();
				}
				// other stanza => not in a term
				else
				{
					in_term = false;
				}
			}
			//********************************************************************************
			//data line
			else if (in_term)
			{
				if (line_wo_spaces.hasPrefix("id:"))
				{
					term.id = line.substr(line.find(':')+1).trim();
				}
				else if (line_wo_spaces.hasPrefix("name:"))
				{
					term.name = line.substr(line.find(':')+1).trim();
				}
				else if (line_wo_spaces.hasPrefix("is_a:"))
				{
					if (line.has('!'))
					{
						String parent_id = line.substr(line.find(':') + 1).prefix('!').trim();
						term.parents.insert(parent_id);

						//check if the parent term name is correct
						String parent_name = line.suffix('!').trim();
						if (!checkName_(parent_id,parent_name)) cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': parent term name '" << parent_name << "' and id '" << parent_id << "' differ." << endl;
					}
					else
					{
						term.parents.insert(line.substr(line.find(':') + 1).trim());
					}
				}
				// brenda tissue special relationships, DRV (derived and part of)
				else if (line_wo_spaces.hasPrefix("relationship:DRV") && name == "brenda")
				{
					if (line.has('!'))
					{
						// e.g. relationship: DRV BTO:0000142 ! brain
						String parent_id = line.substr(line.find("DRV") + 4).prefix(':') + ":" + line.suffix(':').prefix('!').trim();
						term.parents.insert(parent_id);

						//check if the parent term name is correct
						String parent_name = line.suffix('!').trim();
						if (!checkName_(parent_id,parent_name)) cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': DRV relationship term name '" << parent_name << "' and id '" << parent_id << "' differ." << endl;
					}
					else
					{
						// e.g. relationship: DRV BTO:0000142
						term.parents.insert(line.substr(line.find("DRV") + 4).prefix(':') + ":" + line.suffix(':').trim());
					}
				}
				else if (line_wo_spaces.hasPrefix("relationship:part_of") && name == "brenda")
				{
					if (line.has('!'))
					{
						String parent_id = line.substr(line.find("part_of") + 8).prefix(':') + ":" + line.suffix(':').prefix('!').trim();
						term.parents.insert(parent_id);

						//check if the parent term name is correct
						String parent_name = line.suffix('!').trim();
						if (!checkName_(parent_id,parent_name)) cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': part_of relationship term name '" << parent_name << "' and id '" << parent_id << "' differ." << endl;
					}
					else
					{
						 term.parents.insert(line.substr(line.find("part_of") + 8).prefix(':') + ":" + line.suffix(':').trim());
					}
				}
				else if (line_wo_spaces.hasPrefix("relationship:has_units"))
				{
					if (line.has('!'))
					{
						String unit_id = line.substr(line.find("has_units")+10).prefix(':') + ":" + line.suffix(':').prefix('!').trim();
						term.units.insert(unit_id);

						//check if the parent term name is correct
						String unit_name = line.suffix('!').trim();
						if (!checkName_(unit_id,unit_name)) cerr << "Warning: while loading term '" << term.id << "' of CV '" << name_ << "': has_units relationship term name '" << unit_name << "' and id '" << unit_id << "' differ." << endl;
					}
					else
					{
						term.units.insert(line.substr(line.find("has_units")+10).prefix(':') + ":" + line.suffix(':').trim());
					}
				}
				else if (line_wo_spaces.hasPrefix("def:"))
				{
					String description = line.substr(line.find('"')+1);
					description.trim();
					description = description.substr(0,description.find('"'));
					description.trim();
					term.description = description;
				}
				else if (line_wo_spaces.hasPrefix("synonym:"))
				{
					String synonym = line.substr(line.find('"')+1);
					synonym.trim();
					synonym = synonym.substr(0,synonym.find('"'));
					synonym.trim();
					term.synonyms.push_back(synonym);
				}
				else if (line_wo_spaces=="is_obsolete:true")
				{
					term.obsolete = true;
				}
				else if (line_wo_spaces.hasPrefix("xref:value-type") || line_wo_spaces.hasPrefix("xref_analog:value-type"))
				{
					line_wo_spaces.remove('\\');
          if (line_wo_spaces.hasSubstring("value-type:xsd:string")) 
					{ 
						term.xref_type = CVTerm::XSD_STRING; 
						continue; 
					}
          if (line_wo_spaces.hasSubstring("value-type:xsd:integer") || line_wo_spaces.hasSubstring("value-type:xsd:int")) 
					{ 
						term.xref_type = CVTerm::XSD_INTEGER; 
						continue; 
					}
          if (line_wo_spaces.hasSubstring("value-type:xsd:decimal") || 
							line_wo_spaces.hasSubstring("value-type:xsd:float") ||
							line_wo_spaces.hasSubstring("value-type:xsd:double")) 
					{ 
						term.xref_type = CVTerm::XSD_DECIMAL; 
						continue; 
					}
          if (line_wo_spaces.hasSubstring("value-type:xsd:negativeInteger")) 
					{ 
						term.xref_type = CVTerm::XSD_NEGATIVE_INTEGER; 
						continue; 
					}
          if (line_wo_spaces.hasSubstring("value-type:xsd:positiveInteger")) 
					{ 
						term.xref_type = CVTerm::XSD_POSITIVE_INTEGER; 
						continue; 
					}
          if (line_wo_spaces.hasSubstring("value-type:xsd:nonNegativeInteger")) 
					{ 
						term.xref_type = CVTerm::XSD_NON_NEGATIVE_INTEGER; 
						continue; 
					}
          if (line_wo_spaces.hasSubstring("value-type:xsd:nonPositiveInteger")) 
					{ 
						term.xref_type = CVTerm::XSD_NON_POSITIVE_INTEGER; 
						continue; 
					}
          if (line_wo_spaces.hasSubstring("value-type:xsd:boolean") || line_wo_spaces.hasSubstring("value-type:xsd:bool")) 
					{ 
						term.xref_type = CVTerm::XSD_BOOLEAN; 
						continue; 
					}
          if (line_wo_spaces.hasSubstring("value-type:xsd:date"))
					{ 
						term.xref_type = CVTerm::XSD_DATE; 
						continue; 
					}
					if (line_wo_spaces.hasSubstring("value-type:xsd:anyURI"))
					{
						term.xref_type = CVTerm::XSD_ANYURI;
						continue;
					}
					cerr << "ControlledVocabulary: OBOFile: unknown xsd type: " << line_wo_spaces << ", ignoring" << endl;
				}
				else if (line_wo_spaces.hasPrefix("xref:binary-data-type") || line_wo_spaces.hasPrefix("xref_analog:binary-data-type"))
				{
					line_wo_spaces.remove('\\');
					//remove description (if present)
					if (line_wo_spaces.has('!'))
					{
						line_wo_spaces = line_wo_spaces.substr(0,line_wo_spaces.find('!'));
					}
					//trim prefix
					line_wo_spaces = line_wo_spaces.substr(22);
					//trim just to be sure
					line_wo_spaces.trim();
					term.xref_binary.push_back(line_wo_spaces);
				}
				else if (line!="")
				{
					term.unparsed.push_back(line);
				}
			}
		}

		if (term.id!="") //store last term
		{
			terms_[term.id] = term;
		}

		// now build all child terms
		for (Map<String, CVTerm>::iterator it = terms_.begin(); it != terms_.end(); ++it)
		{
			//cerr << it->first << endl;
			for (set<String>::const_iterator pit = terms_[it->first].parents.begin(); pit != terms_[it->first].parents.end(); ++pit)
			{
				//cerr << "Parent: " << *pit << endl;
				terms_[*pit].children.insert(it->first);
			}
		}
	}
		
	const ControlledVocabulary::CVTerm& ControlledVocabulary::getTerm(const String& id) const
	{
		Map<String, CVTerm>::const_iterator it = terms_.find(id);
		if (it==terms_.end())
		{
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid CV identifier!",id);
		}
		return it->second;
	}	

	const Map<String, ControlledVocabulary::CVTerm>& ControlledVocabulary::getTerms() const
	{
		return terms_;
	}

	void ControlledVocabulary::getAllChildTerms(set<String>& terms, const String& parent) const
	{
		//cerr << "Parent: " << parent << endl;
		const set<String>& children = getTerm(parent).children;
		for (set<String>::const_iterator it = children.begin(); it != children.end(); ++it)
		{
			terms.insert(*it);
			//TODO This is not save for cyclic graphs. Are they allowd in CVs?
			getAllChildTerms(terms, *it);
		}
	}

	bool ControlledVocabulary::exists(const String& id) const
	{
		return terms_.has(id);
	}	
	
	bool ControlledVocabulary::isChildOf(const String& child, const String& parent) const
	{
		//cout << "CHECK child:" << child << " parent: " << parent << endl;
		const CVTerm& ch = getTerm(child);
		
		for (set<String>::const_iterator it = ch.parents.begin(); it != ch.parents.end(); ++it)
		{
			//cout << "Parent: " << ch.parents[i] << endl;
			
			//check if it is a direct parent
			if (*it == parent)
			{
				return true;
			}
			//check if it is an indirect parent
			else if (isChildOf(*it, parent))
			{
				return true;
			}
		}
		
		return false;
	}

	std::ostream& operator << (std::ostream& os, const ControlledVocabulary& cv)
	{
		for (Map<String, ControlledVocabulary::CVTerm>::const_iterator it = cv.terms_.begin(); it!=cv.terms_.end(); ++it)
		{
			os << "[Term]" << endl;
			os << "id: '" << it->second.id << "'" <<endl;
			os << "name: '" << it->second.name <<  "'" << endl;
			for (set<String>::const_iterator it2 = it->second.parents.begin(); it2!= it->second.parents.end(); ++it2)
			{
				cout << "is_a: '" << *it2 <<  "'" <<endl;
			}
		}
		return os;
	}
	
	const String& ControlledVocabulary::name() const
	{
		return name_;
	}
	
	bool ControlledVocabulary::checkName_(const String& id, const String& name, bool ignore_case)
	{
		if (!exists(id)) return true;
		
		String parent_name = name;
		String real_parent_name = getTerm(id).name;
		
		if (ignore_case)
		{
			parent_name.toLower();
			real_parent_name.toLower();
		}	

		return (real_parent_name==parent_name);
	}
	
} // namespace OpenMS

