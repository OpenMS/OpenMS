// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
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
	
	void ControlledVocabulary::loadFromOBO(const String& name, const String& filename) throw (Exception::FileNotFound, Exception::ParseError)
	{
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
			
			if (line!="") //skip empty lines
			{
				if (line_wo_spaces.toLower()=="[term]") //new term
				{
					if (term.id!="") //store last term
					{
						terms_[term.id] = term;
						term.id="";
						term.name="";
						term.parents.clear();
						term.obsolete=false;
					}
					else if (term.name!="")
					{
						cout << "Warning: Dropping term without identifier!" << endl;
					}
				}
				//new id line
				else if (line_wo_spaces.hasPrefix("id:"))
				{
					term.id = line.substr(line.find(':')+1).trim();
				}
				else if (line_wo_spaces.hasPrefix("name:"))
				{
					term.name = line.substr(line.find(':')+1).trim();
				}
				else if (line_wo_spaces.hasPrefix("is_a:"))
				{
					term.parents.push_back(line.substr(line.find(':')+1).prefix('!').trim());
				}
				else if (line_wo_spaces=="is_obsolete:true")
				{
					term.obsolete = true;
				}
				else
				{
					//cout << "Ignored line: '" << line << "'" << endl;
				}
			}
		 }

		if (term.id!="") //store last term
		{
			terms_[term.id] = term;
		}
	}
		
	const ControlledVocabulary::CVTerm& ControlledVocabulary::getTerm(const String& id) const throw (Exception::InvalidValue)
	{
		map<String, CVTerm>::const_iterator it = terms_.find(id);
		if (it==terms_.end())
		{
			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Inalid CV identifier!",id);
		}
		return it->second;
	}	

	bool ControlledVocabulary::exists(const String& id) const
	{
		if (terms_.find(id)==terms_.end())
		{
			return false;
		}
		return true;
	}	
	
	bool ControlledVocabulary::isChildOf(const String& child, const String& parent) const throw (Exception::InvalidValue)
	{
		//cout << "CHECK child:" << child << " parent: " << parent << endl;
		const CVTerm& ch = getTerm(child);
		
		for (UInt i=0; i<ch.parents.size(); ++i)
		{
			//cout << "Parent: " << ch.parents[i] << endl;
			
			//check if it is a direct parent
			if (ch.parents[i]==parent)
			{
				return true;
			}
			//check if it is an indirect parent
			else if (isChildOf(ch.parents[i],parent))
			{
				return true;
			}
		}
		
		return false;
	}

	std::ostream& operator << (std::ostream& os, const ControlledVocabulary& cv)
	{
		for (map<String, ControlledVocabulary::CVTerm>::const_iterator it = cv.terms_.begin(); it!=cv.terms_.end(); ++it)
		{
			os << "[Term]" << endl;
			os << "id: '" << it->second.id << "'" <<endl;
			os << "name: '" << it->second.name <<  "'" << endl;
			for (vector<String>::const_iterator it2 = it->second.parents.begin(); it2!= it->second.parents.end(); ++it2)
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
	
} // namespace OpenMS

