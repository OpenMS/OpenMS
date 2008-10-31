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

#include <OpenMS/FORMAT/SemanticValidator.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/CVMappings.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

using namespace xercesc;
using namespace std;

  
//TODO:
// - use of term in wrong location
// - missing terms (MUST, warning for SHOULD)
// - wrong combination of terms (XOR)
// - wrong number of repeats
//   - (isRepeatable)
//   - the same term twice

namespace OpenMS 
{

	SemanticValidator::SemanticValidator(const CVMappings& mapping, const ControlledVocabulary& cv)
		: XMLHandler("", 0),
			XMLFile(),
			mapping_(mapping),
			cv_(cv),
			open_tags_(),
			valid_(true),
			cv_tag_("cvParam"),
			accession_att_("accession"),
			name_att_("name"),
			value_att_("value")
	{
	}
	
	SemanticValidator::~SemanticValidator()
	{
	}

	void SemanticValidator::setTag(const String& tag)
	{
		cv_tag_ = tag;
	}
	
	void SemanticValidator::setAccessionAttribute(const String& accession)
	{
		accession_att_ = accession;
	}
	
	void SemanticValidator::setNameAttribute(const String& name)
	{
		name_att_ = name;
	}
	
	void SemanticValidator::setValueAttribute(const String& value)
	{
		value_att_ = value;
	}

  bool SemanticValidator::validate(const String& filename, ValidationOutput& output)
  {
		//try to open file
		if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
    
    //initialize
    output = output_ = ValidationOutput();
    valid_ = true;
    
    //parse
		file_ = filename;
		parse_(filename, this);
		
		//set output
		if (!valid_) output = output_;
		
		return valid_;
  }

  void SemanticValidator::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
  {
    String tag = sm_.convert(qname);
    
    //TODO Check if all required CVs are loaded => Excpetion if not
    
    if (tag==cv_tag_)
    {
    	//extract path
	    String path;
	    path.implode(open_tags_.begin(), open_tags_.end(),"/");
	    path = "/"+path;
	    
	    //extract accession, name and value
	    String accession = attributeAsString_(attributes, accession_att_.c_str());
	    String name = attributeAsString_(attributes, name_att_.c_str());
	    String value;
	    optionalAttributeAsString_(value, attributes, value_att_.c_str());
			
			//construct location
			ValiationLocation location;
			location.path=path;
			location.accession=accession;
			location.name=name;
			location.value=value;
			
			//check if the term is unknown
			if (!cv_.exists(accession))
			{
				valid_ = false;
				output_.unknown_terms.push_back(location);
			}
			else
			{
				const ControlledVocabulary::CVTerm& term = cv_.getTerm(accession);
				if (term.obsolete)
				{
					valid_ = false;
					output_.obsolete_terms.push_back(location);
				}
			}
			
    }
    
    open_tags_.push_back(tag);
  }
	

	void SemanticValidator::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    String tag = sm_.convert(qname);
    
    open_tags_.pop_back();
  }

  void SemanticValidator::characters(const XMLCh* const /*chars*/, const unsigned int /*length*/)
  {
  	//nothing to do here
  }
  					 
} // namespace OpenMS
