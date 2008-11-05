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

#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>

using namespace xercesc;
using namespace std;

namespace OpenMS 
{
  namespace Internal
  {
		MzMLValidator::MzMLValidator(const CVMappings& mapping, const ControlledVocabulary& cv)
			: SemanticValidator(mapping, cv),
				cv_values_()
		{
			//find and store those cv terms, that can have a value
  		for (Map<String,ControlledVocabulary::CVTerm>::const_iterator it=cv_.getTerms().begin(); it!=cv_.getTerms().end(); ++it)
  		{
				for (UInt i=0; i<it->second.unparsed.size(); ++i)
				{
					if (it->second.unparsed[i].hasSubstring("value-type:xsd\\:int"))
					{
						cv_values_[it->first] = DataValue::INT_VALUE;
					}
					else if (it->second.unparsed[i].hasSubstring("value-type:xsd\\:float"))
					{
						cv_values_[it->first] = DataValue::DOUBLE_VALUE;
					}
					else if (it->second.unparsed[i].hasSubstring("value-type:xsd\\:string"))
					{
						cv_values_[it->first] = DataValue::STRING_VALUE;
					}
				}
  		}
		}
		
		MzMLValidator::~MzMLValidator()
		{
		}
		
		//This method needed to be reimplemented to
		// - check CV term values
		// - handle referenceableParamGroups
	  void MzMLValidator::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
	  {
	    String tag = sm_.convert(qname);
	    String parent_tag;
	    if (open_tags_.size()>0) parent_tag = open_tags_.back();
	    String path = getPath_() + "/" + cv_tag_ + "/@" + accession_att_;
	    open_tags_.push_back(tag);

	    if (tag=="referenceableParamGroup")
	    {
				current_id_ = attributeAsString_(attributes,"id");
	  	}
	  	else if (tag=="referenceableParamGroupRef")
	    {
	    	const std::vector<CVTerm>& terms = param_groups_[attributeAsString_(attributes,"ref")];
				for (UInt i=0; i< terms.size(); ++i)
				{
					handleTerm_(path, terms[i]);
				}
	  	}
	    else if (tag==cv_tag_)
	    {
		    //extract accession, name and value
		    CVTerm parsed_term;
		    getCVTerm_(attributes, parsed_term);
	
				//check if the term is unknown
				if (!cv_.exists(parsed_term.accession))
				{
					warnings_.push_back(String("Unknown CV term: '") + parsed_term.accession + " - " + parsed_term.name + "' at element '" + getPath_(1) + "'");
					return;
				}
				
				//check if the term is obsolete
				if (cv_.getTerm(parsed_term.accession).obsolete)
				{
					warnings_.push_back(String("Obsolete CV term: '") + parsed_term.accession + " - " + parsed_term.name + "' at element '" + getPath_(1) + "'");
				}
				
				//values used in wrong places and wrong value types
				if (parsed_term.value!="")
				{
					if (!cv_values_.has(parsed_term.accession))
					{
						errors_.push_back(String("CV term should not have a value: '") + parsed_term.accession + " - " + parsed_term.name + "' (value: '" + parsed_term.value + "') at element '" + getPath_(1) + "'");
	    		}
					else
					{
						if (cv_values_[parsed_term.accession]==DataValue::INT_VALUE)
						{
							try
							{
								parsed_term.value.toInt();
							}
							catch(Exception::ConversionError&)
							{
								errors_.push_back(String("CV term should have an integer value: '") + parsed_term.accession + " - " + parsed_term.name + "' (value: '" + parsed_term.value + "') at element '" + getPath_(1) + "'");
							}
						}
						else if (cv_values_[parsed_term.accession]==DataValue::DOUBLE_VALUE)
						{
							try
							{
								parsed_term.value.toDouble();
							}
							catch(Exception::ConversionError&)
							{
								errors_.push_back(String("CV term should have a floating-point value: '") + parsed_term.accession + " - " + parsed_term.name + "' (value: '" + parsed_term.value + "') at element '" + getPath_(1) + "'");
							}
						}
					}
				}
				//no value, although there should be a numerical value
				else if (cv_values_.has(parsed_term.accession) && (cv_values_[parsed_term.accession]==DataValue::INT_VALUE || cv_values_[parsed_term.accession]==DataValue::DOUBLE_VALUE))
				{
					errors_.push_back(String("CV term should have a numerical value: '") + parsed_term.accession + " - " + parsed_term.name + "' (value: '" + parsed_term.value + "') at element '" + getPath_(1) + "'");
    		}
				
				//actual handling of the term
	    	if (parent_tag=="referenceableParamGroup")
	    	{
	    		param_groups_[current_id_].push_back(parsed_term);
	    	}
	    	else
	    	{				
					handleTerm_(path, parsed_term);
				}
			}
	  }
		
		//reimplemented in order to remove the "indexmzML" tag from the front (if present)
		String MzMLValidator::getPath_(UInt remove_from_end) const
		{
			String path;
			if (open_tags_.size()!=0 && open_tags_.front()=="indexedmzML")
			{
				path.implode(open_tags_.begin()+1, open_tags_.end()-remove_from_end,"/");
			}
			else
			{
				path.implode(open_tags_.begin(), open_tags_.end()-remove_from_end,"/");
			}
			path = String("/") + path;
			return path;
		}

	} // namespace Internal
} // namespace OpenMS
