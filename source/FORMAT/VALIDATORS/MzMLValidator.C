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
// $Maintainer: Marc Sturm $
// $Authors: $
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
				binary_data_array_(),
				binary_data_type_()
		{
			setCheckUnits(true);
		}
		
		MzMLValidator::~MzMLValidator()
		{
		}
		
		//This method needed to be reimplemented to
		// - check CV term values
		// - handle referenceableParamGroups
		// - check if binaryDataArray name and type match
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
				for (Size i=0; i< terms.size(); ++i)
				{
					handleTerm_(path, terms[i]);
				}
	  	}
	  	else if (tag=="binaryDataArray")
	    {
				binary_data_array_ = "";
				binary_data_type_ = "";
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
				path.concatenate(open_tags_.begin()+1, open_tags_.end()-remove_from_end,"/");
			}
			else
			{
				path.concatenate(open_tags_.begin(), open_tags_.end()-remove_from_end,"/");
			}
			path = String("/") + path;
			return path;
		}
		
		//reimplemented to 
		// - catch non-PSI CVs
		// - check if binaryDataArray name and type match
		void MzMLValidator::handleTerm_(const String& path, const CVTerm& parsed_term) 
		{
			//some CVs cannot be validated because they use 'part_of' which spoils the inheritance
			if (parsed_term.accession.hasPrefix("GO:")) return;
			if (parsed_term.accession.hasPrefix("BTO:")) return;
			
			//check binary data array terms
			if (path.hasSuffix("/binaryDataArray/cvParam/@accession"))
			{
				//binary data array
				if (cv_.isChildOf(parsed_term.accession,"MS:1000513"))
				{
					binary_data_array_ = parsed_term.accession;
				}
				//binary data type
				if (cv_.isChildOf(parsed_term.accession,"MS:1000518"))
				{
					binary_data_type_ = parsed_term.accession;
				}
				//if both are parsed, check if they match
				if (binary_data_type_!="" && binary_data_array_!="")
				{
					if (!cv_.getTerm(binary_data_array_).xref_binary.contains(binary_data_type_))
					{
						errors_.push_back(String("Binary data array of type '") + binary_data_array_ + " ! " + cv_.getTerm(binary_data_array_).name + "' cannot have the value type '" + binary_data_type_ + " ! " + cv_.getTerm(binary_data_type_).name + "'.");
					}
				}
			}
			
			SemanticValidator::handleTerm_(path,parsed_term);
		}

	} // namespace Internal
} // namespace OpenMS
