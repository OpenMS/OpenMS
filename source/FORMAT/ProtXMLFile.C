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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ProtXMLFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>
#include <fstream>
#include <limits>

using namespace std;

namespace OpenMS 
{

	ProtXMLFile::ProtXMLFile()
		: XMLHandler("","1.2"),
			XMLFile("/SCHEMAS/protXML_v4.xsd","4.0"),
			last_meta_(0),
			document_id_()
	{
	}

  void ProtXMLFile::load(const String& filename,  vector<ProteinIdentification>& protein_ids, vector<PeptideIdentification>& peptide_ids)
  {
  	//Filename for error messages in XMLHandler
  	file_ = filename;
  	
  	protein_ids.clear();
  	peptide_ids.clear();
  	
  	prot_ids_ = &protein_ids;
  	pep_ids_ = &peptide_ids;
		
		parse_(filename,this);
    
		resetMembers_();
  }
  					 
  void ProtXMLFile::store(const String& /*filename*/, const vector<ProteinIdentification>& /*protein_ids*/, const vector<PeptideIdentification>& /*peptide_ids*/, const String& /*document_id*/)
  {
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		// resetMembers_();
  }


	/// reset members
	void ProtXMLFile::resetMembers_()
	{
    prot_ids_ = 0;
		pep_ids_ = 0;
		last_meta_ = 0;
		parameters_.clear();
		param_ = ProteinIdentification::SearchParameters();
		String id_ = "";
		prot_id_ = ProteinIdentification();
		pep_id_ = PeptideIdentification();
		prot_hit_ = ProteinHit();
		pep_hit_ = PeptideHit();
		proteinid_to_accession_.clear();
	}
  
	void ProtXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{		
		String tag = sm_.convert(qname);
		
		if (tag =="protein_group")
		{
			p_protein_group_ =  attributeAsDouble_(attributes,"probability");
		}
		else if (tag =="protein")
		{
			protein_p_ =  attributeAsDouble_(attributes,"probability");
			protein_name_ = attributeAsString_(attributes,"protein_name");
		}
		else if (tag =="peptide")
		{
			peptide_seq_ = attributeAsString_(attributes,"peptide_sequence");
			peptide_nspp_ = attributeAsDouble_(attributes,"nsp_adjusted_probability");
		}

	}
	
	void ProtXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		String tag = sm_.convert(qname);
		
		//START
		if (tag =="ProtXML")
		{
		}

	}

} // namespace OpenMS
