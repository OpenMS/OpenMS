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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XTandemXMLHandler.h>
#include <xercesc/sax2/Attributes.hpp>

using namespace std;
using namespace xercesc;

// TODO test all indices if >= 0!!!

namespace OpenMS
{
	namespace Internal
	{
  
  XTandemXMLHandler::XTandemXMLHandler(ProteinIdentification& protein_id,
								  									 map<UInt, vector<PeptideHit> >& peptide_hits, 
																		 map<UInt, String>& descriptions,
      								 							 const String& filename) :
    XMLHandler(filename),
    protein_id_(protein_id),
    peptide_hits_(peptide_hits),
		descriptions_(descriptions),
		actual_protein_id_(),
		actual_charge_(0),
		tag_()
  {
  	
  }
   
  XTandemXMLHandler::~XTandemXMLHandler()
  {
    
  }
  
  void XTandemXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
	{

		tag_ = String(sm_.convert(qname));
		//cerr << "start_local_name: " << String(sm_.convert(local_name)) << endl;
		//cerr << "start_uri: " << String(sm_.convert(uri)) << endl;

		//cerr << "startElement tag_: " << tag_ << endl;

		
		if (tag_ == "domain")
		{
			PeptideHit hit;
			hit.metaRegistry().registerName("E-Value", "E-Value of Hit");

			// get hyperscore
			double hyperscore(String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("hyperscore"))))).toDouble());
			hit.setScore(hyperscore);

			// get sequence of peptide
			String seq(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("seq")))));
			hit.setSequence(seq);

			// get amino acid before
			String pre(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("pre")))));
			if (pre.size() != 0)
			{
				hit.setAABefore(pre[pre.size() - 1]);
			}

			// get amino acid after
			String post(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("post")))));
			if (post.size() != 0)
			{
				hit.setAAAfter(post[0]);
			}
			
			// get expectation value
			double expect(String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("expect"))))).toDouble());
			hit.setMetaValue("E-Value", expect);

			// get precursor m/z
			double mh(String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("mh"))))).toDouble());
			hit.setMetaValue("MZ", mh);
			
			// spectrum id
			String id_string(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("id")))));
			vector<String> split;
			id_string.split('.', split);
			UInt id(split[0].toInt());
			actual_id_ = id;
			
			// add the actual protein accession
			hit.addProteinAccession(actual_protein_id_);
			hit.setCharge(actual_charge_);
			
			peptide_hits_[id].push_back(hit);
			return;
		}


		if (tag_ == "group")
		{
			Int index = attributes.getIndex(sm_.convert("z"));
			if (index >= 0)
			{
				actual_charge_ = String(sm_.convert(attributes.getValue(index))).toInt();
			}
		}

		/*
		if (tag_ == "note")
		{
			Int index = attributes.getIndex(sm_.convert("label"));
			if (index >= 0)
			{
				String label = String(sm_.convert(attributes.getValue(index)));
				if (label == "Description")
				{
					descriptions_[actual_id_] = ((String) sm_.convert(chars)).trim();
				}
			}
		}
		*/

		if (tag_ == "protein")
		{
			ProteinHit hit;
			String accession(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("label")))));
			actual_protein_id_ = accession;

			if (accessions_.find(accession) == accessions_.end())
			{
				accessions_.insert(accession);
				hit.setAccession(accession);

				double score(String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("expect"))))).toDouble());
				hit.setScore(score);

				//actual_protein_id_ = (UInt)String(sm_.convert(attributes.getValue(attributes.getIndex(sm_.convert("id"))))).toInt();

				protein_id_.insertHit(hit);
			}
		}

	}
	  
  void XTandemXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const /*qname*/)
 	{
		tag_ = "";
 	} 

  void XTandemXMLHandler::characters(const XMLCh* const chars, const unsigned int /*length*/)
  {
   	if (tag_ == "note")
	 	{
			String description = ((String) sm_.convert(chars)).trim();
			if (description.size() == 34)
			{
      	descriptions_[actual_id_] = description;
			}
		}
	}

	} // namespace Internal
} // namespace OpenMS
