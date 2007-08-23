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

#include <OpenMS/FORMAT/HANDLERS/OMSSAXMLHandler.h>
#include <xercesc/sax2/Attributes.hpp>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
	namespace Internal
	{
  
  OMSSAXMLHandler::OMSSAXMLHandler(ProteinIdentification& protein_identification,
								  									 vector<PeptideIdentification>& peptide_identifications, 
      								 							 const String& filename) :
    XMLHandler(filename),
    protein_identification_(protein_identification),
    peptide_identifications_(peptide_identifications),
    actual_peptide_hit_(),
		actual_peptide_id_(),
		tag_()
  {
  	
  }
   
  OMSSAXMLHandler::~OMSSAXMLHandler()
  {
    
  }
  
  void OMSSAXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& /*attributes*/)
	{
		tag_ = String(sm_.convert(qname)).trim();

		if (tag_ == "MSHitSet_number")
		{
			return;
		}

		if (tag_ == "MSPepHit")
		{
			return;
		}
		if (tag_ == "MSHits")
		{
			return;
		}

		if (tag_ == "MSHitSet")
		{
			return;
		}
	}
	  
  void OMSSAXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
 	{
 		tag_ = String(sm_.convert(qname)).trim();

		// protein hits (MSPepHits) are handled in ::characters(...)

		// end of peptide hit 
 		if (tag_ == "MSHits")
		{
			actual_peptide_id_.insertHit(actual_peptide_hit_);
			actual_peptide_hit_ = PeptideHit();
			tag_ = "";
			return;
		}

		// end of peptide id
		if (tag_ == "MSHitSet")
		{
			peptide_identifications_.push_back(actual_peptide_id_);
			actual_peptide_id_.assignRanks();
			actual_peptide_id_ = PeptideIdentification();
		}

		tag_ = "";
 	} 

  void OMSSAXMLHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
  {
		String value = ((String)sm_.convert(chars)).trim();
		// MSPepHit section
    // <MSPepHit_start>0</MSPepHit_start>
    // <MSPepHit_stop>8</MSPepHit_stop>
    // <MSPepHit_accession>6599</MSPepHit_accession>
    // <MSPepHit_defline>CRHU2 carbonate dehydratase (EC 4.2.1.1) II [validated] - human</MSPepHit_defline>
    // <MSPepHit_protlength>260</MSPepHit_protlength>
    // <MSPepHit_oid>6599</MSPepHit_oid>
		if (tag_ == "MSPepHit_start")
		{
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_stop")
		{
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_accession")
		{
      actual_peptide_hit_.addProteinAccession(value);
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_defline")
		{
			// TODO add defline to ProteinHit?
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_protlength")
		{
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_oid")
		{
      tag_ = "";
      return;
		}

		// MSHits section
		// <MSHits_evalue>0.00336753988893542</MSHits_evalue>
		// <MSHits_pvalue>1.30819399070598e-08</MSHits_pvalue>
		// <MSHits_charge>1</MSHits_charge>
		// <MSHits_pepstring>MSHHWGYGK</MSHits_pepstring>
    // <MSHits_mass>1101492</MSHits_mass>
    // <MSHits_pepstart></MSHits_pepstart>
    // <MSHits_pepstop>H</MSHits_pepstop>
    // <MSHits_theomass>1101484</MSHits_theomass>
		if (tag_ == "MSHits_evalue")
		{
			actual_peptide_hit_.setScore(value.toDouble());
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_charge")
		{
			actual_peptide_hit_.setCharge(value.toInt());
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pvalue")
		{
			// TODO extra field?
			//actual_peptide_hit_.setScore(value.toDouble());
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pepstring")
		{
			actual_peptide_hit_.setSequence(value.trim());
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_mass")
		{
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pepstart")
		{
			if (value != "")
			{
				actual_peptide_hit_.setAABefore(value[0]);
			}
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pepstop")
		{
			if (value != "")
			{
				actual_peptide_hit_.setAAAfter(value[0]);
			}
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_theomass")
		{
			
			tag_ = "";
			return;
		}
	}

	} // namespace Internal
} // namespace OpenMS

