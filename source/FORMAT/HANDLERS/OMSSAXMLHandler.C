// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
								  									 vector<PeptideIdentification>& id_data, 
      								 							 const String& filename) :
    XMLHandler(filename),
    protein_identification_(protein_identification),
    id_data_(id_data),
    actual_protein_hit_(),
    actual_peptide_hit_(),/*
    peptide_identification_index_(0),*/
		tag_()/*,
		date_() */       
  {
  	
  }
   
  OMSSAXMLHandler::~OMSSAXMLHandler()
  {
    
  }
  
  void OMSSAXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& /*attributes*/)
	{
		tag_ = String(XMLString::transcode(qname));
	
		if (tag_ == "MSPepHit")
		{
			//tag_ = "";
			return;
		}
		if (tag_ == "MSHits")
		{
			//actual_peptide_hit_.clear();
			actual_peptide_hit_ = PeptideHit();
			//tag_ = "";
			return;
		}

		if (tag_ == "MSHitSet")
		{
			//tag_ = "";
			return;
		}
	}
	  
  void OMSSAXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
 	{
 		tag_ = String(XMLString::transcode(qname)).trim();
 		if (tag_ == "MSHits")
		{
			///cerr << "MSHits " << actual_peptide_hit_.getSequence() << endl;
			
			//actual_peptide_hit_.setScoreType("OMSSA");
			//IdentificationData id_dat;
			//id_dat.id = ProteinIdentification();
			//id_dat.id.getHits().push_back(actual_peptide_hit_);
			//id_data_.push_back(id_dat);
			if (id_data_.size() == 0)
			{
				//IdentificationData id_dat;
				//id_dat.id = ProteinIdentification();
				//id_dat.id.getHits().push_back(actual_peptide_hit_);
				//id_data_.push_back(id_dat);
			}
			else
			{
				//id_data_[0].id.getHits().push_back(actual_peptide_hit_);
			}
			tag_ = "";
			return;
		}
		tag_ = "";
 	} 

  void OMSSAXMLHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
  {
		//cerr << tag_ << " '" << XMLString::transcode(chars) << "'" << endl;
		// MSPepHit section
    // <MSPepHit_start>0</MSPepHit_start>
    // <MSPepHit_stop>8</MSPepHit_stop>
    // <MSPepHit_accession>6599</MSPepHit_accession>
    // <MSPepHit_defline>CRHU2 carbonate dehydratase (EC 4.2.1.1) II [validated] - human</MSPepHit_defline>
    // <MSPepHit_protlength>260</MSPepHit_protlength>
    // <MSPepHit_oid>6599</MSPepHit_oid>
		if (tag_ == "MSPepHit_start")
		{
      // TODO
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_stop")
		{
      // TODO
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_accession")
		{
      // TODO
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_defline")
		{
      // TODO
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_protlength")
		{
      // TODO
      tag_ = "";
      return;
		}
		if (tag_ == "MSPepHit_oid")
		{
      // TODO
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
			// TODO
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_charge")
		{
			//cerr << tag_ << " '" << XMLString::transcode(chars) << "'" << endl;
			actual_peptide_hit_.setCharge(((String) XMLString::transcode(chars)).trim().toInt());
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pvalue")
		{
			//cerr << tag_ << " '" << XMLString::transcode(chars) << "'" << endl;
			actual_peptide_hit_.setScore(((String) XMLString::transcode(chars)).trim().toDouble());
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pepstring")
		{
			//cerr << tag_ << " '" << XMLString::transcode(chars) << "'" << endl;
			actual_peptide_hit_.setSequence(((String)XMLString::transcode(chars)).trim());
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
			// TODO
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pepstop")
		{
			// TODO
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_theomass")
		{
			// TODO
			tag_ = "";
			return;
		}
	}

	} // namespace Internal
} // namespace OpenMS
