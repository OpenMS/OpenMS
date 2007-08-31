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

#include <OpenMS/FORMAT/HANDLERS/XTandemInfileXMLHandler.h>
#include <xercesc/sax2/Attributes.hpp>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
	namespace Internal
	{
  
  XTandemInfileXMLHandler::XTandemInfileXMLHandler(const String& filename, XTandemInfile* infile) :
    XMLHandler(filename),
		infile_(infile),
		tag_()/*,
		date_() */       
  {
  	
  }
   
  XTandemInfileXMLHandler::~XTandemInfileXMLHandler()
  {
    
  }
  
  void XTandemInfileXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
	{

		tag_ = String(sm_.convert(qname));
		
		if (tag_ == "domain")
		{
			return;
		}

	}
	  
  void XTandemInfileXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
 	{
 		tag_ = String(sm_.convert(qname)).trim();
 		if (tag_ == "domain")
		{
			tag_ = "";
			return;
		}
		tag_ = "";
 	} 

  void XTandemInfileXMLHandler::characters(const XMLCh* const chars, const unsigned int /*length*/)
  {
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
		//cerr << "tag_: " << tag_ << ", text: " <<  sm_.convert(chars) << endl;
		if (tag_ == "MSHits_evalue")
		{
			// TODO
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_charge")
		{
			//cerr << tag_ << " '" << sm_.convert(chars) << "'" << endl;
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pvalue")
		{
			//cerr << tag_ << " '" << sm_.convert(chars) << "'" << endl;
			tag_ = "";
			return;
		}
		if (tag_ == "MSHits_pepstring")
		{
			//cerr << tag_ << " '" << sm_.convert(chars) << "'" << endl;
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
