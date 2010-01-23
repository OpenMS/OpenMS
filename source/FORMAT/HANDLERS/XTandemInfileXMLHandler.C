// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/XTandemInfileXMLHandler.h>
#include <xercesc/sax2/Attributes.hpp>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
	namespace Internal
	{
  
  XTandemInfileXMLHandler::XTandemInfileXMLHandler(const String& filename, vector<XTandemInfileNote>& notes, XTandemInfile* infile) :
    XMLHandler(filename,""),
		notes_(notes),
		infile_(infile)
  {
  	
  }
   
  XTandemInfileXMLHandler::~XTandemInfileXMLHandler()
  {
    
  }
  
  void XTandemInfileXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
	{

		tag_ = String(sm_.convert(qname));
		
		if (tag_ == "note")
		{
			int type_idx = attributes.getIndex(sm_.convert("type"));
			int label_idx = attributes.getIndex(sm_.convert("label"));

			if (type_idx != -1)
			{
				actual_note_.note_type = String(sm_.convert(attributes.getValue(type_idx)));
			}
			if (label_idx != -1)
			{
				actual_note_.note_label = String(sm_.convert(attributes.getValue(label_idx)));
			}
			return;
		}

	}
	  
  void XTandemInfileXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
 	{
 		tag_ = String(sm_.convert(qname)).trim();
 		if (tag_ == "note")
		{
			return;
		}
 	} 

  void XTandemInfileXMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
  {
		String value = ((String)sm_.convert(chars)).trim();
		if (tag_ == "note")
		{
      actual_note_.note_value = value;
			notes_.push_back(actual_note_);
			actual_note_ = XTandemInfileNote();
      return;
		}
	}

	} // namespace Internal
} // namespace OpenMS
