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

namespace OpenMS 
{

	SemanticValidator::SemanticValidator()
		: XMLHandler("", 0),
			XMLFile()
	{
	  	
	}
	
	SemanticValidator::~SemanticValidator()
	{
	}
	
  bool SemanticValidator::validate(const String& filename, const CVMappings& mapping, const ControlledVocabulary& cv)
  {
		
		return true;
  }

  void SemanticValidator::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
  {

    String tag = sm_.convert(qname);

		if (tag == "")
		{
		}
  }
	

	void SemanticValidator::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    String tag = sm_.convert(qname);

		if (tag == "")
		{
		}
  }

  void SemanticValidator::characters(const XMLCh* const /*chars*/, const unsigned int /*length*/)
  {
  }
  					 
} // namespace OpenMS
