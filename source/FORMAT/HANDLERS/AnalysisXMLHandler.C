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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/AnalysisXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{
	namespace Internal
	{

  AnalysisXMLHandler::AnalysisXMLHandler(const Identification& id, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger)
  {
  	cv_.loadFromOBO("PI",File::find("/CV/psi-pi.obo"));
  }

  AnalysisXMLHandler::AnalysisXMLHandler(Identification& id, const String& filename, const String& version, const ProgressLogger& logger)
		: XMLHandler(filename, version),
    	logger_(logger)
  {
  	cv_.loadFromOBO("PI",File::find("/CV/psi-pi.obo"));
  }	

	AnalysisXMLHandler::~AnalysisXMLHandler()
	{
	}

	void AnalysisXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		String tag = sm_.convert(qname);
		if (tag == "AnalysisXML")
		{
			// TODO handle version
			return;
		}

		if (tag == "SpectrumIdentificationList")
		{
			// This corresponds to a new "PeptideIdentification List"
			
			return;
		}

		if (tag == "SpectrumIdentificationResult")
		{
			// This corresponds to a new "PeptideIdentification"

			return;
		}

		if (tag == "SpectrumIdentificationItem")
		{
			// This corresponds to a new "PeptideHit"

			return;
		}

	}

	void AnalysisXMLHandler::characters(const XMLCh* const chars, const XMLSize_t length)
	{
	}

	void AnalysisXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
	}
	
	void AnalysisXMLHandler::writeTo(std::ostream& os)
	{
	}


	} //namespace Internal
} // namespace OpenMS


