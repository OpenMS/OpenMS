// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzXMLHandler.h>

namespace OpenMS
{
				/*
	namespace Internal
	{
		void MzXMLHandler::initStaticMembers_()
  	{
  		static bool init(false);
    	if (!init)
    	{
    		s_value_ = xercesc::XMLString::transcode("value");
      s_count_ = xercesc::XMLString::transcode("scanCount");
      s_type_ = xercesc::XMLString::transcode("type");
      s_name_ = xercesc::XMLString::transcode("name");
      s_version_ = xercesc::XMLString::transcode("version");
      s_filename_ = xercesc::XMLString::transcode("fileName");
      s_filetype_ = xercesc::XMLString::transcode("fileType");
      s_filesha1_ = xercesc::XMLString::transcode("fileSha1");
      s_completiontime_ = xercesc::XMLString::transcode("completionTime");
    	s_precision_ = xercesc::XMLString::transcode("precision");
    	s_byteorder_ = xercesc::XMLString::transcode("byteOrder");
    	s_pairorder_ = xercesc::XMLString::transcode("pairOrder");
    	s_precursorintensity_ = xercesc::XMLString::transcode("precursorIntensity");
    	s_precursorcharge_ = xercesc::XMLString::transcode("precursorCharge");
    	s_windowwideness_ = xercesc::XMLString::transcode("windowWideness");
    	s_mslevel_ = xercesc::XMLString::transcode("msLevel");
    	s_peakscount_ = xercesc::XMLString::transcode("peaksCount");
    	s_polarity_ = xercesc::XMLString::transcode("polarity");
    	s_scantype_ = xercesc::XMLString::transcode("scanType");
    	s_retentiontime_ = xercesc::XMLString::transcode("retentionTime");
    	s_collisionenergy_ = xercesc::XMLString::transcode("collisionEnergy");
    	s_startmz_ = xercesc::XMLString::transcode("startMz");
    	s_endmz_ = xercesc::XMLString::transcode("endMz");
    	s_first_ = xercesc::XMLString::transcode("first");
    	s_last_ = xercesc::XMLString::transcode("last");
    	s_phone_ = xercesc::XMLString::transcode("phone");
    	s_email_ = xercesc::XMLString::transcode("email");
    	s_uri_ = xercesc::XMLString::transcode("URI");
    	s_num_ = xercesc::XMLString::transcode("num");
    	s_intensitycutoff_ = xercesc::XMLString::transcode("intensityCutoff");
    	s_centroided_ = xercesc::XMLString::transcode("centroided");
    	s_deisotoped_ = xercesc::XMLString::transcode("deisotoped");
    	s_chargedeconvoluted_ = xercesc::XMLString::transcode("chargeDeconvoluted");

    	init = true;
    	}
  		return;
		}
	}
	*/
}
