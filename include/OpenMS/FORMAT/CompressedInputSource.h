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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_COMPRESSEDINPUTSOURCE_H
#define OPENMS_FORMAT_COMPRESSEDINPUTSOURCE_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <xercesc/sax/InputSource.hpp>

namespace OpenMS
{
	/**
		@brief This class is based on xercesc::LocalFileInputSource
	*/
	class OPENMS_DLLAPI CompressedInputSource 
		:	public xercesc::InputSource
	{
		public:
			///Constructor
			CompressedInputSource(const   String& file_path,const char* header , xercesc::MemoryManager* const manager = xercesc::XMLPlatformUtils::fgMemoryManager);
			///Constructor
			CompressedInputSource(const   XMLCh* const file_path,const char* header, xercesc::MemoryManager* const manager = xercesc::XMLPlatformUtils::fgMemoryManager);
		  ///Constructor
		  virtual ~CompressedInputSource();
		 
		 /**
		  	@brief Depending on the header in the Constructor a Bzip2InputStream or a GzipInputStream object is returned
			  @note InputSource interface implementation
			*/
		  virtual xercesc::BinInputStream* makeStream() const;
		  
		private:
			 char head_[2];
		  //not implemented
		  CompressedInputSource();
		  CompressedInputSource(const CompressedInputSource& source);
		  CompressedInputSource& operator=(const CompressedInputSource& source);
	};
	
} // namespace OpenMS

#endif // OPENMS_FORMAT_COMPRESSEDINPUTSOURCE_H
