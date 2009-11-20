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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_GZIPINPUTSTREAM_H
#define OPENMS_FORMAT_GZIPINPUTSTREAM_H

#include <xercesc/util/BinInputStream.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <OpenMS/FORMAT/GzipIfstream.h>


namespace OpenMS
{
	class String;
	/**
		@brief Implements the BinInputStream class of the xerces-c library in order to read gzip compressed XML files.
		
	*/
	class OPENMS_DLLAPI GzipInputStream
		:	public xercesc::BinInputStream
	{
		public:
			///Constructor
			GzipInputStream(const   String& file_name);

   		GzipInputStream(const   char* const     file_name);	 
   		
   		
   		///Destructor
   		virtual ~GzipInputStream();
   		///returns true if file is open
   		 bool getIsOpen() const;
    	/**
    	@brief returns the current position in the file
    	@note Implementation of the xerces-c input stream interface
    	*/
    	virtual XMLFilePos curPos() const;
			
			/**
				@brief writes bytes into buffer from file
				@note Implementation of the xerces-c input stream interface
				
				@param to_fill is the buffer which is written to
				@param max_to_read is the size of the buffer
				
				@return returns the number of bytes which were actually read
			
			*/
   	 	virtual XMLSize_t readBytes(XMLByte* const  to_fill, const XMLSize_t max_to_read);
			/**
				@brief returns 0
    		@note Implementation of the xerces-c input stream interface				
			*/
	    virtual const XMLCh* getContentType() const;

   		
    private:
    ///pointer to an compression stream
    	GzipIfstream* 	gzip_;
    	///current index of the actual file
    	XMLSize_t       file_current_index_;
    	
    	//not implemented
    	GzipInputStream();
    	GzipInputStream(const GzipInputStream& stream);
    	GzipInputStream& operator=(const GzipInputStream& stream);
	};
	
	inline XMLFilePos GzipInputStream::curPos() const
	{
    return file_current_index_;
	}
	
	inline bool GzipInputStream::getIsOpen() const
	{
			return gzip_->isOpen();
	}
} // namespace OpenMS

#endif // OPENMS_FORMAT_GZIPInputStream_H
