// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FORMAT_GZIPIFSTREAM_H
#define	OPENMS_FORMAT_GZIPIFSTREAM_H

#include <OpenMS/config.h>
#include <zlib.h>

namespace OpenMS
{
/**
	@brief Decompresses files which are compressed in the gzip format (*.gzip)
*/
	class OPENMS_DLLAPI GzipIfstream 
	{
		public: 
			///Default Constructor
			GzipIfstream();
			/// Detailed constructor with filename
			GzipIfstream(const char * filename);
			///Destructor
			virtual ~GzipIfstream();
			
			/**
					@brief reads n bytes from the gzip compressed file into buffer s
					
					@param s will be filled with bytes
					@param n is the size of the buffer s
					@return the number of actually read bytes. If it is 0 the end of the file was reached and the stream is closed
					
					@exception Exception::ConversionError is thrown if decompression fails
					@exception Exception::IllegalArgument is thrwon if no file for decompression is given. This can happen even happen if a file was already open but read until the end.
			*/
			size_t read(char* s, size_t n);
			
			/**
				@brief indicates whether the read function can be used safely
				
				@return true if end of file was reached. Otherwise false.
			*/
			bool streamEnd() const;
			
			/**
				@brief returns whether a file is open. 
			*/
			bool isOpen() const;
			
			/**
				@brief opens a file for reading (decompression)
				@note any previous open files will be closed first!
			*/
			void open(const char* filename);
			
			/**
				@brief closes current file.
			*/
			void close();
	
			/*
				@brief updates crc32 check sum whether the buffer is corrupted
				@note if this function is used it has to be called after every call of function read
				@param s the buffer which will be checked
				@param n the size of the buffer
			*	
			//void updateCRC32(const char* s,const size_t n);
			
			*
				@brief	checks if data is corrupted after crc32 was computed
				@note 	it can only be used if updateCRC32 was called after every call of function read
				@return true if the buffer and hence the file is corrupted; no decompression is possible
			*
			//bool isCorrupted();
			
			//unsigned long Crc32_ComputeBuf( unsigned long inCrc32, const void *buf,
      //                                 size_t bufLen );*/
			
		protected:

			///a gzFile object(void*) . Necessary for decompression
			gzFile gzfile_;
			///counts the last read bufffer
			int     n_buffer_;
			///saves the last returned error by the read function
			int     gzerror_;
			///true if end of file is reached
			bool stream_at_end_;
			
			//needed if one wants to know whetther file is okay
			//unsigned long original_crc;
			//needed if one wants to know whetther file is okay			
    	//		unsigned long crc;
			
			///not implemented
			GzipIfstream(const GzipIfstream& bzip2);
			GzipIfstream& operator=(const GzipIfstream& bzip2);
	};
	
	inline bool GzipIfstream::isOpen() const
	{
		return (gzfile_ != NULL);
	}
	
	inline bool GzipIfstream::streamEnd() const
	{
		return stream_at_end_;
	}
/*	inline bool GzipIfstream::isCorrupted()
	{
		std::cout<<"CRC"<<crc<<std::endl;
		return (crc != original_crc);
	}*/	

} //namespace OpenMS
#endif //OPENMS_FORMAT_GZIPIFSTREAM_H
