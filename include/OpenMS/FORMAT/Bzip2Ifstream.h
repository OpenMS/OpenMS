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

#ifndef OPENMS_FORMAT_BZIP2IFSTREAM_H
#define	OPENMS_FORMAT_BZIP2IFSTREAM_H

#include <OpenMS/config.h>
#include <bzlib.h>
#include <istream>

namespace OpenMS
{
/**
	@brief Decompresses files which are compressed in the bzip2 format (*.bz2)
*/
	class OPENMS_DLLAPI Bzip2Ifstream 
	{
		public: 
			///Default Constructor
			Bzip2Ifstream();
			/// Detailed constructor with filename
			Bzip2Ifstream(const char * filename);
			///Destructor
			virtual ~Bzip2Ifstream();
			
			/**
					@brief reads n bytes from the bzip2 compressed file into buffer s
					
					@param s will be filled with bytes
					@param n is the size of the buffer s
					@return the number of actually read bytes. If it is 0 the end of the file was reached and the stream is closed
					
					@exception Exception::ConversionError is thrown if decompression fails
					@exception Exception::IllegalArgument is thrown if no file for decompression is given. This can happen even happen if a file was already open but read until the end.
					@note it is undefined what will happen if parameter n is bigger than the length of the buffer
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
		protected:
			/// pointer to a FILE object. Necessary for opening the file
			FILE*   file_;
			/// a pointer to a BZFILE object. Necessary for decompression
			BZFILE* bzip2file_;
			///counts the last read buffer
			size_t     n_buffer_;
			///saves the last returned error by the read function
			int     bzerror_;
			///true if end of file is reached
			bool stream_at_end_;
			
			//not implemented
			Bzip2Ifstream(const Bzip2Ifstream& bzip2);
			Bzip2Ifstream& operator=(const Bzip2Ifstream& bzip2);
	};
	
	//return bzip2file???!!!!????
	inline bool Bzip2Ifstream::isOpen() const
	{
		return (file_ != NULL);
	}
	
	inline bool Bzip2Ifstream::streamEnd() const
	{
		return stream_at_end_;
	}

} //namespace OpenMS
#endif //OPENMS_FORMAT_BZIP2IFSTREAM_H
