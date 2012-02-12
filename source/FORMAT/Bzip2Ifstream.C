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
#include <iostream>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <cstdlib>
using namespace std;
namespace OpenMS
{
	Bzip2Ifstream::Bzip2Ifstream(const char * filename) : n_buffer_(0),stream_at_end_(false)
	{
		file_ = fopen( filename, "rb" ); //read binary: always open in binary mode because windows and mac open in text mode
		
		//aborting, ahhh!
		if( !file_ ) 
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
		bzip2file_ = BZ2_bzReadOpen ( &bzerror_, file_, 0, 0, NULL, 0 );
		if ( bzerror_ != BZ_OK ) 
		{
	  	close();
	  	throw Exception::ConversionError(__FILE__,__LINE__,__PRETTY_FUNCTION__,"bzip2 compression failed: ");
		}
	}
	
	Bzip2Ifstream::Bzip2Ifstream()
		: file_(NULL),bzip2file_(NULL),n_buffer_(0),bzerror_(0),stream_at_end_(true)
	{
	}
	
	Bzip2Ifstream::~Bzip2Ifstream()
	{
		close();
	}
	
	size_t Bzip2Ifstream::read(char* s, size_t n)
	{
		if(bzip2file_ != NULL)
		{
			bzerror_ = BZ_OK;
  		n_buffer_ = BZ2_bzRead ( &bzerror_, bzip2file_, s, (unsigned int)n/* size of buf */ );		
	  	if(bzerror_ == BZ_OK) 
	 		{
    			return n_buffer_;
  		}
			else if(bzerror_ != BZ_STREAM_END) 
			{
   			close();
   			throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,"	","bzip2 compression failed: ");
			} 
			else 
			{
   			close();
   			return n_buffer_;
			}
		}
		else
		{
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"no file for decompression initialized");
		}
	}
	
	void Bzip2Ifstream::open(const char* filename)
	{
		close();
		file_ = fopen( filename, "rb" ); //read binary: always open in binary mode because windows and mac open in text mode
		
		//aborting, ahhh!
		if( !file_ ) 
		{
			throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		
		bzip2file_ = BZ2_bzReadOpen ( &bzerror_, file_, 0, 0, NULL, 0 );
		if ( bzerror_ != BZ_OK ) 
		{
	  	close();
	  	throw Exception::ConversionError(__FILE__,__LINE__,__PRETTY_FUNCTION__,"bzip2 compression failed: ");
		}
		stream_at_end_ = false;
	}
	
	void Bzip2Ifstream::close()
	{
		if(bzip2file_ != NULL)
		{
			BZ2_bzReadClose(&bzerror_,bzip2file_);
		}
		if(file_ != NULL)
		{
			fclose(file_);
		}
		file_ = NULL;
		bzip2file_ = NULL;
		stream_at_end_ = true;
	}	

} //namespace OpenMS
