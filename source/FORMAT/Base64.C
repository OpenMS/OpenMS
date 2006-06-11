// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: Base64.C,v 1.9 2006/03/28 08:03:32 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Base64.h>
#include <iostream>

namespace OpenMS
{

	const char Base64::encoder_[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	const char Base64::decoder_[] = "|$$$}rstuvwxyz{$$$$$$$>?@ABCDEFGHIJKLMNOPQRSTUVW$$$$$$XYZ[\\]^_`abcdefghijklmnopq";

	Base64::Base64(): in_buffer_(0), out_buffer_(0), in_length_(0), out_length_(0)
	{}

	Base64::~Base64()
  {
		if (out_length_!=0)	delete [] out_buffer_;
	}

	Size Base64::getOutputBufferSize()
	{
		return out_length_;
	}

	void Base64::setOutputBufferSize(Size s)
	{
		if (out_length_!=0) delete [] out_buffer_;
		out_length_ = int(ceil(s/3.0))*3+1;  // add 1 for '\0'
		try{
				out_buffer_ = (char*) new char[out_length_];
				out_length_--;
		}catch( std::bad_alloc)
		{
			throw Exception::OutOfMemory(__FILE__, __LINE__, __PRETTY_FUNCTION__, out_length_);
		}
	}

	//--------- 32 bit precision -------------------------------
	float* Base64::decodeFloatCorrected(const char* src, const Size size)
	{
		decode(src,size);
		Size n = out_length_/4;

		u_int32_t tmp;
		for (UnsignedInt i=0; i<n; i++)
		{
			//byte order correction
			tmp = ntohl( (u_int32_t) ((u_int32_t *)out_buffer_)[i]);
			((float*)out_buffer_)[i] = *((float*)&tmp);
		}
		return (float*) out_buffer_;
	}

	float* Base64::decodeFloat(const char* src, const Size size)
	{
		return (float*) decode(src,size);
	}

	float* Base64::getFloatBuffer(const Size size)
	{
		// increase buffer if necessary
		if (in_length_!=0 && 4*size > in_length_) delete [] in_buffer_;
		in_length_ = 4*size+1;
		try{
				in_buffer_ = (char*) new char[in_length_];
				in_buffer_[in_length_-1] = '\0';
				--in_length_;
		}catch( std::bad_alloc)
		{
			throw Exception::OutOfMemory(__FILE__, __LINE__, __PRETTY_FUNCTION__, in_length_);
		}
		return (float*)in_buffer_;
	}

	char* Base64::encodeFloatCorrected()
	{
		Size n = in_length_/4;
		u_int32_t tmp;

		for (UnsignedInt i=0; i<n; i++)
		{
			//byte order correction
			tmp = htonl( (u_int32_t) ((u_int32_t *)in_buffer_)[i]);
			((float*)in_buffer_)[i] = *((float*)&tmp);
		}
		return encode(in_buffer_,in_length_);

	}

	char* Base64::encodeFloat()
	{
		return encode(in_buffer_,in_length_);
	}

	//--------- 64 bit precision -------------------------------
	double* Base64::decodeDoubleCorrected(const char* src, const Size size)
	{
		decode(src,size*2);
		Size n = out_length_/8;

		u_int32_t tmp1, tmp2;
		for (UnsignedInt i=0; i<n; i++)
		{
			//byte order correction
			tmp1 = ntohl( (u_int32_t) ((u_int32_t *)out_buffer_)[2*i]);
			tmp2 = ntohl( (u_int32_t) ((u_int32_t *)out_buffer_)[2*i+1]);
			((u_int32_t*)out_buffer_)[2*i] 		= tmp2;
			((u_int32_t*)out_buffer_)[2*i+1] 	= tmp1;
		}
		return (double*) out_buffer_;
	}

	double* Base64::decodeDouble(const char* src, const Size size)
	{
		return (double*) decode(src,size);
	}

	double* Base64::getDoubleBuffer(const Size size)
	{
		// increase buffer if necessary
		if (in_length_!=0 && 8*size > in_length_) delete [] in_buffer_;
		in_length_ = 8*size+1;
		try{
				in_buffer_ = (char*) new char[in_length_];
				in_buffer_[in_length_-1] = '\0';
				--in_length_;				
		}catch( std::bad_alloc)
		{
			throw Exception::OutOfMemory(__FILE__, __LINE__, __PRETTY_FUNCTION__, in_length_);
		}
		return (double*)in_buffer_;
	}

	char* Base64::encodeDoubleCorrected()
	{
		Size n = in_length_/8;
		u_int32_t tmp1;
		u_int32_t tmp2;

		for (UnsignedInt i=0; i<n; i++)
		{
			//byte order correction
			tmp1 = htonl( (u_int32_t) ((u_int32_t *)in_buffer_)[2*i]);
			tmp2 = htonl( (u_int32_t) ((u_int32_t *)in_buffer_)[2*i+1]);
			((u_int32_t *)in_buffer_)[2*i] 		= tmp2;
			((u_int32_t *)in_buffer_)[2*i+1] 	= tmp1;
		}
		return encode(in_buffer_,in_length_);
	}

	char* Base64::encodeDouble()
	{
		return encode(in_buffer_,in_length_);
	}

}

