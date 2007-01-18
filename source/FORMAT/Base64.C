// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

	const char Base64::encoder_[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	const char Base64::decoder_[] = "|$$$}rstuvwxyz{$$$$$$$>?@ABCDEFGHIJKLMNOPQRSTUVW$$$$$$XYZ[\\]^_`abcdefghijklmnopq";

	Base64::Base64()
		: in_buffer_(0), 
			out_buffer_(0), 
			in_length_(0), 
			out_length_(0),
			ibuffer_size_(0)
	{
		
	}

	Base64::~Base64()
  {
		delete [] out_buffer_;
		delete [] in_buffer_;
	}

	Size Base64::getOutputBufferSize()
	{
		return out_length_;
	}

	void Base64::setOutputBufferSize(Size s)
	{
		//cout << "adpting output buffer " << out_length_ << " -> " << s << endl;
		delete [] out_buffer_;
		out_length_ = int(ceil(s/3.0))*3+1;  // add 1 for '\0'
		try
		{
			out_buffer_ = (char*) new char[out_length_];
			out_length_--;
		}
		catch(bad_alloc)
		{
			throw Exception::OutOfMemory(__FILE__, __LINE__, __PRETTY_FUNCTION__, out_length_);
		}
	}

	float* Base64::decodeFloatCorrected(const char* src, Size size)
	{
		decode(src,size);

		u_int32_t tmp;
		for (UnsignedInt i=0; i<out_length_/4; i++)
		{
			//byte order correction
			tmp = ntohl( (u_int32_t) ((u_int32_t *)out_buffer_)[i]);
			((float*)out_buffer_)[i] = *((float*)&tmp);
		}
		return (float*) out_buffer_;
	}

	float* Base64::decodeFloat(const char* src, Size size)
	{
		return (float*) decode(src,size);
	}

	float* Base64::getFloatBuffer(Size size)
	{
		in_length_ = 4*size;
		adaptInputBuffer_();
		return (float*)in_buffer_;
	}

	double* Base64::getDoubleBuffer(Size size)
	{
		in_length_ = 8*size;
		adaptInputBuffer_();
		return (double*)in_buffer_;
	}

	void Base64::adaptInputBuffer_()
	{
		// increase buffer if necessary
		if (in_length_ > ibuffer_size_)
		{
			//cout << "adpting input buffer " << ibuffer_size_ << " -> " << in_length_ << endl;
			ibuffer_size_ = in_length_+1;
			delete [] in_buffer_;
			try
			{
				in_buffer_ = (char*) new char[ibuffer_size_];
			}
			catch(bad_alloc)
			{
				throw Exception::OutOfMemory(__FILE__, __LINE__, __PRETTY_FUNCTION__, in_length_);
			}
			in_buffer_[in_length_] = '\0';
		}
	}


	char* Base64::encodeFloatCorrected()
	{
		u_int32_t tmp;

		for (UnsignedInt i=0; i<in_length_/4; i++)
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

	double* Base64::decodeDoubleCorrected(const char* src, Size size)
	{
		decode(src,size*2);

		u_int32_t tmp1, tmp2;
		for (UnsignedInt i=0; i<out_length_/8; i++)
		{
			//byte order correction
			tmp1 = ntohl( (u_int32_t) ((u_int32_t *)out_buffer_)[2*i]);
			tmp2 = ntohl( (u_int32_t) ((u_int32_t *)out_buffer_)[2*i+1]);
			((u_int32_t*)out_buffer_)[2*i] 		= tmp2;
			((u_int32_t*)out_buffer_)[2*i+1] 	= tmp1;
		}
		return (double*) out_buffer_;
	}

	double* Base64::decodeDouble(const char* src, Size size)
	{
		return (double*) decode(src,size);
	}

	char* Base64::encodeDoubleCorrected()
	{
		u_int32_t tmp1;
		u_int32_t tmp2;

		for (UnsignedInt i=0; i<in_length_/8; i++)
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

	char* Base64::encode(const char* src, Size size)
	{
		UnsignedInt padding = 0;
		if (size%3 == 2) padding=1; 
		if (size%3 == 1) padding=2;
		
		Size dest_size = (size+padding)/3*4;
		if (dest_size > out_length_) setOutputBufferSize(dest_size);

		if (size<3)
		{
			if (out_length_==0) setOutputBufferSize(1);
			out_buffer_[0] = '\0';
			return out_buffer_;
		}

		register unsigned char a;
		register unsigned char b;

		UnsignedInt j = 0;
		UnsignedInt i = 0;
		for (i=0; i<size-3; i+=3)  
		{
			// encode 3 Byte to 4 Base64-Chars
			a = src[i];
			b = src[i+1];
			out_buffer_[j++] = encoder_[a>>2];
			out_buffer_[j++] = encoder_[((a&3)<<4) | (b>>4)];
			a = src[i+2];
			out_buffer_[j++] = encoder_[((b&15)<<2) | (a>>6)];
			out_buffer_[j++] = encoder_[a&63];
		}

		// encode last 3 Byte (fill missing bits with 0)
		a = src[i++];
		out_buffer_[j++] = encoder_[a>>2];

		if (padding>1)
		{
			out_buffer_[j++] = encoder_[(a&3)<<4];
			out_buffer_[j++] = '=';
			out_buffer_[j++] = '=';
		}
		else if (padding)
		{
			b = src[i++];
			out_buffer_[j++] = encoder_[((a&3)<<4) | (b>>4)];
			out_buffer_[j++] = encoder_[(b&15)<<2];
			out_buffer_[j++] = '=';				
		} 
		else
		{
			b = src[i++];
			out_buffer_[j++] = encoder_[((a&3)<<4) | (b>>4)];
			a = src[i++];
			out_buffer_[j++] = encoder_[((b&15)<<2) | (a>>6)];
			out_buffer_[j++] = encoder_[a&63];
		}

		// terminate string 
		out_buffer_[j++] = '\0';
		
		return out_buffer_;
	}

	char* Base64::decode(const char* src, Size size)
	{
		Size src_size = size;
		// remove last one or two '=' if contained
		int padding = 0;
		if (src[src_size-1] == '=') padding++;
		if (src[src_size-2] == '=') padding++;
		src_size -= padding;

		// increase buffer if necessary			
		Size dest_size = int(ceil(src_size/4.0))*3;
		if (dest_size > out_length_) setOutputBufferSize(dest_size);

		register UnsignedInt a;
		register UnsignedInt b;

		UnsignedInt j = 0;
		for (UnsignedInt i=0; i<src_size; i+=4)
		{
			// decode 4 Base64-Chars to 3 Byte
			a = decoder_[(int)src[i]-43]-62;
			b = decoder_[(int)src[i+1]-43]-62;
			out_buffer_[j++] = (unsigned char) ((a<<2) | (b>>4));
			a = decoder_[(int)src[i+2]-43]-62;
			out_buffer_[j++] = (unsigned char) (((b&15)<<4) | (a>>2));
			b = decoder_[(int)src[i+3]-43]-62;
			out_buffer_[j++] = (unsigned char) (((a&3)<<6) | b);
		}
		// terminate string depending on padding 
		out_buffer_[dest_size-padding] = '\0';
		
		return out_buffer_;
	}

}

