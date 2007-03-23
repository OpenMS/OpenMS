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

#if defined OPENMS_BIG_ENDIAN
#define OPENMS_IS_BIG_ENDIAN true
#else
#define OPENMS_IS_BIG_ENDIAN false
#endif

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

	void Base64::encode(std::vector<Real>& in, ByteOrder to_byte_order, std::string& out)
	{
		bool convert = false;
		out = string();

		if ((OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::LITTLEENDIAN) ||
			(!OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BIGENDIAN))
		{
			convert = true;
		}
	
		//changeme
		UInt element_size = 4;
		UInt size = element_size * in.size();
 		UInt padding = 0;
		if (size%3 == 2) padding=1; 
		if (size%3 == 1) padding=2;
		
		register unsigned char a;
		register unsigned char b;

/*
	Inline documentation:
	3 bytes are encoded to 4 base64 chars after the following method:
	|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0|7 6 5 4 3 2 1 0| 3 Bytes
  |5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0|5 4 3 2 1 0| 4 "Bytes"

	Each resulting byte is assigned a printable character based upon its
	bit value (0..(2^6)-1), according to a special encoding character list.
	See top.

	On LITTLE endian machines, a vector<Real> looks like this:

	Memory Bytes   0 1 2 3 4 5 6 7 ... 
	Byte of Real   4 3 2 1 4 3 2 1 ...
	Vector element 1       2       ...

	On BIG endian machines, a vector<Real> looks like this:

	Memory Bytes   0 1 2 3 4 5 6 7 ... 
	Byte of Real   1 2 3 4 1 2 3 4 ...
	Vector element 1       2       ...

	When encoding to the same ENDIAN as the host byte order, we don't need to
	do any conversion.

	When encoding to the other ENDIAN, we need to
	do the conversion as follows (for every value):
  - mirror the byte orders
  - go on as usual

	To encode a vector<Real>, we need bytewise access to the original vector.
	We accomplish that by interpreting the Real as char[] and accessing it
	by index. This index is also a method for converting ENDIAN methods.
	- Starting index =
		same endian:      0, increment 1
		different endian: element_size-1, increment -1
*/

		UInt i = 0;
		UInt pos = 0;				// position in vector
		UInt offset = 0;		// offset in Real
		int inc = 1;				// increment

		if (convert == false) inc = 1;
		else inc = -1;

		for (i=0; i<size-3; i+=3)  
		{
			pos = i / element_size;

			if (convert == false)
			{
				offset = i % element_size;			// same endian
			}
			else
			{
				offset = (element_size - 1) - (i % element_size);		// other endian
			}

//			printf ("pos %d, offset %d\n", pos, offset);
//			printf ("i= %d, read %d\n", i, ((char*) &(in[pos]))[offset]);
			
			// encode 3 Byte to 4 Base64-Chars
			// a = byte at position i
			a = ((char*) &(in[pos]))[offset];

			// b = byte at position i+1
			pos = (i+1) / element_size;
			offset = (offset+inc) % element_size;
//			printf ("i+1: pos %d, offset %d\n", pos, offset);
			b = ((char*) &(in[pos]))[offset];
			out.push_back(encoder_[a>>2]);
			out.push_back(encoder_[((a&3)<<4) | (b>>4)]);

			// a = byte at position i + 2
			pos = (i+2) / element_size;
			offset = (offset+inc) % element_size;
//			printf ("i+2: pos %d, offset %d\n", pos, offset);
			a = ((char*) &(in[pos]))[offset];

			out.push_back(encoder_[((b&15)<<2) | (a>>6)]);
			out.push_back(encoder_[a&63]);
		}

		// encode last 3 Byte (fill missing bits with 0)
		pos = i / element_size;
		offset = (offset+inc) % element_size;
//		printf ("i: %d, last byte: pos %d, offset %d\n", i, pos, offset);
		a = ((char*) &(in[pos]))[offset];
		out.push_back(encoder_[a>>2]);

		if (padding == 2)
		{
/*
			One overlapping byte in input (for example 4 bytes = 1 real => 8 chars)
			last sequence:
			8 7 6 5 4 3 2 1|0 0 0 0 0 0 0 0|0 0 0 0 0 0 0 0
			6 5 4 3 2 1|6 5 4 3 2 1|6 5 4 3 2 1|6 5 4 3 2 1
			X           X           =           =
*/
			out.push_back(encoder_[(a&3)<<4]);
			out.push_back('=');
			out.push_back('=');
		}
		else if (padding)
		{
/*
			Two overlapping bytes in input (for example 8 bytes = 2 reals => 12 chars)
			last sequence:
			8 7 6 5 4 3 2 1|8 7 6 5 4 3 2 1|0 0 0 0 0 0 0 0
			6 5 4 3 2 1|6 5 4 3 2 1|6 5 4 3 2 1|6 5 4 3 2 1
			X           X           X           =
*/

			i++;
			pos = i / element_size;
			offset = (offset+inc) % element_size;
			b = ((char*) &(in[pos]))[offset];
			out.push_back(encoder_[((a&3)<<4) | (b>>4)]);
			out.push_back(encoder_[(b&15)<<2]);
			out.push_back('=');				
		} 
		else
		{
			i++;
			pos = i / element_size;
			offset = (offset+inc) % element_size;
			b = ((char*) &(in[pos]))[offset];
			out.push_back(encoder_[((a&3)<<4) | (b>>4)]);

			i++;
			pos = i / element_size;
			offset = (offset+inc) % element_size;
			a = ((char*) &(in[pos]))[offset];
			out.push_back(encoder_[((b&15)<<2) | (a>>6)]);
			out.push_back(encoder_[a&63]);
		}
}


/**
   @brief Encodes a float vector to a Base64 String
     @note @p in will be emtpy after this method
*/
void Base64::encode(std::vector<DoubleReal>& in, ByteOrder to_byte_order, std::string& out)
{
}

	/// Decodes a Base64 string to a float vector
	void Base64::decode(const std::string& in, ByteOrder from_byte_order, std::vector<Real>& out)
	{
		UInt src_size = in.size();
		// last one or two '=' are skipped if contained
		int padding = 0;
		if (in[src_size-1] == '=') padding++;
		if (in[src_size-2] == '=') padding++;

		src_size -= padding;

		register UInt a;
		register UInt b;

		UInt offset = 0;
		bool convert = false;
		int inc = 1;
		UInt written = 0;

		// changeme
		out = std::vector<Real>();
		UInt element_size = 4;

		// enough for either float or double
		char element[8] = "\x00\x00\x00\x00\x00\x00\x00";

		if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::LITTLEENDIAN) ||
			(!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BIGENDIAN))
		{
			convert = true;
			offset = (element_size - 1);		// other endian
			inc = -1;
		}
		else
		{
			offset = 0;
			inc = 1;
		}

		// sort all read bytes correctly into a char[4] (double) or
		// char[8] (Real) and push_back when necessary.
		for (UInt i=0; i<src_size; i+=4)
		{
//		printf ("start: i=%d, offset %d\n", i, offset);

			// decode 4 Base64-Chars to 3 Byte
			a = decoder_[(int)in[i]-43]-62;
			b = decoder_[(int)in[i+1]-43]-62;
			if (i+1 >= src_size) b=0;
			element[offset] = (unsigned char) ((a<<2) | (b>>4));
			written++;
//		printf ("1: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % 4 == 0)
			{
				// changeme
				out.push_back(((float*)element)[0]);
				strcpy(element, "");
			}

			a = decoder_[(int)in[i+2]-43]-62;
			if (i+2 >= src_size) a=0;
			element[offset] = (unsigned char) (((b&15)<<4) | (a>>2));
			written++;
//		printf ("2: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % 4 == 0)
			{
				// debug: output float in binary,
/*
				for (int sl = 0; sl != 4; sl++)
				{
					if (element[sl] & 128) cout << "1"; else cout << "0";
					if (element[sl] & 64) cout << "1"; else cout << "0";
					if (element[sl] & 32) cout << "1"; else cout << "0";
					if (element[sl] & 16) cout << "1"; else cout << "0";
					if (element[sl] & 8) cout << "1"; else cout << "0";
					if (element[sl] & 4) cout << "1"; else cout << "0";
					if (element[sl] & 2) cout << "1"; else cout << "0";
					if (element[sl] & 1) cout << "1"; else cout << "0";
					cout << " ";
				}
*/

				// changeme
				out.push_back(((float*)element)[0]);
				strcpy(element, "");
			}

			b = decoder_[(int)in[i+3]-43]-62;
			if (i+3 >= src_size) b=0;
			element[offset] = (unsigned char) (((a&3)<<6) | b);
			written++;
//		printf ("3: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % 4 == 0)
			{
				// changeme
				out.push_back(((float*)element)[0]);
				strcpy(element, "");
			}
		}
	}

	/// Decodes a Base64 string to a double vector
	void Base64::decode(const std::string& in, ByteOrder from_byte_order, std::vector<DoubleReal>& out)
	{
	}

	UInt Base64::getOutputBufferSize()
	{
		return out_length_;
	}

	void Base64::setOutputBufferSize(UInt s)
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

	float* Base64::decodeFloatCorrected(const char* src, UInt size)
	{
		decode(src,size);

		u_int32_t tmp;
		for (UInt i=0; i<out_length_/4; i++)
		{
			//byte order correction
			tmp = ntohl( (u_int32_t) ((u_int32_t *)out_buffer_)[i]);
			((float*)out_buffer_)[i] = *((float*)&tmp);
		}
		return (float*) out_buffer_;
	}

	float* Base64::decodeFloat(const char* src, UInt size)
	{
		return (float*) decode(src,size);
	}

	float* Base64::getFloatBuffer(UInt size)
	{
		in_length_ = 4*size;
		adaptInputBuffer_();
		return (float*)in_buffer_;
	}

	double* Base64::getDoubleBuffer(UInt size)
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

		for (UInt i=0; i<in_length_/4; i++)
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

	double* Base64::decodeDoubleCorrected(const char* src, UInt size)
	{
		decode(src,size*2);

		u_int32_t tmp1, tmp2;
		for (UInt i=0; i<out_length_/8; i++)
		{
			//byte order correction
			tmp1 = ntohl( (u_int32_t) ((u_int32_t *)out_buffer_)[2*i]);
			tmp2 = ntohl( (u_int32_t) ((u_int32_t *)out_buffer_)[2*i+1]);
			((u_int32_t*)out_buffer_)[2*i] 		= tmp2;
			((u_int32_t*)out_buffer_)[2*i+1] 	= tmp1;
		}
		return (double*) out_buffer_;
	}

	double* Base64::decodeDouble(const char* src, UInt size)
	{
//		#if defined OPENMS_BIG_ENDIAN
		return (double*) decode(src,size);
//		#else
//		return decodeDoubleCorrected(src,size);
//		#endif
	}

	char* Base64::encodeDoubleCorrected()
	{
		u_int32_t tmp1;
		u_int32_t tmp2;

		for (UInt i=0; i<in_length_/8; i++)
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

	char* Base64::encode(const char* src, UInt size)
	{
		UInt padding = 0;
		if (size%3 == 2) padding=1; 
		if (size%3 == 1) padding=2;
		
		UInt dest_size = (size+padding)/3*4;
		if (dest_size > out_length_) setOutputBufferSize(dest_size);

		if (size<3)
		{
			if (out_length_==0) setOutputBufferSize(1);
			out_buffer_[0] = '\0';
			return out_buffer_;
		}

		register unsigned char a;
		register unsigned char b;

		UInt j = 0;
		UInt i = 0;
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

	char* Base64::decode(const char* src, UInt size)
	{
		UInt src_size = size;
		// remove last one or two '=' if contained
		int padding = 0;
		if (src[src_size-1] == '=') padding++;
		if (src[src_size-2] == '=') padding++;
		src_size -= padding;

		// increase buffer if necessary			
		UInt dest_size = int(ceil(src_size/4.0))*3;
		if (dest_size > out_length_) setOutputBufferSize(dest_size);

		register UInt a;
		register UInt b;

		UInt j = 0;
		for (UInt i=0; i<src_size; i+=4)
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

