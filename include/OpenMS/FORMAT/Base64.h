// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FORMAT_BASE64_H
#define OPENMS_FORMAT_BASE64_H

#ifndef OPENMS_IS_BIG_ENDIAN
#if defined OPENMS_BIG_ENDIAN
#define OPENMS_IS_BIG_ENDIAN true
#else
#define OPENMS_IS_BIG_ENDIAN false
#endif
#endif

#include <OpenMS/CONCEPT/Types.h>

#include <cmath>
#include <string>
#include <vector>
#include <cstring>

namespace OpenMS
{
  /**
  	@brief Class to encode and decode Base64
		
		Base64 supports two precisions: 32 Bit (float) and 64 Bit (double).
		
  */
  class OPENMS_DLLAPI Base64
  {
  public:
    /// default constructor
    Base64();

    /// Destructor
    virtual ~Base64();

		enum ByteOrder{BYTEORDER_BIGENDIAN,BYTEORDER_LITTLEENDIAN};

		/**
		   @brief Encodes a vector of floating point numbers to a Base64 String
		     @note @p in will be emtpy after this method
		*/

		template <typename FromType>
		void encode(std::vector<FromType>& in, ByteOrder to_byte_order, std::string& out);

		/**
		   @brief Decodes a Base64 string to a vector of floating point numbers
		     @note @p in will be emtpy after this method
		*/

		template <typename ToType>
		void decode(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out);

	private:
		static const char encoder_[];
		static const char decoder_[];
  };


	template <typename FromType>
	void Base64::encode(std::vector<FromType>& in, ByteOrder to_byte_order, std::string& out)
	{
		bool convert = false;
		out.clear();
		if (in.size() == 0) return;
		
		if ((OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BYTEORDER_LITTLEENDIAN) ||
			(!OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BYTEORDER_BIGENDIAN))
		{
			convert = true;
		}
	
		UInt element_size = sizeof (FromType);
		UInt size = element_size * in.size();
 		UInt padding = 0;
		if (size%3 == 2) padding=1; 
		if (size%3 == 1) padding=2;
		
		register unsigned char a;
		register unsigned char b;

		//reserve enough space in the output string
		out.reserve((UInt)(std::ceil((3.0*size)/4.0)+8.0));

/*
	Inline documentation:

  If you want to understand the following, this link might be crucial:

  http://babbage.cs.qc.edu/IEEE-754/32bit.html (and some links at the bottom of this document)
  (online converter - *extremely* helpful for debugging)

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

	/// Decodes a Base64 string to a float vector
	template <typename ToType>
	void Base64::decode(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out)
	{
		out.clear();
		if (in == "") return;
	
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

		UInt element_size = sizeof(ToType);

		// enough for either float or double
		char element[8] = "\x00\x00\x00\x00\x00\x00\x00";

		if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) ||
			(!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
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

		//reserve enough space in the output vector
		out.reserve((UInt)(std::ceil((4.0*src_size)/3.0)+6.0));

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

			if (written % sizeof(ToType) == 0)
			{
				out.push_back(((ToType*)element)[0]);
				strcpy(element, "");
			}

			a = decoder_[(int)in[i+2]-43]-62;
			if (i+2 >= src_size) a=0;
			element[offset] = (unsigned char) (((b&15)<<4) | (a>>2));
			written++;
//		printf ("2: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % sizeof(ToType) == 0)
			{
				// debug: output float in binary
/*
				for (int sl = 0; sl != sizeof(ToType); sl++)
				{
					for (int sl2 = 128; sl2 >= 1; sl2 /= 2)
					{
						std::cout << (element[sl] & sl2);
					}
					std::cout << " ";
				}
*/
				out.push_back(((ToType*)element)[0]);
				strcpy(element, "");
			}

			b = decoder_[(int)in[i+3]-43]-62;
			if (i+3 >= src_size) b=0;
			element[offset] = (unsigned char) (((a&3)<<6) | b);
			written++;
//		printf ("3: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % sizeof(ToType) == 0)
			{
				out.push_back(((ToType*)element)[0]);
				strcpy(element, "");
			}
		}
	}

} //namespace OpenMS

#endif /* OPENMS_FORMAT_BASE64_H */
