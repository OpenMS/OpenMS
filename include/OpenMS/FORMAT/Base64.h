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
// $Maintainer: Marc Sturm $
// $Authors: $
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
#include <OpenMS/DATASTRUCTURES/String.h>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <vector>
#include <iostream>


#include <QtCore/QString>
namespace OpenMS
{
  /**
  	@brief Class to encode and decode Base64
		
		Base64 supports two precisions: 32 bit (float) and 64 bit (double).
  */
  class OPENMS_DLLAPI Base64
  {
  	
	  public:
	  
	    /// default constructor
	    Base64();
	
	    /// Destructor
	    virtual ~Base64();
			
			/// Byte order type
			enum ByteOrder
			{
				BYTEORDER_BIGENDIAN,		///< Big endian type
				BYTEORDER_LITTLEENDIAN	///< Little endian type
			};
			
			/**
				@brief Encodes a vector of floating point numbers to a Base64 string
				
				You can specify the byte order of the output and if it is zlib-compressed.
				
				@note @p in will be emtpy after this method
			*/
			template <typename FromType>
			void encode(std::vector<FromType>& in, ByteOrder to_byte_order, std::string& out, bool zlib_compression= false);
	
			/**
				@brief Decodes a Base64 string to a vector of floating point numbers
	
				You have to specify the byte order of the input and if it is zlib-compressed.
	
				@note @p in will be emtpy after this method
			*/
			template <typename ToType>
			void decode(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out, bool zlib_compression = false);
			
		private:
			
			static const char encoder_[];
			static const char old_decoder_[];
			char decoder_[256];
			bool table_initialized_;
			
			/**
				@brief Decodes a Base64 string to a vector of floating point numbers
				
				@note this function is not able to decode compresed data
			*/
			template <typename ToType>
			void decode_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out);
		 
			///Decodes a compressed Base64 string to a vector of floating point numbers 
			template <typename ToType>
			void decodeCompressed_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out,bool zlib_compression );
			
			///Endianizes a 32 bit type from big endian to litte endian and vice versa
			inline Int32 endianize32_(Int32 n);
	
			///Endianizes a 64 bit type from  big endian to litte endian and vice versa
			inline Int64 endianize64_(Int64 n);
  };

	inline Int32 Base64::endianize32_(Int32 n)
	{
		return ((n&0xff)<<24) | ((n&0xff00)<<8) | ((n&0xff0000)>>8) | ((n&0xff000000)>>24);
	}
	
	inline Int64 Base64::endianize64_(Int64 n)
	{
		return ((n&0x00000000000000ffll)<<56) | 
		((n&0x000000000000ff00ll)<<40) | 
		((n&0x0000000000ff0000ll)<<24) | 
		((n&0x00000000ff000000ll)<<8)  |
		((n&0x000000ff00000000ll)>>8)  | 
		((n&0x0000ff0000000000ll)>>24) |
		((n&0x00ff000000000000ll)>>40) | 
		((n&0xff00000000000000ll)>>56);
	}
	
	template <typename FromType>
	void Base64::encode(std::vector<FromType>& in, ByteOrder to_byte_order, std::string& out, bool zlib_compression )
	{
		out.clear();
		if (in.size() == 0) return;

		//initialize
		const Size element_size = sizeof(FromType);
		const Size input_bytes = element_size * in.size();
		
		//Change endianness if necessary
		if ((OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BYTEORDER_BIGENDIAN))
		{
			if (element_size == 4)
			{
				for (Size i = 0; i < in.size(); ++i)
				{
					Int32 tmp = endianize32_(reinterpret_cast<Int32&>(in[i]));
					in[i] = reinterpret_cast<FromType&>(tmp);
				}
			}
			else if (element_size == 8)
			{
				for (Size i = 0; i < in.size(); ++i)
				{
					Int64 tmp = endianize64_(reinterpret_cast<Int64&>(in[i]));
					in[i] = reinterpret_cast<FromType&>(tmp);
				}
			}
		}
		
		//encode with compression (use Qt because of zlib support)
		if (zlib_compression)
		{
			QByteArray original = QByteArray::fromRawData(reinterpret_cast<const char*>(&in[0]), (int) input_bytes);
			QByteArray compressed = qCompress((uchar*)original.data(),original.size());
			QByteArray extern_compressed = compressed.right(compressed.size() - 4);			
			QByteArray base64_compressed = extern_compressed.toBase64();
	
			out = QString(base64_compressed).toStdString();
		}
		//encode without compression
		else
		{
			out.resize((Size)ceil(input_bytes/3.) * 4); //resize output array in order to have enough space for all characters
			Byte* to = reinterpret_cast<Byte*>(&out[0]);
			Byte* it = reinterpret_cast<Byte*>(&in[0]);
			Byte* end = it + input_bytes;
			Size written = 0;

			while (it!=end)
			{
				Int int_24bit = 0;
				Int padding_count = 0;

				// construct 24-bit integer from 3 bytes
				for (Size i=0; i<3; i++)
				{
					if (it!=end)
					{
						int_24bit |= *it++<<((2-i)*8);
					}
					else
					{
						padding_count++;
					}
				}

				// write out 4 characters
				for (Int i=3; i>=0; i--)
				{
					to[i] = encoder_[int_24bit & 0x3F];
					int_24bit >>= 6;
				}

				// fixup for padding
				if (padding_count > 0) to[3] = '=';
				if (padding_count > 1) to[2] = '=';

				to += 4;
				written += 4;
			}

			out.resize(written); //no more space is needed
		}
	}
	
	template<typename ToType>
	void Base64::decode(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out, bool zlib_compression)
	{
		if(zlib_compression)
		{
			decodeCompressed_(in,from_byte_order,out,zlib_compression);
		}
		else
		{
			decode_(in,from_byte_order,out);
		}
	}
	
	template <typename ToType>
	void Base64::decodeCompressed_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out, bool zlib_compression)
	{
		out.clear();
		if (in == "") return;
		
		Size src_size = in.size();
		
		void* byteBuffer;
		Size bufferSize;		
		std::vector<unsigned char> binary;
		const Size element_size = sizeof(ToType);
		
		String decompressed;
		if (zlib_compression)
		{
			QByteArray herewego = QByteArray::fromRawData(in.c_str(), (int) in.size());
			QByteArray bazip = QByteArray::fromBase64(herewego);
			QByteArray czip;
			czip.resize(4);
			czip[0] = (bazip.size() & 0xff000000) >> 24;
			czip[1] = (bazip.size() & 0x00ff0000) >> 16;
			czip[2] = (bazip.size() & 0x0000ff00) >> 8;
			czip[3] = (bazip.size()& 0x000000ff);
			czip += bazip;
			QByteArray base64_uncompressed = qUncompress(czip);
			
			if(base64_uncompressed.isEmpty())
			{
				throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Decompression error?");
			}
			decompressed.resize(base64_uncompressed.size());
			
			//decompressed  =  QString(base64_uncompressed).toStdString();  // Unfortunately this method doesnt write all characters to the string always
			for(Int i = 0 ; i < base64_uncompressed.size(); ++i)
			{
				decompressed[i] = base64_uncompressed[i];
			}
			byteBuffer = reinterpret_cast<void*>(&decompressed[0]);
			bufferSize = decompressed.size();
		}
		else
		{
		  binary.resize((Size)ceil(src_size/4.) * 3);

			if(!table_initialized_)
			{
				for (UInt i=0; i<64; i++)
				{
					decoder_[static_cast<int>(encoder_[i])] = static_cast<char>(i);
				}
				table_initialized_ = true;
			}

			Byte* it = (Byte*)&in[0];
			Byte* end = it + src_size;
			Byte* result = (Byte*)&binary[0];
			Size written = 0;

			while (it!=end)
			{
				Int int_24bit = 0;
				Int padding_count = 0;

				// construct 24-bit integer from 4 characters
				for (Int i=0; i<4 && it!=end; i++, it++)
				{
					if (*it != '=')
					{
						int_24bit |= decoder_[*it]<<((3-i)*6);
					}
					else
					{
						padding_count++;
					}
				}

				// write out bytes
				for (Int i=0; i<3-padding_count; i++)
				{
					Byte temp = static_cast<Byte>(int_24bit>>((2-i)*8));
					*result++ = temp;
					int_24bit ^= temp<<((2-i)*8);
					written++;
				}

			}
				
			while(!zlib_compression && written%element_size!= 0)
			{
				Byte temp = 0;

				*result++ = temp;
				written++;
			}
			binary.resize(written);
	
			byteBuffer = &binary[0];
			bufferSize = written;			
			
		}
		
		//change endianness if necessary
		if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
		{
			if(element_size==4)
			{
				const Real* floatBuffer = reinterpret_cast<const Real*>(byteBuffer);

				if (bufferSize % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");

				Size floatCount = bufferSize / element_size;

				out.resize(floatCount);

				UInt i = 0;
				Int32* p = reinterpret_cast<Int32*> (byteBuffer);
				while(i < floatCount)
				{
					*p = endianize32_(*p);
					++p;
					++i;
				}
				out.assign(floatBuffer,floatBuffer+floatCount);	

			}
			else if(element_size==8)
			{
				const DoubleReal* floatBuffer = reinterpret_cast<const DoubleReal*>(byteBuffer);

				if (bufferSize % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");

				Size floatCount = bufferSize / element_size;

				out.resize(floatCount);

				UInt i = 0;
				Int64* p = reinterpret_cast<Int64*> (byteBuffer);
				while(i < floatCount)
				{
					*p = endianize64_(*p);
					++p;
					++i;
				}
				out.assign(floatBuffer,floatBuffer+floatCount);	
			}			
		}
		else
		{
			if(element_size==4)
			{
				const Real* floatBuffer = reinterpret_cast<const Real*>(byteBuffer);
				if (bufferSize % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");

				Size floatCount = bufferSize / element_size;

				out.resize(floatCount);

				std::copy(floatBuffer, floatBuffer+floatCount, out.begin());
			}
			else if(element_size==8)
			{
				const DoubleReal* floatBuffer = reinterpret_cast<const DoubleReal*>(byteBuffer);

				if (bufferSize % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");

				Size floatCount = bufferSize / element_size;

				out.resize(floatCount);

				std::copy(floatBuffer, floatBuffer+floatCount, out.begin());					
			}					
		}
		
	}
	
	template <typename ToType>
	void Base64::decode_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out)
	{
		out.clear();
		if (in == "") return;
	
		Size src_size = in.size();
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

		const Size element_size = sizeof(ToType);

		// enough for either float or double
		char element[8] = "\x00\x00\x00\x00\x00\x00\x00";

		if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
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
		for (Size i=0; i<src_size; i+=4)
		{
//		printf ("start: i=%d, offset %d\n", i, offset);

			// decode 4 Base64-Chars to 3 Byte
			a = old_decoder_[(int)in[i]-43]-62;
			b = old_decoder_[(int)in[i+1]-43]-62;
			if (i+1 >= src_size) b=0;
			element[offset] = (unsigned char) ((a<<2) | (b>>4));
			written++;
//		printf ("1: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % element_size == 0)
			{
				out.push_back(((ToType*)element)[0]);
				strcpy(element, "");
			}

			a = old_decoder_[(int)in[i+2]-43]-62;
			if (i+2 >= src_size) a=0;
			element[offset] = (unsigned char) (((b&15)<<4) | (a>>2));
			written++;
//		printf ("2: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % element_size == 0)
			{
				// debug: output float in binary
/*
				for (int sl = 0; sl != element_size; sl++)
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

			b = old_decoder_[(int)in[i+3]-43]-62;
			if (i+3 >= src_size) b=0;
			element[offset] = (unsigned char) (((a&3)<<6) | b);
			written++;
//		printf ("3: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % element_size == 0)
			{
				out.push_back(((ToType*)element)[0]);
				strcpy(element, "");
			}
		}
	}	

} //namespace OpenMS

#endif /* OPENMS_FORMAT_BASE64_H */
