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
			
			enum DataType
			{
				INTEGER,///< Integer type
				FLOAT,	///< Floating point type
				STRING ///< String type
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
			void decode(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out, bool zlib_compression = false,DataType data_type= FLOAT);

			/**
				@brief Encodes a vector of strings to a Base64 string
				
				You can specify  zlib-compresseion-
				
				@note @p in will be emtpy after this method
			*/
			void encodeStrings(std::vector<String>& in, std::string& out, bool zlib_compression= false);
			
			/**
				@brief Decodes a Base64 string to a vector of (null-terminated) strings
	
				You have to specify whether the Base64 string is zlib-compressed.
	
				@note @p in will be emtpy after this method
			*/		
			void decodeStrings(const std::string& in, std::vector<String>& out, bool zlib_compression = false);
		private:
			
			///Internal class needed for type-punning
			union Reinterpreter64_
			{
				DoubleReal f;
				Int64 i;
			};

			///Internal class needed for type-punning
			union Reinterpreter32_
			{
				Real f;
				Int32 i;
			};

			static const char encoder_[];
			static const char decoder_[];
			/**
				@brief Decodes a Base64 string to a vector of floating point numbers
				
				@note this function is not able to decode compresed data
			*/
			template <typename ToType>
			void decodeUncompressed_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out, DataType data_type);
		 
			///Decodes a compressed Base64 string to a vector of floating point numbers 
			template <typename ToType>
			void decodeCompressed_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out,DataType data_type);

			///Decodes a compressed Base64 string to a vector of floating point numbers representing integers
			template <typename ToType>
			void decodeUncompressedInteger_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out);			
  };
			///Endianizes a 32 bit type from big endian to litte endian and vice versa
	inline Int32 endianize32(Int32& n)
	{
		return ((n&0xff)<<24) | ((n&0xff00)<<8) | ((n&0xff0000)>>8) | ((n&0xff000000)>>24);
	}
	///Endianizes a 64 bit type from  big endian to litte endian and vice versa
	inline Int64 endianize64(Int64& n)
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
					Reinterpreter32_ tmp;
					tmp.f = in[i];
					tmp.i = endianize32(tmp.i);
					in[i] = tmp.f;
				}
			}
			else if (element_size == 8)
			{
				for (Size i = 0; i < in.size(); ++i)
				{
					Reinterpreter64_ tmp;
					tmp.f = in[i];
					tmp.i = endianize64(tmp.i);
					in[i] = tmp.f;
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
	void Base64::decode(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out, bool zlib_compression, DataType data_type)
	{
		if(zlib_compression)
		{
			decodeCompressed_(in,from_byte_order,out,data_type);
		}
		else
		{
			decodeUncompressed_(in,from_byte_order,out,data_type);
		}
	}
	
	template <typename ToType>
	void Base64::decodeCompressed_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out,DataType data_type)
	{
		out.clear();
		if (in == "") return;
		
		void* byte_buffer;
		Size buffer_size;		
		std::vector<unsigned char> binary;
		const Size element_size = sizeof(ToType);
		
		String decompressed;

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
			
			std::copy(base64_uncompressed.begin(),base64_uncompressed.end(),decompressed.begin());

			byte_buffer = reinterpret_cast<void*>(&decompressed[0]);
			buffer_size = decompressed.size();
		
		//change endianness if necessary
		if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) || (!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
		{
			if(element_size==4 && data_type == FLOAT)
			{
				const Real* float_buffer = reinterpret_cast<const Real*>(byte_buffer);
				if (buffer_size % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");
				Size float_count = buffer_size / element_size;
				Int32* p = reinterpret_cast<Int32*> (byte_buffer);
				std::transform(p,p+float_count,p,endianize32);
				out.assign(float_buffer,float_buffer+float_count);	
			}
			else if(element_size==4 && data_type == INTEGER)
			{
				const Int32* float_buffer = reinterpret_cast<const Int32*>(byte_buffer);
				if (buffer_size % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");
				Size float_count = buffer_size / element_size;
				Int32* p = reinterpret_cast<Int32*> (byte_buffer);
				std::transform(p,p+float_count,p,endianize32);

				out.resize(float_count);
				// do NOT use assign here, as it will give a lot of type conversion warnings on VS compiler	
				for (Size i=0;i<float_count;++i)
				{
					out[i]=(ToType) *float_buffer;
					++float_buffer;
				}
			}
			else if(element_size==8 && data_type == FLOAT)
			{
				const DoubleReal* float_buffer = reinterpret_cast<const DoubleReal*>(byte_buffer);

				if (buffer_size % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");

				Size float_count = buffer_size / element_size;

				Int64* p = reinterpret_cast<Int64*> (byte_buffer);
				std::transform(p,p+float_count,p,endianize64);				
				
				out.resize(float_count);
				// do NOT use assign here, as it will give a lot of type conversion warnings on VS compiler
				for (Size i=0;i<float_count;++i)
				{
					out[i]=(ToType) *float_buffer;
					++float_buffer;
				}
			}
			else if(element_size==8 && data_type == INTEGER)
			{
				const Int64* float_buffer = reinterpret_cast<const Int64*>(byte_buffer);

				if (buffer_size % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");

				Size float_count = buffer_size / element_size;

				Int64* p = reinterpret_cast<Int64*> (byte_buffer);
				std::transform(p,p+float_count,p,endianize64);				
				
				out.resize(float_count);
				// do NOT use assign here, as it will give a lot of type conversion warnings on VS compiler
				for (Size i=0;i<float_count;++i)
				{
					out[i]=(ToType) *float_buffer;
					++float_buffer;
				}
				
			}				
		}
		else
		{
			if(element_size==4 && data_type == FLOAT)
			{
				const Real* float_buffer = reinterpret_cast<const Real*>(byte_buffer);
				if (buffer_size % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount while decoding?");

				Size float_count = buffer_size / element_size;
				out.assign(float_buffer, float_buffer+float_count);
			}
			else if(element_size==4 && data_type == INTEGER)
			{
				const Int* float_buffer = reinterpret_cast<const Int*>(byte_buffer);
				if (buffer_size % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount while decoding?");

				Size float_count = buffer_size / element_size;
				out.resize(float_count);
				// do NOT use assign here, as it will give a lot of type conversion warnings on VS compiler			
				for (Size i=0;i<float_count;++i)
				{
					out[i]=(ToType) *float_buffer;
					++float_buffer;
				}
			}
			else if(element_size==8 && data_type == FLOAT)
			{
				const DoubleReal* float_buffer = reinterpret_cast<const DoubleReal*>(byte_buffer);

				if (buffer_size % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount while decoding?");

				Size float_count = buffer_size / element_size;
				out.resize(float_count);
				// do NOT use assign here, as it will give a lot of type conversion warnings on VS compiler			
				for (Size i=0;i<float_count;++i)
				{
					out[i]=(ToType) *float_buffer;
					++float_buffer;
				}
			}	
			else if(element_size==8 && data_type == INTEGER)
			{
				const Int64* float_buffer = reinterpret_cast<const Int64*>(byte_buffer);

				if (buffer_size % element_size != 0) throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount while decoding?");

				Size float_count = buffer_size / element_size;
				out.resize(float_count);
				// do NOT use assign here, as it will give a lot of type conversion warnings on VS compiler			
				for (Size i=0;i<float_count;++i)
				{
					out[i]=(ToType) *float_buffer;
					++float_buffer;
				}
			}
		}
		
	}
	
	template <typename ToType>
	void Base64::decodeUncompressed_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out,DataType data_type)
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
			a = decoder_[(int)in[i]-43]-62;
			b = decoder_[(int)in[i+1]-43]-62;
			if (i+1 >= src_size) b=0;
			element[offset] = (unsigned char) ((a<<2) | (b>>4));
			written++;
//		printf ("1: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % element_size == 0)
			{
				if(data_type == FLOAT)
				{
					ToType* to_type = reinterpret_cast<ToType*>(&element[0]);
					out.push_back((*to_type));
				}
				else
				{
					ToType float_value;
					if(element_size==4)
					{
					Int32* value = reinterpret_cast<Int32*>(&element[0]);
					float_value = (ToType)*value;
					}
					else
					{
						Int64* value = reinterpret_cast<Int64*>(&element[0]);
						float_value = (ToType)*value;
					}
					out.push_back(float_value);			
				}
				strcpy(element, "");
			}

			a = decoder_[(int)in[i+2]-43]-62;
			if (i+2 >= src_size) a=0;
			element[offset] = (unsigned char) (((b&15)<<4) | (a>>2));
			written++;
//		printf ("2: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % element_size == 0)
			{
				// debug: output float in binary
/*
				for (unsigned int sl = 0; sl != element_size; sl++)
				{
					for (int sl2 = 128; sl2 >= 1; sl2 /= 2)
					{
						std::cout << (element[sl] & sl2);
					}
					std::cout << " ";
				}
*/
				if(data_type == FLOAT)
				{
					ToType* to_type = reinterpret_cast<ToType*>(&element[0]);
					out.push_back((*to_type));
				}
				else
				{
					ToType float_value;
					if(element_size==4)
					{
					Int32* value = reinterpret_cast<Int32*>(&element[0]);
					float_value = (ToType)*value;
					}
					else
					{
						Int64* value = reinterpret_cast<Int64*>(&element[0]);
						float_value = (ToType)*value;
					}
					out.push_back(float_value);				
				}
				strcpy(element, "");
			}

			b = decoder_[(int)in[i+3]-43]-62;
			if (i+3 >= src_size) b=0;
			element[offset] = (unsigned char) (((a&3)<<6) | b);
			written++;
//		printf ("3: i=%d, offset %d, wrote %d\n", i, offset, element[offset]);
			offset = (offset + inc) % element_size;

			if (written % element_size == 0)
			{
				if(data_type == FLOAT)
				{
					ToType* to_type = reinterpret_cast<ToType*>(&element[0]);
					out.push_back((*to_type));
				}
				else
				{
					ToType float_value;
					if(element_size==4)
					{
					Int32* value = reinterpret_cast<Int32*>(&element[0]);
					float_value = (ToType)*value;
					}
					else
					{
						Int64* value = reinterpret_cast<Int64*>(&element[0]);
						float_value = (ToType)*value;
					}
					out.push_back(float_value);	
				}
				strcpy(element, "");
			}
		}
	}	
	
} //namespace OpenMS

#endif /* OPENMS_FORMAT_BASE64_H */
