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
#include <string>
#include <vector>
#include <cstring>


#include <QtCore/QString>
namespace OpenMS
{
  /**
  	@brief Class to encode and decode Base64
		
		Base64 supports two precisions: 32 Bit (float) and 64 Bit (double).
		
		@todo Check assign code / replace unsigned long long (HiWi)
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
		void encode(std::vector<FromType>& in, ByteOrder to_byte_order, std::string& out, bool zlib_compression= false);

		/**
		   @brief Decodes a Base64 string to a vector of floating point numbers
		     @note @p in will be emtpy after this method
		*/

		template <typename ToType>
		void decode(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out,bool zlib_compression = false);
		

		
		inline unsigned int endianize32(unsigned int n);


		inline unsigned long long endianize64(unsigned long long n);
		
		union endian64
		{
			DoubleReal value;
			char bytes[sizeof(DoubleReal)];
		};

		inline double SwapBytes(DoubleReal& value,unsigned int size)
		{
			endian64 in;
			endian64 out;
	
			in.value = value;
	
			for (unsigned int i = 0; i < size / 2; ++i)
			{
				out.bytes[i] = in.bytes[size - 1 - i];
				out.bytes[size - 1 - i] = in.bytes[i];
			}
	
			return out.value;
		}

		union endian32
		{
			Real value;
			char bytes[sizeof(Real)];
		};
		
		inline double SwapBytes(Real& value,unsigned int size);

	
	private:
		static const char encoder_[];
		static const char old_decoder_[];
		char decoder_[256];
		bool tableInitialized_;
		//this is the old decode function which doesn't suport compression but is faster.
		template <typename ToType>
	  void decode_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out);
	 //new decode function does support compression and no compression but is only used if compression needed because of performance
			template <typename ToType>
		  void newdecode_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out,bool zlib_compression = false);
  };

	inline unsigned int Base64::endianize32(unsigned int n)
	{
		return ((n&0xff)<<24) | ((n&0xff00)<<8) | ((n&0xff0000)>>8) | ((n&0xff000000)>>24);
		
	}
	
	inline unsigned long long Base64::endianize64(unsigned long long n)
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
	
	inline double Base64::SwapBytes(Real& value,unsigned int size)
	{
			endian32 in;
			endian32 out;
			in.value = value;
		
			for (unsigned int i = 0; i < size / 2; ++i)
			{
				out.bytes[i] = in.bytes[size - 1 - i];
				out.bytes[size - 1 - i] = in.bytes[i];
			}
	
			return out.value;
		}
	template <typename FromType>
	void Base64::encode(std::vector<FromType>& in, ByteOrder to_byte_order, std::string& out, bool zlib_compression )
	{
		out.clear();
		if (in.size() == 0) return;

		
		void* byteBuffer = reinterpret_cast<void*>(& in[0]);
				
		Size element_size = sizeof (FromType);
		Size size = element_size * in.size();
			
		
		if ((OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BYTEORDER_LITTLEENDIAN) ||
			(!OPENMS_IS_BIG_ENDIAN && to_byte_order == Base64::BYTEORDER_BIGENDIAN))
		{
			if(element_size == sizeof(Real))
			{
				//unsigned int* p = reinterpret_cast<unsigned int*>(&in[0]);
				UInt i = 0;
				while(i < in.size())
				{
					in[i] = SwapBytes(in[i],sizeof(in[i]));
					i++;
				}
				//while (p != p+in.size())
				//{*p++ = endianize32(*p++);  }
				//OpenMS::transform(p, p+in.size(), p , endianize32);
			}
			else
			{
				//unsigned long long* p = reinterpret_cast<unsigned long long*>(&in[0]);
				UInt i=0;
				while (i < in.size())
				{in[i] = SwapBytes(in[i],sizeof(in[i]));i++;  }			
			}
		}
		
		if (zlib_compression)
		{
			QByteArray original = QByteArray::fromRawData((const char*)byteBuffer, (int) size);
			QByteArray compressed = qCompress((uchar*)original.data(),original.size());
			QByteArray extern_compressed = compressed.right(compressed.size() - 4);			
			QByteArray base64_compressed = extern_compressed.toBase64();
	
			out = QString(base64_compressed).toStdString();
		}
		else
		{
			out.resize((Size)ceil(size/3.) * 4); //resize  out  in order to have enough space for all characters
			char* to = &out[0];
			Byte* it = (Byte*)byteBuffer;
			Byte* end = it + size;
			Size written = 0;

			while (it!=end)
			{
				Int int24bit = 0;
				Int paddingCount = 0;

				// construct 24-bit integer from 3 bytes
				for (UInt i=0; i<3; i++)
				{
					if (it!=end)
						int24bit |= *it++<<((2-i)*8);
					else
						paddingCount++;
				}

				// write out 4 characters
				for (Int i=3; i>=0; i--)
				{
					to[i] = encoder_[int24bit & 0x3F];
					int24bit >>= 6;
				}

				// fixup for padding
				if (paddingCount > 0)
				to[3] = '=';
				if (paddingCount > 1)
				to[2] = '=';

				to += 4;
				written += 4;
			}

			out.resize(written);//no more space is needed
		}
	}
	
	template<typename ToType>
	void Base64::decode(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out, bool zlib_compression)
	{
		if(zlib_compression)
		{
			newdecode_(in,from_byte_order,out,zlib_compression);
		}
		else
		{
			decode_(in,from_byte_order,out);
		}
	}
	
	template <typename ToType>
	void Base64::newdecode_(const std::string& in, ByteOrder from_byte_order, std::vector<ToType>& out,bool zlib_compression)
	{
		if(zlib_compression)
		out.clear();
		if (in == "") return;
		Size src_size = in.size();
		
		void* byteBuffer;
		Size bufferSize;		
		std::vector<unsigned char> binary;
		UInt element_size = sizeof(ToType);
		
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

			if(!tableInitialized_)
			{
				for (UInt i=0; i<64; i++)
				{
					decoder_[static_cast<int>(encoder_[i])] = static_cast<char>(i);
				}
				tableInitialized_ = true;
			}

			Byte* it = (Byte*)&in[0];
			Byte* end = it + src_size;
			Byte* result = (Byte*)&binary[0];
			Size written = 0;

			while (it!=end)
			{
				Int int24bit = 0;
				Int paddingCount = 0;

				// construct 24-bit integer from 4 characters
				for (Int i=0; i<4 && it!=end; i++, it++)
				{
					if (*it != '=')
					{
						int24bit |= decoder_[*it]<<((3-i)*6);
					}
					else
					{
						paddingCount++;
					}
				}

				// write out bytes
				for (Int i=0; i<3-paddingCount; i++)
				{
					Byte temp = static_cast<Byte>(int24bit>>((2-i)*8));
					*result++ = temp;
					int24bit ^= temp<<((2-i)*8);
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
		if ((OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_LITTLEENDIAN) ||
			(!OPENMS_IS_BIG_ENDIAN && from_byte_order == Base64::BYTEORDER_BIGENDIAN))
		{
			if(element_size == sizeof(Real))
			{

				const Real* floatBuffer = reinterpret_cast<const Real*>(byteBuffer);

				if (bufferSize % sizeof(Real) != 0) 
					throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");

				Size floatCount = bufferSize / sizeof(Real);

				out.resize(floatCount);

				UInt i = 0;
				unsigned int* p = reinterpret_cast<unsigned int*> (byteBuffer);
				while(i < floatCount)
				{
					*p = endianize32(*p);
					++p;
					++i;
				}
				out.assign(floatBuffer,floatBuffer+floatCount);	

			}
			else
			{
				const DoubleReal* floatBuffer = reinterpret_cast<const DoubleReal*>(byteBuffer);

				if (bufferSize % sizeof(DoubleReal) != 0) 
					throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");

				Size floatCount = bufferSize / sizeof(DoubleReal);

				out.resize(floatCount);

				UInt i = 0;
				unsigned long long* p = reinterpret_cast<unsigned long long*> (byteBuffer);
				while(i < floatCount)
				{
					*p = endianize64(*p);
					++p;
					++i;
				}
				out.assign(floatBuffer,floatBuffer+floatCount);	
			}			
		}
		else
		{
			if(element_size == sizeof(Real))
			{
				
				const Real* floatBuffer = reinterpret_cast<const Real*>(byteBuffer);
				if (bufferSize % sizeof(Real) != 0) 
					throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");

				Size floatCount = bufferSize / sizeof(Real);

				out.resize(floatCount);

				std::copy(floatBuffer, floatBuffer+floatCount, out.begin());
			}
			else
			{
				const DoubleReal* floatBuffer = reinterpret_cast<const DoubleReal*>(byteBuffer);

				if (bufferSize % sizeof(DoubleReal) != 0) 
					throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Bad BufferCount?");

				Size floatCount = bufferSize / sizeof(DoubleReal);

				out.resize(floatCount);

				std::copy(floatBuffer, floatBuffer+floatCount, out.begin());					
			}					
		}
		
	}
	
	/// Decodes a Base64 string to a float vector
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

			if (written % sizeof(ToType) == 0)
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

			b = old_decoder_[(int)in[i+3]-43]-62;
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
