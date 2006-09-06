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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_BASE64_H
#define OPENMS_FORMAT_BASE64_H

#include <cmath>
#include <arpa/inet.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include<iostream>

namespace OpenMS
{
  /**
  	@brief Class to encode and decode Base64
  	
  	@ingroup Format
		
		Use the encoding-functions only after filling the internal buffer. To get the
		buffer, call get***Buffer(). Furthermore, do not use the same or different
		encoding functions multiple times without refilling the internal buffer
		inbetween.
		Base64 supports two precisions 32 Bit (float) and 64 Bit (double).
		<br>
		Base64 supports endian-conversion only for little-endian machines:
		to handle big endian data call decode***Corrected() and encode***Corrected().
  */
  class Base64
  {
  public:
    /// default constructor
    Base64();

    /// destructor
    virtual ~Base64();

		/// return size of output buffer
		Size getOutputBufferSize();

		/** @brief set size of output buffer

			 	Given size adapts automaticly to multiple of 3 <br>
				use this function if you know the maximal used buffer size in advance
				to avoid successive buffer allocation
		*/
		void setOutputBufferSize(Size s);


    /// decode given Base64-String of size @p size
		inline char* decode(const char* src, const Size size)
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

    /// encode given Base64-String of size @p size
		inline char* encode(const char* src, const Size size)
		{
			if (size<3)
			{
				out_buffer_[0] = '\0';
				return out_buffer_;
			}
			UnsignedInt padding = 0;
			if (size%3 == 2) padding=1; 
			if (size%3 == 1) padding=2;
			
			Size dest_size = (size+padding)/3*4;
			if (dest_size > out_length_) setOutputBufferSize(dest_size);

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

		/**@name Handling of 32 bit real type */
    //@{

    /// decode given Base64-String of size @p size to array of floats,
		///	each float corrected from network to host byte order
		float* decodeFloatCorrected(const char* src, const Size size);

    /// decode given Base64-String of size @p size to array of floats
		float* decodeFloat(const char* src, const Size size);

    /// return internal input buffer to fill with @p size floats
		float* getFloatBuffer(const Size size);

    /// encode internal input buffer (fill with getFloatBuffer()) to Base64-String
		/// after conversion to network byte order
		char* encodeFloatCorrected();

    /// encode internal input buffer (fill with getFloatBuffer()) to Base64-String
		char* encodeFloat();

		//@}

		/**@name Handling of 64 bit real type */
    //@{

    /// decode given Base64-String of size @p size to array of doubles,
		///	each double corrected from network to host byte order
		double* decodeDoubleCorrected(const char* src, const Size size);

    /// decode given Base64-String of size @p size to array of doubles
		double* decodeDouble(const char* src, const Size size);

    /// return internal input buffer to fill with @p size doubles
		double* getDoubleBuffer(const Size size);

    /// encode internal input buffer (fill with getDoubleBuffer()) to Base64-String
		/// after conversion to network byte order
		char* encodeDoubleCorrected();

    /// encode internal input buffer (fill with getDoubleBuffer()) to Base64-String
		char* encodeDouble();

		//@}

  private:
		/// input buffer
		char* in_buffer_;
		/// output buffer
		char* out_buffer_;
		/// actual length of input buffer in Bytes
		Size in_length_;
		/// actual length of output buffer in Bytes
		Size out_length_;

		static const char encoder_[];
		static const char decoder_[];
  };



} //namespace OpenMS

#endif /* OPENMS_FORMAT_BASE64_H */
