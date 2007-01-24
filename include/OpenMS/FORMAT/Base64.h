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

#ifndef OPENMS_FORMAT_BASE64_H
#define OPENMS_FORMAT_BASE64_H

#include <cmath>
#include <arpa/inet.h>
#include <OpenMS/CONCEPT/Types.h>

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
		
		Base64 supports endian-conversion only for little-endian machines:
		to handle big endian data call decode***Corrected() and encode***Corrected().
		
		@todo Speed up by writing to MSExperiment directly (Marc)
  */
  class Base64
  {
  public:
    /// default constructor
    Base64();

    /// Destructor
    virtual ~Base64();

		/// Return size of output buffer
		Size getOutputBufferSize();

		/** 
			@brief Set size of output buffer.

		 	Given size adapts automaticly to multiple of 3.
		 	
			Use this function if you know the maximal used buffer size in advance
			to avoid successive buffer allocation.
		*/
		void setOutputBufferSize(Size s);


    /// Decode given Base64-String of size @p size
		char* decode(const char* src, Size size);
		
    /// Dncode given Base64-String of size @p size
		char* encode(const char* src, Size size);

		/**@name Handling of 32 bit real type */
    //@{

    /// Decode given Base64-String of size @p size to array of floats, each float corrected from network to host byte order
		float* decodeFloatCorrected(const char* src, Size size);

    /// decode given Base64-String of size @p size to array of floats
		float* decodeFloat(const char* src, Size size);

    /// return internal input buffer to fill with @p size floats
		float* getFloatBuffer(Size size);

    /// encode internal input buffer (fill with getFloatBuffer()) to Base64-String after conversion to network byte order
		char* encodeFloatCorrected();

    /// encode internal input buffer (fill with getFloatBuffer()) to Base64-String
		char* encodeFloat();
		//@}
		
		
		/**@name Handling of 64 bit real type */
    //@{

    /// decode given Base64-String of size @p size to array of doubles, each double corrected from network to host byte order
		double* decodeDoubleCorrected(const char* src, Size size);

    /// decode given Base64-String of size @p size to array of doubles
		double* decodeDouble(const char* src, Size size);

    /// return internal input buffer to fill with @p size doubles
		double* getDoubleBuffer(Size size);

    /// encode internal input buffer (fill with getDoubleBuffer()) to Base64-String after conversion to network byte order
		char* encodeDoubleCorrected();

    /// encode internal input buffer (fill with getDoubleBuffer()) to Base64-String
		char* encodeDouble();
		//@}

  private:
		/// input buffer
		char* in_buffer_;
		/// output buffer
		char* out_buffer_;
		/// Length of input buffer
		Size in_length_;
		/// Length of output buffer
		Size out_length_;
		/// Length of input buffer
		Size ibuffer_size_;
		
		/// Adapts input buffer size to new value of in_length_.
		void adaptInputBuffer_();
		
		static const char encoder_[];
		static const char decoder_[];
  };



} //namespace OpenMS

#endif /* OPENMS_FORMAT_BASE64_H */
