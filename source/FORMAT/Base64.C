// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <QtCore/QList>
#include <QtCore/QString>

using namespace std;

namespace OpenMS
{

	const char Base64::encoder_[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	const char Base64::decoder_[] = "|$$$}rstuvwxyz{$$$$$$$>?@ABCDEFGHIJKLMNOPQRSTUVW$$$$$$XYZ[\\]^_`abcdefghijklmnopq";

	Base64::Base64()
	{
	}

	Base64::~Base64()
  {
	}

	void Base64::encodeStrings(std::vector<String>& in,String& out, bool zlib_compression)
	{
		out.clear();
		if (in.empty()) return;
		std::string str;
		std::string compressed;
		Byte* it;
		Byte* end;
		for(Size i = 0 ; i < in.size(); ++i )
		{
			str = str.append(in[i]);
			str.push_back('\0');
		}
		
		if(zlib_compression)
		{
			unsigned long sourceLen =	(unsigned long)str.size();
			unsigned long compressed_length = //compressBound((unsigned long)str.size());
					sourceLen + (sourceLen >> 12) + (sourceLen >> 14) + 11; // taken from zlib's compress.c, as we cannot use compressBound*
			
			int zlib_error;
			do
			{
      	compressed.resize(compressed_length);
      	zlib_error = compress(reinterpret_cast<Bytef *>(&compressed[0]),&compressed_length , reinterpret_cast<Bytef*>(&str[0]),(unsigned long) str.size());
       
        switch (zlib_error) 
        {
        	case Z_MEM_ERROR:
          	throw Exception::OutOfMemory(__FILE__,__LINE__,__PRETTY_FUNCTION__,compressed_length);
            break;
        	case Z_BUF_ERROR:
            compressed_length *= 2;
     		}
    	}while (zlib_error == Z_BUF_ERROR);
			
			if(zlib_error != Z_OK)
			{
				throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Compression error?");
			}
			
			it = reinterpret_cast<Byte*>(&compressed[0]);
			end =it + compressed_length;
			out.resize((Size)ceil(compressed_length/3.)*4); //resize output array in order to have enough space for all characters
		}
		else
		{
			out.resize((Size)ceil(str.size()/3.) * 4); //resize output array in order to have enough space for all characters
			it = reinterpret_cast<Byte*>(&str[0]);
			end = it + str.size();		
		}
			Byte* to = reinterpret_cast<Byte*>(&out[0]);			
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
	
	void Base64::decodeStrings(const String& in, std::vector<String>& out, bool zlib_compression)
	{
		out.clear();
		if (in == "") return;
		QByteArray base64_uncompressed;
		QByteArray herewego = QByteArray::fromRawData(in.c_str(), (int) in.size());
		base64_uncompressed = QByteArray::fromBase64(herewego);
		if(zlib_compression)
		{		
			QByteArray czip;
			czip.resize(4);
			czip[0] = (base64_uncompressed.size() & 0xff000000) >> 24;
			czip[1] = (base64_uncompressed.size() & 0x00ff0000) >> 16;
			czip[2] = (base64_uncompressed.size() & 0x0000ff00) >> 8;
			czip[3] = (base64_uncompressed.size()& 0x000000ff);
			czip += base64_uncompressed;
			base64_uncompressed = qUncompress(czip);

			if(base64_uncompressed.isEmpty())
			{
				throw Exception::ConversionError (__FILE__,__LINE__,__PRETTY_FUNCTION__,"Decompression error?");
			}
		}
		QList<QByteArray> null_strings = base64_uncompressed.split('\0');
		for(QList<QByteArray>::iterator it = null_strings.begin(); it < null_strings.end();++it)
		{
			if(!it->isEmpty())
			{
				out.push_back(QString(*it).toStdString());
			}
		}
	}	
	

}//end OpenMS

