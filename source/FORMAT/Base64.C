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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <QtCore/QList>

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
		if (in.size() == 0) return;
		std::string str;
		for(Size i = 0 ; i < in.size(); ++i )
		{
			str = str.append(in[i]);
			str.push_back('\0');
		}

		QByteArray original = QByteArray::fromRawData(str.c_str(),(int) str.size());
		QByteArray base64_compressed;
		if (zlib_compression)
		{
			QByteArray compressed = qCompress((uchar*)original.data(),original.size());
			QByteArray extern_compressed = compressed.right(compressed.size() - 4);			
			base64_compressed = extern_compressed.toBase64();
		}
		//encode without compression
		else
		{
			base64_compressed = original.toBase64();
		}	
		out = QString(base64_compressed).toStdString();
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
	
}

