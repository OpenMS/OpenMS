// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/Bzip2InputStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace xercesc;

namespace OpenMS
{
	Bzip2InputStream::Bzip2InputStream(const   String&   file_name)
	:bzip2_(new Bzip2Ifstream(file_name.c_str())),file_current_index_(0)
	{
	}

	Bzip2InputStream::Bzip2InputStream(const   char* file_name)
	:bzip2_(new Bzip2Ifstream(file_name)),file_current_index_(0)
  {
  }
	
/*	Bzip2InputStream::Bzip2InputStream()
	:bzip2_(NULL)
	{
	
	}*/
	
	Bzip2InputStream::~Bzip2InputStream()
	{
		delete bzip2_;
	}
	XMLSize_t Bzip2InputStream::readBytes(XMLByte* const  to_fill, const XMLSize_t  max_to_read)
	{
    //  Figure out whether we can really read. 
   if(bzip2_->streamEnd())
   {
   	return 0;
   }
   
   unsigned char* fill_it = static_cast<unsigned char*>(to_fill);
   XMLSize_t actual_read = (XMLSize_t) bzip2_->read((char*)fill_it, static_cast<const size_t>(max_to_read));
   file_current_index_ += actual_read;
   return actual_read;
	}

	const XMLCh* Bzip2InputStream::getContentType() const
	{
    return 0;
	}	
	
} // namespace OpenMS
