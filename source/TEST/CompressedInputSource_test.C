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
// $Maintainer: David Wojnar$
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/CompressedInputSource.h>
#include <OpenMS/FORMAT/GzipInputStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
using namespace OpenMS;


///////////////////////////

START_TEST(CompressedInputSource, "$Id$")

CompressedInputSource* ptr = 0;
CompressedInputSource* nullPointer = 0;
START_SECTION(CompressedInputSource(const   String& file_path, const char * header,xercesc::MemoryManager* const manager = xercesc::XMLPlatformUtils::fgMemoryManager))
	char header[2];
	header[0] = 'B';
	header[1] = 'Z';
	ptr = new CompressedInputSource(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"),header);
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION
START_SECTION((~CompressedInputSource()))
	delete ptr;
END_SECTION

START_SECTION(CompressedInputSource(const XMLCh *const file_path, const char *header, xercesc::MemoryManager *const manager=xercesc::XMLPlatformUtils::fgMemoryManager))
		char header[2];
	header[0] = 'B';
	header[1] = 'Z';
	String filename(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	ptr = new CompressedInputSource(Internal::StringManager().convert(filename.c_str()),header);
	TEST_NOT_EQUAL(ptr, nullPointer)
	delete ptr;
END_SECTION


START_SECTION(virtual xercesc::BinInputStream* makeStream() const)
	char header[2];
	header[0] = 'B';
	header[1] = 'Z';
	CompressedInputSource source(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist"),header);
	TEST_EXCEPTION(Exception::FileNotFound,source.makeStream())
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
