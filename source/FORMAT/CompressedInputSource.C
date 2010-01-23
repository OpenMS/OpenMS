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
#include <OpenMS/FORMAT/CompressedInputSource.h>
#include <OpenMS/FORMAT/GzipInputStream.h>
#include <OpenMS/FORMAT/Bzip2InputStream.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <xercesc/internal/MemoryManagerImpl.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLUniDefs.hpp>


using namespace xercesc;
namespace OpenMS
{

	CompressedInputSource::CompressedInputSource(const String& file_path,const char* header, MemoryManager* const manager)
   : xercesc::InputSource(manager)
	{
    	if(sizeof(header)/sizeof(char)  > 1)
    	{
    		head_[0] = header[0];
    		head_[1] = header[1];
    	}
    	else
    	{
    		head_[0] = '\0';
    		head_[1] = '\0';
    	}
    	//
    	//  If the path is relative, then complete it acording to the current
    	//  working directory rules of the current platform. Else, just take
    	//  it as is.
    	//
    	Internal::StringManager strman;
    	XMLCh* file = strman.convert(file_path.c_str());
    	if (xercesc::XMLPlatformUtils::isRelative(file, manager))
    	{
        XMLCh* curDir = xercesc::XMLPlatformUtils::getCurrentDirectory(manager);

        XMLSize_t curDirLen = XMLString::stringLen(curDir);
        XMLSize_t filePathLen = XMLString::stringLen(file);
        XMLCh* fullDir = (XMLCh*) manager->allocate
        (
            (curDirLen + filePathLen + 2) * sizeof(XMLCh)
        );//new XMLCh [ curDirLen + filePathLen + 2];

        XMLString::copyString(fullDir, curDir);
        fullDir[curDirLen] = chForwardSlash;
        XMLString::copyString(&fullDir[curDirLen+1], file);
        
        XMLPlatformUtils::removeDotSlash(fullDir, manager);
        XMLPlatformUtils::removeDotDotSlash(fullDir, manager);

        setSystemId(fullDir);

        manager->deallocate(curDir);//delete [] curDir;
        manager->deallocate(fullDir);//delete [] fullDir;
    	}
     	else
    	{
        XMLCh* tmpBuf = XMLString::replicate(file, manager);
        XMLPlatformUtils::removeDotSlash(tmpBuf, manager);
        setSystemId(tmpBuf);
      	manager->deallocate(tmpBuf);//delete [] tmpBuf;
  	  }
	}

	CompressedInputSource::CompressedInputSource(const XMLCh* const file,const char* header, MemoryManager* const manager)
   : xercesc::InputSource(manager)
	{
    	if(sizeof(header)/sizeof(char) > 1)
    	{
    		head_[0] = header[0];
    		head_[1] = header[1];
    	}
    	else
    	{
    		head_[0] = '\0';
    		head_[1] = '\0';
    	}
    	//
    	//  If the path is relative, then complete it acording to the current
    	//  working directory rules of the current platform. Else, just take
    	//  it as is.
    	//
    	if (xercesc::XMLPlatformUtils::isRelative(file, manager))
    	{
        XMLCh* curDir = xercesc::XMLPlatformUtils::getCurrentDirectory(manager);

        XMLSize_t curDirLen = XMLString::stringLen(curDir);
        XMLSize_t filePathLen = XMLString::stringLen(file);
        XMLCh* fullDir = (XMLCh*) manager->allocate
        (
            (curDirLen + filePathLen + 2) * sizeof(XMLCh)
        );//new XMLCh [ curDirLen + filePathLen + 2];

        XMLString::copyString(fullDir, curDir);
        fullDir[curDirLen] = chForwardSlash;
        XMLString::copyString(&fullDir[curDirLen+1], file);
        
        XMLPlatformUtils::removeDotSlash(fullDir, manager);
        XMLPlatformUtils::removeDotDotSlash(fullDir, manager);

        setSystemId(fullDir);

        manager->deallocate(curDir);//delete [] curDir;
        manager->deallocate(fullDir);//delete [] fullDir;
    	}
     	else
    	{
        XMLCh* tmpBuf = XMLString::replicate(file, manager);
        XMLPlatformUtils::removeDotSlash(tmpBuf, manager);
        setSystemId(tmpBuf);
      	manager->deallocate(tmpBuf);//delete [] tmpBuf;
  	  }
	}

	CompressedInputSource::~CompressedInputSource()
	{
	}
	

BinInputStream* CompressedInputSource::makeStream() const
{
		if(head_[0] == 'B' && head_[1] =='Z' )
    {
    	Bzip2InputStream* retStrm = new Bzip2InputStream(Internal::StringManager().convert(getSystemId()));
    	    if (!retStrm->getIsOpen())
    {
       delete retStrm;
        return 0;
    }
    return retStrm;
    }
    else /* 	(bz[0] == g1 && bz[1] == g2), where char g1 = 0x1f and char g2 = 0x8b */
    {
    	GzipInputStream* retStrm = new GzipInputStream(Internal::StringManager().convert(getSystemId()));
     if(!retStrm->getIsOpen())
  	  {
 	      delete retStrm;
  	     return 0;
	    }
   	 return retStrm;
    }

}	
	
	
} // namespace OpenMS
