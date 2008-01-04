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
// $Maintainer: Chris Bielow  $
// --------------------------------------------------------------------------

#ifndef OPENMS_SYSTEM_EXTERNALALLOCATORUNIQUE_H
#define OPENMS_SYSTEM_EXTERNALALLOCATORUNIQUE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>



namespace OpenMS
{

  /**
    @brief Unique data of External allocator. This data will be held by a boost-sharedPtr.
    
    @ingroup System
  */
  class ExternalAllocatorUnique {
  
    protected:

      /// name of temporary file
      String filename_;
      /// size of temporary file
      off64_t filesize_;     
      /// next byte position in file where the next mapping is scheduled
      off64_t nextfree_;
      
      /// filehandle to swap file
      #ifdef OPENMS_WINDOWSPLATFORM
        HANDLE mmap_handle_;
      #else
        int mmap_handle_;      
      #endif
      
      /// just for informational purposes: how many bytes are mapped
      off64_t totalmappingsize_;
      
      
      /// freed blocks which shall be reused before opening a new one (used as a stack) (TODO: use it!)
      std::vector < std::pair <off64_t, off64_t > > freeblock_;


    public:  
    
    /* constructors and destructor
      */
         
      // "real" ctor that we will be using
      ExternalAllocatorUnique(const String &filename, const off64_t &filesize)
      :
        filename_(filename),
        filesize_(filesize),
        nextfree_(0),
        totalmappingsize_(0),
        freeblock_()
      {
        #ifdef DEBUG_ALLOC      
        std::cout << "--- 2-tuple Ctor called \n file:: " << filename_ << " size:: " << filesize_ << std::endl;
        #endif          
        
        String unique_filename = filename;
        // if file exists, another mapping is probably in place there.. we do not want to override that
        while (File::exists(unique_filename))
        {
          unique_filename = filename + std::rand();
        }
        filename_ = unique_filename;
        
        // handle to swap file (create swap file as well)
        mmap_handle_ = File::getSwapFileHandle(filename_, filesize_, true);
                
      }
      
      /* copy C'tor */
      ExternalAllocatorUnique(const ExternalAllocatorUnique& rhs) throw() 
      :
        filename_(rhs.filename_),
        filesize_(rhs.filesize_),
        nextfree_(rhs.nextfree_),
        mmap_handle_(rhs.mmap_handle_),
        totalmappingsize_(rhs.totalmappingsize_),
        freeblock_(rhs.freeblock_)
      {
        #ifdef DEBUG_ALLOC      
        std::cerr << "--- Copy Ctor called with nextfree_ " << nextfree_ << "\n";
        #endif          
      }
      
      
      /* destructor */
      ~ExternalAllocatorUnique() throw() 
      {
        #ifdef DEBUG_ALLOC      
        std::cerr << "--- ~ Destructor called \n";
        #endif          
        
        File::closeSwapFileHandle(mmap_handle_);
        
        if (!File::remove(filename_))
        {
          #ifdef DEBUG_ALLOC
            std::cerr << "Warning: Deletion of file" << filename_ << "failed!" << std::endl;
          #endif
        }
      
      }
      
      
      // accessors
      
      const String& getFilename() const
      {
        return filename_;
      }
      void setFilename(const String& filename)
      {
        filename_ = filename;
      }
      
      const off64_t& getFilesize() const
      {
        return filesize_;
      }      
      void setFilesize(const off64_t& filesize)
      {
        filesize_ = filesize;
      }      

      //nextfree_
      const off64_t& getNextfree() const
      {
        return nextfree_;
      }      
      void setNextfree(const off64_t& nextfree)
      {
        nextfree_ = nextfree;
      }       

      // totalmappingsize_      
      const off64_t& getTotalmappingsize() const
      {
        return totalmappingsize_;
      }      
      void setTotalmappingsize(const off64_t& totalmappingsize)
      {
        totalmappingsize_ = totalmappingsize;
      }  
            
      #ifdef OPENMS_WINDOWSPLATFORM
	      const HANDLE& getMmapHandle() const
	      {
	        return mmap_handle_;
	      }      
      #else
	      const int& getMmapHandle() const
	      {
	        return mmap_handle_;
	      }      
      #endif

            
      
  }; //end class

    
}
#endif //OPENMS_SYSTEM_EXTERNALALLOCATORUNIQUE_H 
