// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_SYSTEM_EXTERNALALLOCATORUNIQUE_H
#define OPENMS_SYSTEM_EXTERNALALLOCATORUNIQUE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Types.h>

#include <vector>
#include <cstdlib>

namespace OpenMS
{

  /**
    @brief Unique settings for an ExternalAllocator. 
		
		When an ExternalAllocator is copied it is necessary to ensure data consistency between the
		copied instances.
		This class is held by a boost sharedPtr within the ExternalAllocator.
    
		It contains the size, name and handle of the swap file.
  */
  class ExternalAllocatorUnique {
  
		private:
			/// do not allow default C'tor
			ExternalAllocatorUnique()
			{
			}
	
    protected:

      /// name of temporary file
      String filename_;
      /// size of temporary file
      Offset64Int filesize_;     
      /// next byte position in file where the next mapping is scheduled
      Offset64Int nextfree_;
      
      /// filehandle to swap file
      #ifdef OPENMS_WINDOWSPLATFORM
        HANDLE mmap_handle_;
      #else
        int mmap_handle_;      
      #endif
      
      /// just for informational purposes: how many bytes are mapped
      Offset64Int totalmappingsize_;
      
      
      // /// freed blocks which shall be reused before opening a new one
      //std::vector < std::pair <Offset64Int, Offset64Int > > freeblock_;

    public:  
    
    /* constructors and destructor
     */
         
      /// C'tor
      ExternalAllocatorUnique(const String &filename, const Offset64Int &filesize)
      :
        filename_(filename),
        filesize_(filesize),
        nextfree_(0),
        totalmappingsize_(0)
				//,freeblock_()
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
      
      /// copy C'tor
      ExternalAllocatorUnique(const ExternalAllocatorUnique& rhs)  
      :
        filename_(rhs.filename_),
        filesize_(rhs.filesize_),
        nextfree_(rhs.nextfree_),
        mmap_handle_(rhs.mmap_handle_),
        totalmappingsize_(rhs.totalmappingsize_)
				//, freeblock_(rhs.freeblock_)
      {
        #ifdef DEBUG_ALLOC      
        std::cerr << "--- Copy Ctor called with nextfree_ " << nextfree_ << "\n";
        #endif          
      }
      
      
      /// D'tor
      ~ExternalAllocatorUnique()  
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
      
      
			/**	@name	read-only accessors */
			//@{
			
      /// get the name of the swap file
      const String& getFilename() const
      {
        return filename_;
      }

			/// get handle to the swap file
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
			//@}
			
			/**	@name	read & write accessors */
			//@{
			/// increase the filesize by @p x bytes
      void advanceFilesize(const Offset64Int& x)
      {
        filesize_ += x;
      }       
			
			/// get the size of the swap file      
      const Offset64Int& getFilesize() const
      {
        return filesize_;
      }      
			
      /// get next free byte position of swap file
      const Offset64Int& getNextfree() const
      {
        return nextfree_;
      }      
			/// advance the next free byte position by @p x bytes
      void advanceNextfree(const Offset64Int& x)
      {
        nextfree_ += x;
      }       

      /// get current number of bytes mapped from swap file into virtual memory      
      const Offset64Int& getTotalmappingsize() const
      {
        return totalmappingsize_;
      }
			/// set new mapping size @p x
      void setTotalmappingsize(const Offset64Int& x)
      {
        totalmappingsize_ = x;
      }  
			//@}
      
			/// determine if a new mapping at the current file position would go beyond EOF
			bool hasFreeSwap(const Offset64Int& bytes_needed)
			{
				return (filesize_ > bytes_needed+nextfree_);
			}
            
      
  }; //end class

    
}
#endif //OPENMS_SYSTEM_EXTERNALALLOCATORUNIQUE_H 
