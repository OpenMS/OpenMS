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
// $Maintainer: Chris Bielow  $
// --------------------------------------------------------------------------

#ifndef OPENMS_SYSTEM_EXTERNALALLOCATOR_H
#define OPENMS_SYSTEM_EXTERNALALLOCATOR_H

#include <limits>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h> // for O_RDWR etc
#include <math.h>  // for ldiv()



#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/ExternalAllocatorUnique.h>
#include <OpenMS/SYSTEM/MemoryMap.h>
#include <boost/shared_ptr.hpp>


namespace OpenMS
{

  /**
    @brief External allocator used in MSExperiment's std::vector to handle virtual memory, mapped to a swap file
		
		@note (Linux) Very slow performance has been observed on Reiser filesystems, when dealing with sparse files.
		@todo evaluate possibility of simply wrapping std:allocator when RAM is sufficient, else use external (Chris)
  */
  template <class T>
  class ExternalAllocator {
  
    protected:

      /// stores the allocator's data and prevent data corruption when copying this allocator
      boost::shared_ptr<ExternalAllocatorUnique> shared_extalloc_;
    
    public:
			
			/// allow other template instances to access private members
			template < typename T_ > friend class ExternalAllocator;
		
      // type definitions
      typedef T        value_type;
      typedef T*       pointer;
      typedef const T* const_pointer;
      typedef T&       reference;
      typedef const T& const_reference;
      typedef std::size_t    size_type;
      typedef std::ptrdiff_t difference_type;

      /// rebind allocator to type U
      template <class U>
      struct rebind 
      {
          typedef ExternalAllocator<U> other;
      };

      /// return address of @p value
      pointer address (reference value) const 
      {
				std::cout << "address of " << value << " is " << &value << "\n";
          return &value;
      }
			/// return address of @p value as const pointer
      const_pointer address (const_reference value) const 
      {
				std::cout << "address of " << value << " is " << &value << " (c)\n";
          return &value;
      }

      /* constructors and destructor
       */
    
      /// C'tor where @p filename specifies the swap file of size @p filesize bytes
      ExternalAllocator(const String& filename = File::getUniqueName(), const Offset64Int &filesize = 1) 
      {
        #ifdef DEBUG_ALLOC      
        std::cout << "<<-->> 2-tuple Ctor called \n";
        #endif          

        // assign our only member
        ExternalAllocatorUnique * ea = new ExternalAllocatorUnique(filename, filesize);
        shared_extalloc_.reset( ea );
          
      }
      
      /// copy C'tor
      ExternalAllocator(const ExternalAllocator& rhs)  
      :
        shared_extalloc_(rhs.shared_extalloc_)
      {
        #ifdef DEBUG_ALLOC      
        std::cerr << "<<-->> Copy Ctor called with nextfree_ " << shared_extalloc_->getNextfree() << "\n";
        #endif          
      }
      
			/// copy C'tor with other template parameter
      template <class U>
      ExternalAllocator (const ExternalAllocator<U>& rhs) 
			:
        shared_extalloc_(rhs.shared_extalloc_)
      {
        #ifdef DEBUG_ALLOC      
        std::cerr << "<<-->> !!strange!! Ctor called \n";
        #endif          
      }
      
      /// D'tor
      ~ExternalAllocator()  
      {
         // we should be good. shared_ptr should call dtor of data-object and delete tmp-file      
      }

      /// return maximum number of elements that can be allocated
      size_type max_size () const  {
          return shared_extalloc_->getFilesize() / sizeof(T);
      }

      
      /// allocate but don't initialize @p num elements of type T
      pointer allocate (size_type num, const void* = 0) 
      {
          pointer map = 0;
          
          // when the default c'tor for a std::vector is called, it will call allocate() for 0 elements
          if (num==0)
          {
            return map;
          }
          
          size_type alloc_bytes = num*sizeof(T);
          // round up to next free page (file location needs to be page-aligned)
          ldiv_t ldiff =  ldiv(  alloc_bytes,  MemoryMap::OpenMS_getFileBlocksize() );

          #ifdef DEBUG_ALLOC
          std::cout << "\n\n in:" << alloc_bytes << "/" << MemoryMap::OpenMS_getFileBlocksize() << "\n out:"<< ldiff.rem << " & " << ldiff.quot << "\n\n";
          #endif          

                              
          size_type block_bytes = (ldiff.rem > 0 ? ldiff.quot + 1 : ldiff.quot) * MemoryMap::OpenMS_getFileBlocksize();
          
          
          #ifdef DEBUG_ALLOC
          std::cout << "file space required: " << block_bytes << " blocksize: " << MemoryMap::OpenMS_getFileBlocksize() << std::endl;
          #endif    
      
               
          /* Now the file is ready to be mmapped */

          
          // print message and allocate memory
          #ifdef DEBUG_ALLOC
          std::cerr << "allocate " << num << " element(s)"
                    << " of size " << sizeof(T) << "(( " << shared_extalloc_->getMmapHandle() << " & " << shared_extalloc_->getNextfree() << "))"
                    << std::endl;
          #endif          

					// check if swap file is big enough
					if (!shared_extalloc_->hasFreeSwap(block_bytes)) 
					{
						//std::cerr << "extending swap in ExternalAllocator from " << shared_extalloc_->getFilesize() << " by " << block_bytes << std::endl;
						// extend swap file
						if (File::extendSparseFile(shared_extalloc_->getMmapHandle(), shared_extalloc_->getFilesize()+block_bytes))
						{
							shared_extalloc_->advanceFilesize(block_bytes);
						}
						else
						{
							std::cerr << "ExternalAllocator: Extending the swap file failed. Maybe you ran out of disk space!" << std::endl;
							return 0;
						}
						
						
					}
					
          map = static_cast<pointer> (MemoryMap::OpenMS_mmap(alloc_bytes, 
                                                             shared_extalloc_->getMmapHandle(), 
                                                             shared_extalloc_->getNextfree()
                                                            )
                                     );
          
          
                                   
          if (map == MAP_FAILED) {
            std::cerr << "MAPPING FAILED:  \n"
											<< " blocksize " << block_bytes << "\n"
                      << " nextfree: " << shared_extalloc_->getNextfree() << "( of allowed " << (shared_extalloc_->getFilesize() - block_bytes) << ")\n"
                      << " totally mapped: " << shared_extalloc_->getTotalmappingsize() << std::endl;
            #ifndef OPENMS_64BIT_ARCHITECTURE
						
						std::cerr << "The most common cause on 32-bit systems (like this one)"
											<< " is lack of virtual address space, which is usually 2-3 GB large. See the 'totally mapped' for information about your system."							
											<< "\nUpdate to a 64-bit OS to circumvent this restriction or use smaller datasets."
											<< std::endl;
											//<< "Alternatively, you forgot to use large file pointers by configuring OpenMS with\n\n"
											//<< "'./configure  --with-cxxflags=\"-D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE\"'\n\n"
											//<< std::endl;
						#endif
						return 0;
          }
          
          #ifdef DEBUG_ALLOC
          std::cerr << " allocated at MEM: " << (void*)map << " <-> filepos: " << shared_extalloc_->getNextfree() << std::endl;
          #endif          
          
          
          // set the file offset where the next mapping will start
          shared_extalloc_->advanceNextfree(block_bytes);
          
          shared_extalloc_->setTotalmappingsize(shared_extalloc_->getTotalmappingsize() + block_bytes);
          
          
          #ifdef DEBUG_ALLOC
          FILE * pFile;
          String file = shared_extalloc_->getFilename() + ".log";
          pFile = fopen (file.c_str(),"a");
          if (pFile!=NULL)
          {
            String s((long long unsigned int)shared_extalloc_->getTotalmappingsize());
            String s2((long long unsigned int)shared_extalloc_->getNextfree());
            s = "totally mapped: " + s + " offset: " + s2 + "\n";
            //std::cerr << "@@@ printing:" << s << std::endl;
            fputs (s.c_str(),pFile);
            fclose (pFile);
          }
          
          std::cerr << " new filepos: " << shared_extalloc_->getNextfree() << std::endl;
          #endif          
          
          return map;
      } // end of allocate

      /// initialize elements of allocated storage @p p with value @p value
      void construct (pointer p, const T& value) {
          // initialize memory with placement new
          //this function is essential because the vector will call this function 
          //whenever it needs to assign a value on an element (e.g. "push_back" or preinitialized constructor)
          new((void*)p)T(value);
      }

      /// destroy elements of initialized storage @p p
      void destroy (pointer /*p*/) {
          //TODO: is that really necessary?! maybe do two versions.. one fast&dangerous, one slow&secure (benchmark!)
          //p->~T();
      }

      /// deallocate storage @p p of deleted elements
      void deallocate (pointer p, size_type num) {
          #ifdef DEBUG_ALLOC
          std::cerr << "deallocate " << num << " element(s)"
                    << " of size " << sizeof(T)
                    << " at: " << (void*)p << std::endl;
          #endif
          
          // round up to the next page size
          size_type alloc_bytes = num*sizeof(T);
          ldiv_t ldiff =  ldiv(  alloc_bytes,  MemoryMap::OpenMS_getFileBlocksize() );
          size_type block_bytes = (ldiff.rem > 0 ? ldiff.quot + 1 : ldiff.quot) * MemoryMap::OpenMS_getFileBlocksize();
          
          
          shared_extalloc_->setTotalmappingsize(shared_extalloc_->getTotalmappingsize() - block_bytes);
          
          //TODO: look at int madvise (void *addr, size_t length, int advice), especially 
          /* advice = POSIX_MADV_DONTNEED: The region is no longer needed. The kernel may free these pages, causing any changes to the pages to be lost, as well as swapped out pages to be discarded. */
          int result = MemoryMap::OpenMS_unmap(p, alloc_bytes);
          if (result == OPENMS_MUNMAP_FAILURE)
          { // this is not fatal, but still severe!
            std::cerr << "Severe WARNING: unable to unmap memory at " << p << std::endl;
          }
      }
      
			/// number of bytes currently mapped from file into virtual memory
  		Offset64Int getMappingSize()
			{
				return shared_extalloc_->getTotalmappingsize();
			}
      
  }; //end class

  /// return that all specializations of this allocator are NOT interchangeable
  template <class T1, class T2>
  bool operator== (const ExternalAllocator<T1>&,
                    const ExternalAllocator<T2>&)  {
      return false;
  }
  template <class T1, class T2>
  bool operator!= (const ExternalAllocator<T1>&,
                    const ExternalAllocator<T2>&)  {
      return true;
  }

  
    
}
#endif //OPENMS_SYSTEM_EXTERNALALLOCATOR_H 
