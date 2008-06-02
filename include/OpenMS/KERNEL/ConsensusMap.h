// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_CONSENSUSMAP_H
#define OPENMS_KERNEL_CONSENSUSMAP_H

#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

namespace OpenMS
{
	/**
    @brief A container for consensus elements.
    
    A ConsensusMap is a container holding 2-dimensional consensus elements (ConsensusFeature)
    which in turn represent combined elements of 2-dimensional experiments.
    The map is implemented as a vector of elements.
    
    The map indices used in the consensus features should be registered in this class.
 		
 		@todo Add method to dump gnuplot files of consensus maps (Clemens)
 		
    @ingroup Kernel
  */
	class ConsensusMap : public DPeakArray<ConsensusFeature>
	{
	  public:
	    /// Base class type
	    typedef DPeakArray<ConsensusFeature > Base;
	  	
	    /// Default onstructor
	    inline ConsensusMap()
	      : Base()
	    {
	    }
	
	    /// Copy constructor
	    inline ConsensusMap(const ConsensusMap& source)
	      : Base(source),
	        filenames_(source.filenames_)
	    {
	    }
	
	    /// Destructor
	    inline ~ConsensusMap()
	    {
	    }
	
	    /// Creates a ConsensusMap with n elements
	    inline ConsensusMap(Base::size_type n) 
	    	: Base(n)
	    {
	    }
	
	    /// Assignment operator
	    ConsensusMap& operator = (const ConsensusMap& source)
	    {
	      if (this==&source) return *this;
	
	      Base::operator=(source);
	      filenames_ = source.filenames_;
	      
	      return *this;
	    }
	
	    /// Non-mutable access to the filenames
	    inline const Map<UInt,String>& getFileNames() const
	    {
	      return filenames_;
	    }
	    
	    /// Set a file name
	    inline void setFileName(UInt index, const String& name)
	    {
	      filenames_[index] = name;
	    }
	    
	    /// Merge overlapping consensus elements
	    void merge(ConsensusMap& new_map);
			
			///Checks if all map identifiers in FeatureHandles are have a filename associated
			bool isValid() const;
			
	  protected:
	  
	    /// Map from index to filenames 
	    Map<UInt,String> filenames_;
  };

  ///Print the contents of a ConsensusMap to a stream.
  std::ostream& operator << (std::ostream& os, const ConsensusMap& cons_map);

} // namespace OpenMS

#endif // OPENMS_KERNEL_CONSENSUSMAP_H
