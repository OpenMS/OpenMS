// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    Flex series file support
// --------------------------------------------------------------------------
//  Copyright (C) 2009 -- Guillaume Belz (guillaume.belz@chu-lyon.fr)
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
// $Maintainer: Guillaume Belz$
// $Authors: Guillaume Belz$
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_FIDHANDLER_H
#define OPENMS_FORMAT_HANDLERS_FIDHANDLER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <fstream>

namespace OpenMS
{
	namespace Internal
	{
    /**
	    @brief Read-only fid File handler for XMass Analysis.
	    
	    fid File contains intensity array. Intensity for each point are coded in 4 bytes integer.
	    
	    @note Do not use this class directly. It is only needed for XMassFile.
    */
    class OPENMS_DLLAPI FidHandler
     : public std::ifstream
    {
      public:
        /**
			    @brief Contructor with filename.

          Open fid File as stream and initialize index.          
          
			    @param filename to fid File.
		    */ 
        FidHandler(const String& filename);
        
        /// Destructor
        virtual ~FidHandler();
        
		    /// Get index of current position (without position moving).
        Size getIndex();
        
        /// Get intensity of current position and move to next position.
        Size getIntensity();
        
      private:
        /// Private default constructor
        FidHandler();
        
        /// Index of position
        Size index_;
    };
	
	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_FIDHANDLER_H

