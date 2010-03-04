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

#ifndef OPENMS_FORMAT_HANDLERS_ACQUSHANDLER_H
#define OPENMS_FORMAT_HANDLERS_ACQUSHANDLER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

namespace OpenMS
{
	namespace Internal
	{
    /**
	    @brief Read-only acqus File handler for XMass Analysis.
	    
	    acqus File contains meta data about calibration (conversion for time to mz ratio), 
	    instrument specification and acquisition method.
	    
	    @note Do not use this class directly. It is only needed for XMassFile.
    */
    class OPENMS_DLLAPI AcqusHandler
    {
      public:
        /**
			    @brief Contructor with filename.

          Open acqus File as stream and import params.
          
			    @param filename to acqus File.
			    
			    @exception Exception::FileNotFound is thrown if the file could not be opened.
			    @exception Exception::ConversionError is thrown if error conversion from String to calibration param.
		    */ 
        AcqusHandler(const String& filename);
        
        /// Destructor
        virtual ~AcqusHandler();

        /// Conversion from index to MZ ratio using internal calibration params
        DoubleReal getPosition(Size index);

        /// Read param as string      
        String getParam(const String& param);

        /// Get size of spectrum
        Size getSize();
        
      private:
        /// Private default constructor
        AcqusHandler();
        
        /// Map for params saving
        Map<String, String> params_;
        
	      /**@name Internal params for calibration */
	      //@{
        DoubleReal dw_;
        Size delay_;
        DoubleReal ml1_;
        DoubleReal ml2_;
        DoubleReal ml3_;
        Size td_;
        //@}
	  };
	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_ACQUSHANDLER_H

