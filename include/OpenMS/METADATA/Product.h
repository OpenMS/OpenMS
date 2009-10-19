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

#ifndef OPENMS_METADATA_PRODUCT_H
#define OPENMS_METADATA_PRODUCT_H

#include <OpenMS/METADATA/CVTermList.h>

namespace OpenMS 
{
	/**
		@brief Product meta information.
		
		This class describes the product isolation window for special scan types, such as MRM.
  
		@ingroup Metadata
	*/  
  class OPENMS_DLLAPI Product
  	: public CVTermList
  {
    
    public:
    
      /// Constructor
      Product();
      /// Copy constructor
      Product(const Product& source);
      /// Destructor
      ~Product();
      
      /// Assignment operator
      Product& operator= (const Product& source);

      /// Equality operator
      bool operator== (const Product& rhs) const;
      /// Equality operator
      bool operator!= (const Product& rhs) const;
			
			/// returns the target m/z
      DoubleReal getMZ() const;
      /// sets the target m/z
      void setMZ(DoubleReal mz);
			
      /// returns the lower offset from the target m/z
      DoubleReal getIsolationWindowLowerOffset() const;
      /// sets the lower offset from the target m/z
      void setIsolationWindowLowerOffset(DoubleReal bound);
      
      /// returns the upper offset from the target m/z
      DoubleReal getIsolationWindowUpperOffset() const;
      /// sets the upper offset from the target m/z
      void setIsolationWindowUpperOffset(DoubleReal bound);
      
    protected:
      
      DoubleReal mz_;
      DoubleReal window_low_;
      DoubleReal window_up_;
    };
} // namespace OpenMS

#endif // OPENMS_METADATA_PRODUCT_H
