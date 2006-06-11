// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: FeaFiModule.h,v 1.11 2006/04/25 12:41:44 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H

#include <OpenMS/CONCEPT/FactoryProduct.h>

namespace OpenMS
{
  class FeaFiTraits;

  /** @brief Class to hold a module of the FeatureFinder algorithm 
      module accesses datastructures using BaseFeaFiTraits.
      
      @ingroup FeatureFinder
    */
  class FeaFiModule 
    : public FactoryProduct
  {

  public:
    /** @brief Inner Classes for Exception handling
      
    	 NoSuccessor-Excpetion if getNext*** or getPrev***-Methods are called on an index 
			that has no successor or predecessor 
		*/
   class NoSuccessor
     : public OpenMS::Exception::Base
     {
     public:
       NoSuccessor(const char* file, int line, const char* function, const UnsignedInt& index) throw();
       
       virtual ~NoSuccessor() throw();
       
     protected:
       UnsignedInt index_;  // index without successor/predecessor
     };

    /// standard constructor 
    FeaFiModule();

    /// copy constructor 
    FeaFiModule(const FeaFiModule& source);

    /// destructor 
    virtual ~FeaFiModule();

    /// assignment operator 
    virtual FeaFiModule& operator = (const FeaFiModule& source);

    /// set FeatureFinder traits 
    void setTraits(FeaFiTraits* traits);

  protected:
   	FeaFiTraits* traits_;
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEAFIMODULE_H
