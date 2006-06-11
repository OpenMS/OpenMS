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
// $Id: FactoryProduct.h,v 1.3 2006/04/07 08:26:29 j-joachim Exp $
// $Author: j-joachim $
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FACTORYPRODUCT_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FACTORYPRODUCT_H

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /** @brief Base class for all classes T whose objects need to be constructed by Factory<T>

      Every derived class T has to implement the static function registerChildren 
      that registers all from T derived classes S at Factory<T>.<br>

      Every from T derived class S has to implement the static function T* create()
      which is going to be registered at Factory<T> 
      as well as the static function const String getName() which returns the name
      the class is registered by.
			The used name has to be assigned to @p name_ in the constructor of class S.
			
			@ingroup Concept
	*/

  class FactoryProduct
  {

  public:
    /// standard constructor
    FactoryProduct();

    /// copy constructor 
    FactoryProduct(const FactoryProduct& source);

    /// destructor
    virtual ~FactoryProduct();

    /// assignment operator
    virtual FactoryProduct& operator = (const FactoryProduct& source);

    /// parameterized constructor
    FactoryProduct(const Param& p);

    /// set parameters
    virtual void setParam(const Param& p);

    /// get parameters (const access)
    virtual const Param& getParam() const;

    /// get parameters
    virtual Param& getParam();

		/// returns the name of this module
    const String& getName() const;

    bool operator == (const FactoryProduct& rhs) const;
		bool operator != (const FactoryProduct& rhs) const;

  protected:
    mutable Param param_;
		Param defaults_;
		bool check_defaults_;
    String name_;

  };

	/// Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const FactoryProduct& prod);

}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FACTORYPRODUCT_H
