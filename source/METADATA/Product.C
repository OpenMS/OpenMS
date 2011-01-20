// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/Product.h>

using namespace std;

namespace OpenMS
{

	Product::Product()
		: CVTermList(),
			mz_(0.0),
			window_low_(0.0),
			window_up_(0.0)
	{
	}
	
	Product::Product(const Product& source):
		CVTermList(source),
	  mz_(source.mz_),
	  window_low_(source.window_low_),
	  window_up_(source.window_up_)
	{
	}
	
	Product::~Product()
	{
	}
	
	Product& Product::operator = (const Product& source)
	{
	  if (&source == this) return *this;
	  
	  CVTermList::operator=(source);
	  mz_ = source.mz_;
	  window_low_ = source.window_low_;
	  window_up_ = source.window_up_;
		
	  return *this;
	}

  bool Product::operator== (const Product& rhs) const
  {
  	return 
	    mz_ == rhs.mz_ &&
	    window_low_ == rhs.window_low_ &&
	    window_up_ == rhs.window_up_ &&
			CVTermList::operator==(rhs)
 		;
  }
  
  bool Product::operator!= (const Product& rhs) const
  {
  	return !(operator==(rhs));
 	}
	
	DoubleReal Product::getMZ() const 
	{
	  return mz_; 
	}
	
	void Product::setMZ(DoubleReal mz)
	{
	  mz_ = mz; 
	}
	
	DoubleReal Product::getIsolationWindowLowerOffset() const
	{
		return window_low_;
	}
	
	void Product::setIsolationWindowLowerOffset(DoubleReal bound)
	{
		window_low_ = bound;
	}
	
	DoubleReal Product::getIsolationWindowUpperOffset() const
	{
		return window_up_;
	}

	void Product::setIsolationWindowUpperOffset(DoubleReal bound)
	{
		window_up_ = bound;
	}

}



