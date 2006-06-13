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
// $Id: SqrtMower.h,v 1.4 2006/04/05 11:18:23 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_SQRTMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_SQRTMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>

namespace OpenMS
{
  /**
  	@brief Scales the intensity of peaks to log(intensity)
  */
  class SqrtMower : public PreprocessingFunctor
  {
  public:
    /// standard constructor
    SqrtMower();

    /// copy constructor
    SqrtMower(const SqrtMower& source);

    /// destructor
    virtual ~SqrtMower();

    /// assignment operator <br>*/
    SqrtMower& operator=(const SqrtMower& source);

    static FactoryProduct* create() { return new SqrtMower();}

    virtual void operator()(MSSpectrum< DPeak<1> >&) const;

    virtual String info() const ;

		static const String getName()
		{
			return "SqrtMower";
		}
  private:
    static const String info_;
  };

}

#endif // OPENMS_FILTERING_TRANSFORMERS_SQRTMOWER_H
