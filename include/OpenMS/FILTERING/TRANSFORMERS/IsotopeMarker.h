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
// $Id: IsotopeMarker.h,v 1.3 2006/04/05 11:18:23 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_ISOTOPEMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_ISOTOPEMARKER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>

namespace OpenMS
{

  /**
  IsotopeMarker marks peak pairs which could represent an ion and its isotope<br>
  
  \param marks times a peak needs to be marked to get marked in the result
  \param mz_variation m/z tolerance
  \param in_variation intensity variation in fraction of the theoretical isotope peak
  */
  class IsotopeMarker
    :public PeakMarker
  {
  public:
     /** @brief standard constructor <br> */
    IsotopeMarker();

    /** @brief copy constructor <br> */
    IsotopeMarker(const IsotopeMarker& source);

    /** @brief destructor <br> */
    ~IsotopeMarker();

    /** @brief assignment operator <br> */
    IsotopeMarker& operator=(const IsotopeMarker& source);

    static FactoryProduct* create() { return new IsotopeMarker();}
    std::map<double,bool> operator()( MSSpectrum< DPeak<1> >&) const;
    String info() const;

		static const String getName()
		{
			return "IsotopeMarker";
		}
		
  private:
    static const String info_;
  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_ISOTOPEMARKER_H
