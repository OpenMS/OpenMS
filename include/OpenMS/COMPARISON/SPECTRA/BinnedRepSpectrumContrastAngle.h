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
// $Id: BinnedRepSpectrumContrastAngle.h,v 1.2 2006/03/29 13:06:18 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_COMPARISON_CLUSTERING_BINNEDREPSPECTRUMCONTRASTANGLE_H
#define OPENMS_COMPARISON_CLUSTERING_BINNEDREPSPECTRUMCONTRASTANGLE_H

#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

namespace OpenMS
{

  /**
  BinnedRepSpectrumContrastAngle calculates the spectral contrast angle between two spectra in bin representation<br>
  the functor does not cover the whole function, the normalization is done by the ClusterRun to allow for more flexibility<br>
  */
  class BinnedRepSpectrumContrastAngle
    : public CompareFunctor
  {
  public:

    /** @brief standard constructor <br> */
    BinnedRepSpectrumContrastAngle();

    /** @brief copy constructor <br> */
    BinnedRepSpectrumContrastAngle(const BinnedRepSpectrumContrastAngle& source);

    /** @brief destructor <br> */
    ~BinnedRepSpectrumContrastAngle();

    /** @brief assignment operator <br> */
    BinnedRepSpectrumContrastAngle& operator=(const BinnedRepSpectrumContrastAngle& source);

    static FactoryProduct* create(){return new BinnedRepSpectrumContrastAngle();}

    double operator()(const ClusterSpectrum& a ,const ClusterSpectrum& b ) const;

		static const String getName()
		{
			return "BinnedRepSpectrumContrastAngle";
		}

    String info() const;

  private:
		
    static const String info_;
  };

}

#endif //OPENMS_COMPARISON_CLUSTERING_BINNEDREPSPECTRUMCONTRASTANGLE_H

