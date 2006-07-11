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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDREPSUMAGREEINGINTENSITIES_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDREPSUMAGREEINGINTENSITIES_H

#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

namespace OpenMS
{

  /**
  BinnedRepSumAgreeingIntensities calculates the sum of agreeing intensities for two spectra in bin representation<br>
  */
  class BinnedRepSumAgreeingIntensities
    : public CompareFunctor
  {
  public:

    /** @brief standard constructor <br> */
    BinnedRepSumAgreeingIntensities();

    /** @brief copy constructor <br> */
    BinnedRepSumAgreeingIntensities(const BinnedRepSumAgreeingIntensities& source);

    /** @brief destructor <br> */
    ~BinnedRepSumAgreeingIntensities();

    /** @brief assignment operator <br> */
    BinnedRepSumAgreeingIntensities& operator=(const BinnedRepSumAgreeingIntensities& source);

    static FactoryProduct* create() { return new BinnedRepSumAgreeingIntensities();}

    double operator()(const ClusterSpectrum& csa ,const ClusterSpectrum& csb )const;

		static const String getName()
		{
			return "BinnedRepSumAgreeingIntensities";
		}

    String info() const;

  private:
    static const String info_;
  };

}

#endif //OPENMS_COMPARISON_SPECTRA_BINNEDREPSUMAGREEINGINTENSITIES_H

