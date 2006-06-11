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
// $Id: BinnedRepSharedPeakCount.h,v 1.2 2006/03/29 13:06:18 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDREPSHAREDPEAKCOUNT_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDREPSHAREDPEAKCOUNT_H

#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

namespace OpenMS
{

  //calculates the Shared Peak Count for 2 BinnedReps
  class BinnedRepSharedPeakCount
    : public CompareFunctor
  {
  public:
    BinnedRepSharedPeakCount();
    BinnedRepSharedPeakCount(const BinnedRepSharedPeakCount& source);
    ~BinnedRepSharedPeakCount();
    BinnedRepSharedPeakCount& operator=(const BinnedRepSharedPeakCount& source);
    static FactoryProduct* create() {return new BinnedRepSharedPeakCount();}
    double operator()(const ClusterSpectrum& csa,const ClusterSpectrum& csb)const;
		static const String getName()
		{
			return "BinnedRepSharedPeakCount";
		}
    String info() const ;
  private:
    static const String info_;
  };

}

#endif //OPENMS_COMPARISON_SPECTRA_BINNEDREPSHAREDPEAKCOUNT_H

