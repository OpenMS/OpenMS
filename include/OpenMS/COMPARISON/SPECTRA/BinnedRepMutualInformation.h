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

#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDREPMUTUALINFORMATION_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDREPMUTUALINFORMATION_H

#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

namespace OpenMS
{

  /**
  BinnedRepMutualInformation calculates the Mutual Information inherent in two spectra in bin representation <br>
  
  \param intervals peak intensity is discretized into intervals
  */
  class BinnedRepMutualInformation
    : public CompareFunctor
  {
  public:

    /** @brief standard constructor <br> */
    BinnedRepMutualInformation();

    /** @brief copy constructor <br> */
    BinnedRepMutualInformation(const BinnedRepMutualInformation& source);

    /** @brief destructor <br> */
    ~BinnedRepMutualInformation();

    /** @brief assignment operator */
    BinnedRepMutualInformation& operator=(const BinnedRepMutualInformation& source);

    static FactoryProduct* create() {return new BinnedRepMutualInformation();}

    double operator()(const ClusterSpectrum& csa,const ClusterSpectrum& csb)const;

		static const String getName()
		{
			return "BinnedRepMutualInformation";
		}

    String info() const ;

  private:
    static const String info_;
  };


}

#endif //OPENMS_COMPARISON_SPECTRA_BINNEDREPMUTUALINFORMATION_H

