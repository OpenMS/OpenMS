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
// $Id: TradSeqQuality.h,v 1.3 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_TRADSEQQUALITY_H
#define OPENMS_COMPARISON_CLUSTERING_TRADSEQQUALITY_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

namespace OpenMS{
  /**
  TradSeqQuality returns a number > 0 if the sequest score are above a certain XCorr and above a certain deltaCN<br>
  
  \param xcorr_1+ min XCorr for charge state 1
  \param xcorr_2+ min XCorr for charge state 2
  \param xcorr_3+ min XCorr for charge state 3
  \param dCn_1+ min deltaCN for charge state 1
  \param dCn_2+ min deltaCN for charge state 2
  \param dCn_3+ min deltaCN for charge state 3
  */
  class TradSeqQuality : public FilterFunctor
  {
  public:
    /** @brief standard constructor <br> */
    TradSeqQuality();

    /** @brief copy constructor <br> */
    TradSeqQuality(const TradSeqQuality& source);

    /** @brief assignment operator <br> */
    TradSeqQuality& operator=(const TradSeqQuality& source );

    /** @brief destructor <br> */
    ~TradSeqQuality();

    static FactoryProduct* create() { return new TradSeqQuality();}

    std::vector<double> operator()(const ClusterSpectrum& spec);

    String info() const;

		static const String getName()
		{
			return "TradSeqQuality";
		}

  private:
    static const String info_;
  };
}
#endif // OPENMS_COMPARISON_CLUSTERING_TRADSEQQUALITY_H
