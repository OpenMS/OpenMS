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
// $Id: KellerQuality.h,v 1.3 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_KELLERQUALITY_H
#define OPENMS_COMPARISON_CLUSTERING_KELLERQUALITY_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <string>

namespace OpenMS{
  /**
  KellerQuality assigns a quality measure to a spectrum based on the linear regression formula from their 2003 Paper<br>
  http://www.systemsbiology.org/PDFs/Keller.Emperical%20statistical%20model.Anal%20Chem.02.pdf<br>
  */
  class KellerQuality : public FilterFunctor
  {
  public:
    /** @brief standard constructor <br> */
    KellerQuality();

    /** @brief copy constructor <br> */
    KellerQuality(const KellerQuality& source);

    /** @brief assignment operator <br> */
    KellerQuality& operator=(const KellerQuality& source );

    /** @brief destructor <br> */
    ~KellerQuality();

    static FactoryProduct* create() { return new KellerQuality();}

    std::vector<double> operator()(const ClusterSpectrum& spec);

    String info() const;

		static const String getName()
		{
			return "KellerQuality";
		}

  private:
    static const String info_;
  };
}
#endif // OPENMS_COMPARISON_CLUSTERING_KELLERQUALITY_H
