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
#ifndef OPENMS_COMPARISON_SPECTRA_COMPAREFUNCTOR_H
#define OPENMS_COMPARISON_SPECTRA_COMPAREFUNCTOR_H

#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

namespace OpenMS
{

  /**
  CompareFunctor classes return a value for a pair of ClusterSpectrum objects<br>
  ideally the value should reflect the similarity of the pair<br>
  similarities of spectra should be > 0<br>
  \param filterwindow
    maximum mass difference for spectra that get similarity > 0
  */
  class CompareFunctor : public FactoryProduct
  {
  public:

    /** @brief standard constructor <br> */
    CompareFunctor();

    /** @brief copy constructor <br> */
    CompareFunctor(const CompareFunctor& source);

    /** @brief destructor <br> */
    virtual ~CompareFunctor();

    /** @brief assignment operator <br> */
    CompareFunctor& operator = ( const CompareFunctor& source );

		static void registerChildren();

    /** @brief function call operator, calculates the similarity <br> */
    virtual double operator()(const ClusterSpectrum&, const ClusterSpectrum&) const = 0;

    /** @brief function call operator, calculates the self similarity <br> */
    virtual double operator()(const ClusterSpectrum& a) const {return (*this)(a,a);}

    /** @brief preliminary check for similarity <br>*/
    virtual double filter(const ClusterSpectrum&,const ClusterSpectrum&) const;

    /** @brief returns type of compared spectrum representation <br>*/
    bool usebins() const {return usebins_;}

  protected:

    /**
    type of spectrum representation the CompareFunctor operates on <br>
      usebin_ = 1 => bin representation <br>
      usebin_ = 0 => stick spectrum representation<br>
    */
    bool usebins_;

  };

}
#endif // OPENMS_COMPARISON_SPECTRA_COMPAREFUNCTOR_H
