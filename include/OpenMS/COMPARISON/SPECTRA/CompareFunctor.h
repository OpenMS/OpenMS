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

#include <OpenMS/CONCEPT/FactoryProduct.h>

namespace OpenMS
{

  /**
		@defgroup Comparison Comparison
	*/

	/**
		@defgroup SpectraComparison Spectra Comparison

		@ingroup Comparison
	*/

	/**
	
		@brief Base class for compare functors of spectra; compare functors returns a similiarity value of two spectra
	
  	CompareFunctor classes return a value for a pair of ClusterSpectrum objects<br>
  	ideally the value should reflect the similarity of the pair<br>
  	similarities of spectra should be > 0<br>
		
  	@param filterwindow
		
    maximum mass difference for spectra that get similarity > 0

		@ingroup SpectraComparison
  */
	class ClusterSpectrum;
	
  class CompareFunctor : public FactoryProduct
  {

  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    CompareFunctor();

    /// copy constructor
    CompareFunctor(const CompareFunctor& source);

    /// destructor
    virtual ~CompareFunctor();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    CompareFunctor& operator = (const CompareFunctor& source);

		/// @todo needed for factory ?
		/// static void registerChildren();

    /// function call operator, calculates the similarity
    virtual double operator () (const ClusterSpectrum&, const ClusterSpectrum&) const { return 0; }

    /// function call operator, calculates the self similarity
    double operator()(const ClusterSpectrum& a) const { return (*this)(a, a); }
		// @}

		// @name Accessors
		// @{
    /// preliminary check for similarity
    double filter(const ClusterSpectrum&, const ClusterSpectrum&) const;

    /// returns type of compared spectrum representation
    bool usebins() const { return usebins_; }
		// @}

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
