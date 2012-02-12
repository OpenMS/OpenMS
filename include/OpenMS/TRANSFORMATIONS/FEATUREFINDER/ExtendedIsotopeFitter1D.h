// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_EXTENDEDISOTOPEFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_EXTENDEDISOTOPEFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MaxLikeliFitter1D.h>

namespace OpenMS
{
	/**
	@brief Extended isotope distribution fitter (1-dim.) approximated using linear interpolation.

	@htmlinclude OpenMS_ExtendedIsotopeFitter1D.parameters
	*/
	class OPENMS_DLLAPI ExtendedIsotopeFitter1D
		: public MaxLikeliFitter1D
	{
	 public:

		/// Default constructor
		ExtendedIsotopeFitter1D();

		/// copy constructor
		ExtendedIsotopeFitter1D(const ExtendedIsotopeFitter1D& source);

		/// destructor
		virtual ~ExtendedIsotopeFitter1D();

		/// assignment operator
		virtual ExtendedIsotopeFitter1D& operator = (const ExtendedIsotopeFitter1D& source);

		/// create new ExtendedIsotopeFitter1D object (function needed by Factory)
		static Fitter1D* create()
		{
			return new ExtendedIsotopeFitter1D();
		}

		/// name of the model (needed by Factory)
		static const String getProductName()
		{
			return "ExtendedIsotopeFitter1D";
		}

		/// return interpolation model
		QualityType fit1d(const RawDataArrayType& range, InterpolationModel*& model);

	 protected:

		/// isotope charge
		CoordinateType charge_;
		/// standard derivation in isotope
		CoordinateType isotope_stdev_;
		/// monoistopic mass
		CoordinateType monoisotopic_mz_;
		/// maximum isotopic rank to be considered
		Int max_isotope_;

		void updateMembers_();
	};
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_EXTENDEDISOTOPEFITTER1D_H
