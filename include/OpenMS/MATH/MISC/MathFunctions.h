// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_MISC_MATHFUNCTIONS_H
#define OPENMS_MATH_MISC_MATHFUNCTIONS_H

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
	/**
		@brief %Math namespace.

		Contains mathemtical auxiliary functions.

		@ingroup Concept
	*/
	namespace Math
	{
		/**
			@brief rounds @p x up to the next decimal power 10 ^ @p decPow

			@verbatim
			e.g.: (123.0 , 1)  => 130
			      (123.0 , 2)  => 200
					  (0.123 ,-2)  => 0.13 ( 10^-2 = 0.01 )
			@endverbatim

			@ingroup MathFunctionsMisc
		*/
		inline static double ceilDecimal(double x, int decPow)
		{
			return (ceil(x/pow(10.0,decPow)))*pow(10.0,decPow); // decimal shift right, ceiling, decimal shift left
		}

		/**
			@brief rounds @p x to the next decimal power 10 ^ @p decPow

			@verbatim
			e.g.: (123.0 , 1)  => 120
			      (123.0 , 2)  => 100
			@endverbatim

			@ingroup MathFunctionsMisc
		*/
		inline static double roundDecimal(double x, int decPow)
		{
			if (x>0) return (floor(0.5+x/pow(10.0,decPow)))*pow(10.0,decPow);
			return -((floor(0.5+fabs(x)/pow(10.0,decPow)))*pow(10.0,decPow));
		}

		/**
			@brief transforms point @p x of interval [left1,right1] into interval [left2,right2]

			@ingroup MathFunctionsMisc
		*/
		inline static double intervalTransformation(double x,double left1,double right1,double left2,double right2)
		{
			return left2 + (x - left1) * (right2 - left2) / (right1 - left1);
		}

		/**
			@brief Transforms a number from linear to log10 scale. Avoids negative logarithms by adding 1.

			@param x The number to transform

			@ingroup MathFunctionsMisc
		*/
		inline double linear2log(double x)
		{
			return log10(x+1); //+1 to avoid negative logarithms
		}

		/**
			@brief Transforms a number from log10 to to linear scale. Subtracts the 1 added by linear2log(double)

			@param x The number to transform

			@ingroup MathFunctionsMisc
		*/
		inline double log2linear(double x)
		{
			return pow(10,x)-1;
		}

		/**
			@brief Returns true if the given interger is odd

			@ingroup MathFunctionsMisc
		*/
		inline bool isOdd(UInt x)
		{
			return ((x & 1)!=0);
		}

		/**
			@brief Rounds the value

			@ingroup MathFunctionsMisc
		*/
		template <typename T>
		T round (T x)
		{
			if (x >= T(0))
			{
				return T(floor(x+T(0.5)));
			}
			else
			{
				return T(ceil(x-T(0.5)));
			}
		}
		
		/**
			@brief Returns if @p a is approximately equal @p, allowing a tolerance of @p tol
		
			@ingroup MathFunctionsMisc
		*/
		inline static bool approximatelyEqual(DoubleReal a, DoubleReal b, DoubleReal tol)
		{
			return std::abs(a - b) <= tol;
		}

	} // namespace Math
} // namespace OpenMS

#endif // OPENMS_MATH_MISC_MATHFUNCTIONS_H
