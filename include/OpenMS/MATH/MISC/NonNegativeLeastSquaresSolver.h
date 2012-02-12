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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_MISC_NONNEGATIVELEASTSQUARESSOLVER_H
#define OPENMS_MATH_MISC_NONNEGATIVELEASTSQUARESSOLVER_H

#include <OpenMS/DATASTRUCTURES/Matrix.h>

namespace OpenMS 
{
	/**
		@brief Wrapper for a non-negative least squares (NNLS) solver.
		
		It solves Ax=b, where x>0 in the least squares sense (i.e. minimum residual)
	*/
	class OPENMS_DLLAPI NonNegativeLeastSquaresSolver
	{
		public:
		
			enum RETURN_STATUS
			{ 
				SOLVED,
				ITERATION_EXCEEDED
			};

			/**
				@brief This is a wrapper for the external nnls library for the non-negative least square problem Ax=b, where x>0
				
				@param A Input matrix A of size mxn
				@param b Input vector (matrix with one column) b of size mx1
				@param x Output vector with non-negative least square solution of size mx1
				@return status of solution (either NonNegativeLeastSquaresSolver::SOLVED, NonNegativeLeastSquaresSolver::ITERATION_EXCEEDED)
				
				@throws Exception::InvalidParameters if Matrix dimensions do not fit
			*/
			static Int solve(const Matrix<double>& A, const Matrix<double>& b, Matrix<double>& x);
	};
	
} // namespace OpenMS

#endif // OPENMS_MATH_MISC_NONNEGATIVELEASTSQUARESSOLVER_H
 
