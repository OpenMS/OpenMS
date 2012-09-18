// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
 
