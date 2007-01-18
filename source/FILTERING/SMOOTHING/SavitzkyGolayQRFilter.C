// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayQRFilter.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

namespace OpenMS
{
	using namespace Math;
	
  SavitzkyGolayQRFilter::SavitzkyGolayQRFilter()
      : SmoothFilter()
  {
  	
  	defaults_.setValue("frame_length",17);
  	defaults_.setValue("polynomial_order",4);
    
    frame_size_=17;
    order_=4;
    
    coeffs_.clear();
    coeffs_.resize(frame_size_*(frame_size_/2+1));
    computeCoeffs_();
  }

  void SavitzkyGolayQRFilter::setParam(Param param) throw (Exception::InvalidValue)
  {
		param.setDefaults(defaults_);
    param.checkDefaults("SavitzkyGolayQRFilter",defaults_);
    
    frame_size_ = (int)param.getValue("frame_length");
    if (!isOdd(frame_size_))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The frame_size has to be an odd integer!", String(frame_size_));
    }

    order_ = (unsigned int)param.getValue("polynomial_order");

    if (frame_size_ <= order_)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The degree of the polynomial has to be less than the frame length.", String(order_));
    }

    coeffs_.clear();
    coeffs_.resize(frame_size_*(frame_size_/2+1));
    computeCoeffs_();
  }

  void SavitzkyGolayQRFilter::setWindowSize(const unsigned int& frame_size)
  {
    frame_size_=frame_size;
    
    coeffs_.clear();
    coeffs_.resize(frame_size_*(frame_size_/2+1));
    computeCoeffs_();
  }

  void SavitzkyGolayQRFilter::setOrder(const unsigned int& order)
  {
    order_=order;
    computeCoeffs_();
  }

  void SavitzkyGolayQRFilter::computeCoeffs_() throw (Exception::InvalidValue)
  {
  	if (!isOdd(frame_size_))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The frame_size has to be an odd integer!", String(frame_size_));
    }

    if (frame_size_ <= order_)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The degree of the polynomial has to be less than the frame length.", String(order_));
    }
    
    int nr, m = frame_size_/2;

    for (int nl=0; nl <= m; ++nl)
    {
      nr=frame_size_-1-nl;

      int imj, ipj, k, kk, mm;
      double fac, sum;
      gsl_vector* b =gsl_vector_calloc((int)(order_+1));
      gsl_vector* tau =gsl_vector_calloc((int)(order_+1));
      gsl_matrix* A = gsl_matrix_calloc((int)(order_+1),(int)(order_+1));

      // compute the "normal equations", which are \f$ A^TA \f$ and \f$ A \f$ is the designmatrix.
      for (ipj = 0; ipj <= 2 * (int)order_; ipj++)
      {
        sum = 0.0;
        if ( ipj == 0 )
        {
          sum = 1.0;
        }
        for (k = 1; k <= nr; k++)
        {
          sum += gsl_pow_int(k,ipj);
        }
        for (k = 1; k <= nl; k++)
        {
          sum += gsl_pow_int((-k),ipj);
        }
        mm = std::min(ipj,2*(int)order_-ipj);
        for (imj = -mm; imj <= mm; imj += 2)
        {
          gsl_matrix_set(A,(ipj+imj)/2,(ipj-imj)/2,sum);
        }
      }

      gsl_vector_set(b,0,1.0);
      gsl_linalg_QR_decomp(A,tau);

      int qr_success = gsl_linalg_QR_solve(A,tau,b,b);

      //get one row of the inverse \f$ (A^TA)^{-1} \f$ by the QR decomposition with only a single backsubstitution.
      OPENMS_PRECONDITION((qr_success != 1),"QR Decomposition of the normal equations is not possible!");

			++qr_success; // suppresses compiler warning 
			
      // we compute the coefficients
      for (k = -nl, kk = frame_size_-1; k <= m; k++, kk--)
      {
        sum = gsl_vector_get(b,0);
        fac = 1.0;
        for (mm = 1; mm <= (int)order_; mm++)
        {
          fac *= k;
          // each Savitzky-Golay coefficient is the dot product of powers of an integer with the inverse matrix row.
          sum += gsl_vector_get(b,mm) * fac;
        }
        coeffs_[nl*frame_size_+kk]=sum;
      }

      gsl_vector_free(b);
      gsl_vector_free(tau);
      gsl_matrix_free(A);
    }
  }
}
