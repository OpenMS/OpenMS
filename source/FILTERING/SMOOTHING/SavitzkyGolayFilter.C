// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

namespace OpenMS
{
  SavitzkyGolayFilter::SavitzkyGolayFilter()
  	: ProgressLogger(),
      DefaultParamHandler("SavitzkyGolayFilter"),
      coeffs_()
  {
    defaults_.setValue("frame_length",11,"The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
    defaults_.setValue("polynomial_order",4,"Order or the polynomial that is fitted.");
    
    defaultsToParam_();
  }


   SavitzkyGolayFilter::~SavitzkyGolayFilter()
  {
  }

	void SavitzkyGolayFilter::updateMembers_()
	{
    frame_size_ = (UInt)param_.getValue("frame_length"); 
    order_ = (UInt)param_.getValue("polynomial_order");
    
    //recalculate coefficients
  	if (!Math::isOdd(frame_size_))
    {
			frame_size_ += 1;
    }
    coeffs_.resize(frame_size_*(frame_size_/2+1));
		
    if (frame_size_ <= order_)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The degree of the polynomial has to be less than the frame length.", String(order_));
    }
      	
    int nr,m =frame_size_/2;
    for (int nl=0; nl<=m; ++nl)
    {
      nr=frame_size_-1-nl;

      int i,j;
      double help;

      gsl_vector* sv = gsl_vector_alloc((int)order_+1);
      gsl_vector* work =gsl_vector_alloc((int)order_+1);
      gsl_matrix* A = gsl_matrix_calloc(frame_size_,(int)order_+1);
      gsl_matrix* V = gsl_matrix_calloc((int)order_+1,(int)order_+1);


      // compute a vandermonde matrix whose columns are powers of the vector [-nL,...,nR]
      for (i=-nl; i<=nr; i++)
      {
        for (j=0; j<=(int)order_; j++)
        {
          gsl_matrix_set(A,i+nl,j,gsl_pow_int (i,j));
        }
      }

      // compute the singular-value decomposition of A
      if (gsl_linalg_SV_decomp(A,V,sv,work)!=1)
      {
        // compute B=V*inv(D)
        for (i=0; i <= (int)order_; ++i)
        {
          gsl_vector_set(sv,i,(gsl_matrix_get(V,0,i)/gsl_vector_get(sv,i)));
        }

        // compute B*transpose(U)*b, where b is the unit vector b=[1 0 ... 0]
        for (i=0; i<(int)frame_size_; ++i)
        {
          help=0;
          for (j=0;j<=(int)order_ ; ++j)
          {
            help+=gsl_vector_get(sv,j)*gsl_matrix_get(A,i,j);
          }
          coeffs_[(nl+1)*frame_size_-i-1]=help;
        }
      }
      gsl_vector_free(sv);
      gsl_vector_free(work);
      gsl_matrix_free(A);
      gsl_matrix_free(V);
    }
  }

}
