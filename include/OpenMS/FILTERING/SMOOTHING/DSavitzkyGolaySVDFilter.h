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
// $Maintainer: Eva Lange  $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_SMOOTHING_DSAVITZKYGOLAYSVDFILTER_H
#define OPENMS_FILTERING_SMOOTHING_DSAVITZKYGOLAYSVDFILTER_H

#include <OpenMS/FILTERING/SMOOTHING/DSmoothFilter.h>

#include <OpenMS/FORMAT/Param.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_pow_int.h>

namespace OpenMS
{
  using namespace Math;

  /**
    @brief Computes the Savitzky-Golay filter coefficients using singular value decomposition

    @note This filter works only for uniform raw data!

    This class represents a Savitzky-Golay lowpass-filter. The idea of the Savitzky-Golay filter
    is to find filtercoefficients that preserve higher moments, which means to approximate the underlying
    function within the moving window by a polynomial of higher order (typically quadratic or quartic).
    Therefore we least-squares fit for each data point a polynomial to all points \f$ f_i \f$ in the window
    and set \f$g_i\f$ to be the value of that polynomial at position \f$ i \f$. This method is superior
    to adjacent averaging because it tends to preserve features of the data such as peak height and width, which
    are usually 'washed out' by adjacent averaging.

    Because of the linearity of the problem, we can reduce the work computing by fitting in advance, for fictious data
    consisiting of all zeros except for a singe 1 and then do the fits on the real data just by taking linear
    combinations. There are a particular sets of filter coefficients \f$ c_n \f$ which accomplish the process of
    polynomial least-squares fit inside a moving window. To get the symmetric coefficient-matrix
    \f$C \in R^{frameSize \times frameSize}\f$ with
    \f[ C= \left[ \begin{array}{cccc} c_{0,0} & c_{0,1} & \cdots & c_{0,frameSize-1} \\
    \vdots  &         &         & \vdots            \\
    c_{frameSize-1,0} & c_{frameSize-1,2} & \ldots  & c_{frameSize-1,frameSize-1} \end{array} \right]\f]
    The first (last) \f$ \frac{frameSize}{2} \f$ rows of \f$ C \f$ we need to smooth the first (last)
    \f$ frameSize \f$ data points of the signal. So we use for the smoothing of the first data point the data point
    itself and the next \f$ frameSize-1 \f$ future points. For the second point we take the first datapoint,
    the data point itself and \f$ frameSize-2 \f$ of rightward data points... .
    We compute the Matrix \f$ C \f$ by solving the underlying least-squares problems with the singular value decomposition.
    Here we demonstrate the computation of the first row of a coefficient-matrix \f$ C \f$ for a Savitzky-Golay Filter
    of order=3 and frameSize=5:
    The design-matrix for the least-squares fit of a linear combination of 3 basis functions to 5 data points is:
    \f[ A=\left[ \begin{array}{ccc} x_0^0 & x_0^1 & x_0^2 \\ \\
    x_1^0 & x_1^1 & x_1^2 \\ \\
    x_2^0 & x_2^1 & x_2^2 \\ \\
    x_3^0 & x_3^1 & x_3^2 \\ \\
    x_4^0 & x_4^1 & x_4^2 \end{array} \right]. \f]
    To smooth the first data point we have to create a design-matrix with \f$ x=[0,\ldots, frameSize-1] \f$.
    Now we have to solve the over-determined set of \f$ frameSize \f$ linear equations
    \f[ Ac=b \f]
    where \f$ b=[1,0,\ldots,0] \f$ represents the fictious data.
    Therefore we decompose \f$ A \f$ in the orthogonal matrices \f$U, V\f$ and the diagonal matrix \f$ D \f$, which includes
    the singular values of \f$ A \f$
    \f[ A=UDV^T. \f]
    Now, it is possible to solve
    \f[ c=A^{-1}b=(VD^{-1}U^T)b.\f]
    The vector \f$c \f$ represents the wanted coefficients.

   @note This filter works only for uniform raw data!
         A polynom order of 4 is recommended.
         The bigger the frame size the smoother the signal (the more detail information get lost!). The frame size corresponds to the number
          of filter coefficients, so the width of the smoothing intervall is given by frame_size*spacing of the raw data.

  */

  template <Size D, typename MapType = MSExperiment< DRawDataPoint<1> > >
  class DSavitzkyGolaySVDFilter : public DSmoothFilter<D, MapType >
  {
  public:
    /** @name Type definitions
     */
    //@{
    ///
    typedef DSmoothFilter<D, MapType> BaseClass;
    typedef typename BaseClass::DimensionDescription DimensionDescription;
    typedef typename BaseClass::RawDataConstIterator RawDataConstIterator;
    typedef typename BaseClass::RawDataIterator RawDataIterator;
    typedef typename BaseClass::RawData RawData;
    ///
    using BaseClass::raw_filtered_;
    using BaseClass::mz_dim_;
    using BaseClass::rt_dim_;
    using BaseClass::coeffs_;
    ///
    //@}


    /** @name Constructors and Destructor
     */
    //@{
    inline DSavitzkyGolaySVDFilter() : DSmoothFilter<D, MapType>()
    {
      frame_size_ = 17;
      order_=4;

      coeffs_.resize(frame_size_*(frame_size_/2+1));
      computeCoeffs_();

      raw_filtered_=0;
    }

    /** Class constructor setting the coefficients in the filter.
     */
    /// Class constructor setting the coefficients in the filter. Please note that the frame size must be odd and that the order of the polynomial is more less than the frame size.
    DSavitzkyGolaySVDFilter(const Param& parameters) throw (Exception::InvalidValue)
    {
      param_ = parameters;

      if (D == 1)
      {
        rt_dim_ = -1;
        mz_dim_ = 0;
      }
      else
        if (D == 2)
        {
          rt_dim_ = DimensionDescription::RT;
          mz_dim_ = DimensionDescription::MZ;
        }


      // if a smoothing parameter is missed in the param object the value should be substituted by a dv value
      DataValue dv = param_.getValue("FrameLength");
      if (dv.isEmpty() || dv.toString() == "") frame_size_ = 17;
      else frame_size_ = (int)dv;

      if (!isOdd(frame_size_))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The frame_size must be an odd integer!", String(frame_size_));
      }


      dv = param_.getValue("PolynomOrder");
      if (dv.isEmpty() || dv.toString() == "") order_ = 4;
      else order_ = (unsigned int)dv;

      if (frame_size_ <= order_)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The degree of the polynomial has to  must be less than the frame length.", String(order_));
      }


      coeffs_.resize(frame_size_*(frame_size_/2+1));
      computeCoeffs_();
    }

    ///
    inline DSavitzkyGolaySVDFilter(const DSavitzkyGolaySVDFilter& s)
        : DSmoothFilter<D, MapType>(s),
        param_(s.param_),
        frame_size_(s.frame_size_),
        order_(s.order_)
    {   }
    ///
    virtual ~DSavitzkyGolaySVDFilter()
    { }
    //@}

    /** @name Assignment
     */
    //@{
    ///
    inline DSavitzkyGolaySVDFilter& operator=(const DSavitzkyGolaySVDFilter& s)
    {
      param_ = s.param_;

      order_=s.order_;
      frame_size_=s.frame_size_;
      raw_filtered_=s.raw_filtered_;
      mz_dim_=s.mz_dim_;
      rt_dim_=s.rt_dim_;
      coeffs_=s.coeffs;

      return *this;
    }
    //@}

    /** Accessors
     */
    //@{
    /// Non-mutable access to the order
    inline const unsigned int& getOrder() const { return order_; }
    /// Mutable access to the order
    inline unsigned int& getOrder() { return order_; }
    /// Mutable access to the order
    inline void setOrder(const unsigned int order)
    {
      order_=order;
      computeCoeffs_();
    }

    /// Non-mutable access to length of the smoothing window
    inline const unsigned int& getWindowSize() const { return frame_size_; }
    /// Mutable access to the length of the window
    inline void setWindowSize(const unsigned int frame_size)
    {
      frame_size_=frame_size;
      coeffs_.resize(frame_size_*(frame_size_/2+1));
      computeCoeffs_();
    }
    /// Non-mutable access to the parameter object
    inline const Param& getParam() const { return param_; }
    /// Mutable access to the parameter object
    inline void setParam(const Param& param) throw (Exception::InvalidValue)
    {
      param_ = param;

      // read the new parameter
      DataValue dv = param_.getValue("FrameLength");
      if (!(dv.isEmpty() || dv.toString() == "")) frame_size_ = (int)dv;

      if (!isOdd(frame_size_))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The frame_size must be an odd integer!", String(frame_size_));
      }

      dv = param_.getValue("PolynomOrder");
      if (!(dv.isEmpty() || dv.toString() == "")) order_ = (unsigned int)dv;

      if (frame_size_ <= order_)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The degree of the polynomial has to  must be less than the frame length.", String(order_));
      }
      coeffs_.resize(frame_size_*(frame_size_/2+1));
      computeCoeffs_();
    }
    //@}
  protected:
    /// Parameter object
    Param param_;
    /// Size of the filter kernel (number of pre-tabulated coefficients)
    unsigned int frame_size_;
    /// The order of the smoothing polynomial.
    unsigned int order_;


    void convolute_(RawDataConstIterator it_begin, RawDataConstIterator it_end, RawDataIterator new_raw_first) throw (Exception::InvalidSize)
    {
      if (distance(it_begin,it_end) <= (int)frame_size_)
      {
        return; // scan does not contain enough points, so we skip it
      }

      int i;
      unsigned int j;
      int mid=(frame_size_/2);
      double help;

      RawDataConstIterator it_forward;
      RawDataConstIterator it_help;

      // compute the transient on
      for (i=0; i <= mid; ++i)
      {
        it_forward=(it_begin-i);
        help=0;

        for (j=0; j < frame_size_; ++j)
        {
          help+=it_forward->getIntensity()*coeffs_[(i+1)*frame_size_-1-j];
          ++it_forward;
        }


        new_raw_first->setPosition(it_begin->getPosition());
        new_raw_first->setIntensity(help);
        ++new_raw_first;
        ++it_begin;
      }

      // compute the steady state output
      it_help=(it_end-mid);
      while (it_begin!=it_help)
      {
        it_forward=(it_begin-mid);
        help=0;

        for (j=0; j < frame_size_; ++j)
        {
          help+=it_forward->getIntensity()*coeffs_[mid*frame_size_+j];
          ++it_forward;
        }


        new_raw_first->setPosition(it_begin->getPosition());
        new_raw_first->setIntensity(help);
        ++new_raw_first;
        ++it_begin;
      }

      // compute the transient off
      for (i=(mid-1); i >= 0; --i)
      {
        it_forward=(it_begin-(frame_size_-i-1));
        help=0;

        for (j=0; j < frame_size_; ++j)
        {
          help+=it_forward->getIntensity()*coeffs_[i*frame_size_+j];
          ++it_forward;
        }


        new_raw_first->setPosition(it_begin->getPosition());
        new_raw_first->setIntensity(help);
        ++new_raw_first;
        ++it_begin;
      }
    }


    /// Compute the coefficient-matrix \f$ C \f$ of the filter.
    void computeCoeffs_()
    {
      int nr,m =frame_size_/2;

      for (int nl=0; nl<=m; ++nl)
      {
        nr=frame_size_-1-nl;

        int i,j;
        double help;

        gsl_vector* sv = gsl_vector_alloc((int)order_+1);
        gsl_vector* work =gsl_vector_alloc((int)order_+1);
        gsl_vector* b =gsl_vector_alloc(frame_size_);
        gsl_vector_set_all(b,1.0);

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
        gsl_vector_free(b);

        gsl_matrix_free(A);
        gsl_matrix_free(V);
      }
    }
  };


} // namespace OpenMS
#endif

