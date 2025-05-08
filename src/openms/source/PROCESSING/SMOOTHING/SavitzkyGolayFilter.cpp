// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <Eigen/Core>
#include <Eigen/SVD>

namespace OpenMS
{
  SavitzkyGolayFilter::SavitzkyGolayFilter() :
    ProgressLogger(),
    DefaultParamHandler("SavitzkyGolayFilter"),
    coeffs_()
  {
    defaults_.setValue("frame_length", 11, "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
    defaults_.setValue("polynomial_order", 4, "Order or the polynomial that is fitted.");

    defaultsToParam_();
  }

  SavitzkyGolayFilter::~SavitzkyGolayFilter() = default;

  void SavitzkyGolayFilter::updateMembers_()
  {
    frame_size_ = (UInt)param_.getValue("frame_length");
    order_ = (UInt)param_.getValue("polynomial_order");

    //recalculate coefficients
    if (!Math::isOdd(frame_size_))
    {
      frame_size_ += 1;
    }
    if (frame_size_ <= order_)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The degree of the polynomial has to be less than the frame length.", String(order_));
    }
    coeffs_.resize(frame_size_ * (frame_size_ / 2 + 1));

    for (int nl = 0; nl <= (int) (frame_size_ / 2); ++nl)
    {
      int nr = frame_size_ - 1 - nl;

      // compute a Vandermonde matrix whose columns are powers of the vector [-nL,...,nR]
      Eigen::MatrixXd A (frame_size_, order_ + 1);
      for (int i = -nl; i <= nr; i++)
      {
        for (int j = 0; j <= static_cast<int>(order_); j++)
        {
          A(i + nl, j) = std::pow((float)i, j); // pow(int, int) is not defined
        }
      }

      // compute the singular-value decomposition of A
      Eigen::JacobiSVD<Eigen::MatrixXd> svd (A, Eigen::ComputeThinU | Eigen::ComputeThinV);

      Eigen::VectorXd B (order_ + 1);
      for (UInt i = 0; i <= order_; ++i)
      {
        B(i) = svd.matrixV()(0, i) / svd.singularValues()(i);
      }

      // compute B*transpose(U)*b, where b is the unit vector b=[1 0 ... 0]
      for (UInt i = 0; i < frame_size_; ++i)
      {
        coeffs_[(nl + 1) * frame_size_ - i - 1] = 0;
        for (UInt j = 0; j <= order_; ++j)
        {
          coeffs_[(nl + 1) * frame_size_ - i - 1] += B(j) * svd.matrixU()(i, j);
        }
      }
    }
  }
}
