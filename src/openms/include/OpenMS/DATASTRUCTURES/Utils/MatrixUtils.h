/*
 * MatrixUtils.h
 *
 *  Created on: Oct 14, 2013
 *      Author: Hans-Christian Ehrlich
 */

#ifndef MATRIXUTILS_H_
#define MATRIXUTILS_H_

#include <OpenMS/DATASTRUCTURES/Matrix.h>

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  /** Matrix utility functions.
   *
   */

  typedef boost::shared_ptr< const Eigen::MatrixXd > EigenMatrixXdPtr;
  typedef boost::shared_ptr< Eigen::MatrixXd > MutableEigenMatrixXdPtr;

  static EigenMatrixXdPtr
  convertOpenMSMatrix2EigenMatrixXd( const Matrix<double>& m )
  {
    MutableEigenMatrixXdPtr em ( new Eigen::MatrixXd(m.rows(), m.cols()) );
    for (unsigned i=0; i<m.rows(); ++i)
    {
      for (unsigned j=0; j<m.cols(); ++j)
      {
        (*em)(i,j) = m(i,j);
      }
    }
    return em;
  }

  static bool
  matrixIsIdentityMatrix(const Matrix<double>& channel_frequency)
  {
    bool is_identity = true;

    for (Matrix<double>::SizeType i = 0; i < channel_frequency.rows(); ++i)
    {
      for (Matrix<double>::SizeType j = 0; j < channel_frequency.rows(); ++j)
      {
        // check if the entries are those of a identity matrix;
        // i==j -> m(i,j) == 1.0 && i!=j -> m(i,j) == 0.0
        if ((i == j && channel_frequency(i, j) != 1.0) || channel_frequency(i, j) != 0.0)
        {
          is_identity = false;
          break;
        }
      }
      // leave outer loop if we have reached the abortion cirteria
      if (!is_identity) break;
    }
    return is_identity;
  }
}//namespace

#endif /* MATRIXUTILS_H_ */
