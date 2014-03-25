// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Christian Ehrlich $
// $Authors: Christian Ehrlich $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_MISC_SPLINE2D_H_
#define OPENMS_MATH_MISC_SPLINE2D_H_


#include <Eigen/Core>
#include <unsupported/Eigen/Splines>

#include <vector>
#include <cassert>
#include <map>


namespace OpenMS {
  /**
      @brief Wrapper for Spline interpolation.
  */
  template<typename ValType = double>
  class Spline2d
  {
  private:
      typedef Eigen::Matrix< ValType, Eigen::Dynamic, Eigen::Dynamic > MatrixT;
      typedef Eigen::Matrix< ValType, Eigen::Dynamic, 1 > VectorT;

  public:
    typedef ValType value_type;

    /** create a spline from two vectors.
     *
     *  One vector holds the x-coordinates the other the y-coordinates.
     *  Coordinates must match by index. Vectors must be the same size.
     *  Vectors must be sorted by x-coordinate. */
    Spline2d(unsigned degree, const std::vector<ValType>& x,
        const std::vector<ValType>& y)
    {
      assert(!x.empty());
      assert(x.size() == y.size());
      size_t num_data_points = x.size();
      MatrixT raw_values (2, num_data_points);
      for (unsigned i=0; i<num_data_points; ++i)
      {
        raw_values(0, i) = x.at(i);
        raw_values(1, i) = y.at(i);
      }
      initialize(degree, raw_values);
    }

    /** create a spline from a std::map */
    Spline2d(unsigned degree, const std::map<ValType, ValType>& m)
    {
      size_t num_data_points = m.size();
      MatrixT raw_values (2, num_data_points);
      typename std::map<ValType, ValType>::const_iterator map_it;
      unsigned colIdx = 0;
      for (map_it = m.begin(); map_it != m.end(); ++map_it, ++colIdx)
      {
        raw_values(0, colIdx) = map_it->first;
        raw_values(1, colIdx) = map_it->second;
      }
      initialize(degree, raw_values);
    }

    /** factory method to create a spline from a Eigen::Matrix */
    Spline2d(unsigned degree, const MatrixT& raw_values)
    {
      initialize(degree, raw_values);
    }

    /** evaluate spline at position x */
    ValType eval(ValType x) const
    {
      return spline_( getNormIndex(x) )(1);
    }
    /** evaluates the spline derivative of the given order */
    ValType derivatives(ValType x, unsigned order) const
    {
      return spline_.derivatives( getNormIndex(x), order )(3);
    }

  private:

    Eigen::Spline<ValType,2> spline_;
    ValType minXCoeff_;
    ValType maxXCoeff_;

    /** calculate the normalized index ([0,1]) for Eigen spline access */
    ValType
    getNormIndex(ValType x) const
    {
      return (x - minXCoeff_) / (maxXCoeff_ - minXCoeff_);
    }

    void
    initialize(unsigned degree, const MatrixT& raw_values)
    {
      minXCoeff_ = raw_values.row(0).minCoeff();
      maxXCoeff_ = raw_values.row(0).maxCoeff();

      // setup spline
      size_t num_data_points = raw_values.row(0).size();
      assert(num_data_points > degree);
      VectorT uvalues (num_data_points);
      for (Eigen::DenseIndex j=0; j<(Eigen::DenseIndex) num_data_points; ++j)
      {
          uvalues(j) = getNormIndex(raw_values(0, j));
      }
      spline_ = Eigen::SplineFitting< Eigen::Spline<ValType,2> >::Interpolate(raw_values, degree, uvalues.transpose());
    }
  };
}//namespace

#endif /* SPLINE2D_H_ */
