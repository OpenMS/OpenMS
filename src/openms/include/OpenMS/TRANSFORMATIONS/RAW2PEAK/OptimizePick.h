// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/KERNEL/Peak1D.h>

// forward decl
namespace Eigen
{
    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    class Matrix;
    using MatrixXd = Matrix<double, -1, -1, 0, -1, -1>;
    using VectorXd = Matrix<double, -1, 1, 0, -1, 1>;
}

#include <vector>

namespace OpenMS
{
  namespace OptimizationFunctions
  {
    /// Profile data vector type
    typedef std::vector<Peak1D> RawDataVector;
    /// Profile data iterator type
    typedef RawDataVector::iterator PeakIterator;

    /** @brief Class for the penalty factors used during the optimization.

        A great deviation (squared deviation) of a peak shape's position or its left or right width parameter can be penalised.
        In each iteration the penalty (for each peak shape) is computed by:
                penalty = penalty_pos * pow(p_position - old_position, 2)
                        + penalty_lwidth * pow(p_width_l - old_width_l, 2)
                        + penalty_rwidth * pow(p_width_r - old_width_r, 2);
    */
    struct OPENMS_DLLAPI PenaltyFactors
    {
      PenaltyFactors() :
        pos(0), lWidth(0), rWidth(0) {}
      PenaltyFactors(const PenaltyFactors & p) :
        pos(p.pos), lWidth(p.lWidth), rWidth(p.rWidth) {}
      inline PenaltyFactors & operator=(const PenaltyFactors & p)
      {
        pos = p.pos;
        lWidth = p.lWidth;
        rWidth = p.rWidth;

        return *this;
      }

      ~PenaltyFactors(){}

      /// Penalty factor for the peak shape's position
      double pos;
      /// Penalty factor for the peak shape's left width parameter
      double lWidth;
      /// Penalty factor for the peak shape's right width parameter
      double rWidth;
    };
  }


  /**
    @brief This class provides the non-linear optimization of the peak parameters.

    Given a vector of peak shapes, this class optimizes all peak shapes parameters using a non-linear optimization.
    For the non-linear optimization we use the Levenberg-Marquardt algorithm provided by the Eigen.
  */
  class OPENMS_DLLAPI OptimizePick
  {
public:

    struct Data
    {
      /// Positions and intensity values of the profile data
      std::vector<double> positions;
      std::vector<double> signal;
      /// This container contains the peak shapes to be optimized
      std::vector<PeakShape> peaks;

      OptimizationFunctions::PenaltyFactors penalties;

    };

    class OptPeakFunctor
    {
    public:
      int inputs() const { return m_inputs; }
      int values() const { return m_values; }

      OptPeakFunctor(unsigned dimensions, unsigned num_data_points, const OptimizePick::Data * data)
      : m_inputs(dimensions), m_values(num_data_points), m_data(data) {}

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec);
      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &J);

    private:
      const int m_inputs, m_values;
      const Data * m_data;
    };

    /// Profile data vector type
    typedef std::vector<Peak1D> RawDataVector;
    /// Profile data iterator type
    typedef RawDataVector::iterator PeakIterator;


    /// Constructor
    OptimizePick() :
      max_iteration_(400)
    {}

    /// Constructor to set the penalty factors, the number of optimization iterations as well as the threshold for the absolute and the relative error.
    OptimizePick(const struct OptimizationFunctions::PenaltyFactors & penalties_,
                 const int max_iteration_);

    /// Destructor
    ~OptimizePick();

    /// Non-mutable access to the penalty factors
    inline const struct OptimizationFunctions::PenaltyFactors & getPenalties() const { return penalties_; }
    /// Mutable access to the penalty factors
    inline struct OptimizationFunctions::PenaltyFactors & getPenalties() { return penalties_; }
    /// Mutable access to the penalty factors
    inline void setPenalties(const struct OptimizationFunctions::PenaltyFactors & penalties) { penalties_ = penalties; }

    /// Non-mutable access to the number of iterations
    inline UInt getNumberIterations() const { return max_iteration_; }
    /// Mutable access to the number of iterations
    inline unsigned int & getNumberIterations() { return max_iteration_; }
    /// Mutable access to the number of iterations
    inline void setNumberIterations(const int max_iteration) { max_iteration_ = max_iteration; }

    /// Start the optimization of the peak shapes peaks. The original peak shapes will be substituted by the optimized peak shapes.
    void optimize(std::vector<PeakShape> & peaks, Data & data);


protected:
    /// Penalty factors
    struct OptimizationFunctions::PenaltyFactors penalties_;

    /// Maximum number of iterations during optimization
    unsigned int max_iteration_;
  };
}

