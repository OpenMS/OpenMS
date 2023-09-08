// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakShape.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

//#define DEBUG_DECONV
#include <vector>

namespace OpenMS
{

  namespace OptimizationFunctions
  {

    /**
         @brief Class for the penalty factors used during the optimization.

         A great deviation (squared deviation) of a peak shape's position or its left or right width parameter can be penalised.
         During the optimization negative heights may occur, they are penalised as well.
    */
    struct OPENMS_DLLAPI PenaltyFactorsIntensity :
      public PenaltyFactors
    {
      PenaltyFactorsIntensity() :
        PenaltyFactors(), height(0){}
      PenaltyFactorsIntensity(const PenaltyFactorsIntensity & p) :
        PenaltyFactors(p), height(p.height) {}
      inline PenaltyFactorsIntensity & operator=(const PenaltyFactorsIntensity & p)
      {
        height = p.height;
        pos = p.pos;
        lWidth = p.lWidth;
        rWidth = p.rWidth;

        return *this;
      }

      ~PenaltyFactorsIntensity(){}

      double height;


    };



  } //namespace OptimizationFunctions

  /**
        @brief This class provides the deconvolution of peak regions using non-linear optimization.

        Given a vector of peak shapes, this class optimizes all peak shapes parameters using a non-linear optimization.
        For the non-linear optimization we use the Levenberg-Marquardt algorithm.
        There are a few constraints for the parameters: the positions are equidistant according to the peptide
        mass rule, e.g. two consecutive isotopic peaks are 1.003/charge away from each other. Besides the
        peaks have all the same left and right width, respectively.

        @htmlinclude OpenMS_OptimizePeakDeconvolution.parameters
    */
  class OPENMS_DLLAPI OptimizePeakDeconvolution :
    public DefaultParamHandler
  {
public:
    /** @name Type definitions
     */
    //@{
    typedef std::vector<Peak1D> RawDataVector;
    typedef RawDataVector::iterator PeakIterator;
    //@}

    /**
         @brief Class containing the data needed for optimization.
    */
    struct Data
    {
      std::vector<PeakShape> peaks;
      std::vector<double> positions;
      std::vector<double> signal;
      OptimizationFunctions::PenaltyFactorsIntensity penalties;
      Int charge;
    };



    /** @name Constructors and Destructor
     */
    //@{
    ///Constructor
    OptimizePeakDeconvolution();

    /// Copy-Constructor
    OptimizePeakDeconvolution(const OptimizePeakDeconvolution & opt) :
      DefaultParamHandler(opt),
      penalties_(opt.penalties_),
      charge_(opt.charge_){}

    ///Destructor
    ~OptimizePeakDeconvolution() override{}
    //@}

    /**	@name Assignment
     */
    //@{
    inline OptimizePeakDeconvolution & operator=(const OptimizePeakDeconvolution & opt)
    {
      DefaultParamHandler::operator=(opt);
      penalties_ = opt.penalties_;
      charge_ = opt.charge_;

      return *this;
    }

    //@}


    /**	Accessors
     */
    //@{
    /// Non-mutable access to the penalty parameter
    inline const OptimizationFunctions::PenaltyFactorsIntensity & getPenalties() const { return penalties_; }
    /// Mutable access to the penalty parameter
    inline void setPenalties(const OptimizationFunctions::PenaltyFactorsIntensity & penalties)
    {
      penalties_ = penalties;
      param_.setValue("penalties:left_width", penalties_.lWidth);
      param_.setValue("penalties:right_width", penalties_.rWidth);
      param_.setValue("penalties:height", penalties_.height);
      param_.setValue("penalties:position", penalties_.pos);
    }

    /// Non-mutable access to the charge
    inline Int getCharge() const { return charge_; }
    /// Mutable access to the charge
    inline void setCharge(const Int charge) { charge_ = charge; }
    //@}


    /// Performs a nonlinear optimization of the peaks that belong to the current isotope pattern
    bool optimize(std::vector<PeakShape> & peaks, Data & data);
    Size getNumberOfPeaks_(Int charge, std::vector<PeakShape> & temp_shapes, Data & data);

protected:
    // Penalty factors for some parameter in the optimization
    OptimizationFunctions::PenaltyFactorsIntensity penalties_;

    /// Charge state of the current isotope pattern
    Int charge_;

    /// distance between two isotopic peaks
    static const double dist_;

    /// A function to determine the number of peaks that lie in the current m/z interval given the distance between the peaks by the current charge state.
    void setNumberOfPeaks_(Data & data, const std::vector<PeakShape> & temp_shapes, Int charge);

    void updateMembers_() override;
  }; // class

} // namespace OpenMS


