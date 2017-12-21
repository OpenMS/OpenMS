// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_OPTIMIZEPEAKDECONVOLUTION_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_OPTIMIZEPEAKDECONVOLUTION_H

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


#endif
