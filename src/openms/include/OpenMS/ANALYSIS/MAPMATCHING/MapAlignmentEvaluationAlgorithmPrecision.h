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
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTEVALUATIONALGORITHMPRECISION_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTEVALUATIONALGORITHMPRECISION_H

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>

namespace OpenMS
{
  /**
      @brief Caap evaluation algorithm to obtain a precision value.

      It evaluates an input consensus map with respect to a ground truth.

      @ingroup MapAlignmentEvaluation
  */
  class OPENMS_DLLAPI MapAlignmentEvaluationAlgorithmPrecision :
    public MapAlignmentEvaluationAlgorithm
  {
public:
    /// Default constructor
    MapAlignmentEvaluationAlgorithmPrecision();

    /// Destructor
    ~MapAlignmentEvaluationAlgorithmPrecision() override;

    /**
        @brief Applies the algorithm
    */
    void evaluate(const ConsensusMap & consensus_map_in, const ConsensusMap & consensus_map_gt, const double & rt_dev, const double & mz_dev, const Peak2D::IntensityType & int_dev, const bool use_charge, double & out) override;

    /// Creates a new instance of this class (for Factory)
    static MapAlignmentEvaluationAlgorithm * create()
    {
      return new MapAlignmentEvaluationAlgorithmPrecision();
    }

    /// Returns the product name (for the Factory)
    static String getProductName()
    {
      return "precision";
    }

private:

    /// Copy constructor intentionally not implemented -> private
    MapAlignmentEvaluationAlgorithmPrecision(const MapAlignmentEvaluationAlgorithmPrecision &);
    /// Assignment operator intentionally not implemented -> private
    MapAlignmentEvaluationAlgorithmPrecision & operator=(const MapAlignmentEvaluationAlgorithmPrecision &);

  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTEVALUATIONALGORITHMPRECISION_H
