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
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMQUANTIFICATION_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMQUANTIFICATION_H

#include <OpenMS/config.h>

//Kernal classes
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MRMFeature.h>

//Analysis classes
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

//Standard library
#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>

namespace OpenMS
{

  /**
    @brief MRMQuantification is a class to handle the relationship between
    runs, transitions, and actual concentrations.
  */
  class OPENMS_DLLAPI MRMQuantification
  public DefaultParamHandler,
  public ProgressLogger
  {

public:    
    //@{
    /// Constructor
    MRMQuantification();

    /// Destructor
    ~MRMQuantification();
    //@}
 
    /**
      @brief This function calculates the ratio between features.

      @param feature_1 feature value of the numerator
      @param feature_2 feature value of the denomenator
      @param feature_name name of the feature to calculate the ratio on
       e.g., peak_apex, peak_area

      @return The ratio.

      @exception Exception::UnableToFit
    */ 
    double calculateRatio(Feature & feature_1, Feature & feature_2, std::string feature_nme);
                   
    /**
      @brief This function calculates the bias of the calibration.

      The bias is defined as the following:
        |actual_concentration - calculated_concentration|/actual_concentration * 100%
      This is in contrast to accuracy, which is defined as the following:
        calculated_concentration/actual_concentration * 100%

      @param actual_concentration the actual concentration of the component
      @param calculated_concentration the calibration curve back calculated concentration 
        of the component

      @return The bias.

      @exception Exception::UnableToFit
    */ 
    double calculateBias(double & actual_concentration, double & calculated_concentration);

    void optimizeCalibrationCurve(std::vector<MSExperiment>);
    
    // members
    std::string run_id;
    std::string transition_id;
    double actual_concentration;

  };

}
#endif // OPENMS_ANALYSIS_OPENSWATH_MRMQUANTIFICATION_H

