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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_RTSIMULATION_H
#define OPENMS_SIMULATION_RTSIMULATION_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS
{
  /**
   @brief Simulates/Predicts retention times for peptides or peptide separation.

   The retention times for the different peptides are determined based on
   a SVM model or are all set to -1 in case of simulations without a HPLC
   column.

   @htmlinclude OpenMS_RTSimulation.parameters

   @ingroup Simulation
  */
  class OPENMS_DLLAPI RTSimulation :
    public DefaultParamHandler
  {


public:
    /** @name Constructors and Destructors
      */
    //@{

    /// Default constructor
    RTSimulation();

    /// Constructor taking a random generator
    explicit RTSimulation(SimTypes::MutableSimRandomNumberGeneratorPtr random_generator);

    /// Copy constructor
    RTSimulation(const RTSimulation& source);

    /// Destructor
    ~RTSimulation() override;
    //@}

    /// Assignment operator
    RTSimulation& operator=(const RTSimulation& source);

    /**
     @brief Predict retention times for given peptide features based for HPLC or CE

     @param features Feature map for which the retention times will be predicted
     */
    void predictRT(SimTypes::FeatureMapSim& features);

    /**
     @brief Set retention times randomly for given contaminants
     */
    void predictContaminantsRT(SimTypes::FeatureMapSim&);

    /**
     @brief Returns true if a RT column was simulated
     */
    bool isRTColumnOn() const;

    /// Wrapper for the SVM RT Prediction (HPLC) using AASequences
    void wrapSVM(std::vector<AASequence>& peptide_sequences, std::vector<double>& predicted_retention_times);

    SimTypes::SimCoordinateType getGradientTime() const;

    /// Size experiment and assign retention time grid
    void createExperiment(SimTypes::MSSimExperiment& experiment);

private:
    /// Set default parameters
    void setDefaultParams_();

    /// Simply set all retention times to -1
    void noRTColumn_(SimTypes::FeatureMapSim&);

    /// smoothes the simulated distortion for the elution profiles with a moving average filter of size 3
    void smoothRTDistortion_(SimTypes::MSSimExperiment& experiment);

    /**
      Wrapper for the Migration time calculation (CE)

      @param features will get modified with metavalue "RT_CE_width_factor",
             describing widening of MT shape.
      @param predicted_retention_times will contain afterwards the predicted retention times.
    */
    void calculateMT_(SimTypes::FeatureMapSim& features, std::vector<double>& predicted_retention_times);

    void getChargeContribution_(Map<String, double>& q_cterm,
                                Map<String, double>& q_nterm,
                                Map<String, double>& q_aa_basic,
                                Map<String, double>& q_aa_acidic);

    // MEMBERS:

    // Name of the svm model file
    OpenMS::String rt_model_file_;

    /// Total gradient time
    SimTypes::SimCoordinateType total_gradient_time_;

    /// gradient ranges

    /// Minimal observed gradient time
    SimTypes::SimCoordinateType gradient_min_;
    /// Maximal observed gradient time
    SimTypes::SimCoordinateType gradient_max_;

    /// bin size in rt dimension
    SimTypes::SimCoordinateType rt_sampling_rate_;

    /// EGH tau value
    double egh_tau_location_;
    /// EGH tau scale parameter of the lorentzian variation
    double egh_tau_scale_;

    /// EGH sigma value
    double egh_variance_location_;
    /// EGH sigma scale parameter of the lorentzian variation
    double egh_variance_scale_;

protected:
    /// Random number generator
    SimTypes::MutableSimRandomNumberGeneratorPtr rnd_gen_;

    /// Synchronize members with param class
    void updateMembers_() override;

  };

}

#endif
