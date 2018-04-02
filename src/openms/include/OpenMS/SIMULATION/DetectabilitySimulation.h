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

#ifndef OPENMS_SIMULATION_DETECTABILITYSIMULATION_H
#define OPENMS_SIMULATION_DETECTABILITYSIMULATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS
{
  /**
   @brief Simulates peptide detectability

   The peptide detectability is predicted based on a support-vector machine. Alternatively
   the detectability can be set to a default value for all peptides if no model is given.

   @htmlinclude OpenMS_DetectabilitySimulation.parameters

   @ingroup Simulation
  */
  class OPENMS_DLLAPI DetectabilitySimulation :
    public DefaultParamHandler
  {

public:
    /** @name Constructors and Destructors
      */
    //@{
    /// Constructor taking a random generator
    DetectabilitySimulation();

    /// Copy constructor
    DetectabilitySimulation(const DetectabilitySimulation& source);

    /// Destructor
    ~DetectabilitySimulation() override;
    //@}

    /// Assignment operator
    DetectabilitySimulation& operator=(const DetectabilitySimulation& source);

    /**
     @brief Filters the given peptide features for detectability

     Based on the provided method (SVM or simple) all peptide features are
     removed that do not have a sufficient peptide detectability.

     @param features Feature map that will be filtered for detectability
     */
    void filterDetectability(SimTypes::FeatureMapSim& features);


    void predictDetectabilities(std::vector<String>& peptides_vector, std::vector<double>& labels,
                                std::vector<double>& detectabilities);
private:
    /// Set default parameters
    void setDefaultParams_();

    /// Synchronize members with param class
    void updateMembers_() override;

    /// Minimum allowed detectability likelihood of a peptide
    double min_detect_;

    /// Name of the svm model file
    OpenMS::String dt_model_file_;

protected:
    /// Filter the feature map using a svm model
    void svmFilter_(SimTypes::FeatureMapSim&);

    /// Do not filter the feature map, just set the detectability to a default value
    void noFilter_(SimTypes::FeatureMapSim&);

  };

}

#endif
