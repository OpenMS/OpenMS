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
// $Authors: Stephan Aiche, Chris Bielow, Sandro Andreotti $
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_RAWTANDEMMSSIGNALSIMULATION_H
#define OPENMS_SIMULATION_RAWTANDEMMSSIGNALSIMULATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS
{

  /**
   @brief Simulates tandem MS signals for a given set of peptides

   Simulates tandem MS signals for a given set of peptides, with charge annotation,
   given detectabilities, predicted retention times and charge values.

   @htmlinclude OpenMS_RawTandemMSSignalSimulation.parameters

   @ingroup Simulation
  */
  class OPENMS_DLLAPI RawTandemMSSignalSimulation :
    public DefaultParamHandler
  {
public:

    /** @name Constructors and Destructors
      */
    //@{

    /// Default constructor (hidden)
    RawTandemMSSignalSimulation();

    /// Constructor taking a random generator
    explicit RawTandemMSSignalSimulation(SimTypes::MutableSimRandomNumberGeneratorPtr rng);

    /// Copy constructor
    RawTandemMSSignalSimulation(const RawTandemMSSignalSimulation& source);

    /// Destructor
    ~RawTandemMSSignalSimulation() override;
    //@}

    RawTandemMSSignalSimulation& operator=(const RawTandemMSSignalSimulation& source);


    /**

     */
    void generateRawTandemSignals(const SimTypes::FeatureMapSim&, SimTypes::MSSimExperiment&, SimTypes::MSSimExperiment&);


protected:

    /// initialize param_ class
    void initParam_();

    void generateMSESpectra_(const SimTypes::FeatureMapSim& features, const SimTypes::MSSimExperiment& experiment, SimTypes::MSSimExperiment& ms2);

    void generatePrecursorSpectra_(const SimTypes::FeatureMapSim& features, const SimTypes::MSSimExperiment& experiment, SimTypes::MSSimExperiment& ms2);

    /// Random number generator
    SimTypes::MutableSimRandomNumberGeneratorPtr  rnd_gen_;
  };

}

#endif
