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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_MSSIM_H
#define OPENMS_SIMULATION_MSSIM_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/SIMULATION/SimTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{


  /**
   @brief Central class for simulation of mass spectrometry experiments

   This implementation is an extended and rewritten version of the concepts and ideas presented in:<br>
   <p>
    Ole Schulz-Trieglaff, Nico Pfeifer, Clemens Gropl, Oliver Kohlbacher, and Knut Reinert.<br>
    LC-MSsim - A simulation software for liquid chromatography mass spectrometry data.<br>
    <em>BMC Bioinformatics</em> 9:423, 2008.
   </p>

   @htmlinclude OpenMS_MSSim.parameters

   @see DetectabilitySimulation
   @see DigestSimulation
   @see IonizationSimulation
   @see RawMSSignalSimulation
   @see RTSimulation

   @ingroup Simulation
  */
  class OPENMS_DLLAPI MSSim :
    public DefaultParamHandler,
    public ProgressLogger
  {

public:
    /** @name Constructors and Destructors
      */
    //@{
    /// Default constructor
    MSSim();

    /// Destructor
    virtual ~MSSim();
    //@}

    /**
     @brief General purpose function to simulate a mass spectrometry run

     @param rnd_gen GSL random number generator which will be passed to the different classes
     @param peptides List of peptides and abundances that will be simulated
     */
    void simulate(const SimRandomNumberGenerator & rnd_gen, SampleChannels & peptides);

    /// Access the simulated experiment
    MSSimExperiment const & getExperiment() const;

    /// Access the simulated features
    FeatureMapSim const & getSimulatedFeatures() const;

    /// Access the charge consensus map of simulated features
    ConsensusMap & getChargeConsensus();

    /// Access the contaminants feature map of simulated features
    FeatureMapSim const & getContaminants() const;

    /// Access the labeling consensus map of simulated features
    ConsensusMap & getLabelingConsensus();

    /// Access the labeling consensus map of simulated features
    MSSimExperiment const & getPeakMap() const;

    /// Returns the default parameters for simulation including the labeling technique with name @p labeling_name
    Param getParameters() const;

protected:
    /// handle global params
    void syncParams_(Param & p, bool to_outer);

    /// Convert a list of peptides with given abundance values into a FeatureMap
    void createFeatureMap_(const SampleProteins & peptides, FeatureMapSim & features, Size map_index);

private:
    /// Synchronize members with param class
    void updateMembers_();

    MSSimExperiment experiment_;

    MSSimExperiment peak_map_;

    FeatureMapSimVector feature_maps_;

    ConsensusMap consensus_map_;

    FeatureMapSim contaminants_map_;

    /// Labeling functionality
    BaseLabeler * labeler_;
  };

}

#endif
