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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_IONIZATIONSIMULATION_H
#define OPENMS_SIMULATION_IONIZATIONSIMULATION_H

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/SIMULATION/SimTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

// STL includes
#include <set>

namespace OpenMS
{

  /**
   @brief Simulates Protein ionization

   Supports ESI and MALDI. The abundance values are distributed among
   the charge states based on a binomial distribution for the ESI and
   based on discrete distribution for MALDI.
     In ESI mode, this class also supports different adduct types in addition to H+
     (e.g. NH4+, K+) which can be specified by the user and influence
     the mass and induce more charge variation.

   @htmlinclude OpenMS_IonizationSimulation.parameters

   @ingroup Simulation
  */
  class OPENMS_DLLAPI IonizationSimulation :
    public DefaultParamHandler,
    public ProgressLogger
  {

public:
    /// possible ionization methods
    typedef enum
    {
      MALDI,
      ESI
    } IonizationType;

    /** @name Constructors and Destructors
      */
    //@{
    ///
    explicit IonizationSimulation(const SimRandomNumberGenerator &);

    /// Copy constructor
    IonizationSimulation(const IonizationSimulation & source);

    /// Destructor
    virtual ~IonizationSimulation();
    //@}

    /// Assignment operator
    IonizationSimulation & operator=(const IonizationSimulation & source);

    /**
     @brief Ionize all peptide features inside the Feature-Map

     Depending on the parameters the passed peptide features are ionized by MALDI
     or by ESI.

     @param features FeatureMap which will be ionized
         @param charge_consensus ConsensusMap which groups children(=charge variants) of input-features
         @param experiment MSSimExperiment map which contains the simulated experiment
     */
    void ionize(FeatureMapSim & features, ConsensusMap & charge_consensus, MSSimExperiment & experiment);

private:
    /// Default constructor
    IonizationSimulation();

    class CompareCmpByEF_;

    /// ionize using ESI
    void ionizeEsi_(FeatureMapSim &, ConsensusMap & charge_consensus);

    /// ionize using MALDI
    void ionizeMaldi_(FeatureMapSim &, ConsensusMap & charge_consensus);

    /// check if feature is within mz bounds of detector
    inline bool isFeatureValid_(const Feature & feature);

    /// set meta values, mz etc after adducts are ready
    void setFeatureProperties_(Feature & f,
                               const DoubleReal & adduct_mass,
                               const String & adduct_formula,
                               const SimChargeType charge,
                               const SimIntensityType new_intensity,
                               const Size parent_index);

    /// set defaults
    void setDefaultParams_();

    /// Synchronize members with param class
    void updateMembers_();

    /**
     @brief counts all basic residues inside the amino acid sequence to give an upper bound on the maximal charge during ESI ionization

     The N-term contributes +1 always. All other ionizable residues (according to param "esi:ionized_residues") in the sequence are summed up.
    */
    UInt countIonizedResidues_(const AASequence &) const;


    // Members //

    /// ESI or MALDI ionization
    IonizationType ionization_type_;

    /*
 @brief List of residues that are counted as basic during execution of countBasicResidues_
*/
    std::set<String> basic_residues_;

    /**
     @brief Probability for the binomial distribution of ESI charge states
     */
    DoubleReal esi_probability_;

    /**
     @brief Discrete distribution of impure charge adducts like Na+, K+, Ca++ etc besides the usual H+
    */
    // important: leave that as vector<double> because gsl expects 'double' and not 'DoubleReal' (which might be something different)
    std::vector<double> esi_impurity_probabilities_;


    /**
     @brief Corresponding table to @p esi_impurity_probabilities_ holding the actual element and its charge
    */
    Adduct::AdductsType esi_adducts_;

    /**
     @brief Maximal charge that any impure adduct from parameter list has
    */
    Size max_adduct_charge_;

    /**
     @brief Preprocessed table of discrete distribution (MALDI charges)
    */
    // important: leave that as vector<double> because gsl expects 'double' and not 'DoubleReal' (which might be something different)
    std::vector<double> maldi_probabilities_;

    /// Maximum m/z detected by mass analyser
    SimCoordinateType maximal_mz_measurement_limit_;
    /// Minimum m/z detected by mass analyser
    SimCoordinateType minimal_mz_measurement_limit_;

protected:
    /// Random number generator
    SimRandomNumberGenerator const * rnd_gen_;
  };

}

#endif
