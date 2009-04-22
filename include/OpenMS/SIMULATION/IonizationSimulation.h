// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_IONIZATIONSIMULATION_H
#define OPENMS_SIMULATION_IONIZATIONSIMULATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// STL includes
#include <set>

namespace OpenMS {

  /**
   @brief 
   @ingroup Simulation
  */
  class OPENMS_DLLAPI IonizationSimulation
    : public DefaultParamHandler
  {

  public: 
    /// possible ionization methods
    typedef enum {
      MALDI,
      ESI
    } IonizationType;
    
    /** @name Constructors and Destructors
      */
    //@{
    /// 
    IonizationSimulation(const gsl_rng* );
    
    /// Copy constructor
    IonizationSimulation(const IonizationSimulation& source);

    /// Destructor
    virtual ~IonizationSimulation();
    //@}

    IonizationSimulation& operator = (const IonizationSimulation& source);

    void ionize(FeatureMap< > &);

  private:  
    /// Default constructor
    IonizationSimulation();
    
    void ionize_esi(FeatureMap< > &);
    
    void ionize_maldi(FeatureMap< > &);

    /// set defaults
    void setDefaultParams_();
    
    /// Synchronize members with param class
		void updateMembers_();        
    
    /// ESI or MALDI ionization
    IonizationType ionization_type;
    
    /**
     @brief counts all basic residues inside the amino acid sequence to give an upper bound on the maximal charge during ESI ionization
    */
    UInt countBasicResidues_(const AASequence & ) const;
    
    /*
     @brief List of residues that are counted as basic during execution of countBasicResidues_
    */
    std::set<String> basic_residues_;
    
    /**
     @brief Probability for the binomial distribution of ESI charge states
     */
    DoubleReal esi_probability;
    
    /**
     @brief List of probabilities for each charge state during MALDI ionization
     */
    DoubleList maldi_probabilities;
    
  protected:
		/// Random number generator
		const gsl_rng* rnd_gen_;
  };

}

#endif
