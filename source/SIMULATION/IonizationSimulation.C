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

#include<OpenMS/SIMULATION/IonizationSimulation.h>

namespace OpenMS {

  IonizationSimulation::IonizationSimulation()
    : DefaultParamHandler("IonizationSimulation")
  {
    setDefaultParams_();
  }

  IonizationSimulation::IonizationSimulation(const IonizationSimulation& source)
    : DefaultParamHandler(source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();  
  }

  IonizationSimulation::IonizationSimulation(const gsl_rng * random_generator)
    : DefaultParamHandler("IonizationSimulation"), rnd_gen_(random_generator)
  {
    setDefaultParams_();
  }
  
  IonizationSimulation& IonizationSimulation::operator = (const IonizationSimulation& source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
    return *this;
  }
  
  IonizationSimulation::~IonizationSimulation()
  {}

  void IonizationSimulation::ionize(FeatureMap< > & features)
  {
    switch (ionization_type) {
      case MALDI:
        ionize_maldi(features);
        break;
      case ESI:
        ionize_esi(features);
        break;
    }
  }

  
  void IonizationSimulation::setDefaultParams_()
  {
    defaults_.setValue("ionization_type", "ESI", "Type of Ionization (MALDI or ESI)");
    defaultsToParam_();
  }
  
  void IonizationSimulation::updateMembers_()
  {
    String type = param_.getValue("ionization_type");
    if(type == "ESI")
    {
      ionization_type = ESI;    
    }
    else if(type == "MALDI")
    {
      ionization_type = MALDI;
    }
    else
    {
      /// unsupported ionization model
      throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IonizationSimulation got invalid Ionization type" + type);
    }
  }
  
  void IonizationSimulation::ionize_esi(FeatureMap< > & features)
  {
    // iterate over all features
    for( Size i = 0 ; i < features.size() ; ++i)
    {
      // sample different charge states
      
      // sample abundances for different charge states
      
      // split feature into different sub-features charge and modified abundance
      features[i].setMetaValue("charge", 1);
    }    
  }
  
  void IonizationSimulation::ionize_maldi(FeatureMap< > & features)
  {
    // set all charges to 1
    for( Size i = 0 ; i < features.size() ; ++i)
    {
      features[i].setMetaValue("charge", 1);
    }
  }
}

