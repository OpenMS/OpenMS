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
// $Maintainer: Chris Bielow$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/RawTandemMSSignalSimulation.h>

namespace OpenMS {

  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation(const gsl_rng * random_generator)
  : DefaultParamHandler("RawTandemMSSignalSimulation"), rnd_gen_(random_generator)
  {
    setDefaultParams_();
    updateMembers_();
  }


  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation()
    : DefaultParamHandler("RawTandemMSSignalSimulation")
  {
    setDefaultParams_();
    updateMembers_();
  }

  RawTandemMSSignalSimulation::RawTandemMSSignalSimulation(const RawTandemMSSignalSimulation& source)
    : DefaultParamHandler(source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
  }

  RawTandemMSSignalSimulation& RawTandemMSSignalSimulation::operator = (const RawTandemMSSignalSimulation& source)
  {
    setParameters( source.getParameters() );
    rnd_gen_ = source.rnd_gen_;
    updateMembers_();
    return *this;
  }

  RawTandemMSSignalSimulation::~RawTandemMSSignalSimulation()
  {}

  void RawTandemMSSignalSimulation::setDefaultParams_()
  {
    // noise params
    // m/z error
    defaults_.setValue("mz:error_mean",0.0,"Average systematic m/z error (Da)");
 
    defaultsToParam_();
  }

  void RawTandemMSSignalSimulation::updateMembers_()
  {
    double tmp    =  param_.getValue("peak_fwhm");
   
  }

  void RawTandemMSSignalSimulation::generateRawSignals(FeatureMapSim & features, MSSimExperiment & experiment)
  {
 
  }
  
}
