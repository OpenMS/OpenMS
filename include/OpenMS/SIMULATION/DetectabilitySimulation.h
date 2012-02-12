// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

   The peptide detectability is predicted based on a support-vector machine. Alternativly
   the detectability can be set to a default value for all peptides if no model is given.

   @htmlinclude OpenMS_DetectabilitySimulation.parameters

   @ingroup Simulation
  */
  class OPENMS_DLLAPI DetectabilitySimulation
    : public DefaultParamHandler
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
    virtual ~DetectabilitySimulation();
    //@}

    /// Assignment operator
    DetectabilitySimulation& operator = (const DetectabilitySimulation& source);

    /**
     @brief Filters the given peptide features for detectibility

     Based on the provided method (SVM or simple) all peptide features are
     removed that do not have a sufficient peptide detectibility.

     @param features Feature map that will be filtered for detectibility
     */
    void filterDetectability(FeatureMapSim & features);


		void predictDetectabilities(std::vector<String>& peptides_vector,std::vector<DoubleReal>& labels,
																std::vector<DoubleReal>& detectabilities);
  private:
    /// Set default parameters
    void setDefaultParams_();

    /// Synchronize members with param class
		void updateMembers_();

    /// Minimum allowed detectability likelihood of a peptide
		DoubleReal min_detect_;

		/// Name of the svm model file
		OpenMS::String dt_model_file_;

  protected:
    /// Filter the feature map using a svm model
    void svmFilter_(FeatureMapSim &);

    /// Do not filter the feature map, just set the detectability to a default value
    void noFilter_(FeatureMapSim &);

  };

}

#endif
