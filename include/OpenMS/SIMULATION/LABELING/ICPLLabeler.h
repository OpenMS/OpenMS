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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_LABELING_ICPLLABELER_H
#define OPENMS_SIMULATION_LABELING_ICPLLABELER_H

#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>

namespace OpenMS
{

  /**
    @brief Simulate ICPL experiments

    Add modified features to MS1 scans. Highly Experimental.

    @htmlinclude OpenMS_ICPLLabeler.parameters
  */
  class OPENMS_DLLAPI ICPLLabeler
    : public BaseLabeler
  {
  public:

    /// default constructor
    ICPLLabeler();

    /// destructor
    virtual ~ICPLLabeler();

    /// create new object (needed by Factory)
    static BaseLabeler* create()
    {
        return new ICPLLabeler();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
        return "ICPL";
    }

    // redeclaration of virtual methods
    void preCheck(Param &param) const;

    void setUpHook(FeatureMapSimVector & /* channels */);

    void postDigestHook(FeatureMapSimVector & /* features_to_simulate */);

    void postRTHook(FeatureMapSimVector & /* features_to_simulate */);

    void postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */);

    void postIonizationHook(FeatureMapSimVector & /* features_to_simulate */);

    void postRawMSHook(FeatureMapSimVector & /* features_to_simulate */);

    void postRawTandemMSHook(FeatureMapSimVector & /* features_to_simulate */, MSSimExperiment & /* simulated map */);

  protected:
    void addModificationToPeptideHit_(Feature& feature, const String& modification) const;

    void addLabelToProteinHits_(FeatureMapSim& features, const String& label) const;

    Feature mergeFeatures_(Feature& labeled_channel_feature, const AASequence& unmodified_sequence, Map<String, Feature>& unlabeled_features_index) const;

    String light_channel_label_;
    String medium_channel_label_;
    String heavy_channel_label_;

    void updateMembers_();

    String getUnmodifiedAASequence_(const Feature& feature, const String& label) const;
  };
} // namespace OpenMS

#endif //#ifndef OPENMS_SIMULATION_LABELING_ICPLLABELER_H
