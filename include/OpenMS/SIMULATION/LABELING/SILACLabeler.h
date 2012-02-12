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

#ifndef OPENMS_SIMULATION_LABELING_SILACLABELER_H
#define OPENMS_SIMULATION_LABELING_SILACLABELER_H

#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>

namespace OpenMS
{


/**
    @brief Simulate SILAC experiments

    Add modified features to MS1 scans.

    @htmlinclude OpenMS_SILACLabeler.parameters
  */
class OPENMS_DLLAPI SILACLabeler
    : public BaseLabeler
{
public:

  /// default constructor
  SILACLabeler();

  /// destructor
  virtual ~SILACLabeler();

  /// create new object (needed by Factory)
  static BaseLabeler* create()
  {
    return new SILACLabeler();
  }

  /// name of the model (needed by Factory)
  static const String getProductName()
  {
    return "SILAC";
  }

  // redeclaration of virtual methods
  void preCheck(Param & /* param */) const;

  void setUpHook(FeatureMapSimVector & /* channels */);

  void postDigestHook(FeatureMapSimVector & /* features_to_simulate */);

  void postRTHook(FeatureMapSimVector & /* features_to_simulate */);

  void postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */);

  void postIonizationHook(FeatureMapSimVector & /* features_to_simulate */);

  void postRawMSHook(FeatureMapSimVector & /* features_to_simulate */);

  void postRawTandemMSHook(FeatureMapSimVector & /* features_to_simulate */, MSSimExperiment & /* simulated map */);

protected:
  void addModificationToPeptideHit_(Feature& feature, const String& modification) const;

  Feature mergeFeatures_(Feature& labeled_channel_feature, const String& unmodified_sequence, Map<String, Feature>& unlabeled_features_index, Int index_channel_id, Int labeled_channel_id) const;

  Feature mergeAllChannelFeatures_(Feature& heavy_channel_feature, const String& unmodified_feature_sequence, Map<String, Feature>& light_channel_feature_index, Map<String, Feature>& medium_channel_feature_index) const;

  String medium_channel_lysine_label_;
  String medium_channel_arginine_label_;

  String heavy_channel_lysine_label_;
  String heavy_channel_arginine_label_;

  bool canModificationBeApplied_(const String& modification_id, const String& aa) const;

  void applyLabelToProteinHit_(FeatureMapSim& protein_hit, const String& arginine_label, const String& lysine_label) const;

  void updateMembers_();

  String getUnmodifiedSequence_(const Feature& feature, const String& arginine_label, const String& lysine_label) const;
};

} // namespace OpenMS

#endif //#ifndef OPENMS_SIMULATION_LABELING_SILACLABELER_H
