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
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <map>
#include <vector>

#include <OpenMS/SIMULATION/LABELING/LabelFreeLabeler.h>
#include <OpenMS/METADATA/ProteinHit.h>

using std::map;
using std::vector;

namespace OpenMS
{
  LabelFreeLabeler::LabelFreeLabeler()
    : BaseLabeler()
  {
  }

  LabelFreeLabeler::~LabelFreeLabeler()
  {
  }

  void LabelFreeLabeler::preCheck(Param & /* param */) const
  {
    // no specific requirements
  }

  // merge all channels into the first one
  // no further influence of the simulation process needed
  void LabelFreeLabeler::setUpHook(FeatureMapSimVector & features)
  {
    if(features.size() == 1) return;
    else
    {
      LOG_INFO << "Merging input FASTA files into one. Intensities will be summed up if duplicates occur.";
      FeatureMapSim final_map = mergeProteinIdentificationsMaps_(features);
      features.clear();
      features.push_back(final_map);
    }
  }

  /// Labeling between digestion and rt simulation
  void LabelFreeLabeler::postDigestHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between RT and Detectability
  void LabelFreeLabeler::postRTHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between Detectability and Ionization
  void LabelFreeLabeler::postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling between Ionization and RawMS
  void LabelFreeLabeler::postIonizationHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  /// Labeling after RawMS
  void LabelFreeLabeler::postRawMSHook(FeatureMapSimVector & /* features_to_simulate */)
  {
  }

  void LabelFreeLabeler::postRawTandemMSHook(FeatureMapSimVector &, MSSimExperiment &)
  {

  }
} // namespace OpenMS
