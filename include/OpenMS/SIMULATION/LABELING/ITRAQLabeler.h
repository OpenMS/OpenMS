// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_LABELING_ITRAQLABELER_H
#define OPENMS_SIMULATION_LABELING_ITRAQLABELER_H

#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS
{

  /**
  @brief Abstract base class for all kinds of labeling techniques
  */
  class OPENMS_DLLAPI ITRAQLabeler
    : public BaseLabeler
  {
  public:

		typedef ItraqConstants::ChannelInfo ChannelInfo;
		typedef ItraqConstants::ChannelMapType ChannelMapType;
		typedef ItraqConstants::IsotopeMatrices IsotopeMatrices;

    /// default constructor
    ITRAQLabeler();

    /// destructor
    virtual ~ITRAQLabeler();

    /// create new object (needed by Factory)
    static BaseLabeler* create()
    {
        return new ITRAQLabeler();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
        return "itraq";
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

    Feature mergeFeatures_(Feature& labeled_channel_feature, const AASequence& unmodified_sequence, std::map<AASequence, Feature>& unlabeled_features_index) const;

    /// Synchronize members with param class
		void updateMembers_();

		/// convert meta information from feature into intensity values for iTRAQ
		Matrix<SimIntensityType> getItraqIntensity_(const Feature & f) const;


		// Members:

		/// set to either ItraqConstants::FOURPLEX or ItraqConstants::EIGHTPLEX
		Int itraq_type_;
		
		/// map the channel-name (eg 114) onto its description and the centroid mass
		/// the channel-name is also the id-string in the mapList section of the ConsensusMap
		ChannelMapType channel_map_;	

		/// Matrixes with isotope correction values (one for each plex-type)
		IsotopeMatrices isotope_corrections_;

  };
} // namespace OpenMS

#endif //#ifndef OPENMS_SIMULATION_LABELING_ITRAQLabeler_H
