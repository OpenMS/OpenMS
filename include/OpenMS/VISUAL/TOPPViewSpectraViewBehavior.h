// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPVIEWSPECTRAVIEWBEHAVIOR_H
#define OPENMS_VISUAL_TOPPVIEWSPECTRAVIEWBEHAVIOR_H

#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <vector>
#include <OpenMS/VISUAL/TOPPViewBehaviorInterface.h>

namespace OpenMS
{
  class TOPPViewBase;

  class TOPPViewSpectraViewBehavior :
      public TOPPViewBehaviorInterface
  {
    Q_OBJECT
    ///@name Type definitions
    //@{
      //Feature map type
      typedef LayerData::FeatureMapType FeatureMapType;
      //Feature map managed type
      typedef LayerData::FeatureMapSharedPtrType FeatureMapSharedPtrType;

      //Consensus feature map type
      typedef LayerData::ConsensusMapType ConsensusMapType;
      //Consensus  map managed type
      typedef LayerData::ConsensusMapSharedPtrType ConsensusMapSharedPtrType;

      //Peak map type
      typedef LayerData::ExperimentType ExperimentType;
      //Main managed data type (experiment)
      typedef LayerData::ExperimentSharedPtrType ExperimentSharedPtrType;
      ///Peak spectrum type
      typedef ExperimentType::SpectrumType SpectrumType;
    //@}

    public:
      /// Construct the behaviour with its parent
      TOPPViewSpectraViewBehavior(TOPPViewBase* parent);

    public slots:
      /// Behavior for showSpectraumAs1D
      virtual void showSpectrumAs1D(int index);

      /// Behavior for activate1DSpectrum
      virtual void activate1DSpectrum(int index);

      /// Behavior for deactivate1DSpectrum
      virtual void deactivate1DSpectrum(int index);

      /// Slot for behavior activation
      virtual void activateBehavior();

      /// Slot for behavior deactivation
      virtual void deactivateBehavior();
    private:
      TOPPViewBase* tv_;
  };
}
#endif // OPENMS_VISUAL_TOPPVIEWSPECTRAVIEWBEHAVIOR_H
