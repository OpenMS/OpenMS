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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPVIEWIDENTIFICATIONVIEWBEHAVIOR_H
#define OPENMS_VISUAL_TOPPVIEWIDENTIFICATIONVIEWBEHAVIOR_H

#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <vector>
#include <OpenMS/VISUAL/TOPPViewBehaviorInterface.h>

namespace OpenMS
{
  class TOPPViewBase;

  class TOPPViewIdentificationViewBehavior :
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
      TOPPViewIdentificationViewBehavior(TOPPViewBase* parent);

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

      void setVisibleArea1D(DoubleReal l, DoubleReal h);

    private:
      /// Adds labels for the provided precursors to the 1D spectrum
      void addPrecursorLabels1D_(const std::vector<Precursor>& pcs);

      /// Removes the precursor labels for from the specified 1D spectrum
      void removeTemporaryAnnotations_(Size spectrum_index);

      /// Adds a theoretical spectrum as set from the preferences dialog for the peptide hit.
      void addTheoreticalSpectrumLayer_(const PeptideHit& ph);

      /// removes all layer with theoretical spectrum generated in identification view
      void removeTheoreticalSpectrumLayer_();

    private:
      TOPPViewBase* tv_;
      /// Used to check which annotation handles have been added automaticaly by the identification view. Ownership
      /// of the AnnotationItems has the Annotation1DContainer 
      std::vector<Annotation1DItem* > temporary_annotations_;
  };
}

#endif // OPENMS_VISUAL_TOPPVIEWIDENTIFICATIONVIEWBEHAVIOR_H
