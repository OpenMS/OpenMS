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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_LAYERDATA_H
#define OPENMS_VISUAL_LAYERDATA_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>
#include <OpenMS/FILTERING/DATAREDUCTION/DataFilters.h>

#include <boost/shared_ptr.hpp>

#include <vector>
#include <bitset>

namespace OpenMS 
{
	/**
    @brief Class that stores the data for one layer
		
		@ingroup SpectrumWidgets
	*/		
  class LayerData
	{
    public:
      /**	@name Type definitions */
      //@{
      /// Dataset types
      enum DataType
      {
        DT_PEAK,		      ///< Spectrum profile or centroided data
        DT_FEATURE,	      ///< Feature data
        DT_CONSENSUS,     ///< Consensus feature data
        DT_CHROMATOGRAM,  ///< Chromatogram data
        DT_IDENT,         ///< Peptide identification data
        DT_UNKNOWN			  ///< Undefined data type indicating an error
      };

      /// Flags that determine which information is shown.
      enum Flags
      {
        F_HULL,        ///< Features: Overall convex hull
        F_HULLS,       ///< Features: Convex hulls of single mass traces
        F_UNASSIGNED,  ///< Features: Unassigned peptide hits
        P_PRECURSORS,  ///< Peaks: Mark precursor peaks of MS/MS scans
        P_PROJECTIONS, ///< Peaks: Show projections
        C_ELEMENTS,    ///< Consensus features: Show elements
        I_PEPTIDEMZ,   ///< Identifications: m/z source
        SIZE_OF_FLAGS
      };

      /// Actual state of each flag
      std::bitset<SIZE_OF_FLAGS> flags;

      /// Label used in visualization
      enum LabelType
      {
        L_NONE,							///< No label is displayed
        L_INDEX,						///< The element number is used
        L_META_LABEL,				///< The 'label' meta information is used
        L_ID,								///< The best peptide hit of the first identification run is used
        L_ID_ALL,						///< All peptide hits of the first identification run are used
        SIZE_OF_LABEL_TYPE
      };

      /// Label names
      static const std::string NamesOfLabelType[SIZE_OF_LABEL_TYPE];

      /// Features
      typedef FeatureMap<> FeatureMapType;

      /// SharedPtr on feature map
      typedef boost::shared_ptr<FeatureMap<> > FeatureMapSharedPtrType;

      /// consensus features
      typedef ConsensusMap ConsensusMapType;

      /// SharedPtr on consensus features
      typedef boost::shared_ptr<ConsensusMap> ConsensusMapSharedPtrType;

      /// Main data type (experiment)
      typedef MSExperiment<> ExperimentType;

      /// SharedPtr on MSExperiment
      typedef boost::shared_ptr<ExperimentType> ExperimentSharedPtrType;

      //@}

      /// Default constructor
      LayerData()
        :	flags(),
          visible(true),
          flipped(false),
          type(DT_UNKNOWN),
          name(),
          filename(),
          peptides(),          
          param(),
          gradient(),
          filters(),
          annotations_1d(),
          modifiable(false),
          modified(false),
          label(L_NONE),
          features(new FeatureMapType()),
          consensus(new ConsensusMapType()),
          peaks(new ExperimentType()),
          current_spectrum_(0)
      {
        annotations_1d.resize(1);
      }

      /// Returns a const reference to the current spectrum (1d view)
      const ExperimentType::SpectrumType& getCurrentSpectrum() const;

      /// Returns a const reference to the current feature data
      const FeatureMapSharedPtrType& getFeatureMap() const
      {
        return features;
      }

      /// Returns a const reference to the current feature data
      FeatureMapSharedPtrType& getFeatureMap()
      {
        return features;
      }

      /// Returns a const reference to the consensus feature data
      const ConsensusMapSharedPtrType& getConsensusMap() const
      {
        return consensus;
      }

      /// Returns current consensus map (mutable)
      ConsensusMapSharedPtrType& getConsensusMap()
      {
        return consensus;
      }

      /// Returns a const reference to the current peak data
      const ExperimentSharedPtrType& getPeakData() const
      {
        return peaks;
      }

      /// Returns a mutable reference to the current peak data
      ExperimentSharedPtrType& getPeakData()
      {
        return peaks;
      }

      /// Returns a const reference to the annotations of the current spectrum (1d view)
      const Annotations1DContainer& getCurrentAnnotations() const
      {
        return annotations_1d[current_spectrum_];
      }

      /// Returns a mutable reference to the annotations of the current spectrum (1d view)
      Annotations1DContainer& getCurrentAnnotations()
      {
        return annotations_1d[current_spectrum_];
      }

      /// Returns a const reference to the annotations of the current spectrum (1d view)
      const Annotations1DContainer& getAnnotations(Size spectrum_index) const
      {
        return annotations_1d[spectrum_index];
      }

      /// Returns a mutable reference to the annotations of the current spectrum (1d view)
      Annotations1DContainer& getAnnotations(Size spectrum_index)
      {
        return annotations_1d[spectrum_index];
      }

      /// Returns a mutable reference to the current spectrum (1d view)
      ExperimentType::SpectrumType& getCurrentSpectrum()
      {        
        return (*peaks)[current_spectrum_];
      }

      /// Get the index of the current spectrum
      Size getCurrentSpectrumIndex() const
      {
        return current_spectrum_;
      }

      /// Set the index of the current spectrum
      void setCurrentSpectrumIndex(Size index)
      {
        current_spectrum_ = index;
      }

      /// if this layer is visible
      bool visible;

      /// if this layer is flipped (1d mirror view)
      bool flipped;

      /// data type (peak or feature data)
      DataType type;

      /// layer name
      String name;

      /// file name of the file the data comes from (if available)
      String filename;

      /// peptide identifications
      std::vector<PeptideIdentification> peptides;

      /// Layer parameters
      Param param;

      /// Gradient for 2D and 3D views
      MultiGradient gradient;

      /// Filters to apply before painting
      DataFilters filters;

      /// Annotations of all spectra of the experiment (1D view)
      std::vector<Annotations1DContainer> annotations_1d;

      /// Flag that indicates if the layer data can be modified (so far used for features only)
      bool modifiable;

      /// Flag that indicates that the layer data was modified since loading it
      bool modified;

      /// Label type
      LabelType label;

    private:
      /// feature data
      FeatureMapSharedPtrType features;

      /// consensus feature data
      ConsensusMapSharedPtrType consensus;

      /// peak data
      ExperimentSharedPtrType peaks;

      /// Index of the current spectrum
      Size current_spectrum_;
	};

  /// Print the contents to a stream.
	OPENMS_GUI_DLLAPI std::ostream& operator << (std::ostream& os, const LayerData& rhs);

} //namespace

#endif
