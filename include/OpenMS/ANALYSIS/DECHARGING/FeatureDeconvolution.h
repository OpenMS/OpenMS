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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DECHARGING_FEATUREDECONVOLUTION_H
#define OPENMS_ANALYSIS_DECHARGING_FEATUREDECONVOLUTION_H

// OpenMS
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>

namespace OpenMS
{

	class Compomer;
	
  /** 
    @brief An algorithm to decharge features (i.e. as found by FeatureFinder).

    @htmlinclude OpenMS_FeatureDeconvolution.parameters	
    
    @ingroup Analysis
  */
  
  class OPENMS_DLLAPI FeatureDeconvolution : public DefaultParamHandler
  {
    public:
    
			enum CHARGEMODE {QFROMFEATURE=1,QHEURISTIC,QALL};
			
      typedef FeatureMap<> FeatureMapType;
      typedef Feature FeatureType;
      typedef DPosition<2> ClusterPointType;
			typedef FeatureMapType::FeatureType::CoordinateType CoordinateType;
			typedef ILPDCWrapper::PairsType PairsType;
			
      /** @name Constructors and Destructor s
      */
      //@{
      /// default constructor
      FeatureDeconvolution();
      
      /// Copy constructor
      FeatureDeconvolution(const FeatureDeconvolution& source);
      
      /// Assignment operator
      FeatureDeconvolution& operator=(const FeatureDeconvolution& source);        

      /// destructor
      virtual ~FeatureDeconvolution();
      //@}    

      /**
				 @brief Compute a zero-charge feature map from a set of charged features
				 
				 Find putative ChargePairs, then score them and hand over to ILP.

				 @param fm_in  Input feature-map
				 @param fm_out Output feature-map (sorted by position and augmented with user params)
				 @param cons_map   [out] Output of grouped features belonging to a charge group
				 @param cons_map_p [out] Output of paired features connected by an egde
			*/
      void compute(const FeatureMapType &fm_in, FeatureMapType &fm_out, ConsensusMap &cons_map, ConsensusMap &cons_map_p);

    protected:

      void updateMembers_();
      
			/**
				@brief 1-sided Compomer for a feature
				
				Holds information on an explicit (with H+) 1-sided Compomer of a feature.
			**/
			struct CmpInfo_;

			/// test if "simple" edges have alternative
			/// (more difficult explanation) supported by neighbouring edges
			/// e.g. (.)   -> (H+) might be augmented to
			///      (Na+) -> (H+Na+)
			void inferMoreEdges_(PairsType& edges, Map<Size, std::set<CmpInfo_> >& feature_adducts);

			/// A function mostly for debugging 
			void printEdgesOfConnectedFeatures_(Size idx_1, Size idx_2, const PairsType& feature_relation);
			
			/**
				@brief returns true if the intensity filter was passed or switched off
				
				Filter for adding an edge only when the two features connected by it, fulfil the
				intensity criterion.
				
			**/
			inline bool intensityFilterPassed_(const Int q1, const Int q2, const Compomer& cmp,const FeatureType& f1,const FeatureType& f2);

			/**
				@brief determines if we should test a putative feature charge
				
				Answer query given the internal status of @em q_try.
				Features with q<=0 always return true.
				
			**/
			bool chargeTestworthy_(const Int feature_charge, const Int putative_charge, const bool other_unchanged) const;

      /// List of adducts used to explain mass differences
      MassExplainer::AdductsType potential_adducts_;
      /// labeling table
      Map<Size, String> map_label_;
      /// labeling table inverse
      Map<String, Size> map_label_inverse_;
			/// status of intensity filter for edges
			bool enable_intensity_filter_;
			/// status of charge discovery
			CHARGEMODE q_try_;

  };
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_DECHARGING_FEATUREDECONVOLUTION_H

