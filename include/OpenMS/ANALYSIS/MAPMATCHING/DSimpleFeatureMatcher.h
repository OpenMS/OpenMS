// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//									 OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//	Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//	This library is free software; you can redistribute it and/or
//	modify it under the terms of the GNU Lesser General Public
//	License as published by the Free Software Foundation; either
//	version 2.1 of the License, or (at your option) any later version.
//
//	This library is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
//	Lesser General Public License for more details.
//
//	You should have received a copy of the GNU Lesser General Public
//	License along with this library; if not, write to the Free Software
//	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA	02111-1307	USA
//
// --------------------------------------------------------------------------
// $Id: DSimpleFeatureMatcher.h,v 1.5 2006/05/02 09:39:07 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DSIMPLEFEATUREMATCHER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DSIMPLEFEATUREMATCHER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DBaseFeatureMatcher.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#if defined OPENMS_DEBUG
#if ! defined V_DSimpleFeatureMatcher
#define V_DSimpleFeatureMatcher(bla) std::cout << bla << std::endl;
#endif
#else
#define V_DSimpleFeatureMatcher(bla)
#endif

namespace OpenMS
{
	
	/**
		 @brief Feature matcher that uses a simple approach.  For each feature, we
		 find a best companion in the other map.  If both point at each other, and
		 the similarity is high enough, they are paired.
	**/
	template <Size D, typename TraitsT = KernelTraits, typename FeatureT = DFeature<D,TraitsT> >
	class DSimpleFeatureMatcher
		: public DBaseFeatureMatcher < D, TraitsT, FeatureT >
	{
	 protected:
		typedef DBaseFeatureMatcher < D, TraitsT, FeatureT > Base;

	 public:

		/** @name Type definitions
		 */
		//@{
		/// Coordinates
		typedef DimensionDescription<DimensionDescriptionTagLCMS> DimensionDescriptionType;
		enum DimensionId
			{
				RT = DimensionDescriptionType::RT,
				MZ = DimensionDescriptionType::MZ
			};
		enum { DIMENSION = D };

		// The base knows it all...
		typedef typename Base::TraitsType             TraitsType;
		typedef typename Base::QualityType            QualityType;
		typedef typename Base::PositionType           PositionType;
		typedef typename Base::IntensityType          IntensityType;
		typedef typename Base::FeatureType            FeatureType;
		typedef typename Base::FeatureMapType         FeatureMapType;
		typedef typename Base::FeaturePairType        FeaturePairType;
		typedef typename Base::FeaturePairVectorType  FeaturePairVectorType;
		typedef typename Base::GridType               GridType;
		//@}

		using Base::setParam;
		using Base::getParam;
		using Base::setGrid;
		using Base::getGrid;
		using Base::setFeatureMap;
		using Base::getFeatureMap;
		using Base::setFeaturePairs;
		using Base::getFeaturePairs;

		//------------------------------------------------------------
		
		///@name Constructors, destructor and assignment
		//@{
		/// Constructor
		DSimpleFeatureMatcher()
			: Base()
		{}
		
		/// Copy constructor
		DSimpleFeatureMatcher(const DSimpleFeatureMatcher& source)
			: Base(source)
		{}
		
		///	 Assignment operator
		DSimpleFeatureMatcher& operator = (DSimpleFeatureMatcher source)
		{
			Base::operator=(source);
			return *this;
		} 
		
		/// Destructor
		~DSimpleFeatureMatcher() {}
		//@}
		
		//------------------------------------------------------------

		/** @name Accesssor methods
		 */
		//@{	
		// inherited
		//@} 
		
		//----------------------------------------------------------------------
		
		/**@brief Runs the algorithm.  Assigns result to what is accessible by
			  getFeaturePairs() and getGrid().
		*/
		virtual void run()
		{
			parseParam_();
			findFeaturePairs_();
			assignTrivialGrid_();
		}

	 protected:

		//----------------------------------------------------------------------
		
		/** @name Methods
		 */
		//@{

		/// Parses the parameters, assigns their values to instance members.
		void parseParam_()
		{
#define V_parseParam_(bla) V_DSimpleFeatureMatcher(bla)
			V_parseParam_("@@@ parseParam_()");

			{
				std::string param_name_prefix = "similarity:diff_exponent:";
				for ( Size dimension = 0; dimension < 2; ++dimension)
				{
					std::string param_name =
						param_name_prefix + DimensionDescriptionType::dimension_name_short[dimension];
					DataValue data_value = getParam().getValue(param_name);
					if ( data_value == DataValue::EMPTY )
					{
						throw Exception::ElementNotFound<std::string>
							(__FILE__,__LINE__,__PRETTY_FUNCTION__,param_name);
					}
					else
					{
						diff_exponent_[dimension] = data_value;
						V_parseParam_(param_name<< ": "<<diff_exponent_[dimension]);
					}
				}
			}

			{
				std::string param_name_prefix = "similarity:diff_intercept:";
				for ( Size dimension = 0; dimension < 2; ++dimension)
				{
					std::string param_name =
						param_name_prefix + DimensionDescriptionType::dimension_name_short[dimension];
					DataValue data_value = getParam().getValue(param_name);
					if ( data_value == DataValue::EMPTY )
					{
						throw Exception::ElementNotFound<std::string>
							(__FILE__,__LINE__,__PRETTY_FUNCTION__,param_name);
					}
					else
					{
						diff_intercept_[dimension] = data_value;
						V_parseParam_(param_name<< ": "<<diff_intercept_[dimension]);
					}
				}
			}

			{
				std::string param_name = "similarity:pair_min_quality";
				DataValue data_value = getParam().getValue(param_name);
				if ( data_value == DataValue::EMPTY )
				{
					throw Exception::ElementNotFound<std::string>
						(__FILE__,__LINE__,__PRETTY_FUNCTION__,param_name);
				}
				else
				{
					pair_min_quality_ = data_value;
					V_parseParam_(param_name<< ": "<<pair_min_quality_);
				}
			}

#undef V_parseParam_
		} // parseParam_


		//----------------------------------------------------------------------

		/// The actual algorithm for finding feature pairs.
		void findFeaturePairs_()
		{
#define V_findFeaturePairs_(bla) V_DSimpleFeatureMatcher(bla)
			V_findFeaturePairs_("@@@ findFeaturePairs_()");


			// For each feature in map 0, find his/her best friend in map 1
			std::vector<Size>        best_companion_index_0(getFeatureMap(0).size(),Size(-1));
			std::vector<QualityType> best_companion_quality_0(getFeatureMap(0).size(),0);
			for ( Size fi0 = 0; fi0 < getFeatureMap(0).size(); ++fi0 )
			{
				QualityType best_quality = -std::numeric_limits<QualityType>::max();
				for ( Size fi1 = 0; fi1 < getFeatureMap(1).size(); ++ fi1 )
				{
					QualityType quality = similarity_( getFeatureMap(0)[fi0], getFeatureMap(1)[fi1] );
					if ( quality > best_quality )
					{
						best_quality = quality;
						best_companion_index_0[fi0] = fi1;
					}
				}
				best_companion_quality_0[fi0] = best_quality;
			}

			// For each feature in map 1, find his/her best friend in map 0
			std::vector<Size>        best_companion_index_1(getFeatureMap(1).size(),Size(-1));
			std::vector<QualityType> best_companion_quality_1(getFeatureMap(1).size(),0);
			for ( Size fi1 = 0; fi1 < getFeatureMap(1).size(); ++fi1 )
			{
				QualityType best_quality = -std::numeric_limits<QualityType>::max();
				for ( Size fi0 = 0; fi0 < getFeatureMap(0).size(); ++ fi0 )
				{
					QualityType quality = similarity_( getFeatureMap(1)[fi1], getFeatureMap(0)[fi0] );
					if ( quality > best_quality )
					{
						best_quality = quality;
						best_companion_index_1[fi1] = fi0;
					}
				}
				best_companion_quality_1[fi1] = best_quality;
			}

			// And if both like each other, they become a pair.
			for ( Size fi0 = 0; fi0 < getFeatureMap(0).size(); ++fi0 )
			{
				// fi0 likes someone ...
				if ( best_companion_quality_0[fi0] > pair_min_quality_ )
				{
					// ... who likes him too ...
					Size best_companion_of_fi0 = best_companion_index_0[fi0];
					if ( best_companion_index_1[best_companion_of_fi0] == fi0 &&
							 best_companion_quality_1[best_companion_of_fi0] > pair_min_quality_
						 )
					{
						getFeaturePairs().push_back
							( FeaturePairType ( getFeatureMap(0)[fi0],
																	getFeatureMap(1)[best_companion_of_fi0],
																	best_companion_quality_0[fi0] + best_companion_quality_1[best_companion_of_fi0]
																)
							);
					}
				}
			}

			// this->dumpFeaturePairs();

#undef V_findFeaturePairs_
		} // findFeaturePairs_

		//----------------------------------------------------------------------

		/**@brief Compute the similarity for a pair of features; larger quality
			 values are better.
			 
			 The returned value should express our confidence that one feature might
			 possibly be matched to the other.

			 The details here are kind of alchemy ...
		*/
		QualityType similarity_ ( FeatureType const & left, FeatureType const & right ) const 
		{
			QualityType right_intensity(right.getIntensity());
			if ( right_intensity == 0 ) return 0;
			QualityType intensity_ratio = left.getIntensity() / right_intensity;
			if ( intensity_ratio > 1. ) intensity_ratio = 1. / intensity_ratio;

			PositionType position_difference = left.getPosition() - right.getPosition();
			for ( Size dimension = 0; dimension < 2; ++dimension )
			{
				// Take the absolute value
				if ( position_difference[dimension] < 0 )
				{
					position_difference[dimension] = -position_difference[dimension];
				}
				// Raise the difference to a (potentially fractional) power
				position_difference[dimension] =
					pow(position_difference[dimension],diff_exponent_[dimension]);
				// Add an absolute number
				position_difference[dimension] += diff_intercept_[dimension];
			}

			return intensity_ratio / position_difference[RT] / position_difference[MZ];
		}

		//----------------------------------------------------------------------

		/**@brief This will assign a trivial grid with only one grid cell whose
			 range is infinite.

			 @internal Alternatively, we could compute a bounding box and assign it
			 to the range of the grid cell.  But that is NEVER correct, since ranges
			 are half-open.  So, how much should we enlarge it, ... ???

		 */
		void assignTrivialGrid_()
		{
#define V_assignTrivialGrid_(bla) V_DSimpleFeatureMatcher(bla)
			V_assignTrivialGrid_("@@@ assignTrivialGrid_()");

			getGrid().clear();
			getGrid().resize(1);
			getGrid()[0].setMin(PositionType::min_negative);
			getGrid()[0].setMax(PositionType::max);
			return;

#undef V_assignTrivialGrid_
		}

		//@} // Methods

		//------------------------------------------------------------

		/** @name Data members
		 */
		//@{	

		/// A parameter for similarity_().
		QualityType diff_exponent_[2];

		/// A parameter for similarity_().
		QualityType diff_intercept_[2];

		/// A parameter for findFeaturePairs_().
		QualityType pair_min_quality_;

		//@}
		
	}; // DSimpleFeatureMatcher
	
} // namespace OpenMS

#endif	// OPENMS_ANALYSIS_MAPMATCHING_DSIMPLEFEATUREMATCHER_H
