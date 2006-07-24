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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DGEOMHASHSHIFTFEATUREMATCHER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DGEOMHASHSHIFTFEATUREMATCHER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DBaseFeatureMatcher.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <fstream>

#if defined OPENMS_DEBUG && ! defined V_DGeomHashShiftFeatureMatcher
#define V_DGeomHashShiftFeatureMatcher(bla) std::cout << bla << std::endl;
#else
#define V_DGeomHashShiftFeatureMatcher(bla)
#endif

namespace OpenMS
{
	
	/**
		 @brief Feature matcher that uses geometric hashing to find a good shift

		 While the feature positions can have D dimensions, only the first two are
		 used to find a shift.
	**/
	template <Size D, typename TraitsT = KernelTraits, typename FeatureT = DFeature<D,TraitsT> >
	class DGeomHashShiftFeatureMatcher
		: public DBaseFeatureMatcher < D, TraitsT, FeatureT >
	{
	 protected:
		typedef DBaseFeatureMatcher < D, TraitsT, FeatureT > Base;

		// We need this to make the intensity bounding box use the intensity type
		// instead of the coordinate type.
		struct IntensityBoundingBoxTraits : TraitsT
		{
			typedef typename TraitsT::IntensityType CoordinateType;
		};

	 public:
		
		/** @brief Nested class to represent a shift.

		    The shift itself is stored as a DPosition.  Also provided is a
		    quality value, with an acompanying comparator.
		*/
		class Shift
		{
		 public:

			Shift()
				: position_(0),
					quality_(0)
			{}

			Shift(Shift const & source)
				: position_(source.position_),
					quality_(source.quality_)
			{}

			Shift & operator= (Shift const & source)
			{
				position_ = source.position_;
				quality_ = source.quality_;
				return *this;
			}
			
			~Shift() {}
			
			typedef TraitsT TraitsType;
			typedef DPosition<D,TraitsType> PositionType;
			typedef typename TraitsType::QualityType QualityType;

			/// Non-mutable access to the data point position (multidimensional)
			const PositionType& getPosition() const { return position_; }
			/// Mutable access to the data point position (multidimensional)
			PositionType& getPosition() { return position_; }
			/// Mutable access to the data point position (multidimensional)
			void setPosition(const PositionType& position) { position_ = position; }

			/// Non-mutable access to the quality
			const QualityType& getQuality() const { return quality_; }
			/// Mutable access to the quality
			QualityType& getQuality() { return quality_; }
			/// Mutable access to the quality
			void setQuality(const QualityType& quality) { quality_ = quality; }

			/*

			not used currently

			/// Compare by getQuality()
			struct QualityLess : std::binary_function < Shift, Shift, bool >
			{
				bool operator () ( Shift const & left, Shift const & right ) const {return ( left.getQuality() < right.getQuality() );}
				bool operator () ( Shift const & left, QualityType const & right ) const {return ( left.getQuality() < right );}
				bool operator () ( QualityType const & left, Shift const & right ) const {return ( left< right.getQuality() );}
				bool operator () ( QualityType const & left, QualityType const & right ) const {return ( left < right );}
			};
			*/

		 protected:
			PositionType position_;
			QualityType quality_;
		};

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
		typedef typename Base::TraitsType TraitsType;
		typedef typename Base::QualityType QualityType;
		typedef typename Base::PositionType PositionType;
		typedef typename Base::IntensityType IntensityType;
		typedef typename Base::GridType GridType;
		typedef typename Base::FeatureType FeatureType;
		typedef typename Base::FeatureMapType FeatureMapType;
		typedef typename Base::FeaturePairType FeaturePairType;
		typedef typename Base::FeaturePairVectorType FeaturePairVectorType;
		typedef DBoundingBox<DIMENSION,TraitsType>  PositionBoundingBoxType;
		typedef DBoundingBox<1,IntensityBoundingBoxTraits> IntensityBoundingBoxType;
		typedef std::vector <Size> FeatureBucketType;
		typedef Matrix < FeatureBucketType > FeatureBucketMatrixType;
		typedef Shift ShiftType;
		typedef Matrix < typename ShiftType::QualityType > ShiftQualityMatrixType;
		typedef Matrix < ShiftType > ShiftMatrixType;
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
		DGeomHashShiftFeatureMatcher()
			: Base()
		{}
		
		/// Copy constructor
		DGeomHashShiftFeatureMatcher(const DGeomHashShiftFeatureMatcher& source)
			: Base(source)
		{}
		
		///	 Assignment operator
		DGeomHashShiftFeatureMatcher& operator = (DGeomHashShiftFeatureMatcher source)
		{
			Base::operator=(source);
			return *this;
		} 
		
		/// Destructor
		~DGeomHashShiftFeatureMatcher() {}
		//@}
		
		//------------------------------------------------------------

		/** @name Accesssor methods
		 */
		//@{	
		// inherited
		//@} 
		
		//----------------------------------------------------------------------
		
		/// estimates the transformation for each grid cell
		virtual void run()
		{
			parseParam_();
			computeFeatureBuckets_();
			computeShiftBuckets_();
			// computeShift_();
			computeFinalShift_();

			this->getFeaturePairs().resize(1);
		}
		
	 protected:

		//----------------------------------------------------------------------
		
		/** @name Methods
		 */
		//@{


		void parseParam_()
		{
#define V_parseParam_(bla) V_DGeomHashShiftFeatureMatcher(bla)
			V_parseParam_("@@@ parseParam_()");

			// Initialize feature_bucket_size_ with values from param_.
			std::string fm_bs = "feature_map:bucket_size:";
			for ( Size dimension = 0; dimension < 2; ++dimension)
			{
				std::string fm_bs_dn = fm_bs + DimensionDescriptionType::dimension_name_short[dimension];
				DataValue data_value = getParam().getValue(fm_bs_dn);
				if ( data_value == DataValue::EMPTY )
				{
					throw Exception::ElementNotFound<std::string>
						(__FILE__,__LINE__,__PRETTY_FUNCTION__,fm_bs_dn);
				}
				else
				{
					feature_bucket_size_[dimension] = data_value;
					V_parseParam_(fm_bs_dn<< ": "<<feature_bucket_size_[dimension]);
				}
			}

			// Initialize shift_bucket_size_ with values from param_.
			std::string tm_bs  = "shift_map:bucket_size:";
			for ( Size dimension = 0; dimension < 2; ++dimension)
			{
				String tm_bs_dn
					= tm_bs
					+ DimensionDescriptionType::dimension_name_short[dimension];
				DataValue data_value = getParam().getValue(tm_bs_dn);
				if ( data_value == DataValue::EMPTY )
				{
					throw Exception::ElementNotFound<std::string>
						(__FILE__,__LINE__,__PRETTY_FUNCTION__,tm_bs_dn);
				}
				else
				{
					shift_bucket_size_[dimension] = data_value;
					V_parseParam_(tm_bs_dn<< ": "<<shift_bucket_size_[dimension]);
				}
			}

			// Initialize feature_bucket_window_ with values from param_.
			std::string tm_fbw = "shift_map:feature_bucket_window:";
			for ( Size dimension = 0; dimension < 2; ++dimension)
			{
				std::string tm_fbw_dn
					= tm_fbw
					+ DimensionDescriptionType::dimension_name_short[dimension];
				DataValue data_value = getParam().getValue(tm_fbw_dn);
				if ( data_value == DataValue::EMPTY )
				{
					throw Exception::ElementNotFound<std::string>
						(__FILE__,__LINE__,__PRETTY_FUNCTION__,tm_fbw_dn);
				}
				else
				{
					feature_bucket_window_[dimension] = data_value;
					V_parseParam_(tm_fbw_dn<< ": "<<feature_bucket_window_[dimension]);
				}
			}

			// Initialize shift_bucket_window_ with values from param_.
			std::string const tm_tbw = "shift_map:shift_bucket_window:";
			for ( Size dimension = 0; dimension < 2; ++dimension)
			{
				std::string tm_tbw_dn
					= tm_tbw
					+ DimensionDescriptionType::dimension_name_short[dimension];
				DataValue data_value = getParam().getValue(tm_tbw_dn);
				if ( data_value == DataValue::EMPTY )
				{
					throw Exception::ElementNotFound<std::string>
						(__FILE__,__LINE__,__PRETTY_FUNCTION__,tm_tbw_dn);
				}
				else
				{
					shift_bucket_window_[dimension] = data_value;
					V_parseParam_(tm_tbw_dn<< ": "<<shift_bucket_window_[dimension]);
				}
			}

#undef V_parseParam_
		}


		/**@brief Fill the buckets with the indices of the corresponding features.
		 */
		void computeFeatureBuckets_()
		{
#define V_computeFeatureBuckets_(bla) V_DGeomHashShiftFeatureMatcher(bla)
			V_computeFeatureBuckets_("@@@ computeFeatureBuckets_()");

			// Shorthands ...
			PositionType & fbs = feature_bucket_size_;

			for ( Size map_index = 0; map_index < 2; ++map_index )
			{
				// Shorthands ...
				V_computeFeatureBuckets_("\n--- map_index: "<<map_index);
				FeatureMapType const     & fm     = getFeatureMap(map_index);
 				PositionBoundingBoxType  & fmpbb  = feature_map_position_bounding_box_[map_index] ;
				IntensityBoundingBoxType & fmibb  = feature_map_intensity_bounding_box_[map_index];

				fmpbb.clear();
				fmibb.clear();

				// Compute the bounding box for the feature map, with respect to
				// position and intensity.
				for ( typename FeatureMapType::ConstIterator fm_iter = fm.begin();
							fm_iter != fm.end();
							++fm_iter
						)
				{
					fmpbb.enlarge(fm_iter->getPosition());
					fmibb.enlarge(fm_iter->getIntensity());
				}
				V_computeFeatureBuckets_("fmpbb: "<<fmpbb<<"fmibb: "<<fmibb);
			}

			// Next we will enlarge each feature_map_position_bounding_box_ such
			// that all buckets will have the same diagonal.  To provide against
			// rounding errors, we allocate one bucket more than needed (in each
			// dimension) and shift the grid by one-half of the difference.
			for ( Size map_index = 0; map_index < 2; ++map_index )
			{
				// Shorthands ...
				V_computeFeatureBuckets_("\n--- map_index: "<<map_index);
				FeatureMapType          const & fm     = getFeatureMap(map_index);
 				PositionBoundingBoxType const & fmpbb  = feature_map_position_bounding_box_[map_index] ;
 				PositionBoundingBoxType       & fmpbbe = feature_map_position_bounding_box_enlarged_[map_index] ;
				FeatureBucketMatrixType       & fb     = feature_bucket_[map_index];

				// Compute num_buckets.  Compute extra margin to make bounding box a
				// multiple of feature buckets.
				PositionType const diagonal = fmpbb.diagonal();
				PositionType diagonal_enlarged;
				V_computeFeatureBuckets_("diagonal: " << diagonal);
				int num_buckets[2];
				for ( Size dimension = 0; dimension < DIMENSION; ++dimension)
				{
					num_buckets[dimension] = int(1.1 + diagonal[dimension]/fbs[dimension]);
					diagonal_enlarged[dimension] = fbs[dimension] * num_buckets[dimension];
				}
				V_computeFeatureBuckets_("num_buckets: "<<num_buckets[0]<<' '<<num_buckets[1]);
				V_computeFeatureBuckets_("diagonal_enlarged: "<<diagonal_enlarged);

				// The extra margin.
				PositionType extra_feature_bucket_size_(diagonal_enlarged-diagonal);
				extra_feature_bucket_size_ /= 2;
				V_computeFeatureBuckets_("efbs: " << extra_feature_bucket_size_);

				// Compute the enlarged feature map bounding box accordingly.
				fmpbbe.clear();
				fmpbbe.enlarge( fmpbb.min() - extra_feature_bucket_size_ );
				fmpbbe.enlarge( fmpbb.max() + extra_feature_bucket_size_ );
				V_computeFeatureBuckets_("fmpbbe: "<<fmpbbe);

				// Resize feature_bucket_[map_index] accordingly.
				fb.resize(num_buckets[0],num_buckets[1]);
				V_computeFeatureBuckets_("rows: "<<fb.rows()<<"  cols: "<<fb.cols());

				// Now, finally, we store the indices of the features in their
				// corresponding buckets.
				PositionType const & fmpbbe_min = fmpbbe.min();
				for ( Size index= 0; index < fm.size(); ++index )
				{
					PositionType position = fm[index].getPosition() - fmpbbe_min;
					fb ( Size(position[0]/fbs[0]), Size(position[1]/fbs[1]) ).push_back(index);
				}
				
				// Optionally, write debug output as specified in param.
				DataValue data_value_dump_feature_buckets = getParam().getValue("debug:dump_feature_buckets");
				if ( data_value_dump_feature_buckets != DataValue::EMPTY )
				{
					std::string   dump_filename = data_value_dump_feature_buckets;
					std::ofstream dump_file(dump_filename.c_str());
					V_computeFeatureBuckets_("### Writing "<<dump_filename);
					dump_file << "# " << dump_filename << " generated " << Date::now() << std::endl;
					dump_file << "# Positions of features in non-empty feature buckets" << std::endl;
					for ( FeatureBucketMatrixType::ConstIterator iter = fb.begin(); iter != fb.end(); ++iter)
					{
						if (iter->empty()) continue;
						std::pair<Size,Size> row_col = fb.indexPair(iter-fb.begin());
						dump_file << row_col.first << ' ' << row_col.second << " #bucket" << std::endl;
						for ( FeatureBucketType::const_iterator viter = iter->begin(); viter != iter->end(); ++viter)
						{ dump_file << fm[*viter].getPosition()[0] <<' '<<fm[*viter].getPosition()[1] << std::endl; }
						dump_file << std::endl;
					}
					dump_file << "# " << dump_filename << " EOF " << Date::now() << std::endl;
				}
			}
			
			return;
#undef V_computeFeatureBuckets_
		} // computeFeatureBuckets_


		//----------------------------------------------------------------------

		/**@brief Fill the buckets of shifts.

		   Note that computeFeatureBuckets_() must have been called before to make
			 this work properly.
		*/
		void computeShiftBuckets_()
		{
#define V_computeShiftBuckets_(bla) V_DGeomHashShiftFeatureMatcher(bla)
			V_computeShiftBuckets_("@@@ computeShiftBuckets_()");

			// Shorthands ...
			ShiftQualityMatrixType & tb      = shift_bucket_;
			PositionType                 & tbs     = shift_bucket_size_;
			PositionBoundingBoxType      & tbb     = shift_bounding_box_ ;
			PositionBoundingBoxType      & tbbe    = shift_bounding_box_enlarged_ ;
			Size                   const (&fbw)[2] = feature_bucket_window_;
			ShiftMatrixType        & tm      = shift_matrix_;
			
			// Compute the bounding box for the shift map
			{
				tbb.clear();
				PositionType fmpbbmindiff =
					feature_map_position_bounding_box_[1].min() -
					feature_map_position_bounding_box_[0].min();
				PositionType windowdiff =
					PositionType ( feature_bucket_size_[0]*(feature_bucket_window_[0]*2+1),
												 feature_bucket_size_[1]*(feature_bucket_window_[1]*2+1) );
				tbb.enlarge( fmpbbmindiff + windowdiff );
				tbb.enlarge( fmpbbmindiff - windowdiff );
			}
			V_computeShiftBuckets_("tbb: "<<tbb);

			// Next we will enlarge each bucket_size_ such that all buckets will
			// have the same diagonal.  To provide against rounding errors, we
			// allocate one bucket more than needed (in each dimension) and shift
			// the grid by one-half.

			PositionType half_of_shift_bucket_size_(tbs);
			half_of_shift_bucket_size_ /= 2;
			V_computeShiftBuckets_("hotbs: " << half_of_shift_bucket_size_);

			// Adjust the enlarged shift map bounding box accordingly.
			tbbe.enlarge( tbb.min() - half_of_shift_bucket_size_ );
			tbbe.enlarge( tbb.max() + half_of_shift_bucket_size_ );
			V_computeShiftBuckets_("tbbe: "<<tbbe);

			// Compute shift_bucket_size_ and num_buckets.
			PositionType diagonal = tbbe.diagonal();
			V_computeShiftBuckets_("diagonal: " << diagonal);
			int num_buckets[2];
			for ( Size dimension = 0; dimension < DIMENSION; ++dimension)
			{
				num_buckets[dimension] = int(diagonal[dimension]/tbs[dimension]);
				tbs[dimension] = diagonal[dimension] / num_buckets[dimension];
			}
			V_computeShiftBuckets_("tbs: "<<tbs);

			// Resize shift_bucket_ accordingly.
			tb.resize(num_buckets[0]+1,num_buckets[1]+1);
			V_computeShiftBuckets_("rows: "<<tb.rows()<<"  cols: "<<tb.cols());

			// Resize shift_matrix_ according to feature_bucket_[1]
			tm.resize(feature_bucket_[1].sizePair());

			// Now we store the shifts for all relevant feature pairs in their
			// corresponding buckets.  Each shift is distributed among its
			// four neighboring "buckets", with weights according to the distances
			// from these corner points.  Note that the outer two loops (over i and
			// j) enumerate the "image" (feature_bucket_[1]), then we search for
			// "pre-images" (feature_bucket_[0}) in the two inner loops (over k and
			// l).  (And of course, finally, we enumerate all feature pairs.)  This
			// way we can associate the shifts vectors to buckets of the
			// image, and when we will later apply it, we will not change the
			// pre-image, which might be a consensus or so.
#define V_computeShiftBuckets_enumeration(bla) V_computeShiftBuckets_(bla)
			PositionType const & tbbe_min = tbbe.min();
			for ( Size j = 0; j < feature_bucket_[1].cols(); ++j )
			{
				// Clear the shift buckets.
				std::fill(tb.begin(),tb.end(),QualityType(0));
				
				for ( Size i = 0; i < feature_bucket_[1].rows(); ++i )
				{
					for ( Size k = std::max ( int(i-fbw[0]), 0 );
								k     <= std::min ( i + fbw[0], feature_bucket_[0].rows()-1 );
								++k )
					{
						for ( Size l = std::max ( int(j-fbw[1]), 0 );
									l     <= std::min ( j + fbw[1], feature_bucket_[0].cols()-1 );
									++l )
						{
							// V_computeShiftBuckets_enumeration("ijkl: "<<i<<' '<<j<<' '<<k<<' '<<l);
							for ( FeatureBucketType::const_iterator fb0_iter = feature_bucket_[0](i,j).begin();
										fb0_iter != feature_bucket_[0](i,j).end();
										++fb0_iter
									)
							{
								for ( FeatureBucketType::const_iterator fb1_iter = feature_bucket_[1](k,l).begin();
											fb1_iter != feature_bucket_[1](k,l).end();
											++fb1_iter
										)
								{
									// Compute the shift corresponding to a pair of features.
									ShiftType shift = shift_( getFeatureMap(0)[*fb0_iter],
																															getFeatureMap(1)[*fb1_iter] );
									V_computeShiftBuckets_enumeration("shift: "<<shift.getPosition());

									PositionType tpwm = shift.getPosition();
									tpwm -= tbbe_min;
									V_computeShiftBuckets_enumeration("trans.pos wrt tbbe_min: "<<shift.getPosition());

									QualityType  const & tq = shift.getQuality();

									// Compute the bucket index (the lowest of the four) for
									// this shift.  Also compute the fractional part of
									// the position within the bucket.
									Size bucket_index[2];
									PositionType bucket_fraction;
									for ( Size dimension = 0; dimension < 2; ++dimension )
									{
										bucket_fraction[dimension] = tpwm[dimension] / tbs[dimension];  // floating point division
										bucket_index[dimension]    = (Size) bucket_fraction[dimension]; // round down (yes we are >= 0)
										bucket_fraction[dimension] -= bucket_index[dimension];          // fractional part
									}
									PositionType bucket_fraction_complement(1,1);
									bucket_fraction_complement -= bucket_fraction;

									// Distribute the quality of the shift among the four neighboring buckets.
									QualityType factor;

									factor = bucket_fraction_complement[0] * bucket_fraction_complement[1];
									tb( bucket_index[0], bucket_index[1] ) += tq * factor;

									factor = bucket_fraction_complement[0] * bucket_fraction[1];
									tb( bucket_index[0], bucket_index[1] + 1 ) += tq * factor;

									factor = bucket_fraction[0] * bucket_fraction_complement[1];
									tb( bucket_index[0] + 1, bucket_index[1] ) += tq * factor;
	
									factor = bucket_fraction[0] * bucket_fraction[1];
									tb( bucket_index[0] + 1, bucket_index[1] + 1 ) += tq * factor;

								} // for fb1_iter
							} // for fb0_iter
						} // for l
					} // for k
				} // for i
				
				// Compute shift for this columns of image buckets.
				ShiftType result = computeShift_();
				shift_matrix_(0,j) = result; // oh, a weird HACK at the moment!
				
				V_computeShiftBuckets_enumeration(/*i<<' '<<*/j<<' '<<result.getPosition()<<' '<<result.getQuality()<< " #result");
					
			} // for j
#undef V_computeShiftBuckets_enumeration
			
			// Optionally, write debug output as specified in param.
			DataValue data_value_dump_shift_buckets = getParam().getValue("debug:dump_shift_buckets");
			if ( data_value_dump_shift_buckets != DataValue::EMPTY )
			{
				std::string   dump_filename = data_value_dump_shift_buckets;
				std::ofstream dump_file(dump_filename.c_str());
				V_computeShiftBuckets_("### Writing "<<dump_filename);
				dump_file << "# " << dump_filename << " generated " << Date::now() << std::endl;
				dump_file << "# Shift buckets: xcoord ycoord quality xindex yindex" << std::endl;

				for ( typename ShiftQualityMatrixType::ConstIterator iter = tb.begin(); iter != tb.end(); ++iter)
				{
					std::pair<Size,Size> row_col = tb.indexPair(iter-tb.begin());
					dump_file << tbbe_min[0] + tbs[0] * row_col.first << ' '
										<< tbbe_min[1] + tbs[1] * row_col.second << ' '
										<< *iter << ' '
										<< row_col.first << ' '
										<< row_col.second
										<< " #tb" << std::endl ;
				}
				dump_file << "# " << dump_filename << " EOF " << Date::now() << std::endl;
			}
			
#undef V_computeShiftBuckets_
		} // computeShiftBuckets_


		//----------------------------------------------------------------------

		/**@brief Compute the final shift.
		 */
		void computeFinalShift_()
		{
#define V_computeFinalShift_(bla) V_DGeomHashShiftFeatureMatcher(bla)
			V_computeFinalShift_("@@@ computeFinalShift_()");

			ShiftType final_shift_;
			for ( typename ShiftMatrixType::ConstIterator iter = shift_matrix_.begin();
						iter != shift_matrix_.end();
						++iter
					)
			{
				typename ShiftType::QualityType quality = iter->getQuality();
				final_shift_.getPosition() += iter->getPosition() * quality;
				final_shift_.getQuality()  += quality;
			}
			final_shift_.getPosition() /= final_shift_.getQuality();
	
			V_computeFinalShift_("final_shift_: "<<final_shift_.getPosition()<<' '<<final_shift_.getQuality());

			Size const & tm_rows = shift_matrix_.rows();
			Size const & tm_cols = shift_matrix_.cols();
			
			std::vector< ShiftType > shift_by_row(tm_rows);
			std::vector< ShiftType > shift_by_col(tm_cols);
			for ( Size row = 0; row < tm_rows; ++row )
			{
				for ( Size col = 0; col < tm_cols; ++col )
				{
					ShiftType const & shift = shift_matrix_(row,col);
					shift_by_row[row].getPosition() += shift.getPosition() * shift.getQuality();
					shift_by_row[row].getQuality() += shift.getQuality();
					shift_by_col[col].getPosition() += shift.getPosition() * shift.getQuality();
					shift_by_col[col].getQuality() += shift.getQuality();
				}
			}

			for ( Size row = 0; row < tm_rows; ++row )
			{
				if ( shift_by_row[row].getQuality() )
					shift_by_row[row].getPosition() /= shift_by_row[row].getQuality();
				V_computeFinalShift_(row<<' '<<shift_by_row[row].getPosition()<<' '<<shift_by_row[row].getQuality()<<" #tr_by_row");
			}

			for ( Size col = 0; col < tm_cols; ++col )
			{
				if ( shift_by_col[col].getQuality() )
					shift_by_col[col].getPosition() /= shift_by_col[col].getQuality();
				V_computeFinalShift_(col<<' '<<shift_by_col[col].getPosition()<<' '<<shift_by_col[col].getQuality()<<" #tr_by_col");
			}

#undef V_computeFinalShift_
		}

		//----------------------------------------------------------------------

		/**@brief Compute the shift.

		    Note that shift_buckets_ must have been calculated before.
		*/
		ShiftType computeShift_()
		{
#define V_computeShift_(bla) // V_DGeomHashShiftFeatureMatcher(bla)
			V_computeShift_("@@@ computeShift_()");

			// Shorthands ...
			ShiftQualityMatrixType const & tb = shift_bucket_;
			PositionType const & tbs = shift_bucket_size_;
			Size const (&tbw)[2] = shift_bucket_window_;

			// Find the transformation bucket with highest impact (quality).
			Size tb_max_element_index = std::max_element(tb.begin(),tb.end()) - tb.begin();
			Size tb_max_indices[2];
			tb_max_indices[0] = tb.rowIndex(tb_max_element_index);
			tb_max_indices[1] = tb.colIndex(tb_max_element_index);
			V_computeShift_("tb_max: "<<tb_max_indices[0]<<' '<<tb_max_indices[1]<<" quality="<<tb(tb_max_indices[0],tb_max_indices[1]));

			// Compute a weighted average of the shifts nearby the tb_max_element.
			ShiftType result; // initially zero
			PositionType const & tbbe_min = shift_bounding_box_enlarged_.min();
			int tb_run_indices[2];
			for ( tb_run_indices[0]  = std::max ( int (tb_max_indices[0] - tbw[0]), 0 );
						tb_run_indices[0] <= std::min ( int (tb_max_indices[0] + tbw[0]), int (tb.rows()) - 1 );
						++tb_run_indices[0]
					)
			{
				for ( tb_run_indices[1]  = std::max ( int (tb_max_indices[1] - tbw[1]), 0 );
							tb_run_indices[1] <= std::min ( int (tb_max_indices[1] + tbw[1]), int (tb.cols()) - 1 );
							++tb_run_indices[1]
						)
				{
					PositionType contribution_position(tbs);
					for ( Size dimension = 0; dimension < 2; ++dimension)
					{ contribution_position[dimension] *= tb_run_indices[dimension]; }
					contribution_position += tbbe_min;
					QualityType contribution_quality = tb( tb_run_indices[0], tb_run_indices[1] );
					result.getQuality() += contribution_quality;
					contribution_position *= contribution_quality;
					result.getPosition() += contribution_position;
				}
			}
			if ( result.getQuality() != 0 )
			{ result.getPosition() /= result.getQuality(); }
			else
			{
				// result.getPosition() is irrelevant anyway
			}
			
			return result;
#undef V_computeShift_
		} // computeShift_

		//----------------------------------------------------------------------

		/**@brief Compute the shift and similarity for a pair of features;
			 larger quality values are better.
			 
			 The returned value should express our confidence that one feature might
			 possibly be matched to the other.

			 Currently this will just calculate the ratio of intensities, either
			 "left/right" or "right/left", such that a value between 0 and 1 is
			 returned.

			 @todo Take the quality of the features themselves into account, i.e. how good they fit to their model.
		*/
		ShiftType shift_ ( FeatureType const & left, FeatureType const & right ) const 
		{
			ShiftType shift;
			shift.setPosition(right.getPosition() - left.getPosition());
			if ( right.getIntensity() == 0 ) shift.setQuality(0);
			QualityType result = left.getIntensity() / right.getIntensity();
			shift.setQuality( result <= 1. ? result : 1. / result );
			return shift;
		}

		//@} // Methods

		//------------------------------------------------------------

		/** @name Data members
		 */
		//@{	

		/// Holds the bounding box of all input features.
		PositionBoundingBoxType  feature_map_position_bounding_box_[2];

		/// Holds the enlarged bounding box for all input features.  It is larger
		/// by about half of a bucket in all directions.
		PositionBoundingBoxType  feature_map_position_bounding_box_enlarged_[2];

		/// Holds a bounding box for the input feature intensities.
		IntensityBoundingBoxType feature_map_intensity_bounding_box_[2];
		
		/// Feature indices are stored in theses buckets.
		FeatureBucketMatrixType feature_bucket_[2];

		/// Diagonal size of each bucket in feature_index_bucket_.
		PositionType feature_bucket_size_;

		/// Shifts are stored (summed up) in these buckets.
		ShiftQualityMatrixType shift_bucket_;

		/// Holds a bounding box for all possible shift vectors.
		PositionBoundingBoxType shift_bounding_box_;

		/// Holds an enlarged bounding box for all shift vectors.  It is
		/// larger by about half of a bucket in all directions.
		PositionBoundingBoxType shift_bounding_box_enlarged_;

		/// Diagonal size of each bucket in shift_bucket_.
		PositionType shift_bucket_size_;

		/// Number of surrounding buckets of feature indices to be considered when
		/// computing shifts.
		Size feature_bucket_window_[2];

		/// Number of surrounding buckets of shift indices to be considered when
		/// computing shifts.
		Size shift_bucket_window_[2];

		/// Matrix of shifts associated with buckets of feature_map_[1].
		ShiftMatrixType shift_matrix_;

		//@}
		
	}; // DGeomHashShiftFeatureMatcher
	
} // namespace OpenMS

#endif	// OPENMS_ANALYSIS_MAPMATCHING_DGEOMHASHSHIFTFEATUREMATCHER_H
