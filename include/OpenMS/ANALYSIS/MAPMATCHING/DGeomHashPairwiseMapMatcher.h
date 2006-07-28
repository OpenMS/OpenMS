// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DPAIRWISEMAPMATCHER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DPAIRWISEMAPMATCHER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DBasePairwiseMapMatcher.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGeomHashShiftSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DSimplePairFinder.h>
#include <OpenMS/KERNEL/DPeakConstReferenceArray.h>


#if defined OPENMS_DEBUG && ! defined DGeomHashPairwiseMapMatcher
#define V_DGeomHashPairwiseMapMatcher(bla) // std::cout << bla << std::endl;
#else
#define V_DGeomHashPairwiseMapMatcher(bla)
#endif


namespace OpenMS
{

  /**
     @brief This class represents a feature matching algorithm.

     It works on two feature maps and computes a vector of corresponding features
     in both maps (given by a feature pairs vector). 
     Therefore it first estimates the shift, which maps one feature map
     onto the other. This is done by the DGeomHashShiftSuperimposer class. 
     Afterwards using the transformation the corresponding features in both maps 
     can be determined using the DSimplePairFinder class.
     
     If a piecewise shift is assumed, the user can define a grid by setting 
     the number of buckets in the RT as well as the MZ dimension.     
     
  **/
  template <Size D, typename Traits = KernelTraits, typename MapT = DFeatureMap< D, DFeature< D, Traits > > >
  class DGeomHashPairwiseMapMatcher : public DBasePairwiseMapMatcher<D,Traits,MapT>
  {
  public:
    typedef DimensionDescription<DimensionDescriptionTagLCMS> DimensionDescriptionType;
    /// Defines the coordinates of peaks / features.
    enum DimensionId
    {
      RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
      MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };
    enum { DIMENSION = D };

    /** @name Type definitions
     */
    //@{
    /// Base
    typedef DBasePairwiseMapMatcher<D,Traits,MapT> Base;

    /// Traits
    typedef Traits TraitsType;

    /// Coordinate
    typedef typename TraitsType::CoordinateType CoordinateType;

    /// Quality
    typedef typename TraitsType::QualityType QualityType;

    /// Position
    typedef DPosition < DIMENSION, TraitsType > PositionType;

    //// Intensity
    typedef typename TraitsType::IntensityType IntensityType;

    /// Container for input features
    typedef MapT FeatureMapType;
    /// Type of features considered here
    typedef typename FeatureMapType::value_type FeatureType;

    /// Type of feature pairs
    typedef DFeaturePair < DIMENSION, FeatureType > FeaturePairType;
    /// Container for generated feature pairs
    typedef DFeaturePairVector < DIMENSION, FeatureType > FeaturePairVectorType;

    typedef DBoundingBox<DIMENSION,TraitsType>  PositionBoundingBoxType;

    typedef DLinearMapping<1, TraitsType> ShiftType;
    typedef DPeakConstReferenceArray<DFeatureMap<2> > PeakConstReferenceMapType;

    using Base::param_;
    using Base::feature_map_;
    using Base::grid_;
    using Base::all_feature_pairs_;
    //@}


    ///@name Constructors, destructor and assignment
    //@{
    /// Constructor
    DGeomHashPairwiseMapMatcher()
        : Base()
    {}

    /// Copy constructor
    DGeomHashPairwiseMapMatcher(const DGeomHashPairwiseMapMatcher& source)
        : Base(source),
        superimposer_(source.superimposer_),
        pair_finder_(source.pair_finder_),
        bounding_box_scene_map_(source.bounding_box_scene_map_),
        box_size_(source.box_size_),
        number_buckets_(source.number_buckets_)
    {}

    ///  Assignment operator
    virtual DGeomHashPairwiseMapMatcher& operator= (DGeomHashPairwiseMapMatcher source)
    {
      param_ = source.param_;
      feature_map_[0] = source.feature_map_[0];
      feature_map_[1] = source.feature_map_[1];
      grid_= source.grid_;
      bounding_box_scene_map_ = source.bounding_box_scene_map_;
      box_size_ = source.box_size_;
      number_buckets_ = source.number_buckets_;
      return *this;
    }

    /// Destructor
    virtual ~DGeomHashPairwiseMapMatcher() {}
    //@}

    /** @name Accesssor methods
     */
    //@{
    /// TODO get and set feature pairs of a given grid cell?, and all feature pairs
    //@}

    /// Estimates the transformation for each grid cell
    virtual void run()
    {
      superimposer_.setParam(param_);
      pair_finder_.setParam(param_);
      parseParam_();

      // compute the bounding boxes of the grid cells and initialize the grid transformation
      initGridTransformation_();

      // assign each feature of the scene map to the grid cells and build a pointer map for each grid cell
      PeakConstReferenceMapType scene_pointer_map(feature_map_[1]->begin(), feature_map_[1]->end());
      Size number_grid_cells = grid_.size();
      std::vector<PeakConstReferenceMapType> scene_grid_maps(number_grid_cells);
      buildGrid_(scene_pointer_map,scene_grid_maps);

      // initialize a pointer map with the features of the first (model or reference) map
      PeakConstReferenceMapType model_pointer_map(feature_map_[0]->begin(), feature_map_[0]->end());

      // compute the matching of each scene's grid cell features and all the features of the model map
      computeMatching_(model_pointer_map,scene_grid_maps);

      String all_feature_pairs_gnuplot_file =
				param_.getValue("debug:all_feature_pairs_gnuplot_file");
			if ( !all_feature_pairs_gnuplot_file.empty() )
			{
				this->dumpFeaturePairs(all_feature_pairs_gnuplot_file);
			}
    }

  protected:

    /** @name Data members
     */
    //@{
    /// This class computes the shift for the best mapping of one feature map to another
    DGeomHashShiftSuperimposer<PeakConstReferenceMapType> superimposer_;
    /// Given the shift, the pair_finder_ searches for corresponding features in the two maps
    DSimplePairFinder<D,Traits,PeakConstReferenceMapType> pair_finder_;
    /// Bounding box of the second map
    PositionBoundingBoxType bounding_box_scene_map_;
    /// Size of the grid cells
    PositionType box_size_;
    /// Number of buckets in each dimension
    PositionType number_buckets_;
    //@}

    /// compute the matching of each grid cell
    void computeMatching_(const PeakConstReferenceMapType& model_map, const std::vector<PeakConstReferenceMapType>& scene_grid_maps)
    {
#define V_computeMatching_(bla) V_DGeomHashPairwiseMapMatcher(bla)

      pair_finder_.setFeatureMap(0, model_map);

			// same model map for all superpositions
      superimposer_.setFeatureMap(0, model_map);

			String shift_buckets_file = param_.getValue("debug:shift_buckets_file");
			String feature_buckets_file = param_.getValue("debug:feature_buckets_file");

			// iterate over all grid cells of the scene map
      for (Size i = 0; i < scene_grid_maps.size(); ++i)
      {

				String algorithm = param_.getValue("algorithm");
				if ( algorithm == "geomhash_shift" )
				{
					V_computeMatching_("DGeomHashPairwiseMapMatcher:  algorithm \"geomhash_shift\", start superimposer");

					superimposer_.setFeatureMap(1, scene_grid_maps[i]);

					// optional debug output: buckets of the shift map
					if ( !shift_buckets_file.empty() )
					{
						superimposer_.getParam().setValue("debug:shift_buckets_file", shift_buckets_file+String(i).fillLeft('0',3));
					}

					// optional debug output: buckets of the feature maps
					if ( !feature_buckets_file.empty() )
					{
						superimposer_.getParam().setValue("debug:feature_buckets_file", feature_buckets_file+String(i).fillLeft('0',3));
					}

					superimposer_.run();

					// ???? copy array to vector -- but the old Grid class will be replaced anyway
					grid_[i].getMappings().resize(2,0);
					for ( Size dim = 0; dim < 2; ++dim )
					{
						ShiftType const& shift = superimposer_.getShift(dim);
						if ( !grid_[i].getMappings()[dim] )
						{
							grid_[i].getMappings()[dim] = new ShiftType;
						}
						*grid_[i].getMappings()[dim] = shift;
						pair_finder_.setTransformation(dim, shift);
					}
				}
				else if ( algorithm == "simple" )
				{
					V_computeMatching_("DGeomHashPairwiseMapMatcher: algorithm \"simple\", skip superimposer");
				}
				else
				{
					V_computeMatching_("DGeomHashPairwiseMapMatcher: algorithm not recognized");
				}

				// all_feature_pairs_.clear();
        pair_finder_.setFeaturePairs(all_feature_pairs_);
        pair_finder_.setFeatureMap(1, scene_grid_maps[i]);

				V_computeMatching_("DGeomHashPairwiseMapMatcher: start pairfinder");
        pair_finder_.run();
				
#if 0 // debug output
        String filename("Pairs_file_map1_gridCell_" + String(i));
        std::cout << "Write " << filename << std::endl;
        pair_finder_.dumpFeaturePairs(filename);
#endif

      }
#undef V_computeMatching_

    }


    /// initialize a peak pointer map for each grid cell of the scene (second) map
    void buildGrid_(const PeakConstReferenceMapType& scene_map, std::vector<PeakConstReferenceMapType>& scene_grid_maps)
    {
#define V_buildGrid_(bla) V_DGeomHashPairwiseMapMatcher(bla)
      V_buildGrid_("DGeomHashPairwiseMapMatcher: buildGrid_(): starting...");

      for (Size i = 0; i < scene_map.size(); ++i)
      {
        CoordinateType x = scene_map[i].getPosition()[0] - bounding_box_scene_map_.min()[0];
        CoordinateType y = scene_map[i].getPosition()[1] - bounding_box_scene_map_.min()[1];
        
        Size grid_index = (int)(x / box_size_[0]) + (int)(y / box_size_[1]) * (int)(number_buckets_[0]);
#if 0 // debug
        std::cout << "Grid index " << grid_index << ' ' <<
					scene_map[i].getPosition()[0] << ' ' << scene_map[i].getPosition()[1] << std::endl;
#endif
        scene_grid_maps[grid_index].push_back(scene_map[i]);
      }

      for (Size i = 0; i < scene_grid_maps.size(); ++i)
      {
        // V_buildGrid_((grid_transformation_.getGridCells())[i]);
#if 1 // debug less verbose
				V_buildGrid_("scene_grid_maps["<<i<<"].size(): "<<scene_grid_maps[i].size()<<'\n');
#else
        for (Size j = 0; j < scene_grid_maps[i].size(); ++j)
        {
          V_buildGrid_(((scene_grid_maps[i])[j]).getPosition());
        }
#endif
      }
#undef V_buildGrid_

    }

    /// Parses the parameters, assigns their values to instance members.
    void parseParam_()
    {
#define V_parseParam_(bla) V_DGeomHashPairwiseMapMatcher(bla)
      V_parseParam_("@@@ parseParam_()");

      {
        /// Check the user defined size of the grid cells
        std::string param_name_prefix = "number_buckets:";
        PositionType number_buckets;
        for ( Size dimension = 0; dimension < 2; ++dimension)
        {
          std::string param_name =
            param_name_prefix + DimensionDescriptionType::dimension_name_short[dimension];
          DataValue data_value = param_.getValue(param_name);
          if ( data_value == DataValue::EMPTY )
          {
            throw Exception::ElementNotFound<std::string>
            (__FILE__,__LINE__,__PRETTY_FUNCTION__,param_name);
          }
          else
          {
            number_buckets_[dimension] = data_value;
            V_parseParam_(param_name<< ": "<<number_buckets_[dimension]);
          }
        }
      }
#undef V_parseParam_

    } // parseParam_

    virtual void initGridTransformation_()
    {
#define V_initGridTransformation_(bla) V_DGeomHashPairwiseMapMatcher(bla)
      V_initGridTransformation_("@@@ initGridTransformation_()");

      // compute the minimal and maximal positions of the second map (the map, which should be transformed)
      for ( typename FeatureMapType::ConstIterator fm_iter = feature_map_[1]->begin();
            fm_iter != feature_map_[1]->end();
            ++fm_iter
          )
      {
        bounding_box_scene_map_.enlarge(fm_iter->getPosition());
      }
      V_initGridTransformation_(bounding_box_scene_map_);

      // compute the grid sizes in each dimension
			// ???? I added the "-0.01" adjustment because otherwise this will almost certainly crash when the right margin point comes by!!  Clemens
			// TODO: find a better way that does not use a magic constant
      for (Size i = 0; i < D; ++i)
      {
        box_size_[i] = (bounding_box_scene_map_.max()[i] - bounding_box_scene_map_.min()[i]) / ( number_buckets_[i] - 0.01 );
      }

      // initialize the grid cells of the grid_
      for (Size x_index = 0; x_index < number_buckets_[0]; ++x_index)
      {
        for (Size y_index = 0; y_index < number_buckets_[1]; ++y_index)
        {
          CoordinateType x_min = (bounding_box_scene_map_.min())[0] + box_size_[0]*x_index;
          CoordinateType x_max = (bounding_box_scene_map_.min())[0] + box_size_[0]*(x_index+1);
          CoordinateType y_min = (bounding_box_scene_map_.min())[1] + box_size_[1]*y_index;
          CoordinateType y_max = (bounding_box_scene_map_.min())[1] + box_size_[1]*(y_index+1);

          DGridCell<D,TraitsType> cell(x_min, y_min, x_max, y_max);
          grid_.push_back(cell);
        }
      }

      // V_initGridTransformation_(grid_);

#undef V_initGridTransformation_

    } // initGridTransformation_

  }; // DGeomHashPairwiseMapMatcher


} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DBASEFEATUREMATCHER_H
