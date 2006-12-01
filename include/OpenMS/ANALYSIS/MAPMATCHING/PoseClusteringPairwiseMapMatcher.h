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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGPAIRWISEMAPMATCHER_H
#define OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGPAIRWISEMAPMATCHER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairwiseMapMatcher.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer_registerChildren.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>
#ifdef CGAL_DEF
#  include <OpenMS/ANALYSIS/MAPMATCHING/DelaunayPairFinder.h>
#endif
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder_registerChildren.h>
#include <OpenMS/KERNEL/DPeakConstReferenceArray.h>
#include <OpenMS/CONCEPT/Factory.h>


#if defined OPENMS_DEBUG && ! defined PoseClusteringPairwiseMapMatcher
#define V_PoseClusteringPairwiseMapMatcher(bla)  std::cout << bla << std::endl;
#else
#define V_PoseClusteringPairwiseMapMatcher(bla)
#endif


namespace OpenMS
{
  /**
     @brief This class represents a feature matching algorithm.

     It works on two feature maps and computes a vector of corresponding features
     in both maps (given by a feature pairs vector). 
     Therefore it first estimates the shift, which maps one feature map
     onto the other. This is done by the PoseClusteringShiftSuperimposer or the PoseClusteringAffineSuperimposer class. 
     Afterwards using the transformation the corresponding features in both maps 
     can be determined using the DSimplePairFinder class.
     
     If a piecewise shift is assumed, the user can define a grid by setting 
     the number of buckets in the RT as well as the MZ dimension.     
     
  **/
  template < typename MapT = DFeatureMap< 2, DFeature< 2, KernelTraits > > >
  class PoseClusteringPairwiseMapMatcher : public BasePairwiseMapMatcher<MapT>
  {
    public:
      typedef DimensionDescription<DimensionDescriptionTagLCMS> DimensionDescriptionType;
      /// Defines the coordinates of peaks / features.
      enum DimensionId
      {
        RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
        MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

      /** Symbolic names for indices of feature maps etc.
            This should make things more understandable and maintainable.
             */
      enum Maps
      {
        MODEL = 0,
        SCENE = 1
    };


      /** @name Type definitions
       */
      //@{
      /// Base
      typedef BasePairwiseMapMatcher<MapT> Base;

      typedef typename Base::PointMapType PointMapType;
      typedef typename Base::FeatureType FeatureType;
      typedef typename Base::TraitsType TraitsType;
      typedef typename Base::PositionType PositionType;
      typedef typename Base::CoordinateType CoordinateType;

      typedef DLinearMapping<1, TraitsType> TransformationType;

      typedef DPeakConstReferenceArray< PointMapType > PeakConstReferenceMapType;

      using Base::param_;
      using Base::feature_map_;
      using Base::grid_;
      using Base::all_feature_pairs_;
      using Base::bounding_box_scene_map_;
      using Base::box_size_;
      using Base::number_buckets_;
      //@}


      ///@name Constructors, destructor and assignment
      //@{
      /// Constructor
      PoseClusteringPairwiseMapMatcher()
          : Base()
      {
        pair_finder_ = 0;
        superimposer_ = 0;
      }

      /// Copy constructor
      PoseClusteringPairwiseMapMatcher(const PoseClusteringPairwiseMapMatcher& source)
          : Base(source)
      {
        pair_finder_ = Factory<BasePairFinder<PeakConstReferenceMapType> >::create((source.pair_finder_)->getName());
        superimposer_ = Factory<BaseSuperimposer<PeakConstReferenceMapType> >::create((source.superimposer_)->getName());
      }

      ///  Assignment operator
      PoseClusteringPairwiseMapMatcher& operator= (const PoseClusteringPairwiseMapMatcher& source)
      {
        if (&source==this)
          return *this;

        Base::operator = (source);
        pair_finder_ = Factory< BasePairFinder<PeakConstReferenceMapType> >::create((source.pair_finder_)->getName());
        superimposer_ = Factory<BaseSuperimposer<PeakConstReferenceMapType> >::create((source.superimposer_)->getName());

        return *this;
      }

      /// Destructor
      virtual ~PoseClusteringPairwiseMapMatcher()
    {}
      //@}

      /// returns an instance of this class
      static BasePairwiseMapMatcher<MapT>* create()
      {
        return new PoseClusteringPairwiseMapMatcher();
      }


      /// returns the name of this module
      static const String getName()
      {
        return "poseclustering_pairwise";
      }


      /** @name Accesssor methods
       */
      //@{
      /// TODO get and set feature pairs of a given grid cell?, and all feature pairs
      //@}


      /// Estimates the transformation for each grid cell
      virtual void run()
      {
        DataValue data_value = param_.getValue("pair_finder");
        pair_finder_ = Factory<BasePairFinder<PeakConstReferenceMapType> >::create(data_value);
        pair_finder_->setParam(param_);

        data_value = param_.getValue("superimposer");
        if (data_value != DataValue::EMPTY)
        {
          superimposer_ = Factory<BaseSuperimposer<PeakConstReferenceMapType> >::create(data_value);
          superimposer_->setParam(param_);
        }

        // compute the bounding boxes of the grid cells and initialize the grid transformation
        if (grid_.size() == 0)
        {
          initGridTransformation(*(feature_map_[SCENE]));
        }

        // assign each feature of the scene map to the grid cells and build a pointer map for each grid cell
        PeakConstReferenceMapType scene_pointer_map(feature_map_[SCENE]->begin(), feature_map_[SCENE]->end());
        Size number_grid_cells = grid_.size();
        std::vector<PeakConstReferenceMapType> scene_grid_maps(number_grid_cells);
        buildGrid_(scene_pointer_map,scene_grid_maps);

        // initialize a pointer map with the features of the first (model or reference) map
        PeakConstReferenceMapType model_pointer_map(feature_map_[MODEL]->begin(), feature_map_[MODEL]->end());

        // compute the matching of each scene's grid cell features and all the features of the model map
        computeMatching_(model_pointer_map,scene_grid_maps);

        //         String all_feature_pairs_gnuplot_file =
        //           param_.getValue("debug:all_feature_pairs_gnuplot_file");
        //         if ( !all_feature_pairs_gnuplot_file.empty() )
        //         {
        //           pair_finder_->dumpFeaturePairs(all_feature_pairs_gnuplot_file);
        //         }
      }

    protected:

      /** @name Data members
       */
      //@{
      /// This class computes the shift for the best mapping of one feature map to another
      BaseSuperimposer<PeakConstReferenceMapType>* superimposer_;
      /// Given the shift, the pair_finder_ searches for corresponding features in the two maps
      BasePairFinder< PeakConstReferenceMapType >* pair_finder_;
      //@}

      /// compute the matching of each grid cell
      void computeMatching_(const PeakConstReferenceMapType& model_map, const std::vector<PeakConstReferenceMapType>& scene_grid_maps)
      {
#define V_computeMatching_(bla) V_PoseClusteringPairwiseMapMatcher(bla)

        pair_finder_->setFeatureMap(MODEL, model_map);

        if (superimposer_ !=0 )
        {
          // same model map for all superpositions
          superimposer_->setFeatureMap(MODEL, model_map);
        }

        String shift_buckets_file = param_.getValue("debug:shift_buckets_file");
        String feature_buckets_file = param_.getValue("debug:feature_buckets_file");

        // iterate over all grid cells of the scene map
        for (Size i = 0; i < scene_grid_maps.size(); ++i)
        {
          String algorithm = param_.getValue("superimposer");
          if ( superimposer_ != 0 )
          {
            V_computeMatching_("PoseClusteringPairwiseMapMatcher:  superimposer \"pose_clustering\", start superimposer");

            superimposer_->setFeatureMap(SCENE, scene_grid_maps[i]);

            //             // optional debug output: buckets of the shift map
            //             if ( !shift_buckets_file.empty() )
            //             {
            //               superimposer_->getParam().setValue("debug:shift_buckets_file", shift_buckets_file+String(i).fillLeft('0',3));
            //             }
            //
            //             // optional debug output: buckets of the feature maps
            //             if ( !feature_buckets_file.empty() )
            //             {
            //               superimposer_->getParam().setValue("debug:feature_buckets_file", feature_buckets_file+String(i).fillLeft('0',3));
            //             }

            superimposer_->run();

            // ???? copy array to vector -- but the old Grid class will be replaced anyway
            grid_[i].getMappings().resize(2,0);
            for ( Size dim = 0; dim < 2; ++dim )
            {
              TransformationType const& trafo = superimposer_->getTransformation(dim);
              if ( !grid_[i].getMappings()[dim] )
              {
                grid_[i].getMappings()[dim] = new TransformationType;
              }
              *grid_[i].getMappings()[dim] = trafo;
              pair_finder_->setTransformation(dim, trafo);
            }
          }
          else
          {
            V_computeMatching_("PoseClusteringPairwiseMapMatcher: algorithm \"simple\", skip superimposer");
          }

          pair_finder_->setFeaturePairs(*all_feature_pairs_);
          pair_finder_->setFeatureMap(SCENE, scene_grid_maps[i]);

          V_computeMatching_("PoseClusteringPairwiseMapMatcher: start pairfinder");
          pair_finder_->run();
        }
#undef V_computeMatching_

      }


      /// initialize a peak pointer map for each grid cell of the scene (second) map
      void buildGrid_(const PeakConstReferenceMapType& scene_map, std::vector<PeakConstReferenceMapType>& scene_grid_maps)
      {
#define V_buildGrid_(bla) V_PoseClusteringPairwiseMapMatcher(bla)
        V_buildGrid_("PoseClusteringPairwiseMapMatcher: buildGrid_(): starting...");

        for (Size i = 0; i < scene_map.size(); ++i)
        {
          CoordinateType x = scene_map[i].getPosition()[RT] - bounding_box_scene_map_.min()[RT];
          CoordinateType y = scene_map[i].getPosition()[MZ] - bounding_box_scene_map_.min()[MZ];

          Size grid_index = (int)(x / box_size_[RT]) + (int)(y / box_size_[MZ]) * (int)(number_buckets_[RT]);
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
#define V_parseParam_(bla) V_PoseClusteringPairwiseMapMatcher(bla)
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


  }
  ; // PoseClusteringPairwiseMapMatcher


} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DBASEFEATUREMATCHER_H
