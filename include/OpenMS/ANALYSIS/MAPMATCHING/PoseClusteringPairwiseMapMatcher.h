// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder_registerChildren.h>
#include <OpenMS/KERNEL/DPeakConstReferenceArray.h>
#include <OpenMS/CONCEPT/Factory.h>


#define V_PoseClusteringPairwiseMapMatcher(bla) // std::cout << bla << std::endl;



namespace OpenMS
{
  /**
     @brief This class represents a point matching algorithm.
     
     It works on two point maps and computes a vector of corresponding points
     in both maps (given by a point pairs vector). 
     A point can be a DPeak, a DFeature, a ConsensusPeak or ConsensusFeature 
     (wheras DFeature is the default element type).

     Therefore it first estimates the transformation, which maps one point map
     onto the other. This is done by a superimposer class 
     (e.g. PoseClusteringShiftSuperimposer or the PoseClusteringAffineSuperimposer).
     Afterwards a pairfinder class (e.g. SimplePairFinder or DelaunayPairFinder) 
     dewarps the data and searches for corresponding points in both maps.
          
     @note If a piecewise transformation is assumed, the user can define a grid by setting 
     the number of buckets in the RT as well as the MZ dimension.   
     Call initGridTransformation() before run()!   
     
     @ingroup MatchingAlgorithm
  */
  template < typename MapT = DFeatureMap< 2, DFeature< 2, KernelTraits > > >
  class PoseClusteringPairwiseMapMatcher 
  	: public BasePairwiseMapMatcher<MapT>
  {
  public:
    typedef DimensionDescription<LCMS_Tag> DimensionDescriptionType;
    /// Defines the coordinates of elements
    enum DimensionId
    {
      RT = DimensionDescription < LCMS_Tag >::RT,
      MZ = DimensionDescription < LCMS_Tag >::MZ
    };
    /** Symbolic names for indices of element maps etc.
          This should make things more understandable and maintainable.
           */
    enum Maps
    {
      MODEL = 0,
      SCENE = 1
    };

    /// Base
    typedef BasePairwiseMapMatcher<MapT> Base;

    typedef typename Base::PointMapType PointMapType;
    typedef typename Base::ElementType ElementType;
    typedef typename Base::TraitsType TraitsType;
    typedef typename Base::PositionType PositionType;
    typedef typename Base::CoordinateType CoordinateType;

    typedef DLinearMapping<1, TraitsType> TransformationType;

    typedef DPeakConstReferenceArray< PointMapType > PeakConstReferenceMapType;

    using Base::param_;
    using Base::defaults_;
    using Base::subsections_;
    using Base::element_map_;
    using Base::grid_;
    using Base::all_element_pairs_;
    using Base::bounding_box_scene_map_;
    using Base::box_size_;
    using Base::number_buckets_;


    /// Constructor
    PoseClusteringPairwiseMapMatcher()
        : Base(),
        	superimposer_(0),
        	pair_finder_(0)
    {
    	setName(getProductName());
    	
      defaults_.setValue("pairfinder:type", "simple");
			defaults_.setValue("superimposer:type", "none");
			subsections_.push_back("debug");
			subsections_.push_back("pairfinder");
			subsections_.push_back("superimposer");
			
			Base::defaultsToParam_();
    }

    /// Copy constructor
    PoseClusteringPairwiseMapMatcher(const PoseClusteringPairwiseMapMatcher& source)
        : Base(source),
        	superimposer_(0),
        	pair_finder_(0)
    {
    	updateMembers_();	
    }

    ///  Assignment operator
    PoseClusteringPairwiseMapMatcher& operator= (const PoseClusteringPairwiseMapMatcher& source)
    {
      if (&source==this) return *this;

      Base::operator = (source);
			
			updateMembers_();	
			
      return *this;
    }

    /// Destructor
    virtual ~PoseClusteringPairwiseMapMatcher()
  	{
  	}

    /// Returns an instance of this class
    static BasePairwiseMapMatcher<MapT>* create()
    {
      return new PoseClusteringPairwiseMapMatcher();
    }

    /// Returns the name of this module
    static const String getProductName()
    {
      return "poseclustering_pairwise";
    }

    /// Estimates the transformation for each grid cell and searches for element pairs.
    virtual void run()
    {
      // compute the bounding boxes of the grid cells and initialize the grid transformation
      if (grid_.size() == 0)
      {
        initGridTransformation(*(element_map_[SCENE]));
      }

      // assign each element of the scene map to the grid cells and build a pointer map for each grid cell
      PeakConstReferenceMapType scene_pointer_map(element_map_[SCENE]->begin(), element_map_[SCENE]->end());
      Size number_grid_cells = grid_.size();
      std::vector<PeakConstReferenceMapType> scene_grid_maps(number_grid_cells);
      buildGrid_(scene_pointer_map,scene_grid_maps);

      // initialize a pointer map with the elements of the first (model or reference) map
      PeakConstReferenceMapType model_pointer_map(element_map_[MODEL]->begin(), element_map_[MODEL]->end());

      // compute the matching of each scene's grid cell elements and all the elements of the model map
      computeMatching_(model_pointer_map,scene_grid_maps);


      //         String all_element_pairs_gnuplot_file =
      //           param_.getValue("debug:all_element_pairs_gnuplot_file");
      //         if ( !all_element_pairs_gnuplot_file.empty() )
      //         {
      //           pair_finder_->dumpFeaturePairs(all_element_pairs_gnuplot_file);
      //         }
    }

  protected:
   virtual void updateMembers_()
    {
      Base::updateMembers_();
			
			//create pairfinder if it does not exist or if it changed
			if (pair_finder_==0 || pair_finder_->getName()!=param_.getValue("pairfinder:type"))
			{
				delete pair_finder_;
      	pair_finder_ = Factory<BasePairFinder<PeakConstReferenceMapType> >::create(param_.getValue("pairfinder:type"));
      }
      //update pairfinder parameters if necessary
      Param param_copy = param_.copy("pairfinder:",true);
      param_copy.remove("type");
      if (!(pair_finder_->getParameters()==param_copy))
      {
      	pair_finder_->setParameters(param_copy);
     	}
     	 
      //update superimposer
      String type = param_.getValue("superimposer:type");
      if (superimposer_==0)
      {
      	if (type != "none")
      	{
      		delete superimposer_;
        	superimposer_ = Factory<BaseSuperimposer<PeakConstReferenceMapType> >::create(type);
      	}
      }
      else
      {
      	if (type == "none")
      	{
      		delete superimposer_;
      	}
      	else if (superimposer_->getName()!=type)
      	{
      		delete superimposer_;
        	superimposer_ = Factory<BaseSuperimposer<PeakConstReferenceMapType> >::create(type);
      	}
      }
      //update superimposer parameters if necessary
      if (superimposer_!=0)
      {
      	Param param_copy = param_.copy("superimposer:",true);	
      	param_copy.remove("type");
	      if (!(superimposer_->getParameters()==param_copy))
	      {
	      	superimposer_->setParameters(param_copy);
	     	}
	  	}
    }
    
    /// This class computes the shift for the best mapping of one element map to another
    BaseSuperimposer<PeakConstReferenceMapType>* superimposer_;
    /// Given the shift, the pair_finder_ searches for corresponding elements in the two maps
    BasePairFinder< PeakConstReferenceMapType >* pair_finder_;

    /// Computes the matching between each grid cell in the scene map and the model map
    void computeMatching_(const PeakConstReferenceMapType& model_map, const std::vector<PeakConstReferenceMapType>& scene_grid_maps)
    {
#define V_computeMatching_(bla) V_PoseClusteringPairwiseMapMatcher(bla)

      pair_finder_->setElementMap(MODEL, model_map);

      if (superimposer_ !=0 )
      {
        // same model map for all superpositions
        superimposer_->setElementMap(MODEL, model_map);
      }

      String shift_buckets_file = param_.getValue("debug:shift_buckets_file");
      String element_buckets_file = param_.getValue("debug:feature_buckets_file");

      // iterate over all grid cells of the scene map
      for (Size i = 0; i < scene_grid_maps.size(); ++i)
      {
        if (scene_grid_maps[i].size() > 0)
        {
          if ( superimposer_ != 0 )
          {
            V_computeMatching_("PoseClusteringPairwiseMapMatcher:  superimposer \"pose_clustering\", start superimposer");

            superimposer_->setElementMap(SCENE, scene_grid_maps[i]);
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


          pair_finder_->setElementPairs(all_element_pairs_);
          pair_finder_->setElementMap(SCENE, scene_grid_maps[i]);

          V_computeMatching_("PoseClusteringPairwiseMapMatcher: start pairfinder");
          pair_finder_->findElementPairs();
        }
      }
#undef V_computeMatching_

    }


    /// Initializes a peak pointer map for each grid cell of the scene (second) map
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

      //       for (Size i = 0; i < scene_grid_maps.size(); ++i)
      //       {
      //         V_buildGrid_("scene_grid_maps["<<i<<"].size(): "<<scene_grid_maps[i].size()<<'\n');
      //
      //         for (Size j = 0; j < scene_grid_maps[i].size(); ++j)
      //         {
      //           V_buildGrid_(((scene_grid_maps[i])[j]).getPosition());
      //         }
      //       }
#undef V_buildGrid_

    }
  }
  ; // PoseClusteringPairwiseMapMatcher
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGPAIRWISEMAPMATCHER_H
