// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                  OpenMS Mass Spectrometry Framework
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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRWISEMAPMATCHER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRWISEMAPMATCHER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>

#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>

#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <utility>
#include <fstream>

namespace OpenMS
{

  /**
     @brief The base class of all pairwise point matching algorithms.

     The base class of all point matching algorithms.
     A point can be a raw data point, a peak or a feature.
     This class defines the basic interface for all point matching
     algorithms.  It works on two point maps and computes a vector of corresponding points
     in both maps (given by a point pairs vector). 
     The point pairs created by the algorithm solve a
     (bipartite) matching problem between two point maps.

     The matching shall minimize the shift in retention time and m/z between
     the two maps after a suitable transformation, yet have large cardinality.

  **/
  template < typename MapT = DFeatureMap< 2, DFeature< 2, KernelTraits > > >
  class BasePairwiseMapMatcher : public FactoryProduct
  {
    public:
      /// Defines the coordinates of peaks / features.
      enum DimensionId
      {
        RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
        MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

      /** @name Type definitions
       */
      //@{
      /// Container for input features
      typedef MapT PointMapType;

      /// Type of features considered here
      typedef typename PointMapType::value_type FeatureType;

      /// Traits type
      typedef typename FeatureType::TraitsType TraitsType;

      /// Type of feature pairs
      typedef DFeaturePair < 2, FeatureType > FeaturePairType;

      /// Container for generated feature pairs
      typedef DFeaturePairVector < 2, FeatureType > FeaturePairVectorType;

      /// Grid
      typedef DGrid<2> GridType;

      /// Position
      typedef DPosition < 2, TraitsType > PositionType;

      ///
      typedef DBoundingBox< 2, TraitsType>  PositionBoundingBoxType;

      /// Coordinate
      typedef typename TraitsType::CoordinateType CoordinateType;
      //@}


      ///@name Constructors, destructor and assignment
      //@{
      /// Constructor
      BasePairwiseMapMatcher()
          : FactoryProduct(),
          param_()
      {
        feature_map_[0] = 0;
        feature_map_[1] = 0;
      }

      /// Copy constructor
      BasePairwiseMapMatcher(const BasePairwiseMapMatcher& source)
          : FactoryProduct(source),
          param_(source.param_),
          bounding_box_scene_map_(source.bounding_box_scene_map_),
          box_size_(source.box_size_),
          number_buckets_(source.number_buckets_)
      {
        feature_map_[0] = source.feature_map_[0];
        feature_map_[1] = source.feature_map_[1];
        grid_ = source.grid_;
      }

      ///  Assignment operator
      BasePairwiseMapMatcher& operator = (const BasePairwiseMapMatcher& source)
      {
        FactoryProduct::operator = (source);

        param_ = source.param_;
        feature_map_[0] = source.feature_map_[0];
        feature_map_[1] = source.feature_map_[1];
        grid_ = source.grid_;
        bounding_box_scene_map_ = source.bounding_box_scene_map_;
        box_size_ = source.box_size_;
        number_buckets_ = source.number_buckets_;
        return *this;
      }

      /// Destructor
      virtual ~BasePairwiseMapMatcher()
      {}
      //@}

      /** @name Accesssor methods
       */
      //@{

      /// Set param
      void setParam(const Param& param)
      {
        param_ = param;
        parseParam_();
      }

      /// Get param (non-mutable)
      const Param& getParam() const
      {
        return param_;
      }


      /// Set feature map
      void setFeatureMap(Size const index, const PointMapType& feature_map)
      {
        feature_map_[index] = &feature_map;
      }

      /// Get feature map (non-mutable)
      const PointMapType& getFeatureMap(Size index) const
      {
        return *feature_map_[index];
      }

      /// Set feature pair list
      void setFeaturePairs(FeaturePairVectorType& feature_pairs)
      {
        all_feature_pairs_ = &feature_pairs;
      }

      /// Get feature pair list (non-mutable)
      const FeaturePairVectorType& getFeaturePairs() const
      {
        return *all_feature_pairs_;
      }

      /// Get grid
      const GridType& getGrid() const
      {
        return grid_;
      }

      /// Get grid
      GridType& getGrid()
      {
        return grid_;
      }

      /// Set grid
      void setGrid(const GridType& grid)
      {
        grid_ = grid;
      }
      //@}

      /// register all derived classes here
      static void registerChildren();

      /// Estimates the transformation for each grid cell
      virtual void run() = 0;

      /// Initializes the grid for the scene map given the number of buckets in rt and mz
      void initGridTransformation(const PointMapType& scene_map)
      {
        // compute the minimal and maximal positions of the second map (the map, which should be transformed)
        for ( typename PointMapType::const_iterator fm_iter = scene_map.begin();
              fm_iter != scene_map.end();
              ++fm_iter
            )
        {
          bounding_box_scene_map_.enlarge(fm_iter->getPosition());
        }

        // compute the grid sizes in each dimension
        // ???? I added the "-0.01" adjustment because otherwise this will almost certainly crash when the right margin point comes by!!  Clemens
        // TODO: find a better way that does not use a magic constant
        for (Size i = 0; i < 2; ++i)
        {
          box_size_[i] = (bounding_box_scene_map_.max()[i] - bounding_box_scene_map_.min()[i]) / ( number_buckets_[i] - 0.01 );
        }

        // initialize the grid cells of the grid_
        for (Size x_index = 0; x_index < number_buckets_[RT]; ++x_index)
        {
          for (Size y_index = 0; y_index < number_buckets_[MZ]; ++y_index)
          {
            CoordinateType x_min = (bounding_box_scene_map_.min())[RT] + box_size_[RT]*x_index;
            CoordinateType x_max = (bounding_box_scene_map_.min())[RT] + box_size_[RT]*(x_index+1);
            CoordinateType y_min = (bounding_box_scene_map_.min())[MZ] + box_size_[MZ]*y_index;
            CoordinateType y_max = (bounding_box_scene_map_.min())[MZ] + box_size_[MZ]*(y_index+1);

            DGridCell<2,TraitsType> cell(x_min, y_min, x_max, y_max);
            grid_.push_back(cell);
          }
        }
      } // initGridTransformation_



      //  int dumpFeaturePairs(const String& filename); // code is below

    protected:

      /** @name Data members
       */
      //@{

      /// Param class containing the parameters for the map matching phase
      Param param_;

      /// Two maps of features to be matched
      PointMapType const * feature_map_[2];

      /// Each element of the vector corresponds to all feature pairs of one gridcell
      FeaturePairVectorType* all_feature_pairs_;

      /// The estimated transformation between the two feature maps
      GridType grid_;

      /// Bounding box of the second map
      PositionBoundingBoxType bounding_box_scene_map_;

      /// Size of the grid cells
      PositionType box_size_;

      /// Number of buckets in each dimension
      PositionType number_buckets_;
      //@}

      /// Parse the parameter
      virtual void parseParam_() = 0;

  }
  ; // BasePairwiseMapMatcher


  //   template < typename MapT >
  //   int BasePairwiseMapMatcher<MapT>::dumpFeaturePairs(const String& filename)
  //   {
  //     // V_dumpFeaturePairs() is used for a few comments about the files being
  //     // written.
  // #define V_dumpFeaturePairs(bla) std::cerr << bla << std::endl;
  //     V_dumpFeaturePairs("### Writing "<<filename);
  //     std::ofstream dump_file(filename.c_str());
  //     dump_file << "# " << filename<< " generated " << Date::now() << std::endl;
  //     dump_file << "# 1:number 2:quality 3:firstRT 4:firstMZ 5:firstIT 6:firstQual 7:secondRT 8:secondMZ 9:secondIT 10:secondQual\n";
  //     for ( Size fp = 0; fp < getFeaturePairs().size(); ++fp )
  //     {
  //       dump_file << fp << ' '
  //       << getFeaturePairs()[fp].getQuality() << ' '
  //       << getFeaturePairs()[fp].getFirst().getPosition()[RT] << ' '
  //       << getFeaturePairs()[fp].getFirst().getPosition()[MZ] << ' '
  //       << getFeaturePairs()[fp].getFirst().getIntensity() << ' '
  //       << getFeaturePairs()[fp].getFirst().getOverallQuality() << ' '
  //       << getFeaturePairs()[fp].getSecond().getPosition()[RT] << ' '
  //       << getFeaturePairs()[fp].getSecond().getPosition()[MZ] << ' '
  //       << getFeaturePairs()[fp].getSecond().getIntensity() << ' '
  //       << getFeaturePairs()[fp].getSecond().getOverallQuality() << ' '
  //       << std::endl;
  //     }
  //     dump_file << "# " << filename << " EOF " << Date::now() << std::endl;
  //     std::string dump_filename_gp = filename + ".gp";
  //     V_dumpFeaturePairs("### Writing "<<dump_filename_gp);
  //     std::ofstream dump_file_gp(dump_filename_gp.c_str());
  //     dump_file_gp << "# " << dump_filename_gp << " generated " << Date::now() << std::endl;
  //     dump_file_gp <<
  //     "# Gnuplot script to view feature pairs\n"
  //     "plot   \"" << filename <<"\" using 3:4 title \"map 1\"\n"
  //     "replot \"" << filename <<"\" using 7:8 title \"map 2\"\n"
  //     "replot \"" << filename <<"\" using 3:4:($7-$3):($8-$4) w vectors nohead title \"pairs\"\n"
  //     ;
  //     dump_file_gp << "# " << dump_filename_gp << " EOF " << Date::now() << std::endl;
  //     V_dumpFeaturePairs("### You can view `"<<filename<<"' using the command line `gnuplot "<<dump_filename_gp<<" -'");
  // #undef V_dumpFeaturePairs
  //
  //     return 0;
  //   }

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DBASEFEATUREMATCHER_H
