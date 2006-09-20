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
// $Id: DDelaunayPairFinder.h,v 1.6 2006/05/02 09:39:07 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DDELAUNAYPAIRPFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DDELAUNAYPAIRPFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DBasePairFinder.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_set_2.h>

#if defined OPENMS_DEBUG && ! defined V_DDelaunayPairFinder
#define V_DDelaunayPairFinder(bla) // std::cout << bla << std::endl;
#else
#define V_DDelaunayPairFinder(bla)
#endif

namespace OpenMS
{

  /**
     @brief This class implements a feature pair finding algorithm.

     This class implements a feature pair finding algorithm.
     It works on two feature maps, a vector of feature pairs, 
     and a transformation defined for the second feature map (if no
     transformation is given, the pairs are found in the two origianl maps).
     To speed up the search for feature pairs it uses the delaunay triangulation
     of CGAL.

     Policy for copy constructor and assignment: grid_, feature_map_, and
     feature_pairs_ are maintained as pointers and taken shallow copies.  But
     param_ is deep.

  **/
  template <Size D, typename Traits = KernelTraits, typename MapT = DFeatureMap< D, DFeature< D, Traits > > >
  class DDelaunayPairFinder : public DBasePairFinder<D,Traits,MapT>
  {
    public:
      typedef DimensionDescription<DimensionDescriptionTagLCMS> DimensionDescriptionType;
      enum DimensionId
      {
        RT = DimensionDescriptionType::RT,
        MZ = DimensionDescriptionType::MZ
    };
      enum { DIMENSION = D };

      typedef DBasePairFinder< D, Traits, MapT > Base;

      // The base knows it all...
      typedef typename Base::TraitsType             TraitsType; 
		  typedef typename Base::QualityType            QualityType;
      typedef typename Base::PositionType           PositionType;
      typedef typename Base::FeatureType            FeatureType;
      typedef typename Base::FeatureMapType         FeatureMapType;
      typedef typename Base::FeaturePairType        FeaturePairType;
      //@}

      using Base::param_;
      using Base::feature_map_;
      using Base::feature_pairs_;
      using Base::transformation_;
      //@}


      ///@name Constructors, destructor and assignment
      //@{
      /// Constructor
      DDelaunayPairFinder()
          : Base(),
          max_pair_distance_(0)
      {}

      /// Copy constructor
      DDelaunayPairFinder(const DDelaunayPairFinder& source)
          : Base(source),
          max_pair_distance_(source.max_pair_distance_)
      {}

      ///  Assignment operator
      virtual DDelaunayPairFinder& operator = (DDelaunayPairFinder source)
      {
        if (&source==this)
          return *this;

        Base::operator = (source);
        max_pair_distance_ = source.max_pair_distance_;
        return *this;
      }

      /// Destructor
      virtual ~DDelaunayPairFinder()
    {}
      //@}

      // nested class, which inherits from the cgal Point_2 class and additionally contains the a reference to
      // the corresponding feature and a unique key
    class Point : public CGAL::Point_2< CGAL::Cartesian<double> >
      {
        public:

          typedef CGAL::Point_2< CGAL::Cartesian<double> > Base;

          inline Point() : Base()
          {
            feature = 0;
            key = 0;
          }

          inline Point(const Base& cgal_point) : Base(cgal_point)
          {
            feature = 0;
            key = 0;
          }

          inline Point(Base::RT hx, Base::RT hy, const FeatureType& f, UnsignedInt k)
              : Base(hx,hy)
          {
            feature = &f;
            key = k;
          }

          inline Point(Base::RT hx, Base::RT hy) : Base(hx,hy)
          {
            feature = 0;
            key = 0;
          }

          ~Point()
          {}

          /// Copy constructor
          Point(const Point& source)
              : Point_2(source),
              key(source.key)
          {
            feature = source.feature;
          }

          ///  Assignment operator
          Point& operator = (const Point& source)
          {
            if (this==&source)
              return *this;

            feature = source.feature;
            key = source.key;
            Base::operator=(source);
            return *this;
          }

          const FeatureType* feature;
          UnsignedInt key;
      };

      // to construct a delaunay triangulation with our Point class we have to write an own
      // geometric traits class and the operator() to generate a Point given a CGAL circle
    class  GeometricTraits : public CGAL::Cartesian<double>
      {
        public:
          typedef Point Point_2;

          class Construct_center_2
          {
              typedef Point   Point_2;
              typedef CGAL::Cartesian<double>::Circle_2  Circle_2;
            public:
              typedef Point_2          result_type;
              typedef CGAL::Arity_tag< 1 >   Arity;

              Point_2
              operator()(const Circle_2& c) const
              {
                return c.center();
              }
          };
      };

      typedef CGAL::Point_set_2< GeometricTraits, CGAL::Triangulation_data_structure_2< CGAL::Triangulation_vertex_base_2< GeometricTraits > > > Point_set_2;
      typedef typename Point_set_2::Vertex_handle Vertex_handle;

      /// Estimates the transformation for each grid cell
      virtual void run()
      {

        V_DDelaunayPairFinder("DDelaunayPairFinder::run(): parse parameters");

        parseParam_();

        V_DDelaunayPairFinder("DDelaunayPairFinder::run(): find feature pairs");

        findFeaturePairs_();
      };

    protected:

      /** @name Data members
       */
      //@{ 
		  /// A parameter for similarity_().
      QualityType diff_intercept_[2];
      // allowed distance from a feature neighbour
      PositionType max_pair_distance_;
      //@}

      /// Parses the parameters, assigns their values to instance members.
      void parseParam_()
      {
#define V_parseParam_(bla) V_DDelaunayPairFinder(bla)
        V_parseParam_("@@@ parseParam_()");
        std::string param_name_prefix = "max_pair_distance:";
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
            max_pair_distance_[dimension] = data_value;
            V_parseParam_(param_name<< ": "<< max_pair_distance_[dimension]);
          }
        }
#undef V_parseParam_

      } // parseParam_


      void findFeaturePairs_()
      {
        parseParam_();

        const FeatureMapType& reference_map = *(feature_map_[0]);
        const FeatureMapType& transformed_map = *(feature_map_[1]);

#define V_findFeaturePairs_(bla) V_DDelaunayPairFinder(bla)

        V_findFeaturePairs_("@@@ findFeaturePairs_()");

        Size n = reference_map.size();

        // Vector to fill the point set for triangulation
        // Penalize a deviation in mz more than in rt: deviation(diff_intercept_[RT]) ~ deviation(diff_intercept_[MZ])
        std::vector< Point > positions_reference_map;
        for (Size i = 0; i < n; ++i)
        {
          positions_reference_map.push_back(Point(reference_map[i].getPosition()[RT],
                                                  reference_map[i].getPosition()[MZ] / (diff_intercept_[MZ] / diff_intercept_[RT]),reference_map[i],i));
        }

        // compute the delaunay triangulation
        Point_set_2 p_set(positions_reference_map.begin(),positions_reference_map.end());

        V_findFeaturePairs_("Translation rt " << transformation_[RT]->getParam());
        V_findFeaturePairs_("Translation mz " << transformation_[MZ]->getParam());

        // Initialize a hash map for the features of reference_map to avoid that features of the reference map occur in several feature pairs
        std::vector< SignedInt > lookup_table(n,-1);
        std::vector< FeaturePairType > all_feature_pairs;

        UnsignedInt index_act_feature_pair = 0;
        // take each point in the first data map and search for its neighbours in the second feature map (within a given (transformed) range)
        for ( Size fi1 = 0; fi1 < transformed_map.size(); ++fi1 )
        {
          // compute the transformed iso-rectangle (upper_left,bottom_left,bottom_right,upper_right) for the range query
          typename TraitsType::RealType rt_pos = transformed_map[fi1].getPosition()[RT];
          typename TraitsType::RealType mz_pos = transformed_map[fi1].getPosition()[MZ] / (diff_intercept_[MZ] / diff_intercept_[RT]);

          V_findFeaturePairs_("Search for two nearest neighbours of " << rt_pos << ' ' << transformed_map[fi1].getPosition()[MZ] );
          transformation_[RT]->apply(rt_pos);
          transformation_[MZ]->apply(mz_pos);
          Point transformed_pos(rt_pos,mz_pos);

          V_findFeaturePairs_("Transformed Position is : " << transformed_pos );

          std::vector< Vertex_handle > resulting_range;
          p_set.nearest_neighbors(transformed_pos,2,std::back_inserter(resulting_range));

          V_findFeaturePairs_("Neighbouring points : ");
          for (typename std::vector< Vertex_handle >::const_iterator it = resulting_range.begin(); it != resulting_range.end(); it++)
          {
            V_findFeaturePairs_((*it)->point());
            V_findFeaturePairs_(*((*it)->point().feature));
          }

          UnsignedInt number_of_neighbours = resulting_range.size();
          if (number_of_neighbours > 0)
          {
            // if the first neighbour is close enough to act_pos and the second_nearest neighbour lies far enough from the nearest neighbour
            Point nearest = resulting_range[0]->point();
            Point second_nearest = resulting_range[1]->point();
            // if the first neighbour is close enough to act_pos and the second_nearest neighbour lies far enough from the nearest neighbour
            if (((fabs(transformed_map[fi1].getPosition()[RT] - nearest.feature->getPosition()[RT])  < max_pair_distance_[RT])
                 &&  (fabs(transformed_map[fi1].getPosition()[MZ] - nearest.feature->getPosition()[MZ])  < max_pair_distance_[MZ]))
                && ((fabs(second_nearest.feature->getPosition()[RT] - nearest.feature->getPosition()[RT])  < max_pair_distance_[RT])
                    ||  (fabs(second_nearest.feature->getPosition()[MZ] - nearest.feature->getPosition()[MZ])  < max_pair_distance_[MZ])))
            {
              all_feature_pairs.push_back(FeaturePairType(*(nearest.feature),transformed_map[fi1]));

              SignedInt feature_key = resulting_range[0]->point().key;
              // if the feature already part of a FeaturePair the value in the lookup_table becomes -2
              if ( lookup_table[feature_key] > -1)
              {
                lookup_table[feature_key] = -2;
              }
              // otherwise if the feature is until now no part of a feature pair,
              // set the value in the lookup_table to the index of the pair in the all_feature_pairs vector
              else
              {
                if ( lookup_table[feature_key] == -1)
                {
                  lookup_table[feature_key] = index_act_feature_pair;
                }
              }
              ++index_act_feature_pair;
            }
          }
        }

        for (Size i = 0; i < n; ++i)
        {
          SignedInt pair_key = lookup_table[i];
          if ( pair_key > -1 )
          {
            feature_pairs_->push_back(all_feature_pairs[pair_key]);
          }
        }

#undef V_findFeaturePairs_

      } // findFeaturePairs_
  }
  ; // DDelaunayPairFinder

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DDelaunayPairFinder_H
