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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRPFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRPFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_set_2.h>

#define V_DelaunayPairFinder(bla) // std::cout << bla << std::endl;
#define V_DelaunayConsenus(bla) // std::cout << bla << std::endl;

namespace OpenMS
{

  /**
     @brief This class implements a feature pair finding algorithm.

     This class implements a point pair finding algorithm.
     It works on two point maps, a vector of feature pairs, 
     and a transformation defined for the second feature map (if no
     transformation is given, the pairs are found in the two origianl maps).
     To speed up the search for feature pairs it uses the delaunay triangulation
     of CGAL.

     Policy for copy constructor and assignment: grid_, feature_map_, and
     feature_pairs_ are maintained as pointers and taken shallow copies.  But
     param_ is deep.

  **/
  template < typename MapT = DFeatureMap< 2, DFeature< 2, KernelTraits > >, typename ConsensusFeatureMapT = DFeatureMap< 2, DFeature< 2, KernelTraits > > >
  class DelaunayPairFinder : public BasePairFinder<MapT>
  {
    public:
      typedef DimensionDescription<DimensionDescriptionTagLCMS> DimensionDescriptionType;
      enum DimensionId
      {
        RT = DimensionDescriptionType::RT,
        MZ = DimensionDescriptionType::MZ
    };
      typedef BasePairFinder< MapT > Base;

      // The base knows it all...
      typedef typename Base::TraitsType             TraitsType;
      typedef typename Base::QualityType            QualityType;
      typedef typename Base::PositionType           PositionType;
      typedef typename Base::IntensityType          IntensityType;
      typedef typename Base::PointType              PointType;
      typedef typename Base::PointMapType           PointMapType;
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
      DelaunayPairFinder()
          : Base()
      {
        max_pair_distance_[0] = 0;
        max_pair_distance_[1] = 0;
      }

      /// Copy constructor
      DelaunayPairFinder(const DelaunayPairFinder& source)
          : Base(source)
      {
        max_pair_distance_[0] = source.max_pair_distance_[0];
        max_pair_distance_[1] = source.max_pair_distance_[1];
      }

      ///  Assignment operator
      virtual DelaunayPairFinder& operator = (DelaunayPairFinder source)
      {
        if (&source==this)
          return *this;

        //Base::operator = (source);  FactoryProduct::operator = (source);
        feature_map_[0] = source.feature_map_[0];
        feature_map_[1] = source.feature_map_[1];
        transformation_[0] = source.transformation_[0];
        transformation_[1] = source.transformation_[1];
        max_pair_distance_[0] = source.max_pair_distance_[0];
        max_pair_distance_[1] = source.max_pair_distance_[1];
        return *this;
      }

      /// Destructor
      virtual ~DelaunayPairFinder()
    {}
      //@}

      /// returns an instance of this class
      static BasePairFinder<PointMapType>* create()
      {
        return new DelaunayPairFinder();
      }

      /// returns the name of this module
      static const String getName()
      {
        return "delaunay";
      }

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

          inline Point(Base::RT hx, Base::RT hy, const PointType& f, UnsignedInt k=0)
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

          const PointType* feature;
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
        V_DelaunayPairFinder("DelaunayPairFinder::run(): parse parameters");

        parseParam_();

        V_DelaunayPairFinder("DelaunayPairFinder::run(): find feature pairs");

        findFeaturePairs();
      };

      /// The actual algorithm for finding feature pairs.
      void findFeaturePairs()
      {
        parseParam_();

        const PointMapType& reference_map = *(feature_map_[0]);
        const PointMapType& transformed_map = *(feature_map_[1]);

#define V_findFeaturePairs_(bla) V_DelaunayPairFinder(bla)

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
        std::vector< std::pair< const PointType*,const PointType*> > all_feature_pairs;

        UnsignedInt index_act_feature_pair = 0;
        // take each point in the first data map and search for its neighbours in the second feature map (within a given (transformed) range)
        for ( Size fi1 = 0; fi1 < transformed_map.size(); ++fi1 )
        {
          // compute the transformed iso-rectangle (upper_left,bottom_left,bottom_right,upper_right) for the range query
          double rt_pos = transformed_map[fi1].getPosition()[RT];
          double mz_pos = transformed_map[fi1].getPosition()[MZ] / (diff_intercept_[MZ] / diff_intercept_[RT]);

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

          // if the first neighbour is close enough to act_pos and the second_nearest neighbour lies far enough from the nearest neighbour
          Point nearest = resulting_range[0]->point();
          Point second_nearest = resulting_range[1]->point();

          if (((fabs(transformed_pos[RT] - nearest.hx())  < max_pair_distance_[RT])
               &&  (fabs(transformed_pos[MZ] - nearest.hy())  < max_pair_distance_[MZ]))
              && ((fabs(second_nearest.hx() - nearest.hx())  > max_pair_distance_[RT])
                  || (fabs(second_nearest.hy() - nearest.hy())  > max_pair_distance_[MZ])))
          {
            all_feature_pairs.push_back(std::pair<const PointType*,const PointType*>(nearest.feature,&transformed_map[fi1]));

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

        for (Size i = 0; i < n; ++i)
        {
          SignedInt pair_key = lookup_table[i];
          if ( pair_key > -1 )
          {
            feature_pairs_->push_back(FeaturePairType(*(all_feature_pairs[pair_key].second),*(all_feature_pairs[pair_key].first)));
          }
        }
#undef V_findFeaturePairs_

      } // findFeaturePairs_


      /// The actual algorithm for consensus map computing.
      template < typename ResultMapType >
      void computeConsensusMap(const PointMapType& first_map, ResultMapType& second_map)
      {

#define V_computeConsenusMap(bla) V_DelaunayConsenus(bla)

        V_computeConsenusMap("@@@ computeConsensusMap()");

        parseParam_();

        // Vector to fill the point set for triangulation
        std::vector< Point > positions_reference_map;
        Size n = first_map.size();
        for (Size i = 0; i < n; ++i)
        {
          positions_reference_map.push_back(Point((double)(first_map[i].getPosition()[RT]),
                                                  (double)(first_map[i].getPosition()[MZ] / (diff_intercept_[MZ] / diff_intercept_[RT])),first_map[i],i));
        }

        // compute the delaunay triangulation
        Point_set_2 p_set(positions_reference_map.begin(),positions_reference_map.end());

        // Initialize a hash map for the features of reference_map to avoid that features of the reference map occur in several feature pairs
        std::vector< SignedInt > lookup_table(n,-1);
        std::vector< std::pair< const PointType*, PointType*> > all_feature_pairs;

        UnsignedInt trans_single = 0;
        UnsignedInt ref_single = 0;
        UnsignedInt pairs = 0;
        UnsignedInt index_act_feature_pair = 0;
        // take each point in the first data map and search for its neighbours in the second feature map (within a given (transformed) range)
        for ( Size fi1 = 0; fi1 < second_map.size(); ++fi1 )
        {
          // compute the transformed iso-rectangle (upper_left,bottom_left,bottom_right,upper_right) for the range query
          double rt_pos = (double)(second_map[fi1].getPosition()[RT]);
          double mz_pos = (double)(second_map[fi1].getPosition()[MZ] / (diff_intercept_[MZ]/diff_intercept_[RT]));

          V_computeConsenusMap("Search for two nearest neighbours of " << rt_pos << ' ' << second_map[fi1].getPosition()[MZ] );
          Point transformed_pos(rt_pos,mz_pos,second_map[fi1]);

          V_computeConsenusMap("Transformed Position is : " << transformed_pos );
          std::vector< Vertex_handle > resulting_range;
          p_set.nearest_neighbors(transformed_pos,2,std::back_inserter(resulting_range));

          V_computeConsenusMap("Neighbouring points : ");
          for (typename std::vector< Vertex_handle >::const_iterator it = resulting_range.begin(); it != resulting_range.end(); it++)
          {
            V_computeConsenusMap((*it)->point());
            V_computeConsenusMap(*((*it)->point().feature));
          }

          Point nearest = resulting_range[0]->point();
          Point second_nearest = resulting_range[1]->point();
          if (((fabs(transformed_pos[RT] - nearest.hx())  < max_pair_distance_[RT])
               &&  (fabs(transformed_pos[MZ] - nearest.hy())  < max_pair_distance_[MZ]))
              && ((fabs(second_nearest.hx() - nearest.hx())  > max_pair_distance_[RT])
                  || (fabs(second_nearest.hy() - nearest.hy())  > max_pair_distance_[MZ])))
          {
            all_feature_pairs.push_back(std::pair<const PointType*,PointType*>(nearest.feature,&(second_map[fi1])));

            SignedInt feature_key = resulting_range[0]->point().key;
            // if the feature a is already part of a FeaturePair (a,b) do:
            //    if (the feature c closer to a than b to a) and (the distance between c and b is > a given threshold) do:
            //    --> push (a,c)
            //    else 
            //    --> the value in the lookup_table becomes -2 because the mapping is not unique
            if ( lookup_table[feature_key] > -1)
            {
              SignedInt pair_key = lookup_table[feature_key];
              const PointType& first_map_a = *(all_feature_pairs[pair_key].first);
              PointType& second_map_b = *(all_feature_pairs[pair_key].second);
              PointType& second_map_c = second_map[fi1];

              V_computeConsenusMap("The element " << first_map_a.getPosition() << " has two element partners \n");
              V_computeConsenusMap(second_map_b.getPosition() << "  and  " << second_map_c.getPosition());

              V_computeConsenusMap("Range " << second_map_b.getPositionRange() << "  and  " << second_map_c.getPositionRange());

              if (second_map_c.getPositionRange().encloses(first_map_a.getPosition())
                  && !second_map_b.getPositionRange().encloses(first_map_a.getPosition()))
              {
                lookup_table[feature_key] = index_act_feature_pair;
              }
              else
              {
                // if second_map_b and second_map_c do not enclose first_map_a
                if (!(second_map_b.getPositionRange().encloses(first_map_a.getPosition())
                      && !second_map_c.getPositionRange().encloses(first_map_a.getPosition())))
                {
                  // check the distance between second_map_b and second_map_c
                  if (fabs(second_map_b.getPosition()[MZ] / (diff_intercept_[MZ]/diff_intercept_[RT])
                           - second_map_c.getPosition()[MZ] / (diff_intercept_[MZ]/diff_intercept_[RT]))
                      > max_pair_distance_[MZ])
                  {
                    // and check which one of the elements lies closer to first_map_a
                    if( sqrt(pow((first_map_a.getPosition()[RT] - second_map_b.getPosition()[RT]), 2)
                             + pow((first_map_a.getPosition()[MZ] - second_map_b.getPosition()[RT]), 2))
                        < sqrt(pow((first_map_a.getPosition()[RT] - second_map_c.getPosition()[RT]), 2)
                               + pow((first_map_a.getPosition()[MZ] - second_map_c.getPosition()[RT]), 2)))
                    {
                      lookup_table[feature_key] = index_act_feature_pair;
                    }
                  }
                  else
                  {
                    lookup_table[feature_key] = -2;
                    ++trans_single;
                  }
                }
              }
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
          // no corresponding feature in reference map
          // add a singleton consensus feature
          else
          {
            ++trans_single;
          }
        }

        std::vector< const PointType* > single_elements_first_map;
        for (Size i = 0; i < n; ++i)
        {
          SignedInt pair_key = lookup_table[i];
          if ( pair_key > -1 )
          {
            IndexTuple< ConsensusFeatureMapT > index_tuple(*((all_feature_pairs[pair_key].first)->begin()));
            (all_feature_pairs[pair_key].second)->insert(index_tuple);
            ++pairs;
          }
          // add a singleton consensus feature
          else
          {
            single_elements_first_map.push_back(positions_reference_map[i].feature);
            ++ref_single;
          }
        }

        UnsignedInt length = single_elements_first_map.size();
        for (UnsignedInt i = 0; i < length; ++i)
        {
          second_map.push_back(*(single_elements_first_map[i]));
        }

        V_computeConsenusMap("SINGLE TRANS: " << trans_single);
        V_computeConsenusMap("SINGLE REF: " << ref_single);
        V_computeConsenusMap("PAIRS: " << pairs);

#undef V_computeConsensusMap

      } // computeConsensusMap


    protected:

      /** @name Data members
       */
      //@{

      /// A parameter for similarity_().
      double diff_intercept_[2];
      // allowed distance from a feature neighbour
      float max_pair_distance_[2];
      //@}

      /// Parses the parameters, assigns their values to instance members.
      void parseParam_()
      {
#define V_parseParam_(bla) V_DelaunayPairFinder(bla)
        V_parseParam_("@@@ parseParam_()");

        String param_name_prefix = "similarity:max_pair_distance:";
        for ( Size dimension = 0; dimension < 2; ++dimension)
        {
          std::string param_name = param_name_prefix + DimensionDescriptionType::dimension_name_short[dimension];
          DataValue data_value = param_.getValue(param_name);
          if ( data_value == DataValue::EMPTY )
          {
            throw Exception::ElementNotFound<std::string>
            (__FILE__,__LINE__,__PRETTY_FUNCTION__,param_name);
          }
          else
          {
            max_pair_distance_[dimension] = (float)data_value;
            V_parseParam_(param_name<< ": "<< diff_intercept_[dimension]);
          }
        }

        param_name_prefix = "similarity:diff_intercept:";
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
            diff_intercept_[dimension] = (double)data_value;
            V_parseParam_(param_name<< ": "<< diff_intercept_[dimension]);
          }
        }
#undef V_parseParam_
      } // parseParam_
  }
  ; // DelaunayPairFinder
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DelaunayPairFinder_H
