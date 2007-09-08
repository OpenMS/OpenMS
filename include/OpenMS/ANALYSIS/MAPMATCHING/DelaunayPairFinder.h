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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/IndexTuple.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_set_2.h>

#define V_DelaunayPairFinder(bla) // std::cout << bla << std::endl;
#define V_DelaunayConsenus(bla) // std::cout << bla << std::endl;

namespace OpenMS
{

  /**
		@brief This class implements an element pair finding algorithm.
		
		This class implements a point pair finding algorithm.
		It offers a method to determine element pairs in two element maps,
		given two point maps and a transformation defined for the second element map (if no
		transformation is given, the pairs are found in the two original maps). 
		The pair finder also offers a method to compute consensus elements given 
		two element maps. This algorithm is similar to the pair finding method as mentioned above,
		but it implies that the scene map is already dewarped.
		     
		To speed up the search for element pairs an consensus elements, the DelaunayPairFinder
		uses the CGAL delaunay triangulation for the nearest neighbour search.
		
		The first template parameter is the type of the consensus map and the second parameter is the type of the element maps.
		
		@note The RT and the MZ dimension are not equivalent, because two elements that differ in RT by 1s (or minute) are 
		more similar than two points that differ in MZ by 1Th. To be able to use the euclidean distance in the nearest neighbour search, 
		we have to transform the elements MZ position m into a new MZ position m'= m / (diff_intercept_RT/diff_intercept_MZ).
		E.g. given diff_intercept_RT=1 and diff_intercept_MZ=0.1 results in 1s difference in RT is similar to 0.1Th difference in MZ.
		 
		@ref DelaunayPairFinder_Parameters are explained on a separate page.  
  */
  template < typename ConsensusMapT = FeatureMap< Feature >, typename ElementMapT = FeatureMap< > >
  class DelaunayPairFinder
        : public BasePairFinder<ConsensusMapT>
  {
  public:

    /** Symbolic names for indices of element maps etc.
        This should make things more understandable and maintainable.
    */
    enum Maps
    {
      MODEL = 0,
      SCENE = 1
    };

    typedef BasePairFinder< ConsensusMapT > Base;

    // The base knows it all...
    typedef typename Base::QualityType            QualityType;
    typedef typename Base::PositionType           PositionType;
    typedef typename Base::IntensityType          IntensityType;
    typedef typename Base::PointType              PointType;
    typedef typename Base::PointMapType           PointMapType;
    typedef typename Base::ElementPairType        ElementPairType;

    using Base::param_;
    using Base::defaults_;
    using Base::element_map_;
    using Base::element_pairs_;
    using Base::transformation_;

    /// Constructor
    DelaunayPairFinder()
        : Base()
    {
    	//set the name for DefaultParamHandler error messages
      Base::setName(getProductName());

      defaults_.setValue("similarity:max_pair_distance:RT",20,"Maximum distance in RT dimension");
      defaults_.setValue("similarity:max_pair_distance:MZ",1,"Maximum distance in m/z dimension");
      defaults_.setValue("similarity:precision:RT",60,"The precision of the element's retention time position.");
      defaults_.setValue("similarity:precision:MZ",0.5,"The precision of the element's m/z position.");
      defaults_.setValue("similarity:diff_intercept:RT",1,"Factor for RT position used to balance the influence of RT and m/z deviations");
      defaults_.setValue("similarity:diff_intercept:MZ",0.01,"Factor for m/z position used to balance the influence of RT and m/z deviations");

      Base::defaultsToParam_();
    }

    /// Destructor
    virtual ~DelaunayPairFinder()
  	{}

    /// Returns an instance of this class
    static BasePairFinder<PointMapType>* create()
    {
      return new DelaunayPairFinder();
    }

    /// Returns the name of this module
    static const String getProductName()
    {
      return "DelaunayPairFinder";
    }

    /// Nested class, which inherits from the cgal Point_2 class and additionally contains the a reference to
    /// the corresponding element and an unique key
  	class Point 
  	: public CGAL::Point_2< CGAL::Cartesian<double> >
    {
    public:

      typedef CGAL::Point_2< CGAL::Cartesian<double> > Base;

      inline Point() : Base()
      {
        element = 0;
        key = 0;
      }

      inline Point(const Base& cgal_point) : Base(cgal_point)
      {
        element = 0;
        key = 0;
      }

      inline Point(Base::RT hx, Base::RT hy, const PointType& f, UInt k=0)
          : Base(hx,hy)
      {
        element = &f;
        key = k;
      }

      inline Point(Base::RT hx, Base::RT hy) : Base(hx,hy)
      {
        element = 0;
        key = 0;
      }

      ~Point()
      {}

      /// Copy constructor
      Point(const Point& source)
          : Point_2(source),
          key(source.key)
      {
        element = source.element;
      }

      ///  Assignment operator
      Point& operator = (const Point& source)
      {
        if (this==&source)
          return *this;

        element = source.element;
        key = source.key;
        Base::operator=(source);
        return *this;
      }

      const PointType* element;
      UInt key;
    };

    /// To construct a delaunay triangulation with our Point class we have to write an own
    /// geometric traits class and the operator() (that generates a Point given a CGAL circle)
  	class  GeometricTraits 
  		: public CGAL::Cartesian<double>
    {
    public:
      typedef Point Point_2;
      ///
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

    /// Get diff intercept
    double getDiffIntercept(UInt dim)
    {
      return diff_intercept_[dim];
    }

    /// Set diff intercept
    void setDiffIntercept(UInt dim, DoubleReal intercept)
    {
         param_.setValue(String("similarity:diff_intercept:") + RawDataPoint2D::shortDimensionName(dim), intercept);
         updateMembers_();   
    }

    /// Get max_pair_distance_
    float getMaxPairDistance(UInt dim)
    {
      if (dim == RawDataPoint2D::RT) return max_pair_distance_[dim];
      else 
        {
          double diff = diff_intercept_[RawDataPoint2D::RT] / diff_intercept_[RawDataPoint2D::MZ];
     
          return max_pair_distance_[dim] / diff;      
        }
    }

    /// Set max_pair_distance_
    void setMaxPairDistance(UInt dim, Real max_pair_distance)
    {
      if (dim == RawDataPoint2D::RT) 
      {
        max_pair_distance_[dim] = max_pair_distance; 
      }
      else 
      {
        double diff = diff_intercept_[RawDataPoint2D::RT] / diff_intercept_[RawDataPoint2D::MZ];
        max_pair_distance_[dim] = max_pair_distance * diff;     
       }
      param_.setValue(String("similarity:max_pair_distance:") + RawDataPoint2D::shortDimensionName(dim), max_pair_distance);
    }

    /// Get precision
    float getPrecision(UInt dim)
    {
      if (dim == RawDataPoint2D::RT) return precision_[dim];
      else 
        {
          double diff = diff_intercept_[RawDataPoint2D::RT] / diff_intercept_[RawDataPoint2D::MZ];
     
          return precision_[dim] / diff;      
        }
      
    }

    /// Set precision
    void setPrecision(UInt dim, Real precision)
    {
      if (dim == RawDataPoint2D::RT) 
      {
        precision_[dim] = precision; 
      }
      else 
      {
        double diff = diff_intercept_[RawDataPoint2D::RT] / diff_intercept_[RawDataPoint2D::MZ];
        precision_[dim] = precision * diff;     
      }
      param_.setValue(String("similarity:precision:") + RawDataPoint2D::shortDimensionName(dim), precision);
    }

    /// The actual algorithm for finding element pairs.
    void findElementPairs()
    {
      const PointMapType& reference_map = *(element_map_[MODEL]);
      const PointMapType& transformed_map = *(element_map_[SCENE]);

#define V_findElementPairs(bla) V_DelaunayPairFinder(bla)

      V_findElementPairs("@@@ findElementPairs_()");

      UInt n = reference_map.size();

      // Vector to fill the point set for triangulation
      // Penalize a deviation in mz more than in rt: deviation(diff_intercept_[RawDataPoint2D::RT]) ~ deviation(diff_intercept_[RawDataPoint2D::MZ])
      std::vector< Point > positions_reference_map;
      for (UInt i = 0; i < n; ++i)
      {
        positions_reference_map.push_back(Point(reference_map[i].getRT(),
                                                reference_map[i].getMZ() / (diff_intercept_[RawDataPoint2D::MZ] / diff_intercept_[RawDataPoint2D::RT]),reference_map[i],i));
      }

      // compute the delaunay triangulation
      Point_set_2 p_set(positions_reference_map.begin(),positions_reference_map.end());

      V_findElementPairs("Translation rt " << transformation_[RawDataPoint2D::RT].getParam());
      V_findElementPairs("Translation mz " << transformation_[RawDataPoint2D::MZ].getParam());

      // Initialize a hash map for the elements of reference_map to avoid that elements of the reference map occur in several element pairs
      std::vector< Int > lookup_table(n,-1);
      std::vector< std::pair< const PointType*,const PointType*> > all_element_pairs;

      UInt index_act_element_pair = 0;
      // take each point in the first data map and search for its neighbours in the second element map (within a given (transformed) range)
      for ( UInt fi1 = 0; fi1 < transformed_map.size(); ++fi1 )
      {
        // compute the transformed iso-rectangle (upper_left,bottom_left,bottom_right,upper_right) for the range query
        double rt_pos = transformed_map[fi1].getRT();
        double mz_pos = transformed_map[fi1].getMZ();

        V_findElementPairs("Search for two nearest neighbours of " << rt_pos << ' ' << transformed_map[fi1].getMZ() );
        transformation_[RawDataPoint2D::RT].apply(rt_pos);
        transformation_[RawDataPoint2D::MZ].apply(mz_pos);

        mz_pos /= (diff_intercept_[RawDataPoint2D::MZ] / diff_intercept_[RawDataPoint2D::RT]);
        Point transformed_pos(rt_pos,mz_pos);

        V_findElementPairs("Transformed Position is : " << transformed_pos );

        std::vector< Vertex_handle > resulting_range;
        p_set.nearest_neighbors(transformed_pos,2,std::back_inserter(resulting_range));

        V_findElementPairs("Neighbouring points : ");
        for (typename std::vector< Vertex_handle >::const_iterator it = resulting_range.begin(); it != resulting_range.end(); it++)
        {
          V_findElementPairs((*it)->point());
          V_findElementPairs(*((*it)->point().element));
        }

        // if the first neighbour is close enough to act_pos and the second_nearest neighbour lies far enough from the nearest neighbour
        Point nearest = resulting_range[0]->point();
        Point second_nearest = resulting_range[1]->point();

        if (((fabs(transformed_pos[RawDataPoint2D::RT] - nearest.hx())  < precision_[RawDataPoint2D::RT])
             &&  (fabs(transformed_pos[RawDataPoint2D::MZ] - nearest.hy())  < precision_[RawDataPoint2D::MZ]))
            && ((fabs(second_nearest.hx() - nearest.hx())  > max_pair_distance_[RawDataPoint2D::RT])
                || (fabs(second_nearest.hy() - nearest.hy())  > max_pair_distance_[RawDataPoint2D::MZ])))
        {
          all_element_pairs.push_back(std::pair<const PointType*,const PointType*>(nearest.element,&transformed_map[fi1]));

          Int element_key = resulting_range[0]->point().key;
          // if the element already part of a ElementPair the value in the lookup_table becomes -2
          if ( lookup_table[element_key] > -1)
          {
            lookup_table[element_key] = -2;
          }
          // otherwise if the element is until now no part of a element pair,
          // set the value in the lookup_table to the index of the pair in the all_element_pairs vector
          else
          {
            if ( lookup_table[element_key] == -1)
            {
              lookup_table[element_key] = index_act_element_pair;
            }
          }
          ++index_act_element_pair;
        }
      }

//       std::ofstream out("pairs.dat",std::ios::out);
      for (UInt i = 0; i < n; ++i)
      {
        Int pair_key = lookup_table[i];
        if ( pair_key > -1 )
        {
          element_pairs_->push_back(ElementPairType(*(all_element_pairs[pair_key].second),*(all_element_pairs[pair_key].first)));
        }
      }
#undef V_findElementPairs

    } // findElementPairs_
    /// The actual algorithm for finding consensus consensus elements.
    /// Elements in the first_map are aligned to elements in the second_map, so the second_map contains the resulting consensus elements.
    template < typename ResultMapType >
    void computeConsensusMap(const PointMapType& first_map, ResultMapType& second_map)
    {
#define V_computeConsensusMap(bla) V_DelaunayConsenus(bla)
      V_computeConsensusMap("@@@ computeConsensusMap()");

      // Vector to fill the point set for triangulation
      std::vector< Point > positions_reference_map;
      UInt n = first_map.size();
      for (UInt i = 0; i < n; ++i)
      {
        positions_reference_map.push_back(Point((double)(first_map[i].getRT()),
                                                (double)(first_map[i].getMZ() / (diff_intercept_[RawDataPoint2D::MZ] / diff_intercept_[RawDataPoint2D::RT])),first_map[i],i));
      }
      StopWatch timer;
      V_computeConsensusMap("Start delaunay triangulation for " << positions_reference_map.size() << " elements");
      // compute the delaunay triangulation
      timer.start();
      Point_set_2 p_set(positions_reference_map.begin(),positions_reference_map.end());
      timer.stop();
      V_computeConsensusMap("End delaunay triangulation after " << timer.getCPUTime() << "s");

      // Initialize a hash map for the elements of reference_map to avoid that elements of the reference map occur in several element pairs
      std::vector< Int > lookup_table(n,-1);
      std::vector< std::pair< const PointType*, PointType*> > all_element_pairs;

      UInt trans_single = 0;
      UInt ref_single = 0;
      UInt pairs = 0;
      UInt index_act_element_pair = 0;
      // take each point in the first data map and search for its neighbours in the second element map (within a given (transformed) range)
      for ( UInt fi1 = 0; fi1 < second_map.size(); ++fi1 )
      {
        // compute the transformed iso-rectangle (upper_left,bottom_left,bottom_right,upper_right) for the range query
        double rt_pos = (double)(second_map[fi1].getRT());
        double mz_pos = (double)(second_map[fi1].getMZ() / (diff_intercept_[RawDataPoint2D::MZ]/diff_intercept_[RawDataPoint2D::RT]));

        V_computeConsensusMap("Search for two nearest neighbours of " << rt_pos << ' ' << second_map[fi1].getMZ() << ' ' << second_map[fi1].getIntensity() );
        Point transformed_pos(rt_pos,mz_pos,second_map[fi1]);

        V_computeConsensusMap("Transformed Position is : " << transformed_pos );
        std::vector< Vertex_handle > resulting_range;
        p_set.nearest_neighbors(transformed_pos,2,std::back_inserter(resulting_range));

        Point nearest;
        Point second_nearest;
        if (resulting_range.size() == 1)
        {
          nearest = resulting_range[0]->point();
          if ((fabs(transformed_pos[RawDataPoint2D::RT] - nearest.hx())  < precision_[RawDataPoint2D::RT])
              &&  (fabs(transformed_pos[RawDataPoint2D::MZ] - nearest.hy())  < precision_[RawDataPoint2D::MZ]))
          {
            all_element_pairs.push_back(std::pair<const PointType*,PointType*>(nearest.element,&(second_map[fi1])));
          }
        }
        else
          if (resulting_range.size() > 1)
          {
            nearest = resulting_range[0]->point();
            second_nearest = resulting_range[1]->point();
            V_computeConsensusMap("Nearest: " << (nearest.element)->getRT() << ' ' << fabs(transformed_pos[RawDataPoint2D::RT] - nearest.hx()) << ' ' 
                << (nearest.element)->getMZ() << ' ' << fabs(transformed_pos[RawDataPoint2D::MZ] - nearest.hy()) << ' '
                                              << (nearest.element)->getIntensity());
                
            V_computeConsensusMap("Second nearest: " << (second_nearest.element)->getRT() << ' ' << second_nearest.hx() << ' ' 
                                                     << (second_nearest.element)->getMZ() << ' ' << second_nearest.hy() << ' '
                                                     << (second_nearest.element)->getIntensity());

            if (((fabs(transformed_pos[RawDataPoint2D::RT] - nearest.hx())  < precision_[RawDataPoint2D::RT])
                 &&  (fabs(transformed_pos[RawDataPoint2D::MZ] - nearest.hy())  < precision_[RawDataPoint2D::MZ]))
                && ((fabs(second_nearest.hx() - nearest.hx())  > max_pair_distance_[RawDataPoint2D::RT])
                    || (fabs(second_nearest.hy() - nearest.hy())  > max_pair_distance_[RawDataPoint2D::MZ])))
            {
              all_element_pairs.push_back(std::pair<const PointType*,PointType*>(nearest.element,&(second_map[fi1])));
              V_computeConsensusMap("Push first: " << *(nearest.element))
              V_computeConsensusMap("Push second: " << second_map[fi1])

              Int element_key = resulting_range[0]->point().key;

              // if the element a is already part of a ElementPair (a,b) do:
              //    if (the element c closer to a than b to a) and (the distance between c and b is > a given threshold) do:
              //    --> push (a,c)
              //    else
              //    --> the value in the lookup_table becomes -2 because the mapping is not unique
              if ( lookup_table[element_key] > -1)
              {
                Int pair_key = lookup_table[element_key];
                const PointType& first_map_a = *(all_element_pairs[pair_key].first);
                PointType& second_map_b = *(all_element_pairs[pair_key].second);
                PointType& second_map_c = second_map[fi1];

                V_computeConsensusMap("The element " << first_map_a.getPosition() << ' ' << first_map_a.getIntensity() << " has two element partners \n");
                V_computeConsensusMap(second_map_b.getPosition() << ' ' << second_map_b.getIntensity() 
                    << "  and  " << second_map_c.getPosition() << ' ' << second_map_c.getIntensity()) ;

                V_computeConsensusMap("Range " << second_map_b.getPositionRange() << "  and  " << second_map_c.getPositionRange());

                if (second_map_c.getPositionRange().encloses(first_map_a.getPosition())
                    && !second_map_b.getPositionRange().encloses(first_map_a.getPosition()))
                {
                  lookup_table[element_key] = index_act_element_pair;
                  V_computeConsensusMap(second_map_c.getPosition() << " and " << first_map_a.getPosition() << " are a pair");
                }
                else
                {
                  // if second_map_b and second_map_c do not enclose first_map_a
                  if (!(second_map_b.getPositionRange().encloses(first_map_a.getPosition())
                        && !second_map_c.getPositionRange().encloses(first_map_a.getPosition())))
                  {
                    V_computeConsensusMap(second_map_b.getPosition() << " and " << first_map_a.getPosition() << " are a pair, but check the distance between c and b");
                    // check the distance between second_map_b and second_map_c
                    if (fabs(second_map_b.getMZ() / (diff_intercept_[RawDataPoint2D::MZ]/diff_intercept_[RawDataPoint2D::RT])
                             - second_map_c.getMZ() / (diff_intercept_[RawDataPoint2D::MZ]/diff_intercept_[RawDataPoint2D::RT]))
                        > max_pair_distance_[RawDataPoint2D::MZ])
                    {
                      V_computeConsensusMap("distance ok");
                      // and check which one of the elements lies closer to first_map_
                      if( sqrt(pow((first_map_a.getRT() - second_map_b.getRT()), 2)
                               + pow((first_map_a.getMZ() - second_map_b.getMZ()), 2))
                          > sqrt(pow((first_map_a.getRT() - second_map_c.getRT()), 2)
                                 + pow((first_map_a.getMZ() - second_map_c.getMZ()), 2)))
                      {
                        lookup_table[element_key] = index_act_element_pair;
                        V_computeConsensusMap("take a and c");
                      }
                    }
                    else
                    {
                      lookup_table[element_key] = -2;
                      ++trans_single;
                    }
                  }
                }
              }
              // otherwise if the element is until now no part of a element pair,
              // set the value in the lookup_table to the index of the pair in the all_element_pairs vector
              else
              {
                if ( lookup_table[element_key] == -1)
                {
                  lookup_table[element_key] = index_act_element_pair;
                }
              }
              ++index_act_element_pair;
            }
            // no corresponding element in reference map
            // add a singleton consensus element
            else
            {
              ++trans_single;
            }
          }
      }
      V_computeConsensusMap("Insert elements ");
      std::vector< const PointType* > single_elements_first_map;
      for (UInt i = 0; i < n; ++i)
      {
        Int pair_key = lookup_table[i];
        if ( pair_key > -1 )
        {
          IndexTuple< ElementMapT > index_tuple(*((all_element_pairs[pair_key].first)->begin()));
          V_computeConsensusMap("First: " << *((all_element_pairs[pair_key].first)))
          V_computeConsensusMap("Second: " << *((all_element_pairs[pair_key].second)))
          (all_element_pairs[pair_key].second)->insert(index_tuple);
          ++pairs;
        }
        // add a singleton consensus element
        else
        {
          single_elements_first_map.push_back(positions_reference_map[i].element);
          ++ref_single;
        }
      }

      UInt length = single_elements_first_map.size();
      for (UInt i = 0; i < length; ++i)
      {
        second_map.push_back(*(single_elements_first_map[i]));
      }

      V_computeConsensusMap("SINGLE TRANS: " << trans_single);
      V_computeConsensusMap("SINGLE REF: " << ref_single);
      V_computeConsensusMap("PAIRS: " << pairs);

#undef V_computeConsensusMap

    } // computeConsensusMap


  protected:
    virtual void updateMembers_()
    {
    	diff_intercept_[RawDataPoint2D::RT] = (double)param_.getValue("similarity:diff_intercept:RT");
      diff_intercept_[RawDataPoint2D::MZ] = (double)param_.getValue("similarity:diff_intercept:MZ");
      
      double diff = diff_intercept_[RawDataPoint2D::RT] / diff_intercept_[RawDataPoint2D::MZ];
      
      max_pair_distance_[RawDataPoint2D::RT] = (float)param_.getValue("similarity:max_pair_distance:RT");
      max_pair_distance_[RawDataPoint2D::MZ] = (float)param_.getValue("similarity:max_pair_distance:MZ") * diff;

      precision_[RawDataPoint2D::RT] = (float)param_.getValue("similarity:precision:RT");
      precision_[RawDataPoint2D::MZ] = (float)param_.getValue("similarity:precision:MZ") * diff;

    }

    /// A parameter for similarity_().
    double diff_intercept_[2];
    /// To uniquely assign an element e1 of the scene map to another element e2 in the model map
    /// all elements in the scene map have to lie at least max_pair_distance_ far from e1 and
    /// all elements in the model map have to lie at least max_pair_distance_ far from e2.
    float max_pair_distance_[2];
    /// Only points that differ not more than precision_ can be assigned as a pair
    float precision_[2];
  }
  ; // DelaunayPairFinder
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H
