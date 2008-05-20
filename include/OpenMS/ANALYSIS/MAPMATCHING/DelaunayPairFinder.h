// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_set_2.h>

#ifdef Debug_DelaunayPairFinder
#define V_(bla) std::cout << __FILE__ ":" << __LINE__ << ": " << bla << std::endl;
#else
#define V_(bla)
#endif

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
			
  The template parameter is the type of the consensus map.
			
  @note The RT and the MZ dimension are not equivalent, because two elements that differ in RT by 1s (or minute) are 
  more similar than two points that differ in MZ by 1Th. To be able to use the euclidean distance in the nearest neighbour search, 
  we have to transform the elements MZ position m into a new MZ position m'= m / (diff_intercept_RT/diff_intercept_MZ).
  E.g. given diff_intercept_RT=1 and diff_intercept_MZ=0.1 results in 1s difference in RT is similar to 0.1Th difference in MZ.
			 
  @ref DelaunayPairFinder_Parameters are explained on a separate page.  

  @todo work out all TODOs in the code

  @ingroup FeatureGrouping
  */
  class DelaunayPairFinder : public BasePairFinder
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

    typedef BasePairFinder Base;

    /// Constructor
    DelaunayPairFinder() : Base()
    {
      //set the name for DefaultParamHandler error messages
      Base::setName(getProductName());

      defaults_.setValue("similarity:max_pair_distance:RT",20.0,"Consider element el_1 of a map map_1 and element el_2  el_2 of a different map map_2. Then el_1 is assigned to el_2 if they are nearest neighbors and if all other elements in map_1 have a distance greater max_pair_distance::RT to el_2 as well as all other elements in map_2 have a distance greater max_pair_distance::RT to el_1.",true);
      defaults_.setValue("similarity:max_pair_distance:MZ",1.0,"Consider element el_1 of a map map_1 and element el_2  el_2 of a different map map_2. Then el_1 is assigned to el_2 if they are nearest neighbors and if all other elements in map_1 have a distance greater max_pair_distance::MZ to el_2 as well as all other elements in map_2 have a distance greater max_pair_distance::MZ to el_1.",true);
      defaults_.setValue("similarity:precision:RT",60.0,"Maximal deviation in RT of two elements in different maps to be considered as a corresponding element.");
      defaults_.setValue("similarity:precision:MZ",0.5,"Maximal deviation in m/z of two elements in different maps to be considered as a corresponding element.");
      defaults_.setValue("similarity:diff_intercept:RT",1.0,"Factor for RT position used to balance the influence of RT and m/z deviations");
      defaults_.setValue("similarity:diff_intercept:MZ",0.01,"Factor for m/z position used to balance the influence of RT and m/z deviations");

      Base::defaultsToParam_();
    }

    /// Destructor
    virtual ~DelaunayPairFinder()
    {}

    /// Returns an instance of this class
    static BasePairFinder* create()
    {
      return new DelaunayPairFinder();
    }

    /// Returns the name of this module
    static const String getProductName()
    {
      return "delaunay";
    }

    /// Nested class, which inherits from the CGAL Point_2 class and
    /// additionally contains a reference to the corresponding element and a
    /// unique key
    struct Point : public CGAL::Point_2< CGAL::Cartesian<double> >
    {
      typedef CGAL::Point_2< CGAL::Cartesian<double> > Base;

      const ConsensusFeature* element;
      UInt key;
			
      /// Default ctor
      inline Point() : Base(), element(0), key(0)
      {
        element = 0;
        key = 0;
      }
			
      /// Ctor from Base class, aka CGAL:Point_2<...>
      inline Point(const Base& cgal_point) : Base(cgal_point), element(0), key(0) {}

      /// Ctor from coordinates, element, and key
      inline Point(double hx, double hy, const ConsensusFeature& element, UInt key) : Base(hx,hy), element(&element), key(key) {}

      /// Ctor from coordinates
      inline Point(double hx, double hy) : Base(hx,hy), element(0), key(0) {}

      /// Dtor
      ~Point() {}

      /// Copy ctor
      Point(const Point& rhs) : Base(rhs), element(rhs.element) , key(rhs.key){}

      ///  Assignment operator
      Point& operator = (const Point& source)
      {
        if (this==&source) return *this;
        Base::operator=(source);
        element = source.element;
        key = source.key;
        return *this;
      }
    };

    /// To construct a Delaunay triangulation with our Point class we have to
    /// write an own geometric traits class and the operator() (that generates
    /// a Point given a CGAL circle)
    class  GeometricTraits : public CGAL::Cartesian<double>
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
    typedef Point_set_2::Vertex_handle Vertex_handle;


    // TODO replace diff_intercept_ by balance_MZ_
    /// Get diff intercept
    double getDiffIntercept(UInt dim)
    {
      return diff_intercept_[dim];
    }

    // TODO replace diff_intercept_ by balance_MZ_
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

    /// documented in base class
    void run(ConsensusMap &result_map)
    {
      V_("@@@ DelaunayPairFinder::run()");

      // Every derived class should set maps_.result_ at the beginning.
      maps_.result_ = &result_map;

			// Check whether map indices are meaningful
			{
				if ( map_index_.model_ < -1 )
				{
					throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
				if ( map_index_.scene_ < -1 )
				{
					throw Exception::OutOfRange(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
			}

      // Check whether input and output are not the same
      {
				if ( maps_.result_ == maps_.model_ )
				{
					// sorry, result must not overwrite model!
					throw Exception::IllegalSelfOperation(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}

				if ( maps_.result_ == maps_.scene_ )
				{
					// sorry, result must not overwrite scene!
					throw Exception::IllegalSelfOperation(__FILE__,__LINE__,__PRETTY_FUNCTION__);
				}
      }

      // Empty output destination
      result_map.clear();
			
      // derived from a literal copy of findElementPairs, TODO: rewrite, improve, clean up
      {

				// TODO Find out whether it is (1) correct and (2) fast if we
				// push_back() the Points into the Delaunay triangulation. Otherwise,
				// use an iterator adapter and construct Point_set_2 p_set from
				// iterator range.
				Point_set_2 p_set;

				for (UInt i = 0; i < getModelMap().size(); ++i)
				{
					// incrementally compute the delaunay triangulation
					p_set.push_back( Point( getModelMap()[i].getRT(), getModelMap()[i].getMZ() * balance_MZ_, getModelMap()[i], i) );
				}

				// We add two outlier points so that we will always have at least
				// three points in the Delaunay triangulation and the nearest neighbor
				// search will succeed.  [NOTE about half_infinity_doublereal:
				// Something like std::numeric_limits<DoubleReal>::max() / 2.  does
				// not work here, because CGAL will fail on a precondition.  But 1E10
				// feels like a really big number, doesn't it?  Found by trial and
				// error.  Clemens, 2008-05-19]
				DoubleReal really_big_doublereal = 1E10;
				Feature outlier_feature_1;
				outlier_feature_1.setRT(-really_big_doublereal);
				outlier_feature_1.setMZ(-really_big_doublereal);
				ConsensusFeature outlier_consensusfeature_1( std::numeric_limits<UInt>::max(), std::numeric_limits<UInt>::max(), outlier_feature_1 );
				Point outlier_point_1( -really_big_doublereal, -really_big_doublereal, outlier_consensusfeature_1, getModelMap().size() );
				p_set.push_back(outlier_point_1);
				Feature outlier_feature_2;
				outlier_feature_2.setRT(really_big_doublereal);
				outlier_feature_2.setMZ(really_big_doublereal);
				ConsensusFeature outlier_consensusfeature_2( std::numeric_limits<UInt>::max(), std::numeric_limits<UInt>::max(), outlier_feature_2 );
				Point outlier_point_2( really_big_doublereal, really_big_doublereal, outlier_consensusfeature_2, getModelMap().size() );
				p_set.push_back(outlier_point_2);


				V_("p_set.number_of_vertices(): " << p_set.number_of_vertices() << "  [includes two dummy outliers]");
				V_("Transformation rt " << transformation_[RawDataPoint2D::RT]);
				V_("Transformation mz " << transformation_[RawDataPoint2D::MZ]);
				V_("balance_MZ_: " << balance_MZ_);

				// Initialize a hash map for the elements of model map to avoid that elements of the reference map occur in several element pairs
				// The semantics for the hashed values is:
				// -1: not touched,
				// -2: cannot assign unambiguously,
				// >=0: index of the matching scene ConsensusFeature
				std::vector<Int> match_model_to_scene( getModelMap().size(), -1 );

				// inverse of match_model_to_scene is match_scene_to_model
				std::vector<Int> match_scene_to_model( getSceneMap().size(), -1);

				UInt current_result_cf_index = 0;

				// TODO: The matching algorithm below is somewhat weird, but on the
				// other hand it does not rely on a Delaunay triangulation for the
				// scene map.  We should try a symmetric version and see if it
				// performs better, besides aesthetic advantages.  (Clemens,
				// 2008-05-19)

				// take each point in the scene map and search for its neighbours in the model map (within a given (transformed) range)
				for ( UInt scene_cf_index = 0; scene_cf_index < getSceneMap().size(); ++scene_cf_index )
				{
					// compute the transformed pos
					double rt_pos = getSceneMap()[scene_cf_index].getRT();
					transformation_[RawDataPoint2D::RT].apply(rt_pos); // TODO
					double mz_pos = getSceneMap()[scene_cf_index].getMZ();
					transformation_[RawDataPoint2D::MZ].apply(mz_pos); // TODO

					mz_pos *= balance_MZ_;
					Point transformed_pos(rt_pos,mz_pos);
					V_("Transformed Position is : " << transformed_pos );

					// Search for nearest and second nearest neighbor.
					std::vector< Vertex_handle > resulting_range;
					p_set.nearest_neighbors(transformed_pos,2,std::back_inserter(resulting_range));
					Point nearest        = resulting_range[0]->point();
					Point second_nearest = resulting_range[1]->point();

					V_("Neighbouring points : ");
					V_("nearest : " << nearest << " key=" << nearest.key << "\n" << *nearest.element );
					V_("second_nearest : " << second_nearest << " key=" << second_nearest.key << "\n" << *second_nearest.element );

					if (
							(
							 (fabs(transformed_pos[RawDataPoint2D::RT] - nearest.hx())  < precision_[RawDataPoint2D::RT]) &&
							 (fabs(transformed_pos[RawDataPoint2D::MZ] - nearest.hy())  < precision_[RawDataPoint2D::MZ])
							)
							&&
							(
							 (fabs(second_nearest.hx() - nearest.hx())  > max_pair_distance_[RawDataPoint2D::RT]) ||
							 (fabs(second_nearest.hy() - nearest.hy())  > max_pair_distance_[RawDataPoint2D::MZ])
							)
						 )
					{
						// if the element already part of a ElementPair the value in the match_model_to_scene becomes -2
						if ( match_model_to_scene[nearest.key] > -1)
						{
							match_model_to_scene[nearest.key] = -2;
						}
						// otherwise if the element is until now no part of a element pair,
						// set the value in the match_model_to_scene to the index of the pair in the result_map vector
						else
						{
							if ( match_model_to_scene[nearest.key] == -1)
							{
								match_model_to_scene[nearest.key] = scene_cf_index;
								match_scene_to_model[scene_cf_index] = nearest.key;
								maps_.result_->push_back(ConsensusFeature());
								if ( map_index_.model_ == -1 )
								{
									maps_.result_->back().insert( *nearest.element );
								}
								else
								{
										maps_.result_->back().insert( map_index_.model_, nearest.key, *nearest.element );
								}
								if ( map_index_.scene_ == -1 )
								{
									maps_.result_->back().insert( getSceneMap()[scene_cf_index] );
								}
								else
								{
									maps_.result_->back().insert( map_index_.scene_, scene_cf_index, getSceneMap()[scene_cf_index] );
								}
								maps_.result_->back().computeConsensus();
								V_("Result " << current_result_cf_index << " : " << maps_.result_->back());
								++current_result_cf_index;
							}
						}
					}
				}

				// restore match_model_to_scene from match_scene_to_model
				for ( UInt scene_cf_index = 0; scene_cf_index < getSceneMap().size(); ++scene_cf_index )
				{
					Int model_cf_index = match_scene_to_model[scene_cf_index];
					if ( model_cf_index >= 0 )
					{
						match_model_to_scene[ model_cf_index ] = scene_cf_index;
					}
				}

				// write out singleton consensus features for unmatched consensus features in model map
				for ( UInt model_cf_index = 0; model_cf_index < getModelMap().size(); ++model_cf_index )
				{
					if ( match_model_to_scene[model_cf_index] < 0 )
					{
						maps_.result_->push_back(ConsensusFeature());
						if ( map_index_.model_ == -1)
						{
							maps_.result_->back().insert( getModelMap()[model_cf_index] );
						}
						else
						{
							maps_.result_->back().insert( map_index_.model_, model_cf_index, getModelMap()[model_cf_index] );
						}
						maps_.result_->back().computeConsensus();
						V_("Result " << current_result_cf_index << " : " << maps_.result_->back());
						V_("match_model_to_scene[model_cf_index]" << match_model_to_scene[model_cf_index]);
						++current_result_cf_index;
					}
				}

				// write out singleton consensus features for unmatched consensus features in scene map
				for ( UInt scene_cf_index = 0; scene_cf_index < getSceneMap().size(); ++scene_cf_index )
				{
					if ( match_scene_to_model[scene_cf_index] < 0 )
					{
						maps_.result_->push_back(ConsensusFeature());
						if ( map_index_.scene_ )
						{
							maps_.result_->back().insert( getSceneMap()[scene_cf_index] );
						}
						else
						{
							maps_.result_->back().insert( map_index_.scene_, scene_cf_index, getSceneMap()[scene_cf_index] );
						}
						maps_.result_->back().computeConsensus();
						V_("Result " << current_result_cf_index << " : " << maps_.result_->back());
						++current_result_cf_index;
					}
				}

      }

      return;
    }


    /// The actual algorithm for finding element pairs.
    void findElementPairs()
    {
      const ConsensusMap& reference_map = *(element_map_[MODEL]);
      const ConsensusMap& transformed_map = *(element_map_[SCENE]);

      V_("@@@ findElementPairs_()");

      UInt n = reference_map.size();

      // Vector to fill the point set for triangulation
      // Penalize a deviation in mz more than in rt: deviation(diff_intercept_[RawDataPoint2D::RT]) ~ deviation(diff_intercept_[RawDataPoint2D::MZ])
      std::vector< Point > positions_reference_map;
      for (UInt i = 0; i < n; ++i)
      {
        positions_reference_map.push_back( Point( reference_map[i].getRT(), reference_map[i].getMZ() * balance_MZ_, reference_map[i], i) );
      }
			
      // compute the delaunay triangulation
      Point_set_2 p_set(positions_reference_map.begin(),positions_reference_map.end());

      V_("Translation rt " << transformation_[RawDataPoint2D::RT]);
      V_("Translation mz " << transformation_[RawDataPoint2D::MZ]);

      // Initialize a hash map for the elements of reference_map to avoid that elements of the reference map occur in several element pairs
      std::vector< Int > lookup_table(n,-1);
      std::vector< std::pair< const ConsensusFeature*,const ConsensusFeature*> > all_element_pairs;

      UInt index_act_element_pair = 0;
      // take each point in the first data map and search for its neighbours in the second element map (within a given (transformed) range)
      for ( UInt fi1 = 0; fi1 < transformed_map.size(); ++fi1 )
      {
        // compute the transformed iso-rectangle (upper_left,bottom_left,bottom_right,upper_right) for the range query
        double rt_pos = transformed_map[fi1].getRT();
        double mz_pos = transformed_map[fi1].getMZ();

        V_("Search for two nearest neighbours of " << rt_pos << ' ' << transformed_map[fi1].getMZ() );
        transformation_[RawDataPoint2D::RT].apply(rt_pos);
        transformation_[RawDataPoint2D::MZ].apply(mz_pos);

        mz_pos /= (diff_intercept_[RawDataPoint2D::MZ] / diff_intercept_[RawDataPoint2D::RT]);
        Point transformed_pos(rt_pos,mz_pos);

        V_("Transformed Position is : " << transformed_pos );

        std::vector< Vertex_handle > resulting_range;
        p_set.nearest_neighbors(transformed_pos,2,std::back_inserter(resulting_range));

        V_("Neighbouring points : ");
        for (std::vector< Vertex_handle >::const_iterator it = resulting_range.begin(); it != resulting_range.end(); it++)
        {
          V_((*it)->point());
          V_(*((*it)->point().element));
        }

        // if the first neighbour is close enough to act_pos and the second_nearest neighbour lies far enough from the nearest neighbour
        Point nearest = resulting_range[0]->point();
        Point second_nearest = resulting_range[1]->point();

        if (((fabs(transformed_pos[RawDataPoint2D::RT] - nearest.hx())  < precision_[RawDataPoint2D::RT])
             &&  (fabs(transformed_pos[RawDataPoint2D::MZ] - nearest.hy())  < precision_[RawDataPoint2D::MZ]))
            && ((fabs(second_nearest.hx() - nearest.hx())  > max_pair_distance_[RawDataPoint2D::RT])
                || (fabs(second_nearest.hy() - nearest.hy())  > max_pair_distance_[RawDataPoint2D::MZ])))
        {
          all_element_pairs.push_back(std::pair<const ConsensusFeature*,const ConsensusFeature*>(nearest.element,&transformed_map[fi1]));

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

      // TODO improve running time (linear search!)
      for (UInt i = 0; i < n; ++i)
      {
        Int pair_key = lookup_table[i];
        if ( pair_key > -1 )
        {
          element_pairs_->push_back(ElementPairType(*(all_element_pairs[pair_key].second),*(all_element_pairs[pair_key].first)));
        }
      }

    }

    /**@brief The actual algorithm for finding consensus elements.  Elements
    in the first_map are aligned to elements in the second_map, so the
    second_map contains the resulting consensus elements.
    */
    void computeConsensusMap(const ConsensusMap& first_map, ConsensusMap& second_map)
    {
      V_("@@@ computeConsensusMap()");

      // Vector to fill the point set for triangulation
      std::vector< Point > positions_reference_map;
      UInt n = first_map.size();
      for (UInt i = 0; i < n; ++i)
      {
        positions_reference_map.push_back( Point(first_map[i].getRT(), first_map[i].getMZ() * balance_MZ_, first_map[i], i) );
      }     

      StopWatch timer;
      V_("Start delaunay triangulation for " << positions_reference_map.size() << " elements");
      // compute the delaunay triangulation
      timer.start();
      Point_set_2 p_set(positions_reference_map.begin(),positions_reference_map.end());
      timer.stop();
      V_("End delaunay triangulation after " << timer.getCPUTime() << "s");

      // Initialize a hash map for the elements of reference_map to avoid that elements of the reference map occur in several element pairs
      std::vector< Int > lookup_table(n,-1);
      std::vector< std::pair< const ConsensusFeature*, ConsensusFeature*> > all_element_pairs;

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

        V_("Search for two nearest neighbours of " << rt_pos << ' ' << second_map[fi1].getMZ() << ' ' << second_map[fi1].getIntensity() );
        Point transformed_pos(rt_pos,mz_pos,second_map[fi1],fi1);

        V_("Transformed Position is : " << transformed_pos );
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
            all_element_pairs.push_back(std::pair<const ConsensusFeature*,ConsensusFeature*>(nearest.element,&(second_map[fi1])));
          }
        }
        else
          if (resulting_range.size() > 1)
          {
            nearest = resulting_range[0]->point();
            second_nearest = resulting_range[1]->point();
            V_("Nearest: " << (nearest.element)->getRT() << ' ' << fabs(transformed_pos[RawDataPoint2D::RT] - nearest.hx()) << ' ' 
							 << (nearest.element)->getMZ() << ' ' << fabs(transformed_pos[RawDataPoint2D::MZ] - nearest.hy()) << ' '
							 << (nearest.element)->getIntensity());
						
            V_("Second nearest: " << (second_nearest.element)->getRT() << ' ' << second_nearest.hx() << ' ' 
							 << (second_nearest.element)->getMZ() << ' ' << second_nearest.hy() << ' '
							 << (second_nearest.element)->getIntensity());
						
            if (
								(
								 (fabs(transformed_pos[RawDataPoint2D::RT] - nearest.hx())  < precision_[RawDataPoint2D::RT]) &&
								 (fabs(transformed_pos[RawDataPoint2D::MZ] - nearest.hy())  < precision_[RawDataPoint2D::MZ])
								)
                &&
								(
								 (fabs(second_nearest.hx() - nearest.hx())  > max_pair_distance_[RawDataPoint2D::RT]) ||
								 (fabs(second_nearest.hy() - nearest.hy())  > max_pair_distance_[RawDataPoint2D::MZ])
								)
							 )
            {
              all_element_pairs.push_back(std::pair<const ConsensusFeature*,ConsensusFeature*>(nearest.element,&(second_map[fi1])));
              V_("Push first: " << *(nearest.element));
							V_("Push second: " << second_map[fi1]);
							
							Int element_key = resulting_range[0]->point().key;

              // if the element a is already part of a ElementPair (a,b) do:
              //    if (the element c closer to a than b to a) and (the distance between c and b is > a given threshold) do:
              //    --> push (a,c)
              //    else
              //    --> the value in the lookup_table becomes -2 because the mapping is not unique
              if ( lookup_table[element_key] > -1)
              {
                Int pair_key = lookup_table[element_key];
                const ConsensusFeature& first_map_a = *(all_element_pairs[pair_key].first);
                ConsensusFeature& second_map_b = *(all_element_pairs[pair_key].second);
                ConsensusFeature& second_map_c = second_map[fi1];

                V_("The element " << first_map_a.getPosition() << ' ' << first_map_a.getIntensity() << " has two element partners \n");
                V_(second_map_b.getPosition() << ' ' << second_map_b.getIntensity() 
									 << "  and  " << second_map_c.getPosition() << ' ' << second_map_c.getIntensity()) ;

                V_("Range " << second_map_b.getPositionRange() << "  and  " << second_map_c.getPositionRange());

                if (second_map_c.getPositionRange().encloses(first_map_a.getPosition())
                    && !second_map_b.getPositionRange().encloses(first_map_a.getPosition()))
                {
                  lookup_table[element_key] = index_act_element_pair;
                  V_(second_map_c.getPosition() << " and " << first_map_a.getPosition() << " are a pair");
                }
                else
                {
                  // if second_map_b and second_map_c do not enclose first_map_a
                  if (!(second_map_b.getPositionRange().encloses(first_map_a.getPosition())
                        && !second_map_c.getPositionRange().encloses(first_map_a.getPosition())))
                  {
                    V_(second_map_b.getPosition() << " and " << first_map_a.getPosition() << " are a pair, but check the distance between c and b");
                    // check the distance between second_map_b and second_map_c
                    if (fabs(second_map_b.getMZ() / (diff_intercept_[RawDataPoint2D::MZ]/diff_intercept_[RawDataPoint2D::RT])
                             - second_map_c.getMZ() / (diff_intercept_[RawDataPoint2D::MZ]/diff_intercept_[RawDataPoint2D::RT]))
                        > max_pair_distance_[RawDataPoint2D::MZ])
                    {
                      V_("distance ok");
                      // and check which one of the elements lies closer to first_map_
                      if( sqrt(pow((first_map_a.getRT() - second_map_b.getRT()), 2)
                               + pow((first_map_a.getMZ() - second_map_b.getMZ()), 2))
                          > sqrt(pow((first_map_a.getRT() - second_map_c.getRT()), 2)
                                 + pow((first_map_a.getMZ() - second_map_c.getMZ()), 2)))
                      {
                        lookup_table[element_key] = index_act_element_pair;
                        V_("take a and c");
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
      V_("Insert elements ");
      std::vector< const ConsensusFeature* > single_elements_first_map;
      for (UInt i = 0; i < n; ++i)
      {
        Int pair_key = lookup_table[i];
        if ( pair_key > -1 )
        {
          FeatureHandle index_tuple(*((all_element_pairs[pair_key].first)->begin()));
					(all_element_pairs[pair_key].second)->insert(index_tuple);
          (all_element_pairs[pair_key].second)->computeConsensus();
          ++pairs;
        }
        // add a singleton consensus element
        else
        {
          single_elements_first_map.push_back(positions_reference_map[i].element);
          ++ref_single;
        }
      }

      for (UInt i = 0; i < single_elements_first_map.size(); ++i)
      {
        second_map.push_back(*(single_elements_first_map[i]));
        second_map.back().computeConsensus();
      }

      V_("SINGLE TRANS: " << trans_single);
      V_("SINGLE REF: " << ref_single);
      V_("PAIRS: " << pairs);

    }


   protected:
    virtual void updateMembers_()
    {
      diff_intercept_[RawDataPoint2D::RT] = (double)param_.getValue("similarity:diff_intercept:RT");
      diff_intercept_[RawDataPoint2D::MZ] = (double)param_.getValue("similarity:diff_intercept:MZ");

      balance_MZ_ = diff_intercept_[RawDataPoint2D::RT] / diff_intercept_[RawDataPoint2D::MZ];

      max_pair_distance_[RawDataPoint2D::RT] = (float)param_.getValue("similarity:max_pair_distance:RT");
      max_pair_distance_[RawDataPoint2D::MZ] = (float)param_.getValue("similarity:max_pair_distance:MZ") * balance_MZ_;

      precision_[RawDataPoint2D::RT] = (float)param_.getValue("similarity:precision:RT");
      precision_[RawDataPoint2D::MZ] = (float)param_.getValue("similarity:precision:MZ") * balance_MZ_;

      return;
    }
		
    // TODO replace diff_intercept_ by balance_MZ_
    /// A parameter for similarity_().
    DoubleReal diff_intercept_[2];

    /// Factor by which MZ has to be rescaled so that differences in MZ and RT are equally significant.
    DoubleReal balance_MZ_;

    /// To uniquely assign an element e1 of the scene map to another element e2 in the model map
    /// all elements in the scene map have to lie at least max_pair_distance_ far from e1 and
    /// all elements in the model map have to lie at least max_pair_distance_ far from e2.
    float max_pair_distance_[2];
    /// Only points that differ not more than precision_ can be assigned as a pair
    float precision_[2];
  }
    ; // DelaunayPairFinder
} // namespace OpenMS

#undef V_

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H
