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
  we have to transform the elements MZ position m into a new MZ position m'= m * internal_mz_scaling;
  E.g. given internal_mz_scaling=10 results in 1s difference in RT being similar to 0.1Th difference in MZ.
			 
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
      defaults_.setValue("similarity:internal_mz_scaling",100.0,"Factor by which MZ has to be rescaled so that differences in MZ and RT are equally significant");
      defaults_.setValue("similarity:internal_mz_scaling",100.0,"Factor by which MZ has to be rescaled so that differences in MZ and RT are equally significant");

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

    /**@brief Nested class, which inherits from the CGAL Point_2 class and
		additionally contains a reference to the corresponding element and a
		unique key

		@todo check which ctors are really needed and make sense
		*/
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
			

			// TODO Find out whether it is (1) correct and (2) fast if we
			// push_back() the Points into the Delaunay triangulation. Otherwise,
			// use an iterator adapter and construct Point_set_2 p_set from
			// iterator range.
			Point_set_2 p_set;

			for (UInt i = 0; i < getModelMap().size(); ++i)
			{
				// incrementally compute the delaunay triangulation
				p_set.push_back( Point( getModelMap()[i].getRT(), getModelMap()[i].getMZ() * internal_mz_scaling_, getModelMap()[i], i) );
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
			V_("internal_mz_scaling_: " << internal_mz_scaling_);

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

				mz_pos *= internal_mz_scaling_;
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
      return;
    }

   protected:
    virtual void updateMembers_()
    {
      internal_mz_scaling_ = (DoubleReal)param_.getValue("similarity:internal_mz_scaling");

      max_pair_distance_[RawDataPoint2D::RT] = (DoubleReal)param_.getValue("similarity:max_pair_distance:RT");
      max_pair_distance_[RawDataPoint2D::MZ] = (DoubleReal)param_.getValue("similarity:max_pair_distance:MZ") * internal_mz_scaling_;

      precision_[RawDataPoint2D::RT] = (DoubleReal)param_.getValue("similarity:precision:RT");
      precision_[RawDataPoint2D::MZ] = (DoubleReal)param_.getValue("similarity:precision:MZ") * internal_mz_scaling_;

      return;
    }
		
    /// Factor by which MZ has to be rescaled so that differences in MZ and RT are equally significant.
    DoubleReal internal_mz_scaling_;

    /// To uniquely assign an element e1 of the scene map to another element e2 in the model map
    /// all elements in the scene map have to lie at least max_pair_distance_ far from e1 and
    /// all elements in the model map have to lie at least max_pair_distance_ far from e2.
    DoubleReal max_pair_distance_[2];

    /// Only points that differ not more than precision_ can be assigned as a pair
    DoubleReal precision_[2];
  }
    ; // DelaunayPairFinder
} // namespace OpenMS

#undef V_

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DELAUNAYPAIRFINDER_H
