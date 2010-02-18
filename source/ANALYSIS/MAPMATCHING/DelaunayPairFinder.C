// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl, Eva Lange $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/DelaunayPairFinder.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>

#include <climits> // damn CGAL uses INT_MAX!

#ifdef _MSC_VER // disable some CGAL warnings that distract from ours
#	pragma warning( push ) // save warning state
#	pragma warning( disable : 4396 )
#	pragma warning( disable : 4267 )
#endif
#include <CGAL/Cartesian.h>
#include <CGAL/Point_set_2.h>
#ifdef _MSC_VER
#	pragma warning( pop )  // restore old warning state
#endif


// #define Debug_DelaunayPairFinder
#ifdef Debug_DelaunayPairFinder
#define V_(bla) std::cout << __FILE__ ":" << __LINE__ << ": " << bla << std::endl;
#else
#define V_(bla)
#endif
#define VV_(bla) V_(""#bla": " << bla)

namespace OpenMS
{
	/**
  	@brief Nested class, which inherits from the CGAL Point_2 class and
  	additionally contains a reference to the corresponding element and a
  	unique key

  	@todo check which ctors are really needed and make sense (Clemens)
  */
  struct DelaunayPairFinder::Point: public CGAL::Point_2< CGAL::Cartesian<double> >
  {
    typedef CGAL::Point_2< CGAL::Cartesian<double> > Base;

    const ConsensusFeature* element;
    SignedSize key;

    /// Default ctor
    inline Point() : Base(), element(0), key(0) {}

    /// Ctor from Base class, aka CGAL:Point_2<...>
    inline Point(const Base& cgal_point) : Base(cgal_point), element(0), key(0) {}

    /// Ctor from coordinates, element, and key
    inline Point(double hx, double hy, const ConsensusFeature& element, SignedSize key) : Base(hx,hy), element(&element), key(key) {}

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
  class  DelaunayPairFinder::GeometricTraits
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

	/**
		@brief Just a silly array of Point of size two which has constructors
		missing in Point[2].  (std::vector<Point[2]> doesn't work.)
	*/
	struct DelaunayPairFinder::PointArray2
	{
		/// default ctor
		PointArray2() {array_[0] = Point(); array_[1] = Point();}
		/// ctor from two Ints
		PointArray2(Point const i0, Point const i1 ) {array_[0] = i0; array_[1] = i1;}
		/// copy ctor
		PointArray2(PointArray2 const& rhs) {array_[0] = rhs.array_[0]; array_[1] = rhs.array_[1];}
		/// assignment op
		PointArray2 & operator=(PointArray2 const& rhs) {array_[0] = rhs.array_[0]; array_[1] = rhs.array_[1]; return *this;}
		/// indexing op
		Point & operator[](UInt const index) {return array_[index];}
		/// indexing op
		Point const& operator[](UInt const index) const {return array_[index];}

	 protected:
		/// the underlying array
		Point array_[2];
	};

  DelaunayPairFinder::DelaunayPairFinder()
  	: Base()
  {
    //set the name for DefaultParamHandler error messages
    Base::setName(getProductName());

    defaults_.setValue("similarity:max_pair_distance:RT", 20.0, "Maximal allowed distance in retention time for a pair to be matched");
    defaults_.setValue("similarity:max_pair_distance:MZ",  1.0, "Maximal allowed distance in mass-to-charge for a pair to be matched");
    defaults_.setValue("similarity:second_nearest_gap",    0.0, "The distance of the second nearest neighbors must be this factor larger than to the nearest");

    Base::defaultsToParam_();
  }

	void DelaunayPairFinder::updateMembers_()
  {
		max_pair_distance_[MZ] = (DoubleReal) param_.getValue("similarity:max_pair_distance:MZ");
		max_pair_distance_[RT] = (DoubleReal) param_.getValue("similarity:max_pair_distance:RT");
		second_nearest_gap_    = (DoubleReal) param_.getValue("similarity:second_nearest_gap");
		internal_mz_scaling_   = max_pair_distance_[RT] / max_pair_distance_[MZ];
		max_squared_distance_  = pow(max_pair_distance_[RT],2);

		V_("@@@ DelaunayPairFinder::updateMembers_()");
		VV_(max_pair_distance_[RT]);
		VV_(max_pair_distance_[MZ]);
		VV_(internal_mz_scaling_);
		VV_(max_squared_distance_);
		VV_(second_nearest_gap_);
  }

  void DelaunayPairFinder::run(const std::vector<ConsensusMap>& input_maps, ConsensusMap &result_map)
  {
  	if (input_maps.size()!=2) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"exactly two input maps required");
		checkIds_(input_maps);

  	typedef CGAL::Point_set_2< DelaunayPairFinder::GeometricTraits, CGAL::Triangulation_data_structure_2< CGAL::Triangulation_vertex_base_2< DelaunayPairFinder::GeometricTraits > > > Point_set_2;
  	typedef Point_set_2::Vertex_handle Vertex_handle;

  	const ConsensusMap* maps_array[2];
  	maps_array[MODEL_] = &(input_maps[0]);
		maps_array[SCENE_] = &(input_maps[1]);

  	startProgress(0,10,"Delaunay pair finder");
		UInt progress = 0;
		setProgress(++progress);

    V_("@@@ DelaunayPairFinder::run()");
		VV_(max_pair_distance_[RT]);
		VV_(max_pair_distance_[MZ]);
		VV_(internal_mz_scaling_);
		VV_(max_squared_distance_);
		VV_(second_nearest_gap_);

		// The delaunay triangulation data structures for model and scene.
		Point_set_2 p_set[2];

		// We will add two outlier points to each p_set so that we will always
		// have at least three points in the Delaunay triangulation and the
		// nearest neighbor search will succeed.  [NOTE about
		// really_big_doublereal: Something like
		// std::numeric_limits<DoubleReal>::max() / 2.  does not work here,
		// because CGAL fails on a precondition.  But 1E10 is a really big
		// number, isn't it?  Found by trial and error.  Clemens, 2008-05-19]
		DoubleReal const really_big_doublereal = 1E10;
		Feature outlier_feature_1;
		outlier_feature_1.setRT(-really_big_doublereal); outlier_feature_1.setMZ(-really_big_doublereal);
		ConsensusFeature const outlier_consensusfeature_1( std::numeric_limits<UInt>::max(), outlier_feature_1, std::numeric_limits<UInt>::max() );
		Point const outlier_point_1( -really_big_doublereal, -really_big_doublereal, outlier_consensusfeature_1, -1 );
		Feature outlier_feature_2;
		outlier_feature_2.setRT(really_big_doublereal); outlier_feature_2.setMZ(really_big_doublereal);
		ConsensusFeature const outlier_consensusfeature_2( std::numeric_limits<UInt>::max(), outlier_feature_2, std::numeric_limits<UInt>::max() );
		Point const outlier_point_2( really_big_doublereal, really_big_doublereal, outlier_consensusfeature_2, -1 );
		PointArray2 const outlier_points( outlier_point_1, outlier_point_2 );

		// do the preprocessing for both input maps
    for ( UInt input = MODEL_; input <= SCENE_; ++ input )
    {
			setProgress(++progress);

			// TODO Find out whether it is (1) correct and (2) fast if we
			// push_back() the Points into the Delaunay triangulation. Otherwise,
			// use an iterator adapter and construct Point_set_2 p_set from an
			// iterator range.
			if ( input == MODEL_ )
			{
				for (Size i = 0; i < input_maps[0].size(); ++i)
				{
					DoubleReal trans_rt = input_maps[0][i].getRT();
					DoubleReal trans_mz = input_maps[0][i].getMZ() * internal_mz_scaling_;
					p_set[MODEL_].push_back( Point( trans_rt, trans_mz, input_maps[0][i], i) );
					V_("MODEL_: trans_rt:"<<trans_rt<<" trans_mz:"<<trans_mz);
				}
			}
			else // input == SCENE_
			{
				for (Size i = 0; i < input_maps[1].size(); ++i)
				{
					DoubleReal trans_rt = input_maps[1][i].getRT();
					//TODO: use offset -- transformation_[Peak2D::RT].apply(trans_rt);
					DoubleReal trans_mz = input_maps[1][i].getMZ() * internal_mz_scaling_;
					p_set[SCENE_].push_back( Point( trans_rt, trans_mz, input_maps[1][i], i) );
					V_("SCENE_: trans_rt:"<<trans_rt<<" trans_mz:"<<trans_mz);
				}
			}

			p_set[input].push_back(outlier_point_1);
			p_set[input].push_back(outlier_point_2);
			V_("p_set[" << input << "].number_of_vertices(): " << p_set[input].number_of_vertices() << "  [includes two dummy outliers]");
    }

    Size p_set_number_of_vertices[2];
    for ( UInt input = MODEL_; input <= SCENE_; ++ input )
    {
        p_set_number_of_vertices[input] = p_set[input].number_of_vertices();
    }

    // Empty output destination
    result_map.clear(true);

		// In this we store the best ([0]) and second best ([1]) neighbours;
		// essentially this is a regular bipartite graph of degree two.
		// E.g. neighbours[MODEL_][i][j] is the j-th best match (index of an
		// element of scene) for the i-th consensus feature in the model.
		std::vector<PointArray2> neighbours[2];

		// For each input, we find nearest and second nearest neighbors in the
		// other map.
		std::vector<Vertex_handle> neighbors_buffer;
    for ( UInt input = MODEL_; input <= SCENE_; ++ input )
    {
			setProgress(++progress);

			UInt const other_input = 1 - input;

			neighbours[input].resize( maps_array[input]->size(), outlier_points );
      Size loop_count = 0;
			for ( Point_set_2::Point_iterator iter = p_set[input].points_begin(); iter != p_set[input].points_end(); ++iter )
			{
        VV_(iter->key);
        if ( loop_count++ == p_set_number_of_vertices[input] ) break; // FIXME this is a hack introduced for Mac using gcc 4.0.1 (Clemens)
        if ( iter->key == -1 ) continue;
				neighbors_buffer.clear();
				p_set[other_input].nearest_neighbors( *iter, 2, std::back_inserter(neighbors_buffer) );
				neighbours[input][iter->key][0] = neighbors_buffer[0]->point();
				neighbours[input][iter->key][1] = neighbors_buffer[1]->point();
				VV_(neighbours[input][iter->key][0].key);
				VV_(neighbours[input][iter->key][1].key);
				V_( iter->x() << " " << iter->y() << " " <<
						neighbours[input][iter->key][0].x() << " " << neighbours[input][iter->key][0].y() << " " <<
						neighbours[input][iter->key][1].x() << " " << neighbours[input][iter->key][1].y() );
			}
		}

    // Initialize a hash map for the elements of model map to avoid that
    // elements of the reference map occur in several element pairs
		//
    // The semantics for the hashed values is:
    // -1: not touched,
    // -2: cannot assign unambiguously, (currently not being used!)
    // >=0: index of the matching scene ConsensusFeature
    std::vector<Int> matches[2];
		matches[MODEL_].resize(input_maps[0].size(),-1);
		matches[SCENE_].resize(input_maps[1].size(),-1);

		V_("max_pair_distance_[RT]: " << max_pair_distance_[RT]);
		V_("max_pair_distance_[MZ]: " << max_pair_distance_[MZ]);
		V_("max_squared_distance_:  " << max_squared_distance_ );
		V_("second_nearest_gap:     " << second_nearest_gap_   );

    UInt current_result_cf_index = 0;

		setProgress(++progress);

    // take each point in the model map and search for its neighbours in the scene map

    Size loop_count = 0;
		for ( Point_set_2::Point_iterator model_iter = p_set[MODEL_].points_begin(); model_iter != p_set[MODEL_].points_end(); ++model_iter )
    {

      if ( loop_count++ == p_set_number_of_vertices[MODEL_] ) break; // FIXME this is a hack introduced for Mac using gcc 4.0.1 (Clemens)

			Int const model_cf_index = model_iter->key;
			VV_(model_cf_index);
			if ( model_cf_index == -1 ) continue;

			Point const & scene_point = neighbours[MODEL_][model_cf_index][0];
			Int const scene_cf_index = scene_point.key;
			VV_(scene_cf_index);
			if ( scene_cf_index == -1 ) continue;

			bool const is_bidirectionally_nearest = neighbours[SCENE_][ scene_cf_index ][0].key == model_cf_index;
			VV_(is_bidirectionally_nearest);

			if ( is_bidirectionally_nearest )
			{
				DoubleReal pair_distance_m_s =
					squaredDistance_( model_iter->x(), model_iter->y(),
														 scene_point.x(), scene_point.y()  );

				VV_(pair_distance_m_s);
				bool const is_close_enough = ( pair_distance_m_s <= max_squared_distance_ );
				VV_(is_close_enough);

				if ( is_close_enough ) /* && is_bidirectionally_nearest */
				{
					Point const & scene_point_second_nearest = neighbours[MODEL_][model_cf_index][1];
					DoubleReal pair_distance_m_s2nd =
						squaredDistance_( model_iter->x(),                model_iter->y(),
															scene_point_second_nearest.x(), scene_point_second_nearest.y() );
					VV_(pair_distance_m_s2nd);
					Point const & model_point_second_nearest = neighbours[SCENE_][scene_cf_index][1];
					VV_(model_point_second_nearest);
					DoubleReal pair_distance_m2nd_s =
						squaredDistance_( model_point_second_nearest.x(), model_point_second_nearest.y(),
															scene_point.x(),                scene_point.y()                 );
					VV_(pair_distance_m2nd_s);
					DoubleReal min_second_pair_distance = pair_distance_m_s * second_nearest_gap_;
					VV_(min_second_pair_distance);
					bool const is_unambiguous =
						( pair_distance_m_s2nd >= min_second_pair_distance &&
							pair_distance_m2nd_s >= min_second_pair_distance    );
					VV_(is_unambiguous);
					if ( is_unambiguous ) /* && is_close_enough && is_bidirectionally_nearest */
					{
						// assign matching pair
						matches[MODEL_][model_cf_index] = scene_cf_index;
						matches[SCENE_][scene_cf_index] = model_cf_index;

						VV_(matches[MODEL_][model_cf_index]);
						VV_(matches[SCENE_][scene_cf_index]);
						VV_(current_result_cf_index);

						// create a consensus feature
						result_map.push_back(ConsensusFeature());
						result_map.back().insert( input_maps[1][scene_cf_index] );
						result_map.back().getPeptideIdentifications().insert(result_map.back().getPeptideIdentifications().end(),input_maps[1][scene_cf_index].getPeptideIdentifications().begin(),input_maps[1][scene_cf_index].getPeptideIdentifications().end());
						result_map.back().insert( input_maps[0][model_cf_index] );
						result_map.back().getPeptideIdentifications().insert(result_map.back().getPeptideIdentifications().end(),input_maps[0][model_cf_index].getPeptideIdentifications().begin(),input_maps[0][model_cf_index].getPeptideIdentifications().end());

						result_map.back().computeConsensus();
						V_("Result " << current_result_cf_index << " : " << result_map.back());
						++current_result_cf_index;
					}
				}
			}
    }

		setProgress(++progress);

    // write out unmatched consensus features in model and scene
    for ( UInt input = MODEL_; input <= SCENE_; ++ input )
    {
			for ( UInt index = 0; index < maps_array[input]->size(); ++ index )
			{
				if ( matches[input][index] < 0 )
				{
					result_map.push_back(ConsensusFeature());
					result_map.back().insert( (*maps_array[input])[index] );
					result_map.back().getPeptideIdentifications().insert(result_map.back().getPeptideIdentifications().end(),(*maps_array[input])[index].getPeptideIdentifications().begin(),(*maps_array[input])[index].getPeptideIdentifications().end());
					result_map.back().computeConsensus();
					V_("Result " << current_result_cf_index << " : " << result_map.back());
					V_("matches["<<input<<"]["<<index<< "]: " << matches[input][index] );
					++current_result_cf_index;
				}
			}
		}

		setProgress(++progress);

		// Acquire statistics about distances
		// TODO: use these statistics to derive aposteriori estimates for optimal matching parameters.
		if ( 1 )
		{
			DoubleReal squared_dist_RT = 0;
			DoubleReal squared_dist_MZ = 0;
			UInt count = 0;
			for ( Int model_feature_index = 0; model_feature_index < (Int) input_maps[0].size(); ++ model_feature_index )
			{
				Int scene_feature_index = matches[0][model_feature_index];
				V_("model_feature_index:" << model_feature_index << "  scene_feature_index: " << scene_feature_index);
				if ( scene_feature_index >= 0 )
				{
					++ count;
					squared_dist_RT += pow( input_maps[1][scene_feature_index].getRT() - input_maps[0][model_feature_index].getRT(), 2);
					// VV_(pow( input_maps[1][scene_feature_index].getRT() - input_maps[0][model_feature_index].getRT(), 2));
					squared_dist_MZ += pow( input_maps[1][scene_feature_index].getMZ() - input_maps[0][model_feature_index].getMZ(), 2);
					// VV_(pow( input_maps[1][scene_feature_index].getMZ() - input_maps[0][model_feature_index].getMZ(), 2));
				}
			}
			if ( count )
			{
				squared_dist_RT /= count;
				squared_dist_MZ /= count;
			}
			DoubleReal avg_dist_RT = sqrt(squared_dist_RT); avg_dist_RT+=0; //TODO
			DoubleReal avg_dist_MZ = sqrt(squared_dist_MZ); avg_dist_MZ+=0;
			VV_(count);
			VV_(avg_dist_RT);
			VV_(avg_dist_MZ);
		}

		setProgress(++progress);

		// Very useful for checking the results, and the ids have no real meaning anyway
		result_map.sortByMZ();

		//Add protein identifications to result map
		result_map.getProteinIdentifications().insert(result_map.getProteinIdentifications().end(),input_maps[1].getProteinIdentifications().begin(), input_maps[1].getProteinIdentifications().end());

		//Add unassigned peptide identifications to result map
		result_map.getUnassignedPeptideIdentifications().insert(result_map.getUnassignedPeptideIdentifications().end(),input_maps[1].getUnassignedPeptideIdentifications().begin(), input_maps[1].getUnassignedPeptideIdentifications().end());

		endProgress();

		return;
  }

}
