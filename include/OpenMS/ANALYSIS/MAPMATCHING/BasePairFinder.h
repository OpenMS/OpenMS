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

#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <utility>
#include <fstream>

namespace OpenMS
{
 
  /**
	@brief The base class of all element pair finding algorithms.
		
	This class defines the basic interface for all element pair finding
	algorithms. It works on two element maps (FeatureMap is the default map
	type) and a transformation defined for the second element map (if no
	transformation is given, the pairs are found in the two original maps).  A
	element can be a DPeak, a DFeature or ConsensusFeature (wheras DFeature is
	the default element type).
		    
	@todo Remove transformation stuff (Marc, Clemens)
	
	@todo Redefine interface: see code in class docu (Marc)
	
	@code
	virtual void setModelMap((Int map_index, const ConsensusMap& map)
	virtual void setSceneMap((Int map_index, const ConsensusMap& map)
	void run(ConsensusMap& result_map);
	template <InputMapType> static void convert(UInt input_map_index, const InputMapType& input, ConsensusMap& output);
	@endcode
  	
  */
  // template < typename MapT = FeatureMap< > >
  class BasePairFinder : public FactoryProduct
  {
  public:

    /// Type of elements considered here
    // typedef ConsensusMap::value_type PointType;

    /// Position
    typedef DPosition<2> PositionType;

    /// Quality
    typedef DoubleReal QualityType;

    //// Intensity
    typedef DoubleReal IntensityType;

    /// Type of element pairs
    typedef ElementPair < ConsensusFeature > ElementPairType;

    /// Container for generated element pairs
    typedef std::vector < ElementPairType > ElementPairVectorType;

    /// Type of estimated transformation
    typedef LinearMapping TransformationType;

    /// Default constructor
    BasePairFinder();

		/// Destructor
    virtual ~BasePairFinder();

		/** Set model map (might do some preprocessing).

		@param map_index If map_index>=0, then run() will use it as a map index
		for the input and store feature handles pointing to the consensus features
		of theq input to the consensus features in the result, which adds another
		level of indirection/nesting.<br> If -1, then run() will "unpack" the
		consensus features from the input, i.e. it will store the feature handles
		contained in the consensus features rather than the consensus features
		themselves in the result.

		@param model_map Consensus map to be used as model.
		*/
		virtual void setModelMap(Int map_index, ConsensusMap const& model_map)
		{
			maps_.model_ = &model_map;
		// Note this is a virtual method. Hence it is strictly forbidden to
		// provide a default value for unpack!  I will bite off your fingers!  We
		// will get a big mess if the default is set differently in derived
		// classes. Clemens 2008-05-18
			map_index_.model_ = map_index;
		}

		/// Get model map
		virtual ConsensusMap const & getModelMap() const
		{
			return *maps_.model_;
		}

		/// Set scene map.  @sa setModelMap()
		virtual void setSceneMap(Int map_index, ConsensusMap const& scene_map)
		{
			maps_.scene_ = &scene_map;
			map_index_.scene_ = map_index;
		}

		/// Get scene map
		virtual ConsensusMap const & getSceneMap() const
		{
			return *maps_.scene_;
		}

		/// Run the algorithm
		virtual void run(ConsensusMap& result_map)
		{
			// Every derived class should set maps_.result_ at the beginning.
			maps_.result_ = &result_map;
			return;
		};

		/**@brief Convert any (random access) container of features to a ConsensusMap.  Each
		ConsensusFeature contains a map index, so this has to be given as well.
		The previous content of output_map is cleared.

		@param input_map_index The index of the input map.
		@param input_map The container to be converted.  (Must support size() and operator[].)
		@param output_map The resulting ConsensusMap.

		*/
		template <typename ContainerT>
		static void convert( UInt const input_map_index, ContainerT const & input_map, ConsensusMap& output_map )
		{
			output_map.clear();
			output_map.reserve(input_map.size());
			for ( UInt element_index = 0; element_index < input_map.size(); ++element_index )
			{
				output_map.push_back( ConsensusFeature( input_map_index, element_index, input_map[element_index] ) );
			}
			return;
		}
		
		/**@brief Similar to convert, but copies only the @p n most intense
		elements from an MSExperiment.

		@param input_map_index The index of the input map.
		@param input_map The input map to be converted.
		@param output_map The resulting ConsensusMap.
		@param n The maximum number of elements to be copied.
		*/
		static void convert( UInt const input_map_index, MSExperiment<> & input_map, ConsensusMap& output_map, UInt n ) // TODO find out what goes wrong in template instantiation (?!!)
		// template <typename PeakT, typename AllocT>
		// static void convert( UInt const input_map_index, MSExperiment<PeakT,AllocT> & input_map, ConsensusMap& output_map, UInt n )
		{
			input_map.updateRanges(1);
			if ( n > input_map.getSize() )
			{
				n = input_map.getSize();
			}
			output_map.clear();
			output_map.reserve(n);
			std::vector<RawDataPoint2D> tmp; // TODO let's see if this will pass the nightly build
			// std::vector<RawDataPoint2D,AllocT> tmp;
			tmp.reserve(input_map.getSize());
			input_map.get2DData(tmp);
			std::partial_sort( tmp.begin(), tmp.begin()+n, tmp.end(), reverseComparator(RawDataPoint2D::IntensityLess()) );
			for ( UInt element_index = 0; element_index < n; ++element_index )
			{
				output_map.push_back( ConsensusFeature( input_map_index, element_index, tmp[element_index] ) );
			}
			return;
		}
		
    /// Set element pair list
    void setElementPairs(ElementPairVectorType& element_pairs)
    {
      element_pairs_ = &element_pairs;
    }

    /// Get element pair list (non-mutable)
    const ElementPairVectorType& getElementPairs() const
    {
      return *element_pairs_;
    }

    /// Set transformation
    void setTransformation(UInt dim, const TransformationType& trafo)
    {
      transformation_[dim] = trafo;
    }

    /// Get transformation
    const TransformationType& getTransformation(UInt dim) const
    {
      return transformation_[dim];
    }

    /// Register all derived classes here
    static void registerChildren();

    /**@brief Write a debug dump of element pairs, e.g. for viewing the result
       with Gnuplot.  The details (e.g. output destination) are controlled by
       param_.  Returns 0 upon success and -1 if no debug dump was requested
       according to param_.  You might invoke this at the end of run() in
       derived classes.

       The base class will treat parameter "debug:dump_element_pairs" as a
       string, namely the filename for the element pair data and append ".gp"
       to that filename for a gnuplot script.

		@todo Adapt this to ConsensusMap. The functionality itself is still quite useful. (Clemens)
     */
    virtual int dumpElementPairs(const String& filename); // code is below

#if 0
    /// Estimates the transformation for each grid cell
    virtual void findElementPairs() = 0;
#endif

  protected:
		
		/** @brief Array of pointers to model and scene map
		
		Normally you will use maps_.model_ etc. to access these.
		The reason why we use an array is because this way algorithms can easily <i>loop</i> over all maps.
		*/
		union
		{
			ConsensusMap * element_map_[3]; ///< @sa Maps_
			struct
			{
				ConsensusMap const * model_; ///< pointer to model map
				ConsensusMap const * scene_; ///< pointer to scene map 
				ConsensusMap * result_; ///< pointer to result map
			} maps_;
		};
		
		/**
		@brief Symbolic names to make usage of element_map_ more understandable and maintainable.
		*/
		// note: RESULT is already #defined in ClassTest.h!
		// note2: You are not allowed to remove or comment this out ;-)
		enum Maps_ { MODEL_ = 0, SCENE_ = 1, RESULT_ = 2 };
		
		/**@brief This tells us the map indices of the model and the scene map or
		whether their consensus features shall be unpacked when they are added to
		the result.
		
		@sa element_map_
		*/
		union
		{
			Int element_map_index_[2]; ///< @sa Maps_
			struct
			{
				Int model_;
				Int scene_;
			} map_index_;
		};

		/// Transformation in rt and mz dimension // TODO remove
    TransformationType transformation_[2];

    /// Vector of pairs of elements that have been identified by the element matcher // TODO remove
    mutable ElementPairVectorType * element_pairs_;

	 private:

    /// Copy constructor intentionally not implemented
    BasePairFinder(const BasePairFinder&);
		
    /// Assignment operator intentionally not implemented
    BasePairFinder & operator=(const BasePairFinder&);
		
	};

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H
