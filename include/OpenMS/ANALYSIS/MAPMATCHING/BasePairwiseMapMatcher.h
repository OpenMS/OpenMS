// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                  OpenMS Mass Spectrometry Framework
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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRWISEMAPMATCHER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRWISEMAPMATCHER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/Grid.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <utility>
#include <fstream>

namespace OpenMS
{

  /**
		@brief The base class of all pairwise point matching algorithms.
		
		This class defines the basic interface for all point matching
		algorithms.  
		It works on two point maps and computes a vector of corresponding points
		in both maps (given by a point pairs vector). 
		A point can be a DPeak, a DFeature, a ConsensusPeak or ConsensusFeature 
		(wheras DFeature is the default element type).
		
		The point pairs created by the algorithm solve a
		(bipartite) matching problem between two point maps.
		Therefore first a transformation is estimated, that maps the one map 
		(the so called scene map) onto the other map (the so called model map).
		Given the transformation correspoinding elements in the two maps are determined.
		
		@note If a piecewise transformation is assumed, the user can define a grid by setting 
		the number of buckets in the RT as well as the MZ dimension. 
		Call initGridTransformation() before run()!
		
		@todo Avoid the "0.01 hack" in initGridTransformation(). (Eva)
		
		@ingroup Analysis
  */
  template < typename MapT = FeatureMap< > >
  class BasePairwiseMapMatcher 
  	: public FactoryProduct
  {
  public:
    /// Container for input elements
    typedef MapT PointMapType;

    /// Type of elements considered here
    typedef typename PointMapType::value_type ElementType;

    /// Type of element pairs
    typedef ElementPair < ElementType > ElementPairType;

    /// Container for generated element pairs
    typedef std::vector < ElementPairType > ElementPairVectorType;

    /// Position
    typedef DPosition < 2 > PositionType;

    ///
    typedef DBoundingBox< 2>  PositionBoundingBoxType;

    /// Coordinate
    typedef DoubleReal CoordinateType;


    /// Constructor
    BasePairwiseMapMatcher()
        : FactoryProduct("BasePairWiseMapMatcher")
    {
      element_map_[0] = 0;
      element_map_[1] = 0;
      defaults_.setValue("number_buckets:RT",1);
      defaults_.setValue("number_buckets:MZ",1);
      
			// no need to call defaultsToParam_() as it is called in the non-abstract children 
    }

    /// Copy constructor
    BasePairwiseMapMatcher(const BasePairwiseMapMatcher& source)
        : FactoryProduct(source),
        all_element_pairs_(source.all_element_pairs_),
        bounding_box_scene_map_(source.bounding_box_scene_map_),
        box_size_(source.box_size_)
    {
      element_map_[0] = source.element_map_[0];
      element_map_[1] = source.element_map_[1];
      grid_ = source.grid_;
    	
    	// no need to call defaultsToParam_() as it is called in the non-abstract children 
    }

    ///  Assignment operator
    BasePairwiseMapMatcher& operator = (const BasePairwiseMapMatcher& source)
    {
      if (&source==this) return *this;

      FactoryProduct::operator = (source);
      	
      element_map_[0] = source.element_map_[0];
      element_map_[1] = source.element_map_[1];
      all_element_pairs_ = source.all_element_pairs_;
      grid_ = source.grid_;
      bounding_box_scene_map_ = source.bounding_box_scene_map_;
      box_size_ = source.box_size_;
      
      // no need to call defaultsToParam_() as it is called in the non-abstract children 
      
      return *this;
    }

    /// Destructor
    virtual ~BasePairwiseMapMatcher()
  	{
  	}
		
    /// Set element map
    void setElementMap(UInt const index, const PointMapType& element_map)
    {
      element_map_[index] = &element_map;
    }

    /// Get element map
    const PointMapType& getElementMap(UInt index) const
    {
      return *element_map_[index];
    }

    /// Get element pair list
    const ElementPairVectorType& getElementPairs() const
    {
      return all_element_pairs_;
    }

    /// Get grid
    const Grid& getGrid() const
    {
      return grid_;
    }

    /// Set @p number of buckets in dimension @p dim
    void setNumberBuckets(UInt dim, UInt number)
    {
      number_buckets_[dim] = number;
			param_.setValue(String("number_buckets:") + RawDataPoint2D::shortDimensionName(dim), (Int)number);
    }

    /// Get number of buckets in dimension index
    UInt getNumberBuckets(UInt index) const
    {
      return number_buckets_[index];
    }

    void clearGrid()
    {
      grid_.clear();
    }

    /// Register all derived classes here
    static void registerChildren();

    /// Determine corresponding elements (element pairs)
    virtual void run() = 0;

    /// Initializes the grid for the scene map given the number of buckets in rt and mz. This method has to be called before run()!
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
      // TODO: find a better way that does not use such a magic constant. As is, the bounding box reported in the output is incorrect!
      for (UInt i = 0; i < 2; ++i)
      {
        box_size_[i] =
          (bounding_box_scene_map_.max()[i] - bounding_box_scene_map_.min()[i]) /
          ( number_buckets_[i] - 0.01 /* <- magic constant */);
      }

      // initialize the grid cells of the grid_
      for (UInt x_index = 0; x_index < number_buckets_[RawDataPoint2D::RT]; ++x_index)
      {
        for (UInt y_index = 0; y_index < number_buckets_[RawDataPoint2D::MZ]; ++y_index)
        {
          CoordinateType x_min = (bounding_box_scene_map_.min())[RawDataPoint2D::RT] + box_size_[RawDataPoint2D::RT]*x_index;
          CoordinateType x_max = (bounding_box_scene_map_.min())[RawDataPoint2D::RT] + box_size_[RawDataPoint2D::RT]*(x_index+1);
          CoordinateType y_min = (bounding_box_scene_map_.min())[RawDataPoint2D::MZ] + box_size_[RawDataPoint2D::MZ]*y_index;
          CoordinateType y_max = (bounding_box_scene_map_.min())[RawDataPoint2D::MZ] + box_size_[RawDataPoint2D::MZ]*(y_index+1);

          grid_.push_back(GridCell(x_min, y_min, x_max, y_max));
        }
      }
    } // initGridTransformation_



    //  int dumpElementPairs(const String& filename); // code is below

  protected:
  	virtual void updateMembers_()
  	{
      for ( UInt dim = 0; dim < 2; ++dim)
      {
				number_buckets_[dim] = param_.getValue(String("number_buckets:") + RawDataPoint2D::shortDimensionName(dim));
      }
  	}
  	
    /// Two maps of elements to be matched
    PointMapType const * element_map_[2];

    /// Each element of the vector corresponds to all element pairs of one gridcell
    ElementPairVectorType all_element_pairs_;

    /// The estimated transformation between the two element maps
    Grid grid_;

    /// Bounding box of the second map
    PositionBoundingBoxType bounding_box_scene_map_;

    /// UInt of the grid cells
    PositionType box_size_;

    /// Number of buckets in each dimension
    UInt number_buckets_[2];
  }
  ; // BasePairwiseMapMatcher


  //   template < typename MapT >
  //   int BasePairwiseMapMatcher<MapT>::dumpElementPairs(const String& filename)
  //   {
  //     // V_dumpElementPairs() is used for a few comments about the files being
  //     // written.
  // #define V_dumpElementPairs(bla) std::cerr << bla << std::endl;
  //     V_dumpElementPairs("### Writing "<<filename);
  //     std::ofstream dump_file(filename.c_str());
  //     dump_file << "# " << filename<< " generated " << Date::now() << std::endl;
  //     dump_file << "# 1:number 2:quality 3:firstRT 4:firstMZ 5:firstIT 6:firstQual 7:secondRT 8:secondMZ 9:secondIT 10:secondQual\n";
  //     for ( UInt fp = 0; fp < getElementPairs().size(); ++fp )
  //     {
  //       dump_file << fp << ' '
  //       << getElementPairs()[fp].getQuality() << ' '
  //       << getElementPairs()[fp].getFirst().getRT() << ' '
  //       << getElementPairs()[fp].getFirst().getMZ() << ' '
  //       << getElementPairs()[fp].getFirst().getIntensity() << ' '
  //       << getElementPairs()[fp].getFirst().getOverallQuality() << ' '
  //       << getElementPairs()[fp].getSecond().getRT() << ' '
  //       << getElementPairs()[fp].getSecond().getMZ() << ' '
  //       << getElementPairs()[fp].getSecond().getIntensity() << ' '
  //       << getElementPairs()[fp].getSecond().getOverallQuality() << ' '
  //       << std::endl;
  //     }
  //     dump_file << "# " << filename << " EOF " << Date::now() << std::endl;
  //     std::string dump_filename_gp = filename + ".gp";
  //     V_dumpElementPairs("### Writing "<<dump_filename_gp);
  //     std::ofstream dump_file_gp(dump_filename_gp.c_str());
  //     dump_file_gp << "# " << dump_filename_gp << " generated " << Date::now() << std::endl;
  //     dump_file_gp <<
  //     "# Gnuplot script to view element pairs\n"
  //     "plot   \"" << filename <<"\" using 3:4 title \"map 1\"\n"
  //     "replot \"" << filename <<"\" using 7:8 title \"map 2\"\n"
  //     "replot \"" << filename <<"\" using 3:4:($7-$3):($8-$4) w vectors nohead title \"pairs\"\n"
  //     ;
  //     dump_file_gp << "# " << dump_filename_gp << " EOF " << Date::now() << std::endl;
  //     V_dumpElementPairs("### You can view `"<<filename<<"' using the command line `gnuplot "<<dump_filename_gp<<" -'");
  // #undef V_dumpElementPairs
  //
  //     return 0;
  //   }

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRWISEMAPMATCHER_H
