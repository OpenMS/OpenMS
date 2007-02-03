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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DLinearMapping.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <utility>
#include <fstream>

namespace OpenMS
{

  /**
     @brief The base class of all element pair finding algorithms.

     This class defines the basic interface for all element pair finding 
     algorithms. It works on two element maps (DFeatureMap is the default map type, 
     but you can also use a pointer map like DPeakConstReferenceArray),
     and a transformation defined for the second element map (if no
     transformation is given, the pairs are found in the two original maps).
     A element can be a DPeak, a DFeature, a ConsensusPeak or ConsensusFeature 
     (wheras DFeature is the default element type).
          
     Policy for copy constructor and assignment: element_map_ is 
     maintained as pointer and taken shallow copy. 
     But param_ is deep.

  */
  template < typename MapT = DFeatureMap< 2, DFeature< 2, KernelTraits > > >
  class BasePairFinder 
  	: public FactoryProduct
  {
  public:
    /// Defines the coordinates of elements.
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

    /// Container for input elements
    typedef MapT PointMapType;

    /// Type of elements considered here
    typedef typename PointMapType::value_type PointType;

    /// Traits type
    typedef typename PointType::TraitsType TraitsType;

    /// Position
    typedef DPosition < 2, TraitsType > PositionType;

    /// Quality
    typedef typename TraitsType::QualityType QualityType;

    //// Intensity
    typedef typename TraitsType::IntensityType IntensityType;

    /// Type of element pairs
    typedef DFeaturePair < 2, PointType > ElementPairType;

    /// Container for generated element pairs
    typedef DFeaturePairVector < 2, PointType > ElementPairVectorType;

    /// Type of estimated transformation
    typedef DLinearMapping< 1, TraitsType > TransformationType;

    /// Constructor
    BasePairFinder()
        : FactoryProduct("BasePairFinder"),
        element_pairs_(0)
    {
      element_map_[MODEL] = 0;
      element_map_[SCENE] = 0;
      transformation_[RT].setSlope(1);
      transformation_[RT].setIntercept(0);
      transformation_[MZ].setSlope(1);
      transformation_[MZ].setIntercept(0);
    }

    /// Copy constructor
    BasePairFinder(const BasePairFinder& source)
        : FactoryProduct(source),
        element_pairs_(source.element_pairs_)
    {
      element_map_[MODEL] = source.element_map_[MODEL];
      element_map_[SCENE] = source.element_map_[SCENE];
      transformation_[RT] = source.transformation_[RT];
      transformation_[MZ] = source.transformation_[MZ];
    }

    ///  Assignment operator
    BasePairFinder& operator = (const BasePairFinder& source)
    {
      if (&source==this) return *this;

      FactoryProduct::operator = (source);
      
      element_map_[MODEL] = source.element_map_[MODEL];
      element_map_[SCENE] = source.element_map_[SCENE];
      element_pairs_ = source.element_pairs_;
      transformation_[RT] = source.transformation_[RT];
      transformation_[MZ] = source.transformation_[MZ];
      
      return *this;
    }

    /// Destructor
    virtual ~BasePairFinder()
		{
		}

    void setElementMap(Size const index, const PointMapType& element_map)
    {
      element_map_[index] = &element_map;
    }

    /// Get element maps by arg (non-mutable)
    const PointMapType& getElementMap(Size index) const
    {
      return *element_map_[index];
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
    void setTransformation(Size dim, const TransformationType& trafo)
    {
      transformation_[dim] = trafo;
    }

    /// Get transformation
    const TransformationType& getTransformation(Size dim) const
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
     */
    virtual int dumpElementPairs(const String& filename); // code is below

    /// Estimates the transformation for each grid cell
    virtual void findElementPairs() = 0;

  protected:
    /// Two maps of elements to be matched
    PointMapType const * element_map_[2];

    /// Transformation in rt and mz dimension
    TransformationType transformation_[2];

    /// Vector of pairs of elements that have been identified by the element matcher
    mutable ElementPairVectorType * element_pairs_;


    /// Given a position element positoon this method computes the grid cell that covers this point.
    SignedInt computeGridCellIndex_(const PositionType& pos, const DGrid<2>& grid) throw (Exception::InvalidValue)
    {
      UnsignedInt index = 0;
      typename DGrid<2>::ConstIterator it = grid.begin();
      while ( it != grid.end() )
      {
        if (it->encloses(pos))
        {
          return index;
        }
        ++it;
        ++index;
      }
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"The position is not contained in any of the grid cells.","") ;
    }
  }
  ; // BasePairFinder


  template <typename MapT >
  int BasePairFinder<MapT>::dumpElementPairs(const String& filename)
  {
    // V_dumpElementPairs() is used for a few comments about the files being
    // written.  We are silent unless output is actually being written, so
    // it is defined here inside the "else" branch.
#define V_dumpElementPairs(bla) std::cerr << bla << std::endl;
    V_dumpElementPairs("### Writing "<<filename);
    std::ofstream dump_file(filename.c_str());
    dump_file << "# " << filename<< " generated " << Date::now() << std::endl;
    dump_file << "# 1:number 2:quality 3:firstRT 4:firstMZ 5:firstIT 6:firstQual 7:secondRT 8:secondMZ 9:secondIT 10:secondQual\n";
    for ( Size fp = 0; fp < getElementPairs().size(); ++fp )
    {
      dump_file << fp << ' '
      << getElementPairs()[fp].getFirst().getPosition()[RT] << ' '
      << getElementPairs()[fp].getFirst().getPosition()[MZ] << ' '
      << getElementPairs()[fp].getFirst().getIntensity() << ' '
      << getElementPairs()[fp].getSecond().getPosition()[RT] << ' '
      << getElementPairs()[fp].getSecond().getPosition()[MZ] << ' '
      << getElementPairs()[fp].getSecond().getIntensity() << ' '
      << std::endl;
    }
    dump_file << "# " << filename << " EOF " << Date::now() << std::endl;
    std::string dump_filename_gp = filename + ".gp";
    V_dumpElementPairs("### Writing "<<dump_filename_gp);
    std::ofstream dump_file_gp(dump_filename_gp.c_str());
    dump_file_gp << "# " << dump_filename_gp << " generated " << Date::now() << std::endl;
    dump_file_gp <<
    "# Gnuplot script to view element pairs\n"
    "plot   \"" << filename <<"\" using 2:3 title \"map 1\"\n"
    "replot \"" << filename <<"\" using 5:6 title \"map 2\"\n"
    "replot \"" << filename <<"\" using 2:3:($5-$2):($6-$3) w vectors nohead title \"pairs\"\n"
    ;
    dump_file_gp << "# " << dump_filename_gp << " EOF " << Date::now() << std::endl;
    V_dumpElementPairs("### You can view `"<<filename<<"' using the command line `gnuplot "<<dump_filename_gp<<" -'");
#undef V_dumpElementPairs

    return 0;
  }
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H
