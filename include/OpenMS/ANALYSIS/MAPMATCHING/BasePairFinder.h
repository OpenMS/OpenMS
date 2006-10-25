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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DLinearMapping.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <utility>
#include <fstream>

namespace OpenMS
{

  /**
     @brief The base class of all feature pair finding algorithms.

     This class defines the basic interface for all feature pair finding 
     algorithms. It works on two feature maps (DFeatureMap< D, Traits, DFeature< D, Traits > >
     is the default map type, but you can also use a pointer map like 
     DPeakConstReferenceArray<2, KernelTraits, DFeatureMap<2> >), a vector of feature
     pairs, and a transformation defined for the second feature map (if no
     transformation is given, the pairs are found in the two original maps).

  **/
  template < typename MapT = DFeatureMap< 2, DFeature< 2, KernelTraits > > >
  class BasePairFinder : public FactoryProduct
  {
    public:
      /// Defines the coordinates of peaks / features.
      enum DimensionId
      {
        RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
        MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

      /** Symbolic names for indices of feature maps etc.
            This should make things more understandable and maintainable.
             */
      enum Maps
      {
        MODEL = 0,
        SCENE = 1
    };


      /** @name Type definitions
       */
      //@{
      /// Container for input features
      typedef MapT PointMapType;

      /// Type of features considered here
      typedef typename PointMapType::value_type PointType;

      /// Traits type
      typedef typename PointType::TraitsType TraitsType;

      /// Position
      typedef DPosition < 2, TraitsType > PositionType;

      /// Quality
      typedef typename TraitsType::QualityType QualityType;

      //// Intensity
      typedef typename TraitsType::IntensityType IntensityType;

      /// Type of feature pairs
      typedef DFeaturePair < 2, PointType > FeaturePairType;

      /// Container for generated feature pairs
      typedef DFeaturePairVector < 2, PointType > FeaturePairVectorType;

      /// Type of estimated transformation
      typedef DBaseMapping< 1, TraitsType > TransformationType;
      //@}


      ///@name Constructors, destructor and assignment
      //@{
      /// Constructor
      BasePairFinder()
          : FactoryProduct(),
          param_(),
          feature_pairs_(0)
      {
        feature_map_[MODEL] = 0;
        feature_map_[SCENE] = 0;
        transformation_[MODEL] = 0;
        transformation_[MZ] = 0;
      }

      /// Copy constructor
      BasePairFinder(const BasePairFinder& source)
          : FactoryProduct(source),
          param_(source.param_),
          feature_pairs_(source.feature_pairs_)
      {
        feature_map_[MODEL] = source.feature_map_[MODEL];
        feature_map_[SCENE] = source.feature_map_[SCENE];
        transformation_[RT] = source.transformation_[RT];
        transformation_[MZ] = source.transformation_[MZ];
      }

      ///  Assignment operator
      virtual BasePairFinder& operator = (BasePairFinder source)
      {
        if (&source==this)
          return *this;

        FactoryProduct::operator = (source);
        feature_map_[MODEL] = source.feature_map_[MODEL];
        feature_map_[SCENE] = source.feature_map_[SCENE];
        transformation_[RT] = source.transformation_[RT];
        transformation_[MZ] = source.transformation_[MZ];
        return *this;
      }

      /// Destructor
      virtual ~BasePairFinder()
    {}
      //@}

      /** @name Accesssor methods
       */
      //@{

      /// Set param class
      void setParam(const Param& param)
      {
        param_ = param;
      }
      /// Get param class
      Param& getParam()
      {
        return param_;
      }
      /// Get param class (non-mutable)
      const Param& getParam() const
      {
        return param_;
      }
      /// Set feature map by arg
      void setFeatureMap(Size const index, const PointMapType& feature_map)
      {
        feature_map_[index] = &feature_map;
      }
      /// Get feature map by arg
      const PointMapType& getFeatureMap(Size index)
      {
        return *feature_map_[index];
      }
      /// Get feature maps by arg (non-mutable)
      const PointMapType& getFeatureMap(Size index) const
      {
        return *feature_map_[index];
      }
      /// Set feature pair list
      void setFeaturePairs(FeaturePairVectorType& feature_pairs)
      {
        feature_pairs_ = &feature_pairs;
      }
      /// Get feature pair list
      FeaturePairVectorType& getFeaturePairs()
      {
        return *feature_pairs_;
      }
      /// Get feature pair list (non-mutable)
      const FeaturePairVectorType getFeaturePairs() const
      {
        return *feature_pairs_;
      }
      /// Set transformation
      void setTransformation(Size const dim, const TransformationType& t)
      {
        transformation_[dim] = &t;
      }
      /// Get transformation
      const TransformationType& getTransformation() const
      {
        return *transformation_;
      }
      //@}

      /// register all derived classes here
      static void registerChildren();

      /**@brief Write a debug dump of feature pairs, e.g. for viewing the result
         with Gnuplot.  The details (e.g. output destination) are controlled by
         param_.  Returns 0 upon success and -1 if no debug dump was requested
         according to param_.  You might invoke this at the end of run() in
         derived classes.

         The base class will treat parameter "debug:dump_feature_pairs" as a
         string, namely the filename for the feature pair data and append ".gp"
         to that filename for a gnuplot script.
       */
      virtual int dumpFeaturePairs(const String& filename); // code is below

      /// Estimates the transformation for each grid cell
      virtual void run()
      {}

    protected:

      /** @name Data members
       */
      //@{
      /// Param class containing the parameters for the map matching phase
      Param param_;

      /// Two maps of features to be matched
      PointMapType const * feature_map_[2];

      /// Vector of pairs of features that have been identified by the feature matcher
      mutable FeaturePairVectorType * feature_pairs_;

      /// If the second map should be transformed before the pair finding process, this transformation should be set
      TransformationType const* transformation_[2];
      //@}

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


  //   template <typename MapT >
  //   int BasePairFinder<MapT>::dumpFeaturePairs(const String& filename)
  //   {
  //     // V_dumpFeaturePairs() is used for a few comments about the files being
  //     // written.  We are silent unless output is actually being written, so
  //     // it is defined here inside the "else" branch.
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

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BasePairFinder_H
