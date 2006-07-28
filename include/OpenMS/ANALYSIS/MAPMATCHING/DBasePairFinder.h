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
// $Id: DBasePairFinder.h,v 1.6 2006/05/02 09:39:07 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_DBASEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_DBASEPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DBaseMapping.h>

#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DimensionDescription.h>

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
  template <Size D, typename Traits = KernelTraits, typename MapT = DFeatureMap< D, DFeature< D, Traits > > >
  class DBasePairFinder
  {
  public:
    /// Defines the coordinates of peaks / features.
    enum DimensionId
    {
      RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
      MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };
    enum { DIMENSION = D };

    /** @name Type definitions
     */
    //@{
    /// Traits
    typedef Traits TraitsType;

    /// Quality
    typedef typename TraitsType::QualityType QualityType;

    /// Position
    typedef DPosition < DIMENSION, TraitsType > PositionType;

    //// Intensity
    typedef typename TraitsType::IntensityType IntensityType;

    /// Container for input features
    typedef MapT FeatureMapType;
    /// Type of features considered here
    typedef typename FeatureMapType::value_type FeatureType;

    /// Type of feature pairs
    typedef DFeaturePair < DIMENSION, FeatureType > FeaturePairType;
    /// Container for generated feature pairs
    typedef DFeaturePairVector < DIMENSION, FeatureType > FeaturePairVectorType;

    /// Type of estimated transformation
    typedef DBaseMapping< 1, TraitsType > TransformationType;
    //@}


    ///@name Constructors, destructor and assignment
    //@{
    /// Constructor
    DBasePairFinder()
        : param_(),
        feature_pairs_(0)
    {
      feature_map_[0] = 0;
      feature_map_[1] = 0;
      transformation_[0] = 0;
      transformation_[1] = 0;
    }

    /// Copy constructor
    DBasePairFinder(const DBasePairFinder& source)
        : param_(source.param_),
        feature_pairs_(source.feature_pairs_)
    {
      feature_map_[0] = source.feature_map_[0];
      feature_map_[1] = source.feature_map_[1];
      transformation_[0] = source.transformation_[0];
      transformation_[1] = source.transformation_[1];
    }

    ///  Assignment operator
    virtual DBasePairFinder& operator = (DBasePairFinder source)
    {
      param_ = source.param_;
      feature_map_[0] = source.feature_map_[0];
      feature_map_[1] = source.feature_map_[1];
      feature_pairs_ = source.feature_pairs_;
      transformation_[0] = source.transformation_[0];
      transformation_[1] = source.transformation_[1];
      return *this;
    }

    /// Destructor
    virtual ~DBasePairFinder() {}
    //@}

    /** @name Accesssor methods
     */
    //@{

    /// Set param class
    void setParam(const Param& param)
    {
      param_ = param;
      // std::cout << param_ << std::endl; // debugging
    }
    /// Get param class
    Param& getParam() { return param_; }
    /// Get param class (non-mutable)
    const Param& getParam() const { return param_; }

    /// Set feature map by arg
    void setFeatureMap(Size const index, const FeatureMapType& feature_map) { feature_map_[index] = &feature_map; }
    /// Get feature map by arg
    const FeatureMapType& getFeatureMap(Size index) { return *feature_map_[index]; }
    /// Get feature maps by arg (non-mutable)
    const FeatureMapType& getFeatureMap(Size index) const { return *feature_map_[index]; }

    /// Set feature pair list
    void setFeaturePairs(FeaturePairVectorType& feature_pairs) { feature_pairs_ = &feature_pairs; }
    /// Get feature pair list
    FeaturePairVectorType& getFeaturePairs() { return *feature_pairs_; }
    /// Get feature pair list (non-mutable)
    const FeaturePairVectorType getFeaturePairs() const { return *feature_pairs_; }

    /// Set transformation
    void setTransformation(Size const dim, const TransformationType& t)
    {
      transformation_[dim] = &t;
    }
    /// Get transformation
    const TransformationType& getTransformation() const { return *transformation_; }


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

    //@}

    /// Estimates the transformation for each grid cell
    virtual void run() {}; // = 0;



  protected:

    /** @name Data members
     */
    //@{
    /// Param class containing the parameters for the map matching phase
    Param param_;

    /// Two maps of features to be matched
    FeatureMapType const * feature_map_[2];

    /// Vector of pairs of features that have been identified by the feature matcher
    mutable FeaturePairVectorType * feature_pairs_;

    /// If the second map should be transformed before the pair finding process, this transformation should be set
    TransformationType const* transformation_[2];
    //@}

  }
  ; // DBasePairFinder


  template <Size D, typename Traits, typename MapT >
  int DBasePairFinder<D,Traits,MapT>::dumpFeaturePairs(const String& filename)
  {
    // V_dumpFeaturePairs() is used for a few comments about the files being
    // written.  We are silent unless output is actually being written, so
    // it is defined here inside the "else" branch.
#define V_dumpFeaturePairs(bla) std::cerr << bla << std::endl;
    V_dumpFeaturePairs("### Writing "<<filename);
    std::ofstream dump_file(filename.c_str());
    dump_file << "# " << filename<< " generated " << Date::now() << std::endl;
    dump_file << "# 1:number 2:quality 3:firstRT 4:firstMZ 5:firstIT 6:firstQual 7:secondRT 8:secondMZ 9:secondIT 10:secondQual\n";
    for ( Size fp = 0; fp < getFeaturePairs().size(); ++fp )
    {
      dump_file << fp << ' '
      << getFeaturePairs()[fp].getQuality() << ' '
      << getFeaturePairs()[fp].getFirst().getPosition()[RT] << ' '
      << getFeaturePairs()[fp].getFirst().getPosition()[MZ] << ' '
      << getFeaturePairs()[fp].getFirst().getIntensity() << ' '
      << getFeaturePairs()[fp].getFirst().getOverallQuality() << ' '
      << getFeaturePairs()[fp].getSecond().getPosition()[RT] << ' '
      << getFeaturePairs()[fp].getSecond().getPosition()[MZ] << ' '
      << getFeaturePairs()[fp].getSecond().getIntensity() << ' '
      << getFeaturePairs()[fp].getSecond().getOverallQuality() << ' '
      << std::endl;
    }
    dump_file << "# " << filename << " EOF " << Date::now() << std::endl;
    std::string dump_filename_gp = filename + ".gp";
    V_dumpFeaturePairs("### Writing "<<dump_filename_gp);
    std::ofstream dump_file_gp(dump_filename_gp.c_str());
    dump_file_gp << "# " << dump_filename_gp << " generated " << Date::now() << std::endl;
    dump_file_gp <<
    "# Gnuplot script to view feature pairs\n"
    "plot   \"" << filename <<"\" using 3:4 title \"map 1\"\n"
    "replot \"" << filename <<"\" using 7:8 title \"map 2\"\n"
    "replot \"" << filename <<"\" using 3:4:($7-$3):($8-$4) w vectors nohead title \"pairs\"\n"
    ;
    dump_file_gp << "# " << dump_filename_gp << " EOF " << Date::now() << std::endl;
    V_dumpFeaturePairs("### You can view `"<<filename<<"' using the command line `gnuplot "<<dump_filename_gp<<" -'");
#undef V_dumpFeaturePairs

    return 0;
  }


} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_DBasePairFinder_H
