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

#ifndef OPENMS_KERNEL_CONSENSUSFEATURE_H
#define OPENMS_KERNEL_CONSENSUSFEATURE_H

#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/Group.h>

namespace OpenMS
{

  /**
    @brief A D-dimensional feature.
    
    A feature represents a subset of peaks in a map.  In general, it
    summarizes all peaks related to a specific peptide or chemical entity and
    thus reduces partitions of the LCMS dataset to a more meaningful entity.
    Peaks can be converted to features through the FeatureFinder.

    Features are usually contained in ConsensusFeatureMap.  Features themselves can
    either contain features again (composite design pattern) or their
    constituent peaks.

    Hierarchical relationships in features (i.e. features containing features
    containing features...) can be used to express group relationships.  For
    example, one might group the features corresponding to an ICAT pair into
    an aggregate ICAT feature.

    Features are derived from DPeak, as they inherit most of their properties.
    In particular, a feature has a position and an intensity (area).  The
    position of a feature is defined as maximum position of the model for the
    retention time dimension and the mass of the monoisotopic peak for the m/z
    dimension.
    
    @ingroup Kernel, Serialization
  */
  template < typename ContainerT = DFeatureMap< 2, DFeature<2, KernelTraits > > >
  class ConsensusFeature :  public DFeature< 2, KernelTraits >,
        public Group< ContainerT >
  {
    public:
      /**
        @name Type definitions
      */
      //@{
      typedef DFeature<2,KernelTraits> ElementType;
      typedef ContainerT ContainerType;

      typedef typename ElementType::TraitsType TraitsType;
      typedef DPosition < 2, TraitsType > PositionType;
      typedef typename TraitsType::IntensityType IntensityType;
      typedef IndexTuple< ContainerType > IndexTuple;
      typedef Group< ContainerType > Group;
      typedef DRange<2, TraitsType> PositionBoundingBoxType;
      typedef DRange<1, TraitsType> IntensityBoundingBoxType;
      //@}


      /** @name Constructors and Destructor
      */
      //@{
      /// Default constructor
      ConsensusFeature()
          : ElementType(),
          Group(),
          position_range_(),
          intensity_range_()
      {}

      ///
      ConsensusFeature(const PositionType& pos)
          : ElementType(),
          Group(),
          position_range_(),
          intensity_range_()
      {
        this->getPosition() = pos;
      }

      /// Constructor for a singleton consensus feature
      ConsensusFeature(const ContainerType& map, const ElementType& feature)
      {
        try
        {
          this->insert(IndexTuple(map,feature));
        }
        catch(Exception::InvalidValue)
        {}

        this->getPosition() = feature.getPosition();
        this->getIntensity() = feature.getIntensity();

        position_range_.setMinMax(feature.getPosition(),feature.getPosition());
        intensity_range_.setMinMax(feature.getIntensity(),feature.getIntensity());
      }

      /// Constructor
      ConsensusFeature(const ContainerType& map_1, const ElementType& feature_1, const ContainerType& map_2, const ElementType& feature_2)
      {
        try
        {
          this->insert(IndexTuple(map_1,feature_1));
          this->insert(IndexTuple(map_2,feature_2));
        }
        catch(Exception::InvalidValue)
        {}

        computeConsensus_();
      }

      /// Constructor
      ConsensusFeature(const ContainerType& map, const ElementType& feature, const ConsensusFeature& c_feature)
      {
        Group::operator=(c_feature);
        this->insert(IndexTuple(map,feature));

        computeConsensus_();
      }

      /// Constructor
      ConsensusFeature(const ConsensusFeature& c_feature_1, const ConsensusFeature& c_feature_2)
      {
        Group::operator=(c_feature_1);

        for (typename Group::iterator it = c_feature_2.begin(); it != c_feature_2.end(); ++it)
        {
          try
          {
            this->insert(*it);
          }
          catch(Exception::InvalidValue)
          {}
        }

        computeConsensus_();
      }


      /// Copy constructor
      inline ConsensusFeature(const ConsensusFeature& source)
          : ElementType(source),
          Group(source),
          position_range_(source.position_range_),
          intensity_range_(source.intensity_range_)
      {}

      /// Assignement operator
      ConsensusFeature& operator=(const ConsensusFeature& source)
      {
        if (&source==this)
          return *this;

        Group::operator=(source);
        ElementType::operator=(source);
        position_range_=source.position_range_;
        intensity_range_=source.intensity_range_;

        return *this;
      }

      /// Destructor
      virtual ~ConsensusFeature()
    {}
      //@}


      void insert(const IndexTuple& tuple) 
      {
        Group::insert(tuple);

        computeConsensus_();
      }

      /// Non-mutable access to the position range
      inline const PositionBoundingBoxType& getPositionRange() const
      {
        return position_range_;
      }
      /// Mutable access to the position range
      inline PositionBoundingBoxType& getPositionRange()
      {
        return position_range_;
      }
      /// Set the position range
      inline void setPositionRange(const PositionBoundingBoxType& p)
      {
        position_range_= p;
      }

      /// Non-mutable access to the intensity range
      inline const IntensityBoundingBoxType& getIntensityRange() const
      {
        return intensity_range_;
      }
      /// Mutable access to the intensity range
      inline IntensityBoundingBoxType& getIntensityRange()
      {
        return intensity_range_;
      }
      /// Set the intensity range
      inline void setIntensityRange(const IntensityBoundingBoxType& i)
      {
        intensity_range_ = i;
      }

      /// Non-mutable access to the combined features
      inline const Group& getFeatures() const
      {
        return *this;
      }
      /// Mutable access to the combined features
      inline Group& getFeatures()
      {
        return *this;
      }
      /// Set the combined features
      inline void setFeatures(const Group& g)
      {
        Group::operator=(g);
      }


    protected:
      PositionBoundingBoxType position_range_;
      IntensityBoundingBoxType intensity_range_;

      // compute the consensus attributes like intensity and position as well as the position and intensity range given by the group elements
      void computeConsensus_()
      {
        unsigned int n = Group::size();
        DPosition<2> sum_position;
        DPosition<2> pos_min(std::numeric_limits< typename DPeak<2,KernelTraits>::CoordinateType>::max());
        DPosition<2> pos_max(std::numeric_limits< typename DPeak<2,KernelTraits>::CoordinateType>::min());
        DPosition<1> sum_intensities = 0;
        DPosition<1> int_min(std::numeric_limits< typename DPeak<2,KernelTraits>::IntensityType>::max());
        DPosition<1> int_max(std::numeric_limits< typename DPeak<2,KernelTraits>::IntensityType>::min());
        for (typename Group::const_iterator it = Group::begin(); it != Group::end(); ++it)
        {
          DPosition<1> act_int = (it->getElement()).getIntensity();
          DPosition<2> act_pos = it->getTransformedPosition();

          if (int_min > act_int)
          {
            int_min = act_int;
          }
          if (int_max < act_int)
          {
            int_max = act_int;
          }

          for (UnsignedInt dim=0; dim < 2; ++dim)
          {
            if (act_pos[dim] > pos_max[dim])
              pos_max[dim] = act_pos[dim];
            if (act_pos[dim] < pos_min[dim])
              pos_min[dim] = act_pos[dim];
          }

          sum_intensities += act_int;
          sum_position += act_pos;
        }

        for (UnsignedInt dim = 0; dim< 2 ; ++dim)
        {
          this->position_[dim] = sum_position[dim] / n;
        }
        this->intensity_ = sum_intensities[0] / n;

        intensity_range_.setMinMax(int_min,int_max);
        position_range_.setMinMax(pos_min,pos_max);
      }
  };

  ///Print the contents to a stream.
  template < typename ContainerT >
  std::ostream& operator << (std::ostream& os, const ConsensusFeature<ContainerT>& cons)
  {
    os << "---------- CONSENSUS ELEMENT BEGIN -----------------\n"
    << "Position: " << cons.getPosition()<< std::endl
    << "Intensity " << cons.getIntensity() << std::endl
    << "Position range " << cons.getPositionRange() << std::endl
    << "Intensity range " << cons.getIntensityRange() << std::endl
    << "Grouped elements: " << std::endl;

    unsigned int i = 1;
    os << "Size " << cons.count() << std::endl;
    for (typename ConsensusFeature<ContainerT>::Group::const_iterator it = cons.begin(); it != cons.end(); ++it, ++i)
    {
      os  << "Element: " << i << '\n'
      << "Transformed Position: " << it->getTransformedPosition() << '\n'
      << "Original Position: " << it->getElement() << '\n';
    }
    os << "---------- CONSENSUS ELEMENT END ----------------- " << std::endl;

    return os;
  }

} // namespace OpenMS

#endif // OPENMS_KERNEL_DFEATURE_H
