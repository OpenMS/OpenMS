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

#ifndef OPENMS_KERNEL_CONSENSUSPEAK_H
#define OPENMS_KERNEL_CONSENSUSPEAK_H

#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/Group.h>

namespace OpenMS
{

  /**
    @brief A 2-dimensional consensus peak.
    
    A consensus peak represents corresponding peaks in multiple peakmaps.
    
    @ingroup Kernel, Serialization
  */
  template < typename ContainerT = DPeakArray< 2, Peak2D > >
  class ConsensusPeak :  public Peak2D,
        public Group< ContainerT >
  {
    public:
      /**
        @name Type definitions
      */
      //@{
      typedef Peak2D ElementType;
      typedef ContainerT ElementContainerType;
      typedef Group< ElementContainerType > Group;

      typedef typename ElementType::TraitsType TraitsType;
      typedef DPosition < 2, TraitsType > PositionType;
      typedef typename TraitsType::IntensityType IntensityType;
      typedef IndexTuple< ElementContainerType > IndexTuple;
      typedef DRange<2, TraitsType> PositionBoundingBoxType;
      typedef DRange<1, TraitsType> IntensityBoundingBoxType;
      //@}


      /** @name Constructors and Destructor
      */
      //@{
      /// Default constructor
      ConsensusPeak()
          : ElementType(),
          Group(),
          position_range_(),
          intensity_range_()
      {}

      ///
      ConsensusPeak(const PositionType& pos, const IntensityType& i)
          : ElementType(),
          Group(),
          position_range_(),
          intensity_range_()
      {
        this->getPosition() = pos;
        this->getIntensity() = i;
      }

      /// Constructor for a singleton consensus peak
      ConsensusPeak(const UnsignedInt& map_index,  const UnsignedInt& peak_index, const ElementType& peak)
      {
        try
        {
          IndexTuple i(map_index,peak_index,peak);
          i.setTransformedPosition(peak.getPosition());
          this->insert(IndexTuple(map_index,peak_index,peak));
        }
        catch(Exception::InvalidValue)
        {}

        this->getPosition() = peak.getPosition();
        this->getIntensity() = peak.getIntensity();

        position_range_.setMinMax(peak.getPosition()
                                  ,peak.getPosition());
        intensity_range_.setMinMax(peak.getIntensity(),peak.getIntensity());
      }

      /// Constructor
      ConsensusPeak(const UnsignedInt& map_1_index, const UnsignedInt& peak_index_1, const ElementType& peak_1,
                    const UnsignedInt& map_2_index, const UnsignedInt& peak_index_2, const ElementType& peak_2)
      {
        try
        {
          IndexTuple i1(map_1_index,peak_index_1, peak_1);
          i1.setTransformedPosition(peak_1.getPosition());
          this->insert(i1);
          IndexTuple i2(map_2_index,peak_index_2, peak_2);
          i2.setTransformedPosition(peak_2.getPosition());
          this->insert(i2);
        }
        catch(Exception::InvalidValue)
        {}

        computeConsensus_();
      }

      /// Constructor
      ConsensusPeak(const UnsignedInt& map_index, const UnsignedInt& peak_index, const ElementType& peak, const ConsensusPeak& c_peak)
      {
        Group::operator=(c_peak);
        IndexTuple i(map_index,peak_index,peak);
        i.setTransformedPosition(peak.getPosition());
        this->insert(IndexTuple(map_index,peak_index,peak));

        computeConsensus_();
      }

      /// Constructor
      ConsensusPeak(const ConsensusPeak& c_peak_1, const ConsensusPeak& c_peak_2)
      {
        Group::operator=(c_peak_1);

        for (typename Group::iterator it = c_peak_2.begin(); it != c_peak_2.end(); ++it)
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
      inline ConsensusPeak(const ConsensusPeak& source)
          : ElementType(source),
          Group(source),
          position_range_(source.position_range_),
          intensity_range_(source.intensity_range_)
      {}

      /// Assignement operator
      ConsensusPeak& operator=(const ConsensusPeak& source)
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
      virtual ~ConsensusPeak()
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

      /// Non-mutable access to the combined peaks
      inline const Group& getPeaks() const
      {
        return *this;
      }
      /// Mutable access to the combined peaks
      inline Group& getPeaks()
      {
        return *this;
      }
      /// Set the combined peaks
      inline void setPeaks(const Group& g)
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
  std::ostream& operator << (std::ostream& os, const ConsensusPeak<ContainerT>& cons)
  {
    os << "---------- CONSENSUS ELEMENT BEGIN -----------------\n";
    os << "Position: " << cons.getPosition()<< std::endl;
    os << "Intensity " << cons.getIntensity() << std::endl;
    os << "Position range " << cons.getPositionRange() << std::endl;
    os << "Intensity range " << cons.getIntensityRange() << std::endl;
    os << "Grouped elements: " << std::endl;

    unsigned int i = 1;
    os << "Size " << cons.count() << std::endl;
    for (typename ConsensusPeak<ContainerT>::Group::const_iterator it = cons.begin(); it != cons.end(); ++it, ++i)
    {
      os  << "Element: " << i << '\n'
      << "Transformed Position: " << it->getTransformedPosition() << '\n'
      << "Original Position: " << it->getElement() << '\n'
      << "Element index " << it->getElementIndex() << '\n'
      << "Map index " << it->getMapIndex() << std::endl;

    }
    os << "---------- CONSENSUS ELEMENT END ----------------- " << std::endl;

    return os;
  }

} // namespace OpenMS

#endif // OPENMS_KERNEL_Dpeak_H
