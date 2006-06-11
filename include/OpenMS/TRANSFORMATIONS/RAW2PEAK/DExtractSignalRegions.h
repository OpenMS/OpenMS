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
// $Id: DExtractSignalRegions.h,v 1.31 2006/05/29 15:52:15 elange Exp $
// $Author: elange $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_DEXTRACTSIGNALREGIONS_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_DEXTRACTSIGNALREGIONS_H

#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/KERNEL/DPeakArrayNonPolymorphic.h>

#include <OpenMS/FORMAT/Param.h>

#include <OpenMS/KERNEL/DimensionDescription.h>

#include <iostream>
#include <vector>
#include <math.h>

//#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING

namespace OpenMS
{
  /**
    @brief Class to decompose raw data into smaller boxes.

    To improve the run time of the peak picker the raw mass spectra are decomposed into
    smaller parts.

    @todo fix and add test
  */

  template <Size D, typename ContainerType = DPeakArrayNonPolymorphic<D,DRawDataPoint<D> > >
  class DExtractSignalRegions
  {
  public:
    /** @name Type definitions
     */
    //@{
    typedef OpenMS::DimensionDescription < DimensionDescriptionTagLCMS > DimensionDescription;
    ///
    typedef typename ContainerType::PeakType PeakType;
    ///
    typedef typename ContainerType::const_iterator ConstIterator;
    ///
    typedef std::vector<ConstIterator> IteratorVector;
    ///
    typedef typename IteratorVector::iterator IteratorIterator;
    ///
    typedef typename IteratorVector::const_iterator IteratorConstIterator;
    //@}


    /** @name Constructors and Destructor
     */
    //@{
    /// Default constructor
    inline DExtractSignalRegions()
        : dalton_per_split_(10)
    {
      if (D == 1)
      {
        rt_dim_ = -1;
        mz_dim_ = 0;
      }
      else
        if (D == 2)
        {
          rt_dim_ = DimensionDescription::RT;
          mz_dim_ = DimensionDescription::MZ;
        }
    }
    ///
    inline DExtractSignalRegions(const Param& parameters)
    {
      param_ = parameters;
      if (D == 1)
      {
        rt_dim_ = -1;
        mz_dim_ = 0;
      }
      else
        if (D == 2)
        {
          rt_dim_ = DimensionDescription::RT;
          mz_dim_ = DimensionDescription::MZ;
        }

      // if a peak picking parameter is missed in the param object the value should be substituted by a default value
      DataValue dv;
      dv = param_.getValue("Split:DaltonPerSplit");
      if (dv.isEmpty() || dv.toString() == "") dalton_per_split_ = 10;
      else dalton_per_split_ = (float)dv;

#ifdef DEBUG_PEAK_PICKING
      std::cout << "dalton per split " << dalton_per_split_
      << std::endl;
#endif

    }

    /// Copy constructor
  inline DExtractSignalRegions(const DExtractSignalRegions&  e) : mz_dim_(e.mz_dim_), rt_dim_(e.rt_dim_), dalton_per_split_(e.dalton_per_split_) {}
    /// Destructor
    virtual ~DExtractSignalRegions()
    {}
    //@}

    /** @name Assignment
     */
    //@{
    inline DExtractSignalRegions& operator=(const DExtractSignalRegions& e)
    {
      mz_dim_=e.mz_dim_;
      rt_dim_=e.rt_dim_;
      dalton_per_split_ = e.dalton_per_split_;
      param_ = e.param_;

      return *this;
    }
    //@}

    /** Accessors
     */
    //@{
    /// Non-mutable access to he mz dimension
    inline const int& getMZdim() const { return mz_dim_; }
    /// Mutable access to the mz dimensin
    inline int& getMZdim() { return mz_dim_; }
    /// Mutable access to the mz dimensin
    inline void setMZdim(const int& mz_dim) { mz_dim_ = mz_dim; }

    /// Non-mutable access to the rt dimension
    inline const int& getRTdim() const { return rt_dim_; }
    /// Mutable access to the rt dimensin
    inline int& getRTdim() { return rt_dim_; }
    /// Mutable access to the rt dimensin
    inline void setRTdim(const int& rt_dim) { rt_dim_ = rt_dim; }

    /// Non-mutable access to the decomposition length
    inline const float& getDaltonPerSplit() const { return dalton_per_split_; }
    /// Mutable access to the decomposition length
    inline float& getDaltonPerSplit() { return dalton_per_split_; }
    /// Mutable access to the decomposition length
    inline void setDaltonPerSplit(const float& dalton_per_split) { dalton_per_split_ = dalton_per_split; }

    /// Non-mutable access to the parameter object
    inline const Param& getParam() const { return param_; }
    /// Mutable access to the parameter object
    inline Param& getParam() { return param_; }
    /// Mutable access to the parameter object
    inline void setParam(const Param& param) { param_ = param; }
    //@}

    /// The filtering method.
    /** Splits the array defined by begin/end into independent areas where peaks are picked independently.
        This reduces the runtime of finding peaks.
        As criterium, a signal is splitted into parts if more than points_for_split
        points are less than maxIntensity_*signalToNoiseThreshold_ (means that a quite large area of the signal is
        unimportant)
        The output is written to splitted_array_ (pairs of beginning/end iterators)
    */
    void splitScan(ConstIterator it_begin, ConstIterator it_end, double noise_level, IteratorVector &splitted_array)
    {
#ifdef DEBUG_PEAK_PICKING
      std::cout << "start splitScan " << dalton_per_split_ << " with noise_level " << noise_level
      << std::endl;
#endif
      splitted_array.clear();

      ConstIterator end_point = it_end-1;
      ConstIterator first=it_begin,last=it_end-1;

      /** First find the first point of interest. **/
      while((first < end_point) && (first->getIntensity() < noise_level))
      {
        first++;
      }

      /** There are no datapoints greater than the noise level **/
      if (first==last)
        return;

      /** Then find the last point of interest. **/
      while((last > first) && (last->getIntensity() < noise_level))
      {
        last--;
      }

      first = ((first-3) < it_begin) ? it_begin : (first-3);
      last = ((last+3) > it_end) ? it_end : (last+3);

      ConstIterator new_end = first+2;
      while(new_end != (last-1))
      {
#ifdef DEBUG_PEAK_PICKING
        std::cout << "while : first " << first->getPosition()[mz_dim_] << " while : new_end " << new_end->getPosition()[mz_dim_]
        <<  " " << ((new_end+1)->getPosition()[mz_dim_] - new_end->getPosition()[mz_dim_])<< std::endl;
#endif
        /**is there a gap between new_end and the next data point, which is greater than one dalton? **/
        if (((new_end+1)->getPosition()[mz_dim_] - new_end->getPosition()[mz_dim_]) > 1)
        {
#ifdef DEBUG_PEAK_PICKING
          std::cout << "gap " << new_end->getPosition()[mz_dim_] << " und first " << first->getPosition()[mz_dim_] << " dist: " << distance(first,new_end)<<std::endl;
#endif
          // are there enough datapoints between first and new_end?
          if (distance(first,(new_end+1)) > 3)
          {
            // cut the signal
            splitted_array.push_back(first);
            splitted_array.push_back(new_end+1);
#ifdef DEBUG_PEAK_PICKING
            std::cout << " find the end x da away: " <<  std::endl
											<< " push first " << first->getPosition()[mz_dim_] 
											<< " push last " << (new_end)->getPosition()[mz_dim_]
											<< " distance " << distance(first,(new_end+1))
											<< std::endl;
#endif
          }
          //else throw away this split
          first = new_end+1;
          ++new_end;
#ifdef DEBUG_PEAK_PICKING
          std::cout << "after gap: first " << first->getPosition()[mz_dim_] << std::endl;
#endif
          continue;
        }

        // compute the moving average
        if((distance(first,new_end) >= 5) && average_(new_end,it_begin,5) < noise_level)
        {
          // are there enough datapoints between first and new_end?
          if (distance(first,new_end-1) > 3)
          {
            // if the signal is homogeneous cut
            splitted_array.push_back(first);
            splitted_array.push_back(new_end-1);
#ifdef DEBUG_PEAK_PICKING
            std::cout << " average find the end x da away: " << std::endl
											<< " push first " << first->getPosition()[mz_dim_] 
											<< " last " << (new_end-2)->getPosition()[mz_dim_]
											<< " distance " << distance(first,(new_end+1))
											<< std::endl;
#endif

          }
          first = new_end-1;
#ifdef DEBUG_PEAK_PICKING
          std::cout << "after : first " << first->getPosition()[mz_dim_] << std::endl;
#endif
          ++new_end;
          continue;
        }


        /** if the split has already  a length of dalton_per_split_ cut search for a minimum **/
        if (((new_end)->getPosition()[mz_dim_] - first->getPosition()[mz_dim_]) > dalton_per_split_)
        {
#ifdef DEBUG_PEAK_PICKING
          std::cout << "dalton_per_split_ reached " << first->getPosition()[mz_dim_]  << " " << (new_end)->getPosition()[mz_dim_] << std::endl;
#endif
          bool found = false;
          int  search_radius = 5;

          /** Search to the left for a minimum **/
          for (int i=0; ((i < search_radius) && (found == false)); i++)
          {
            ConstIterator left_temp, right_temp;
            left_temp = new_end - i;

            if ((average_(left_temp-1,it_begin, 2) < average_(left_temp-2,it_begin, 2))
                &&(average_(left_temp-1,it_begin,2) < average_(left_temp,it_begin,2))
                && left_temp->getIntensity() < noise_level)
            {
              splitted_array.push_back(first);
              splitted_array.push_back(left_temp);
#ifdef DEBUG_PEAK_PICKING
              std::cout << "minimum found " << std::endl
												<< " push " << first->getPosition()[mz_dim_] 
												<< " and " << (left_temp-1)->getPosition()[mz_dim_] << std::endl;
#endif
              new_end = left_temp;
              first = new_end;
              found = true;
              break;
            }
          }
        }
        /** If there is no minimum to be found inside one search_radius,
        *  go on
        */
        ++new_end;
      }

      if (distance(first, new_end) > 3)
      {
#ifdef DEBUG_PEAK_PICKING
        std::cout << " push the last split: first " << first->getPosition()[mz_dim_] << std::endl
        << "push last " << (new_end+1)->getPosition()[mz_dim_]
        << "distance " << distance(first,(new_end+2)) << std::endl;
#endif
        splitted_array.push_back(first);
        splitted_array.push_back(new_end+1);
      }

#ifdef DEBUG_PEAK_PICKING
      std::cout << "end splitScan " << dalton_per_split_
      << std::endl;
#endif

    }

  protected:
    ///
    // MZ dimension
    int mz_dim_;
    ///
    // RT dimension
    int rt_dim_;
    ///
    // length of the boxes the signal is decomposed to
    float dalton_per_split_;
    ///
    // Parameter object
    Param param_;
    //
    inline double average_(ConstIterator act_pos, ConstIterator first, int number)
    {
      int i = 0;
      int k = 1;
      double mean = act_pos->getIntensity();
      for (i = 1; i < number; ++i)
      {
        if ((act_pos-i) >= first)
        {
          mean+=(act_pos-i)->getIntensity();
          ++k;
        }
      }
      return (mean / k) ;
    }
  };

}// namespace OpenMS

#endif
