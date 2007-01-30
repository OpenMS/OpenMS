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

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_EXTRACTSIGNALREGIONS_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_EXTRACTSIGNALREGIONS_H

#include <OpenMS/FORMAT/Param.h>

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
  
  @todo use DefaultParamHandler (Eva)
*/

class ExtractSignalRegions
{
public:
    /// Default constructor
    inline ExtractSignalRegions()
            : dalton_per_split_(10)
    {}

    /// Constructor given a param object
    inline ExtractSignalRegions(const Param& parameters)
    {
        param_ = parameters;

        // if a peak picking parameter is missed in the param object the value should be substituted by a default value
        DataValue dv;
        dv = param_.getValue("Split:DaltonPerSplit");
        if (dv.isEmpty() || dv.toString() == "")
            dalton_per_split_ = 10;
        else
            dalton_per_split_ = (float)dv;

#ifdef DEBUG_PEAK_PICKING

        std::cout << "dalton per split " << dalton_per_split_
        << std::endl;
#endif

    }

    /// Copy constructor
    inline ExtractSignalRegions(const ExtractSignalRegions&  e)
            : dalton_per_split_(e.dalton_per_split_)
    {}

    /// Destructor
    virtual ~ExtractSignalRegions()
    {}


    /// Assignment operator
    inline ExtractSignalRegions& operator=(const ExtractSignalRegions& e)
    {
        // take care of self assignments
        if (this == &e)
        {
            return *this;
        }

        dalton_per_split_ = e.dalton_per_split_;
        param_ = e.param_;

        return *this;
    }

    /// Non-mutable access to the decomposition length
    inline const float& getDaltonPerSplit() const
    {
        return dalton_per_split_;
    }
    /// Mutable access to the decomposition length
    inline float& getDaltonPerSplit()
    {
        return dalton_per_split_;
    }
    /// Mutable access to the decomposition length
    inline void setDaltonPerSplit(const float& dalton_per_split)
    {
        dalton_per_split_ = dalton_per_split;
    }

    /// Non-mutable access to the parameter object
    inline const Param& getParam() const
    {
        return param_;
    }
    /// Mutable access to the parameter object
    inline Param& getParam()
    {
        return param_;
    }
    /// Mutable access to the parameter object
    inline void setParam(const Param& param)
    {
        param_ = param;
    }


    /** @brief Splits the array defined by begin/end into independent areas.
     
        This can for example reduces the runtime of the peak picking algorithm.
        As criterium, a signal is splitted into parts if more than points_for_split
        points are less than maxIntensity_*signalToNoiseThreshold_ (means that a quite large area of the signal is
        unimportant)
        The begin and end iterators of each split are written into an iterator vector.
    */
    template < typename InputPeakIterator >
    void splitScan(InputPeakIterator it_begin, InputPeakIterator it_end, double noise_level, std::vector<InputPeakIterator>& splitted_array)
    {
#ifdef DEBUG_PEAK_PICKING
        std::cout << "start splitScan " << dalton_per_split_ << " with noise_level " << noise_level
        << std::endl;
#endif

        splitted_array.clear();

        InputPeakIterator end_point = it_end-1;
        InputPeakIterator first=it_begin,last=it_end-1;

        // Search (from both ends) for the first raw data points which exeed the noise level
        while((first < end_point) && (first->getIntensity() < noise_level))
        {
            first++;
        }

        if (first==last)
            return;


        while((last > first) && (last->getIntensity() < noise_level))
        {
            last--;
        }

        first = ((first-10) < it_begin) ? it_begin : (first-10);
        last = ((last+10) > it_end) ? it_end : (last+10);

        // if the split contains no data points with intensity greater than noise, discard it
        InputPeakIterator greater_noise = it_end;
        InputPeakIterator new_end = first+2;
        while(new_end != (last-1))
        {
#ifdef DEBUG_PEAK_PICKING
            std::cout << "while : first " << first->getPos() << " while : new_end " << new_end->getPos()
            <<  " " << ((new_end+1)->getPos() - new_end->getPos())<< std::endl;
#endif

            if (new_end->getIntensity() > noise_level)
            {
                greater_noise = new_end;
            }

            // is there a gap between new_end and the next data point, which is greater than one dalton?
            if (((new_end+1)->getPos() - new_end->getPos()) > 1)
            {
#ifdef DEBUG_PEAK_PICKING
                std::cout << "gap " << new_end->getPos() << " und first " << first->getPos() << " dist: " << distance(first,new_end)<<std::endl;
#endif
                // and are there enough datapoints between first and new_end?
                if ((distance(first,(new_end+1)) > 3) && (greater_noise != it_end))
                {
                    // cut the signal
                    splitted_array.push_back(first);
                    splitted_array.push_back(new_end+1);
#ifdef DEBUG_PEAK_PICKING

                    std::cout << " find the end x da away: " <<  std::endl
                    << " push first " << first->getPos()
                    << " push last " << (new_end)->getPos()
                    << " distance " << distance(first,(new_end+1))
                    << std::endl;
#endif

                    greater_noise = it_end;

                }
                //else throw away this split
                first = new_end+1;
                ++new_end;
#ifdef DEBUG_PEAK_PICKING

                std::cout << "after gap: first " << first->getPos() << std::endl;
#endif

                continue;
            }

            // if the split has already  a length of dalton_per_split_ search for a minimum and cut
            if (((new_end)->getPos() - first->getPos()) > dalton_per_split_)
            {
#ifdef DEBUG_PEAK_PICKING
                std::cout << "dalton_per_split_ reached " << first->getPos()  << " " << (new_end)->getPos() << std::endl;
#endif

                int  search_radius = 5;

                // if the split contains no data point with intensity greater than noise, discard the split
                if (greater_noise == it_end)
                {
                    first = new_end;
                    std::cout << "No greater_noise " << std::endl;
                }
                else
                {
                    // Search to the left for a minimum
                    for (int i=0; i < search_radius; i++)
                    {
                        InputPeakIterator left_temp, right_temp;
                        left_temp = new_end - i;

                        if (((average_(left_temp-1,it_begin, 2) < average_(left_temp-2,it_begin, 2))
                                &&(average_(left_temp-1,it_begin,2) < average_(left_temp,it_begin,2))
                                && left_temp->getIntensity() < noise_level) && (greater_noise <= left_temp))
                        {
                            // if the signal is still falling to the left go on until a minimum is found
                            while ((left_temp-2)->getIntensity() < (left_temp-1)->getIntensity())
                            {
                                --left_temp;
                            }

                            splitted_array.push_back(first);
                            splitted_array.push_back(left_temp);
#ifdef DEBUG_PEAK_PICKING

                            std::cout << "minimum found " << std::endl
                            << " push " << first->getPos()
                            << " and " << (left_temp-1)->getPos() << std::endl;
#endif

                            new_end = left_temp;
                            first = new_end;
                            greater_noise = it_end;

                            break;
                        }
                    }
                }
            }
            ++new_end;
        }

        // does the last split contain enough data points and does the split contain any data value with an intensity greater then noise?
        if ((distance(first, new_end) > 3) && (greater_noise != it_end))
        {
#ifdef DEBUG_PEAK_PICKING
            std::cout << " push the last split: first " << first->getPos() << std::endl
            << "push last " << (new_end)->getPos()
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
    /// Length of the regions the signal is decomposed to
    float dalton_per_split_;
    ///
    /// Parameter object
    Param param_;

    /// compute the average of the act_pos intensity and the intensities of its left number data points
    template < typename InputPeakIterator >
    inline double average_(InputPeakIterator act_pos, InputPeakIterator first, int number)
    {
        int i = 0;
        int k = 1;
        double mean = act_pos->getIntensity();
        for (i = 1; ((i < number) && (act_pos - i) >= first) ; ++i)
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
