// -*- mode: C++; tab-width: 2; -*-
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

#ifndef OPENMS_FILTERING_BASELINE_MORPHFILTER_H
#define OPENMS_FILTERING_BASELINE_MORPHFILTER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <iostream>

namespace OpenMS
{
  /**  
  	@brief This class is the baseclass of morphological filters.
  
		The basic idea of a morphological filter is to inhibit selected signal structures.
		Such structures could be noise or some irrelevant signal structures like the baseline.
		A morphological filter is an increasing operator and has the feature of idempotenz,
		which means that structures which should be received will not be modified by more
		applications of the filter.
		This class provides the basic morphological operations Erosion and Dilatation with a
		structuring element (a flat line) of length frameSize_.
		
		Erosion and dilatation are implemented using van Herk's method.
		 
		@htmlinclude OpenMS_MorphFilter.parameters
*/

  class OPENMS_DLLAPI MorphFilter 
  	: public DefaultParamHandler, 
  		public ProgressLogger
  {
	  public:
	    /// Constructor
	    inline MorphFilter()
	      : DefaultParamHandler("MorphFilter"),
	      	struc_size_(3)
	    {
	      defaults_.setValue("struc_elem_length",10.0,"Length of the structuring element. Should be wider than the expected peak width.");
	    
	      defaultsToParam_();
	    }
	
	    /// Destructor
	    virtual ~MorphFilter()
	    {
	    }
	    
	    /** Van Herk's method of the Dilatation. The algorithm requires 3 min/max comparisons for every data point
	        independent from the length of the structuring element.
	        Basic idea of the dilatation is "Does the structuring element touch a given set?". The value of a data point
	        \f$ x \f$ in the signal \f$ s \f$ after a dilatation is the maximal data point in a window which is represented by
	        the structuring element \f$ B \f$, when the \f$ B \f$'s point of reference is at \f$ x \f$:
	        \f[ [\delta_B(s)](x)=max_{b \in B} s(x+b). \f]
	        \image html Dilation.png "Dilation with a structuring element of length 3"
	    */
	    template < typename InputPeakIterator, typename OutputPeakContainer >
	    void dilatation(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& result, int l)
	    {
	      //--------------van Herk's method of the dilatation --------------------
	      result.resize(distance(first,last));
	
	      int middle = l/2;
	      int i,k,m,n;
	      int length = distance(first,last);
	
	      DoubleReal *g = new DoubleReal[l];
	      DoubleReal *h = new DoubleReal[l];
	      k=length-(length%l)-1;
	
	      calcGDilatation_<InputPeakIterator>(first,last,l,g,true);
	      calcHDilatation_<InputPeakIterator>(first,first+l-1,l,h,true);
	
	      typename OutputPeakContainer::iterator it = result.begin();
	      for (i = 0; i < middle; ++i)
	      {
	        it->setIntensity(g[i+middle]);
	        it->setPosition(first->getPosition());
	        ++first;
	        ++it;
	      }
	
	      m = l-1;
	      n = 0;
	      for (i = middle; i<(length-middle); ++i)
	      {
	        if ((i%l)==(middle+1))
	        {
	          if (i==k)
	          {
					    calcGDilatation_<InputPeakIterator>((first+middle),last,l,g,false);
	          }
	          else
	          {
							calcGDilatation_<InputPeakIterator>((first+middle),last,l,g,true);
	          }
	          m=0;
	        }
	        if ((i%l)==middle && (i>middle))
	        {
	          if (i>k)
	          {
							calcHDilatation_<InputPeakIterator>(first,last,l,h,false); 
	          }
	          else
	          {
							calcHDilatation_<InputPeakIterator>((first-middle),(first+middle),l,h,true);
	          }
	          n=0;
	        }
	
	        it->setIntensity(std::max(g[m],h[n]));
	        it->getPosition() = first->getPosition();	
	        ++it;
	        ++first;
	        ++m;
	        ++n;
	      }
	      
				DoubleReal last_int = (first-1)->getIntensity();
				for (i=0; i<middle; ++i)
	      {
	        it->setIntensity(last_int);
	        it->getPosition() = first->getPosition();
	        ++it;
	        ++first;
	      }
	
	      delete [] g;
	      delete [] h;
	    }
	
	
	
	    /**Van Herk's method of the Erosion. The algorithm requires 3 min/max comparisons for every data point
	       independent from the length of the structuring element.
	       Basic idea of the erosion is "Does the structuring element fit completely in a given set?". The value of a data point
	       \f$ x \f$ in the signal \f$s \f$ after an erosion is the minimal data point in a window which is represented by the
	       structuring element \f$ B\f$, when the \f$ B\f$'s point of reference is at \f$ x \f$:
	       \f[ [\epsilon_B(s)](x)=min_{b \in B} s(x+b). \f]
	       \image html Erosion.png "Erosion with a structuring element of length 3"
	    */
	    template < typename InputPeakIterator, typename OutputPeakContainer >
	    void erosion(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& result, int l)
	    {
	      //-------------- van Herk's method of the erosion --------------------
	      result.resize(distance(first,last));
	
	      int middle=l/2;
	      int i,k,m,n;
	      int length=distance(first,last);
	
	      DoubleReal *g = new DoubleReal[l];
	      DoubleReal *h = new DoubleReal[l];
	      k=length-(length%l)-1;
	
	      calcGErosion_<InputPeakIterator>(first,last,l,g,true);
	      calcHErosion_<InputPeakIterator>(first+l-1,l,h,true);
	
	      typename OutputPeakContainer::iterator it = result.begin();
	      for (i=0; i<middle; ++i)
	      {
	        it->setIntensity(0);
	        it->getPosition() = first->getPosition();
	        ++it;
	        ++first;
	      }
	
	      m = l-1;
	      n = 0;
	      for (i=middle; i<(length-middle); ++i)
	      {
	        if ((i%l)==(middle+1))
	        {
	          if (i==k)
	          {
	            calcGErosion_<InputPeakIterator>((first+middle),last,l,g,false);
	          }
	          else
	          {
	            calcGErosion_<InputPeakIterator>((first+middle),last,l,g,true);
	          }
	          m=0;
	        }
	        if ((i%l)==middle && (i>middle) )
	        {
	          if (i>k)
	          {
	            calcHErosion_<InputPeakIterator>((first+middle),l,h,false);
	          }
	          else
	          {
	            calcHErosion_<InputPeakIterator>((first+middle),l,h,true);
	          }
	          n=0;
	        }
	
	        it->setIntensity(std::min(g[m],h[n]));
	        it->setPosition(first->getPosition());		
	        ++it;
	        ++first;
	        ++m;
	        ++n;
	      }  
	
	      for (i=0; i<middle; ++i)
	      {
	        it->setIntensity(0);
	        it->getPosition() = first->getPosition();
	        ++it;
	        ++first;
	      }
	
	      delete [] g;
	      delete [] h;
	    }
	
	
	
	
	  protected:
	    ///The length of the structuring element.
	    DoubleReal struc_size_;
	   
	    
	    virtual void updateMembers_()
	    {
	      struc_size_ = (DoubleReal)param_.getValue("struc_elem_length"); 
	    }
	
	    /// Subtracted the intensities of all data points in [first, last] from the intensities in result
	    template < typename InputPeakIterator, typename OutputPeakContainer >
	    inline void minusIntensities_(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& result)
	    {
	      typename OutputPeakContainer::iterator it = result.begin();
	      while (first != last)
	      {
	        it->setIntensity(std::max(0.0f, first->getIntensity() - it->getIntensity()));
	        ++first;
	        ++it;
	      }
	    }
	
	    /// Compute the auxiliary fields g and h for the erosion
	    template < typename InputPeakIterator >
	    void calcGErosion_(InputPeakIterator first, InputPeakIterator last, int l, DoubleReal* g, bool b)
	    {
	      int i,j;
	
	      if (b)
	      {
	        for (j=0; j<l; ++j)
	        {
	          if (first < last)
	          {
	            if (j==0)
	            {
	              g[j]=first->getIntensity();
	            }
	            else
	            {
	              g[j]=std::min((DoubleReal)(first->getIntensity()),g[j-1]);
	            }
	            ++first;
	          }
	          else
	          {
	            break;
	          }
	        }
	      }
	      else
	      {
	        j=0;
	        while (first!=last)
	        {
	          if (j==0)
	          {
	            g[j]=first->getIntensity();
	          }
	          else
	          {
	            g[j]=std::min((DoubleReal)(first->getIntensity()),g[j-1]);
	          }
	          ++first;
	          ++j;
	        }
	
	        for (i=j; i<l; ++i)
	        {
	          g[i]=0;
	        }
	      }
	    }
	
	
	    template < typename InputPeakIterator >
	    void calcHErosion_(InputPeakIterator first, int l, DoubleReal* h, bool b)
	    {
	      int j;
	      if (b)
	      {
	        for (j=l-1; j>=0; --j)
	        {
	          if (j==(l-1))
	          {
	            h[j]=first->getIntensity();
	          }
	          else
	          {
	            h[j]=std::min((DoubleReal)(first->getIntensity()),h[j+1]);
	          }
	          if (j > 0)
	          {
	            --first;
	          }
	        }
	      }
	      else
	      {
	        for (j=0;j<l;++j)
	        {
	          h[j]=0;
	        }
	      }
	    }
	
	
	    /// Compute the auxiliary fields g and h for the dilatation
	    template < typename InputPeakIterator >
	    void calcGDilatation_(InputPeakIterator first, InputPeakIterator last, int l, DoubleReal* g, bool b)
	    {
	      int i,j;
	
	      if (b)
	      {
	        for (j=0; j<l; ++j)
	        {
	          if (first < last)
	          {
	            if (j==0)
	            {
	              g[j]=first->getIntensity();
	            }
	            else
	            {
	              g[j]=std::max((DoubleReal)(first->getIntensity()),g[j-1]);
	            }
	            ++first;
	          }
	          else
	          {
	            break;
	          }
	        }
	      }
	      else
	      {
	        j=0;
	        while (first!=last)
	        {
	          if (j==0)
	          {
	            g[j]=first->getIntensity();
	          }
	          else
	          {
	            g[j]=std::max((DoubleReal)(first->getIntensity()),g[j-1]);
	          }
	          ++first;
	          ++j;
	        }
	        for (i=j; i<l; ++i)
	        {
	          g[i]=g[j-1];
	        }
	      }
	    }
	
	
	    template < typename InputPeakIterator >
	    void calcHDilatation_(InputPeakIterator first, InputPeakIterator last, int l, DoubleReal* h, bool b)
	    {
	      int j;
	
	      if (b)
	      {
	        for (j=l-1; j>=0; --j)
	        {
	          if (j==(l-1))
	          {
	            h[j]=last->getIntensity();
	          }
	          else
	          {
	            h[j]=std::max((DoubleReal)(last->getIntensity()),h[j+1]);
	          }
	          if (j > 0)
	          {
	            --last;
	          }
	        }
	      }
	      else
	      {
	        j=(last-first)-1;
	        h[j]=(--last)->getIntensity();
	        while (last!=first)
	        {
	          --j;
	          h[j]=std::max((DoubleReal)(first->getIntensity()),h[j+1]);;
	          --last;
	        }
	      }
	    }

  };

}

#endif
