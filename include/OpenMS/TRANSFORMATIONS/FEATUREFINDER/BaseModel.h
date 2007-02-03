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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H

#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/KERNEL/DPeakArray.h>

namespace OpenMS
{
  /** 
  	@brief Abstract base class for all D-dimensional models.
		
		Every derived class has to implement the static functions
    "T* create()" and "const String getProductName()" (see FactoryProduct for details)
		
		@ingroup FeatureFinder
   */
  template <UnsignedInt D, typename Traits = KernelTraits>
    class BaseModel
    : public FactoryProduct
    {

      public:
	      typedef int Flag;
	      typedef std::vector<Flag> Flags;
	
	      typedef typename DPeak<D,Traits>::IntensityType IntensityType;
		  	typedef typename DPeak<D,Traits>::CoordinateType CoordinateType;
	      typedef DPosition<D,Traits> PositionType;
	      typedef DPeak<D,Traits> PeakType;
		  	typedef DPeakArray<D, DPeak<D,Traits> > SamplesType;
	
	
	      /// Default constructor. 
	      BaseModel()
					: FactoryProduct("BaseModel")
				{
					defaults_.setValue("cutoff",0.0);
				}
	
	      /// copy constructor 
	      BaseModel(const BaseModel& source)
					: FactoryProduct(source),
						cut_off_(source.cut_off_)
				{
				}
	
	      /// Destructor 
	      virtual ~BaseModel()
	      {	
	      }
	
	      /// assignment operator
	      virtual BaseModel& operator = (const BaseModel& source)
				{
					if (&source == this) return *this;
					
					FactoryProduct::operator = (source);
					cut_off_ = source.cut_off_;
					
					return *this;
				}
	
	      /// register all derived classes here
	      static void registerChildren();
				
	      /// acess model predicted intensity at position @p pos
	      virtual IntensityType getIntensity(const PositionType& pos) const=0;
	      
	      /// check if position @p pos is part of the model regarding the models cut-off.
	      virtual bool isContained(const PositionType& pos) const
				{
					return getIntensity(pos) >= cut_off_;
				}
	
	      /// set DPeaks intensity to model predicted intensity.
	      virtual void  fillIntensity(PeakType& peak) const
				{
					peak.setIntensity( getIntensity(peak.getPosition()) );
				}
	
	      //// set DPeaks intensity to model predicted intensity.
				template <class PeakIterator>
	      void  fillIntensities(PeakIterator beg, PeakIterator end) const
				{
		    	for (PeakIterator it=beg; it!=end; it++)
					{
						fillIntensity(*it);
					}
				}
	
				/// get cutoff value
				virtual IntensityType getCutOff() const	
				{	
					return cut_off_;	
				}
	
				///	set cutoff value
				virtual void setCutOff(IntensityType cut_off)
				{
					cut_off_ = cut_off;
					param_.setValue("cutoff",(double)cut_off_);
				}
				
				/// get reasonable set of samples from the model (i.e. for printing)
				virtual void getSamples(SamplesType& cont) const =0;
	
				/// fill stream with reasonable set of samples from the model (i.e. for printing)
				virtual void getSamples(std::ostream& os)
				{
					SamplesType samples;
					getSamples(samples);
					for (typename SamplesType::ConstIterator it=samples.begin();it!=samples.end(); ++it)
					{
						os << *it << std::endl;
					}
				}
	
			protected:
				IntensityType cut_off_;
	
				virtual void updateMembers_()
				{
					cut_off_ = (double)param_.getValue("cutoff");
				}
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H
