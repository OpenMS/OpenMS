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
// $Id: BaseModel.h,v 1.19 2006/04/25 12:41:44 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------


#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H


#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/KERNEL/DPeakArray.h>
#include <vector>

namespace OpenMS
{
  /** @brief Abstract base class for all D-dimensional models.
		
			Every derived class has to implement the static functions
      "T* create()" and "const String getName()" (see FactoryProduct for details)

			Every derived class should implement a getParm method to allow lazy parameter
			setting i.e. each time the getParam is called the param_ object is updated rather than
			when the set methods are called allowing fast resetting of the model during fitting.
			
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
      typedef DPosition<D,Traits> PositionType;
      typedef DPeak<D,Traits> PeakType;
			typedef DPeakArray<D, DPeak<D,Traits> > SamplesType;


      /// standard constructor. 
      BaseModel()
				: FactoryProduct()
			{
				setCutOff(0.0);
				defaults_.setValue("cutoff",0.0);
			}

      /// copy constructor 
      BaseModel(const BaseModel& source)
				: FactoryProduct(source)
			{
				setCutOff(source.cut_off_);
			}

      /// Destructor 
      virtual ~BaseModel(){}

      /// assignment operator
      virtual BaseModel& operator = (const BaseModel& source)
			{
				FactoryProduct::operator = (source);
				setCutOff(source.cut_off_);
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
			virtual IntensityType getCutOff() const	{	return cut_off_;	}


			///	set cutoff value
			virtual void setCutOff(IntensityType cut_off)
			{
				cut_off_ = cut_off;
			}

			virtual void setParam(const Param& param)
			{
				FactoryProduct::setParam(param);
				cut_off_ = static_cast<double>(param_.getValue("cutoff"));
			}

    	/// get parameters (const access)
    	virtual const Param& getParam() const
			{
				param_.setValue("cutoff",(double)cut_off_);
				return FactoryProduct::getParam();
			}

    	/// get parameters
    	virtual Param& getParam()
			{
				param_.setValue("cutoff",(double)cut_off_);
				return FactoryProduct::getParam();
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
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H
