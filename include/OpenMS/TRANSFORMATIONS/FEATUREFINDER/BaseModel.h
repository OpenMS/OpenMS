// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASEMODEL_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/DPeak.h>

namespace OpenMS
{

  /**
	@brief Abstract base class for all D-dimensional models.

	Every derived class has to implement the static functions
	"T* create()" and "const String getProductName()" (see DefaultParamHandler for details)
  */
  template <UInt D>
	class BaseModel
    : public DefaultParamHandler
	{

	 public:

		typedef DoubleReal IntensityType;
		typedef DoubleReal CoordinateType;
		typedef DPosition<D> PositionType;
		typedef typename DPeak<D>::Type PeakType;
		typedef std::vector<PeakType> SamplesType;


		/// Default constructor.
		BaseModel()
			: DefaultParamHandler("BaseModel")
		{
			defaults_.setValue("cutoff",0.0,"Low intensity cutoff of the model.  Peaks below this intensity are not considered part of the model.");
		}

		/// copy constructor
		BaseModel(const BaseModel& source)
			: DefaultParamHandler(source),
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

			DefaultParamHandler::operator = (source);
			cut_off_ = source.cut_off_;

			return *this;
		}

		/// register all derived classes here (implemented in file BaseModel_impl.h)
		static void registerChildren();

		/// access model predicted intensity at position @p pos
		virtual IntensityType getIntensity(const PositionType& pos) const=0;

		/// check if position @p pos is part of the model regarding the models cut-off.
		virtual bool isContained(const PositionType& pos) const
		{
			return getIntensity(pos) >= cut_off_;
		}

		/**@brief Convenience function to set the intensity of a peak to the
		predicted intensity at its current position, calling virtual void
		getIntensity().
		*/
		template <typename PeakType>
		void fillIntensity(PeakType& peak) const
		{
			peak.setIntensity(getIntensity(peak.getPosition()));
		}

		/**@brief Convenience function that applies fillIntensity() to an iterator
		range.
		*/
		template <class PeakIterator>
		void  fillIntensities(PeakIterator begin, PeakIterator end) const
		{
			for ( PeakIterator it = begin; it != end; ++it )
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
			param_.setValue("cutoff",cut_off_);
		}

		/// get reasonable set of samples from the model (i.e. for printing)
		virtual void getSamples(SamplesType& cont) const =0;

		/// fill stream with reasonable set of samples from the model (i.e. for printing)
		virtual void getSamples(std::ostream& os)
		{
			SamplesType samples;
			getSamples(samples);
			for (typename SamplesType::const_iterator it=samples.begin();it!=samples.end(); ++it)
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
