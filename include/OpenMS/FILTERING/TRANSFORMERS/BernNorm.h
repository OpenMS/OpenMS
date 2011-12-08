// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_BERNNORM_H
#define OPENMS_FILTERING_TRANSFORMERS_BERNNORM_H

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <map>

namespace OpenMS
{
  /**
  	@brief BernNorm scales the peaks by ranking them and then scaling them according to rank.

  	For exact formula look in Bioinformatics, Aug 2004; 20: i49 - i54

		@improvement read paper and try to confirm implementation (andreas)

		@htmlinclude OpenMS_BernNorm.parameters

		@ingroup SpectraPreprocessers
  */
	class OPENMS_DLLAPI BernNorm
		: public DefaultParamHandler
	{
	public:

	// @name Constructors and Destructors
	//@{
	/// default constructor
	BernNorm();

	/// copy constructor
	BernNorm(const BernNorm& source);

	/// destructor
	virtual ~BernNorm();
	//@}

	// @name Operators
	// @{
	/// assignment operator
	BernNorm& operator=(const BernNorm& source);
	//@}

	// @name Accessors
	// @{

	///
	template <typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)
	{
		typedef typename SpectrumType::Iterator Iterator;
		typedef typename SpectrumType::ConstIterator ConstIterator;

		c1_ = (DoubleReal)param_.getValue("C1");
		c2_ = (DoubleReal)param_.getValue("C2");
		th_ = (DoubleReal)param_.getValue("threshold");

		spectrum.sortByPosition();

		// find highest peak and ranking
		double maxint = 0;
		std::map<double, Size> peakranks;
		for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
		{
		peakranks[it->getIntensity()] = 0;
		if (it->getIntensity() > maxint)
		{
			maxint = it->getIntensity();
		}
		}
		UInt rank = 0;
		for (std::map<double, Size>::reverse_iterator mit = peakranks.rbegin(); mit != peakranks.rend(); ++mit)
		{
			mit->second = ++rank;
		}

		// find maxmz i.e. significant (> threshold * maxpeak) peak with highest m/z
		double maxmz = 0;
		for (SignedSize i = spectrum.size() -1 ; i >= 0 ; --i)
		{
			if (spectrum[i].getIntensity() > maxint * th_)
			{
				maxmz = spectrum[i].getMZ();
				break;
			}
		}

		// rank
		for (Iterator it = spectrum.begin() ; it != spectrum.end(); )
		{
			double newint = c1_ - (c2_ / maxmz) * peakranks[it->getIntensity()];
			if (newint < 0)
			{
				it = spectrum.erase(it);
			}
			else
			{
				it->setIntensity(newint);
				++it;
			}
		}
		return;
	}

	void filterPeakSpectrum(PeakSpectrum& spectrum);

	void filterPeakMap(PeakMap& exp);
	//TODO reimplement DefaultParamHandler::updateMembers_()

	private:
		DoubleReal c1_;
		DoubleReal c2_;
		DoubleReal th_;

		// @}

	};

} // namespace OpenMS

#endif //OPENMS_FILTERING_TRANSFORMERS_BERNNORM_H
