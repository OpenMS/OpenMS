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
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#include <string>
#include <vector>
#include <map>
#include <math.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>

namespace OpenMS
{

	using namespace std;

	BackgroundIntensityBin& BackgroundIntensityBin::operator =(const BackgroundIntensityBin& bib)
	{
		if (this == &bib)
			return *this;

		this->mzCoord_ = bib.mzCoord_;
		this->trCoord_ = bib.trCoord_;
		this->zCoord_ = bib.zCoord_;
		this->mean_ = bib.mean_;
		this->intensityMap_ = bib.intensityMap_;
		this->intensityHist_ = bib.intensityHist_;

		return *this;
	}

	BackgroundIntensityBin::BackgroundIntensityBin(const BackgroundIntensityBin& bib)
	{
		this->mzCoord_ = bib.mzCoord_;
		this->trCoord_ = bib.trCoord_;
		this->zCoord_ = bib.zCoord_;
		this->mean_ = bib.mean_;
		this->intensityMap_ = bib.intensityMap_;
		this->intensityHist_ = bib.intensityHist_;
	}

	BackgroundIntensityBin::BackgroundIntensityBin(double mz, double tr)
	{
		mzCoord_ = mz;
		trCoord_ = tr;
		zCoord_ = -1;
	}

// check if a peak belongs to this intenity bin
	bool BackgroundIntensityBin::checkBelonging(MSPeak* peak)
	{

		// check charge state:
		if (zCoord_ != -1)
		{
			if (peak->get_charge_state() != zCoord_)
			{
				return false;
			}
		}

		double TR_BINS = SuperHirnParameters::instance()->getBackgroundIntensityBinsTR();

		// check tr:
		double tr = peak->get_retention_time();
		if ((tr < (trCoord_ - TR_BINS / 2.0)) || (tr > (trCoord_ + TR_BINS / 2.0)))
		{
			return false;
		}

		double MZ_BINS = SuperHirnParameters::instance()->getBackgroundIntensityBinsMZ();

		double mz = peak->get_MZ();

		if ((mz < (mzCoord_ - MZ_BINS / 2.0)) || (mz > (mzCoord_ + MZ_BINS / 2.0)))
		{
			return false;
		}

		addIntensity(peak->get_intensity());
		peak = NULL;
		return true;
	}

	void BackgroundIntensityBin::addMSPeak(MSPeak* peak)
	{
		addIntensity(peak->get_intensity());
		peak = NULL;
	}

	void BackgroundIntensityBin::addIntensity(double intens)
	{
		intensityMap_.push_back(intens);
	}

// copied from simple_math
	double simple_math_WEIGHTED_AVERAGE(map<double, double>* IN)
	{

		double AVERAGE = 0;
		double TOT_WEIGHT = 0;

		if (IN->size() > 1)
		{

			map<double, double>::iterator START = IN->begin();
			while (START != IN->end())
			{
				TOT_WEIGHT += (*START).second;
				AVERAGE += ((*START).first * (*START).second );
				START++;
			}

			return AVERAGE / TOT_WEIGHT;
		}
		else
		{
			return (*IN->begin()).first;
		}
	}

// process collected intensities in the map
	void BackgroundIntensityBin::processIntensities()
	{

		computeIntensityHist();

		if (!intensityHist_.empty())
		{
			mean_ = simple_math_WEIGHTED_AVERAGE(&intensityHist_);
		}
		else
		{
			mean_ = 0;
		}
	}

// compute an intensity histogram
	void BackgroundIntensityBin::computeIntensityHist()
	{

		double constraint = SuperHirnParameters::instance()->getBackgroundIntensityBinsIntens();

		// insert into the histogram map
		vector<double>::iterator P = intensityMap_.begin();
		while (P != intensityMap_.end())
		{

			// intensity to bin:
			double intens = (*P);

			// find a key:
			map<double, double>::iterator F = intensityHist_.lower_bound(intens);
			if (F != intensityHist_.end())
			{

				// check this one:
				map<double, double>::iterator check = F;
				double mainLow = fabs(check->first - intens);
				double deltaHigh = 1000000;
				if (check != intensityHist_.begin())
				{
					check--;
					deltaHigh = fabs(check->first - intens);
					if (mainLow > deltaHigh)
					{
						mainLow = deltaHigh;
						F = check;
					}
				}
				if (mainLow > constraint)
				{
					F = intensityHist_.end();
				}
				else
				{
					F->second += 1.0;
				}
			}

			if (F == intensityHist_.end())
			{
				intensityHist_.insert(make_pair(intens, 1.0));
			}

			P++;
		}

		// filter out bins of only 1 counts:
		map<double, double>::iterator F = intensityHist_.begin();
		while (F != intensityHist_.end())
		{

			if (F->second == SuperHirnParameters::instance()->getBackgroundIntensityBinsMinBinCount())
			{
				intensityHist_.erase(F++);
			}
			else
			{
				F++;
			}
		}
	}

	BackgroundIntensityBin::~BackgroundIntensityBin()
	{
		intensityMap_.clear();
	}

}
