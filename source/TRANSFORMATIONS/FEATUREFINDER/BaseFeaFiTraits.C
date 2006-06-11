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
// $Id: BaseFeaFiTraits.C,v 1.17 2006/03/08 21:32:32 j-joachim Exp $
// $Author: j-joachim $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseFeaFiTraits.h>

// all from BaseFeaFiTraits derived classes
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleFeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedFeaFiTraits.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FastFeaFiTraits.h>

namespace OpenMS
{
	void BaseFeaFiTraits::registerChildren()
	{
		Factory<BaseFeaFiTraits>::registerProduct(SimpleFeaFiTraits::getName(), &SimpleFeaFiTraits::create);
		Factory<BaseFeaFiTraits>::registerProduct(ExtendedFeaFiTraits::getName(), &ExtendedFeaFiTraits::create);
	}

	BaseFeaFiTraits::BaseFeaFiTraits() 
		: FactoryProduct(), debug_(0),
		  instance_(1) 
	{
		debug_stream_ = &std::cout;
	}

	BaseFeaFiTraits::BaseFeaFiTraits(const BaseFeaFiTraits& source)
		: FactoryProduct(source), debug_(0),
		  debug_stream_(source.debug_stream_), instance_(1)
	{	}

	BaseFeaFiTraits::~BaseFeaFiTraits() {}

	BaseFeaFiTraits& BaseFeaFiTraits::operator = (const BaseFeaFiTraits& source)
	{
		FactoryProduct::operator = (source);
		debug_									 = source.debug_;
		debug_stream_            = source.debug_stream_;
		instance_                = source.instance_;
		return *this;
	}

	void BaseFeaFiTraits::setSeeders(std::vector<BaseSeeder*> seeders)
	{
		std::vector<BaseSeeder*>::iterator it = seeders.begin();
		while (it != seeders.end()) { (*it++)->setTraits(this); }
		seeders_ = seeders;
		for (UnsignedInt i = 0; i<seeders.size(); i++)
			seeders[i]->setTraits(this);
	}

	void BaseFeaFiTraits::setExtenders(std::vector<BaseExtender*> extenders)
	{
			std::vector<BaseExtender*>:: iterator it = extenders.begin();
			while (it != extenders.end()) { (*it++)->setTraits(this); }
      extenders_ = extenders;
	}
	
	void BaseFeaFiTraits::setDebugStream(std::ostream* os) { debug_stream_ = os; }
	 
	std::ostream* BaseFeaFiTraits::getDebugStream() { return debug_stream_; }
		
	void BaseFeaFiTraits::setDebugLevel(UnsignedInt lvl) { debug_ = lvl; }
  
  const UnsignedInt& BaseFeaFiTraits::getDebugLevel() const { return debug_; }
  
  UnsignedInt& BaseFeaFiTraits::getDebugLevel() { return debug_; }
  
  void BaseFeaFiTraits::setInstanceId(String num) { instance_ = num; }
  
  const String& BaseFeaFiTraits::getInstanceId() const  { return instance_; }
   
  String& BaseFeaFiTraits::getInstanceId() { 	return instance_; }
	
	void BaseFeaFiTraits::setFitters(std::vector<BaseModelFitter*> fitters)
	{
		std::vector<BaseModelFitter*>:: iterator it = fitters.begin();
		while (it != fitters.end()) { (*it++)->setTraits(this); }
    fitters_ = fitters;
	}
		
	const BaseFeaFiTraits::ConvexHullType& BaseFeaFiTraits::calculateConvexHull(const IndexSet& set)
	{
		const double PRECISION = 0.0001;
		convex_hull_.clear();
		if (set.size()<3)
			return convex_hull_;

		// keep track of already in hull included peaks to avoid unnecessary computations of triangle area
		std::map<UnsignedInt, bool> isIncluded;

		CoordinateType min_mz = std::numeric_limits<CoordinateType>::max();
		IndexSet::const_iterator min = set.begin();

		// Find peak with minimal mz to start wrapping
		for (IndexSet::const_iterator it = set.begin(); it!=set.end(); ++it)
		{
			if (getPeakMz(*it) < min_mz)
			{
				min_mz = getPeakMz(*it);
				min = it;
			}
			isIncluded[*it] = false;
		}
		convex_hull_.push_back( getPeak(*min).getPosition() );

		// Hull peaks denoting current hull line
		IndexSet::const_iterator hull_peak1 = min;
		IndexSet::const_iterator start = set.begin();
		if (start==min) ++start;  // don't start at "min" because of while-condition
		IndexSet::const_iterator hull_peak2 = start;

		while (hull_peak2!=min)
		{
			bool found_any = false;
			for (IndexSet::const_iterator it = set.begin(); it!=set.end(); ++it)
			{
				// skip if already used
				if (isIncluded[*it] || it==hull_peak1 || it==hull_peak2) continue;

				found_any = true;
				// "it" lies to the right of the line [hull_peak1,hull_peak2]
				double area = triangleArea_(hull_peak1,hull_peak2,it);
				if (area>-PRECISION)
				{
					// area almost 0 -> collinear points
					// -> avoid consecutive peaks with equal mz or rt coordinate
					if (fabs(area)<PRECISION)
					{
						double mz1 = getPeakMz(*hull_peak1);
						double mz2 = getPeakMz(*hull_peak2);
						double mz3 = getPeakMz(*it);
						double rt1 = getPeakRt(*hull_peak1);
						double rt2 = getPeakRt(*hull_peak2);
						double rt3 = getPeakRt(*it);
			      if ( 	( fabs(mz2-mz3)<PRECISION && fabs(rt2-rt1) > fabs(rt3-rt1) )
								||( fabs(rt2-rt3)<PRECISION && fabs(mz2-mz1) > fabs(mz3-mz1) ))
						{
							isIncluded[*it] = true;
							continue;
						}
					}
					hull_peak2 = it;  // "it" becomes new hull peak
				}
			}

			if (!found_any){
				hull_peak2 = min; // no available peaks anymore
				continue;
			}

			if (hull_peak2 == min) continue;  // finish loop
			isIncluded[*hull_peak2] = true;

			// continue wrapping
			hull_peak1 = hull_peak2;
			// hull_peak2 satisfies the contition: all peaks lie to the left of [hull_peak1,hull_peak2]
			convex_hull_.push_back( getPeak(*hull_peak2).getPosition() );


			start = set.begin();
			if (start==min) ++start;  // don't start at "min" because of while-condition
			hull_peak2 = start;
		}

		return convex_hull_;
	}

}
