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
// $Id: PeakDensityFilter.C,v 1.4 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/PeakDensityFilter.h>
#include <cmath>

using namespace std;

namespace OpenMS
{
  //const String PeakDensityFilter::info_ = "Scores the PeakDensity according to some Function";
  
  PeakDensityFilter::PeakDensityFilter()
    :FilterFunctor()
  {
		name_ = PeakDensityFilter::getName();
  }

  PeakDensityFilter::PeakDensityFilter(const PeakDensityFilter& source )
    :FilterFunctor(source)
  {
		name_ = source.getName();
  }
    
  PeakDensityFilter& PeakDensityFilter::operator = (const PeakDensityFilter& source)
  {
    FilterFunctor::operator = (source);
    return *this;
  }
  
  PeakDensityFilter::~PeakDensityFilter()
  {
  }
/*
  String PeakDensityFilter::info() const
  {
    return info_;
  }

  vector<double> PeakDensityFilter::operator()(const ClusterSpectrum& cspec)
  {
    int nrpeaks = cspec.getSpec().getContainer().size();
    double size = cspec.getParentMass();
    double density = nrpeaks/size;
    vector<double> result;
    result.push_back(density);
    return result;
  }
*/
  double PeakDensityFilter::gaussfit_(double x)
  {
    //parameters obtained by a fitting
    //values of a good sample (high sequest scores)
    //versus a gauss function (using gnuplot)
    double b = 0.0934731;
    double c = 0.0458471;
    return exp(-((x-b)*(x-b))/(c*c));
  }
}
