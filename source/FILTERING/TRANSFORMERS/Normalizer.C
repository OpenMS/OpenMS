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
// $Id: Normalizer.C,v 1.5 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>

#include <cmath>

using namespace std;
namespace OpenMS
{
  Normalizer::Normalizer()
    : PreprocessingFunctor()
  {
		name_ = Normalizer::getName();
    defaults_.setValue("windows", 10); // only applicable for max normalizing, run twice for other methods
		param_ = defaults_;
  }

  Normalizer::Normalizer(const Normalizer& source)
    : PreprocessingFunctor(source)
  {
		name_ = source.getName();
  }

  Normalizer::~Normalizer()
  {
  }

  Normalizer& Normalizer::operator=(const Normalizer& source)
  {
    PreprocessingFunctor::operator=(source);
    return *this;
  }

/*
  void Normalizer::operator()(MSSpectrum< DPeak<1> >& spec) const
  {
    vector<double> max = vector<double>((unsigned int)param_.getValue("windows"));
    double minmz = 10000;
    double maxmz = 0;
    for (MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end();++it )
    {
      if ( it->getPosition()[0] < minmz ) minmz = it->getPosition()[0];
      else if ( it->getPosition()[0] > maxmz ) maxmz = it->getPosition()[0];
    }
    for (MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end();++it )
    {
      uint pos = (uint) ( ( it->getPosition()[0] - minmz ) / ( maxmz - minmz ) *  max.size() );
      if ( pos == max.size() ) --pos;
      if ( max.at(pos) < it->getIntensity() ) max.at(pos) = it->getIntensity(); 
    }
    for (MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end();++it )
    {
      uint pos = (uint) ( ( it->getPosition()[0] - minmz ) / ( maxmz - minmz ) * max.size() );
      if ( pos == max.size() ) --pos;
      it->setIntensity(it->getIntensity()/max[pos]);
    }
  }
*/
}
