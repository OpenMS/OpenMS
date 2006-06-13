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
// $Id: MarkerMower.h,v 1.5 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_MARKERMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_MARKERMOWER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <vector>

namespace OpenMS
{
  /**
  	@brief MarkerMower uses PeakMarker to find peaks, those that are not marked get removed<br>
  */
  class MarkerMower : public PreprocessingFunctor
  {
  public:
    /// standard constructor
    MarkerMower();

    /// copy constructor
    MarkerMower(const MarkerMower& source);

    /// destructor
    ~MarkerMower();

    /// assignment operator
    MarkerMower& operator=(const MarkerMower& source);

    static FactoryProduct* create() { return new MarkerMower();}
    void operator()(MSSpectrum< DPeak<1> >&) const;
    String info() const;

		static const String getName()
		{
			return "MowerMarker";
		}

    /// insert new Marker
    void insertmarker(PeakMarker*);
  private:
    static const String info_;
    /**
    used Markers
    */
    std::vector<PeakMarker*> markers_;
  };
}
#endif // OPENMS_COMPARISON_CLUSTERING_MARKERMOWER_H
