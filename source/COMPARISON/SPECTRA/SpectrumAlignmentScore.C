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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <math.h>

using namespace std;

namespace OpenMS
{
  SpectrumAlignmentScore::SpectrumAlignmentScore()
    : PeakSpectrumCompareFunctor()
  {
		setName(SpectrumAlignmentScore::getProductName());
		defaults_.setValue("epsilon", 0.3, "Defines the absolut error of the mass spectrometer");
		defaultsToParam_();
  }

  SpectrumAlignmentScore::SpectrumAlignmentScore(const SpectrumAlignmentScore& source)
    : PeakSpectrumCompareFunctor(source)
  {
  }

  SpectrumAlignmentScore::~SpectrumAlignmentScore()
  {
  }

  SpectrumAlignmentScore& SpectrumAlignmentScore::operator = (const SpectrumAlignmentScore& source)
  {
		if (this != &source)
		{
    	PeakSpectrumCompareFunctor::operator = (source);
		}
    return *this;
  }

	double SpectrumAlignmentScore::operator () (const PeakSpectrum& spec) const
	{
		return operator () (spec, spec);
	}

  double SpectrumAlignmentScore::operator () (const PeakSpectrum& s1, const PeakSpectrum& s2) const
  {			
		const double epsilon = (double)param_.getValue("epsilon");

		SpectrumAlignment aligner;
		Param p;
		p.setValue("epsilon", epsilon);
		aligner.setParameters(p);

		vector<pair<UInt, UInt> > alignment;
		aligner.getSpectrumAlignment(alignment, s1, s2);

		double score(0), sum(0), sum1(0), sum2(0);
		for (PeakSpectrum::ConstIterator it1 = s1.begin(); it1 != s1.end(); ++it1)
		{
			sum1 += it1->getIntensity() * it1->getIntensity();
		}
		
		for (PeakSpectrum::ConstIterator it1 = s2.begin(); it1 != s2.end(); ++it1)
		{
			sum2 += it1->getIntensity() * it1->getIntensity();
		}
		
		for (vector<pair<UInt, UInt> >::const_iterator it = alignment.begin(); it != alignment.end(); ++it)
		{
			sum += sqrt(s1.getContainer()[it->first].getIntensity() * s2.getContainer()[it->second].getIntensity());
		}

    score = sum / (sqrt(sum1 * sum2));

    return score;
	}

}
