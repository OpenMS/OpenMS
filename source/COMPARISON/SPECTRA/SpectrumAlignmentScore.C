// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <cmath>

using namespace std;

namespace OpenMS
{
  SpectrumAlignmentScore::SpectrumAlignmentScore()
    : PeakSpectrumCompareFunctor()
  {
		setName(SpectrumAlignmentScore::getProductName());
		//defaults_.setValue("epsilon", 0.3, "Defines the absolut error of the mass spectrometer");
		defaults_.setValue("tolerance", 0.3, "Defines the absolut (in Da) or relative (in ppm) tolerance");
		defaults_.setValue("is_relative_tolerance", "false", "if true, the tolerance value is interpreted as ppm");
		defaults_.setValidStrings("is_relative_tolerance", StringList::create("true,false"));
		defaults_.setValue("use_linear_factor", "false", "if true, the intensities are weighted with the relative m/z difference");
		defaults_.setValidStrings("use_linear_factor", StringList::create("true,false"));
		defaults_.setValue("use_gaussian_factor", "false", "if true, the intensities are weighted with the relative m/z difference using a gaussian");
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
		const double tolerance = (double)param_.getValue("tolerance");
		bool is_relative_tolerance = param_.getValue("is_relative_tolerance").toBool();
		bool use_linear_factor = param_.getValue("use_linear_factor").toBool();
		bool use_gaussian_factor = param_.getValue("use_gaussian_factor").toBool();

		if (use_linear_factor && use_gaussian_factor)
		{
			cerr << "Warning: SpectrumAlignmentScore, use either 'use_linear_factor' or 'use_gaussian_factor'!" << endl;
		}

		SpectrumAlignment aligner;
		Param p;
		p.setValue("tolerance", tolerance);
		p.setValue("is_relative_tolerance", (String)param_.getValue("is_relative_tolerance"));
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
			//double factor(0.0);
			//factor = (epsilon - fabs(s1[it->first].getPosition()[0] - s2[it->second].getPosition()[0])) / epsilon;
			double mz_tolerance(tolerance);

			if (is_relative_tolerance)
			{
				mz_tolerance = mz_tolerance * s1[it->first].getPosition()[0] / 10e6;
			}
	
			double mz_difference(fabs(s1[it->first].getPosition()[0] - s2[it->second].getPosition()[0]));
			double factor = 1.0;
			
			if (use_linear_factor || use_gaussian_factor)
			{
				factor = getFactor_(mz_tolerance, mz_difference, use_gaussian_factor);
			}
			sum += sqrt(s1[it->first].getIntensity() * s2[it->second].getIntensity() * factor);
		}

    score = sum / (sqrt(sum1 * sum2));

    return score;
	}

	double SpectrumAlignmentScore::getFactor_(double mz_tolerance, double mz_difference, bool is_gaussian) const
	{
		double factor(0.0);

		if (is_gaussian)
		{
			throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			// to be implemented
		}
		else
		{
			factor = (mz_tolerance - mz_difference) / mz_tolerance;
		}
		return factor;
	}

}
