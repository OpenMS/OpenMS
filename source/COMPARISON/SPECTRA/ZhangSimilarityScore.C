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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>

#include <boost/math/special_functions/erf.hpp>

#include <cmath>

using namespace std;

namespace OpenMS
{
  ZhangSimilarityScore::ZhangSimilarityScore()
    : PeakSpectrumCompareFunctor()
  {
		setName(ZhangSimilarityScore::getProductName());
		defaults_.setValue("tolerance", 0.2, "defines the absolut (in Da) or relative (in ppm) tolerance");
		defaults_.setValue("is_relative_tolerance", "false", "If set to true, the tolerance is interpreted as relative");
		defaults_.setValidStrings("is_relative_tolerance", StringList::create("true,false"));
		defaults_.setValue("use_linear_factor", "false", "if true, the intensities are weighted with the relative m/z difference");
    defaults_.setValidStrings("use_linear_factor", StringList::create("true,false"));
		defaults_.setValue("use_gaussian_factor", "false", "if true, the intensities are weighted with the relative m/z difference using a gaussian");
    defaults_.setValidStrings("use_gaussian_factor", StringList::create("true,false"));
		defaultsToParam_();
  }

  ZhangSimilarityScore::ZhangSimilarityScore(const ZhangSimilarityScore& source)
    : PeakSpectrumCompareFunctor(source)
  {
  }

  ZhangSimilarityScore::~ZhangSimilarityScore()
  {
  }

  ZhangSimilarityScore& ZhangSimilarityScore::operator = (const ZhangSimilarityScore& source)
  {
		if (this != &source)
		{
    	PeakSpectrumCompareFunctor::operator = (source);
		}
    return *this;
  }

	double ZhangSimilarityScore::operator () (const PeakSpectrum& spec) const
	{
		return operator () (spec, spec);
	}
	
  double ZhangSimilarityScore::operator () (const PeakSpectrum& s1, const PeakSpectrum& s2) const
  {
		const double tolerance = (double)param_.getValue("tolerance");
		bool use_linear_factor = param_.getValue("use_linear_factor").toBool();
		bool use_gaussian_factor = param_.getValue("use_gaussian_factor").toBool();
    double score(0), sum(0), sum1(0), sum2(0)/*, squared_sum1(0), squared_sum2(0)*/;

    for (PeakSpectrum::ConstIterator it1 = s1.begin(); it1 != s1.end(); ++it1)
    {
      sum1 += it1->getIntensity();
			/*
      for (PeakSpectrum::ConstIterator it2 = s1.begin(); it2 != s1.end(); ++it2)
      {
        if (abs(it1->getPosition()[0] - it2->getPosition()[0]) <= 2 * tolerance)
        {
          squared_sum1 += it1->getIntensity() * it2->getIntensity();
        }
      }*/
    }

/*
		UInt i_left(0);
		for (Size i = 0; i != s1.size(); ++i)
		{
			sum1 += s1[i].getIntensity();
			for (Size j = i_left; j != s1.size(); ++j)
			{
				double pos1(s1[i].getPosition()[0]), pos2(s1[j].getPosition()[0]);
				if (abs(pos1 - pos2) <= 2 * tolerance)
				{
					squared_sum1 += s1[i].getIntensity() * s1[j].getIntensity();
				}
				else
				{
					if (pos2 > pos1)
					{
						break;
					}
					else
					{
						i_left = i;
					}
				}
			}
		}*/

/*
    i_left = 0;
    for (Size i = 0; i != s2.size(); ++i)
    {
      sum2 += s2[i].getIntensity();
      for (Size j = i_left; j != s2.size(); ++j)
      {
        double pos1(s2[i].getPosition()[0]), pos2(s2[j].getPosition()[0]);
        if (abs(pos1 - pos2) <= 2 * tolerance)
        {
          squared_sum1 += s2[i].getIntensity() * s2[j].getIntensity();
        }
        else
        {
          if (pos2 > pos1)
          {
            break;
          }
          else
          {
            i_left = i;
          }
        }
      }
    }*/

    for (PeakSpectrum::ConstIterator it1 = s2.begin(); it1 != s2.end(); ++it1)
    {
      sum2 += it1->getIntensity();
			/*
      for (PeakSpectrum::ConstIterator it2 = s2.begin(); it2 != s2.end(); ++it2)
      {
        if (abs(it1->getPosition()[0] - it2->getPosition()[0]) <= 2 * tolerance)
        {
          squared_sum2 += it1->getIntensity() * it2->getIntensity();
        }
      }
			*/
    }

		Size j_left(0);
		for (Size i = 0; i != s1.size(); ++i)
		{
			for (Size j = j_left; j != s2.size(); ++j)
			{
				double pos1(s1[i].getMZ()), pos2(s2[j].getMZ());
				if (fabs(pos1 - pos2) < tolerance)
				{
					//double factor((tolerance - fabs(pos1 - pos2)) / tolerance);
					double factor = 1.0;
					
					if (use_linear_factor || use_gaussian_factor)
					{
						factor = getFactor_(tolerance, fabs(pos1 - pos2), use_gaussian_factor);
					}
					sum += sqrt(s1[i].getIntensity() * s2[j].getIntensity() * factor);
				}
				else
				{
					if (pos2 > pos1)
					{
						break;
					}
					else
					{
						j_left = j;
					}
				}
			}
		}


		/*
    for (PeakSpectrum::ConstIterator it1 = s1.begin(); it1 != s1.end(); ++it1)
    {
      for (PeakSpectrum::ConstIterator it2 = s2.begin(); it2 != s2.end(); ++it2)
      {
        if (abs(it1->getPosition()[0] - it2->getPosition()[0]) <= 2 * tolerance)
        {
          sum += sqrt(it1->getIntensity() * it2->getIntensity());
        }
      }
    }*/

    score = sum / (sqrt(sum1 * sum2));

    return score;
	
	}


	double ZhangSimilarityScore::getFactor_(double mz_tolerance, double mz_difference, bool is_gaussian) const
  {
    double factor(0.0);

    if (is_gaussian)
    {
			static const DoubleReal denominator = mz_tolerance * 3.0 * sqrt(2.0);
			factor = boost::math::erfc(mz_difference / denominator);
			//cerr << "Factor: " << factor << " " << mz_tolerance << " " << mz_difference << endl;
    }
    else
    {
      factor = (mz_tolerance - mz_difference) / mz_tolerance;
    }
    return factor;
  }

	
}
