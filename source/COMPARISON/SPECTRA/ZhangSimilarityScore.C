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

#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  ZhangSimilarityScore::ZhangSimilarityScore()
    : PeakSpectrumCompareFunctor()
  {
		setName(ZhangSimilarityScore::getProductName());
		defaults_.setValue("epsilon", 0.2);
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
		const double epsilon = (double)param_.getValue("epsilon");
  	//const double epsilon(0.4);
    //const double epsilon(0.2);
    //const double c(0.0004);
    double score(0), sum(0), sum1(0), sum2(0)/*, squared_sum1(0), squared_sum2(0)*/;

    for (PeakSpectrum::ConstIterator it1 = s1.begin(); it1 != s1.end(); ++it1)
    {
      sum1 += it1->getIntensity();
			/*
      for (PeakSpectrum::ConstIterator it2 = s1.begin(); it2 != s1.end(); ++it2)
      {
        if (abs(it1->getPosition()[0] - it2->getPosition()[0]) <= 2 * epsilon)
        {
          squared_sum1 += it1->getIntensity() * it2->getIntensity();
        }
      }*/
    }

/*
		Size i_left(0);
		for (Size i = 0; i != s1.getContainer().size(); ++i)
		{
			sum1 += s1.getContainer()[i].getIntensity();
			for (Size j = i_left; j != s1.getContainer().size(); ++j)
			{
				double pos1(s1.getContainer()[i].getPosition()[0]), pos2(s1.getContainer()[j].getPosition()[0]);
				if (abs(pos1 - pos2) <= 2 * epsilon)
				{
					squared_sum1 += s1.getContainer()[i].getIntensity() * s1.getContainer()[j].getIntensity();
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
    for (Size i = 0; i != s2.getContainer().size(); ++i)
    {
      sum2 += s2.getContainer()[i].getIntensity();
      for (Size j = i_left; j != s2.getContainer().size(); ++j)
      {
        double pos1(s2.getContainer()[i].getPosition()[0]), pos2(s2.getContainer()[j].getPosition()[0]);
        if (abs(pos1 - pos2) <= 2 * epsilon)
        {
          squared_sum1 += s2.getContainer()[i].getIntensity() * s2.getContainer()[j].getIntensity();
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
        if (abs(it1->getPosition()[0] - it2->getPosition()[0]) <= 2 * epsilon)
        {
          squared_sum2 += it1->getIntensity() * it2->getIntensity();
        }
      }
			*/
    }

		Size j_left(0);
		for (Size i = 0; i != s1.getContainer().size(); ++i)
		{
			for (Size j = j_left; j != s2.getContainer().size(); ++j)
			{
				double pos1(s1.getContainer()[i].getPosition()[0]), pos2(s2.getContainer()[j].getPosition()[0]);
				if (abs(pos1 - pos2) <= 2 * epsilon)
				{
					sum += sqrt(s1.getContainer()[i].getIntensity() * s2.getContainer()[j].getIntensity());
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
        if (abs(it1->getPosition()[0] - it2->getPosition()[0]) <= 2 * epsilon)
        {
          sum += sqrt(it1->getIntensity() * it2->getIntensity());
        }
      }
    }*/

    score = sum / (sqrt(sum1 * sum2));

    return score;
	
	}

}
