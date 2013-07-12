// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
  ZhangSimilarityScore::ZhangSimilarityScore() :
    PeakSpectrumCompareFunctor()
  {
    setName(ZhangSimilarityScore::getProductName());
    defaults_.setValue("tolerance", 0.2, "defines the absolute (in Da) or relative (in ppm) tolerance");
    defaults_.setValue("is_relative_tolerance", "false", "If set to true, the tolerance is interpreted as relative");
    defaults_.setValidStrings("is_relative_tolerance", StringList::create("true,false"));
    defaults_.setValue("use_linear_factor", "false", "if true, the intensities are weighted with the relative m/z difference");
    defaults_.setValidStrings("use_linear_factor", StringList::create("true,false"));
    defaults_.setValue("use_gaussian_factor", "false", "if true, the intensities are weighted with the relative m/z difference using a gaussian");
    defaults_.setValidStrings("use_gaussian_factor", StringList::create("true,false"));
    defaultsToParam_();
  }

  ZhangSimilarityScore::ZhangSimilarityScore(const ZhangSimilarityScore & source) :
    PeakSpectrumCompareFunctor(source)
  {
  }

  ZhangSimilarityScore::~ZhangSimilarityScore()
  {
  }

  ZhangSimilarityScore & ZhangSimilarityScore::operator=(const ZhangSimilarityScore & source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  double ZhangSimilarityScore::operator()(const PeakSpectrum & spec) const
  {
    return operator()(spec, spec);
  }

  double ZhangSimilarityScore::operator()(const PeakSpectrum & s1, const PeakSpectrum & s2) const
  {
    const double tolerance = (double)param_.getValue("tolerance");
    bool use_linear_factor = param_.getValue("use_linear_factor").toBool();
    bool use_gaussian_factor = param_.getValue("use_gaussian_factor").toBool();
    double score(0), sum(0), sum1(0), sum2(0) /*, squared_sum1(0), squared_sum2(0)*/;

    // TODO remove parameter 
    if (param_.getValue("is_relative_tolerance").toBool() )
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

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
