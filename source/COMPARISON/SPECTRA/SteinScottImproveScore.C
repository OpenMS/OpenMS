// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Erhan Kenar $
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/SteinScottImproveScore.h>


namespace OpenMS
{
  /// default constructor
  SteinScottImproveScore::SteinScottImproveScore() :
    PeakSpectrumCompareFunctor()
  {
    setName(SteinScottImproveScore::getProductName());
    defaults_.setValue("tolerance", 0.2, "defines the absolut error of the mass spectrometer");
    defaults_.setValue("threshold", 0.2, "if the calculated score is smaller than the threshold, a zero is given back");
    defaultsToParam_();
  }

  /// copy constructor
  SteinScottImproveScore::SteinScottImproveScore(const SteinScottImproveScore & source) :
    PeakSpectrumCompareFunctor(source)
  {
  }

  /// destructor
  SteinScottImproveScore::~SteinScottImproveScore()
  {
  }

  /// assignment operator
  SteinScottImproveScore & SteinScottImproveScore::operator=(const SteinScottImproveScore & source)
  {
    if (this != &source)
    {
      PeakSpectrumCompareFunctor::operator=(source);
    }
    return *this;
  }

  /**
  @brief Similarity pairwise score itself

  This function return the similarity score of itself based on SteinScott.

  @param spec  const PeakSpectrum Spectrum 1
@see SteinScottImproveScore()
  */
  double SteinScottImproveScore::operator()(const PeakSpectrum & spec) const
  {
    return operator()(spec, spec);
  }

  /**
  @brief Similarity pairwise score

  This function return the similarity score of two Spectrums based on SteinScott.

      @param s1  const PeakSpectrum Spectrum 1
  @param s2  const PeakSpectrum Spectrum 2
  @see SteinScottImproveScore()
  */
  double SteinScottImproveScore::operator()(const PeakSpectrum & s1, const PeakSpectrum & s2) const
  {
    const double epsilon = (double)param_.getValue("tolerance");
    const double constant = epsilon / 10000;

    //const double c(0.0004);
    double score(0), sum(0), sum1(0), sum2(0), sum3(0), sum4(0);
    /* std::cout << s1 << std::endl;
    std::cout << std::endl;
    std::cout << s2 << std::endl;*/

    for (PeakSpectrum::ConstIterator it1 = s1.begin(); it1 != s1.end(); ++it1)
    {
      double temp = it1->getIntensity();
      sum1 += temp * temp;
      sum3 += temp;
    }

    for (PeakSpectrum::ConstIterator it1 = s2.begin(); it1 != s2.end(); ++it1)
    {
      double temp = it1->getIntensity();
      sum2 += temp * temp;
      sum4 += temp;
    }
    double z = constant * (sum3 * sum4);
    Size j_left(0);
    for (Size i = 0; i != s1.size(); ++i)
    {
      for (Size j = j_left; j != s2.size(); ++j)
      {
        double pos1(s1[i].getMZ()), pos2(s2[j].getMZ());
        if (std::abs(pos1 - pos2) <= 2 * epsilon)
        {
          sum += s1[i].getIntensity() * s2[j].getIntensity();
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
    //std::cout<< sum << " Sum " << z << " z " << std::endl;
    score = (sum - z) / (std::sqrt((sum1 * sum2)));
    // std::cout<<score<< " score" << std::endl;
    if (score < (Real)param_.getValue("threshold"))
      score = 0;
    return score;
  }

}
