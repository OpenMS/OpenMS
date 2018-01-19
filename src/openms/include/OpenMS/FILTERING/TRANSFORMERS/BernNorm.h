// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_BERNNORM_H
#define OPENMS_FILTERING_TRANSFORMERS_BERNNORM_H

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <map>

namespace OpenMS
{
  /**
    @brief BernNorm scales the peaks by ranking them and then scaling them according to rank.

    For exact formula look in Bioinformatics, Aug 2004; 20: i49 - i54

        @improvement read paper and try to confirm implementation (Andreas)

        @htmlinclude OpenMS_BernNorm.parameters

        @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI BernNorm :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    //@{
    /// default constructor
    BernNorm();

    /// copy constructor
    BernNorm(const BernNorm & source);

    /// destructor
    ~BernNorm() override;
    //@}

    // @name Operators
    // @{
    /// assignment operator
    BernNorm & operator=(const BernNorm & source);
    //@}

    // @name Accessors
    // @{

    ///
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      typedef typename SpectrumType::Iterator Iterator;
      typedef typename SpectrumType::ConstIterator ConstIterator;

      c1_ = (double)param_.getValue("C1");
      c2_ = (double)param_.getValue("C2");
      th_ = (double)param_.getValue("threshold");

      spectrum.sortByPosition();

      // find highest peak and ranking
      double maxint = 0;
      std::map<double, Size> peakranks;
      for (ConstIterator it = spectrum.begin(); it != spectrum.end(); ++it)
      {
        peakranks[it->getIntensity()] = 0;
        if (it->getIntensity() > maxint)
        {
          maxint = it->getIntensity();
        }
      }
      UInt rank = 0;
      for (std::map<double, Size>::reverse_iterator mit = peakranks.rbegin(); mit != peakranks.rend(); ++mit)
      {
        mit->second = ++rank;
      }

      // find maxmz i.e. significant (> threshold * maxpeak) peak with highest m/z
      double maxmz = 0;
      for (SignedSize i = spectrum.size() - 1; i >= 0; --i)
      {
        if (spectrum[i].getIntensity() > maxint * th_)
        {
          maxmz = spectrum[i].getMZ();
          break;
        }
      }

      // rank
      for (Iterator it = spectrum.begin(); it != spectrum.end(); )
      {
        double newint = c1_ - (c2_ / maxmz) * peakranks[it->getIntensity()];
        if (newint < 0)
        {
          it = spectrum.erase(it);
        }
        else
        {
          it->setIntensity(newint);
          ++it;
        }
      }
      return;
    }

    void filterPeakSpectrum(PeakSpectrum & spectrum);

    void filterPeakMap(PeakMap & exp);
    //TODO reimplement DefaultParamHandler::updateMembers_()

private:
    double c1_;
    double c2_;
    double th_;

    // @}

  };

} // namespace OpenMS

#endif //OPENMS_FILTERING_TRANSFORMERS_BERNNORM_H
