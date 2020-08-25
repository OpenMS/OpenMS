// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataAggregatingConsumer.h>

#include <OpenMS/FILTERING/TRANSFORMERS/LinearResamplerAlign.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>

namespace OpenMS
{

  MSDataAggregatingConsumer::~MSDataAggregatingConsumer()
  {
    // flush remaining spectra
    if (!s_list.empty())
    {
      MSSpectrum tmps = SpectrumAddition::addUpSpectra(s_list, -1, true);
      tmps.SpectrumSettings::operator=(s_list[0]); // copy over SpectrumSettings of first spectrum
      tmps.setName( s_list[0].getName() );
      tmps.setRT( s_list[0].getRT() );
      tmps.setDriftTime( s_list[0].getDriftTime() );
      tmps.setDriftTimeUnit( s_list[0].getDriftTimeUnit() );
      tmps.setMSLevel( s_list[0].getMSLevel() );

      next_consumer_->consumeSpectrum(tmps);
    }
  }

  void MSDataAggregatingConsumer::consumeSpectrum(SpectrumType & s)
  {
    // aggregate by RT
    double RT = s.getRT();

    if (rt_initialized_ && std::fabs(RT - previous_rt_) < 1e-5)
    {
      // need to aggregate spectrum
      s_list.push_back(s);
    }
    else
    {
      // consume the previous list
      if (rt_initialized_ && !s_list.empty())
      {
        MSSpectrum tmps = SpectrumAddition::addUpSpectra(s_list, -1, true);
        tmps.SpectrumSettings::operator=(s_list[0]); // copy over SpectrumSettings of first spectrum
        tmps.setName( s_list[0].getName() );
        tmps.setRT( s_list[0].getRT() );
        tmps.setDriftTime( s_list[0].getDriftTime() );
        tmps.setDriftTimeUnit( s_list[0].getDriftTimeUnit() );
        tmps.setMSLevel( s_list[0].getMSLevel() );

        next_consumer_->consumeSpectrum(tmps);
      }

      // start new spectrum list
      int expected_size = s_list.size();
      s_list.clear();
      s_list.reserve(expected_size);
      s_list.push_back(s);
    }

    previous_rt_ = RT;
    rt_initialized_ = true;
  }

  void MSDataAggregatingConsumer::consumeChromatogram(ChromatogramType & c)
  {
    // NOP
    next_consumer_->consumeChromatogram(c);
  }


} // namespace OpenMS
