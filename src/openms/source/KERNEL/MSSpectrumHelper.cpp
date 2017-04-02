// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSSpectrumHelper.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <algorithm>

namespace OpenMS
{
  MSSpectrum<RichPeak1D> MSSpectrumHelper::clone(const MSSpectrum<Peak1D>& p)
  {
    MSSpectrum<RichPeak1D> rp;
    rp.retention_time_ = p.retention_time_;
    rp.ms_level_ = p.ms_level_;
    rp.name_ = p.name_;
    rp.float_data_arrays_.clear();

    //ugly stuff from MSSpectrum
    for (MSSpectrum<Peak1D>::FloatDataArrays::const_iterator fit = p.float_data_arrays_.begin(); fit!= p.float_data_arrays_.end(); ++fit)
    {
      MSSpectrum<RichPeak1D>::FloatDataArray f;
      f.MetaInfoInterface::operator=(*fit);
      std::copy(fit->begin(),fit->end(),f.begin());
      rp.float_data_arrays_.push_back(f);
    }
    for (MSSpectrum<Peak1D>::StringDataArrays::const_iterator fit = p.string_data_arrays_.begin(); fit!= p.string_data_arrays_.end(); ++fit)
    {
      MSSpectrum<RichPeak1D>::StringDataArray f;
      f.MetaInfoInterface::operator=(*fit);
      std::copy(fit->begin(),fit->end(),f.begin());
      rp.string_data_arrays_.push_back(f);
    }
    for (MSSpectrum<Peak1D>::IntegerDataArrays::const_iterator fit = p.integer_data_arrays_.begin(); fit!= p.integer_data_arrays_.end(); ++fit)
    {
      MSSpectrum<RichPeak1D>::IntegerDataArray f;
      f.MetaInfoInterface::operator=(*fit);
      std::copy(fit->begin(),fit->end(),f.begin());
      rp.integer_data_arrays_.push_back(f);
    }

    //grab SpectrumSettings
    rp.SpectrumSettings::operator=(p);
    //grab MetaInfoInterface
    rp.MetaInfoInterface::operator=(p);

    //copy cast peaks to riches
    std::copy(p.begin(), p.end(), std::back_inserter(rp));
    return rp;
  }


} // namespace OpenMS
