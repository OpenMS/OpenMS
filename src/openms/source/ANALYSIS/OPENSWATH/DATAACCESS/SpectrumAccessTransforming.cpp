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

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>

namespace OpenMS
{

  SpectrumAccessTransforming::SpectrumAccessTransforming(OpenSwath::SpectrumAccessPtr sptr) :
    sptr_(sptr)
  {}

  SpectrumAccessTransforming::~SpectrumAccessTransforming() {}

  size_t SpectrumAccessTransforming::getNrChromatograms() const
  {
    return sptr_->getNrChromatograms();
  }


  OpenSwath::SpectrumPtr SpectrumAccessTransforming::getSpectrumById(int id)
  {
    return sptr_->getSpectrumById(id);
  }

  OpenSwath::SpectrumMeta SpectrumAccessTransforming::getSpectrumMetaById(int id) const
  {
    return sptr_->getSpectrumMetaById(id);
  }

  std::vector<std::size_t> SpectrumAccessTransforming::getSpectraByRT(double RT, double deltaRT) const
  {
    return sptr_->getSpectraByRT(RT, deltaRT);
  }

  size_t SpectrumAccessTransforming::getNrSpectra() const
  {
    return sptr_->getNrSpectra();
  }

  OpenSwath::ChromatogramPtr SpectrumAccessTransforming::getChromatogramById(int id)
  {
    return sptr_->getChromatogramById(id);
  }

  std::string SpectrumAccessTransforming::getChromatogramNativeID(int id) const
  {
    return sptr_->getChromatogramNativeID(id);
  }

}

