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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>

using namespace std;

namespace OpenMS
{
  ParentPeakMower::ParentPeakMower() :
    DefaultParamHandler("ParentPeakMower")
  {
    defaults_.setValue("window_size", 2.0, "The size of the m/z window where the peaks are removed, +/- window_size.");
    defaults_.setValue("default_charge", 2, "If the precursor has no charge set, the default charge is assumed.");
    defaults_.setValue("clean_all_charge_states", 1, "Set to 1 if precursor ions of all possible charge states should be removed.", ListUtils::create<String>("advanced"));
    defaults_.setValue("consider_NH3_loss", 1, "Whether NH3 loss peaks from the precursor should be removed.");
    defaults_.setValue("consider_H2O_loss", 1, "Whether H2O loss peaks from the precursor should be removed.");
    defaults_.setValue("reduce_by_factor", 0, "Reduce the intensities of the precursor and related ions by a given factor (set 'set_to_zero' to 0).", ListUtils::create<String>("advanced"));
    defaults_.setValue("factor", 1000.0, "Factor which is used to reduce the intensities if 'reduce_by_factor' is selected.", ListUtils::create<String>("advanced"));
    defaults_.setValue("set_to_zero", 1, "Reduce the intensities of the precursor and related ions to zero.", ListUtils::create<String>("advanced"));
    defaultsToParam_();
  }

  ParentPeakMower::~ParentPeakMower()
  {
  }

  ParentPeakMower::ParentPeakMower(const ParentPeakMower & source) :
    DefaultParamHandler(source)
  {
  }

  ParentPeakMower & ParentPeakMower::operator=(const ParentPeakMower & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void ParentPeakMower::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    filterSpectrum(spectrum);
  }

  void ParentPeakMower::filterPeakMap(PeakMap & exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

}
