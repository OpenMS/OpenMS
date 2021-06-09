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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/IONMOBILITY/MSRunIMSplitter.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <map>

namespace OpenMS
{
  std::vector<PeakMap> MSRunIMSplitter::splitByFAIMSCV(PeakMap&& exp)
  {
    std::vector<PeakMap> split_peakmap;

    // TODO test with any random PeakMap without FAIMS data.
    // What breaks, how should it break?
    std::set<double> CVs = FAIMSHelper::getCompensationVoltages(exp);

    if (CVs.empty())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not FAIMS data!");
    }

    // create map to easily turn a CV value into a PeakMap index
    std::map<double, size_t> cv2index;
    size_t counter(0);
    for (double cv : CVs)
    {
      cv2index[cv] = counter;
      counter++;
    }

    // make as many PeakMaps as there are different CVs and fill their Meta Data
    split_peakmap.resize(CVs.size());
    for (auto it = split_peakmap.begin(); it != split_peakmap.end(); ++it)
    {
      it->getExperimentalSettings() = exp.getExperimentalSettings();
    }

    // fill up the PeakMaps by moving spectra from the input PeakMap
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      split_peakmap[cv2index[it->getDriftTime()]].addSpectrum(std::move(*it));
    }

    return split_peakmap;
  }

}  //end namespace OpenMS
