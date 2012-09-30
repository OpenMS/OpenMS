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
// $Maintainer: Hannes Roest ,  Witold Wolski$
// $Authors:  Hannes Roest , Witold Wolski$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSCached.h>

namespace OpenMS{

  OpenSwath::SpectrumAccessPtr SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(OpenMS::MSExperiment<OpenMS::Peak1D> & exp)
  {
    bool is_cached = false;
    for (std::size_t i = 0; i < exp.size(); ++i)
      {
        for (std::size_t j = 0; j < exp[i].getDataProcessing().size(); j++)
          {
            const OpenMS::DataProcessing & dp = exp[i].getDataProcessing()[j];
            if (dp.metaValueExists("cached_data"))
              {
                is_cached = true;
              }
          }
      }
    for (std::size_t i = 0; i < exp.getChromatograms().size(); ++i)
      {
        for (std::size_t j = 0; j < exp.getChromatograms()[i].getDataProcessing().size(); j++)
          {
            const OpenMS::DataProcessing & dp = exp.getChromatograms()[i].getDataProcessing()[j];
            if (dp.metaValueExists("cached_data"))
              {
                is_cached = true;
              }
          }
      }
    if (is_cached)
      {
        OpenSwath::SpectrumAccessPtr experiment(new OpenMS::SpectrumAccessOpenMSCached(exp.getLoadedFilePath()));
        return experiment;
      }
    else
      {
        OpenSwath::SpectrumAccessPtr experiment(new OpenMS::SpectrumAccessOpenMS(exp));
        return experiment;
      }
  }

}//end Namespace
