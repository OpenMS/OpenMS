// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors: Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <boost/numeric/conversion/cast.hpp>
#include <cmath>
#include <string>
#include <vector>


namespace OpenSwath
{
  ISpectrumAccess::~ISpectrumAccess()
  {
  }

  SpectrumSequence ISpectrumAccess::getMultipleSpectra(double RT, int nr_spectra_to_fetch)
  {
    std::vector<std::size_t> indices = getSpectraByRT(RT, 0.0);
    SpectrumSequence all_spectra;

    if (indices.empty() )
    {
      return all_spectra;
    }
    int closest_idx = boost::numeric_cast<int>(indices[0]);
    if (indices[0] != 0 &&
        std::fabs(getSpectrumMetaById(boost::numeric_cast<int>(indices[0]) - 1).RT - RT) <
        std::fabs(getSpectrumMetaById(boost::numeric_cast<int>(indices[0])).RT - RT))
    {
      closest_idx--;
    }

    all_spectra.push_back(getSpectrumById(closest_idx));

    int nrSpectra = (int) getNrSpectra();
    for (int i = 1; i <= nr_spectra_to_fetch / 2; i++) // cast to int is intended!
    {
      if (closest_idx - i >= 0)
      {
        all_spectra.push_back(getSpectrumById(closest_idx - i));
      }
      if (closest_idx + i < nrSpectra)
      {
        all_spectra.push_back(getSpectrumById(closest_idx + i));
      }
    }

    return all_spectra;
  }


  SpectrumSequence ISpectrumAccess::getMultipleSpectra(double RT, int nr_spectra_to_fetch, double drift_start, double drift_end)
  {
    std::vector<std::size_t> indices = getSpectraByRT(RT, 0.0);
    SpectrumSequence all_spectra;

    if (indices.empty() )
    {
      return all_spectra;
    }
    int closest_idx = boost::numeric_cast<int>(indices[0]);
    if (indices[0] != 0 &&
        std::fabs(getSpectrumMetaById(boost::numeric_cast<int>(indices[0]) - 1).RT - RT) <
        std::fabs(getSpectrumMetaById(boost::numeric_cast<int>(indices[0])).RT - RT))
    {
      closest_idx--;
    }

    all_spectra.push_back(getSpectrumById(closest_idx, drift_start, drift_end));

    int nrSpectra = (int) getNrSpectra();
    for (int i = 1; i <= nr_spectra_to_fetch / 2; i++) // cast to int is intended!
    {
      if (closest_idx - i >= 0)
      {
        all_spectra.push_back(getSpectrumById(closest_idx - i, drift_start, drift_end));
      }
      if (closest_idx + i < nrSpectra)
      {
        all_spectra.push_back(getSpectrumById(closest_idx + i, drift_start, drift_end));
      }
    }

    return all_spectra;
  }


  SpectrumPtr ISpectrumAccess::getSpectrumById(int id, double drift_start, double drift_end)
  {
    // first fetch the spectrum
    OpenSwath::SpectrumPtr spectrum = getSpectrumById(id);

    // then filter by drift
    return ISpectrumAccess::filterByDrift(spectrum, drift_start, drift_end);
  }
}
