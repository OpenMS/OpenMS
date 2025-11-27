// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
