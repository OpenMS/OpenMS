// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <boost/shared_ptr.hpp>
#include <string>
#include <vector>

namespace OpenSwath
{

  /**
    @brief The interface of a mass spectrometry experiment.
  */
  class OPENSWATHALGO_DLLAPI ISpectrumAccess
  {
public:
    /// Destructor
    virtual ~ISpectrumAccess();

    /**
      @brief Light clone operator to produce a copy for concurrent read access.

      This function guarantees to produce a copy of the underlying object that
      provides thread-safe concurrent read access to the underlying data. It
      should be implemented with minimal copy-overhead to make this operation
      as fast as possible.

      To use this function, each thread should call this function to produce an
      individual copy on which it can operate.

    */
    virtual boost::shared_ptr<ISpectrumAccess> lightClone() const = 0;

    /// Return a pointer to a spectrum at the given id
    virtual SpectrumPtr getSpectrumById(int id) = 0;
    /// Return a vector of ids of spectra that are within RT +/- deltaRT
    virtual std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const = 0;
    /// Returns the number of spectra available
    virtual size_t getNrSpectra() const = 0;
    /// Returns the meta information for a spectrum
    virtual SpectrumMeta getSpectrumMetaById(int id) const = 0;

    /// Return a pointer to a chromatogram at the given id
    virtual ChromatogramPtr getChromatogramById(int id) = 0;
    /// Returns the number of chromatograms available
    virtual std::size_t getNrChromatograms() const = 0;
    /// Returns the native id of the chromatogram at the given id
    virtual std::string getChromatogramNativeID(int id) const = 0;
  };

  typedef boost::shared_ptr<ISpectrumAccess> SpectrumAccessPtr;
}

