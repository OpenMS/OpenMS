// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  /**
   * @brief An implementation of the Spectrum Access interface using SQL files
   *
   * The interface takes an MzMLSqliteHandler object to access spectra and
   * chromatograms from a sqlite file (sqMass). Currently access to individual
   * spectra and chromatograms are not supported due to large overhead of
   * opening a DB connection and performing a single query.
   *
   * Instead, the users should use getAllSpectra which returns all available
   * spectra together.
   *
   * The interface allows to be constructed in a way as to only provide access
   * to a subset of spectra / chromatograms by supplying a set of indices which
   * are then used to provide a transparent interface to any consumer who will
   * get access to the described subset of spectra / chromatograms. This can be
   * useful to provide a specific interface to MS1 or MS2  spectra only or to
   * different DIA / SWATH-MS windows.
   *
   * Parallel access is supported through this interface as it is read-only and
   * sqlite3 supports multiple parallel read threads as long as they use a
   * different db connection.
   *
   * Sample usage:
   *
   *
   * @code
   *   // Obtain swath_map with boundaries first
   *   std::vector<int> indices = sql_mass_reader.readSpectraForWindow(swath_map);
   *   OpenMS::Internal::MzMLSqliteHandler handler(file);
   *   OpenSwath::SpectrumAccessPtr sptr(new OpenMS::SpectrumAccessSqMass(handler, indices));
   *   swath_maps[k].sptr = sptr;
   * @endcode
   *
   *
  */
  class OPENMS_DLLAPI SpectrumAccessSqMass :
    public OpenSwath::ISpectrumAccess
  {

public:
    typedef OpenMS::MSSpectrum MSSpectrumType;
    typedef OpenMS::MSChromatogram MSChromatogramType;

    /// Constructor
    SpectrumAccessSqMass(const OpenMS::Internal::MzMLSqliteHandler& handler);

    SpectrumAccessSqMass(const OpenMS::Internal::MzMLSqliteHandler& handler, const std::vector<int> & indices);

    SpectrumAccessSqMass(const SpectrumAccessSqMass& sp, const std::vector<int>& indices);

    /// Destructor
    ~SpectrumAccessSqMass() override;

    /// Copy constructor
    SpectrumAccessSqMass(const SpectrumAccessSqMass & rhs);

    /// Light clone operator (actual data will not get copied)
    boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;

    OpenSwath::SpectrumPtr getSpectrumById(int /* id */) override;

    OpenSwath::SpectrumMeta getSpectrumMetaById(int /* id */) const override;

    /// Load all spectra from the underlying sqMass file into memory
    void getAllSpectra(std::vector< OpenSwath::SpectrumPtr > & spectra, std::vector< OpenSwath::SpectrumMeta > & spectra_meta) const;

    std::vector<std::size_t> getSpectraByRT(double /* RT */, double /* deltaRT */) const override;

    size_t getNrSpectra() const override;

    OpenSwath::ChromatogramPtr getChromatogramById(int /* id */) override;

    size_t getNrChromatograms() const override;

    std::string getChromatogramNativeID(int /* id */) const override;

private:

    /// Access to underlying sqMass file
    OpenMS::Internal::MzMLSqliteHandler handler_;
    /// Optional subset of spectral indices
    std::vector<int> sidx_;
  };
} //end namespace OpenMS



