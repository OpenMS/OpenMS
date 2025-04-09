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

#include <OpenMS/FORMAT/CachedMzML.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <fstream>

namespace OpenMS
{

  /**
    @brief An implementation of the Spectrum Access interface using on-disk caching

    This class implements the OpenSWATH Spectrum Access interface
    (ISpectrumAccess) using the CachedmzML class which is able to read and
    write a cached mzML file.

    @note This implementation is @a not thread-safe since it keeps internally a
    single file access pointer which it moves when accessing a specific
    data item. The caller is responsible to ensure that access is performed
    atomically.

  */
  class OPENMS_DLLAPI SpectrumAccessOpenMSCached :
    public OpenSwath::ISpectrumAccess,
    public OpenMS::CachedmzML
  {

public:
    typedef OpenMS::PeakMap MSExperimentType;
    typedef OpenMS::MSSpectrum MSSpectrumType;

    /**
      @brief Constructor, opens the file stream

      @param filename The filename of the .mzML file (it is assumed a second
      file .mzML.cached exists).

      @throws Exception::FileNotFound is thrown if the file is not found
      @throws Exception::ParseError is thrown if the file cannot be parsed
    */
    explicit SpectrumAccessOpenMSCached(const String& filename);

    /**
      @brief Destructor
    */
    ~SpectrumAccessOpenMSCached() override;

    /// Copy constructor
    SpectrumAccessOpenMSCached(const SpectrumAccessOpenMSCached & rhs);

    /// Light clone operator (actual data will not get copied)
    boost::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override;

    OpenSwath::SpectrumPtr getSpectrumById(int id) override;

    OpenSwath::SpectrumMeta getSpectrumMetaById(int id) const override;

    std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const override;

    size_t getNrSpectra() const override;

    SpectrumSettings getSpectraMetaInfo(int id) const;

    OpenSwath::ChromatogramPtr getChromatogramById(int id) override;

    size_t getNrChromatograms() const override;

    ChromatogramSettings getChromatogramMetaInfo(int id) const;

    std::string getChromatogramNativeID(int id) const override;
  };

} //end namespace

