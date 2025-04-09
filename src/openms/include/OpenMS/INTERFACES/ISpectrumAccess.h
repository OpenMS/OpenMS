// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/DataStructures.h>

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

namespace OpenMS
{
namespace Interfaces
{

  /**
    @brief The interface of read-access to a list of spectra
  */
  class OPENMS_DLLAPI ISpectraReader
  {
public:
    virtual ~ISpectraReader() {}
    /// Return a pointer to a spectrum at the given id
    virtual SpectrumPtr getSpectrumById(int id) const = 0;
    /// Return a pointer to a spectrum at the given string id
    virtual SpectrumPtr getSpectrumById(const std::string& id) const = 0;
    /// Return a vector of ids of spectra that are within RT +/- deltaRT
    virtual std::vector<std::size_t> getSpectraByRT(double RT, double deltaRT) const = 0;
    /// Returns the number of spectra available
    virtual size_t getNrSpectra() const = 0;
    /// Returns the meta information for a spectrum
    virtual SpectrumMetaPtr getSpectrumMetaById(int id) const = 0;

    /*
     * Do we need an Iterator here?
     * We would have to provide our own iterator wrapper class because we don't
     * know whether all the spectra are loaded at any given timepoint
    typedef SpectrumPtr const ConstSpectraIterator;
    ConstSpectraIterator beginSpectra() const;
    ConstSpectraIterator endSpectra() const;
    */
  };
  typedef boost::shared_ptr<ISpectraReader> SpectraReaderPtr;


  /**
    @brief The interface of read-access to a list of chromatograms
  */
  class OPENMS_DLLAPI IChromatogramsReader
  {
public:
    virtual ~IChromatogramsReader() {}
    /// Return a pointer to a chromatogram at the given id
    virtual ChromatogramPtr getChromatogramById(int id) const = 0;
    /// Return a pointer to a chromatogram at the given string id
    virtual ChromatogramPtr getChromatogramById(const std::string& id) const = 0;
    /// Return a vector of ids of chromatograms that are within mz +/- deltaMz
    virtual std::vector<std::size_t> getChromatogramByPrecursorMZ(double mz, double deltaMZ) const = 0;
    /// Returns the number of chromatograms available
    virtual std::size_t getNrChromatograms() const = 0;
    /// Returns the meta information for a chromatogram
    virtual ChromatogramMetaPtr getChromatogramMetaById(int id) const = 0;

    /*
     * Do we need an Iterator here?
     * We would have to provide our own iterator wrapper class because we don't
     * know whether all the chromatograms are loaded at any given timepoint
    ConstChromatogramIterator beginChromatograms() const;
    ConstChromatogramIterator endChromatograms() const;
    */
  };
  typedef boost::shared_ptr<IChromatogramsReader> ChromatogramsReaderPtr;


  class OPENMS_DLLAPI ISpectraWriter
  {
public:
    virtual ~ISpectraWriter() {}
    /// Append a spectrum to the end
    virtual void appendSpectrum(SpectrumPtr spectrum, bool write_through=false) = 0;
    /// write all cached data to disk
    virtual void flush() = 0;
  };
  typedef boost::shared_ptr<ISpectraWriter> SpectraWriterPtr;


  class OPENMS_DLLAPI IChromatogramsWriter
  {
public:
    virtual ~IChromatogramsWriter() {}
    /// Append a chromatogram to the end
    virtual void appendChromatogram(ChromatogramPtr chromatogram, bool write_through=false) = 0;
    /// write all cached data to disk
    virtual void flush() = 0;
  };
  typedef boost::shared_ptr<IChromatogramsWriter> ChromatogramsWriterPtr;

} //end namespace Interfaces
} //end namespace OpenMS
