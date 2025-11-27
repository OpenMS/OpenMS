// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/INTERFACES/DataStructures.h>
#include <OpenMS/INTERFACES/ISpectrumAccess.h>

namespace OpenMS::Interfaces
{

  class OPENMS_DLLAPI MockISpectraReader :
    public ISpectraReader
  {
public:
    MockISpectraReader() = default;
    ~MockISpectraReader() override = default;
    /// Return a pointer to a spectrum at the given id
    SpectrumPtr getSpectrumById(int /* id */) const override
    {
      SpectrumPtr spectrum(new Spectrum);
      return spectrum;
    }
    /// Return a pointer to a spectrum at the given string id
    SpectrumPtr getSpectrumById(const std::string& /* id */) const override
    {
      SpectrumPtr spectrum(new Spectrum);
      return spectrum;
    }
    /// Return a vector of ids of spectra that are within RT +/- deltaRT
    std::vector<std::size_t> getSpectraByRT(double /* RT */, double /* deltaRT */) const override
    {
      return std::vector<std::size_t>();
    }
    /// Returns the number of spectra available
    size_t getNrSpectra() const override
    {
      return 0;
    }
    /// Returns the meta information for a spectrum
    SpectrumMetaPtr getSpectrumMetaById(int /* id */) const override
    {
      SpectrumMetaPtr spectrum_meta(new SpectrumMeta);
      return spectrum_meta;
    }
  };
  // create an instance of the mock object to test
  MockISpectraReader test_mock_spectra_reader;

  class OPENMS_DLLAPI MockIChromatogramsReader :
    public IChromatogramsReader
  {
public:
    MockIChromatogramsReader() = default;
    ~MockIChromatogramsReader() override = default;
    /// Return a pointer to a chromatogram at the given id
    ChromatogramPtr getChromatogramById(int /* id */) const override
    {
      ChromatogramPtr chromatogram(new Chromatogram);
      return chromatogram;
    }
    /// Return a pointer to a chromatogram at the given string id
    ChromatogramPtr getChromatogramById(const std::string& /* id */) const override
    {
      ChromatogramPtr chromatogram(new Chromatogram);
      return chromatogram;
    }
    /// Return a vector of ids of spectra that are within RT +/- deltaRT
    std::vector<std::size_t> getChromatogramByPrecursorMZ(double /* mz */, double /* deltaMZ */) const override
    {
      return std::vector<std::size_t>();
    }
    /// Returns the number of spectra available
    size_t getNrChromatograms() const override
    {
      return 0;
    }
    /// Returns the meta information for a chromatogram
    ChromatogramMetaPtr getChromatogramMetaById(int /* id */) const override
    {
      ChromatogramMetaPtr chromatogram_meta(new ChromatogramMeta);
      return chromatogram_meta;
    }
  };
  // create an instance of the mock object to test
  MockIChromatogramsReader test_mock_chromatograms_reader;

  class OPENMS_DLLAPI MockISpectraWriter :
    public ISpectraWriter
  {
public:
    MockISpectraWriter() = default;
    ~MockISpectraWriter() override = default;
    /// Append a spectrum to the end
    void appendSpectrum(SpectrumPtr /* spectrum */, bool /* write_through*/) override
    {
      // do nothing
    }
    /// write all cached data to disk
    void flush() override
    {
      // do nothing
    }
  };
  // create an instance of the mock object to test
  MockISpectraWriter test_mock_specrta_writer;

  class OPENMS_DLLAPI MockIChromatogramsWriter :
    public IChromatogramsWriter
  {
public:
    MockIChromatogramsWriter() = default;
    ~MockIChromatogramsWriter() override = default;
    /// Append a chromatogram to the end
    void appendChromatogram(ChromatogramPtr /* chromatogram */, bool /* write_through */) override
    {
      // do nothing
    }
    /// write all cached data to disk
    void flush() override
    {
      // do nothing
    }
  };
  // create an instance of the mock object to test
  MockIChromatogramsWriter test_mock_chromatograms_writer;

}
