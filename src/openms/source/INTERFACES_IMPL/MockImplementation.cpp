// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/INTERFACES/DataStructures.h>
#include <OpenMS/INTERFACES/ISpectrumAccess.h>

namespace OpenMS
{

/**
  @brief Mock implementations of the interfaces (empty ones)
*/
namespace Interfaces
{

  class OPENMS_DLLAPI MockISpectraReader :
    public ISpectraReader
  {
public:
    MockISpectraReader() {}
    ~MockISpectraReader() override {}
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
    MockIChromatogramsReader() {}
    ~MockIChromatogramsReader() override {}
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
    MockISpectraWriter() {}
    ~MockISpectraWriter() override {}
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
    MockIChromatogramsWriter() {}
    ~MockIChromatogramsWriter() override {}
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
}
