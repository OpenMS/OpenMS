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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>

#include <OpenMS/CONCEPT/Macros.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

namespace OpenMS
{
  MSDataCachedConsumer::MSDataCachedConsumer(const String& filename, bool clearData) :
    ofs_(filename.c_str(), std::ios::binary),
    clearData_(clearData),
    spectra_written_(0),
    chromatograms_written_(0)
  {
    int file_identifier = CACHED_MZML_FILE_IDENTIFIER;
    ofs_.write((char*)&file_identifier, sizeof(file_identifier));
  }

  MSDataCachedConsumer::~MSDataCachedConsumer()
  {
    // Write size of file (to the end of the file)
    ofs_.write((char*)&spectra_written_, sizeof(spectra_written_));
    ofs_.write((char*)&chromatograms_written_, sizeof(chromatograms_written_));

    // Close file stream: close() _should_ call flush() but it might not in
    // all cases. To be sure call flush() first.
    ofs_.flush();
    ofs_.close();
  }

  void MSDataCachedConsumer::consumeSpectrum(SpectrumType & s)
  {
    if (chromatograms_written_ > 0)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Cannot write spectra after writing chromatograms.");
    }
    writeSpectrum_(s, ofs_);
    spectra_written_++;

    // Clear all spectral data including all float/int data arrays (but not string arrays)
    if (clearData_)
    {
      s.clear(false);
      s.setFloatDataArrays({});
      s.setIntegerDataArrays({});
    }
    OPENMS_POSTCONDITION( (!clearData_ || s.empty() ), "clearData implies spectrum is empty")
    OPENMS_POSTCONDITION( (!clearData_ || s.getFloatDataArrays().empty() ), "clearData implies spectrum is empty")
    OPENMS_POSTCONDITION( (!clearData_ || s.getIntegerDataArrays().empty() ), "clearData implies spectrum is empty")
  }

  void MSDataCachedConsumer::consumeChromatogram(ChromatogramType & c)
  {
    writeChromatogram_(c, ofs_);
    chromatograms_written_++;

    // Clear all chromatogram data including all float/int data arrays (but not string arrays)
    if (clearData_)
    {
      c.clear(false);
      c.setFloatDataArrays({});
      c.setIntegerDataArrays({});
    }
    OPENMS_POSTCONDITION( (!clearData_ || c.empty() ), "clearData implies chromatogram is empty")
    OPENMS_POSTCONDITION( (!clearData_ || c.getFloatDataArrays().empty() ), "clearData implies chromatogram is empty")
    OPENMS_POSTCONDITION( (!clearData_ || c.getIntegerDataArrays().empty() ), "clearData implies chromatogram is empty")
  }

} // namespace OpenMS
