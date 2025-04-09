// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
