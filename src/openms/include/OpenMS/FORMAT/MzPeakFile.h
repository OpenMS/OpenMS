// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <string>

namespace OpenMS
{
  /**
    @brief mzPeak file reader and writer (HUPO-PSI mzPeak: Parquet tables in a ZIP archive or a directory).

    Reads an mzPeak run into an MSExperiment or streams it into an Interfaces::IMSDataConsumer (as MzMLFile::transform()
    does), and writes an MSExperiment as an mzPeak archive. The format is implemented by the vendored mzPeak C++ library
    (src/openms/thirdparty/mzpeak, see its README.OpenMS.md for the library commit and the format version), which is
    available when OpenMS is built with OPENMS_ENABLE_MZPEAK (the WITH_MZPEAK define); without it, load(), transform()
    and store() throw Exception::NotImplemented.

    An mzPeak run is a ZIP archive ('.mzpeak' file) with an index ('mzpeak_index.json') and Parquet members, or a
    directory with the same members; isMzPeak() tells the two from other files and directories.

    Reading: spectra (m/z, intensity, per-peak ion mobility arrays), their retention time (in seconds, converted from the
    minutes the format stores), MS level, native ID, polarity, centroid/profile representation, scan windows, ion
    mobility (drift time with its unit), precursors (isolation window, selected ion m/z, charge, intensity, ion mobility,
    activation methods and energy) and chromatograms; run-level metadata (run start time, source files, instrument
    model). The scan-level PeakFileOptions (MS levels, RT range, precursor m/z range) are applied before the peaks are
    decoded; the peak-level options (m/z, intensity, ion mobility ranges) after. The spectra are decoded in parallel with
    the OpenMP threads of the calling thread (one reader per thread, as the library's thread-safety model requires).

    Writing: m/z, intensity, representation, MS level, retention time, polarity, native ID, precursors (isolation window,
    selected ion, charge, intensity) and the spectrum's ion mobility; run-level metadata (run start time, source files
    with their checksum, file format and native ID format). Not written (yet): per-peak ion mobility arrays,
    activation methods, scan windows, chromatograms and instrument metadata.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MzPeakFile :
    public ProgressLogger
  {
public:
    MzPeakFile();
    ~MzPeakFile() override;

    /// Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();
    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;
    /// set options for loading/storing
    void setOptions(const PeakFileOptions&);

    /**
      @brief Loads a run from an mzPeak archive or directory into @p exp

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
      @exception Exception::NotImplemented is thrown if OpenMS was built without mzPeak support
    */
    void load(const std::string& filename, PeakMap& exp);

    /**
      @brief Streams the spectra (and chromatograms) of an mzPeak archive or directory into @p consumer

      The consumer gets the experimental settings and the expected sizes first, then the spectra in file order,
      decoded in batches with the OpenMP threads of the calling thread, then the chromatograms.

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
      @exception Exception::NotImplemented is thrown if OpenMS was built without mzPeak support
    */
    void transform(const std::string& filename, Interfaces::IMSDataConsumer* consumer);

    /**
      @brief Stores a run as an mzPeak archive (ZIP) at @p filename

      @exception Exception::UnableToCreateFile is thrown if the file could not be created
      @exception Exception::NotImplemented is thrown if OpenMS was built without mzPeak support
    */
    void store(const std::string& filename, const PeakMap& exp) const;

    /**
      @brief Whether @p path looks like an mzPeak run: a directory holding 'mzpeak_index.json', or a file that is a ZIP
      archive (whose members are checked when it is opened). False for anything that does not exist.
    */
    static bool isMzPeak(const std::string& path);

    /// Whether OpenMS was built with mzPeak support (OPENMS_ENABLE_MZPEAK)
    static bool isSupported();

private:
    PeakFileOptions options_;
  };

} // namespace OpenMS
