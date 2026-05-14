// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// forward declarations
struct sqlite3;

namespace OpenMS
{

  namespace Internal
  {

    /**
        @brief Read-only accessor for SWATH/DIA spectrum indices stored in an sqMass file.

        Exposes three lookups against a single sqMass file:
          - @ref readSwathWindows -- the SWATH window definitions
            (precursor isolation centre and lower/upper offsets),
          - @ref readMS1Spectra -- the spectrum indices of all MS1 scans,
          - @ref readSpectraForWindow -- the spectrum indices that belong
            to a given SWATH window, identified by its centre m/z.

        Internal helper class used by the OpenSWATH I/O stack; not part of
        the stable public API.

        @ingroup FileIO
    */
    class OPENMS_DLLAPI MzMLSqliteSwathHandler
    {

public:

      /**
          @brief Construct from the path of an sqMass file.

          The file is not opened by the constructor; each accessor opens its
          own connection on demand.

          @param[in] filename Path of the sqMass file to read.
      */
      MzMLSqliteSwathHandler(const String& filename) :
        filename_(filename)
      {}

      /**
          @brief Read the SWATH window boundaries from the file.

          Returns one @c SwathMap entry per distinct precursor isolation
          centre present at MS level 2; @c center, @c lower and @c upper
          are filled. Other @c SwathMap fields are left at their default
          values.

          @return SWATH window definitions in the file's natural order.

          @throws Exception::SqlOperationFailed if the sqMass file cannot
                                                be opened.
          @throws Exception::IllegalArgument    if preparing the SQL query
                                                fails.
      */
      std::vector<OpenSwath::SwathMap> readSwathWindows();

      /**
          @brief Read the indices of every MS1 spectrum.

          @return Spectrum indices for the MS1 scans in the file's natural
                  order.

          @throws Exception::SqlOperationFailed if the sqMass file cannot
                                                be opened.
          @throws Exception::IllegalArgument    if preparing the SQL query
                                                fails.
      */
      std::vector<int> readMS1Spectra();

      /**
          @brief Read the indices of every spectrum that belongs to a given
                 SWATH window.

          The window is identified by @p swath_map's @c center value;
          spectra whose precursor isolation centre lies within a small fixed
          tolerance (@c 0.01 m/z) of @c swath_map.center are returned.

          @note Only @c swath_map.center is consulted; the @c lower and
                @c upper fields are ignored.

          @param[in] swath_map SWATH window to look up; only its @c center
                               is used.
          @return Spectrum indices that fall into the window, in the file's
                  natural order.

          @throws Exception::SqlOperationFailed if the sqMass file cannot
                                                be opened.
          @throws Exception::IllegalArgument    if preparing the SQL query
                                                fails.
      */
      std::vector<int> readSpectraForWindow(const OpenSwath::SwathMap & swath_map);

protected:

      /// Path of the sqMass file passed in at construction time.
      String filename_;

      /// Next spectrum id used when appending spectra to the database; lets multiple append passes coexist on the same file.
      Int spec_id_;

      /// Next chromatogram id used when appending chromatograms to the database; lets multiple append passes coexist on the same file.
      Int chrom_id_;

    };


  }   // namespace Internal
} // namespace OpenMS

