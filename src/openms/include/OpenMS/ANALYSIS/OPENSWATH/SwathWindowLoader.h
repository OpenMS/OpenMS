// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest$
// $Authors: Hannes Roest$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <vector>

namespace OpenMS
{

  /**
    @brief Reads a SWATH-window definition file and annotates SwathMaps with it.

    The file format is tab- or space-delimited with two columns per row:
\verbatim
window_lower window_upper
400 425
425 450
...
\endverbatim

    The first line is treated as a header only if its first token cannot
    be parsed as a number; otherwise the first line is read as data. The
    reader does not require any specific header text and accepts files
    without a header row.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI SwathWindowLoader
  {

  public:

    /**
      @brief Overwrite the lower/upper precursor boundaries of @p swath_maps with the windows read from @p filename.

      Walks @p swath_maps in order, skipping entries whose @c ms1 flag is
      set, and replaces each remaining entry's @c lower and @c upper with
      the corresponding row of the SWATH-window file. The number of
      non-MS1 entries in @p swath_maps must equal the number of rows in
      @p filename.

      @note The non-MS1 entries in @p swath_maps and the rows of the file
            are matched by position, so both must be in the same order
            (usually low-to-high m/z). When @p do_sort is @c true,
            @p swath_maps is sorted by its @c upper field before matching.

      @param[in]     filename   Path of the SWATH-window file.
      @param[in,out] swath_maps SWATH maps to annotate; their @c lower /
                                @c upper fields are overwritten in place
                                (and the vector may be reordered when
                                @p do_sort is @c true).
      @param[in]     do_sort    If @c true, sort @p swath_maps by ascending
                                @c upper before matching against the file.
      @param[in]     force      If @c true, a file window that extends
                                beyond the corresponding map's data
                                boundaries is logged as an error and
                                applied anyway; if @c false the same
                                situation throws.

      @throws Exception::IllegalArgument if the number of non-MS1 entries
                                         in @p swath_maps does not equal
                                         the number of rows in
                                         @p filename, or if a file window
                                         extends beyond the corresponding
                                         map's @c lower / @c upper and
                                         @p force is @c false.
      @throws Exception::InvalidValue    propagated from @ref readSwathWindows
                                         when a row has @c lower >= @c upper.
    */
    static void annotateSwathMapsFromFile(const std::string& filename,
                                          std::vector<OpenSwath::SwathMap>& swath_maps,
                                          bool do_sort,
                                          bool force);

    /**
      @brief Read a SWATH-window definition file into two parallel vectors.

      The file format is two whitespace-separated columns
      (@c window_lower @c window_upper) per row:
\verbatim
window_lower window_upper
400 425
425 450
...
\endverbatim

      A header row is detected automatically by trying to parse the first
      token of the first line as a number; if parsing fails, the line is
      treated as a header and skipped, otherwise the first line is read as
      data. No specific header text is required.

      @param[in]  filename         Path of the SWATH-window file. A
                                   missing or unreadable file is not
                                   reported as an exception; the output
                                   vectors are left empty.
      @param[out] swath_prec_lower Lower boundary of each window, in file
                                   order.
      @param[out] swath_prec_upper Upper boundary of each window, in file
                                   order; same length as
                                   @p swath_prec_lower.

      @throws Exception::InvalidValue if a row has @c lower >= @c upper.
    */
    static void readSwathWindows(const std::string& filename,
                                 std::vector<double>& swath_prec_lower,
                                 std::vector<double>& swath_prec_upper);
  };
}

