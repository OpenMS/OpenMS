// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:	Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Lukas Zimmermann $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/TypeAliases.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <map>
#include <set>

namespace OpenMS
{
  class ExperimentalDesign;
  class TextFile; 
  /**
  @brief Load an experimental design from a TSV file

  The format -- both the one-table and the two-table variant -- its columns, the rules a design
  has to satisfy, and worked examples are documented on OpenMS::ExperimentalDesign. In short:

  - TAB-separated; cells are whitespace trimmed and lines starting with a hash character
    (a comment) are ignored.
  - The variant is auto-detected. A file is read as two-table as soon as some line carries a
    @c Sample column header but no @c Fraction_Group column header (that line is the header of
    the sample section, which a blank line separates from the MS file section above it).
    Otherwise it is read as one-table.
  - Mandatory columns of the MS file section are @c Fraction_Group, @c Fraction and
    @c Spectra_Filepath; @c Label (default 1) and @c Sample are optional. In the one-table
    format any further column is read as sample metadata, whereas the file section of a
    two-table design rejects unknown columns.
  - A relative @c Spectra_Filepath is resolved against the directory of the design file first,
    then against the current working directory.

  @ingroup Format
  */

  class OPENMS_DLLAPI ExperimentalDesignFile
  {
    public:

    /**
      @brief Loads an experimental design from a tabular separated file

      @param[in] tsv_file Path of the design file
      @param[in] require_spectra_files If true, every @c Spectra_Filepath must resolve to an
                 existing file; otherwise unresolvable paths are kept as written
      @throw Exception::ParseError on a malformed table, a missing mandatory column, or -- with
             @p require_spectra_files -- a spectra file that does not exist
      @throw Exception::InvalidValue if the fraction groups are not consecutive starting at 1, or
             a (path, label) combination is ambiguous
      @throw Exception::MissingInformation if a (fraction group, fraction, label) combination
             repeats, or a label-free (fraction group, label) maps to more than one sample
    */
    static ExperimentalDesign load(const std::string& tsv_file, bool require_spectra_files);

    /// Loads an experimental design from an already loaded or generated, tabular file.
    /// @p filename is only used to name the source in error messages.
    /// @see load(const std::string&, bool) for the exceptions thrown
    static ExperimentalDesign load(const TextFile& text_file, const bool require_spectra_file, std::string filename);

    private:
    static bool isOneTableFile_(const TextFile& text_file);
    static ExperimentalDesign parseOneTableFile_(const TextFile& text_file, const std::string& tsv_file, bool require_spectra_file);
    static ExperimentalDesign parseTwoTableFile_(const TextFile& text_file, const std::string& tsv_file, bool require_spectra_file);
    /// Reads header line of File and Sample section, checks for the existence of required headers
    /// and maps the column name to its position
    static void parseHeader_(
      const StringList &header,
      const std::string &filename,
      std::map <std::string, Size> &column_map,
      const std::set <std::string> &required,
      const std::set <std::string> &optional,
      bool allow_other_header);

    /// Throws @class ParseError with @p filename and @p message if @p test is false.
    static void parseErrorIf_(const bool test, const std::string &filename, const std::string &message);
  };
}

