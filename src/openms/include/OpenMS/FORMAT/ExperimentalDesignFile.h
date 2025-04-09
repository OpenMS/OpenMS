// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:	Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Lukas Zimmermann $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>
#include <set>

namespace OpenMS
{
  class ExperimentalDesign;
  class TextFile; 
  /**
  @brief Load an experimental design from a TSV file. (see ExperimentalDesign for details on the supported format)

  @ingroup Format
  */

  class OPENMS_DLLAPI ExperimentalDesignFile
  {
    public:

    /// Loads an experimental design from a tabular separated file
    static ExperimentalDesign load(const String& tsv_file, bool require_spectra_files);

    /// Loads an experimental design from an already loaded or generated, tabular file
    static ExperimentalDesign load(const TextFile& text_file, const bool require_spectra_file, String filename);

    private:
    static bool isOneTableFile_(const TextFile& text_file);
    static ExperimentalDesign parseOneTableFile_(const TextFile& text_file, const String& tsv_file, bool require_spectra_file);
    static ExperimentalDesign parseTwoTableFile_(const TextFile& text_file, const String& tsv_file, bool require_spectra_file);
    /// Reads header line of File and Sample section, checks for the existence of required headers
    /// and maps the column name to its position
    static void parseHeader_(
      const StringList &header,
      const String &filename,
      std::map <String, Size> &column_map,
      const std::set <String> &required,
      const std::set <String> &optional,
      bool allow_other_header);

    /// Throws @class ParseError with @p filename and @p message if @p test is false.
    static void parseErrorIf_(const bool test, const String &filename, const String &message);
  };
}

