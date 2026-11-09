// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------

#pragma once

#include <map>
#include <string>
#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Private helpers shared by the Inspect and Sequest input/output file classes.

      Not installed. InspectInfile/SequestInfile::handlePTMs and
      InspectOutfile/SequestOutfile::getACAndACType forward here.
    */
    namespace SearchEngineFileHelpers
    {
      /**
        @brief Parse a modification line into @p PTMname_residues_mass_type.

        See InspectInfile::handlePTMs / SequestInfile::handlePTMs for the format.
        @p function is the name of the calling member function; it is reported by the thrown exceptions.

        @throw Exception::FileNotFound if @p modifications_filename is needed but does not exist
        @throw Exception::FileNotReadable if @p modifications_filename is needed but not readable
        @throw Exception::ParseError if a modification cannot be parsed (@p PTMname_residues_mass_type is cleared)
      */
      void handlePTMs(const std::string& modification_line,
                      const std::string& modifications_filename,
                      bool monoisotopic,
                      std::map<std::string, std::vector<std::string> >& PTMname_residues_mass_type,
                      const char* function);

      /// Retrieve the accession and the accession type of a protein from a (FASTA header) line.
      void getACAndACType(std::string line, std::string& accession, std::string& accession_type);
    }
  }
}
