// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>

#include <string>
#include <vector>

namespace OpenMS
{
  class ResidueModification;

  /**
    @brief Used to load XML files from unimod.org files

    @ingroup FileIO
  */
  class OPENMS_DLLAPI UnimodXMLFile :
    public Internal::XMLFile
  {
public:

    /// Default constructor
    UnimodXMLFile();

    /// Destructor
    ~UnimodXMLFile() override;
    /**
      @brief loads data from unimod.xml file

          @param[in] filename the filename were the unimod xml file should be read from
          @param[in] modifications the modifications which are read from the file
          @throw FileNotFound is thrown if the file could not be found
          @throw ParseError is thrown if the given file could not be parsed

      @ingroup FileIO
    */
    void load(const std::string & filename, std::vector<ResidueModification *> & modifications);

private:

    /**
      @brief Parses a large UniMod file in slices on several threads; false if it did not (the caller parses it as a whole)

      The \<umod:modifications\> section of @p file is cut at end tags of \<umod:mod\> elements into a few slices; each
      is parsed as the file with the other slices left out, and the modifications are appended to @p modifications in
      the order of the file. A small file, a section that may hold such end tags as text (comments, CDATA sections,
      processing instructions), a specificity position the handler warns about and any error in a slice leave
      @p modifications unchanged and return false, so that messages and errors are those of the whole file.
    */
    bool loadInSlices_(const std::string & file, std::vector<ResidueModification *> & modifications);

    ///Not implemented
    UnimodXMLFile(const UnimodXMLFile & rhs);
    ///Not implemented
    UnimodXMLFile & operator=(const UnimodXMLFile & rhs);

  };

} // namespace OpenMS

