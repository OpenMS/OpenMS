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
      is parsed with its own handler as the file with the other slices left out, and the modifications are appended to
      @p modifications in the order of the file. Only for a regular file whose text the handler provably reads the same
      way in slices: an XML declaration and comments before the root element only, no comments, CDATA sections,
      processing instructions or document type declaration after that, nothing the handler acts on outside the
      section, and in it complete \<umod:mod\> elements with specificities, neutral losses, deltas and elements nested
      as in unimod.xml and positions the handler knows. A file that is not like that, any error in a slice (whose
      messages are not written) or any failure (e.g. memory) leaves @p modifications unchanged and returns false, so
      that the result, messages and errors are those of the whole file.
    */
    bool loadInSlices_(const std::string & file, std::vector<ResidueModification *> & modifications);

    ///Not implemented
    UnimodXMLFile(const UnimodXMLFile & rhs);
    ///Not implemented
    UnimodXMLFile & operator=(const UnimodXMLFile & rhs);

  };

} // namespace OpenMS

