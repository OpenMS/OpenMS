// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>

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

          @param filename the filename were the unimod xml file should be read from
          @param modifications the modifications which are read from the file
          @throw FileNotFound is thrown if the file could not be found
          @throw ParseError is thrown if the given file could not be parsed

      @ingroup FileIO
    */
    void load(const String & filename, std::vector<ResidueModification *> & modifications);

private:

    ///Not implemented
    UnimodXMLFile(const UnimodXMLFile & rhs);
    ///Not implemented
    UnimodXMLFile & operator=(const UnimodXMLFile & rhs);

  };

} // namespace OpenMS

