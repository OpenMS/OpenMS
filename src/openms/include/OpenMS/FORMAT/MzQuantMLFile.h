// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/METADATA/MSQuantifications.h>

namespace OpenMS
{
  /**
      @brief File adapter for MzQuantML files

      If a critical error occurs due to the missing functionality, Exception::NotImplemented is thrown.

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MzQuantMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    MzQuantMLFile();
    ///Destructor
    ~MzQuantMLFile() override;

    /**
        @brief Loads a map from a MzQuantML file.

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String & filename, MSQuantifications & msq);

    /**
        @brief Stores a map in a MzQuantML file.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const MSQuantifications & cmsq) const;

    //~ this is not overwritten, XMLFile version works fine
    //~ bool isValid(const String& filename, std::ostream& os = std::cerr);

    /**
        @brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

        @param filename File name of the file to be checked.
        @param errors Errors during the validation are returned in this output parameter.
        @param warnings Warnings during the validation are returned in this output parameter.

        @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings);

private:

  };

} // namespace OpenMS

