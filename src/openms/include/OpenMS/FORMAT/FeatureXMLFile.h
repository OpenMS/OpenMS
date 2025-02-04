// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <iosfwd>

namespace OpenMS
{
  class Feature;
  class FeatureMap;

  /**
    @brief This class provides Input/Output functionality for feature maps

    A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

    @todo Take care that unique ids are assigned properly by TOPP tools before
    calling FeatureXMLFile::store().  There will be a message on OPENMS_LOG_INFO but
    we will make no attempt to fix the problem in this class.  (all developers)

    @ingroup FileIO
  */
  class OPENMS_DLLAPI FeatureXMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {

public:

    /** @name Constructors and Destructor */
    //@{
    ///Default constructor
    FeatureXMLFile();
    ///Destructor
    ~FeatureXMLFile() override;
    //@}

    /**
        @brief loads the file with name @p filename into @p map and calls updateRanges().

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, FeatureMap& feature_map);

    Size loadSize(const String& filename);

    /**
        @brief stores the map @p feature_map in file with name @p filename.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename, const FeatureMap& feature_map);

    /// Mutable access to the options for loading/storing
    FeatureFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const FeatureFileOptions& getOptions() const;

    /// setter for options for loading/storing
    void setOptions(const FeatureFileOptions&);

protected:

    /// Options that can be set
    FeatureFileOptions options_;

  };

} // namespace OpenMS

