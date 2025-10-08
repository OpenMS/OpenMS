// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>

#include <vector>

namespace OpenMS
{

  /**
    @brief Used to load and store TransformationXML files

    This class is used to load and store documents that implement the schema of
    TransformationXML files.

    A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

    @ingroup FileIO
  */
  class OPENMS_DLLAPI TransformationXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Constructor
    TransformationXMLFile();

    /**
    @brief Loads the transformation from an TransformationXML file

    The information is read in and the information is stored in the
    corresponding variables

    @exception Exception::FileNotFound is thrown if the file could not be opened
    @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String& filename, TransformationDescription& transformation, bool fit_model=true);

    /**
    @brief Stores the data in an TransformationXML file

    The data is read in and stored in the file named 'filename'.

    @exception Exception::UnableToCreateFile is thrown if the file could not be created
    @exception Exception::IllegalArgument is thrown if unsupported parameter types have been set
    */
    void store(const String& filename, const TransformationDescription& transformation);

protected:
    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    /// @name Members for use during loading data
    //@{
    /// Param to fill in during parse
    Param params_;
    /// Data vector
    TransformationDescription::DataPoints data_;
    /// Model type
    String model_type_;
    //@}

  };

} // namespace OpenMS

