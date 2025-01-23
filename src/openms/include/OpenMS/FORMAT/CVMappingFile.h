// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>

#include <vector>

namespace OpenMS
{
  class String;

  /**
    @brief Used to load CvMapping files

    This file contains the mapping of CV terms to the schema, which
    is used by PSI standard formats to semantically validate files.

    A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

    @ingroup FileIO
  */
  class OPENMS_DLLAPI CVMappingFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Default constructor
    CVMappingFile();

    /// Destructor
    ~CVMappingFile() override;

    /**
        @brief loads CvMappings from the given file

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
        @param strip_namespaces if enable, namespace definitions of the paths are eliminated, e.g. 'pf:cvParam' -> 'cvParam'
        @param cv_mappings  The CVMappings instance in which the rules, cvs and other content from the file should be stored
        @param filename  The filename to read from
    */
    void load(const String& filename, CVMappings& cv_mappings, bool strip_namespaces = false);

protected:

    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    // Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    // Docu in base class
    void characters(const XMLCh* const chars, const XMLSize_t /*length*/) override;

private:

    ///Not implemented
    CVMappingFile(const CVMappingFile& rhs);

    ///Not implemented
    CVMappingFile& operator=(const CVMappingFile& rhs);

    String tag_;

    bool strip_namespaces_;

    CVMappingRule actual_rule_;

    std::vector<CVMappingRule> rules_;

    std::vector<CVReference> cv_references_;

  };

} // namespace OpenMS

