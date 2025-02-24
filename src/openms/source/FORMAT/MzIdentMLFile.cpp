// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch, Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/VALIDATORS/MzIdentMLValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/HANDLERS/MzIdentMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/MzIdentMLDOMHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>

namespace OpenMS
{

  MzIdentMLFile::MzIdentMLFile() :
    XMLFile("/SCHEMAS/mzIdentML1.1.0.xsd", "1.1.0")
  {
  }

  MzIdentMLFile::~MzIdentMLFile() = default;

  void MzIdentMLFile::load(const String& filename, std::vector<ProteinIdentification>& poid, std::vector<PeptideIdentification>& peid)
  {
    Internal::MzIdentMLDOMHandler handler(poid, peid, schema_version_, *this);
    handler.readMzIdentMLFile(filename);
  }

  void MzIdentMLFile::store(const String& filename, const std::vector<ProteinIdentification>& poid, const std::vector<PeptideIdentification>& peid) const
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::MZIDENTML))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::MZIDENTML) + "'");
    }

    Internal::MzIdentMLHandler handler(poid, peid, filename, schema_version_, *this);
    save_(filename, &handler);
//    Internal::MzIdentMLDOMHandler handler(poid, peid, schema_version_, *this);
//    handler.writeMzIdentMLFile(filename);
  }

  bool MzIdentMLFile::isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings)
  {
    // load mapping
    CVMappings mapping;
    CVMappingFile().load(File::find("/MAPPING/mzIdentML-mapping.xml"), mapping);

    // validate
    Internal::MzIdentMLValidator v(mapping, ControlledVocabulary::getPSIMSCV());
    bool result = v.validate(filename, errors, warnings);

    return result;
  }

} // namespace OpenMS
