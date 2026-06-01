// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch, Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/VALIDATORS/MzIdentMLValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/HANDLERS/MzIdentMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/MzIdentMLDOMHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/SYSTEM/File.h>


namespace OpenMS
{

  MzIdentMLFile::MzIdentMLFile() :
    XMLFile("/SCHEMAS/mzIdentML1.3.0.xsd", "1.3.0")
  {
  }

  MzIdentMLFile::~MzIdentMLFile() = default;

  void MzIdentMLFile::load(const String& filename, std::vector<ProteinIdentification>& poid, PeptideIdentificationList& peid)
  {
    Internal::MzIdentMLDOMHandler handler(poid, peid, schema_version_, *this);
    handler.readMzIdentMLFile(filename);
  }

  void MzIdentMLFile::store(const String& filename, const std::vector<ProteinIdentification>& poid, const PeptideIdentificationList& peid) const
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

  String MzIdentMLFile::detectVersion(const String& filename) const
  {
    // Read the first few lines of the header and look for the mzIdentML version. The version is
    // declared both as a 'version="x.y.z"' attribute on the <MzIdentML> root and in the target
    // namespace ('.../mzIdentML/x.y'). We prefer the explicit version attribute and fall back to
    // the namespace; if neither maps to a bundled schema we keep the adapter default (1.3.0).
    TextFile file(filename, true, 15);
    String header;
    header.concatenate(file.begin(), file.end(), " ");

    // candidate versions we ship schemas for, newest first
    static const std::vector<String> known_versions = {"1.3.0", "1.2.0", "1.1.0", "1.0.0"};

    // 1) explicit version attribute, e.g. version="1.1.0"
    for (const String& v : known_versions)
    {
      if (header.hasSubstring(String("version=\"") + v + "\""))
      {
        return v;
      }
    }
    // 2) target namespace, e.g. http://psidev.info/psi/pi/mzIdentML/1.1
    for (const String& v : known_versions)
    {
      // namespace carries only major.minor (e.g. "1.1"), so match on that prefix
      String ns_version = v.prefix(v.rfind('.')); // "1.1.0" -> "1.1"
      if (header.hasSubstring(String("mzIdentML/") + ns_version))
      {
        return v;
      }
    }
    // 3) fall back to the adapter default
    return schema_version_;
  }

  bool MzIdentMLFile::isValid(const String& filename, std::ostream& os, String& used_version)
  {
    used_version = detectVersion(filename);
    // map the detected version to its schema; fall back to the default schema if not bundled
    String schema = "/SCHEMAS/mzIdentML" + used_version + ".xsd";
    if (!File::exists(File::getOpenMSDataPath() + schema))
    {
      schema = schema_location_;
      used_version = schema_version_;
    }
    return XMLValidator().isValid(filename, File::find(schema), os);
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
