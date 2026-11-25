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
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/FORMAT/HANDLERS/MzIdentMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/MzIdentMLDOMHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/METADATA/ProteinRunReferences.h>


namespace OpenMS
{

  MzIdentMLFile::MzIdentMLFile() :
    XMLFile("/SCHEMAS/mzIdentML1.3.0.xsd", "1.3.0")
  {
  }

  MzIdentMLFile::~MzIdentMLFile() = default;

  void MzIdentMLFile::load(const std::string& filename, std::vector<ProteinIdentification>& poid, PeptideIdentificationList& peid)
  {
    // replace (not accumulate): the DOM handler only appends, so clear the outputs first to match the
    // behavior of the sibling identification loaders (e.g. IdXMLFile::load) and FileHandler::loadIdentifications.
    poid.clear();
    peid.clear();
    Internal::MzIdentMLDOMHandler handler(poid, peid, schema_version_, *this);
    handler.readMzIdentMLFile(filename);
  }

  void MzIdentMLFile::store(const std::string& filename, const std::vector<ProteinIdentification>& poid, const PeptideIdentificationList& peid) const
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::MZIDENTML))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::MZIDENTML) + "'");
    }

    // every peptide identification needs its protein identification run: fail before the file is opened
    ProteinRunReferences::check(poid, peid);

    Internal::MzIdentMLHandler handler(poid, peid, filename, schema_version_, *this);
    save_(filename, &handler);
//    Internal::MzIdentMLDOMHandler handler(poid, peid, schema_version_, *this);
//    handler.writeMzIdentMLFile(filename);
  }

  namespace
  {
    /// Reads the version attribute and the namespace of the root element, then stops parsing
    class RootVersionHandler final : public Internal::XMLHandler
    {
    public:
      explicit RootVersionHandler(const std::string& filename) :
        XMLHandler(filename, "")
      {
      }

      void onStartElement(const char16_t* qname, const Internal::XMLAttributes& attributes) override
      {
        attributes.valueUTF8("version", version);
        // the default namespace, or the one of the root's prefix (e.g. <mzid:MzIdentML xmlns:mzid="...">)
        const std::string name = sm_.convert(qname);
        const std::string::size_type colon = name.find(':');
        const std::string xmlns_attribute = colon == std::string::npos ? std::string("xmlns") : "xmlns:" + name.substr(0, colon);
        attributes.valueUTF8(xmlns_attribute.c_str(), xmlns);
        throw EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }

      std::string version; ///< version attribute of the root (empty if absent)
      std::string xmlns; ///< namespace of the root (empty if not declared)
    };
  }

  std::string MzIdentMLFile::detectVersion(const std::string& filename) const
  {
    // The mzIdentML version is declared by the root element, both as its 'version="x.y.z"' attribute and in its
    // namespace ('.../mzIdentML/x.y'). Only the root counts: other elements have version attributes of their own
    // (e.g. <AnalysisSoftware version="1.3.0">). We prefer the version attribute and fall back to the namespace;
    // if neither maps to a bundled schema we keep the adapter default (1.3.0).
    RootVersionHandler root(filename);
    try
    {
      MzIdentMLFile reader; // parse_() is not const
      reader.parse_(filename, &root); // stops at the root element
    }
    catch (const Exception::ParseError&)
    {
      // not well-formed up to the root element: use the default
    }

    // candidate versions we ship schemas for, newest first
    static const std::vector<std::string> known_versions = {"1.3.0", "1.2.0", "1.1.0", "1.0.0"};

    // 1) explicit version attribute, e.g. version="1.1.0"
    for (const std::string& v : known_versions)
    {
      if (root.version == v)
      {
        return v;
      }
    }
    // 2) target namespace, e.g. http://psidev.info/psi/pi/mzIdentML/1.1
    for (const std::string& v : known_versions)
    {
      // namespace carries only major.minor (e.g. "1.1"), so match on that prefix
      std::string ns_version = StringUtils::prefix(v, v.rfind('.')); // "1.1.0" -> "1.1"
      if (StringUtils::hasSubstring(root.xmlns, std::string("mzIdentML/") + ns_version))
      {
        return v;
      }
    }
    // 3) fall back to the adapter default
    return schema_version_;
  }

  bool MzIdentMLFile::isValid(const std::string& filename, std::ostream& os, std::string& used_version)
  {
    used_version = detectVersion(filename);
    // map the detected version to its schema; fall back to the default schema if not bundled
    std::string schema = "/SCHEMAS/mzIdentML" + used_version + ".xsd";
    if (!File::exists(File::getOpenMSDataPath() + schema))
    {
      schema = schema_location_;
      used_version = schema_version_;
    }
    return XMLValidator().isValid(filename, File::find(schema), os);
  }

  bool MzIdentMLFile::isSemanticallyValid(const std::string& filename, StringList& errors, StringList& warnings)
  {
    // load mapping
    CVMappings mapping;
    CVMappingFile().load(File::find("/MAPPING/mzIdentML-mapping.xml"), mapping);
    // mzIdentML 1.2 made some terms required for protein grouping: check them only for files of version 1.2 and later
    const std::string version = detectVersion(filename);
    if (version != "1.0.0" && version != "1.1.0")
    {
      CVMappings mapping_1_2;
      CVMappingFile().load(File::find("/MAPPING/mzIdentML-1.2-mapping.xml"), mapping_1_2);
      for (const CVMappingRule& rule : mapping_1_2.getMappingRules())
      {
        mapping.addMappingRule(rule);
      }
    }

    // validate
    Internal::MzIdentMLValidator v(mapping, ControlledVocabulary::getPSIMSCV());
    bool result = v.validate(filename, errors, warnings);

    return result;
  }

} // namespace OpenMS
