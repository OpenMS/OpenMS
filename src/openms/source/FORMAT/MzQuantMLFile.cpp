// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/HANDLERS/MzQuantMLHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/MzQuantMLValidator.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{

  MzQuantMLFile::MzQuantMLFile() :
    XMLFile("/SCHEMAS/mzQuantML_1_0_0-rc2", "1.0.0")
  {
  }

  MzQuantMLFile::~MzQuantMLFile() = default;

  void MzQuantMLFile::load(const String & filename, MSQuantifications & msq)
  {
    Internal::MzQuantMLHandler handler(msq, filename, schema_version_, *this);
    parse_(filename, &handler);
  }

  void MzQuantMLFile::store(const String & filename, const MSQuantifications & cmsq) const
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::MZQUANTML))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::MZIDENTML) + "'");
    }

    Internal::MzQuantMLHandler handler(cmsq, filename, schema_version_, *this);
    save_(filename, &handler);
  }

  bool MzQuantMLFile::isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings)
  {
    //load mapping
    CVMappings mapping;
    CVMappingFile().load(File::find("/MAPPING/mzQuantML-mapping_1.0.0-rc2-general.xml"), mapping);

    //load cvs
    ControlledVocabulary cv;
    cv.loadFromOBO("MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("PATO", File::find("/CV/quality.obo"));
    cv.loadFromOBO("UO", File::find("/CV/unit.obo"));
    cv.loadFromOBO("BTO", File::find("/CV/brenda.obo"));
    cv.loadFromOBO("GO", File::find("/CV/goslim_goa.obo"));

    //validate TODO
    Internal::MzQuantMLValidator v(mapping, cv);
    bool result = v.validate(filename, errors, warnings);

    return result;
  }

} // namespace OpenMS
