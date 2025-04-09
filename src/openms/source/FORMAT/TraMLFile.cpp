// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TraMLFile.h>

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/VALIDATORS/TraMLValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/HANDLERS/TraMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{

  TraMLFile::TraMLFile() :
    XMLFile("/SCHEMAS/TraML1.0.0.xsd", "1.0.0")
  {
  }

  TraMLFile::~TraMLFile() = default;

  void TraMLFile::load(const String & filename, TargetedExperiment & exp)
  {
    Internal::TraMLHandler handler(exp, filename, schema_version_, *this);
    parse_(filename, &handler);
  }

  void TraMLFile::store(const String & filename, const TargetedExperiment & exp) const
  {
    Internal::TraMLHandler handler(exp, filename, schema_version_, *this);
    save_(filename, &handler);
  }

  bool TraMLFile::isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings)
  {
    //load mapping
    CVMappings mapping;
    CVMappingFile().load(File::find("/MAPPING/TraML-mapping.xml"), mapping);

    //load cvs
    ControlledVocabulary cv;
    cv.loadFromOBO("MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("UO", File::find("/CV/unit.obo"));

    //validate
    Internal::TraMLValidator v(mapping, cv);
    bool result = v.validate(filename, errors, warnings);

    return result;
  }

} // namespace OpenMS

