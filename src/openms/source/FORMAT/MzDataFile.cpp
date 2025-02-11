// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/VALIDATORS/MzDataValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

namespace OpenMS
{
  MzDataFile::MzDataFile() :
    XMLFile("/SCHEMAS/mzData_1_05.xsd", "1.05"),
    options_()
  {
  }

  MzDataFile::~MzDataFile() = default;

  PeakFileOptions & MzDataFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions & MzDataFile::getOptions() const
  {
    return options_;
  }

  void MzDataFile::setOptions(const PeakFileOptions & options)
  {
      options_ = options;
  }

  bool MzDataFile::isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings)
  {
    //load mapping
    CVMappings mapping;
    CVMappingFile().load(File::find("/MAPPING/mzdata-mapping.xml"), mapping);

    //load cvs
    ControlledVocabulary cv;
    cv.loadFromOBO("PSI", File::find("/CV/psi-mzdata.obo"));

    //validate
    Internal::MzDataValidator v(mapping, cv);
    bool result = v.validate(filename, errors, warnings);

    return result;
  }

  void MzDataFile::load(const String & filename, PeakMap & map)
  {
    map.reset();

    //set DocumentIdentifier
    map.setLoadedFileType(filename);
    map.setLoadedFilePath(filename);

    Internal::MzDataHandler handler(map, filename, schema_version_, *this);
    handler.setOptions(options_);
    parse_(filename, &handler);
  }

  void MzDataFile::store(const String & filename, const PeakMap & map) const
  {
    Internal::MzDataHandler handler(map, filename, schema_version_, *this);
    handler.setOptions(options_);
    save_(filename, &handler);
  }


} // namespace OpenMS
