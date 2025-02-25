// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/FeatureXMLHandler.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <fstream>

using namespace std;

namespace OpenMS
{
  FeatureXMLFile::FeatureXMLFile() :
    Internal::XMLFile("/SCHEMAS/FeatureXML_1_9.xsd", "1.9")
  {
  }

  FeatureXMLFile::~FeatureXMLFile() = default;

  Size FeatureXMLFile::loadSize(const String& filename)
  {
    FeatureMap dummy;
    Internal::FeatureXMLHandler handler(dummy, filename);
    handler.setOptions(options_);
    handler.setSizeOnly(true);
    handler.setLogType(getLogType());
    parse_(filename, &handler);

    return handler.getSize();
  }

  void FeatureXMLFile::load(const String& filename, FeatureMap& feature_map)
  {
    feature_map.clear(true);
    //set DocumentIdentifier
    feature_map.setLoadedFileType(filename);
    feature_map.setLoadedFilePath(filename);

    Internal::FeatureXMLHandler handler(feature_map, filename);
    handler.setOptions(options_);
    handler.setLogType(getLogType());
    parse_(filename, &handler);

    // !!! Hack: set feature FWHM from meta info entries as
    // long as featureXML doesn't support a width entry.
    // See also hack in BaseFeature::setWidth().
    for (auto& feature : feature_map)
    {
      if (feature.metaValueExists("FWHM"))
      {
        feature.setWidth((double)feature.getMetaValue("FWHM"));
      }
    }

    // put ranges into defined state
    feature_map.updateRanges();
  }

  void FeatureXMLFile::store(const String& filename, const FeatureMap& feature_map)
  {

    if (!FileHandler::hasValidExtension(filename, FileTypes::FEATUREXML))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::FEATUREXML) + "'");
    }

    if (Size invalid_unique_ids = feature_map.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId))
    {

      // TODO Take care *outside* that this does not happen.
      // We can detect this here but it is too late to fix the problem;
      // there is no straightforward action to be taken in all cases.
      // Note also that we are given a const reference.
      OPENMS_LOG_INFO << String("FeatureXMLHandler::store():  found ") + invalid_unique_ids + " invalid unique ids" << std::endl;
    }

    // This will throw if the unique ids are not unique,
    // so we never create bad files in this respect.
    try
    {
      feature_map.updateUniqueIdToIndex();
    }
    catch (Exception::Postcondition& e)
    {
      OPENMS_LOG_FATAL_ERROR << e.getName() << ' ' << e.what() << std::endl;
      throw;
    }

    Internal::FeatureXMLHandler handler(feature_map, filename);
    handler.setOptions(options_);
    handler.setLogType(getLogType());
    save_(filename, &handler);
  }

  FeatureFileOptions& FeatureXMLFile::getOptions()
  {
    return options_;
  }

  const FeatureFileOptions& FeatureXMLFile::getOptions() const
  {
    return options_;
  }

  void FeatureXMLFile::setOptions(const FeatureFileOptions& options)
  {
    options_ = options;
  }


}
