// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzXMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/MzXMLHandler.h>

using namespace std;

namespace OpenMS
{
  MzXMLFile::MzXMLFile() :
    XMLFile("/SCHEMAS/mzXML_idx_3.1.xsd", "3.1")
  {
  }

  MzXMLFile::~MzXMLFile() = default;

  PeakFileOptions & MzXMLFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions & MzXMLFile::getOptions() const
  {
    return options_;
  }

  void MzXMLFile::setOptions(const PeakFileOptions & options)
  {
      options_ = options;
  }

  void MzXMLFile::load(const String & filename, MapType & map)
  {
    map.reset();

    //set DocumentIdentifier
    map.setLoadedFileType(filename);
    map.setLoadedFilePath(filename);

    Internal::MzXMLHandler handler(map, filename, schema_version_, *this);
    handler.setOptions(options_);
    parse_(filename, &handler);
  }

  void MzXMLFile::store(const String & filename, const MapType & map) const
  {
    Internal::MzXMLHandler handler(map, filename, schema_version_, *this);
    handler.setOptions(options_);
    save_(filename, &handler);
  }

  void MzXMLFile::transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count)
  {
    // First pass through the file -> get the meta-data and hand it to the consumer
    transformFirstPass_(filename_in, consumer, skip_full_count);
    
    // Second pass through the data, now read the spectra!
    {
      MapType dummy;
      Internal::MzXMLHandler handler(dummy, filename_in, getVersion(), *this);
      handler.setOptions(options_);
      handler.setMSDataConsumer(consumer);
      parse_(filename_in, &handler);
    }
  }

  void MzXMLFile::transform(const String& filename_in, Interfaces::IMSDataConsumer * consumer, MapType& map, bool skip_full_count)
  {
    // First pass through the file -> get the meta-data and hand it to the consumer
    transformFirstPass_(filename_in, consumer, skip_full_count);

    // Second pass through the data, now read the spectra!
    {
      PeakFileOptions tmp_options(options_);
      Internal::MzXMLHandler handler(map, filename_in, getVersion(), *this);
      tmp_options.setAlwaysAppendData(true);
      handler.setOptions(tmp_options);
      handler.setMSDataConsumer(consumer);

      parse_(filename_in, &handler);
    }
  }

  void MzXMLFile::transformFirstPass_(const String& filename_in, Interfaces::IMSDataConsumer * consumer, bool skip_full_count)
  {
    // Create temporary objects and counters
    PeakFileOptions tmp_options(options_);
    Size scount = 0, ccount = 0;
    MapType experimental_settings;
    Internal::MzXMLHandler handler(experimental_settings, filename_in, getVersion(), *this);

    // set temporary options for handler
    tmp_options.setMetadataOnly( skip_full_count );
    handler.setOptions(tmp_options);
    handler.setLoadDetail(Internal::XMLHandler::LD_RAWCOUNTS);
    parse_(filename_in, &handler);

    // After parsing, collect information
    scount = handler.getScanCount();
    consumer->setExpectedSize(scount, ccount);
    consumer->setExperimentalSettings(experimental_settings);
  }

} // namespace OpenMS

