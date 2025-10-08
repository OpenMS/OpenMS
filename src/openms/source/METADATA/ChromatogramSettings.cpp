// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ChromatogramSettings.h>

#include <OpenMS/CONCEPT/Helpers.h>
#include <boost/iterator/indirect_iterator.hpp> // for equality

using namespace std;

namespace OpenMS
{

  // keep this in sync with enum ChromatogramType
  const char * const ChromatogramSettings::ChromatogramNames[] = {"mass chromatogram", "total ion current chromatogram", "selected ion current chromatogram" ,"base peak chromatogram",
                                                                  "selected ion monitoring chromatogram" ,"selected reaction monitoring chromatogram" ,"electromagnetic radiation chromatogram",
                                                                  "absorption chromatogram", "emission chromatogram", "unknown chromatogram"}; // last entry should be "unknown", since this is the default in FileInfo.cpp

  ChromatogramSettings::ChromatogramSettings() :
    MetaInfoInterface(),
    native_id_(),
    comment_(),
    instrument_settings_(),
    source_file_(),
    acquisition_info_(),
    precursor_(),
    product_(),
    data_processing_(),
    type_(MASS_CHROMATOGRAM)
  {
  }

  ChromatogramSettings::~ChromatogramSettings() = default;

  bool ChromatogramSettings::operator==(const ChromatogramSettings & rhs) const
  {
    return MetaInfoInterface::operator==(rhs) &&
           native_id_ == rhs.native_id_ &&
           comment_ == rhs.comment_ &&
           instrument_settings_ == rhs.instrument_settings_ &&
           acquisition_info_ == rhs.acquisition_info_ &&
           source_file_ == rhs.source_file_ &&
           precursor_ == rhs.precursor_ &&
           product_ == rhs.product_ &&
           // We are not interested whether the pointers are equal but whether
           // the contents are equal
           ( data_processing_.size() == rhs.data_processing_.size() &&
           std::equal( boost::make_indirect_iterator(data_processing_.begin()),
                       boost::make_indirect_iterator(data_processing_.end()),
                       boost::make_indirect_iterator(rhs.data_processing_.begin()) ) ) &&
           type_ == rhs.type_;
  }

  bool ChromatogramSettings::operator!=(const ChromatogramSettings & rhs) const
  {
    return !(operator==(rhs));
  }

  const String & ChromatogramSettings::getComment() const
  {
    return comment_;
  }

  void ChromatogramSettings::setComment(const String & comment)
  {
    comment_ = comment;
  }

  const InstrumentSettings & ChromatogramSettings::getInstrumentSettings() const
  {
    return instrument_settings_;
  }

  InstrumentSettings & ChromatogramSettings::getInstrumentSettings()
  {
    return instrument_settings_;
  }

  void ChromatogramSettings::setInstrumentSettings(const InstrumentSettings & instrument_settings)
  {
    instrument_settings_ = instrument_settings;
  }

  const AcquisitionInfo & ChromatogramSettings::getAcquisitionInfo() const
  {
    return acquisition_info_;
  }

  AcquisitionInfo & ChromatogramSettings::getAcquisitionInfo()
  {
    return acquisition_info_;
  }

  void ChromatogramSettings::setAcquisitionInfo(const AcquisitionInfo & acquisition_info)
  {
    acquisition_info_ = acquisition_info;
  }

  const SourceFile & ChromatogramSettings::getSourceFile() const
  {
    return source_file_;
  }

  SourceFile & ChromatogramSettings::getSourceFile()
  {
    return source_file_;
  }

  void ChromatogramSettings::setSourceFile(const SourceFile & source_file)
  {
    source_file_ = source_file;
  }

  const Precursor & ChromatogramSettings::getPrecursor() const
  {
    return precursor_;
  }

  Precursor & ChromatogramSettings::getPrecursor()
  {
    return precursor_;
  }

  void ChromatogramSettings::setPrecursor(const Precursor & precursor)
  {
    precursor_ = precursor;
  }

  const Product & ChromatogramSettings::getProduct() const
  {
    return product_;
  }

  Product & ChromatogramSettings::getProduct()
  {
    return product_;
  }

  void ChromatogramSettings::setProduct(const Product & product)
  {
    product_ = product;
  }

  std::ostream & operator<<(std::ostream & os, const ChromatogramSettings & /*spec*/)
  {
    os << "-- CHROMATOGRAMSETTINGS BEGIN --" << std::endl;
    os << "-- CHROMATOGRAMSETTINGS END --" << std::endl;
    return os;
  }

  const String & ChromatogramSettings::getNativeID() const
  {
    return native_id_;
  }

  void ChromatogramSettings::setNativeID(const String & native_id)
  {
    native_id_ = native_id;
  }

  ChromatogramSettings::ChromatogramType ChromatogramSettings::getChromatogramType() const
  {
    return type_;
  }

  void ChromatogramSettings::setChromatogramType(ChromatogramType type)
  {
    type_ = type;
  }

  void ChromatogramSettings::setDataProcessing(const std::vector< DataProcessingPtr > & data_processing)
  {
    data_processing_ = data_processing;
  }

  std::vector< DataProcessingPtr > & ChromatogramSettings::getDataProcessing()
  {
    return data_processing_;
  }

  const std::vector< boost::shared_ptr<const DataProcessing > > ChromatogramSettings::getDataProcessing() const 
  {
    return OpenMS::Helpers::constifyPointerVector(data_processing_);
  }

}
