// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/SpectrumSettings.h>

#include <OpenMS/CONCEPT/Helpers.h>
#include <boost/iterator/indirect_iterator.hpp> // for equality

using namespace std;

namespace OpenMS
{

  const std::string SpectrumSettings::NamesOfSpectrumType[] = {"Unknown", "Centroid", "Profile"};

  SpectrumSettings::SpectrumSettings() :
    MetaInfoInterface(),
    type_(UNKNOWN),
    native_id_(),
    comment_(),
    instrument_settings_(),
    source_file_(),
    acquisition_info_(),
    precursors_(),
    products_(),
    identification_(),
    data_processing_()
  {
  }

  SpectrumSettings::~SpectrumSettings() = default;

  bool SpectrumSettings::operator==(const SpectrumSettings & rhs) const
  {
    return MetaInfoInterface::operator==(rhs) &&
           type_ == rhs.type_ &&
           native_id_ == rhs.native_id_ &&
           comment_ == rhs.comment_ &&
           instrument_settings_ == rhs.instrument_settings_ &&
           acquisition_info_ == rhs.acquisition_info_ &&
           source_file_ == rhs.source_file_ &&
           precursors_ == rhs.precursors_ &&
           products_ == rhs.products_ &&
           identification_ == rhs.identification_ &&
           ( data_processing_.size() == rhs.data_processing_.size() &&
           std::equal(data_processing_.begin(),
                      data_processing_.end(),
                      rhs.data_processing_.begin(), 
                      OpenMS::Helpers::cmpPtrSafe<DataProcessingPtr>) );
  }

  bool SpectrumSettings::operator!=(const SpectrumSettings & rhs) const
  {
    return !(operator==(rhs));
  }

  void SpectrumSettings::unify(const SpectrumSettings & rhs)
  {
    // append metavalues (overwrite when already present)
    std::vector<UInt> keys;
    rhs.getKeys(keys);
    for (Size i = 0; i < keys.size(); ++i)
    {
      setMetaValue(keys[i], rhs.getMetaValue(keys[i]));
    }

    if (type_ != rhs.type_)
    {
      type_ = UNKNOWN;                       // only keep if both are equal
    }
    //native_id_ == rhs.native_id_ // keep
    comment_ += rhs.comment_;        // append
    //instrument_settings_ == rhs.instrument_settings_  // keep
    //acquisition_info_ == rhs.acquisition_info_
    //source_file_ == rhs.source_file_ &&
    precursors_.insert(precursors_.end(), rhs.precursors_.begin(), rhs.precursors_.end());
    products_.insert(products_.end(), rhs.products_.begin(), rhs.products_.end());
    identification_.insert(identification_.end(), rhs.identification_.begin(), rhs.identification_.end());
    data_processing_.insert(data_processing_.end(), rhs.data_processing_.begin(), rhs.data_processing_.end());
  }

  SpectrumSettings::SpectrumType SpectrumSettings::getType() const
  {
    return type_;
  }

  void SpectrumSettings::setType(SpectrumSettings::SpectrumType type)
  {
    type_ = type;
  }

  const String & SpectrumSettings::getComment() const
  {
    return comment_;
  }

  void SpectrumSettings::setComment(const String & comment)
  {
    comment_ = comment;
  }

  const InstrumentSettings & SpectrumSettings::getInstrumentSettings() const
  {
    return instrument_settings_;
  }

  InstrumentSettings & SpectrumSettings::getInstrumentSettings()
  {
    return instrument_settings_;
  }

  void SpectrumSettings::setInstrumentSettings(const InstrumentSettings & instrument_settings)
  {
    instrument_settings_ = instrument_settings;
  }

  const AcquisitionInfo & SpectrumSettings::getAcquisitionInfo() const
  {
    return acquisition_info_;
  }

  AcquisitionInfo & SpectrumSettings::getAcquisitionInfo()
  {
    return acquisition_info_;
  }

  void SpectrumSettings::setAcquisitionInfo(const AcquisitionInfo & acquisition_info)
  {
    acquisition_info_ = acquisition_info;
  }

  const SourceFile & SpectrumSettings::getSourceFile() const
  {
    return source_file_;
  }

  SourceFile & SpectrumSettings::getSourceFile()
  {
    return source_file_;
  }

  void SpectrumSettings::setSourceFile(const SourceFile & source_file)
  {
    source_file_ = source_file;
  }

  const vector<Precursor> & SpectrumSettings::getPrecursors() const
  {
    return precursors_;
  }

  vector<Precursor> & SpectrumSettings::getPrecursors()
  {
    return precursors_;
  }

  void SpectrumSettings::setPrecursors(const vector<Precursor> & precursors)
  {
    precursors_ = precursors;
  }

  const vector<Product> & SpectrumSettings::getProducts() const
  {
    return products_;
  }

  vector<Product> & SpectrumSettings::getProducts()
  {
    return products_;
  }

  void SpectrumSettings::setProducts(const vector<Product> & products)
  {
    products_ = products;
  }

  std::ostream & operator<<(std::ostream & os, const SpectrumSettings & /*spec*/)
  {
    os << "-- SPECTRUMSETTINGS BEGIN --" << std::endl;
    os << "-- SPECTRUMSETTINGS END --" << std::endl;
    return os;
  }

  const std::vector<PeptideIdentification> & SpectrumSettings::getPeptideIdentifications() const
  {
    return identification_;
  }

  std::vector<PeptideIdentification> & SpectrumSettings::getPeptideIdentifications()
  {
    return identification_;
  }

  void SpectrumSettings::setPeptideIdentifications(const std::vector<PeptideIdentification> & identification)
  {
    identification_ = identification;
  }

  const String & SpectrumSettings::getNativeID() const
  {
    return native_id_;
  }

  void SpectrumSettings::setNativeID(const String & native_id)
  {
    native_id_ = native_id;
  }

  void SpectrumSettings::setDataProcessing(const std::vector< DataProcessingPtr > & data_processing)
  {
    data_processing_ = data_processing;
  }

  std::vector< DataProcessingPtr > & SpectrumSettings::getDataProcessing()
  {
    return data_processing_;
  }

  const std::vector< boost::shared_ptr<const DataProcessing > > SpectrumSettings::getDataProcessing() const
  {
    return OpenMS::Helpers::constifyPointerVector(data_processing_);
  }

}

