// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

using namespace std;

namespace OpenMS
{
  PeakFileOptions::PeakFileOptions() :
    metadata_only_(false),
    force_maxquant_compatibility_(false),
    force_tpp_compatibility_(false),
    write_supplemental_data_(true),
    has_rt_range_(false),
    has_mz_range_(false),
    has_intensity_range_(false),
    mz_32_bit_(false),
    int_32_bit_(true),
    rt_range_(),
    mz_range_(),
    intensity_range_(),
    ms_levels_(),
    zlib_compression_(false),
    always_append_data_(false),
    skip_xml_checks_(false),
    sort_spectra_by_mz_(true),
    sort_chromatograms_by_rt_(true),
    fill_data_(true),
    write_index_(true),
    np_config_mz_(),
    np_config_int_(),
    np_config_fda_(),
    maximal_data_pool_size_(100),
    precursor_mz_selected_ion_(true)
  {
  }

  PeakFileOptions::PeakFileOptions(const PeakFileOptions& options) = default;

  PeakFileOptions::~PeakFileOptions() = default;

  void PeakFileOptions::setMetadataOnly(bool only)
  {
    metadata_only_ = only;
  }

  bool PeakFileOptions::getMetadataOnly() const
  {
    return metadata_only_;
  }

  void PeakFileOptions::setForceMQCompatability(bool forceMQ)
  {
    force_maxquant_compatibility_ = forceMQ;
  }

  bool PeakFileOptions::getForceMQCompatability() const
  {
    return force_maxquant_compatibility_;
  }

  void PeakFileOptions::setForceTPPCompatability(bool forceTPP)
  {
    force_tpp_compatibility_ = forceTPP;
  }

  bool PeakFileOptions::getForceTPPCompatability() const
  {
    return force_tpp_compatibility_;
  }

  void PeakFileOptions::setWriteSupplementalData(bool write)
  {
    write_supplemental_data_ = write;
  }

  bool PeakFileOptions::getWriteSupplementalData() const
  {
    return write_supplemental_data_;
  }

  void PeakFileOptions::setRTRange(const DRange<1>& range)
  {
    rt_range_ = range;
    has_rt_range_ = !rt_range_.isEmpty();
  }

  bool PeakFileOptions::hasRTRange() const
  {
    return has_rt_range_;
  }

  const DRange<1>& PeakFileOptions::getRTRange() const
  {
    return rt_range_;
  }

  void PeakFileOptions::setMZRange(const DRange<1>& range)
  {
    mz_range_ = range;
    has_mz_range_ = true;
  }

  bool PeakFileOptions::hasMZRange() const
  {
    return has_mz_range_;
  }

  const DRange<1>& PeakFileOptions::getMZRange() const
  {
    return mz_range_;
  }

  void PeakFileOptions::setIntensityRange(const DRange<1>& range)
  {
    intensity_range_ = range;
    has_intensity_range_ = true;
  }

  bool PeakFileOptions::hasIntensityRange() const
  {
    return has_intensity_range_;
  }

  const DRange<1>& PeakFileOptions::getIntensityRange() const
  {
    return intensity_range_;
  }

  void PeakFileOptions::setMSLevels(const vector<Int>& levels)
  {
    ms_levels_ = levels;
  }

  void PeakFileOptions::addMSLevel(int level)
  {
    ms_levels_.push_back(level);
  }

  void PeakFileOptions::clearMSLevels()
  {
    ms_levels_.clear();
  }

  bool PeakFileOptions::hasMSLevels() const
  {
    return !ms_levels_.empty();
  }

  bool PeakFileOptions::containsMSLevel(int level) const
  {
    return find(ms_levels_.begin(), ms_levels_.end(), level) != ms_levels_.end();
  }

  const vector<Int>& PeakFileOptions::getMSLevels() const
  {
    return ms_levels_;
  }

  void PeakFileOptions::setCompression(bool compress)
  {
    zlib_compression_ = compress;
  }

  bool PeakFileOptions::getCompression() const
  {
    return zlib_compression_;
  }

  bool PeakFileOptions::getAlwaysAppendData() const
  {
    return always_append_data_;
  }

  void PeakFileOptions::setAlwaysAppendData(bool always_append_data)
  {
    always_append_data_ = always_append_data;
  }

  bool PeakFileOptions::getFillData() const
  {
    return fill_data_;
  }

  void PeakFileOptions::setSkipXMLChecks(bool skip)
  {
    skip_xml_checks_ = skip;
  }

  bool PeakFileOptions::getSkipXMLChecks() const
  {
    return skip_xml_checks_;
  }

  void PeakFileOptions::setSortSpectraByMZ(bool sort)
  {
    sort_spectra_by_mz_ = sort;
  }

  bool PeakFileOptions::getSortSpectraByMZ() const
  {
    return sort_spectra_by_mz_;
  }

  void PeakFileOptions::setSortChromatogramsByRT(bool sort)
  {
    sort_chromatograms_by_rt_ = sort;
  }

  bool PeakFileOptions::getSortChromatogramsByRT() const
  {
    return sort_chromatograms_by_rt_;
  }

  void PeakFileOptions::setFillData(bool fill_data)
  {
    fill_data_ = fill_data;
  }

  void PeakFileOptions::setMz32Bit(bool mz_32_bit)
  {
    mz_32_bit_ = mz_32_bit;
  }

  bool PeakFileOptions::getMz32Bit() const
  {
    return mz_32_bit_;
  }

  void PeakFileOptions::setIntensity32Bit(bool int_32_bit)
  {
    int_32_bit_ = int_32_bit;
  }

  bool PeakFileOptions::getIntensity32Bit() const
  {
    return int_32_bit_;
  }

  bool PeakFileOptions::getWriteIndex() const
  {
    return write_index_;
  }

  void PeakFileOptions::setWriteIndex(bool write_index)
  {
    write_index_ = write_index;
  }

  MSNumpressCoder::NumpressConfig PeakFileOptions::getNumpressConfigurationMassTime() const
  {
    return np_config_mz_;
  }

  void PeakFileOptions::setNumpressConfigurationMassTime(MSNumpressCoder::NumpressConfig config)
  {
    if (config.np_compression == MSNumpressCoder::SLOF || config.np_compression == MSNumpressCoder::PIC)
    {
      std::cerr << "Warning, compression of m/z or time dimension with pic or slof algorithms can lead to data loss" << std::endl;
    }
    np_config_mz_ = config;
  }

  MSNumpressCoder::NumpressConfig PeakFileOptions::getNumpressConfigurationIntensity() const
  {
    return np_config_int_;
  }

  void PeakFileOptions::setNumpressConfigurationIntensity(MSNumpressCoder::NumpressConfig config)
  {
    np_config_int_ = config;
  }

  MSNumpressCoder::NumpressConfig PeakFileOptions::getNumpressConfigurationFloatDataArray() const
  {
    return np_config_fda_;
  }

  void PeakFileOptions::setNumpressConfigurationFloatDataArray(MSNumpressCoder::NumpressConfig config)
  {
    np_config_fda_ = config;
  }

  Size PeakFileOptions::getMaxDataPoolSize() const
  {
    return maximal_data_pool_size_;
  }

  void PeakFileOptions::setMaxDataPoolSize(Size size)
  {
    maximal_data_pool_size_ = size;
  }

  bool PeakFileOptions::getPrecursorMZSelectedIon() const
  {
    return precursor_mz_selected_ion_;
  }

  void PeakFileOptions::setPrecursorMZSelectedIon(bool choice)
  {
    precursor_mz_selected_ion_ = choice;
  }

  bool PeakFileOptions::hasFilters() const
  {
    return (has_rt_range_ || hasMSLevels());
  }

} // namespace OpenMS
