// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

#include <algorithm>
#include <iostream>

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
    size_only_(false),
    always_append_data_(false),
    skip_xml_checks_(false),
    sort_spectra_by_mz_(true),
    sort_chromatograms_by_rt_(true),
    fill_data_(true),
    write_index_(true),
    np_config_mz_(),
    np_config_int_(),
    maximal_data_pool_size_(100)
  {
  }

  PeakFileOptions::PeakFileOptions(const PeakFileOptions& options) :
    metadata_only_(options.metadata_only_),
    force_maxquant_compatibility_(options.force_maxquant_compatibility_),
    force_tpp_compatibility_(options.force_tpp_compatibility_),
    write_supplemental_data_(options.write_supplemental_data_),
    has_rt_range_(options.has_rt_range_),
    has_mz_range_(options.has_mz_range_),
    has_intensity_range_(options.has_intensity_range_),
    mz_32_bit_(options.mz_32_bit_),
    int_32_bit_(options.int_32_bit_),
    rt_range_(options.rt_range_),
    mz_range_(options.mz_range_),
    intensity_range_(options.intensity_range_),
    ms_levels_(options.ms_levels_),
    zlib_compression_(options.zlib_compression_),
    size_only_(options.size_only_),
    always_append_data_(options.always_append_data_),
    skip_xml_checks_(options.skip_xml_checks_),
    sort_spectra_by_mz_(options.sort_spectra_by_mz_),
    sort_chromatograms_by_rt_(options.sort_chromatograms_by_rt_),
    fill_data_(options.fill_data_),
    write_index_(options.write_index_),
    np_config_mz_(options.np_config_mz_),
    np_config_int_(options.np_config_int_),
    maximal_data_pool_size_(options.maximal_data_pool_size_)
  {
  }

  PeakFileOptions::~PeakFileOptions()
  {
  }

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
    has_rt_range_ = true;
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

  bool PeakFileOptions::getSizeOnly() const
  {
    return size_only_;
  }

  void PeakFileOptions::setSizeOnly(bool size_only)
  {
    size_only_ = size_only;
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

  Size PeakFileOptions::getMaxDataPoolSize() const
  {
    return maximal_data_pool_size_;
  }

  void PeakFileOptions::setMaxDataPoolSize(Size size)
  {
    maximal_data_pool_size_ = size;
  }

} // namespace OpenMS
