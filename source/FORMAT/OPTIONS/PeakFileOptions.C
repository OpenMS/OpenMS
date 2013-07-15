// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{
  PeakFileOptions::PeakFileOptions() :
    metadata_only_(false),
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
    always_append_data_(false)
  {
  }

  PeakFileOptions::PeakFileOptions(const PeakFileOptions & options) :
    metadata_only_(options.metadata_only_),
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
    always_append_data_(options.always_append_data_)
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

  void PeakFileOptions::setWriteSupplementalData(bool write)
  {
    write_supplemental_data_ = write;
  }

  bool PeakFileOptions::getWriteSupplementalData() const
  {
    return write_supplemental_data_;
  }

  void PeakFileOptions::setRTRange(const DRange<1> & range)
  {
    rt_range_ = range;
    has_rt_range_ = true;
  }

  bool PeakFileOptions::hasRTRange() const
  {
    return has_rt_range_;
  }

  const DRange<1> & PeakFileOptions::getRTRange() const
  {
    return rt_range_;
  }

  void PeakFileOptions::setMZRange(const DRange<1> & range)
  {
    mz_range_ = range;
    has_mz_range_ = true;
  }

  bool PeakFileOptions::hasMZRange() const
  {
    return has_mz_range_;
  }

  const DRange<1> & PeakFileOptions::getMZRange() const
  {
    return mz_range_;
  }

  void PeakFileOptions::setIntensityRange(const DRange<1> & range)
  {
    intensity_range_ = range;
    has_intensity_range_ = true;
  }

  bool PeakFileOptions::hasIntensityRange() const
  {
    return has_intensity_range_;
  }

  const DRange<1> & PeakFileOptions::getIntensityRange() const
  {
    return intensity_range_;
  }

  void PeakFileOptions::setMSLevels(const vector<Int> & levels)
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

  const vector<Int> & PeakFileOptions::getMSLevels() const
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

} // namespace OpenMS
