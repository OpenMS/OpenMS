/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "mzpeak/data/array_index.h"
#include "mzpeak/data/encoding.h"
#include "mzpeak/data/signals.h"
#include "mzpeak/util/slice.h"
#include "mzpeak/wavelength_spectrum_metadata.h"

namespace MzPeak {

// Forward declaration.
class WavelengthSpectra;

/**
 * Access to a single wavelength spectrum in an MzPeak file.
 *
 * A wavelength spectrum is structurally identical to a chromatogram: a point
 * table holding a primary axis array (the wavelength array, MS:1000617) and an
 * intensity array (MS:1000515).
 */
class WavelengthSpectrum final {
public:
  using wavelength_type = double;
  using intensity_type = float;

  // Special members are implicit (rule of zero): copyable AND movable.

  /**
   * Wavelength values (the MS:1000617 wavelength array).
   */
  const std::vector<wavelength_type>& wavelength() const;

  /**
   * Intensity values.
   */
  const std::vector<intensity_type>& intensity() const;

  /**
   * Unit CURIE of the values returned by @ref wavelength (UO:0000018 is
   * nanometres), or empty when the file does not say.
   */
  const std::string& wavelength_unit() const;

  /**
   * Unit CURIE of the values returned by @ref intensity, or empty when the file
   * does not say.  UV/Vis intensities are frequently absorbance units rather
   * than detector counts, and the two are not interchangeable.
   */
  const std::string& intensity_unit() const;

  /**
   * Descriptive metadata for this spectrum (id, acquisition time, observed
   * wavelength range).  Returns an empty record when the file carries no
   * wavelength metadata table.
   */
  const WavelengthSpectrumMetadata& metadata() const;

protected:
  friend class WavelengthSpectra;

  WavelengthSpectrum(
      uint64_t index,
      std::shared_ptr<Data::Signals>,
      const std::vector<Data::ArrayIndex::Dimension>&,
      std::unique_ptr<Util::Slice>,
      std::shared_ptr<const std::map<uint64_t, WavelengthSpectrumMetadata>> = {});

private:
  using decoder_type = Data::Encoding::Decoder<double>;

  uint64_t index_;
  decoder_type decoder_;
  std::vector<wavelength_type> wavelength_;
  std::vector<intensity_type> intensity_;
  std::string wavelength_unit_;
  std::string intensity_unit_;
  std::shared_ptr<const std::map<uint64_t, WavelengthSpectrumMetadata>> md_map_;
};

} // namespace MzPeak
