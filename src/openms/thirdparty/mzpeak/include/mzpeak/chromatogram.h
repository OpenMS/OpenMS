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

#include "mzpeak/chromatogram_metadata.h"
#include "mzpeak/data/array_index.h"
#include "mzpeak/data/encoding.h"
#include "mzpeak/data/signals.h"
#include "mzpeak/util/slice.h"

namespace MzPeak {

// Forward declaration.
class Chromatograms;

/**
 * Access to a single chromatogram in an MzPeak file.
 */
class Chromatogram final {
public:
  using time_type = double;
  using intensity_type = float;

  // Special members are implicit (rule of zero): copyable AND movable.

  /**
   * Time values, in SECONDS.
   *
   * mzPeak stores chromatogram time in minutes (UO:0000031) in every file seen
   * so far, and the specification only recommends rather than requires it.  The
   * stored unit is read from the array index and converted here, so this agrees
   * with EicPoint::time and SpectrumMetadata::retention_time, which are also
   * seconds.  Returning the raw stored number instead made the same library
   * report two chromatogram time bases differing by 60x.
   *
   * @throws ParquetError when the file declares a time unit that is neither
   * minutes nor seconds -- converting it would be a guess.
   */
  const std::vector<time_type>& time() const;

  /**
   * Intensity values.
   *
   * @note A file may carry several intensity arrays in different units --
   * detector counts and absorbance units both appear in the bundled `has_uv`
   * fixture -- and this returns whichever one belongs to this chromatogram.
   * Ask @ref intensity_unit which it was.
   */
  const std::vector<intensity_type>& intensity() const;

  /**
   * Unit CURIE of the values returned by @ref intensity, or empty when the file
   * does not say.
   *
   * A file may store several intensity arrays in DIFFERENT units and give each
   * chromatogram one of them -- the bundled `has_uv` fixture holds a TIC in
   * detector counts (MS:1000131) and a DAD trace in absorbance units
   * (UO:0000269).  Both decode correctly; only this tells them apart.
   */
  const std::string& intensity_unit() const;

  /**
   * Descriptive metadata for this chromatogram (id, type, polarity,
   * precursors).  Returns an empty record when the file carries no
   * chromatogram metadata table.
   */
  const ChromatogramMetadata& metadata() const;

protected:
  friend class Chromatograms;

  Chromatogram(uint64_t index,
               std::shared_ptr<Data::Signals>,
               const std::vector<Data::ArrayIndex::Dimension>&,
               std::unique_ptr<Util::Slice>,
               std::shared_ptr<const std::map<uint64_t, ChromatogramMetadata>> = {});

private:
  using decoder_type = Data::Encoding::Decoder<double>;

  uint64_t index_;
  decoder_type decoder_;
  std::vector<time_type> time_;
  std::vector<intensity_type> intensity_;
  std::string intensity_unit_;
  std::shared_ptr<const std::map<uint64_t, ChromatogramMetadata>> md_map_;
};

} // namespace MzPeak
