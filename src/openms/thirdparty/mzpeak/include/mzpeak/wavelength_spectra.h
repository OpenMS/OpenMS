/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <optional>
#include <string>

#include "mzpeak/util/enumerable_proxy.h"
#include "mzpeak/wavelength_spectrum.h"

namespace MzPeak::Data {
class Signals;
}

namespace MzPeak {

/**
 * Access all wavelength spectra in a MzPeak file.
 */
class WavelengthSpectra final : public Util::EnumerableProxy<WavelengthSpectrum> {
public:
  /// Default constructor — empty (no wavelength spectrum table present).
  WavelengthSpectra();

  /// Destructor.
  ~WavelengthSpectra() = default;

  /// Low-level constructor.
  /// @p count is ignored when zero and the signal table has rows; see
  /// Chromatograms for why a declared count is not trusted on its own.
  explicit WavelengthSpectra(
      std::unique_ptr<Data::Signals>,
      std::optional<std::size_t> count = std::nullopt,
      std::map<uint64_t, WavelengthSpectrumMetadata> metadata = {});

  /// Resolve a native id to its index, or nullopt when unknown.
  std::optional<std::size_t> index_for_id(const std::string& id) const;

  /// Fetch the spectrum with the given native id.
  /// @throws ParquetError when the id is unknown.
  WavelengthSpectrum by_id(const std::string& id) const;

private:
  std::shared_ptr<Data::Signals> data_;

  // Shared into every WavelengthSpectrum this collection hands out.
  std::shared_ptr<const std::map<uint64_t, WavelengthSpectrumMetadata>> md_map_;

  // Native id -> index, built from md_map_ at construction.
  std::map<std::string, std::size_t> id_to_index_;

  // Const for the same reason as Chromatograms::fetch.
  WavelengthSpectrum fetch_(uint64_t) const;

  /// EnumerableProxy's fetch is a non-const pure virtual; ours is const so the
  /// const public API (by_id, range queries) still works.  Forward.
  WavelengthSpectrum fetch(uint64_t n) override { return fetch_(n); }
};

} // namespace MzPeak
