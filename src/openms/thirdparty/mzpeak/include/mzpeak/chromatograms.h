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

#include "mzpeak/chromatogram.h"
#include "mzpeak/util/enumerable_proxy.h"

namespace MzPeak::Data {
class Signals;
}

namespace MzPeak {

/**
 * Access all chromatograms in a MzPeak file.
 */
class Chromatograms final : public Util::EnumerableProxy<Chromatogram> {
public:
  /// Default constructor — empty (no chromatogram table present).
  Chromatograms();

  /// Destructor.
  ~Chromatograms() = default;

  /// Low-level constructor.
  ///
  /// @p count is the entity count the file declares.  It is IGNORED when it is
  /// zero and the signal table has rows: the reference writer emits zero for
  /// several count keys, and trusting one would report an empty run rather
  /// than raising anything.
  explicit Chromatograms(std::unique_ptr<Data::Signals>,
                         std::optional<std::size_t> count = std::nullopt,
                         std::map<uint64_t, ChromatogramMetadata> metadata = {});

  /**
   * Resolve a native chromatogram id to its index, or nullopt when unknown.
   * The FIRST chromatogram carrying an id wins; ids are not guaranteed unique.
   */
  std::optional<std::size_t> index_for_id(const std::string& id) const;

  /**
   * Fetch the chromatogram with the given native id.
   * @throws ParquetError when the id is unknown.
   */
  Chromatogram by_id(const std::string& id) const;

private:
  std::shared_ptr<Data::Signals> data_;

  // Shared into every Chromatogram this collection hands out.
  std::shared_ptr<const std::map<uint64_t, ChromatogramMetadata>> md_map_;

  // Native id -> index, built from md_map_ at construction.
  std::map<std::string, std::size_t> id_to_index_;

  // Const because it mutates no Chromatograms state, which lets by_id be
  // const without casting the constness away.
  Chromatogram fetch_(uint64_t) const;

  /// EnumerableProxy's fetch is a non-const pure virtual; ours is const so the
  /// const public API (by_id, range queries) still works.  Forward.
  Chromatogram fetch(uint64_t n) override { return fetch_(n); }
};

} // namespace MzPeak
