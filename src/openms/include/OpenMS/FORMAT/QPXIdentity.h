// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/Types.h>

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace OpenMS
{

/**
  @brief Opaque QPX row identities (@c feature_id, @c psm_id, @c pg_id)

  Every QPX view carries a mandatory @c int64 identity column that is its primary key. The
  value is opaque -- it is not meant to be parsed or reversed, only compared -- and is derived
  deterministically from a footer-declared @c identity_composite of ordinary columns.

  <b>These functions reproduce @c qpx.core.data.identity byte for byte.</b> That is a hard
  requirement, not a nicety. @c qpxc re-derives @c feature_id for identified rows from the
  declared composite, but never rewrites the @c psm.feature_id that points at it. An OpenMS
  collection whose ids were invented rather than derived would therefore turn into one with
  dangling cross-references the moment it was converted -- which qpx's own dataset validation
  reports as @c dangling_feature_id.

  The derivation is BLAKE2b truncated to 8 bytes over the canonical JSON encoding of the
  composite, read as a big-endian signed 64-bit integer. Negative values are normal and carry
  no meaning.

  @note Identity is meaningful <b>within a file only</b>. Two QPX files must not be joined on
        @c feature_id alone, and a test reference must not pin id values across files: the
        feature composite contains @c rt and @c observed_mz, so one ULP of platform drift
        changes the whole id rather than one digit of it.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
namespace QPXIdentity
{
  /**
    @brief The feature&harr;PSM edge of one QPX collection, as @c psm_id &rarr; @c feature_id

    Collected by ConsensusMapArrowExport::exportToArrow() / exportToParquet() /
    exportToParquetStreaming() through their @c out_links parameter, then handed to the psm
    exporter, which fills @c psm.feature_id from it. Keying on the <i>id</i> rather than on a
    pointer or a row index is what keeps the two directions consistent: both views derive the key
    from the same four persisted columns, so they agree by construction rather than by both
    happening to walk the identifications in the same order.
  */
  using FeatureLinks = std::unordered_map<Int64, Int64>;

  /// A composite component that is JSON @c null (a null column value)
  struct Null {};

  /**
    @brief A @c float32 composite component

    Distinct from @c double so a caller cannot pass an unrounded value by accident: QPX hashes
    the value <i>as persisted</i>, so a column stored as @c float32 must be narrowed to
    @c float before it is encoded, whatever precision it was computed at.
  */
  struct Float32
  {
    float value;
  };

  /// One component of an identity composite, in its persisted Arrow type
  using Value = std::variant<Null,                      ///< JSON @c null
                             std::string,               ///< JSON string
                             Int64,                     ///< JSON integer
                             Float32,                   ///< JSON number
                             std::vector<Int32>,        ///< JSON array of integers (@c scan)
                             std::vector<std::string>>; ///< JSON array of strings (@c grouped_runs)

  /**
    @brief Render one @c float32 the way Python's @c repr writes it inside JSON

    QPX encodes composites with @c json.dumps, which formats floats with @c float.__repr__:
    the shortest decimal that round-trips, switching to exponential notation when the decimal
    point would fall at position <tt>&lt;= -4</tt> or <tt>&gt; 16</tt>, and appending @c ".0"
    to a value that would otherwise look like an integer. C++'s own @c general format picks
    between fixed and scientific by a different rule, so it cannot be used directly.

    @param[in] value The value as persisted; it is widened to @c double first, exactly as
               Arrow's @c cast + @c to_pylist does on the reading side
    @return e.g. @c "0.1", @c "1000000000000000.0", @c "1e+16", @c "9.999e-05", @c "-0.0"
  */
  OPENMS_DLLAPI std::string formatFloat(float value);

  /**
    @brief Canonical JSON encoding of an identity composite

    Compact separators, ASCII-escaped strings -- the exact bytes qpx hashes.

    @param[in] values The composite components, in the order the composite declares
    @param[in] unordered_list_indices Positions whose list value is a <i>set</i> and is sorted
               before encoding. QPX uses this for @c grouped_runs; ordered lists such as
               @c scan keep their order.
    @return The encoded composite, e.g. <tt>["BSA1_F1",[2311],"PEPTIDER",2]</tt>
  */
  OPENMS_DLLAPI std::string canonical(const std::vector<Value>& values,
                                      const std::vector<size_t>& unordered_list_indices = {});

  /**
    @brief Derive an opaque signed 64-bit identity from a composite

    @param[in] values The composite components, in the order the composite declares
    @param[in] unordered_list_indices Set-valued positions, see canonical()
    @return The identity; uniqueness is not guaranteed by construction but by the primary-key
            uniqueness check, which reports a collision rather than losing a row to it
  */
  OPENMS_DLLAPI Int64 deriveId(const std::vector<Value>& values,
                               const std::vector<size_t>& unordered_list_indices = {});

  // -------------------------------------------------------------------------------------
  // Per-view identities
  //
  // The composite and its ORDER are declared once, here, and every exporter goes through
  // these functions. The order is part of the hash, so a call site assembling its own vector
  // would be one argument swap away from ids that no other reader can reproduce.
  //
  // These are the composites qpxc hard-codes for OpenMS `-out_qpx` input
  // (qpx/converters/openms/converter.py), so a conversion re-derives the same values.
  // -------------------------------------------------------------------------------------

  /// Footer @c identity_composite for the feature view
  OPENMS_DLLAPI extern const char* const FEATURE_COMPOSITE;
  /// Footer @c identity_composite for the psm view
  OPENMS_DLLAPI extern const char* const PSM_COMPOSITE;
  /// Footer @c identity_composite for the pg view
  OPENMS_DLLAPI extern const char* const PG_COMPOSITE;

  /**
    @brief @c feature_id from a feature row's persisted column values

    Composite: <tt>(run_file_name, peptidoform, charge, rt, scan, observed_mz)</tt>. The floats
    are what make an <i>unidentified</i> feature row identifiable at all -- such a row has an
    empty peptidoform and an empty scan, and its run and charge are shared with dozens of
    others.

    @param[in] run_file_name Bare run name (no directory, no extension)
    @param[in] peptidoform ProForma notation; empty for an unidentified feature
    @param[in] charge Charge state as persisted (@c int16)
    @param[in] rt Retention time in seconds, or @c nullopt when the column is null
    @param[in] scan Scan components; empty for an unidentified feature
    @param[in] observed_mz Experimental m/z
  */
  OPENMS_DLLAPI Int64 featureId(const std::string& run_file_name,
                                const std::string& peptidoform,
                                Int64 charge,
                                std::optional<float> rt,
                                const std::vector<Int32>& scan,
                                float observed_mz);

  /**
    @brief @c psm_id from a psm row's persisted column values

    Composite: <tt>(run_file_name, scan, peptidoform, charge)</tt> -- no floats, so a psm
    identity is stable across platforms.

    @param[in] run_file_name Bare run name (no directory, no extension)
    @param[in] scan Scan components, in order
    @param[in] peptidoform ProForma notation
    @param[in] charge Precursor charge as persisted (@c int16)
  */
  OPENMS_DLLAPI Int64 psmId(const std::string& run_file_name,
                            const std::vector<Int32>& scan,
                            const std::string& peptidoform,
                            Int64 charge);

  /**
    @brief @c pg_id from a protein-group row's persisted column values

    Composite: <tt>(anchor_protein, grouped_runs, label)</tt>, with @c grouped_runs sorted --
    it is the set of raw files aggregated into one quantity, and its order is not part of the
    group's identity.

    @param[in] anchor_protein Leading protein of the group
    @param[in] grouped_runs Raw files of this quantification unit, in any order
    @param[in] label Channel label, or @c nullopt for an identification-only group
  */
  OPENMS_DLLAPI Int64 pgId(const std::string& anchor_protein,
                           const std::vector<std::string>& grouped_runs,
                           const std::optional<std::string>& label);

} // namespace QPXIdentity

} // namespace OpenMS
