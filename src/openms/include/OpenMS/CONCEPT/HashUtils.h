// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <type_traits>

namespace OpenMS
{
  /**
   * @brief Hash utilities for OpenMS classes.
   *
   * This header provides hash functions that are:
   * - Consistent: Equal inputs produce equal hashes
   * - Fast: Uses efficient FNV-1a algorithm
   * - Low collision: Uses proper hash combining with golden ratio mixing
   *
   * @note These hash functions satisfy the requirements for std::unordered_map/set:
   *       if a == b, then hash(a) == hash(b). When used with stable identifiers
   *       (strings, not pointers), hashes are reproducible across process runs
   *       on the same platform. Not suitable for cross-platform persistent storage
   *       due to endianness and size_t differences.
   *
   * @ingroup Concept
   */

  /**
   * @brief FNV-1a hash for a byte sequence.
   *
   * FNV-1a is a non-cryptographic hash function that is:
   * - Fast: Only XOR and multiply per byte
   * - Simple: Easy to implement and verify
   * - Good distribution: Well-tested for hash table use
   *
   * @param data Pointer to data bytes
   * @param size Number of bytes
   * @return Hash value (size_t)
   */
  inline std::size_t fnv1a_hash_bytes(const void* data, std::size_t size) noexcept
  {
    const auto* bytes = static_cast<const unsigned char*>(data);
    // FNV-1a 64-bit constants
    std::size_t hash = 14695981039346656037ull; // FNV offset basis
    for (std::size_t i = 0; i < size; ++i)
    {
      hash ^= bytes[i];
      hash *= 1099511628211ull; // FNV prime
    }
    return hash;
  }

  /**
   * @brief FNV-1a hash for a string.
   *
   * Consistent alternative to std::hash<std::string> which may vary
   * between standard library implementations.
   *
   * @param s String to hash
   * @return Hash value
   */
  inline std::size_t fnv1a_hash_string(const std::string& s) noexcept
  {
    return fnv1a_hash_bytes(s.data(), s.size());
  }

  /**
   * @brief Combine a hash value with additional data using golden ratio mixing.
   *
   * This function combines two hash values using the well-known formula from Boost:
   *   seed ^= hash + golden_ratio + (seed << 6) + (seed >> 2)
   *
   * The constant 0x9e3779b97f4a7c15 is derived from the 64-bit golden ratio
   * (2^64 / phi) and provides good bit mixing properties on 64-bit systems.
   *
   * @param seed The existing hash value (modified in place)
   * @param value The new hash value to combine
   */
  inline void hash_combine(std::size_t& seed, std::size_t value) noexcept
  {
    seed ^= value + 0x9e3779b97f4a7c15 + (seed << 6) + (seed >> 2);
  }

  /**
   * @brief Hash for an integer type.
   *
   * Uses FNV-1a on the bytes of the integer for good distribution.
   *
   * @note This function hashes the platform's native byte representation.
   *       Results may differ on big-endian vs little-endian systems.
   *       This is acceptable for in-process hash containers but not for
   *       cross-platform persistent storage. See header documentation.
   *
   * @tparam T Must be an integral type (int, long, size_t, etc.)
   * @param[in] value Integer value to hash
   * @return Hash value
   */
  template<typename T>
  inline std::size_t hash_int(T value) noexcept
  {
    static_assert(std::is_integral_v<T>, "hash_int requires an integral type");
    return fnv1a_hash_bytes(&value, sizeof(T));
  }

  /**
   * @brief Hash for a character.
   *
   * @param c Character to hash
   * @return Hash value (same as fnv1a of single byte)
   */
  inline std::size_t hash_char(char c) noexcept
  {
    // Apply FNV-1a algorithm for single character
    std::size_t hash = 14695981039346656037ull;
    hash ^= static_cast<unsigned char>(c);
    hash *= 1099511628211ull;
    return hash;
  }

  /**
   * @brief Hash for a floating point type (float or double).
   *
   * Hashes the raw bytes of the floating point value. This is consistent
   * with the default operator== for floating point types.
   *
   * @note Normalizes -0.0 to +0.0 since they compare equal but have different
   *       bit patterns. This ensures the hash contract is satisfied.
   *
   * @tparam T Must be a floating point type (float, double)
   * @param value Floating point value to hash
   * @return Hash value
   */
  template<typename T>
  inline std::size_t hash_float(T value) noexcept
  {
    static_assert(std::is_floating_point_v<T>, "hash_float requires a floating point type");
    // Normalize -0.0 to +0.0 (they compare equal but have different bit patterns)
    if (value == T(0)) value = T(0);
    return fnv1a_hash_bytes(&value, sizeof(T));
  }

} // namespace OpenMS
