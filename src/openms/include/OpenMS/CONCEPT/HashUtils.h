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

namespace OpenMS
{
  /**
   * @brief Portable hash utilities for OpenMS classes.
   *
   * This header provides hash functions that are:
   * - Portable: Same input produces same hash across platforms/compilers
   * - Fast: Uses efficient FNV-1a algorithm
   * - Low collision: Uses proper hash combining with golden ratio mixing
   *
   * @note std::hash<std::string> is implementation-defined and NOT portable.
   *       Use fnv1a_hash_string() or fnv1a_hash_bytes() for portable hashing.
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
   * @return Portable 64-bit hash value
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
   * Portable alternative to std::hash<std::string>.
   *
   * @param s String to hash
   * @return Portable hash value
   */
  inline std::size_t fnv1a_hash_string(const std::string& s) noexcept
  {
    return fnv1a_hash_bytes(s.data(), s.size());
  }

  /**
   * @brief Combine a hash value with additional data using golden ratio mixing.
   *
   * This function combines two hash values using the well-known formula from Boost:
   *   seed ^= hash + 0x9e3779b9 + (seed << 6) + (seed >> 2)
   *
   * The constant 0x9e3779b9 is derived from the golden ratio and provides
   * good bit mixing properties.
   *
   * @param seed The existing hash value (modified in place)
   * @param value The new hash value to combine
   */
  inline void hash_combine(std::size_t& seed, std::size_t value) noexcept
  {
    seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  /**
   * @brief Portable hash for an integer type.
   *
   * Uses bit mixing to ensure good distribution for integer values.
   *
   * @param value Integer value to hash
   * @return Hash value
   */
  template<typename T>
  inline std::size_t hash_int(T value) noexcept
  {
    // Use FNV-1a on the bytes of the integer for portability
    return fnv1a_hash_bytes(&value, sizeof(T));
  }

  /**
   * @brief Portable hash for a character.
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

} // namespace OpenMS
