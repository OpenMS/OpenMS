// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow 
// $Authors: Chris Bielow
// --------------------------------------------------------------------------
#pragma once
#include <OpenMS/SYSTEM/SIMDe.h>
/**
   @brief Wrapper around simde__m128i for better readability and portability
   
   Provides overloaded operators and utility functions for common SIMD operations.
 */
struct Vec128i {
  simde__m128i val;

  Vec128i() = default;
  Vec128i(simde__m128i v) : val(v) {}
  Vec128i(const Vec128i& other) = default;

  // Conversion to simde__m128i
  operator simde__m128i() const { return val; }

  // Bitwise AND
  Vec128i operator&(const Vec128i& other) const {
    return Vec128i(simde_mm_and_si128(val, other.val));
  }

  Vec128i& operator&=(const Vec128i& other) {
    val = simde_mm_and_si128(val, other.val);
    return *this;
  }

  // Bitwise OR
  Vec128i operator|(const Vec128i& other) const {
    return Vec128i(simde_mm_or_si128(val, other.val));
  }

  Vec128i& operator|=(const Vec128i& other) {
    val = simde_mm_or_si128(val, other.val);
    return *this;
  }

  // Bitwise XOR
  Vec128i operator^(const Vec128i& other) const {
    return Vec128i(simde_mm_xor_si128(val, other.val));
  }

  Vec128i& operator^=(const Vec128i& other) {
    val = simde_mm_xor_si128(val, other.val);
    return *this;
  }

  /// Arithmetic addition on 8-bit lanes (for byte-wise operations)
  Vec128i operator+(const Vec128i& other) const {
    return Vec128i(simde_mm_add_epi8(val, other.val));
  }

  /// Arithmetic subtraction on 8-bit lanes (for byte-wise operations)
  Vec128i operator-(const Vec128i& other) const {
    return Vec128i(simde_mm_sub_epi8(val, other.val));
  }
  /// Left shift by constant on 32-bit lanes (for Base64 bit manipulation)
  Vec128i slli_epi32(int imm) const {
    return Vec128i(simde_mm_slli_epi32(val, imm));
  }

  /// Right shift by constant on 32-bit lanes (for Base64 bit manipulation)
  Vec128i srli_epi32(int imm) const {
    return Vec128i(simde_mm_srli_epi32(val, imm));
  }

  /// Shuffle bytes according to mask (requires SSSE3 support)
  /// @param mask Control mask for byte shuffling
  /// @return New Vec128i with shuffled bytes
  Vec128i shuffle_epi8(const Vec128i& mask) const {
    return Vec128i(simde_mm_shuffle_epi8(val, mask.val));
  }

  /// Compare 8-bit lanes for less-than condition
  /// @param a First operand
  /// @param b Second operand  
  /// @return Comparison mask (0xFF for true lanes, 0x00 for false)
  static Vec128i cmplt_epi8(const Vec128i& a, const Vec128i& b) {
    return Vec128i(simde_mm_cmplt_epi8(a.val, b.val));
  }

  /// Compare 8-bit lanes for equality
  /// @param a First operand
  /// @param b Second operand
  /// @return Comparison mask (0xFF for equal lanes, 0x00 for unequal)
  static Vec128i cmpeq_epi8(const Vec128i& a, const Vec128i& b) {
    return Vec128i(simde_mm_cmpeq_epi8(a.val, b.val));
  }

  /// Bitwise ANDNOT operation (a AND NOT b)
  /// @param a First operand
  /// @param b Second operand (to be inverted before AND)
  /// @return Result of (a AND NOT b)
  static Vec128i andnot(const Vec128i& a, const Vec128i& b) {
    return Vec128i(simde_mm_andnot_si128(a.val, b.val));
  }
};
