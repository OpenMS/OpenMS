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
 * @brief Wrapper around simde__m128i for better readability and portability
 * 
 * Provides overloaded operators and utility functions for common SIMD operations.
 */
struct Vec128i {
  simde__m128i val;

  Vec128i() = default;
  Vec128i(simde__m128i v) : val(v) {}
  Vec128i(const Vec128i& other) = default;
  // simde__m128i& raw() { return val; }
  // const simde__m128i& raw() const { return val; }

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

    // Arithmetic Add
  Vec128i operator+(const Vec128i& other) const {
    return Vec128i(simde_mm_add_epi8(val, other.val));
  }

  // Arithmetic Sub
  Vec128i operator-(const Vec128i& other) const {
    return Vec128i(simde_mm_sub_epi8(val, other.val));
  }
    // Left shift by constant
  Vec128i slli_epi32(int imm) const {
    return Vec128i(simde_mm_slli_epi32(val, imm));
  }

  // Right shift by constant
  Vec128i srli_epi32(int imm) const {
    return Vec128i(simde_mm_srli_epi32(val, imm));
  }

  // Shuffle (requires SSSE3)
  Vec128i shuffle_epi8(const Vec128i& mask) const {
    return Vec128i(simde_mm_shuffle_epi8(val, mask.val));
  }

  // Comparison (less than, returns mask)
  static Vec128i cmplt_epi8(const Vec128i& a, const Vec128i& b) {
    return Vec128i(simde_mm_cmplt_epi8(a.val, b.val));
  }

  // Equality
  static Vec128i cmpeq_epi8(const Vec128i& a, const Vec128i& b) {
    return Vec128i(simde_mm_cmpeq_epi8(a.val, b.val));
  }
   // Bitwise ANDNOT
  static Vec128i andnot(const Vec128i& a, const Vec128i& b) {
    return Vec128i(simde_mm_andnot_si128(a.val, b.val));
  }
  // Optional: shift and other operators can be added similarly
};
