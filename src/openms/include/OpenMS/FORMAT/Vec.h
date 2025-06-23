#pragma once
#include <simde/x86/sse2.h>  // or the lowest needed header

class Vec128 {
public:
  simde__m128i val;

  Vec128() = default;
  Vec128(simde__m128i v) : val(v) {}
  Vec128(const Vec128& other) = default;
  simde__m128i& raw() { return val; }
  const simde__m128i& raw() const { return val; }

  // Conversion to simde__m128i
  operator simde__m128i() const { return val; }

  // Bitwise AND
  Vec128 operator&(const Vec128& other) const {
    return Vec128(simde_mm_and_si128(val, other.val));
  }

  Vec128& operator&=(const Vec128& other) {
    val = simde_mm_and_si128(val, other.val);
    return *this;
  }

  // Bitwise OR
  Vec128 operator|(const Vec128& other) const {
    return Vec128(simde_mm_or_si128(val, other.val));
  }

  Vec128& operator|=(const Vec128& other) {
    val = simde_mm_or_si128(val, other.val);
    return *this;
  }

  // Bitwise XOR
  Vec128 operator^(const Vec128& other) const {
    return Vec128(simde_mm_xor_si128(val, other.val));
  }

  Vec128& operator^=(const Vec128& other) {
    val = simde_mm_xor_si128(val, other.val);
    return *this;
  }

    // Arithmetic Add
  Vec128 operator+(const Vec128& other) const {
    return Vec128(simde_mm_add_epi8(val, other.val));
  }

  // Arithmetic Sub
  Vec128 operator-(const Vec128& other) const {
    return Vec128(simde_mm_sub_epi8(val, other.val));
  }
    // Left shift by constant
  Vec128 slli_epi32(int imm) const {
    return Vec128(simde_mm_slli_epi32(val, imm));
  }

  // Right shift by constant
  Vec128 srli_epi32(int imm) const {
    return Vec128(simde_mm_srli_epi32(val, imm));
  }

  // Shuffle (requires SSSE3)
  Vec128 shuffle_epi8(const Vec128& mask) const {
    return Vec128(simde_mm_shuffle_epi8(val, mask.val));
  }

  // Comparison (less than, returns mask)
  static Vec128 cmplt_epi8(const Vec128& a, const Vec128& b) {
    return Vec128(simde_mm_cmplt_epi8(a.val, b.val));
  }

  // Equality
  static Vec128 cmpeq_epi8(const Vec128& a, const Vec128& b) {
    return Vec128(simde_mm_cmpeq_epi8(a.val, b.val));
  }

  // And-not
  static Vec128 andnot_si128(const Vec128& a, const Vec128& b) {
    return Vec128(simde_mm_andnot_si128(a.val, b.val));
  }
   // Bitwise ANDNOT
  static Vec128 andnot(const Vec128& a, const Vec128& b) {
    return Vec128(simde_mm_andnot_si128(a.val, b.val));
  }
  // Optional: shift and other operators can be added similarly
};
