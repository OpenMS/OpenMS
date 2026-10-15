// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/QPXIdentity.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstring>
#include <string_view>

namespace OpenMS
{
namespace QPXIdentity
{

namespace
{
  // ---------------------------------------------------------------------------------------
  // BLAKE2b (RFC 7693), unkeyed, digest length 8
  //
  // Vendored rather than pulled in as a dependency: the only thing OpenMS needs from it is a
  // single small digest of a few dozen bytes, and the algorithm is fixed by the QPX spec, so
  // there is nothing to keep up to date. Note that BLAKE2b-64 is NOT BLAKE2b-512 truncated --
  // the digest length goes into the parameter block and therefore into the initial state.
  // ---------------------------------------------------------------------------------------

  constexpr std::array<UInt64, 8> BLAKE2B_IV = {
    0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
    0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL, 0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL};

  constexpr Byte BLAKE2B_SIGMA[12][16] = {
    { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15},
    {14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3},
    {11,  8, 12,  0,  5,  2, 15, 13, 10, 14,  3,  6,  7,  1,  9,  4},
    { 7,  9,  3,  1, 13, 12, 11, 14,  2,  6,  5, 10,  4,  0, 15,  8},
    { 9,  0,  5,  7,  2,  4, 10, 15, 14,  1, 11, 12,  6,  8,  3, 13},
    { 2, 12,  6, 10,  0, 11,  8,  3,  4, 13,  7,  5, 15, 14,  1,  9},
    {12,  5,  1, 15, 14, 13,  4, 10,  0,  7,  6,  3,  9,  2,  8, 11},
    {13, 11,  7, 14, 12,  1,  3,  9,  5,  0, 15,  4,  8,  6,  2, 10},
    { 6, 15, 14,  9, 11,  3,  0,  8, 12,  2, 13,  7,  1,  4, 10,  5},
    {10,  2,  8,  4,  7,  6,  1,  5, 15, 11,  9, 14,  3, 12, 13,  0},
    { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15},
    {14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3}};

  inline UInt64 rotr64(UInt64 x, unsigned n)
  {
    return (x >> n) | (x << (64 - n));
  }

  inline void mix(UInt64& a, UInt64& b, UInt64& c, UInt64& d, UInt64 x, UInt64 y)
  {
    a = a + b + x;
    d = rotr64(d ^ a, 32);
    c = c + d;
    b = rotr64(b ^ c, 24);
    a = a + b + y;
    d = rotr64(d ^ a, 16);
    c = c + d;
    b = rotr64(b ^ c, 63);
  }

  /// Compress one 128-byte block; @p counter is the number of bytes hashed INCLUDING this block
  void compress(std::array<UInt64, 8>& h, const Byte block[128], UInt64 counter, bool last)
  {
    UInt64 m[16];
    for (int i = 0; i < 16; ++i)
    {
      // Little-endian load, spelled out so the result does not depend on host byte order.
      UInt64 word = 0;
      for (int b = 7; b >= 0; --b) { word = (word << 8) | block[i * 8 + b]; }
      m[i] = word;
    }

    UInt64 v[16];
    for (int i = 0; i < 8; ++i) { v[i] = h[i]; }
    for (int i = 0; i < 8; ++i) { v[8 + i] = BLAKE2B_IV[i]; }
    v[12] ^= counter;
    // v[13] ^= counter >> 64 -- the high half of the 128-bit counter, always 0 for the inputs
    // here (a composite is at most a few kilobytes).
    if (last) { v[14] = ~v[14]; }

    for (int r = 0; r < 12; ++r)
    {
      const Byte* s = BLAKE2B_SIGMA[r];
      mix(v[0], v[4], v[8],  v[12], m[s[0]],  m[s[1]]);
      mix(v[1], v[5], v[9],  v[13], m[s[2]],  m[s[3]]);
      mix(v[2], v[6], v[10], v[14], m[s[4]],  m[s[5]]);
      mix(v[3], v[7], v[11], v[15], m[s[6]],  m[s[7]]);
      mix(v[0], v[5], v[10], v[15], m[s[8]],  m[s[9]]);
      mix(v[1], v[6], v[11], v[12], m[s[10]], m[s[11]]);
      mix(v[2], v[7], v[8],  v[13], m[s[12]], m[s[13]]);
      mix(v[3], v[4], v[9],  v[14], m[s[14]], m[s[15]]);
    }

    for (int i = 0; i < 8; ++i) { h[i] ^= v[i] ^ v[8 + i]; }
  }

  /// BLAKE2b with digest length 8, returned as the signed big-endian reading of those 8 bytes
  Int64 blake2b64(const std::string& message)
  {
    std::array<UInt64, 8> h = BLAKE2B_IV;
    h[0] ^= 0x01010000ULL ^ 8ULL; // parameter block: no key (kk=0), digest length nn=8

    const Byte* data = reinterpret_cast<const Byte*>(message.data());
    const size_t length = message.size();

    // Every block but the last is compressed with last=false. The last one is zero-padded and
    // carries the total length; an empty message is one all-zero block with counter 0.
    size_t offset = 0;
    while (length - offset > 128)
    {
      offset += 128;
      compress(h, data + offset - 128, offset, false);
    }
    Byte final_block[128] = {};
    std::memcpy(final_block, data + offset, length - offset);
    compress(h, final_block, length, true);

    // The digest is the first 8 bytes of the little-endian serialization of the state, i.e.
    // h[0] little-endian. Python reads those bytes big-endian and signed, which is the byte
    // reversal of h[0] reinterpreted as two's complement.
    UInt64 be = 0;
    for (int i = 0; i < 8; ++i) { be = (be << 8) | static_cast<Byte>(h[0] >> (8 * i)); }
    Int64 result;
    std::memcpy(&result, &be, sizeof(result)); // two's complement reinterpret, no UB
    return result;
  }

  // ---------------------------------------------------------------------------------------
  // Canonical JSON (json.dumps(separators=(",", ":"), sort_keys=True, ensure_ascii=True))
  // ---------------------------------------------------------------------------------------

  void appendEscapedCodepoint(std::string& out, UInt32 cp)
  {
    static const char* const HEX = "0123456789abcdef"; // Python formats \uXXXX in lower case
    out += "\\u";
    out += HEX[(cp >> 12) & 0xF];
    out += HEX[(cp >> 8) & 0xF];
    out += HEX[(cp >> 4) & 0xF];
    out += HEX[cp & 0xF];
  }

  /// Decode one UTF-8 sequence at @p pos; false (and @p length untouched) if it is malformed
  bool decodeUtf8(const std::string& text, size_t pos, UInt32& cp, size_t& length)
  {
    const auto byte = [&](size_t i) { return static_cast<Byte>(text[i]); };
    const Byte lead = byte(pos);
    size_t needed;
    if (lead < 0x80)         { cp = lead;         needed = 0; }
    else if (lead < 0xC2)    { return false; }    // continuation byte, or an overlong 2-byte lead
    else if (lead < 0xE0)    { cp = lead & 0x1F;  needed = 1; }
    else if (lead < 0xF0)    { cp = lead & 0x0F;  needed = 2; }
    else if (lead < 0xF5)    { cp = lead & 0x07;  needed = 3; }
    else                     { return false; }

    if (pos + needed >= text.size()) { return false; }
    for (size_t i = 1; i <= needed; ++i)
    {
      const Byte continuation = byte(pos + i);
      if ((continuation & 0xC0) != 0x80) { return false; }
      cp = (cp << 6) | (continuation & 0x3F);
    }
    // Reject overlong forms and surrogates, which would encode differently than Python does.
    if ((needed == 2 && cp < 0x800) || (needed == 3 && cp < 0x10000)) { return false; }
    if (cp >= 0xD800 && cp <= 0xDFFF) { return false; }
    if (cp > 0x10FFFF) { return false; }
    length = needed + 1;
    return true;
  }

  void encodeString(const std::string& text, std::string& out)
  {
    out += '"';
    size_t i = 0;
    while (i < text.size())
    {
      const Byte c = static_cast<Byte>(text[i]);
      switch (c)
      {
        case '"':  out += "\\\""; ++i; continue;
        case '\\': out += "\\\\"; ++i; continue;
        case '\n': out += "\\n";  ++i; continue;
        case '\r': out += "\\r";  ++i; continue;
        case '\t': out += "\\t";  ++i; continue;
        case '\b': out += "\\b";  ++i; continue;
        case '\f': out += "\\f";  ++i; continue;
        default: break;
      }
      if (c >= 0x20 && c <= 0x7E) // printable ASCII passes through; everything else is escaped
      {
        out += static_cast<char>(c);
        ++i;
      }
      else if (c < 0x20)
      {
        appendEscapedCodepoint(out, c);
        ++i;
      }
      else
      {
        UInt32 cp = 0;
        size_t length = 0;
        if (decodeUtf8(text, i, cp, length))
        {
          if (cp <= 0xFFFF)
          {
            appendEscapedCodepoint(out, cp);
          }
          else // outside the BMP: a surrogate pair, exactly as Python's ensure_ascii writes it
          {
            const UInt32 rest = cp - 0x10000;
            appendEscapedCodepoint(out, 0xD800 + (rest >> 10));
            appendEscapedCodepoint(out, 0xDC00 + (rest & 0x3FF));
          }
          i += length;
        }
        else
        {
          // Arrow validates its utf8 columns, so this is unreachable from a well-formed table.
          // Escape the offending byte rather than substitute U+FFFD: an identity must stay a
          // function of the input, even when the input is broken.
          appendEscapedCodepoint(out, c);
          ++i;
        }
      }
    }
    out += '"';
  }

  void encodeValue(const Value& value, std::string& out)
  {
    std::visit(
      [&out](const auto& held)
      {
        using T = std::decay_t<decltype(held)>;
        if constexpr (std::is_same_v<T, Null>)
        {
          out += "null";
        }
        else if constexpr (std::is_same_v<T, std::string>)
        {
          encodeString(held, out);
        }
        else if constexpr (std::is_same_v<T, Int64>)
        {
          out += std::to_string(held);
        }
        else if constexpr (std::is_same_v<T, Float32>)
        {
          out += formatFloat(held.value);
        }
        else if constexpr (std::is_same_v<T, std::vector<Int32>>)
        {
          out += '[';
          for (size_t i = 0; i < held.size(); ++i)
          {
            if (i != 0) { out += ','; }
            out += std::to_string(held[i]);
          }
          out += ']';
        }
        else // std::vector<std::string>
        {
          out += '[';
          for (size_t i = 0; i < held.size(); ++i)
          {
            if (i != 0) { out += ','; }
            encodeString(held[i], out);
          }
          out += ']';
        }
      },
      value);
  }
} // namespace

const char* const FEATURE_COMPOSITE = "run_file_name,peptidoform,charge,rt,scan,observed_mz";
const char* const PSM_COMPOSITE     = "run_file_name,scan,peptidoform,charge";
const char* const PG_COMPOSITE      = "pg_accessions,grouped_runs,label";

std::string formatFloat(float value)
{
  const double magnitude_source = static_cast<double>(value);
  if (std::isnan(magnitude_source)) { return "NaN"; }
  if (std::isinf(magnitude_source)) { return magnitude_source < 0 ? "-Infinity" : "Infinity"; }

  std::string out;
  if (std::signbit(magnitude_source)) { out += '-'; } // -0.0 keeps its sign, as repr() does
  const double magnitude = std::abs(magnitude_source);
  if (magnitude == 0.0)
  {
    out += "0.0";
    return out;
  }

  // Shortest round-tripping digits, obtained in scientific form so the digit string and the
  // decimal exponent can be read off directly. std::chars_format::general would pick between
  // fixed and scientific by C's rule, not Python's, so it cannot be used here.
  char buffer[64];
  const auto conversion =
    std::to_chars(buffer, buffer + sizeof(buffer), magnitude, std::chars_format::scientific);
  const std::string_view scientific(buffer, static_cast<size_t>(conversion.ptr - buffer));
  const size_t exponent_marker = scientific.find('e');

  std::string digits(scientific.substr(0, exponent_marker));
  digits.erase(std::remove(digits.begin(), digits.end(), '.'), digits.end());
  // std::to_chars writes the exponent printf-style, with an explicit sign; std::from_chars does
  // not accept a leading '+', so step over it rather than silently parse nothing.
  size_t exponent_start = exponent_marker + 1;
  if (exponent_start < scientific.size() && scientific[exponent_start] == '+') { ++exponent_start; }
  int exponent = 0;
  std::from_chars(scientific.data() + exponent_start,
                  scientific.data() + scientific.size(), exponent);

  // decpt is where the decimal point falls relative to the digit string, i.e. the value is
  // 0.<digits> * 10^decpt. This is the quantity Python's format_float_short() switches on.
  const int decimal_point = exponent + 1;
  const int digit_count = static_cast<int>(digits.size());

  if (decimal_point <= -4 || decimal_point > 16)
  {
    out += digits[0];
    if (digit_count > 1)
    {
      out += '.';
      out.append(digits, 1, std::string::npos);
    }
    out += 'e';
    const int printed = decimal_point - 1;
    out += (printed < 0 ? '-' : '+');
    const std::string magnitude_digits = std::to_string(std::abs(printed));
    if (magnitude_digits.size() < 2) { out += '0'; } // Python pads the exponent to two digits
    out += magnitude_digits;
  }
  else if (decimal_point <= 0)
  {
    out += "0.";
    out.append(static_cast<size_t>(-decimal_point), '0');
    out += digits;
  }
  else if (decimal_point >= digit_count)
  {
    out += digits;
    out.append(static_cast<size_t>(decimal_point - digit_count), '0');
    out += ".0"; // Py_DTSF_ADD_DOT_0: a float never prints as a bare integer
  }
  else
  {
    out.append(digits, 0, static_cast<size_t>(decimal_point));
    out += '.';
    out.append(digits, static_cast<size_t>(decimal_point), std::string::npos);
  }
  return out;
}

std::string canonical(const std::vector<Value>& values,
                      const std::vector<size_t>& unordered_list_indices)
{
  std::string out;
  out += '[';
  for (size_t i = 0; i < values.size(); ++i)
  {
    if (i != 0) { out += ','; }

    const bool unordered = std::find(unordered_list_indices.begin(),
                                     unordered_list_indices.end(), i)
                           != unordered_list_indices.end();
    if (unordered && std::holds_alternative<std::vector<std::string>>(values[i]))
    {
      // A set-valued field. qpx builds a dict keyed by each element's JSON encoding and then
      // walks the sorted keys (qpx/core/data/identity.py), so equal elements collapse to one
      // and the survivors are ordered by their encoded form -- a set, not just a sorted list.
      // Sorting without deduplicating agrees with qpx on every list that happens to have no
      // repeat and disagrees on exactly the ones that do.
      const auto& raw = std::get<std::vector<std::string>>(values[i]);
      std::vector<std::pair<std::string, const std::string*>> keyed;
      keyed.reserve(raw.size());
      for (const auto& element : raw)
      {
        std::string encoded;
        encodeString(element, encoded);
        keyed.emplace_back(std::move(encoded), &element);
      }
      std::sort(keyed.begin(), keyed.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
      keyed.erase(std::unique(keyed.begin(), keyed.end(),
                              [](const auto& a, const auto& b) { return a.first == b.first; }),
                  keyed.end());

      std::vector<std::string> as_set;
      as_set.reserve(keyed.size());
      for (const auto& entry : keyed) { as_set.push_back(*entry.second); }
      encodeValue(Value{std::move(as_set)}, out);
    }
    else
    {
      encodeValue(values[i], out);
    }
  }
  out += ']';
  return out;
}

Int64 deriveId(const std::vector<Value>& values, const std::vector<size_t>& unordered_list_indices)
{
  return blake2b64(canonical(values, unordered_list_indices));
}

Int64 featureId(const std::string& run_file_name,
                const std::string& peptidoform,
                Int64 charge,
                std::optional<float> rt,
                const std::vector<Int32>& scan,
                float observed_mz)
{
  return deriveId({Value{run_file_name},
                   Value{peptidoform},
                   Value{charge},
                   rt.has_value() ? Value{Float32{*rt}} : Value{Null{}},
                   Value{scan},
                   Value{Float32{observed_mz}}});
}

Int64 psmId(const std::string& run_file_name,
            const std::vector<Int32>& scan,
            const std::string& peptidoform,
            Int64 charge)
{
  return deriveId({Value{run_file_name}, Value{scan}, Value{peptidoform}, Value{charge}});
}

Int64 pgId(const std::vector<std::string>& pg_accessions,
           const std::vector<std::string>& grouped_runs,
           const std::optional<std::string>& label)
{
  // Both lists are sets. qpx picks the unordered positions by field name -- `grouped_runs` and
  // `pg_accessions` -- rather than by index (writers/base.py), so the two happen to be positions
  // 0 and 1 of this composite.
  return deriveId({Value{pg_accessions},
                   Value{grouped_runs},
                   label.has_value() ? Value{*label} : Value{Null{}}},
                  /*unordered_list_indices=*/{0, 1});
}

} // namespace QPXIdentity
} // namespace OpenMS
