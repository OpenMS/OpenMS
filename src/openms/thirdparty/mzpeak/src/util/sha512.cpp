/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/util/sha512.h"

#include <bit>
#include <cstring>

namespace MzPeak::Util {

namespace {

/******************************************************************************/
// FIPS 180-4 section 4.2.3: the first 64 bits of the fractional parts of the
// cube roots of the first eighty primes.
constexpr std::uint64_t K[80] = {
    0x428a2f98d728ae22ULL, 0x7137449123ef65cdULL, 0xb5c0fbcfec4d3b2fULL,
    0xe9b5dba58189dbbcULL, 0x3956c25bf348b538ULL, 0x59f111f1b605d019ULL,
    0x923f82a4af194f9bULL, 0xab1c5ed5da6d8118ULL, 0xd807aa98a3030242ULL,
    0x12835b0145706fbeULL, 0x243185be4ee4b28cULL, 0x550c7dc3d5ffb4e2ULL,
    0x72be5d74f27b896fULL, 0x80deb1fe3b1696b1ULL, 0x9bdc06a725c71235ULL,
    0xc19bf174cf692694ULL, 0xe49b69c19ef14ad2ULL, 0xefbe4786384f25e3ULL,
    0x0fc19dc68b8cd5b5ULL, 0x240ca1cc77ac9c65ULL, 0x2de92c6f592b0275ULL,
    0x4a7484aa6ea6e483ULL, 0x5cb0a9dcbd41fbd4ULL, 0x76f988da831153b5ULL,
    0x983e5152ee66dfabULL, 0xa831c66d2db43210ULL, 0xb00327c898fb213fULL,
    0xbf597fc7beef0ee4ULL, 0xc6e00bf33da88fc2ULL, 0xd5a79147930aa725ULL,
    0x06ca6351e003826fULL, 0x142929670a0e6e70ULL, 0x27b70a8546d22ffcULL,
    0x2e1b21385c26c926ULL, 0x4d2c6dfc5ac42aedULL, 0x53380d139d95b3dfULL,
    0x650a73548baf63deULL, 0x766a0abb3c77b2a8ULL, 0x81c2c92e47edaee6ULL,
    0x92722c851482353bULL, 0xa2bfe8a14cf10364ULL, 0xa81a664bbc423001ULL,
    0xc24b8b70d0f89791ULL, 0xc76c51a30654be30ULL, 0xd192e819d6ef5218ULL,
    0xd69906245565a910ULL, 0xf40e35855771202aULL, 0x106aa07032bbd1b8ULL,
    0x19a4c116b8d2d0c8ULL, 0x1e376c085141ab53ULL, 0x2748774cdf8eeb99ULL,
    0x34b0bcb5e19b48a8ULL, 0x391c0cb3c5c95a63ULL, 0x4ed8aa4ae3418acbULL,
    0x5b9cca4f7763e373ULL, 0x682e6ff3d6b2b8a3ULL, 0x748f82ee5defb2fcULL,
    0x78a5636f43172f60ULL, 0x84c87814a1f0ab72ULL, 0x8cc702081a6439ecULL,
    0x90befffa23631e28ULL, 0xa4506cebde82bde9ULL, 0xbef9a3f7b2c67915ULL,
    0xc67178f2e372532bULL, 0xca273eceea26619cULL, 0xd186b8c721c0c207ULL,
    0xeada7dd6cde0eb1eULL, 0xf57d4f7fee6ed178ULL, 0x06f067aa72176fbaULL,
    0x0a637dc5a2c898a6ULL, 0x113f9804bef90daeULL, 0x1b710b35131c471bULL,
    0x28db77f523047d84ULL, 0x32caab7b40c72493ULL, 0x3c9ebe0a15c9bebcULL,
    0x431d67c49c100d4cULL, 0x4cc5d4becb3e42b6ULL, 0x597f299cfc657e2aULL,
    0x5fcb6fab3ad6faecULL, 0x6c44198c4a475817ULL};

/******************************************************************************/
std::uint64_t ch(std::uint64_t x, std::uint64_t y, std::uint64_t z)
{
  return (x & y) ^ (~x & z);
}

std::uint64_t maj(std::uint64_t x, std::uint64_t y, std::uint64_t z)
{
  return (x & y) ^ (x & z) ^ (y & z);
}

std::uint64_t big_sigma0(std::uint64_t x)
{
  return std::rotr(x, 28) ^ std::rotr(x, 34) ^ std::rotr(x, 39);
}

std::uint64_t big_sigma1(std::uint64_t x)
{
  return std::rotr(x, 14) ^ std::rotr(x, 18) ^ std::rotr(x, 41);
}

std::uint64_t small_sigma0(std::uint64_t x)
{
  return std::rotr(x, 1) ^ std::rotr(x, 8) ^ (x >> 7);
}

std::uint64_t small_sigma1(std::uint64_t x)
{
  return std::rotr(x, 19) ^ std::rotr(x, 61) ^ (x >> 6);
}

/// Big-endian load, written out rather than byte-swapped so the code is
/// correct on either endianness without a conditional.
std::uint64_t load_be64(const unsigned char* p)
{
  return (static_cast<std::uint64_t>(p[0]) << 56) |
         (static_cast<std::uint64_t>(p[1]) << 48) |
         (static_cast<std::uint64_t>(p[2]) << 40) |
         (static_cast<std::uint64_t>(p[3]) << 32) |
         (static_cast<std::uint64_t>(p[4]) << 24) |
         (static_cast<std::uint64_t>(p[5]) << 16) |
         (static_cast<std::uint64_t>(p[6]) << 8) | static_cast<std::uint64_t>(p[7]);
}

} // namespace

/******************************************************************************/
Sha512::Sha512()
    : state_{0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL, 0x3c6ef372fe94f82bULL,
             0xa54ff53a5f1d36f1ULL, 0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
             0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL}
    , length_(0)
    , buffer_{}
    , buffered_(0)
{
}

/******************************************************************************/
void Sha512::compress(const unsigned char* block)
{
  std::uint64_t w[80];
  for (int t = 0; t < 16; ++t)
    w[t] = load_be64(block + t * 8);
  for (int t = 16; t < 80; ++t) {
    w[t] = small_sigma1(w[t - 2]) + w[t - 7] + small_sigma0(w[t - 15]) + w[t - 16];
  }

  std::uint64_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
  std::uint64_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];

  for (int t = 0; t < 80; ++t) {
    const std::uint64_t t1 = h + big_sigma1(e) + ch(e, f, g) + K[t] + w[t];
    const std::uint64_t t2 = big_sigma0(a) + maj(a, b, c);
    h = g;
    g = f;
    f = e;
    e = d + t1;
    d = c;
    c = b;
    b = a;
    a = t1 + t2;
  }

  state_[0] += a;
  state_[1] += b;
  state_[2] += c;
  state_[3] += d;
  state_[4] += e;
  state_[5] += f;
  state_[6] += g;
  state_[7] += h;
}

/******************************************************************************/
void Sha512::update(const void* data, std::size_t size)
{
  const unsigned char* p = static_cast<const unsigned char*>(data);
  length_ += size;

  // Top up a partial block first, so the result depends only on the byte
  // sequence and not on how the caller chunked it.
  if (buffered_ > 0) {
    const std::size_t take = std::min(sizeof(buffer_) - buffered_, size);
    std::memcpy(buffer_ + buffered_, p, take);
    buffered_ += take;
    p += take;
    size -= take;
    if (buffered_ < sizeof(buffer_)) return;
    compress(buffer_);
    buffered_ = 0;
  }

  while (size >= sizeof(buffer_)) {
    compress(p);
    p += sizeof(buffer_);
    size -= sizeof(buffer_);
  }

  if (size > 0) {
    std::memcpy(buffer_, p, size);
    buffered_ = size;
  }
}

/******************************************************************************/
std::string Sha512::hex_digest()
{
  // FIPS 180-4 section 5.1.2: a 0x80 byte, then zeros until 112 mod 128, then
  // the message length in BITS as a 128-bit big-endian integer.  The high 64
  // bits are always zero here: they would need a message of 2 exabytes.
  const std::uint64_t bits = length_ * 8;

  // buffered_ is always < 128, so pad_len lands in [1, 128] and never zero.
  unsigned char pad[128] = {};
  pad[0] = 0x80;
  const std::size_t pad_len =
      (buffered_ < 112) ? (112 - buffered_) : (240 - buffered_);

  unsigned char tail[16] = {};
  for (int i = 0; i < 8; ++i) {
    tail[15 - i] = static_cast<unsigned char>((bits >> (8 * i)) & 0xff);
  }

  update(pad, pad_len);
  update(tail, sizeof(tail));

  static const char* digits = "0123456789abcdef";
  std::string out;
  out.reserve(128);
  for (const std::uint64_t word : state_) {
    for (int i = 7; i >= 0; --i) {
      const unsigned char byte =
          static_cast<unsigned char>((word >> (8 * i)) & 0xff);
      out.push_back(digits[byte >> 4]);
      out.push_back(digits[byte & 0x0f]);
    }
  }
  return out;
}

/******************************************************************************/
std::string sha512_hex(std::string_view bytes)
{
  Sha512 hash;
  hash.update(bytes.data(), bytes.size());
  return hash.hex_digest();
}

} // namespace MzPeak::Util
