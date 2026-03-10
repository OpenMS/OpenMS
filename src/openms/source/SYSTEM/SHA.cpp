// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/SHA.h>

#include <algorithm>
#include <cstring>

namespace OpenMS
{

SHA::SHA()
{
  reset();
}

void SHA::reset()
{
  h_[0] = 0x67452301u;
  h_[1] = 0xEFCDAB89u;
  h_[2] = 0x98BADCFEu;
  h_[3] = 0x10325476u;
  h_[4] = 0xC3D2E1F0u;

  block_byte_index_ = 0;
  total_bits_ = 0;
}

void SHA::process_byte(unsigned char byte)
{
  process_byte_impl(byte);
  total_bits_ += 8;
}

void SHA::process_byte_impl(unsigned char byte)
{
  block_[block_byte_index_] = byte;
  ++block_byte_index_;

  if (block_byte_index_ == 64) {
    block_byte_index_ = 0;
    process_block();
  }
}

void SHA::process_bytes(void const* buffer, std::size_t byte_count)
{
  const auto* p = static_cast<const unsigned char*>(buffer);

  // update bit count
  total_bits_ += static_cast<std::uint64_t>(byte_count) * 8;

  // process in chunks filling the internal block
  while (byte_count > 0)
  {
    std::size_t n = std::min(byte_count, std::size_t(64) - block_byte_index_);
    std::memcpy(block_ + block_byte_index_, p, n);
    block_byte_index_ += n;
    p += n;
    byte_count -= n;
    if (block_byte_index_ == 64)
    {
      block_byte_index_ = 0;
      process_block();
    }
  }
}

void SHA::process_block()
{
  std::uint32_t w[80];

  for (std::size_t i = 0; i < 16; ++i) {
    w[i]  = (static_cast<std::uint32_t>(block_[i*4 + 0]) << 24);
    w[i] |= (static_cast<std::uint32_t>(block_[i*4 + 1]) << 16);
    w[i] |= (static_cast<std::uint32_t>(block_[i*4 + 2]) << 8);
    w[i] |= (static_cast<std::uint32_t>(block_[i*4 + 3]));
  }

  for (std::size_t i = 16; i < 80; ++i) {
    w[i] = left_rotate((w[i-3] ^ w[i-8] ^ w[i-14] ^ w[i-16]), 1);
  }

  std::uint32_t a = h_[0];
  std::uint32_t b = h_[1];
  std::uint32_t c = h_[2];
  std::uint32_t d = h_[3];
  std::uint32_t e = h_[4];

  for (std::size_t i = 0; i < 80; ++i) {
    std::uint32_t f;
    std::uint32_t k;

    if (i < 20) {
      f = (b & c) | (~b & d);
      k = 0x5A827999u;
    } else if (i < 40) {
      f = b ^ c ^ d;
      k = 0x6ED9EBA1u;
    } else if (i < 60) {
      f = (b & c) | (b & d) | (c & d);
      k = 0x8F1BBCDCu;
    } else {
      f = b ^ c ^ d;
      k = 0xCA62C1D6u;
    }

    const std::uint32_t temp = left_rotate(a, 5) + f + e + k + w[i];
    e = d;
    d = c;
    c = left_rotate(b, 30);
    b = a;
    a = temp;
  }

  h_[0] += a;
  h_[1] += b;
  h_[2] += c;
  h_[3] += d;
  h_[4] += e;
}

void SHA::get_digest(digest_type& digest)
{
  // append the bit '1' to the message
  process_byte_impl(0x80);

  // append k bits '0', where k is the minimum number >= 0
  // such that the resulting message length is congruent to 56 (mod 64)
  // check if there is enough space for padding and bit_count
  if (block_byte_index_ > 56) {
    // finish this block
    while (block_byte_index_ != 0) {
      process_byte_impl(0);
    }

    // one more block
    while (block_byte_index_ < 56) {
      process_byte_impl(0);
    }
  } else {
    while (block_byte_index_ < 56) {
      process_byte_impl(0);
    }
  }

  // append length of message (before pre-processing)
  // as a 64-bit big-endian integer
  process_byte_impl(static_cast<unsigned char>((total_bits_ >> 56) & 0xFF));
  process_byte_impl(static_cast<unsigned char>((total_bits_ >> 48) & 0xFF));
  process_byte_impl(static_cast<unsigned char>((total_bits_ >> 40) & 0xFF));
  process_byte_impl(static_cast<unsigned char>((total_bits_ >> 32) & 0xFF));
  process_byte_impl(static_cast<unsigned char>((total_bits_ >> 24) & 0xFF));
  process_byte_impl(static_cast<unsigned char>((total_bits_ >> 16) & 0xFF));
  process_byte_impl(static_cast<unsigned char>((total_bits_ >> 8)  & 0xFF));
  process_byte_impl(static_cast<unsigned char>((total_bits_)       & 0xFF));

  // process the final block (if needed)
  // If block_byte_index_ != 0 here, process_block() will be invoked by process_byte_impl above.

  // get final digest (big-endian words)
  digest[0] = h_[0];
  digest[1] = h_[1];
  digest[2] = h_[2];
  digest[3] = h_[3];
  digest[4] = h_[4];
}

} // namespace OpenMS