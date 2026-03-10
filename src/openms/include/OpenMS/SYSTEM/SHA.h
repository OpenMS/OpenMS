// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <cstddef>
#include <cstdint>

namespace OpenMS
{

/**
  A portable SHA-1 implementation inspired by Boost's boost::uuids::detail::sha1.
  The interface mirrors the parts we use:
   - digest_type is 5 32-bit words
   - reset, process_byte(s), get_digest
*/
class OPENMS_DLLAPI SHA
{
public:
  typedef std::uint32_t(digest_type)[5];

  SHA();

  void reset();

  void process_byte(unsigned char byte);
  void process_bytes(void const* buffer, std::size_t byte_count);

  void get_digest(digest_type& digest);

private:
  static inline std::uint32_t left_rotate(std::uint32_t x, std::size_t n)
  {
    return (x << n) | (x >> (32 - n));
  }

  void process_block();
  void process_byte_impl(unsigned char byte);

private:
  std::uint32_t h_[5];

  unsigned char block_[64];

  std::size_t block_byte_index_;
  std::size_t bit_count_low;
  std::size_t bit_count_high;
};

static_assert(sizeof(unsigned char) * 8 == 8, "OpenMS::SHA requires 8-bit unsigned char");
static_assert(sizeof(std::uint32_t) * 8 == 32, "OpenMS::SHA requires 32-bit uint32_t");

} // namespace OpenMS