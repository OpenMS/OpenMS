// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "Sha1Helpers.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace OpenMS::Internal
{
  namespace // SHA-1 implementation per RFC 3174
  {
    inline uint32_t sha1_left_rotate(uint32_t x, uint32_t n)
    {
      return (x << n) | (x >> (32 - n));
    }

    struct SHA1
    {
      uint32_t h[5] = {0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476, 0xC3D2E1F0};
      uint64_t total_bytes = 0;
      uint8_t block[64] = {};
      size_t block_len = 0;

      void process_block()
      {
        uint32_t w[80];
        for (int i = 0; i < 16; ++i)
        {
          w[i] = (uint32_t(block[i * 4]) << 24) | (uint32_t(block[i * 4 + 1]) << 16) |
                 (uint32_t(block[i * 4 + 2]) << 8) | uint32_t(block[i * 4 + 3]);
        }
        for (int i = 16; i < 80; ++i)
        {
          w[i] = sha1_left_rotate(w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16], 1);
        }

        uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4];
        for (int i = 0; i < 80; ++i)
        {
          uint32_t f, k;
          if (i < 20)      { f = (b & c) | (~b & d);           k = 0x5A827999; }
          else if (i < 40) { f = b ^ c ^ d;                    k = 0x6ED9EBA1; }
          else if (i < 60) { f = (b & c) | (b & d) | (c & d); k = 0x8F1BBCDC; }
          else              { f = b ^ c ^ d;                    k = 0xCA62C1D6; }
          uint32_t temp = sha1_left_rotate(a, 5) + f + e + k + w[i];
          e = d; d = c; c = sha1_left_rotate(b, 30); b = a; a = temp;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e;
      }

      void update(const void* data, size_t len)
      {
        auto* p = static_cast<const uint8_t*>(data);
        total_bytes += len;
        while (len > 0)
        {
          size_t n = std::min(len, size_t(64) - block_len);
          std::memcpy(block + block_len, p, n);
          block_len += n;
          p += n;
          len -= n;
          if (block_len == 64)
          {
            process_block();
            block_len = 0;
          }
        }
      }

      void finalize(uint32_t digest[5])
      {
        uint64_t total_bits = total_bytes * 8;
        uint8_t pad = 0x80;
        update(&pad, 1);
        pad = 0;
        while (block_len != 56)
        {
          update(&pad, 1);
        }
        uint8_t len_be[8];
        for (int i = 7; i >= 0; --i)
        {
          len_be[i] = uint8_t(total_bits);
          total_bits >>= 8;
        }
        update(len_be, 8);
        std::memcpy(digest, h, sizeof(h));
      }
    };

    std::string sha1ToHexString_(const uint32_t digest[5])
    {
      std::ostringstream result;
      for (int i = 0; i < 5; ++i)
      {
        result << std::hex << std::setfill('0') << std::setw(8) << digest[i];
      }
      return result.str();
    }

    bool appendFileToSha1_(SHA1& sha, const std::string& filename)
    {
      std::ifstream file{std::filesystem::path{std::string(filename)}, std::ios::binary};
      if (!file.is_open())
      {
        return false;
      }
      char buffer[8192];
      while (file.read(buffer, sizeof(buffer)) || file.gcount() > 0)
      {
        sha.update(buffer, static_cast<std::size_t>(file.gcount()));
      }
      return true;
    }

  } // anonymous namespace

  std::string sha1FileHex(const std::string& filename)
  {
    SHA1 sha;
    if (!appendFileToSha1_(sha, filename))
    {
      return "";
    }
    uint32_t digest[5];
    sha.finalize(digest);
    return sha1ToHexString_(digest);
  }

} // namespace OpenMS::Internal
