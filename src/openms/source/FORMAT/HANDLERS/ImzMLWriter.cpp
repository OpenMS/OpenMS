// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/ImzMLWriter.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Instrument.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <random>
#include <sstream>
#include <unordered_set>
#include <vector>

namespace OpenMS
{
namespace Internal
{

namespace
{
  /// imzML reserves the first 16 bytes of the .ibd for the dataset UUID.
  constexpr uint64_t IBD_UUID_HEADER_BYTES = 16;

  struct PixelCoord
  {
    uint32_t x {1};
    uint32_t y {1};
    uint32_t z {1};
  };

  struct SpectrumWritePlan
  {
    PixelCoord pixel;
    uint64_t mz_offset {0};
    uint64_t mz_count {0};
    uint64_t mz_encoded {0};
    uint64_t int_offset {0};
    uint64_t int_count {0};
    uint64_t int_encoded {0};
    bool share_mz {false};
  };

  String inferIbdPath_(const String& imzml_path)
  {
    String lower = imzml_path;
    lower.toLower();
    if (lower.hasSuffix(".imzml"))
    {
      return imzml_path.substr(0, imzml_path.size() - 6) + ".ibd";
    }
    return imzml_path + ".ibd";
  }

  PixelCoord readPixelCoord_(const MSSpectrum& spec)
  {
    PixelCoord p;
    p.x = static_cast<uint32_t>(static_cast<Int>(spec.getMetaValue("imzml:x")));
    p.y = static_cast<uint32_t>(static_cast<Int>(spec.getMetaValue("imzml:y")));
    if (spec.metaValueExists("imzml:z"))
    {
      p.z = static_cast<uint32_t>(static_cast<Int>(spec.getMetaValue("imzml:z")));
    }
    return p;
  }

  void validatePixelMetadataForStore_(const MSExperiment& exp)
  {
    struct PixelKey
    {
      uint32_t x {0};
      uint32_t y {0};
      uint32_t z {1};

      bool operator==(const PixelKey& other) const noexcept
      {
        return x == other.x && y == other.y && z == other.z;
      }
    };

    struct PixelKeyHash
    {
      std::size_t operator()(const PixelKey& key) const noexcept
      {
        const std::size_t hx = std::hash<uint32_t>()(key.x);
        const std::size_t hy = std::hash<uint32_t>()(key.y);
        const std::size_t hz = std::hash<uint32_t>()(key.z);
        return hx ^ (hy << 1) ^ (hz << 2);
      }
    };

    std::unordered_set<PixelKey, PixelKeyHash> seen;
    seen.reserve(exp.size());

    for (Size i = 0; i < exp.size(); ++i)
    {
      const MSSpectrum& spec = exp[i];
      if (!spec.metaValueExists("imzml:x") || !spec.metaValueExists("imzml:y"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            String("spectrum ") + i + " missing imzml:x/y MetaValues required for imzML export");
      }

      const Int x_imz = spec.getMetaValue("imzml:x");
      const Int y_imz = spec.getMetaValue("imzml:y");
      if (x_imz < 1 || y_imz < 1)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "imzML pixel coordinates must be >= 1",
                                      String("(") + x_imz + "," + y_imz + ") at spectrum " + i);
      }

      Int z_imz = 1;
      if (spec.metaValueExists("imzml:z"))
      {
        z_imz = spec.getMetaValue("imzml:z");
      }
      if (z_imz < 1)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "imzML pixel z coordinate must be >= 1",
                                      String("z=") + z_imz + " at spectrum " + i);
      }

      const PixelKey key {static_cast<uint32_t>(x_imz),
                          static_cast<uint32_t>(y_imz),
                          static_cast<uint32_t>(z_imz)};
      if (!seen.emplace(key).second)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Duplicate pixel coordinate",
                                      String("(") + x_imz + "," + y_imz + "," + z_imz + ") at spectrum " + i);
      }
    }
  }

  ImzMLMeta extractMeta_(const MSExperiment& exp)
  {
    ImzMLMeta meta;
    if (exp.metaValueExists("imzml:imaging_mode"))
    {
      meta.imaging_mode = exp.getMetaValue("imzml:imaging_mode").toString();
    }
    if (exp.metaValueExists("imzml:max_count_x"))
    {
      meta.max_count_x = static_cast<uint32_t>(static_cast<UInt>(exp.getMetaValue("imzml:max_count_x")));
    }
    if (exp.metaValueExists("imzml:max_count_y"))
    {
      meta.max_count_y = static_cast<uint32_t>(static_cast<UInt>(exp.getMetaValue("imzml:max_count_y")));
    }
    if (exp.metaValueExists("imzml:max_count_z"))
    {
      meta.max_count_z = static_cast<uint32_t>(static_cast<UInt>(exp.getMetaValue("imzml:max_count_z")));
    }
    if (exp.metaValueExists("imzml:pixel_size_x"))
    {
      meta.pixel_size_x = static_cast<double>(exp.getMetaValue("imzml:pixel_size_x"));
    }
    if (exp.metaValueExists("imzml:pixel_size_y"))
    {
      meta.pixel_size_y = static_cast<double>(exp.getMetaValue("imzml:pixel_size_y"));
    }
    if (exp.metaValueExists("imzml:max_dim_x"))
    {
      meta.max_dim_x = static_cast<double>(exp.getMetaValue("imzml:max_dim_x"));
    }
    if (exp.metaValueExists("imzml:max_dim_y"))
    {
      meta.max_dim_y = static_cast<double>(exp.getMetaValue("imzml:max_dim_y"));
    }
    if (exp.metaValueExists("imzml:uuid"))
    {
      meta.uuid = exp.getMetaValue("imzml:uuid").toString();
    }
    if (exp.metaValueExists("imzml:scan_pattern"))
    {
      meta.scan_pattern = exp.getMetaValue("imzml:scan_pattern").toString();
    }
    if (exp.metaValueExists("imzml:scan_direction"))
    {
      meta.scan_direction = exp.getMetaValue("imzml:scan_direction").toString();
    }
    if (exp.metaValueExists("imzml:line_scan_direction"))
    {
      meta.line_scan_direction = exp.getMetaValue("imzml:line_scan_direction").toString();
    }
    if (exp.metaValueExists("imzml:polarity"))
    {
      meta.polarity = exp.getMetaValue("imzml:polarity").toString();
    }
    return meta;
  }

  String instrumentModelForExport_(const MSExperiment& exp)
  {
    const Instrument& inst = exp.getInstrument();
    if (!inst.getModel().empty())
    {
      return inst.getModel();
    }
    if (!inst.getName().empty())
    {
      return inst.getName();
    }
    return "OpenMS export";
  }

  void applyStoreOptions_(MSExperiment& exp, const PeakFileOptions& options)
  {
    if (options.hasMSLevels())
    {
      exp.getSpectra().erase(
        std::remove_if(exp.getSpectra().begin(), exp.getSpectra().end(),
                       [&options](const MSSpectrum& s) { return !options.containsMSLevel(s.getMSLevel()); }),
        exp.getSpectra().end());
    }

    if (options.hasRTRange())
    {
      const auto& rt_range = options.getRTRange();
      exp.getSpectra().erase(
        std::remove_if(exp.getSpectra().begin(), exp.getSpectra().end(),
                       [&rt_range](const MSSpectrum& s) { return !rt_range.encloses(DPosition<1>(s.getRT())); }),
        exp.getSpectra().end());
    }

    if (options.hasPrecursorMZRange())
    {
      const auto& prec_range = options.getPrecursorMZRange();
      exp.getSpectra().erase(
        std::remove_if(exp.getSpectra().begin(), exp.getSpectra().end(),
                       [&prec_range](const MSSpectrum& s)
                       {
                         if (s.getMSLevel() <= 1) { return false; }
                         if (s.getPrecursors().empty()) { return false; }
                         return !prec_range.encloses(DPosition<1>(s.getPrecursors()[0].getMZ()));
                       }),
        exp.getSpectra().end());
    }

    if (options.hasMZRange() || options.hasIntensityRange())
    {
      const bool filter_mz = options.hasMZRange();
      const bool filter_int = options.hasIntensityRange();
      const auto& mz_range = options.getMZRange();
      const auto& int_range = options.getIntensityRange();

      std::vector<Size> kept;
      for (auto& spectrum : exp.getSpectra())
      {
        if (options.getSortSpectraByMZ() && !spectrum.empty() && !spectrum.isSorted())
        {
          spectrum.sortByPosition();
        }

        kept.clear();
        kept.reserve(spectrum.size());
        for (Size i = 0; i < spectrum.size(); ++i)
        {
          const Peak1D& p = spectrum[i];
          if (filter_mz && !mz_range.encloses(DPosition<1>(p.getMZ()))) { continue; }
          if (filter_int && !int_range.encloses(DPosition<1>(p.getIntensity()))) { continue; }
          kept.push_back(i);
        }
        if (kept.size() != spectrum.size())
        {
          spectrum.select(kept);
        }
      }
    }
    else if (options.getSortSpectraByMZ())
    {
      for (auto& spectrum : exp.getSpectra())
      {
        if (!spectrum.empty() && !spectrum.isSorted())
        {
          spectrum.sortByPosition();
        }
      }
    }

    if (options.getMetadataOnly())
    {
      for (auto& spectrum : exp.getSpectra())
      {
        spectrum.clear(false);
      }
    }

    exp.updateRanges();
  }

  uint64_t elementByteSize_(const bool use_32_bit)
  {
    return use_32_bit ? 4ULL : 8ULL;
  }

  void writeMzArray_(FILE* ibd,
                     const std::vector<double>& mz,
                     const bool mz_32_bit,
                     const String& ibd_path)
  {
    if (mz_32_bit)
    {
      ImzMLBinaryIO::writeMzAsFloat32(ibd, mz, ibd_path);
    }
    else
    {
      ImzMLBinaryIO::writeMzAsFloat64(ibd, mz, ibd_path);
    }
  }

  void writeIntArray_(FILE* ibd,
                      const std::vector<float>& intensities,
                      const bool int_32_bit,
                      const String& ibd_path)
  {
    if (int_32_bit)
    {
      ImzMLBinaryIO::writeFloat32Array(ibd, intensities.data(), intensities.size(), ibd_path);
    }
    else
    {
      std::vector<double> tmp(intensities.size());
      for (Size i = 0; i < intensities.size(); ++i)
      {
        tmp[i] = static_cast<double>(intensities[i]);
      }
      ImzMLBinaryIO::writeFloat64Array(ibd, tmp.data(), tmp.size(), ibd_path);
    }
  }

  uint32_t sha1LeftRotate_(uint32_t value, uint32_t bits)
  {
    return (value << bits) | (value >> (32 - bits));
  }

  bool uuidStringToBytes_(const String& uuid, unsigned char out[16])
  {
    String hex = uuid;
    hex.substitute("-", "");
    if (hex.size() != 32)
    {
      return false;
    }
    for (int i = 0; i < 16; ++i)
    {
      const String byte_str = hex.substr(static_cast<Size>(i * 2), 2);
      char* end = nullptr;
      const unsigned long value = std::strtoul(byte_str.c_str(), &end, 16);
      if (end == byte_str.c_str() || value > 255)
      {
        return false;
      }
      out[i] = static_cast<unsigned char>(value);
    }
    return true;
  }

  String uuidBytesToString_(const unsigned char bytes[16])
  {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (int i = 0; i < 16; ++i)
    {
      if (i == 4 || i == 6 || i == 8 || i == 10)
      {
        oss << '-';
      }
      oss << std::setw(2) << static_cast<unsigned>(bytes[i]);
    }
    return oss.str();
  }

  void ensureUuidBytes_(ImzMLMeta& meta, unsigned char out[16])
  {
    if (!meta.uuid.empty() && uuidStringToBytes_(meta.uuid, out))
    {
      return;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);
    for (int i = 0; i < 16; ++i)
    {
      out[i] = static_cast<unsigned char>(dist(gen));
    }
    out[6] = static_cast<unsigned char>((out[6] & 0x0f) | 0x40);
    out[8] = static_cast<unsigned char>((out[8] & 0x3f) | 0x80);
    meta.uuid = uuidBytesToString_(out);
  }

  void writeIbdUuidHeader_(FILE* ibd, const unsigned char uuid[16], const String& ibd_path)
  {
    if (fwrite(uuid, 1, IBD_UUID_HEADER_BYTES, ibd) != IBD_UUID_HEADER_BYTES)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path,
                                  "failed to write imzML .ibd UUID header");
    }
  }

  void updateGridFromPixels_(ImzMLMeta& meta, const std::vector<SpectrumWritePlan>& plans)
  {
    for (const auto& plan : plans)
    {
      if (plan.pixel.x > meta.max_count_x)
      {
        meta.max_count_x = plan.pixel.x;
      }
      if (plan.pixel.y > meta.max_count_y)
      {
        meta.max_count_y = plan.pixel.y;
      }
      if (plan.pixel.z > meta.max_count_z)
      {
        meta.max_count_z = plan.pixel.z;
      }
    }
    if (meta.max_count_z == 0)
    {
      meta.max_count_z = 1;
    }
    if (meta.pixel_size_x > 0 && meta.max_count_x > 0)
    {
      meta.max_dim_x = meta.pixel_size_x * meta.max_count_x / 1000.0;
    }
    if (meta.pixel_size_y > 0 && meta.max_count_y > 0)
    {
      meta.max_dim_y = meta.pixel_size_y * meta.max_count_y / 1000.0;
    }
  }

  bool spectraShareMz_(const MSExperiment& exp, const double tolerance = 1e-5)
  {
    if (exp.empty())
    {
      return false;
    }
    Size ref = exp.size();
    for (Size i = 0; i < exp.size(); ++i)
    {
      if (!exp[i].empty())
      {
        ref = i;
        break;
      }
    }
    if (ref >= exp.size() || exp[ref].empty())
    {
      return false;
    }
    const MSSpectrum& reference = exp[ref];
    for (Size i = 0; i < exp.size(); ++i)
    {
      if (i == ref)
      {
        continue;
      }
      const MSSpectrum& s = exp[i];
      if (s.empty())
      {
        continue;
      }
      if (s.size() != reference.size())
      {
        return false;
      }
      for (Size p = 0; p < s.size(); ++p)
      {
        if (std::fabs(s[p].getMZ() - reference[p].getMZ()) > tolerance)
        {
          return false;
        }
      }
    }
    return true;
  }

  bool isContinuousMode_(const MSExperiment& exp, const ImzMLMeta& meta)
  {
    if (meta.imaging_mode == "continuous")
    {
      if (!spectraShareMz_(exp))
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                           "continuous imzML requires all non-empty spectra to share the same m/z axis and array length");
      }
      return true;
    }
    if (meta.imaging_mode == "processed")
    {
      return false;
    }
    return spectraShareMz_(exp);
  }

  std::vector<float> intensitiesAsFloat_(const MSSpectrum& spec)
  {
    std::vector<float> out(spec.size());
    for (Size i = 0; i < spec.size(); ++i)
    {
      out[i] = spec[i].getIntensity();
    }
    return out;
  }

  std::vector<double> mzAsDouble_(const MSSpectrum& spec)
  {
    std::vector<double> out(spec.size());
    for (Size i = 0; i < spec.size(); ++i)
    {
      out[i] = spec[i].getMZ();
    }
    return out;
  }

  struct UniqueFile_
  {
    FILE* fp {nullptr};

    explicit UniqueFile_(FILE* file) : fp(file) {}

    ~UniqueFile_()
    {
      if (fp != nullptr)
      {
        fclose(fp);
      }
    }

    UniqueFile_(const UniqueFile_&) = delete;
    UniqueFile_& operator=(const UniqueFile_&) = delete;

    FILE* get() const { return fp; }
  };

  void md5ProcessBlock_(uint32_t& a0, uint32_t& b0, uint32_t& c0, uint32_t& d0, const uint8_t block[64])
  {
    auto leftrotate = [](uint32_t x, uint32_t c) { return (x << c) | (x >> (32 - c)); };

    uint32_t M[16];
    for (int i = 0; i < 16; ++i)
    {
      M[i] = uint32_t(block[i * 4])
           | (uint32_t(block[i * 4 + 1]) << 8)
           | (uint32_t(block[i * 4 + 2]) << 16)
           | (uint32_t(block[i * 4 + 3]) << 24);
    }
    uint32_t A = a0, B = b0, C = c0, D = d0;
    auto round = [&](uint32_t& a, uint32_t b, uint32_t c, uint32_t d, uint32_t k, uint32_t s, uint32_t t) {
      a = b + leftrotate(a + ((b & c) | (~b & d)) + k + t, s);
    };
    auto round2 = [&](uint32_t& a, uint32_t b, uint32_t c, uint32_t d, uint32_t k, uint32_t s, uint32_t t) {
      a = b + leftrotate(a + ((b & d) | (c & ~d)) + k + t, s);
    };
    auto round3 = [&](uint32_t& a, uint32_t b, uint32_t c, uint32_t d, uint32_t k, uint32_t s, uint32_t t) {
      a = b + leftrotate(a + (b ^ c ^ d) + k + t, s);
    };
    auto round4 = [&](uint32_t& a, uint32_t b, uint32_t c, uint32_t d, uint32_t k, uint32_t s, uint32_t t) {
      a = b + leftrotate(a + (c ^ (b | ~d)) + k + t, s);
    };

    round (A, B, C, D, M[0],  7, 0xd76aa478); round (D, A, B, C, M[1],  12, 0xe8c7b756);
    round (C, D, A, B, M[2],  17, 0x242070db); round (B, C, D, A, M[3],  22, 0xc1bdceee);
    round (A, B, C, D, M[4],  7, 0xf57c0faf); round (D, A, B, C, M[5],  12, 0x4787c62a);
    round (C, D, A, B, M[6],  17, 0xa8304613); round (B, C, D, A, M[7],  22, 0xfd469501);
    round (A, B, C, D, M[8],  7, 0x698098d8); round (D, A, B, C, M[9],  12, 0x8b44f7af);
    round (C, D, A, B, M[10], 17, 0xffff5bb1); round (B, C, D, A, M[11], 22, 0x895cd7be);
    round (A, B, C, D, M[12], 7, 0x6b901122); round (D, A, B, C, M[13], 12, 0xfd987193);
    round (C, D, A, B, M[14], 17, 0xa679438e); round (B, C, D, A, M[15], 22, 0x49b40821);

    round2(A, B, C, D, M[1],  5, 0xf61e2562); round2(D, A, B, C, M[6],  9, 0xc040b340);
    round2(C, D, A, B, M[11], 14, 0x265e5a51); round2(B, C, D, A, M[0],  20, 0xe9b6c7aa);
    round2(A, B, C, D, M[5],  5, 0xd62f105d); round2(D, A, B, C, M[10], 9, 0x02441453);
    round2(C, D, A, B, M[15], 14, 0xd8a1e681); round2(B, C, D, A, M[4],  20, 0xe7d3fbc8);
    round2(A, B, C, D, M[9],  5, 0x21e1cde6); round2(D, A, B, C, M[14], 9, 0xc33707d6);
    round2(C, D, A, B, M[3],  14, 0xf4d50d87); round2(B, C, D, A, M[8],  20, 0x455a14ed);
    round2(A, B, C, D, M[13], 5, 0xa9e3e905); round2(D, A, B, C, M[2],  9, 0xfcefa3f8);
    round2(C, D, A, B, M[7],  14, 0x676f02d9); round2(B, C, D, A, M[12], 20, 0x8d2a4c8a);

    round3(A, B, C, D, M[5],  4, 0xfffa3942); round3(D, A, B, C, M[8],  11, 0x8771f681);
    round3(C, D, A, B, M[11], 16, 0x6d9d6122); round3(B, C, D, A, M[14], 23, 0xfde5380c);
    round3(A, B, C, D, M[1],  4, 0xa4beea44); round3(D, A, B, C, M[4],  11, 0x4bdecfa9);
    round3(C, D, A, B, M[7],  16, 0xf6bb4b60); round3(B, C, D, A, M[10], 23, 0xbebfbc70);
    round3(A, B, C, D, M[13], 4, 0x289b7ec6); round3(D, A, B, C, M[0],  11, 0xeaa127fa);
    round3(C, D, A, B, M[3],  16, 0xd4ef3085); round3(B, C, D, A, M[6],  23, 0x04881d05);
    round3(A, B, C, D, M[9],  4, 0xd9d4d039); round3(D, A, B, C, M[12], 11, 0xe6db99e5);
    round3(C, D, A, B, M[15], 16, 0x1fa27cf8); round3(B, C, D, A, M[2],  23, 0xc4ac5665);

    round4(A, B, C, D, M[0],  6, 0xf4292244); round4(D, A, B, C, M[7],  10, 0x432aff97);
    round4(C, D, A, B, M[14], 15, 0xab9423a7); round4(B, C, D, A, M[5],  21, 0xfc93a039);
    round4(A, B, C, D, M[12], 6, 0x655b59c3); round4(D, A, B, C, M[3],  10, 0x8f0ccc92);
    round4(C, D, A, B, M[10], 15, 0xffeff47d); round4(B, C, D, A, M[1],  21, 0x85845dd1);
    round4(A, B, C, D, M[8],  6, 0x6fa87e4f); round4(D, A, B, C, M[15], 10, 0xfe2ce6e0);
    round4(C, D, A, B, M[6],  15, 0xa3014314); round4(B, C, D, A, M[13], 21, 0x4e0811a1);
    round4(A, B, C, D, M[4],  6, 0xf7537e82); round4(D, A, B, C, M[11], 10, 0xbd3af235);
    round4(C, D, A, B, M[2],  15, 0x2ad7d2bb); round4(B, C, D, A, M[9],  21, 0xeb86d391);

    a0 += A;
    b0 += B;
    c0 += C;
    d0 += D;
  }

  String md5DigestToHex_(uint32_t a0, uint32_t b0, uint32_t c0, uint32_t d0)
  {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (uint32_t v : {a0, b0, c0, d0})
    {
      for (int i = 0; i < 4; ++i)
      {
        oss << std::setw(2) << ((v >> (8 * i)) & 0xff);
      }
    }
    return oss.str();
  }

  // RFC 1321 MD5 (hex digest), streaming over file contents.
  String md5Hex_(const String& file_path)
  {
    std::ifstream in(file_path.c_str(), std::ios::binary);
    if (!in)
    {
      return "";
    }

    uint32_t a0 = 0x67452301, b0 = 0xefcdab89, c0 = 0x98badcfe, d0 = 0x10325476;
    uint8_t block[64] {};
    size_t block_len = 0;
    uint64_t total_bytes = 0;
    char read_buf[64 * 1024];

    auto consumeBytes = [&](const uint8_t* data, size_t len) {
      size_t idx = 0;
      while (idx < len)
      {
        const size_t take = std::min(len - idx, 64 - block_len);
        std::memcpy(block + block_len, data + idx, take);
        block_len += take;
        idx += take;
        if (block_len == 64)
        {
          md5ProcessBlock_(a0, b0, c0, d0, block);
          block_len = 0;
        }
      }
    };

    while (in.read(read_buf, sizeof(read_buf)) || in.gcount() > 0)
    {
      const size_t chunk = static_cast<size_t>(in.gcount());
      total_bytes += chunk;
      consumeBytes(reinterpret_cast<const uint8_t*>(read_buf), chunk);
    }

    const uint64_t orig_bits = total_bytes * 8;
    block[block_len++] = 0x80;
    if (block_len > 56)
    {
      while (block_len < 64)
      {
        block[block_len++] = 0;
      }
      md5ProcessBlock_(a0, b0, c0, d0, block);
      block_len = 0;
    }
    while (block_len < 56)
    {
      block[block_len++] = 0;
    }
    for (int i = 0; i < 8; ++i)
    {
      block[block_len++] = static_cast<uint8_t>((orig_bits >> (8 * i)) & 0xff);
    }
    md5ProcessBlock_(a0, b0, c0, d0, block);

    return md5DigestToHex_(a0, b0, c0, d0);
  }

  String sha1DigestToHex_(const uint32_t digest[5])
  {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (int i = 0; i < 5; ++i)
    {
      oss << std::setw(8) << digest[i];
    }
    return oss.str();
  }

  String sha1Hex_(const String& file_path)
  {
    std::ifstream in(file_path.c_str(), std::ios::binary);
    if (!in)
    {
      return "";
    }

    uint32_t h[5] = {0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476, 0xC3D2E1F0};
    uint64_t total_bytes = 0;
    uint8_t block[64] {};
    size_t block_len = 0;
    char read_buf[64 * 1024];

    auto processBlock = [&](const uint8_t chunk[64]) {
      uint32_t w[80];
      for (int i = 0; i < 16; ++i)
      {
        w[i] = (uint32_t(chunk[i * 4]) << 24) | (uint32_t(chunk[i * 4 + 1]) << 16)
               | (uint32_t(chunk[i * 4 + 2]) << 8) | uint32_t(chunk[i * 4 + 3]);
      }
      for (int i = 16; i < 80; ++i)
      {
        w[i] = sha1LeftRotate_(w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16], 1);
      }

      uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4];
      for (int i = 0; i < 80; ++i)
      {
        uint32_t f = 0;
        uint32_t k = 0;
        if (i < 20)      { f = (b & c) | (~b & d);           k = 0x5A827999; }
        else if (i < 40) { f = b ^ c ^ d;                    k = 0x6ED9EBA1; }
        else if (i < 60) { f = (b & c) | (b & d) | (c & d); k = 0x8F1BBCDC; }
        else              { f = b ^ c ^ d;                    k = 0xCA62C1D6; }
        const uint32_t temp = sha1LeftRotate_(a, 5) + f + e + k + w[i];
        e = d;
        d = c;
        c = sha1LeftRotate_(b, 30);
        b = a;
        a = temp;
      }
      h[0] += a;
      h[1] += b;
      h[2] += c;
      h[3] += d;
      h[4] += e;
    };

    auto consumeBytes = [&](const uint8_t* data, size_t len) {
      size_t idx = 0;
      while (idx < len)
      {
        const size_t take = std::min(len - idx, 64 - block_len);
        std::memcpy(block + block_len, data + idx, take);
        block_len += take;
        idx += take;
        if (block_len == 64)
        {
          processBlock(block);
          block_len = 0;
        }
      }
    };

    while (in.read(read_buf, sizeof(read_buf)) || in.gcount() > 0)
    {
      const size_t chunk = static_cast<size_t>(in.gcount());
      total_bytes += chunk;
      consumeBytes(reinterpret_cast<const uint8_t*>(read_buf), chunk);
    }

    const uint64_t orig_bits = total_bytes * 8;
    block[block_len++] = 0x80;
    if (block_len > 56)
    {
      while (block_len < 64)
      {
        block[block_len++] = 0;
      }
      processBlock(block);
      block_len = 0;
    }
    while (block_len < 56)
    {
      block[block_len++] = 0;
    }
    for (int i = 7; i >= 0; --i)
    {
      block[block_len++] = static_cast<uint8_t>((orig_bits >> (8 * i)) & 0xff);
    }
    processBlock(block);

    return sha1DigestToHex_(h);
  }

  void writeCvParam_(std::ostream& os,
                     const String& cv_ref,
                     const String& accession,
                     const String& name,
                     const String& value = "",
                     const String& unit_attrs = "")
  {
    os << "<cvParam cvRef=\"" << cv_ref << "\" accession=\"" << accession << "\" name=\"" << name << "\"";
    if (!value.empty())
    {
      os << " value=\"" << XMLHandler::writeXMLEscape(value) << "\"";
    }
    if (!unit_attrs.empty())
    {
      os << " " << unit_attrs;
    }
    os << "/>";
  }

  void writeImsGeometryCvParams_(std::ostream& os, const ImzMLMeta& meta)
  {
    if (meta.scan_pattern == "top down")
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000401", "top down");
      os << "\n";
    }
    else if (meta.scan_pattern == "bottom up")
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000402", "bottom up");
      os << "\n";
    }

    if (meta.scan_direction == "flyback")
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000413", "flyback");
      os << "\n";
    }
    else if (meta.scan_direction == "meander")
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000412", "meander");
      os << "\n";
    }
    else if (meta.scan_direction == "horizontal")
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000480", "horizontal");
      os << "\n";
    }
    else if (meta.scan_direction == "vertical")
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000481", "vertical");
      os << "\n";
    }

    if (meta.line_scan_direction == "left-right")
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000491", "left-right");
      os << "\n";
    }
    else if (meta.line_scan_direction == "right-left")
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000492", "right-left");
      os << "\n";
    }
  }

  void writeExternalBinaryArray_(std::ostream& os,
                                 const String& ref_group,
                                 uint64_t offset,
                                 uint64_t count,
                                 uint64_t encoded,
                                 int indent)
  {
    const String pad(indent, '\t');
    os << pad << "<binaryDataArray encodedLength=\"0\">\n";
    os << pad << "\t<referenceableParamGroupRef ref=\"" << ref_group << "\"/>\n";
    os << pad << "\t";
    writeCvParam_(os, "IMS", "IMS:1000102", "external offset", String(offset));
    os << "\n" << pad << "\t";
    writeCvParam_(os, "IMS", "IMS:1000103", "external array length", String(count));
    os << "\n" << pad << "\t";
    writeCvParam_(os, "IMS", "IMS:1000104", "external encoded length", String(encoded));
    os << "\n";
    os << pad << "\t<binary/>\n";
    os << pad << "</binaryDataArray>\n";
  }

  void writeImzMLXml_(const String& imzml_path,
                      const MSExperiment& exp,
                      const ImzMLMeta& meta,
                      const std::vector<SpectrumWritePlan>& plans,
                      const String& ibd_md5,
                      const String& ibd_sha1,
                      const String& instrument_model,
                      const bool continuous,
                      const bool mz_32_bit,
                      const bool int_32_bit)
  {
    std::ofstream os(imzml_path.c_str(), std::ios::binary);
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, imzml_path);
    }

    const String uuid = meta.uuid;
    const String mode_acc = continuous ? "IMS:1000030" : "IMS:1000031";
    const String mode_name = continuous ? "continuous" : "processed";
    const String mz_precision_acc = mz_32_bit ? "MS:1000521" : "MS:1000523";
    const String mz_precision_name = mz_32_bit ? "32-bit float" : "64-bit float";
    const String int_precision_acc = int_32_bit ? "MS:1000521" : "MS:1000523";
    const String int_precision_name = int_32_bit ? "32-bit float" : "64-bit float";

    os << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
    os << "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" version=\"1.1.0\">\n";
    os << "\t<cvList count=\"3\">\n";
    os << "\t\t<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" version=\"4.1.30\" URI=\"https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo\"/>\n";
    os << "\t\t<cv id=\"IMS\" fullName=\"Imaging MS Ontology\" version=\"1.1.0\" URI=\"https://raw.githubusercontent.com/imzML/imzML/master/imagingMS.obo\"/>\n";
    os << "\t\t<cv id=\"UO\" fullName=\"Unit Ontology\" version=\"09:04:2014\" URI=\"https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo\"/>\n";
    os << "\t</cvList>\n";
    os << "\t<fileDescription>\n";
    os << "\t\t<fileContent>\n";
    os << "\t\t\t";
    writeCvParam_(os, "IMS", "IMS:1000080", "universally unique identifier", uuid);
    os << "\n";
    if (!ibd_sha1.empty())
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000091", "ibd SHA-1", ibd_sha1);
      os << "\n";
    }
    if (!ibd_md5.empty())
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000090", "ibd MD5", ibd_md5);
      os << "\n";
    }
    os << "\t\t\t";
    writeCvParam_(os, "IMS", mode_acc, mode_name);
    os << "\n";
    if (meta.polarity == "positive")
    {
      os << "\t\t\t";
      writeCvParam_(os, "MS", "MS:1000130", "positive scan");
      os << "\n";
    }
    else if (meta.polarity == "negative")
    {
      os << "\t\t\t";
      writeCvParam_(os, "MS", "MS:1000129", "negative scan");
      os << "\n";
    }
    os << "\t\t\t";
    writeCvParam_(os, "MS", "MS:1000294", "mass spectrum");
    os << "\n";
    os << "\t\t</fileContent>\n";
    os << "\t</fileDescription>\n";
    os << "\t<referenceableParamGroupList count=\"2\">\n";
    os << "\t\t<referenceableParamGroup id=\"mzArray\">\n";
    os << "\t\t\t";
    writeCvParam_(os, "MS", "MS:1000514", "m/z array");
    os << "\n\t\t\t";
    writeCvParam_(os, "MS", mz_precision_acc, mz_precision_name);
    os << "\n\t\t\t";
    writeCvParam_(os, "IMS", "IMS:1000101", "external data", "true");
    os << "\n";
    os << "\t\t</referenceableParamGroup>\n";
    os << "\t\t<referenceableParamGroup id=\"intensityArray\">\n";
    os << "\t\t\t";
    writeCvParam_(os, "MS", "MS:1000515", "intensity array");
    os << "\n\t\t\t";
    writeCvParam_(os, "MS", int_precision_acc, int_precision_name);
    os << "\n\t\t\t";
    writeCvParam_(os, "IMS", "IMS:1000101", "external data", "true");
    os << "\n";
    os << "\t\t</referenceableParamGroup>\n";
    os << "\t</referenceableParamGroupList>\n";
    os << "\t<softwareList count=\"1\">\n";
    os << "\t\t<software id=\"sw1\" version=\"" << VersionInfo::getVersion() << "\">\n";
    os << "\t\t\t";
    writeCvParam_(os, "MS", "MS:1000799", "custom unreleased software tool", "OpenMS");
    os << "\n";
    os << "\t\t</software>\n";
    os << "\t</softwareList>\n";
    os << "\t<scanSettingsList count=\"1\">\n";
    os << "\t\t<scanSettings id=\"scanSettings1\">\n";
    if (meta.max_count_x > 0)
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000042", "max count of pixels x", String(meta.max_count_x));
      os << "\n";
    }
    if (meta.max_count_y > 0)
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000043", "max count of pixels y", String(meta.max_count_y));
      os << "\n";
    }
    if (meta.max_dim_x > 0)
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000044", "max dimension x", String(meta.max_dim_x),
                    "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"");
      os << "\n";
    }
    if (meta.max_dim_y > 0)
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000045", "max dimension y", String(meta.max_dim_y),
                    "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"");
      os << "\n";
    }
    if (meta.pixel_size_x > 0)
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000046", "pixel size x", String(meta.pixel_size_x),
                    "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"");
      os << "\n";
    }
    if (meta.pixel_size_y > 0)
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000047", "pixel size y", String(meta.pixel_size_y),
                    "unitCvRef=\"UO\" unitAccession=\"UO:0000017\" unitName=\"micrometer\"");
      os << "\n";
    }
    writeImsGeometryCvParams_(os, meta);
    os << "\t\t</scanSettings>\n";
    os << "\t</scanSettingsList>\n";
    os << "\t<instrumentConfigurationList count=\"1\">\n";
    os << "\t\t<instrumentConfiguration id=\"IC1\">\n";
    os << "\t\t\t";
    writeCvParam_(os, "MS", "MS:1000031", "instrument model", instrument_model);
    os << "\n";
    os << "\t\t</instrumentConfiguration>\n";
    os << "\t</instrumentConfigurationList>\n";
    os << "\t<dataProcessingList count=\"1\">\n";
    os << "\t\t<dataProcessing id=\"dp1\">\n";
    os << "\t\t\t<processingMethod order=\"1\" softwareRef=\"sw1\">\n";
    os << "\t\t\t\t";
    writeCvParam_(os, "MS", "MS:1000544", "Conversion to imzML");
    os << "\n";
    os << "\t\t\t</processingMethod>\n";
    os << "\t\t</dataProcessing>\n";
    os << "\t</dataProcessingList>\n";
    os << "\t<run id=\"ru_0\" defaultInstrumentConfigurationRef=\"IC1\">\n";
    os << "\t\t<spectrumList count=\"" << exp.size() << "\" defaultDataProcessingRef=\"dp1\">\n";

    for (Size i = 0; i < exp.size(); ++i)
    {
      const MSSpectrum& spec = exp[i];
      const SpectrumWritePlan& plan = plans[i];
      String native_id = spec.getNativeID();
      if (native_id.empty())
      {
        native_id = String("spectrum=") + (i + 1);
      }

      os << "\t\t\t<spectrum index=\"" << i << "\" id=\"" << XMLHandler::writeXMLEscape(native_id)
         << "\" defaultArrayLength=\"" << spec.size() << "\">\n";
      if (continuous)
      {
        os << "\t\t\t\t<referenceableParamGroupRef ref=\"mzArray\"/>\n";
      }
      if (spec.getMSLevel() != 0)
      {
        os << "\t\t\t\t";
        writeCvParam_(os, "MS", "MS:1000511", "ms level", String(spec.getMSLevel()));
        os << "\n";
      }
      os << "\t\t\t\t<scanList count=\"1\">\n";
      os << "\t\t\t\t\t";
      writeCvParam_(os, "MS", "MS:1000795", "no combination");
      os << "\n";
      os << "\t\t\t\t\t<scan instrumentConfigurationRef=\"IC1\">\n";
      os << "\t\t\t\t\t\t";
      writeCvParam_(os, "MS", "MS:1000016", "scan start time", String(spec.getRT()),
                    "unitAccession=\"UO:0000010\" unitName=\"second\" unitCvRef=\"UO\"");
      os << "\n\t\t\t\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000050", "position x", String(plan.pixel.x));
      os << "\n\t\t\t\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000051", "position y", String(plan.pixel.y));
      os << "\n\t\t\t\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000052", "position z", String(plan.pixel.z));
      os << "\n";
      os << "\t\t\t\t\t</scan>\n";
      os << "\t\t\t\t</scanList>\n";
      os << "\t\t\t\t<binaryDataArrayList count=\"2\">\n";
      writeExternalBinaryArray_(os, "mzArray", plan.mz_offset, plan.mz_count, plan.mz_encoded, 4);
      writeExternalBinaryArray_(os, "intensityArray", plan.int_offset, plan.int_count, plan.int_encoded, 4);
      os << "\t\t\t\t</binaryDataArrayList>\n";
      os << "\t\t\t</spectrum>\n";
    }

    os << "\t\t</spectrumList>\n";
    os << "\t</run>\n";
    os << "</mzML>\n";
  }

} // namespace

void ImzMLWriter::store(const String& imzml_path,
                          const MSExperiment& exp,
                          const PeakFileOptions& options,
                          ProgressLogger& logger)
{
  MSExperiment work = exp;
  applyStoreOptions_(work, options);

  if (work.empty())
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Cannot store empty MSExperiment as imzML (all spectra removed by PeakFileOptions filters?)");
  }

  validatePixelMetadataForStore_(work);

  ImzMLMeta meta = extractMeta_(work);
  const bool continuous = isContinuousMode_(work, meta);
  meta.imaging_mode = continuous ? "continuous" : "processed";

  const bool mz_32_bit = options.getMz32Bit();
  const bool int_32_bit = options.getIntensity32Bit();
  meta.mz_data_type = mz_32_bit ? "float32" : "float64";
  meta.int_data_type = int_32_bit ? "float32" : "float64";
  const uint64_t mz_elem_bytes = elementByteSize_(mz_32_bit);
  const uint64_t int_elem_bytes = elementByteSize_(int_32_bit);

  std::vector<SpectrumWritePlan> plans(work.size());
  for (Size i = 0; i < work.size(); ++i)
  {
    plans[i].pixel = readPixelCoord_(work[i]);
    plans[i].int_count = work[i].size();
    plans[i].int_encoded = plans[i].int_count * int_elem_bytes;
    plans[i].mz_count = work[i].size();
    plans[i].mz_encoded = plans[i].mz_count * mz_elem_bytes;
  }
  updateGridFromPixels_(meta, plans);

  const String ibd_path = inferIbdPath_(imzml_path);
  const String instrument_model = instrumentModelForExport_(work);
  logger.startProgress(0, work.size() + 2, "storing imzML file");

  {
    UniqueFile_ ibd(fopen(ibd_path.c_str(), "wb"));
    if (!ibd.get())
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path);
    }

    unsigned char uuid_bytes[16] {};
    ensureUuidBytes_(meta, uuid_bytes);
    writeIbdUuidHeader_(ibd.get(), uuid_bytes, ibd_path);

    if (continuous)
    {
      Size ref = 0;
      for (Size i = 0; i < work.size(); ++i)
      {
        if (!work[i].empty())
        {
          ref = i;
          break;
        }
      }
      const std::vector<double> shared_mz = mzAsDouble_(work[ref]);
      const uint64_t mz_bytes = shared_mz.size() * mz_elem_bytes;

      writeMzArray_(ibd.get(), shared_mz, mz_32_bit, ibd_path);
      uint64_t current = IBD_UUID_HEADER_BYTES + mz_bytes;
      for (Size i = 0; i < work.size(); ++i)
      {
        logger.setProgress(static_cast<Size>(i + 1));
        plans[i].share_mz = true;
        plans[i].mz_offset = IBD_UUID_HEADER_BYTES;
        plans[i].mz_count = shared_mz.size();
        plans[i].mz_encoded = mz_bytes;
        plans[i].int_offset = current;
        plans[i].int_count = work[i].size();
        plans[i].int_encoded = plans[i].int_count * int_elem_bytes;

        const std::vector<float> intensities = intensitiesAsFloat_(work[i]);
        writeIntArray_(ibd.get(), intensities, int_32_bit, ibd_path);
        current += plans[i].int_encoded;
      }
    }
    else
    {
      uint64_t offset = IBD_UUID_HEADER_BYTES;
      for (Size i = 0; i < work.size(); ++i)
      {
        logger.setProgress(i + 1);
        plans[i].mz_offset = offset;
        plans[i].mz_count = work[i].size();
        plans[i].mz_encoded = plans[i].mz_count * mz_elem_bytes;
        offset += plans[i].mz_encoded;
        plans[i].int_offset = offset;
        plans[i].int_count = work[i].size();
        plans[i].int_encoded = plans[i].int_count * int_elem_bytes;
        offset += plans[i].int_encoded;

        writeMzArray_(ibd.get(), mzAsDouble_(work[i]), mz_32_bit, ibd_path);
        const std::vector<float> intensities = intensitiesAsFloat_(work[i]);
        writeIntArray_(ibd.get(), intensities, int_32_bit, ibd_path);
      }
    }

    if (fflush(ibd.get()) != 0)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path,
                                  "failed to flush imzML .ibd file");
    }
  }

  logger.setProgress(work.size() + 1);

  const String ibd_md5 = md5Hex_(ibd_path);
  const String ibd_sha1 = sha1Hex_(ibd_path);
  writeImzMLXml_(imzml_path, work, meta, plans, ibd_md5, ibd_sha1, instrument_model, continuous,
                 mz_32_bit, int_32_bit);
  logger.endProgress();
}

} // namespace Internal
} // namespace OpenMS
