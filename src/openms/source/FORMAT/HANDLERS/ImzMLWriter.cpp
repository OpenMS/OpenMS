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
#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

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
    if (exp.metaValueExists("imzml:uuid"))
    {
      meta.uuid = exp.getMetaValue("imzml:uuid").toString();
    }
    return meta;
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
                      bool continuous)
  {
    std::ofstream os(imzml_path.c_str(), std::ios::binary);
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, imzml_path);
    }

    const String uuid = meta.uuid;
    const String mode_acc = continuous ? "IMS:1000030" : "IMS:1000031";
    const String mode_name = continuous ? "continuous" : "processed";

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
    if (!ibd_md5.empty())
    {
      os << "\t\t\t";
      writeCvParam_(os, "IMS", "IMS:1000090", "ibd MD5", ibd_md5);
      os << "\n";
    }
    os << "\t\t\t";
    writeCvParam_(os, "IMS", mode_acc, mode_name);
    os << "\n";
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
    writeCvParam_(os, "MS", "MS:1000521", "32-bit float");
    os << "\n\t\t\t";
    writeCvParam_(os, "IMS", "IMS:1000101", "external data", "true");
    os << "\n";
    os << "\t\t</referenceableParamGroup>\n";
    os << "\t\t<referenceableParamGroup id=\"intensityArray\">\n";
    os << "\t\t\t";
    writeCvParam_(os, "MS", "MS:1000515", "intensity array");
    os << "\n\t\t\t";
    writeCvParam_(os, "MS", "MS:1000521", "32-bit float");
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
    os << "\t\t</scanSettings>\n";
    os << "\t</scanSettingsList>\n";
    os << "\t<instrumentConfigurationList count=\"1\">\n";
    os << "\t\t<instrumentConfiguration id=\"IC1\">\n";
    os << "\t\t\t";
    writeCvParam_(os, "MS", "MS:1000031", "instrument model", "OpenMS export");
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
      os << "\t\t\t\t<scanList count=\"1\">\n";
      os << "\t\t\t\t\t";
      writeCvParam_(os, "MS", "MS:1000795", "no combination");
      os << "\n";
      os << "\t\t\t\t\t<scan instrumentConfigurationRef=\"IC1\">\n";
      os << "\t\t\t\t\t\t";
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
                          const PeakFileOptions& /*options*/,
                          ProgressLogger& logger)
{
  if (exp.empty())
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Cannot store empty MSExperiment as imzML");
  }

  validatePixelMetadataForStore_(exp);

  ImzMLMeta meta = extractMeta_(exp);
  const bool continuous = isContinuousMode_(exp, meta);
  meta.imaging_mode = continuous ? "continuous" : "processed";

  std::vector<SpectrumWritePlan> plans(exp.size());
  for (Size i = 0; i < exp.size(); ++i)
  {
    plans[i].pixel = readPixelCoord_(exp[i]);
    plans[i].int_count = exp[i].size();
    plans[i].int_encoded = plans[i].int_count * 4;
    plans[i].mz_count = exp[i].size();
    plans[i].mz_encoded = plans[i].mz_count * 4;
  }
  updateGridFromPixels_(meta, plans);

  const String ibd_path = inferIbdPath_(imzml_path);
  logger.startProgress(0, exp.size() + 2, "storing imzML file");

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
      for (Size i = 0; i < exp.size(); ++i)
      {
        if (!exp[i].empty())
        {
          ref = i;
          break;
        }
      }
      const std::vector<double> shared_mz = mzAsDouble_(exp[ref]);
      const uint64_t mz_bytes = shared_mz.size() * 4;

      ImzMLBinaryIO::writeMzAsFloat32(ibd.get(), shared_mz, ibd_path);
      uint64_t current = IBD_UUID_HEADER_BYTES + mz_bytes;
      for (Size i = 0; i < exp.size(); ++i)
      {
        logger.setProgress(static_cast<Size>(i + 1));
        plans[i].share_mz = true;
        plans[i].mz_offset = IBD_UUID_HEADER_BYTES;
        plans[i].mz_count = shared_mz.size();
        plans[i].mz_encoded = mz_bytes;
        plans[i].int_offset = current;
        plans[i].int_count = exp[i].size();
        plans[i].int_encoded = plans[i].int_count * 4;

        const std::vector<float> intensities = intensitiesAsFloat_(exp[i]);
        ImzMLBinaryIO::writeFloat32Array(ibd.get(), intensities.data(), intensities.size(), ibd_path);
        current += plans[i].int_encoded;
      }
    }
    else
    {
      uint64_t offset = IBD_UUID_HEADER_BYTES;
      for (Size i = 0; i < exp.size(); ++i)
      {
        logger.setProgress(i + 1);
        plans[i].mz_offset = offset;
        plans[i].mz_count = exp[i].size();
        plans[i].mz_encoded = plans[i].mz_count * 4;
        offset += plans[i].mz_encoded;
        plans[i].int_offset = offset;
        plans[i].int_count = exp[i].size();
        plans[i].int_encoded = plans[i].int_count * 4;
        offset += plans[i].int_encoded;

        ImzMLBinaryIO::writeMzAsFloat32(ibd.get(), mzAsDouble_(exp[i]), ibd_path);
        const std::vector<float> intensities = intensitiesAsFloat_(exp[i]);
        ImzMLBinaryIO::writeFloat32Array(ibd.get(), intensities.data(), intensities.size(), ibd_path);
      }
    }

    if (fflush(ibd.get()) != 0)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path,
                                  "failed to flush imzML .ibd file");
    }
  }

  logger.setProgress(exp.size() + 1);

  const String ibd_md5 = md5Hex_(ibd_path);
  writeImzMLXml_(imzml_path, exp, meta, plans, ibd_md5, continuous);
  logger.endProgress();
}

} // namespace Internal
} // namespace OpenMS
