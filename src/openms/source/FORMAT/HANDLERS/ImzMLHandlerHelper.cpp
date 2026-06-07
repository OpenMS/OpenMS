// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Aditya Sarna $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/ImzMLHandlerHelper.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>

#include <cstdint>
#include <cstring>
#include <limits>

#ifndef OPENMS_IS_BIG_ENDIAN
#if defined(OPENMS_BIG_ENDIAN)
#define OPENMS_IS_BIG_ENDIAN 1
#else
#define OPENMS_IS_BIG_ENDIAN 0
#endif
#endif

namespace OpenMS
{

namespace
{
  void throwReadError_(const String& ibd_path, const String& detail)
  {
    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, ibd_path, detail);
  }

  /// Guard against malformed imzML metadata requesting huge vector allocations.
  static constexpr uint64_t MAX_IBD_ARRAY_ELEMENTS = 100'000'000ULL;

  void validateCount_(const uint64_t count, const String& ibd_path, const char* context)
  {
    if (count > MAX_IBD_ARRAY_ELEMENTS)
    {
      throwReadError_(ibd_path,
                      String(context) + ": element count " + count + " exceeds limit of "
                      + MAX_IBD_ARRAY_ELEMENTS);
    }
    if (count > static_cast<uint64_t>(std::numeric_limits<Size>::max()))
    {
      throwReadError_(ibd_path, String(context) + ": element count exceeds platform limit");
    }
  }

  void seekIbd_(FILE* ibd, const uint64_t offset, const String& ibd_path, const char* context)
  {
    if (offset > static_cast<uint64_t>(std::numeric_limits<off_t>::max()))
    {
      throwReadError_(ibd_path,
                      String(context) + ": byte offset " + offset + " exceeds platform seek limit");
    }
    if (fseeko(ibd, static_cast<off_t>(offset), SEEK_SET) != 0)
    {
      throwReadError_(ibd_path, String(context) + ": failed to seek in .ibd");
    }
  }

  inline uint32_t swapU32_(const uint32_t v)
  {
    return ((v & 0x000000ffU) << 24) | ((v & 0x0000ff00U) << 8) | ((v & 0x00ff0000U) >> 8)
           | ((v & 0xff000000U) >> 24);
  }

  inline uint64_t swapU64_(const uint64_t v)
  {
    return ((v >> 56) & 0x00000000000000FFULL) | ((v >> 40) & 0x000000000000FF00ULL)
           | ((v >> 24) & 0x0000000000FF0000ULL) | ((v >> 8) & 0x00000000FF000000ULL)
           | ((v << 8) & 0x000000FF00000000ULL) | ((v << 24) & 0x0000FF0000000000ULL)
           | ((v << 40) & 0x00FF000000000000ULL) | ((v << 56) & 0xFF00000000000000ULL);
  }

  /// imzML .ibd arrays are little-endian (imzML 1.1.0); byte-swap on big-endian hosts.
  template<typename T>
  void decodeLittleEndian_(T* data, const Size count)
  {
#if OPENMS_IS_BIG_ENDIAN
    if constexpr (sizeof(T) == 4)
    {
      auto* raw = reinterpret_cast<uint32_t*>(data);
      for (Size i = 0; i < count; ++i)
      {
        raw[i] = swapU32_(raw[i]);
      }
    }
    else if constexpr (sizeof(T) == 8)
    {
      auto* raw = reinterpret_cast<uint64_t*>(data);
      for (Size i = 0; i < count; ++i)
      {
        raw[i] = swapU64_(raw[i]);
      }
    }
#endif
  }
} // namespace

void ImzMLBinaryIO::readMzArray(FILE* ibd,
                                const uint64_t offset,
                                const uint64_t count,
                                const ImzMLSpectrumIndex::DataType dt,
                                std::vector<double>& out,
                                const String& ibd_path)
{
  out.clear();
  if (!ibd || count == 0)
  {
    return;
  }

  validateCount_(count, ibd_path, "m/z array");

  seekIbd_(ibd, offset, ibd_path, "m/z array");

  out.resize(count);

  switch (dt)
  {
    case ImzMLSpectrumIndex::DataType::FLOAT32:
    {
      std::vector<float> tmp(count);
      if (fread(tmp.data(), 4, count, ibd) != count)
      {
        throwReadError_(ibd_path, "failed to read float32 m/z array from .ibd");
      }
      decodeLittleEndian_(tmp.data(), count);
      for (Size i = 0; i < count; ++i)
      {
        out[i] = tmp[i];
      }
      break;
    }
    case ImzMLSpectrumIndex::DataType::FLOAT64:
      if (fread(out.data(), 8, count, ibd) != count)
      {
        throwReadError_(ibd_path, "failed to read float64 m/z array from .ibd");
      }
      decodeLittleEndian_(out.data(), count);
      break;
    case ImzMLSpectrumIndex::DataType::INT32:
    {
      std::vector<int32_t> tmp(count);
      if (fread(tmp.data(), 4, count, ibd) != count)
      {
        throwReadError_(ibd_path, "failed to read int32 m/z array from .ibd");
      }
      decodeLittleEndian_(tmp.data(), count);
      for (Size i = 0; i < count; ++i)
      {
        out[i] = tmp[i];
      }
      break;
    }
    case ImzMLSpectrumIndex::DataType::INT64:
    {
      std::vector<int64_t> tmp(count);
      if (fread(tmp.data(), 8, count, ibd) != count)
      {
        throwReadError_(ibd_path, "failed to read int64 m/z array from .ibd");
      }
      decodeLittleEndian_(tmp.data(), count);
      for (Size i = 0; i < count; ++i)
      {
        out[i] = static_cast<double>(tmp[i]);
      }
      break;
    }
    default:
      throwReadError_(ibd_path, "unsupported m/z array data type in .ibd");
  }
}

void ImzMLBinaryIO::readIntArray(FILE* ibd,
                                 const uint64_t offset,
                                 const uint64_t count,
                                 const ImzMLSpectrumIndex::DataType dt,
                                 std::vector<float>& out,
                                 const String& ibd_path)
{
  out.clear();
  if (!ibd || count == 0)
  {
    return;
  }

  validateCount_(count, ibd_path, "intensity array");

  seekIbd_(ibd, offset, ibd_path, "intensity array");

  out.resize(count);

  switch (dt)
  {
    case ImzMLSpectrumIndex::DataType::FLOAT32:
      if (fread(out.data(), 4, count, ibd) != count)
      {
        throwReadError_(ibd_path, "failed to read float32 intensity array from .ibd");
      }
      decodeLittleEndian_(out.data(), count);
      break;
    case ImzMLSpectrumIndex::DataType::FLOAT64:
    {
      std::vector<double> tmp(count);
      if (fread(tmp.data(), 8, count, ibd) != count)
      {
        throwReadError_(ibd_path, "failed to read float64 intensity array from .ibd");
      }
      decodeLittleEndian_(tmp.data(), count);
      for (Size i = 0; i < count; ++i)
      {
        out[i] = static_cast<float>(tmp[i]);
      }
      break;
    }
    case ImzMLSpectrumIndex::DataType::INT32:
    {
      std::vector<int32_t> tmp(count);
      if (fread(tmp.data(), 4, count, ibd) != count)
      {
        throwReadError_(ibd_path, "failed to read int32 intensity array from .ibd");
      }
      decodeLittleEndian_(tmp.data(), count);
      for (Size i = 0; i < count; ++i)
      {
        out[i] = static_cast<float>(tmp[i]);
      }
      break;
    }
    case ImzMLSpectrumIndex::DataType::INT64:
    {
      std::vector<int64_t> tmp(count);
      if (fread(tmp.data(), 8, count, ibd) != count)
      {
        throwReadError_(ibd_path, "failed to read int64 intensity array from .ibd");
      }
      decodeLittleEndian_(tmp.data(), count);
      for (Size i = 0; i < count; ++i)
      {
        out[i] = static_cast<float>(tmp[i]);
      }
      break;
    }
    default:
      throwReadError_(ibd_path, "unsupported intensity array data type in .ibd");
  }
}

void ImzMLBinaryIO::writeFloat32Array(FILE* ibd,
                                      const float* data,
                                      const uint64_t count,
                                      const String& ibd_path)
{
  if (!ibd)
  {
    throwReadError_(ibd_path, "invalid .ibd handle (null FILE*) in float32 writer");
  }
  if (count == 0)
  {
    return;
  }
  validateCount_(count, ibd_path, "float32 array write");
#if OPENMS_IS_BIG_ENDIAN
  std::vector<float> le_data(static_cast<size_t>(count));
  std::memcpy(le_data.data(), data, static_cast<size_t>(count) * sizeof(float));
  decodeLittleEndian_(le_data.data(), static_cast<Size>(count));
  if (fwrite(le_data.data(), 4, static_cast<size_t>(count), ibd) != static_cast<size_t>(count))
  {
    throwReadError_(ibd_path, "failed to write float32 array to .ibd");
  }
#else
  if (fwrite(data, 4, static_cast<size_t>(count), ibd) != static_cast<size_t>(count))
  {
    throwReadError_(ibd_path, "failed to write float32 array to .ibd");
  }
#endif
}

void ImzMLBinaryIO::writeMzAsFloat32(FILE* ibd,
                                     const std::vector<double>& mz,
                                     const String& ibd_path)
{
  if (mz.empty())
  {
    return;
  }
  std::vector<float> tmp(mz.size());
  for (Size i = 0; i < mz.size(); ++i)
  {
    tmp[i] = static_cast<float>(mz[i]);
  }
  writeFloat32Array(ibd, tmp.data(), tmp.size(), ibd_path);
}

void ImzMLBinaryIO::writeFloat64Array(FILE* ibd,
                                      const double* data,
                                      const uint64_t count,
                                      const String& ibd_path)
{
  if (!ibd)
  {
    throwReadError_(ibd_path, "invalid .ibd handle (null FILE*) in float64 writer");
  }
  if (count == 0)
  {
    return;
  }
  validateCount_(count, ibd_path, "float64 array write");
#if OPENMS_IS_BIG_ENDIAN
  std::vector<double> le_data(static_cast<size_t>(count));
  std::memcpy(le_data.data(), data, static_cast<size_t>(count) * sizeof(double));
  decodeLittleEndian_(le_data.data(), static_cast<Size>(count));
  if (fwrite(le_data.data(), 8, static_cast<size_t>(count), ibd) != static_cast<size_t>(count))
  {
    throwReadError_(ibd_path, "failed to write float64 array to .ibd");
  }
#else
  if (fwrite(data, 8, static_cast<size_t>(count), ibd) != static_cast<size_t>(count))
  {
    throwReadError_(ibd_path, "failed to write float64 array to .ibd");
  }
#endif
}

void ImzMLBinaryIO::writeMzAsFloat64(FILE* ibd,
                                     const std::vector<double>& mz,
                                     const String& ibd_path)
{
  if (mz.empty())
  {
    return;
  }
  writeFloat64Array(ibd, mz.data(), mz.size(), ibd_path);
}

} // namespace OpenMS
