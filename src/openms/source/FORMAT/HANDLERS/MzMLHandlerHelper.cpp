// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzMLHandlerHelper.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/ZstdCompression.h>

#include <cstring>

namespace OpenMS::Internal
{

  namespace
  {
    /// Convert a numeric array to little-endian bytes, byte-shuffle and zstd-compress them and encode the result in Base64
    template <typename T>
    void encodeZstdNumeric(std::vector<T>& in, std::string& out)
    {
      out.clear();
      if (in.empty())
      {
        return;
      }
      if constexpr (OPENMS_IS_BIG_ENDIAN)
      {
        invertEndianess<sizeof(T)>(in.data(), in.size());
      }
      std::string compressed;
      ZstdCompression::encode(in.data(), in.size() * sizeof(T), ZstdCompression::ByteTransform::BYTE_SHUFFLE, sizeof(T), compressed);
      Base64::encodeStrings({compressed}, out, false, false);
    }

    /// Decode the Base64 string of a zstd-compressed array (with the transform given in @p bindata) into a numeric array
    template <typename T>
    void decodeZstdNumeric(const MzMLHandlerHelper::BinaryData& bindata, std::vector<T>& out)
    {
      out.clear();
      std::string compressed;
      Base64::decodeSingleString(bindata.base64, compressed, false);
      std::string decoded;
      ZstdCompression::decode(compressed.data(), compressed.size(), bindata.zstd_transform, sizeof(T), decoded);
      if (decoded.size() % sizeof(T) != 0)
      {
        throw Exception::ConversionError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "zstd-compressed binary data array '" + bindata.meta.getName() + "' has a size of " +
                                         std::to_string(decoded.size()) + " bytes, which is not a multiple of the element size " + std::to_string(sizeof(T)) + ".");
      }
      out.resize(decoded.size() / sizeof(T));
      if (!out.empty())
      {
        std::memcpy(out.data(), decoded.data(), decoded.size());
      }
      if constexpr (OPENMS_IS_BIG_ENDIAN)
      {
        invertEndianess<sizeof(T)>(out.data(), out.size());
      }
    }

    /// Decode the Base64 string of a zstd-compressed array of null-terminated strings
    void decodeZstdStrings(const MzMLHandlerHelper::BinaryData& bindata, std::vector<std::string>& out)
    {
      out.clear();
      std::string compressed;
      Base64::decodeSingleString(bindata.base64, compressed, false);
      std::string decoded;
      ZstdCompression::decode(compressed.data(), compressed.size(), bindata.zstd_transform, 1, decoded);
      // split at null bytes (skipping empty strings, identical to Base64::decodeStrings)
      size_t start = 0;
      while (start < decoded.size())
      {
        size_t end = decoded.find('\0', start);
        if (end == std::string::npos)
        {
          end = decoded.size();
        }
        if (end > start)
        {
          out.emplace_back(decoded, start, end - start);
        }
        start = end + 1;
      }
    }

    /// Decode the Base64 string of a numpress + zstd compressed array
    void decodeZstdNumpress(const MzMLHandlerHelper::BinaryData& bindata, std::vector<double>& out)
    {
      std::string compressed;
      Base64::decodeSingleString(bindata.base64, compressed, false);
      std::string numpressed;
      ZstdCompression::uncompressData(compressed.data(), compressed.size(), numpressed);
      MSNumpressCoder::NumpressConfig config;
      config.np_compression = bindata.np_compression;
      MSNumpressCoder().decodeNPRaw(numpressed, out, config);
    }
  }

    void MzMLHandlerHelper::warning(int mode, const std::string & msg, UInt line, UInt column)
    {
      std::string error_message;
      if (mode == 0)
      {
        error_message =std::string("While loading '") + "': " + msg;
      }
      else if (mode == 1)
      {
        error_message =std::string("While storing '") + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message +=std::string("( in line ") + line + " column " + column + ")";
      }
      OPENMS_LOG_WARN << error_message << std::endl;
    }

  std::string MzMLHandlerHelper::getCompressionTerm_(const PeakFileOptions& opt, MSNumpressCoder::NumpressConfig np, const std::string& indent, bool use_numpress, bool zstd_byte_shuffle)
  {
    if (opt.getZstdCompression())
    {
      if (np.np_compression == MSNumpressCoder::NONE || !use_numpress)
      {
        if (zstd_byte_shuffle)
        {
          return indent + R"(<cvParam cvRef="MS" accession="MS:1003781" name="byte-shuffled zstd compression" />)";
        }
        return indent + R"(<cvParam cvRef="MS" accession="MS:1003780" name="zstd compression" />)";
      }
      else if (np.np_compression == MSNumpressCoder::LINEAR)
      {
        return indent + R"(<cvParam cvRef="MS" accession="MS:1003783" name="MS-Numpress linear prediction compression followed by zstd compression" />)";
      }
      else if (np.np_compression == MSNumpressCoder::PIC)
      {
        return indent + R"(<cvParam cvRef="MS" accession="MS:1003784" name="MS-Numpress positive integer compression followed by zstd compression" />)";
      }
      else if (np.np_compression == MSNumpressCoder::SLOF)
      {
        return indent + R"(<cvParam cvRef="MS" accession="MS:1003785" name="MS-Numpress short logged float compression followed by zstd compression" />)";
      }
    }
    else if (opt.getCompression())
    {
      if (np.np_compression == MSNumpressCoder::NONE || !use_numpress)
      {
        return indent + R"(<cvParam cvRef="MS" accession="MS:1000574" name="zlib compression" />)";
      }
      else if (np.np_compression == MSNumpressCoder::LINEAR)
      {
        return indent + R"(<cvParam cvRef="MS" accession="MS:1002746" name="MS-Numpress linear prediction compression followed by zlib compression" />)";
      }
      else if (np.np_compression == MSNumpressCoder::PIC)
      {
        return indent + R"(<cvParam cvRef="MS" accession="MS:1002747" name="MS-Numpress positive integer compression followed by zlib compression" />)";
      }
      else if (np.np_compression == MSNumpressCoder::SLOF)
      {
        return indent + R"(<cvParam cvRef="MS" accession="MS:1002748" name="MS-Numpress short logged float compression followed by zlib compression" />)";
      }
    }
    else
    {
      if (np.np_compression == MSNumpressCoder::NONE || !use_numpress)
      {
        // default
        return indent + R"(<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />)";
      }
      else if (np.np_compression == MSNumpressCoder::LINEAR)
      {
        return indent + R"(<cvParam cvRef="MS" accession="MS:1002312" name="MS-Numpress linear prediction compression" />)";
      }
      else if (np.np_compression == MSNumpressCoder::PIC)
      {
        return indent + R"(<cvParam cvRef="MS" accession="MS:1002313" name="MS-Numpress positive integer compression" />)";
      }
      else if (np.np_compression == MSNumpressCoder::SLOF)
      {
        return indent + R"(<cvParam cvRef="MS" accession="MS:1002314" name="MS-Numpress short logged float compression" />)";
      }
    }
    // default
    return indent + R"(<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />)";
  }

  void MzMLHandlerHelper::encodeNumericArray(std::vector<float>& in, const PeakFileOptions& opt, std::string& out)
  {
    if (opt.getZstdCompression())
    {
      encodeZstdNumeric(in, out);
    }
    else
    {
      Base64::encode(in, Base64::BYTEORDER_LITTLEENDIAN, out, opt.getCompression());
    }
  }

  void MzMLHandlerHelper::encodeNumericArray(std::vector<double>& in, const PeakFileOptions& opt, std::string& out)
  {
    if (opt.getZstdCompression())
    {
      encodeZstdNumeric(in, out);
    }
    else
    {
      Base64::encode(in, Base64::BYTEORDER_LITTLEENDIAN, out, opt.getCompression());
    }
  }

  void MzMLHandlerHelper::encodeNumericArray(std::vector<Int32>& in, const PeakFileOptions& opt, std::string& out)
  {
    if (opt.getZstdCompression())
    {
      encodeZstdNumeric(in, out);
    }
    else
    {
      Base64::encodeIntegers(in, Base64::BYTEORDER_LITTLEENDIAN, out, opt.getCompression());
    }
  }

  void MzMLHandlerHelper::encodeNumericArray(std::vector<Int64>& in, const PeakFileOptions& opt, std::string& out)
  {
    if (opt.getZstdCompression())
    {
      encodeZstdNumeric(in, out);
    }
    else
    {
      Base64::encodeIntegers(in, Base64::BYTEORDER_LITTLEENDIAN, out, opt.getCompression());
    }
  }

  void MzMLHandlerHelper::encodeStringArray(const std::vector<std::string>& in, const PeakFileOptions& opt, std::string& out)
  {
    if (opt.getZstdCompression())
    {
      out.clear();
      if (in.empty())
      {
        return;
      }
      std::string raw;
      for (const auto& str : in)
      {
        raw.append(str);
        raw.push_back('\0');
      }
      std::string compressed;
      ZstdCompression::compressData(raw.data(), raw.size(), compressed);
      Base64::encodeStrings({compressed}, out, false, false);
    }
    else
    {
      Base64::encodeStrings(in, out, opt.getCompression());
    }
  }

  void MzMLHandlerHelper::encodeNumpressArray(const std::vector<double>& in, const PeakFileOptions& opt, const MSNumpressCoder::NumpressConfig& config, std::string& out)
  {
    if (opt.getZstdCompression())
    {
      out.clear();
      std::string numpressed;
      MSNumpressCoder().encodeNPRaw(in, numpressed, config);
      if (numpressed.empty())
      {
        return; // numpress failed (or empty input)
      }
      std::string compressed;
      ZstdCompression::compressData(numpressed.data(), numpressed.size(), compressed);
      Base64::encodeStrings({compressed}, out, false, false);
    }
    else
    {
      MSNumpressCoder().encodeNP(in, out, opt.getCompression(), config);
    }
  }

  void MzMLHandlerHelper::encodeNumpressArray(const std::vector<float>& in, const PeakFileOptions& opt, const MSNumpressCoder::NumpressConfig& config, std::string& out)
  {
    encodeNumpressArray(std::vector<double>(in.begin(), in.end()), opt, config, out);
  }

  void MzMLHandlerHelper::writeFooter_(std::ostream& os,
                                       const PeakFileOptions& options_, 
                                       const std::vector< std::pair<std::string, Int64> > & spectra_offsets,
                                       const std::vector< std::pair<std::string, Int64> > & chromatograms_offsets)
  {
    os << "\t</run>\n";
    os << "</mzML>";

    if (options_.getWriteIndex())
    {
      int indexlists = (int) !spectra_offsets.empty() + (int) !chromatograms_offsets.empty();

      Int64 indexlistoffset = os.tellp();
      os << "\n";
      // NOTE: indexList is required, so we need to write one 
      // NOTE: the spectra and chromatogram ids are user-supplied, so better XML-escape them!
      os << "<indexList count=\"" << indexlists << "\">\n";
      if (!spectra_offsets.empty())
      {
        os << "\t<index name=\"spectrum\">\n";
        for (Size i = 0; i < spectra_offsets.size(); i++)
        {
          os << "\t\t<offset idRef=\"" << XMLHandler::writeXMLEscape(spectra_offsets[i].first) << "\">" << spectra_offsets[i].second << "</offset>\n";
        }
        os << "\t</index>\n";
      }
      if (!chromatograms_offsets.empty())
      {
        os << "\t<index name=\"chromatogram\">\n";
        for (Size i = 0; i < chromatograms_offsets.size(); i++)
        {
          os << "\t\t<offset idRef=\"" << XMLHandler::writeXMLEscape(chromatograms_offsets[i].first) << "\">" << chromatograms_offsets[i].second << "</offset>\n";
        }
        os << "\t</index>\n";
      }
      if (indexlists == 0)
      {
        // dummy: at least one index subelement is required by the standard,
        // and at least one offset element is required so we need to handle
        // the case where no spectra/chromatograms are present.
        os << "\t<index name=\"dummy\">\n";
        os << "\t\t<offset idRef=\"dummy\">-1</offset>\n";
        os << "\t</index>\n";
      }
      os << "</indexList>\n";
      os << "<indexListOffset>" << indexlistoffset << "</indexListOffset>\n";
      os << "<fileChecksum>";

      // TODO calculate checksum here:
      // SHA-1 checksum from beginning of file to end of 'fileChecksum' open tag.
      std::string sha1_checksum = "0";
      os << sha1_checksum << "</fileChecksum>\n";

      os << "</indexedmzML>";
    }
  }

  void MzMLHandlerHelper::decodeBase64Arrays(std::vector<BinaryData>& data, const bool skipXMLCheck)
  {
    // decode all base64 arrays
    for (auto& bindata : data)
    {
      // remove whitespaces from binary data
      // this should not be necessary, but line breaks inside the base64 data are unfortunately no exception
      if (!skipXMLCheck)
      {
        StringUtils::removeWhitespaces(bindata.base64);
      }

      // Catch proteowizard invalid conversion where 
      // (i) no data type is set 
      // (ii) data type is set to integer for pic compression
      //
      // Since numpress arrays are always 64 bit and decode to double arrays,
      // this should be safe. However, we cannot generally assume that DT_NONE
      // means that we are dealing with a 64 bit float type. 
      if (bindata.np_compression != MSNumpressCoder::NONE && 
          bindata.data_type == BinaryData::DT_NONE)
      {
        MzMLHandlerHelper::warning(0,std::string("Invalid mzML format: Numpress-compressed binary data array '") + 
            bindata.meta.getName() + "' has no child term of MS:1000518 (binary data type) set. Assuming 64 bit float data type.");
        bindata.data_type = BinaryData::DT_FLOAT;
        bindata.precision = BinaryData::PRE_64;
      }
      if (bindata.np_compression == MSNumpressCoder::PIC && 
          bindata.data_type == BinaryData::DT_INT)
      {
        bindata.data_type = BinaryData::DT_FLOAT;
        bindata.precision = BinaryData::PRE_64;
      }

      // decode data and check if the length of the decoded data matches the expected length
      if (bindata.data_type == BinaryData::DT_FLOAT)
      {
        if (bindata.np_compression != MSNumpressCoder::NONE)
        {
          // If its numpress, we don't distinguish 32 / 64 bit as the numpress
          // decoder always works with 64 bit (takes std::vector<double>)
          if (bindata.zstd_compression)
          {
            decodeZstdNumpress(bindata, bindata.floats_64);
          }
          else
          {
            MSNumpressCoder::NumpressConfig config;
            config.np_compression = bindata.np_compression;
            MSNumpressCoder().decodeNP(bindata.base64, bindata.floats_64,  bindata.compression, config);
          }

          // Next, ensure that we only look at the float array even if the
          // mzML tags say 32 bit data (I am looking at you, proteowizard)
          bindata.precision = BinaryData::PRE_64;
        }
        else if (bindata.precision == BinaryData::PRE_64)
        {
          if (bindata.zstd_compression)
          {
            decodeZstdNumeric(bindata, bindata.floats_64);
          }
          else
          {
            Base64::decode(bindata.base64, Base64::BYTEORDER_LITTLEENDIAN, bindata.floats_64, bindata.compression);
          }
          if (bindata.size != bindata.floats_64.size())
          {
            MzMLHandlerHelper::warning(0,std::string("Float binary data array '") + bindata.meta.getName() + 
                "' has length " + bindata.floats_64.size() + ", but should have length " + bindata.size + ".");
            bindata.size = bindata.floats_64.size();
          }
        }
        else if (bindata.precision == BinaryData::PRE_32)
        {
          if (bindata.zstd_compression)
          {
            decodeZstdNumeric(bindata, bindata.floats_32);
          }
          else
          {
            Base64::decode(bindata.base64, Base64::BYTEORDER_LITTLEENDIAN, bindata.floats_32, bindata.compression);
          }
          if (bindata.size != bindata.floats_32.size())
          {
            MzMLHandlerHelper::warning(0,std::string("Float binary data array '") + bindata.meta.getName() + 
                "' has length " + bindata.floats_32.size() + ", but should have length " + bindata.size + ".");
            bindata.size = bindata.floats_32.size();
          }
        }

        // check for unit multiplier and correct our units (e.g. seconds vs minutes)
        double unit_multiplier = bindata.unit_multiplier;
        if (unit_multiplier != 1.0 && bindata.precision == BinaryData::PRE_64)
        {
          for (auto& it : bindata.floats_64)
          {
            it = it * unit_multiplier;
          }
        }
        else if (unit_multiplier != 1.0 && bindata.precision == BinaryData::PRE_32)
        {
          for (auto& it : bindata.floats_32)
          {
            it = it * unit_multiplier;
          }
        }
      }
      else if (bindata.data_type == BinaryData::DT_INT)
      {
        if (bindata.precision == BinaryData::PRE_64)
        {
          if (bindata.zstd_compression)
          {
            decodeZstdNumeric(bindata, bindata.ints_64);
          }
          else
          {
            Base64::decodeIntegers(bindata.base64, Base64::BYTEORDER_LITTLEENDIAN, bindata.ints_64, bindata.compression);
          }
          if (bindata.size != bindata.ints_64.size())
          {
            MzMLHandlerHelper::warning(0,std::string("Integer binary data array '") + bindata.meta.getName() + 
                "' has length " + bindata.ints_64.size() + ", but should have length " + bindata.size + ".");
            bindata.size = bindata.ints_64.size();
          }
        }
        else if (bindata.precision == BinaryData::PRE_32)
        {
          if (bindata.zstd_compression)
          {
            decodeZstdNumeric(bindata, bindata.ints_32);
          }
          else
          {
            Base64::decodeIntegers(bindata.base64, Base64::BYTEORDER_LITTLEENDIAN, bindata.ints_32, bindata.compression);
          }
          if (bindata.size != bindata.ints_32.size())
          {
            MzMLHandlerHelper::warning(0,std::string("Integer binary data array '") + bindata.meta.getName() + 
                "' has length " + bindata.ints_32.size() + ", but should have length " + bindata.size + ".");
            bindata.size = bindata.ints_32.size();
          }
        }
      }
      else if (bindata.data_type == BinaryData::DT_STRING)
      {
        if (bindata.zstd_compression)
        {
          decodeZstdStrings(bindata, bindata.decoded_char);
        }
        else
        {
          Base64::decodeStrings(bindata.base64, bindata.decoded_char, bindata.compression);
        }
        if (bindata.size != bindata.decoded_char.size())
        {
          MzMLHandlerHelper::warning(0,std::string("std::string binary data array '") + bindata.meta.getName() + 
              "' has length " + bindata.decoded_char.size() + ", but should have length " + bindata.size + ".");
          bindata.size = bindata.decoded_char.size();
        }
      }
      else 
      {
        // TODO throw error?
        MzMLHandlerHelper::warning(0,std::string("Invalid mzML format: Binary data array '") + bindata.meta.getName() + 
            "' has no child term of MS:1000518 (binary data type) set. Cannot automatically deduce data type.");
      }
    }

  }

  void MzMLHandlerHelper::computeDataProperties_(const std::vector<BinaryData>& data, bool& precision_64, SignedSize& index, const std::string& index_name)
  {
    SignedSize i(0);
    for (auto const&  bindata : data)
    {
      if (bindata.meta.getName() == index_name)
      {
        index = i;
        precision_64 = (bindata.precision == BinaryData::PRE_64);
        return;
      }
      ++i;
    }
  }

  bool MzMLHandlerHelper::handleBinaryDataArrayCVParam(std::vector<BinaryData>& data,
                                                       const std::string& accession,
                                                       const std::string& value,
                                                       const std::string& name,
                                                       const std::string& unit_accession)
  {
    bool is_default_array = (accession == "MS:1000514" || accession == "MS:1000515" || accession == "MS:1000595");

    // zstd compression (optionally preceded by a byte transform or numpress compression)
    auto setZstd = [&data](ZstdCompression::ByteTransform transform, MSNumpressCoder::NumpressCompression np)
    {
      data.back().compression = false;
      data.back().zstd_compression = true;
      data.back().zstd_transform = transform;
      data.back().np_compression = np;
    };

    // store unit accession for non-default arrays
    if (!unit_accession.empty() && !is_default_array)
    {
      data.back().meta.setMetaValue("unit_accession", unit_accession);
    }

    //MS:1000518 ! binary data type
    if (accession == "MS:1000523") //64-bit float
    {
      data.back().precision = BinaryData::PRE_64;
      data.back().data_type = BinaryData::DT_FLOAT;
    }
    else if (accession == "MS:1000521") //32-bit float
    {
      data.back().precision = BinaryData::PRE_32;
      data.back().data_type = BinaryData::DT_FLOAT;
    }
    else if (accession == "MS:1000519") //32-bit integer
    {
      data.back().precision = BinaryData::PRE_32;
      data.back().data_type = BinaryData::DT_INT;
    }
    else if (accession == "MS:1000522") //64-bit integer
    {
      data.back().precision = BinaryData::PRE_64;
      data.back().data_type = BinaryData::DT_INT;
    }
    else if (accession == "MS:1001479")
    {
      data.back().precision = BinaryData::PRE_NONE;
      data.back().data_type = BinaryData::DT_STRING;
    }
    //MS:1000513 ! binary data array
    else if (accession == "MS:1000786") // non-standard binary data array (with name as value)
    {
      data.back().meta.setName(value);
    }
    //MS:1000572 ! binary data compression type
    else if (accession == "MS:1000574") //zlib compression
    {
      data.back().compression = true;
      data.back().zstd_compression = false;
    }
    else if (accession == "MS:1002312") //numpress compression: linear
    {
      data.back().np_compression = MSNumpressCoder::LINEAR;
      data.back().zstd_compression = false;
    }
    else if (accession == "MS:1002313") //numpress compression: pic
    {
      data.back().np_compression = MSNumpressCoder::PIC;
      data.back().zstd_compression = false;
    }
    else if (accession == "MS:1002314") //numpress compression: slof
    {
      data.back().np_compression = MSNumpressCoder::SLOF;
      data.back().zstd_compression = false;
    }
    else if (accession == "MS:1002746") //numpress compression: linear + zlib
    {
      data.back().np_compression = MSNumpressCoder::LINEAR;
      data.back().compression = true;
      data.back().zstd_compression = false;
    }
    else if (accession == "MS:1002747") //numpress compression: pic + zlib
    {
      data.back().np_compression = MSNumpressCoder::PIC;
      data.back().compression = true;
      data.back().zstd_compression = false;
    }
    else if (accession == "MS:1002748") //numpress compression: slof + zlib
    {
      data.back().np_compression = MSNumpressCoder::SLOF;
      data.back().compression = true;
      data.back().zstd_compression = false;
    }
    else if (accession == "MS:1003780") // zstd compression
    {
      setZstd(ZstdCompression::ByteTransform::NONE, MSNumpressCoder::NONE);
    }
    else if (accession == "MS:1003781") // byte-shuffled zstd compression
    {
      setZstd(ZstdCompression::ByteTransform::BYTE_SHUFFLE, MSNumpressCoder::NONE);
    }
    else if (accession == "MS:1003782") // dictionary-encoded zstd compression
    {
      setZstd(ZstdCompression::ByteTransform::DICTIONARY, MSNumpressCoder::NONE);
    }
    else if (accession == "MS:1003783") // numpress compression: linear + zstd
    {
      setZstd(ZstdCompression::ByteTransform::NONE, MSNumpressCoder::LINEAR);
    }
    else if (accession == "MS:1003784") // numpress compression: pic + zstd
    {
      setZstd(ZstdCompression::ByteTransform::NONE, MSNumpressCoder::PIC);
    }
    else if (accession == "MS:1003785") // numpress compression: slof + zstd
    {
      setZstd(ZstdCompression::ByteTransform::NONE, MSNumpressCoder::SLOF);
    }
    else if (accession == "MS:1000576") // no compression
    {
      data.back().compression = false;
      data.back().zstd_compression = false;
      data.back().np_compression = MSNumpressCoder::NONE;
    }
    else if (is_default_array) // handle m/z, intensity, rt
    {
      data.back().meta.setName(name);

      // time array is given in minutes instead of seconds, we need to convert
      if (accession == "MS:1000595" && unit_accession == "UO:0000031")
      {
        data.back().unit_multiplier = 60.0;
      }
    }
    else
    {
      // CV term not identified
      return false;
    }

    // CV term found
    return true;
  }


} // namespace OpenMS // namespace Internal
