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

namespace OpenMS::Internal
{

    void MzMLHandlerHelper::warning(int mode, const String & msg, UInt line, UInt column)
    {
      String error_message;
      if (mode == 0)
      {
        error_message =  String("While loading '") + "': " + msg;
      }
      else if (mode == 1)
      {
        error_message =  String("While storing '") + "': " + msg;
      }
      if (line != 0 || column != 0)
      {
        error_message += String("( in line ") + line + " column " + column + ")";
      }
      OPENMS_LOG_WARN << error_message << std::endl;
    }

  String MzMLHandlerHelper::getCompressionTerm_(const PeakFileOptions& opt, MSNumpressCoder::NumpressConfig np, const String& indent, bool use_numpress)
  {
    if (opt.getCompression())
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
      String sha1_checksum = "0";
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
        bindata.base64.removeWhitespaces();
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
        MzMLHandlerHelper::warning(0, String("Invalid mzML format: Numpress-compressed binary data array '") + 
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
          MSNumpressCoder::NumpressConfig config;
          config.np_compression = bindata.np_compression;
          MSNumpressCoder().decodeNP(bindata.base64, bindata.floats_64,  bindata.compression, config);

          // Next, ensure that we only look at the float array even if the
          // mzML tags say 32 bit data (I am looking at you, proteowizard)
          bindata.precision = BinaryData::PRE_64;
        }
        else if (bindata.precision == BinaryData::PRE_64)
        {
          Base64::decode(bindata.base64, Base64::BYTEORDER_LITTLEENDIAN, bindata.floats_64, bindata.compression);
          if (bindata.size != bindata.floats_64.size())
          {
            MzMLHandlerHelper::warning(0, String("Float binary data array '") + bindata.meta.getName() + 
                "' has length " + bindata.floats_64.size() + ", but should have length " + bindata.size + ".");
            bindata.size = bindata.floats_64.size();
          }
        }
        else if (bindata.precision == BinaryData::PRE_32)
        {
          Base64::decode(bindata.base64, Base64::BYTEORDER_LITTLEENDIAN, bindata.floats_32, bindata.compression);
          if (bindata.size != bindata.floats_32.size())
          {
            MzMLHandlerHelper::warning(0, String("Float binary data array '") + bindata.meta.getName() + 
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
          Base64::decodeIntegers(bindata.base64, Base64::BYTEORDER_LITTLEENDIAN, bindata.ints_64, bindata.compression);
          if (bindata.size != bindata.ints_64.size())
          {
            MzMLHandlerHelper::warning(0, String("Integer binary data array '") + bindata.meta.getName() + 
                "' has length " + bindata.ints_64.size() + ", but should have length " + bindata.size + ".");
            bindata.size = bindata.ints_64.size();
          }
        }
        else if (bindata.precision == BinaryData::PRE_32)
        {
          Base64::decodeIntegers(bindata.base64, Base64::BYTEORDER_LITTLEENDIAN, bindata.ints_32, bindata.compression);
          if (bindata.size != bindata.ints_32.size())
          {
            MzMLHandlerHelper::warning(0, String("Integer binary data array '") + bindata.meta.getName() + 
                "' has length " + bindata.ints_32.size() + ", but should have length " + bindata.size + ".");
            bindata.size = bindata.ints_32.size();
          }
        }
      }
      else if (bindata.data_type == BinaryData::DT_STRING)
      {
        Base64::decodeStrings(bindata.base64, bindata.decoded_char, bindata.compression);
        if (bindata.size != bindata.decoded_char.size())
        {
          MzMLHandlerHelper::warning(0, String("String binary data array '") + bindata.meta.getName() + 
              "' has length " + bindata.decoded_char.size() + ", but should have length " + bindata.size + ".");
          bindata.size = bindata.decoded_char.size();
        }
      }
      else 
      {
        // TODO throw error?
        MzMLHandlerHelper::warning(0, String("Invalid mzML format: Binary data array '") + bindata.meta.getName() + 
            "' has no child term of MS:1000518 (binary data type) set. Cannot automatically deduce data type.");
      }
    }

  }

  void MzMLHandlerHelper::computeDataProperties_(const std::vector<BinaryData>& data, bool& precision_64, SignedSize& index, const String& index_name)
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
                                                       const String& accession,
                                                       const String& value,
                                                       const String& name,
                                                       const String& unit_accession)
  {
    bool is_default_array = (accession == "MS:1000514" || accession == "MS:1000515" || accession == "MS:1000595");

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
    }
    else if (accession == "MS:1002312") //numpress compression: linear
    {
      data.back().np_compression = MSNumpressCoder::LINEAR;
    }
    else if (accession == "MS:1002313") //numpress compression: pic
    {
      data.back().np_compression = MSNumpressCoder::PIC;
    }
    else if (accession == "MS:1002314") //numpress compression: slof
    {
      data.back().np_compression = MSNumpressCoder::SLOF;
    }
    else if (accession == "MS:1002746") //numpress compression: linear + zlib
    {
      data.back().np_compression = MSNumpressCoder::LINEAR;
      data.back().compression = true;
    }
    else if (accession == "MS:1002747") //numpress compression: pic + zlib
    {
      data.back().np_compression = MSNumpressCoder::PIC;
      data.back().compression = true;
    }
    else if (accession == "MS:1002748") //numpress compression: slof + zlib
    {
      data.back().np_compression = MSNumpressCoder::SLOF;
      data.back().compression = true;
    }
    else if (accession == "MS:1000576") // no compression
    {
      data.back().compression = false;
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
