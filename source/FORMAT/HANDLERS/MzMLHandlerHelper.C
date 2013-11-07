// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzMLHandlerHelper.h>

#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
  namespace Internal
  {
    void MzMLHandlerHelper::warning(int mode, const String & msg, UInt line, UInt column)
    {
      String error_message_;
      if (mode == 0)
        error_message_ =  String("While loading '") + "': " + msg;
      else if (mode == 1)
        error_message_ =  String("While storing '") + "': " + msg;
      if (line != 0 || column != 0)
        error_message_ += String("( in line ") + line + " column " + column + ")";
      LOG_WARN << error_message_ << std::endl;
    }

  String MzMLHandlerHelper::getCompressionTerm_(const PeakFileOptions& opt, MSNumpressCoder::NumpressConfig np, bool use_numpress)
  {
    if (np.np_compression != MSNumpressCoder::NONE && opt.getCompression() )
    {
      // TODO check if zlib AND numpress are allowed at the same time by the standard ... 
      // It is technically possible but
      //
      // MUST supply a *child* term of MS:1000572 (binary data compression type) only once
      //
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Cannot have numpress and zlib compression at the same time", "numpress, zlib");
    }

    if (np.np_compression == MSNumpressCoder::NONE || ! use_numpress)
    {
      if (opt.getCompression())
      {
        return "<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\" />";
      }
      else
      {
        return "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" />";
      }
    }
    else if (np.np_compression == MSNumpressCoder::LINEAR)
    {
      return "<cvParam cvRef=\"MS\" accession=\"MS:1002312\" name=\"MS-Numpress linear prediction compression\" />";
    }
    else if (np.np_compression == MSNumpressCoder::PIC)
    {
      return "<cvParam cvRef=\"MS\" accession=\"MS:1002313\" name=\"MS-Numpress linear prediction compression\" />";
    }
    else if (np.np_compression == MSNumpressCoder::SLOF)
    {
      return "<cvParam cvRef=\"MS\" accession=\"MS:1002314\" name=\"MS-Numpress short logged float compression\" />";
    }
    else
    {
      // default
      return "<cvParam cvRef=\"MS\" accession=\"MS:1000576\" name=\"no compression\" />";
    }
  }

  void MzMLHandlerHelper::writeFooter_(std::ostream& os, const PeakFileOptions& options_, 
    std::vector< std::pair<std::string, long> > & spectra_offsets,
    std::vector< std::pair<std::string, long> > & chromatograms_offsets)
  {
    os << "\t</run>\n";
    os << "</mzML>";

    if (options_.getWriteIndex())
    {
      int indexlists = (int) !spectra_offsets.empty() + (int) !chromatograms_offsets.empty();

      long indexlistoffset = os.tellp();
      os << "\n";
      // NOTE: indexList is required, so we need to write one 
      os << "  <indexList count=\"" << indexlists << "\">\n";
      if (!spectra_offsets.empty())
      {
        os << "    <index name=\"spectrum\">\n";
        for (Size i = 0; i < spectra_offsets.size(); i++)
        {
          os << "      <offset idRef=\"" << spectra_offsets[i].first << "\">" << spectra_offsets[i].second << "</offset>\n";
        }
        os << "    </index>\n";
      }
      if (!chromatograms_offsets.empty())
      {
        os << "    <index name=\"chromatogram\">\n";
        for (Size i = 0; i < chromatograms_offsets.size(); i++)
        {
          os << "      <offset idRef=\"" << chromatograms_offsets[i].first << "\">" << chromatograms_offsets[i].second << "</offset>\n";
        }
        os << "    </index>\n";
      }
      if (indexlists == 0)
      {
        // dummy: at least one index subelement is required by the standard,
        // and at least one offset element is required so we need to handle
        // the case where no spectra/chromatograms are present.
        os << "    <index name=\"dummy\">\n";
          os << "      <offset idRef=\"dummy\">-1</offset>\n";
        os << "    </index>\n";
      }
      os << "  </indexList>\n";
      os << "  <indexListOffset>" << indexlistoffset << "</indexListOffset>\n";
      os << "<fileChecksum>";

      // TODO calculate checksum here:
      //  SHA-1 checksum from beginning of file to end of 'fileChecksum' open tag.
      String sha1_checksum = "0";
      os << sha1_checksum << "</fileChecksum>\n";

      os << "</indexedmzML>";
    }
  }

  void MzMLHandlerHelper::decodeBase64Arrays(std::vector<BinaryData>& data_)
  {
    /// Decoder/Encoder for Base64-data in MzML
    Base64 decoder_;

    //decode all base64 arrays
    for (Size i = 0; i < data_.size(); i++)
    {
      //remove whitespaces from binary data
      //this should not be necessary, but linebreaks inside the base64 data are unfortunately no exception
      data_[i].base64.removeWhitespaces();

      //decode data and check if the length of the decoded data matches the expected length
      if (data_[i].data_type == BinaryData::DT_FLOAT)
      {
        if (data_[i].np_compression != MSNumpressCoder::NONE)
        {
          // If its numpress, we don't care about 32 / 64 bit but the numpress
          // decoder expects std::vector<double> which are 64 bit.
          MSNumpressCoder::NumpressConfig config;
          config.np_compression = data_[i].np_compression;
          MSNumpressCoder().decodeNP(data_[i].base64, data_[i].floats_64,  data_[i].compression, config);
        }
        else if (data_[i].precision == BinaryData::PRE_64)
        {
          decoder_.decode(data_[i].base64, Base64::BYTEORDER_LITTLEENDIAN, data_[i].floats_64, data_[i].compression);
          if (data_[i].size != data_[i].floats_64.size())
          {
            MzMLHandlerHelper::warning(0, String("Float binary data array '") + data_[i].meta.getName() + /* "' of spectrum '" + spec_.getNativeID() + */ "' has length " + data_[i].floats_64.size() + ", but should have length " + data_[i].size + ".");
            data_[i].size = data_[i].floats_64.size();
          }
        }
        else if (data_[i].precision == BinaryData::PRE_32)
        {
          decoder_.decode(data_[i].base64, Base64::BYTEORDER_LITTLEENDIAN, data_[i].floats_32, data_[i].compression);
          if (data_[i].size != data_[i].floats_32.size())
          {
            MzMLHandlerHelper::warning(0, String("Float binary data array '") + data_[i].meta.getName() + /* "' of spectrum '" + spec_.getNativeID() + */ "' has length " + data_[i].floats_32.size() + ", but should have length " + data_[i].size + ".");
            data_[i].size = data_[i].floats_32.size();
          }
        }
      }
      else if (data_[i].data_type == BinaryData::DT_INT)
      {
        if (data_[i].precision == BinaryData::PRE_64)
        {
          decoder_.decodeIntegers(data_[i].base64, Base64::BYTEORDER_LITTLEENDIAN, data_[i].ints_64, data_[i].compression);
          if (data_[i].size != data_[i].ints_64.size())
          {
            MzMLHandlerHelper::warning(0, String("Integer binary data array '") + data_[i].meta.getName() + /* "' of spectrum '" + spec_.getNativeID() + */ "' has length " + data_[i].ints_64.size() + ", but should have length " + data_[i].size + ".");
            data_[i].size = data_[i].ints_64.size();
          }
        }
        else if (data_[i].precision == BinaryData::PRE_32)
        {
          decoder_.decodeIntegers(data_[i].base64, Base64::BYTEORDER_LITTLEENDIAN, data_[i].ints_32, data_[i].compression);
          if (data_[i].size != data_[i].ints_32.size())
          {
            MzMLHandlerHelper::warning(0, String("Integer binary data array '") + data_[i].meta.getName() + /* "' of spectrum '" + spec_.getNativeID() + */ "' has length " + data_[i].ints_32.size() + ", but should have length " + data_[i].size + ".");
            data_[i].size = data_[i].ints_32.size();
          }
        }
      }
      else if (data_[i].data_type == BinaryData::DT_STRING)
      {
        decoder_.decodeStrings(data_[i].base64, data_[i].decoded_char, data_[i].compression);
        if (data_[i].size != data_[i].decoded_char.size())
        {
          MzMLHandlerHelper::warning(0, String("String binary data array '") + data_[i].meta.getName() + /* "' of spectrum '" + spec_.getNativeID() +  */"' has length " + data_[i].decoded_char.size() + ", but should have length " + data_[i].size + ".");
          data_[i].size = data_[i].decoded_char.size();
        }
      }
    }

  }

  void MzMLHandlerHelper::computeDataProperties_(std::vector<BinaryData>& data_, bool& precision_64, SignedSize& index, String index_name)
  {
    for (Size i = 0; i < data_.size(); i++)
    {
      if (data_[i].meta.getName() == index_name)
      {
        index = i;
        precision_64 = (data_[i].precision == BinaryData::PRE_64);
      }
    }
  }

  bool MzMLHandlerHelper::handleBinaryDataArrayCVParam(std::vector<BinaryData>& data_,
    const String& accession, const String& value, const String& name)
  {
    //MS:1000518 ! binary data type
    if (accession == "MS:1000523") //64-bit float
    {
      data_.back().precision = BinaryData::PRE_64;
      data_.back().data_type = BinaryData::DT_FLOAT;
    }
    else if (accession == "MS:1000521") //32-bit float
    {
      data_.back().precision = BinaryData::PRE_32;
      data_.back().data_type = BinaryData::DT_FLOAT;
    }
    else if (accession == "MS:1000519") //32-bit integer
    {
      data_.back().precision = BinaryData::PRE_32;
      data_.back().data_type = BinaryData::DT_INT;
    }
    else if (accession == "MS:1000522") //64-bit integer
    {
      data_.back().precision = BinaryData::PRE_64;
      data_.back().data_type = BinaryData::DT_INT;
    }
    else if (accession == "MS:1001479")
    {
      data_.back().precision = BinaryData::PRE_NONE;
      data_.back().data_type = BinaryData::DT_STRING;
    }
    //MS:1000513 ! binary data array
    else if (accession == "MS:1000786") // non-standard binary data array (with name as value)
    {
      data_.back().meta.setName(value);
    }
    //MS:1000572 ! binary data compression type
    else if (accession == "MS:1000574") //zlib compression
    {
      data_.back().compression = true;
    }
    else if (accession == "MS:1002312") //numpress compression: linear (proposed CV term)
    {
      data_.back().np_compression = MSNumpressCoder::LINEAR;
    }
    else if (accession == "MS:1002313") //numpress compression: pic (proposed CV term)
    {
      data_.back().np_compression = MSNumpressCoder::PIC;
    }
    else if (accession == "MS:1002314") //numpress compression: slof (proposed CV term)
    {
      data_.back().np_compression = MSNumpressCoder::SLOF;
    }
    else if (accession == "MS:1000576") // no compression
    {
      data_.back().compression = false;
      data_.back().np_compression = MSNumpressCoder::NONE;
    }
    else if (accession == "MS:1000514" || accession == "MS:1000515" || accession == "MS:1000595")    // handle m/z, intensity, rt
    {
      data_.back().meta.setName(name);
    }
    else
    {
      // CV term not identified
      return false;
    }

    // CV term found
    return true;
  }


  }
} // namespace OpenMS
