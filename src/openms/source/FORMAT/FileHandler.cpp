// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/EDTAFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/MS2File.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/MSPGenericFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/MzQCFile.h>
#include <OpenMS/FORMAT/OMSSAXMLFile.h>
#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/FORMAT/ProtXMLFile.h>
#include <OpenMS/FORMAT/QcMLFile.h>
#include <OpenMS/FORMAT/SqMassFile.h>
#include <OpenMS/FORMAT/XMassFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/FORMAT/PSMArrowIO.h>
#include <OpenMS/FORMAT/FeatureMapArrowIO.h>
#include <OpenMS/FORMAT/ConsensusMapArrowIO.h>

#include <OpenMS/FORMAT/MsInspectFile.h>
#include <OpenMS/FORMAT/SpecArrayFile.h>
#include <OpenMS/FORMAT/KroenikFile.h>

#include <OpenMS/KERNEL/ChromatogramTools.h>

#include <OpenMS/FORMAT/GzipIfstream.h>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>
#include <OpenMS/FORMAT/ZipIfstream.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>

#ifdef WITH_OPENTIMS
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#endif

#ifdef WITH_THERMO_RAW
#include <OpenMS/FORMAT/ThermoRawFile.h>
#endif

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace std;

namespace OpenMS
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
  } // anonymous namespace

  std::string allowedToString_(vector<FileTypes::Type> types)
  {
    std::string aStrings;
    for (auto i : types)
    {
      if (i != FileTypes::SIZE_OF_TYPE)
      {
        aStrings +=  ", " + FileTypes::typeToName(i);
      }
    }
    return aStrings;
  }

  FileTypes::Type FileHandler::getType(const std::string& filename)
  {
    // Strip trailing directory separators (important for directory-based formats
    // like .d where shell tab-completion appends '/')
    std::string normalized = filename;
    while (StringUtils::hasSuffix(normalized, "/") || StringUtils::hasSuffix(normalized, "\\"))
    {
      normalized = StringUtils::prefix(normalized, normalized.size() - 1);
    }

    FileTypes::Type type = getTypeByFileName(normalized);

    // Directory-based formats: validate marker files before returning.
    // Must happen before getTypeByContent() which would fail on directories.
#ifdef WITH_OPENTIMS
    if (type == FileTypes::BRUKER_TDF)
    {
      if (File::exists(normalized + "/analysis.tdf") || File::exists(normalized + "/analysis.tdf_bin"))
      {
        return FileTypes::BRUKER_TDF;
      }
      // Check for .d.zip: a ZIP archive containing a Bruker .d directory
      if (File::exists(normalized) && !File::isDirectory(normalized) && StringUtils::hasSuffix(StringUtils::toLowered(normalized), ".zip"))
      {
        return FileTypes::BRUKER_TDF;
      }
      return FileTypes::UNKNOWN; // .d suffix but not a TDF directory
    }
#else
    if (type == FileTypes::BRUKER_TDF)
    {
      return FileTypes::UNKNOWN; // .d format requires WITH_OPENTIMS
    }
#endif

    // Note: .raw is always recognized as FileTypes::RAW, even without WITH_THERMO_RAW.
    // The in-process reader (loadExperiment below) is #ifdef WITH_THERMO_RAW guarded,
    // but FileConverter can still convert .raw via the external ThermoRawFileParser
    // (mono) path, which is available on all platforms. Callers that try to load a
    // .raw in-process without WITH_THERMO_RAW get a clear "type is not supported"
    // ParseError from loadExperiment's default case.

    if (type == FileTypes::UNKNOWN)
    {
      type = getTypeByContent(normalized);
    }
    return type;
  }

  FileTypes::Type FileHandler::getTypeByFileName(const std::string& filename)
  {
    std::string basename = File::basename(filename), tmp;
    // special rules for "double extensions":
    if (StringUtils::hasSuffix(basename, ".pep.xml"))
    {
      return FileTypes::PEPXML;
    }
    if (StringUtils::hasSuffix(basename, ".prot.xml"))
    {
      return FileTypes::PROTXML;
    }
    if (StringUtils::hasSuffix(basename, ".xquest.xml"))
    {
      return FileTypes::XQUESTXML;
    }
    if (StringUtils::hasSuffix(basename, ".spec.xml"))
    {
      return FileTypes::SPECXML;
    }
    try
    {
      tmp = StringUtils::suffix(basename, '.');
    }
    // no '.' => unknown type
    catch (Exception::ElementNotFound&)
    {
      // last chance, Bruker fid file
      if (basename == "fid")
      {
        return FileTypes::XMASS;
      }
      return FileTypes::UNKNOWN;
    }
    StringUtils::toUpper(tmp);
    if (tmp == "BZ2" || tmp == "GZ" || tmp == "ZIP")
    {
      // do not use getTypeByContent() here, as this is deadly for output files!
      return getTypeByFileName(StringUtils::prefix(filename, filename.size() - tmp.size() - 1)); // check name without compression suffix (e.g. bla.mzML.gz --> bla.mzML)
    }

    return FileTypes::nameToType(tmp);
  }

  bool FileHandler::hasValidExtension(const std::string& filename, const FileTypes::Type type)
  {
    FileTypes::Type ft = FileHandler::getTypeByFileName(filename);
    return (ft == type || ft == FileTypes::UNKNOWN);
  }

  std::string FileHandler::stripExtension(const std::string& filename)
  {
    if (!StringUtils::has(filename, '.'))
    {
      return filename;
    }
    // we don't just search for the last '.' and remove the suffix, because this could be wrong, e.g. bla.mzML.gz would become bla.mzML
    auto type = getTypeByFileName(filename);
    auto s_type = FileTypes::typeToName(type);
    size_t pos = StringUtils::toLowered(filename).rfind(StringUtils::toLowered(s_type)); // search backwards in entire string, because we could search for 'mzML' and have 'mzML.gz'
    if (pos == string::npos) // file type was FileTypes::UNKNOWN and we did not find '.unknown' as ending
    {
      size_t ext_pos = filename.rfind('.');
      size_t dir_sep = filename.find_last_of("/\\"); // look for '/' or '\'
      if (dir_sep != string::npos && dir_sep > ext_pos) // we found a directory separator after the last '.', e.g. '/my.dotted.dir/filename'! Ouch!
      { // do not strip anything, because there is no extension to strip
        return filename;
      }
      return StringUtils::prefix(filename, ext_pos);
    }
    return StringUtils::prefix(filename, pos - 1); // strip the '.' as well
  }

  std::string FileHandler::swapExtension(const std::string& filename, const FileTypes::Type new_type)
  {
    return stripExtension(filename) + "." + FileTypes::typeToName(new_type);
  }

  bool FileHandler::isSupported(FileTypes::Type type)
  {
    if (type == FileTypes::UNKNOWN || type == FileTypes::SIZE_OF_TYPE)
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  FileTypes::Type FileHandler::getConsistentOutputfileType(const std::string& output_filename, const std::string& requested_type)
  {
    FileTypes::Type t_file = getTypeByFileName(output_filename);
    FileTypes::Type t_req = FileTypes::nameToType(requested_type);
    // both UNKNOWN
    if (t_file == FileTypes::Type::UNKNOWN && t_req == FileTypes::Type::UNKNOWN) 
    {
      OPENMS_LOG_ERROR << "Type of '" << output_filename << "' and requested output type '" << requested_type << "' are both unknown." << std::endl;
      return FileTypes::Type::UNKNOWN;
    }
    // or inconsistent (while both are known)
    if ((t_file != t_req) && (t_file != FileTypes::Type::UNKNOWN) + (t_req != FileTypes::Type::UNKNOWN) == 2)
    {
      OPENMS_LOG_ERROR << "Type of '" << output_filename << "' and requested output type '" << requested_type << "' are inconsistent." << std::endl;
      return FileTypes::Type::UNKNOWN;
    }

    if (t_file != FileTypes::Type::UNKNOWN)
    {
      return t_file;
    }
    else
    {
      return t_req;
    }
  }

  FileTypes::Type FileHandler::getTypeByContent(const std::string& filename)
  {
    std::string first_line;
    std::string two_five;
    std::string all_simple;

    // only the first five lines will be set for compressed files
    // so far, compression is only supported for XML files
    vector<std::string> complete_file;

    // test whether the file is compressed (bzip2, gzip, or zip)
    ifstream compressed_file(filename.c_str());
    char bz[4] = {};
    compressed_file.read(bz, 4);
    char g1 = 0x1f;
    char g2 = 0;
    g2 |= 1 << 7;
    g2 |= 1 << 3;
    g2 |= 1 << 1;
    g2 |= 1 << 0;
    compressed_file.close();
    if (bz[0] == 'B' && bz[1] == 'Z') // bzip2
    {
      Bzip2Ifstream bzip2_file(filename.c_str());

      // read in 1024 bytes (keep last byte for zero to end string)
      char buffer[1024];
      size_t bytes_read = bzip2_file.read(buffer, 1024-1);
      buffer[bytes_read] = '\0';

      // get first five lines
      std::string buffer_str(buffer);
      vector<std::string> split;
      StringUtils::split(buffer_str, '\n', split);
      split.resize(5);

      first_line = split[0];
      two_five = split[1] + ' ' + split[2] + ' ' + split[3] + ' ' + split[4];
      all_simple = first_line + ' ' + two_five;
      complete_file = split;
    }
    else if (bz[0] == g1 && bz[1] == g2) // gzip
    {
      GzipIfstream gzip_file(filename.c_str());

      // read in 1024 bytes (keep last byte for zero to end string)
      char buffer[1024];
      size_t bytes_read = gzip_file.read(buffer, 1024-1);
      buffer[bytes_read] = '\0';

      // get first five lines
      std::string buffer_str(buffer);
      vector<std::string> split;
      StringUtils::split(buffer_str, '\n', split);
      split.resize(5);

      first_line = split[0];
      two_five = split[1] + ' ' + split[2] + ' ' + split[3] + ' ' + split[4];
      all_simple = first_line + ' ' + two_five;
      complete_file = split;
    }
    else if (bz[0] == 'P' && bz[1] == 'K' && bz[2] == 0x03 && bz[3] == 0x04) // ZIP local file header
    {
      ZipIfstream zip_file(filename.c_str());

      // read in 1024 bytes (keep last byte for zero to end string)
      char buffer[1024];
      size_t bytes_read = zip_file.read(buffer, 1024-1);
      buffer[bytes_read] = '\0';

      // get first five lines
      std::string buffer_str(buffer);
      vector<std::string> split;
      StringUtils::split(buffer_str, '\n', split);
      split.resize(5);

      first_line = split[0];
      two_five = split[1] + ' ' + split[2] + ' ' + split[3] + ' ' + split[4];
      all_simple = first_line + ' ' + two_five;
      complete_file = split;
    }
    else // uncompressed
    {
      //load first 5 lines
      TextFile file(filename, true, 5);
      TextFile::ConstIterator file_it = file.begin();

      // file could be empty
      if (file_it == file.end())
      {
        two_five = " ";
        all_simple = " ";
        first_line = " ";
      }
      else
      {
        // concat elements 2 to 5
        two_five = "";
        ++file_it;
        for (int i = 1; i < 5; ++i)
        {
          if (file_it != file.end())
          {
            two_five += *file_it;
            ++file_it;
          }
          else
          {
            two_five += "";
          }
          two_five += " ";
        }

        // remove trailing space
        two_five = StringUtils::chop(two_five, 1);
        StringUtils::substitute(two_five, '\t', ' ');
        all_simple = *(file.begin()) + ' ' + two_five;
        first_line = *(file.begin());
      }

      complete_file.insert(complete_file.end(), file.begin(), file.end());
    }
    //std::cerr << "\n Line1:\n" << first_line << "\nLine2-5:\n" << two_five << "\nall:\n" << all_simple << "\n\n";


    //mzXML (all lines)
    if (StringUtils::hasSubstring(all_simple, "<mzXML"))
    {
      return FileTypes::MZXML;
    }
    //mzData (all lines)
    if (StringUtils::hasSubstring(all_simple, "<mzData"))
    { 
      return FileTypes::MZDATA;
    }
    //mzML (all lines)
    if (StringUtils::hasSubstring(all_simple, "<mzML"))
    {
      return FileTypes::MZML;
    }
    //"analysisXML" aka. mzid (all lines)
    if (StringUtils::hasSubstring(all_simple, "<MzIdentML"))
    {
      return FileTypes::MZIDENTML;
    }
    //subject to change!
    if (StringUtils::hasSubstring(all_simple, "<MzQualityMLType"))
    {
      return FileTypes::QCML;
    }
    //pepXML (all lines)
    if (StringUtils::hasSubstring(all_simple, "xmlns=\"http://regis-web.systemsbiology.net/pepXML\""))
    {
      return FileTypes::PEPXML;
    }
    //protXML (all lines)
    if (StringUtils::hasSubstring(all_simple, "xmlns=\"http://regis-web.systemsbiology.net/protXML\""))
    {
      return FileTypes::PROTXML;
    }
    //feature map (all lines)
    if (StringUtils::hasSubstring(all_simple, "<featureMap"))
    {
      return FileTypes::FEATUREXML;
    }
    //idXML (all lines)
    if (StringUtils::hasSubstring(all_simple, "<IdXML"))
    {
      return FileTypes::IDXML;
    }
    //consensusXML (all lines)
    if (StringUtils::hasSubstring(all_simple, "<consensusXML"))
    {
      return FileTypes::CONSENSUSXML;
    }
    //TOPPAS (all lines)
    if (StringUtils::hasSubstring(all_simple, "<PARAMETERS") && StringUtils::hasSubstring(all_simple, "<NODE name=\"info\"") && StringUtils::hasSubstring(all_simple, "<ITEM name=\"num_vertices\""))
    {
      return FileTypes::TOPPAS;
    }
    //INI (all lines) (must be AFTER TOPPAS) - as this is less restrictive
    if (StringUtils::hasSubstring(all_simple, "<PARAMETERS"))
    {
      return FileTypes::INI;
    }
    //TrafoXML (all lines)
    if (StringUtils::hasSubstring(all_simple, "<TrafoXML"))
    {
      return FileTypes::TRANSFORMATIONXML;
    }
    //GelML (all lines)
    if (StringUtils::hasSubstring(all_simple, "<GelML"))
    {
      return FileTypes::GELML;
    }
    //traML (all lines)
    if (StringUtils::hasSubstring(all_simple, "<TraML"))
    {
      return FileTypes::TRAML;
    }
    //OMSSAXML file
    if (StringUtils::hasSubstring(all_simple, "<MSResponse"))
    {
      return FileTypes::OMSSAXML;
    }
    //MASCOTXML file
    if (StringUtils::hasSubstring(all_simple, "<mascot_search_results"))
    {
      return FileTypes::MASCOTXML;
    }
    if (StringUtils::hasPrefix(all_simple, "{"))
    {
      return FileTypes::JSON;
    }
    //FASTA file
    // .. check this fairly early on, because other file formats might be less specific
    {
      Size i = 0;
      Size bigger_than = 0;
      while (i < complete_file.size())
      {
        StringUtils::trim(complete_file[i]);
        if (StringUtils::hasPrefix(complete_file[i], ">"))
        {
          ++bigger_than;
          ++i;
        }
        else if (StringUtils::hasPrefix(complete_file[i], "#"))
        {
          ++i;
        }
        else
        {
          break;
        }
      }
      if (bigger_than > 0)
      {
        return FileTypes::FASTA;
      }
    }

    // PNG file (to be really correct, the first eight bytes of the file would
    // have to be checked; see e.g. the Wikipedia article)
    if (StringUtils::substr(first_line, 1, 3) == "PNG")
    {
      return FileTypes::PNG;
    }
    //MSP (all lines)
    for (Size i = 0; i != complete_file.size(); ++i)
    {
      if (StringUtils::hasPrefix(complete_file[i], "Name: ") && StringUtils::hasSubstring(complete_file[i], "/"))
      {
        return FileTypes::MSP;
      }
      if (StringUtils::hasPrefix(complete_file[i], "Num peaks: "))
      {
        return FileTypes::MSP;
      }
    }

    //tokenize lines 2-5
    vector<std::string> parts;
    StringUtils::split(two_five, ' ', parts);

    //DTA
    if (parts.size() == 8)
    {
      bool conversion_error = false;
      try
      {
        for (Size i = 0; i < 8; ++i)
        {
          StringUtils::toFloat(parts[i]);
        }
      }
      catch ( Exception::ConversionError& )
      {
        conversion_error = true;
      }
      if (!conversion_error)
      {
        return FileTypes::DTA;
      }
    }

    //DTA2D
    if (parts.size() == 12)
    {
      bool conversion_error = false;
      try
      {
        for (Size i = 0; i < 12; ++i)
        {
          StringUtils::toFloat(parts[i]);
        }
      }
      catch ( Exception::ConversionError& )
      {
        conversion_error = true;
      }
      if (!conversion_error)
      {
        return FileTypes::DTA2D;
      }
    }

    // MGF (Mascot Generic Format)
    if (StringUtils::hasSubstring(two_five, "BEGIN IONS"))
    {
      return FileTypes::MGF;
    }
    else
    {
      for (Size i = 0; i != complete_file.size(); ++i)
      {
        if (StringUtils::trim(complete_file[i]) == "FORMAT=Mascot generic" || StringUtils::trim(complete_file[i]) == "BEGIN IONS")
        {
          return FileTypes::MGF;
        }
      }
    }

    // MS2 file format
    if (StringUtils::hasSubstring(all_simple, "CreationDate"))
    {
      if (!all_simple.empty() && all_simple[0] == 'H')
      {
        return FileTypes::MS2;
      }
    }

    // mzTab file format
    for (Size i = 0; i != complete_file.size(); ++i) {
        if (StringUtils::hasSubstring(complete_file[i], "MTD\tmzTab-version")) {
            return FileTypes::MZTAB;
        }
    }

    // msInspect file (.tsv)
    for (Size i = 0; i != complete_file.size(); ++i)
    {
      if (StringUtils::hasSubstring(complete_file[i], "scan\ttime\tmz\taccurateMZ\tmass\tintensity\tcharge\tchargeStates\tkl\tbackground\tmedian\tpeaks\tscanFirst\tscanLast\tscanCount\ttotalIntensity\tsumSquaresDist\tdescription"))
      {
        return FileTypes::TSV;
      }
    }

    // specArray file (.pepList)
    if (StringUtils::hasSubstring(first_line, "       m/z\t     rt(min)\t       snr\t      charge\t   intensity"))
    {
      return FileTypes::PEPLIST;
    }

    // hardkloer file (.hardkloer)
    /**
      NOT IMPLEMENTED YET
      if (StringUtils::hasSubstring(first_line, "File	First Scan	Last Scan	Num of Scans	Charge	Monoisotopic Mass	Base Isotope Peak	Best Intensity	Summed Intensity	First RTime	Last RTime	Best RTime	Best Correlation	Modifications"))
    {
        return FileTypes::HARDKLOER;
    }
    **/

    // kroenik file (.kroenik)
    if (StringUtils::hasSubstring(first_line, "File\tFirst Scan\tLast Scan\tNum of Scans\tCharge\tMonoisotopic Mass\tBase Isotope Peak\tBest Intensity\tSummed Intensity\tFirst RTime\tLast RTime\tBest RTime\tBest Correlation\tModifications"))
    {
      return FileTypes::KROENIK;
    }

    // Percolator tab-delimited output (PSM level, .psms)
    if (StringUtils::hasPrefix(first_line, "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds"))
    {
      return FileTypes::PSMS;
    }

    // EDTA file
    // hard to tell... so we don't even try...

    return FileTypes::UNKNOWN;
  }

  PeakFileOptions& FileHandler::getOptions()
  {
    return options_;
  }

  const PeakFileOptions& FileHandler::getOptions() const
  {
    return options_;
  }

  void FileHandler::setOptions(const PeakFileOptions& options)
  {
    options_ = options;
  }

  FeatureFileOptions& FileHandler::getFeatOptions()
  {
    return f_options_;
  }

  const FeatureFileOptions& FileHandler::getFeatOptions() const
  {
    return f_options_;
  }

  void FileHandler::setFeatOptions(const FeatureFileOptions& f_options)
  {
    f_options_ = f_options;
  }

  std::string FileHandler::computeFileHash(const std::string& filename)
  {
    std::ifstream file{std::filesystem::path{std::string(filename)}, std::ios::binary};
    if (!file.is_open())
    {
      return "";
    }
    SHA1 sha;
    char buffer[8192];
    while (file.read(buffer, sizeof(buffer)) || file.gcount() > 0)
    {
      sha.update(buffer, static_cast<std::size_t>(file.gcount()));
    }
    uint32_t digest[5];
    sha.finalize(digest);

    std::ostringstream result;
    for (int i = 0; i < 5; ++i)
    {
      result << std::hex << std::setfill('0') << std::setw(8) << digest[i];
    }
    return result.str();
  }

  void FileHandler::loadSpectrum(const std::string& filename, MSSpectrum& spec, const std::vector<FileTypes::Type> allowed_types)
  {
    // determine file type
    FileTypes::Type type = getType(filename);
    // If we have a restricted set of file types check that we match them
    if (!allowed_types.empty())
    {    
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for loading a spectrum. Allowed types are: " + allowedToString_(allowed_types));
      }
    }
    switch (type)
    {
      case FileTypes::DTA: 
      {
        DTAFile().load(filename, spec);
      }
      break;

      case FileTypes::XMASS: 
      {
        XMassFile().load(filename, spec);
      }
      break;

      default: 
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) +" is not supported for loading a spectrum");
      }
    }
  }

  void FileHandler::storeSpectrum(const std::string& filename, MSSpectrum& spec, const std::vector<FileTypes::Type> allowed_types)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for storing an spectrum. Allowed types are: " + allowedToString_(allowed_types));
      }
    }
    switch (type)
    {
      case FileTypes::DTA: 
      {
        DTAFile().store(filename, spec);
      }
      break;

      case FileTypes::XMASS: 
      {
        XMassFile().store(filename, spec);
      }
      break;
      
      default: 
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type is not supported for loading experiments");
      }
    }
  }

  void FileHandler::loadExperiment(const std::string& filename, PeakMap& exp, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log, const bool rewrite_source_file,
                                   const bool compute_hash)
  {
    // setting the flag for hash recomputation only works if source file entries are rewritten
    OPENMS_PRECONDITION(rewrite_source_file || !compute_hash, "Can't compute hash if no SourceFile written");

    // determine file type
    FileTypes::Type type = getType(filename);
    // If we have a restricted set of file types check that we match them
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for loading an experiment. Allowed types are: " + allowedToString_(allowed_types));
      }
    }

    // load right file
    switch (type)
    {
      case FileTypes::DTA: 
      {
        exp.reset();
        exp.resize(1);
        DTAFile().load(filename, exp[0]);
      }
      break;

      case FileTypes::DTA2D: 
      {
        DTA2DFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        f.load(filename, exp);
      }
      break;

      case FileTypes::MZXML: 
      {
        MzXMLFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        f.load(filename, exp);
      }
      break;

      case FileTypes::MZDATA: 
      {
        MzDataFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        f.load(filename, exp);
      }
      break;

      case FileTypes::MZML: 
      {
        MzMLFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        f.load(filename, exp);
        ChromatogramTools().convertSpectraToChromatograms<PeakMap>(exp, true);
      }
      break;

      case FileTypes::MGF: 
      {
        MascotGenericFile f;
        f.setLogType(log);
        f.load(filename, exp);
      }
      break;

      case FileTypes::MS2: 
      {
        MS2File f;
        f.setLogType(log);
        f.load(filename, exp);
      }
      break;

      case FileTypes::SQMASS: 
      {
        SqMassFile().load(filename, exp);
      }
      break;

      case FileTypes::XMASS: 
      {
        exp.reset();
        exp.resize(1);
        XMassFile().load(filename, exp[0]);
        XMassFile().importExperimentalSettings(filename, exp);
      }
      break;

      case FileTypes::MSP:
      {
        MSPGenericFile().load(filename, exp);
      }
      break;

#ifdef WITH_OPENTIMS
      case FileTypes::BRUKER_TDF:
      {
        // If the input is a .d.zip archive, extract to a temp directory first.
        std::unique_ptr<File::TempDir> temp_dir;
        std::string load_path = filename;
        if (!File::isDirectory(filename) && StringUtils::hasSuffix(StringUtils::toLowered(filename), ".zip"))
        {
          load_path = ZipArchiveFile::unzipDirectory(filename, temp_dir);
          // Find the .d directory inside the extracted archive (may be nested)
          bool found_d = false;
          for (const auto& entry : std::filesystem::recursive_directory_iterator(std::string(load_path)))
          {
            if (entry.is_directory() && entry.path().extension() == ".d")
            {
              load_path = entry.path().string();
              found_d = true;
              break;
            }
          }
          if (!found_d)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              filename, "ZIP archive does not contain a .d directory");
          }
        }
        BrukerTimsFile f;
        f.setLogType(log);
        f.load(load_path, exp);

        // BrukerTimsFile loads everything; apply PeakFileOptions filters post-load.
        applyPostLoadOptions_(exp);
      }
      break;
#endif

#ifdef WITH_THERMO_RAW
      case FileTypes::RAW:
      {
        ThermoRawFile f;
        f.setLogType(log);
        f.load(filename, exp);

        // ThermoRawFile loads everything; apply PeakFileOptions filters post-load.
        applyPostLoadOptions_(exp);
      }
      break;
#endif

      default:
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type is not supported for loading experiments");
      }
    }
    if (rewrite_source_file)
    {
      SourceFile src_file;
      if (exp.getSourceFiles().empty()) // copy settings like native ID format
      {
        OPENMS_LOG_WARN << "No source file annotated." << endl;
      }
      else
      {
        if (exp.getSourceFiles().size() > 1)
        {
          OPENMS_LOG_WARN << "Expecting a single source file in mzML. Found " << exp.getSourceFiles().size() << " will take only first one for rewriting." << endl;
        }
        src_file = exp.getSourceFiles()[0];
      }

      // Normalize directory paths: strip trailing slashes for basename and URI
      std::string normalized_name = filename;
      while (StringUtils::hasSuffix(normalized_name, "/") || StringUtils::hasSuffix(normalized_name, "\\"))
        normalized_name = StringUtils::prefix(normalized_name, normalized_name.size() - 1);
      src_file.setNameOfFile(File::basename(normalized_name));
      std::string path_to_file = File::path(File::absolutePath(normalized_name)); // convert to absolute path and strip file name

      // make sure we end up with at most 3 forward slashes
      std::string uri = StringUtils::hasPrefix(path_to_file, "/") ? std::string("file://") + path_to_file : std::string("file:///") + path_to_file;
      src_file.setPathToFile(uri);
      // this is more complicated since the data formats allowed by mzML are very verbose.
      // this is prone to changing CV's... our writer will fall back to a default if the name given here is invalid.
      src_file.setFileType(FileTypes::typeToMZML(type));

      if (compute_hash && type != FileTypes::BRUKER_TDF)
      {
        src_file.setChecksum(computeFileHash(filename), SourceFile::ChecksumType::SHA1);
      }

      exp.getSourceFiles().clear();
      exp.getSourceFiles().push_back(src_file);
    }
  }

  void FileHandler::applyPostLoadOptions_(MSExperiment& exp) const
  {
    // Filter by MS level
    if (options_.hasMSLevels())
    {
      exp.getSpectra().erase(
        std::remove_if(exp.getSpectra().begin(), exp.getSpectra().end(),
          [this](const MSSpectrum& s) { return !options_.containsMSLevel(s.getMSLevel()); }),
        exp.getSpectra().end());
    }

    // Filter by RT range
    if (options_.hasRTRange())
    {
      const auto& rt_range = options_.getRTRange();
      exp.getSpectra().erase(
        std::remove_if(exp.getSpectra().begin(), exp.getSpectra().end(),
          [&rt_range](const MSSpectrum& s) { return !rt_range.encloses(DPosition<1>(s.getRT())); }),
        exp.getSpectra().end());
    }

    // Filter by precursor m/z range (for MSn spectra)
    if (options_.hasPrecursorMZRange())
    {
      const auto& prec_range = options_.getPrecursorMZRange();
      exp.getSpectra().erase(
        std::remove_if(exp.getSpectra().begin(), exp.getSpectra().end(),
          [&prec_range](const MSSpectrum& s)
          {
            if (s.getMSLevel() <= 1) { return false; } // keep MS1
            if (s.getPrecursors().empty()) { return false; }
            return !prec_range.encloses(DPosition<1>(s.getPrecursors()[0].getMZ()));
          }),
        exp.getSpectra().end());
    }

    // Filter peaks by m/z range and intensity range within each spectrum.
    // Use MSSpectrum::select() rather than erase(): select() rebuilds the parallel data
    // arrays so the per-peak ion-mobility FloatDataArray carried by Bruker .d frames (and any
    // other float/string/integer array) stays aligned with the peaks. A raw erase() would
    // shrink only the peak vector, leaving the IM array too long -> corrupt mzML output or a
    // Precondition throw on the next sortByPosition()/select().
    if (options_.hasMZRange() || options_.hasIntensityRange())
    {
      const bool filter_mz = options_.hasMZRange();
      const bool filter_int = options_.hasIntensityRange();
      const auto& mz_range = options_.getMZRange();
      const auto& int_range = options_.getIntensityRange();

      std::vector<Size> kept;
      for (auto& spectrum : exp.getSpectra())
      {
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

    // Metadata only: strip peak data
    if (options_.getMetadataOnly())
    {
      for (auto& spectrum : exp.getSpectra())
      {
        spectrum.clear(false); // clear data but keep metadata
      }
    }

    exp.updateRanges();
  }

  void FileHandler::storeExperiment(const std::string& filename, const PeakMap& exp, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for storing an experiment. Allowed types are: " + allowedToString_(allowed_types));
      }
    }

    // load right file
    switch (type)
    {
      case FileTypes::DTA: 
      {
        DTAFile().store(filename, exp[0]);
      }
      break;

      case FileTypes::DTA2D: 
      {
        DTA2DFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        f.store(filename, exp);
      }
      break;

      case FileTypes::MGF: 
      {
        MascotGenericFile f;
        f.setLogType(log);
        f.store(filename, exp);
      }
      break;

      case FileTypes::MSP: 
      {
        MSPGenericFile f;
        // TODO add support for parameters
        f.store(filename, exp);
      }
      break;

      case FileTypes::MZXML: 
      {
        MzXMLFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        if (!exp.getChromatograms().empty())
        {
          PeakMap exp2 = exp;
          ChromatogramTools().convertChromatogramsToSpectra<PeakMap>(exp2);
          f.store(filename, exp2);
        }
        else
        {
          f.store(filename, exp);
        }
      }
      break;

      case FileTypes::SQMASS: 
      {
        SqMassFile f;
        // f.setConfig()
        f.store(filename, exp);
      }
      break;

      case FileTypes::MZDATA: 
      {
        MzDataFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        if (!exp.getChromatograms().empty())
        {
          PeakMap exp2 = exp;
          ChromatogramTools().convertChromatogramsToSpectra<PeakMap>(exp2);
          f.store(filename, exp2);
        }
        else
        {
          f.store(filename, exp);
        }
      }
      break;

      case FileTypes::MZML: 
      {
        MzMLFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        f.store(filename, exp);
      }
      break;

      default: 
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for storing experiments");
      }
    }
  }

  void FileHandler::loadFeatures(const std::string& filename, FeatureMap& map, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    // determine file type
    FileTypes::Type type = getType(filename);
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for loading features. Allowed types are: " + allowedToString_(allowed_types));
      }
    }

    // load right file
    switch (type)
    {
      case FileTypes::FEATUREXML: 
      {
        FeatureXMLFile f;
        f.setLogType(log);
        f.getOptions() = f_options_;
        f.load(filename, map);
      }
      break;

      case FileTypes::TSV: 
      {
        MsInspectFile().load(filename, map);
      }
      break;

      case FileTypes::PEPLIST: 
      {
        SpecArrayFile().load(filename, map);
      }
      break;

      case FileTypes::KROENIK: 
      {
        KroenikFile().load(filename, map);
      }
      break;

      case FileTypes::OMS:
      {
        OMSFile f;
        f.setLogType(log);
        f.load(filename, map);
      }
      break;

      case FileTypes::FEATUREPARQUET:
      {
        if (!FeatureMapArrowIO::importFromParquet(filename, map))
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                      "FeatureMapArrowIO::importFromParquet failed");
        }
      }
      break;

      default:
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,"type: " + FileTypes::typeToName(type) + " is not supported for loading features");
      }
    }
  }

  void FileHandler::storeFeatures(const std::string& filename, const FeatureMap& map, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }

    // If we have a restricted set of file types check that we match them
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for storing features. Allowed types are: " + allowedToString_(allowed_types));
      }
    }

    //store right file
    switch (type)
    {
      case FileTypes::FEATUREXML:
      {
        FeatureXMLFile f;
        f.setLogType(log);
        f.getOptions() = f_options_;
        f.store(filename, map);
      }
      break;

      case FileTypes::EDTA:
      {
        EDTAFile f;
        f.store(filename, map);
      }
      break;

      case FileTypes::TSV:
      {
        MsInspectFile().store(filename, map);
      }
      break;

      case FileTypes::OMS:
      {
        OMSFile f;
        f.setLogType(log);
        f.store(filename, map);
      }
      break;

      case FileTypes::PEPLIST:
      {
        SpecArrayFile().store(filename, map);
      }
      break;

      case FileTypes::KROENIK:
      {
        KroenikFile().store(filename, map);
      }
      break;

      case FileTypes::FEATUREPARQUET:
      {
        if (!FeatureMapArrowIO::exportToParquet(map, filename))
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                              "FeatureMapArrowIO::exportToParquet failed");
        }
      }
      break;

      default:
      {
          throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for storing features");
      }
    }
  }

  void FileHandler::loadConsensusFeatures(const std::string& filename, ConsensusMap& map, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {

    //determine file type
    FileTypes::Type type = getType(filename);

    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for loading consensus features, Allowed types are: " + allowedToString_(allowed_types));
      }
    }

    switch (type)
    {
      case FileTypes::CONSENSUSXML:
      {
        ConsensusXMLFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        f.load(filename, map);
      }
      break;

      case FileTypes::EDTA:
      {
        EDTAFile f;
        f.load(filename, map);
      }
      break;

      case FileTypes::OMS:
      {
        OMSFile f;
        f.setLogType(log);
        f.load(filename, map);
      }
      break;

      case FileTypes::CONSENSUSPARQUET:
      {
        if (!ConsensusMapArrowIO::importFromParquet(filename, map))
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                      "ConsensusMapArrowIO::importFromParquet failed");
        }
      }
      break;

      default:
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not  supported for loading consensus features");
      }
    }
  }

  void FileHandler::storeConsensusFeatures(const std::string& filename, const ConsensusMap& map,  const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for storing an Consensus Features. Allowed types are: " + allowedToString_(allowed_types));
      }
    }
    switch (type)
    {
      case FileTypes::CONSENSUSXML:
      {
        ConsensusXMLFile f;
        f.setLogType(log);
        f.store(filename, map);
      }
      break;

      case FileTypes::EDTA:
      {
        EDTAFile f;
        f.store(filename, map);
      }
      break;

      case FileTypes::OMS:
      {
        OMSFile f;
        f.setLogType(log);
        f.store(filename, map);
      }
      break;

      case FileTypes::CONSENSUSPARQUET:
      {
        if (!ConsensusMapArrowIO::exportToParquet(map, filename))
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                              "ConsensusMapArrowIO::exportToParquet failed");
        }
      }
      break;

      default:
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for storing consensus features");
      }
    }
  }

  void FileHandler::loadIdentifications(const std::string& filename, std::vector<ProteinIdentification>& additional_proteins, PeptideIdentificationList& additional_peptides, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    
    //determine file type
    FileTypes::Type type = getType(filename);
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for loading identifications, Allowed types are: " + allowedToString_(allowed_types));
      }
    }

    switch (type)
    {
      case FileTypes::IDXML:
      {
        IdXMLFile f;
        f.setLogType(log);
        f.load(filename, additional_proteins, additional_peptides);
      }
      break;

      case FileTypes::MZIDENTML:
      {
        MzIdentMLFile f;
        f.setLogType(log);
        f.load(filename, additional_proteins, additional_peptides);
      }
      break;

      case FileTypes::OMS:
      {
        OMSFile f;
        f.setLogType(log);
        IdentificationData idd;
        f.load(filename, idd);
        IdentificationDataConverter::exportIDs(idd, additional_proteins, additional_peptides);
      }
      break;

      case FileTypes::XQUESTXML:
      {
      XQuestResultXMLFile f;
      f.setLogType(log);
      f.load(filename, additional_peptides, additional_proteins);
      }
      break;


      case FileTypes::OMSSAXML:
      {
        additional_proteins.push_back(ProteinIdentification());
        OMSSAXMLFile().load(filename, additional_proteins[0],
                            additional_peptides, true, true);
      }
      break;
      
      /*case FileTypes::MASCOTXML:
      {
        OPENMS_LOG_ERROR << "File " << filename << " Loading Identifications is not yet supported for MASCOTXML files" << endl;
        return false;
      }*/

      case FileTypes::PROTXML:
      {
        additional_proteins.push_back(ProteinIdentification());
        additional_peptides.push_back(PeptideIdentification());
        ProtXMLFile().load(filename, additional_proteins.back(), additional_peptides.back());
      }
      break;

      case FileTypes::IDPARQUET:
      {
        if (!PSMArrowIO::importFromParquet(filename, additional_proteins, additional_peptides))
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                      "PSMArrowIO::importFromParquet failed");
        }
      }
      break;

      default:
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for loading identifications");
    }   
  }

  void FileHandler::storeIdentifications(const std::string& filename, const std::vector<ProteinIdentification>& additional_proteins, const PeptideIdentificationList& additional_peptides, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for storing identifications. Allowed types are: " + allowedToString_(allowed_types));
      }
    }

    switch (type)
    {
      case FileTypes::IDXML:
      {
        IdXMLFile f;
        f.setLogType(log);
        f.store(filename, additional_proteins, additional_peptides);
      }
      break;

      case FileTypes::MZIDENTML:
      {
        MzIdentMLFile f;
        f.setLogType(log);
        f.store(filename, additional_proteins, additional_peptides);
      }
      break;

      case FileTypes::OMS:
      {
        OMSFile f;
        f.setLogType(log);
        IdentificationData idd;
        IdentificationDataConverter::importIDs(idd, additional_proteins, additional_peptides);
        f.store(filename, idd);
      }
    break;

      case FileTypes::XQUESTXML:
      {
      XQuestResultXMLFile f;
      f.setLogType(log);
      f.store(filename, additional_proteins, additional_peptides);
      }
      break;

      case FileTypes::IDPARQUET:
      {
        if (!PSMArrowIO::exportToParquet(additional_proteins, additional_peptides, filename))
        {
          throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,
                                              "PSMArrowIO::exportToParquet failed");
        }
      }
      break;

      default:
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for storing Identifications");
      }
    }
  }

  void FileHandler::loadTransitions(const std::string& filename,TargetedExperiment& library, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    //determine file type
    FileTypes::Type type = getType(filename);
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for loading transitions, Allowed types are: " + allowedToString_(allowed_types));
      }
    }

    switch (type)
    {
      case FileTypes::TRAML:
      {
        TraMLFile f;
        f.setLogType(log);
        f.load(filename, library);
      }
      break;

      default:
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,"type: " + FileTypes::typeToName(type) + " is not supported for loading transitions");
      }
    }
  }

  void FileHandler::storeTransitions(const std::string& filename, const TargetedExperiment& library, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for storing transitions. Allowed types are: " + allowedToString_(allowed_types));
      }
    }
    switch (type)
    {
      case FileTypes::TRAML:
      {
        TraMLFile f;
        f.setLogType(log);
        f.store(filename, library);
      }
      break;

      default:
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for storing transitions"); 
      }
    }
  }

  void FileHandler::loadTransformations(const std::string& filename, TransformationDescription& map, bool fit_model, const std::vector<FileTypes::Type> allowed_types)
  {
    //determine file type
    FileTypes::Type type = getType(filename);
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for loading transformations, Allowed types are: " + allowedToString_(allowed_types));
      }
    }
    switch (type)
    {
      case FileTypes::TRANSFORMATIONXML:
      {
        TransformationXMLFile().load(filename, map, fit_model);
      }
      break;
      
      default:
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,"type: " + FileTypes::typeToName(type) + " is not supported for loading transformations");
      }
    }
  }

  void FileHandler::storeTransformations(const std::string& filename, const TransformationDescription& map,  const std::vector<FileTypes::Type> allowed_types)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for storing transformations. Allowed types are: " + allowedToString_(allowed_types));
      }
    }
    
    switch (type)
    {
      case FileTypes::TRANSFORMATIONXML:
      {
        TransformationXMLFile().store(filename, map);
      }
      break;

      default:
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for storing transformations");
      }
    }
  }

  void FileHandler::storeQC(const std::string& input_file,
               const std::string& filename,
               const MSExperiment& exp,
               const FeatureMap& feature_map,
               std::vector<ProteinIdentification>& prot_ids,
               PeptideIdentificationList& pep_ids,
               const ConsensusMap& consensus_map,
               const std::string& contact_name,
               const std::string& contact_address,
               const std::string& description,
               const std::string& label,
               const bool remove_duplicate_features,
               const std::vector<FileTypes::Type> allowed_types
             )
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (!allowed_types.empty())
    {
      if (!FileTypeList(allowed_types).contains(type))
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not allowed for storing QC data. Allowed types are: " + allowedToString_(allowed_types));
      }
    }

    switch (type)
    {
      case FileTypes::QCML:
      {
      QcMLFile qcmlfile;
      qcmlfile.collectQCData(prot_ids, pep_ids, feature_map,
                    consensus_map, input_file, remove_duplicate_features, exp);
      qcmlfile.store(filename);
      }
      break;
      
      case FileTypes::MZQC:
      {
      MzQCFile().store(input_file, filename, exp, contact_name, contact_address,
                     description, label, feature_map, prot_ids, pep_ids);
      }
      break;
      
      default:
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for storing QC data");
      }
    }
  }

} // namespace OpenMS
