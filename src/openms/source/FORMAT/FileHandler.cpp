// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>

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

#include <OpenMS/FORMAT/MsInspectFile.h>
#include <OpenMS/FORMAT/SpecArrayFile.h>
#include <OpenMS/FORMAT/KroenikFile.h>

#include <OpenMS/KERNEL/ChromatogramTools.h>

#include <OpenMS/FORMAT/GzipIfstream.h>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>

#include <QtCore/QFile>
#include <QtCore/QCryptographicHash>

using namespace std;

namespace OpenMS
{
  String allowedToString_(vector<FileTypes::Type> types)
  {
    String aStrings;
    for (auto i : types)
    {
      if (i != FileTypes::SIZE_OF_TYPE)
      {
        aStrings +=  ", " + FileTypes::typeToName(i);
      }
    }
    return aStrings;
  }

  FileTypes::Type FileHandler::getType(const String& filename)
  {
    FileTypes::Type type = getTypeByFileName(filename);
    if (type == FileTypes::UNKNOWN)
    {
      type = getTypeByContent(filename);
    }
    return type;
  }

  FileTypes::Type FileHandler::getTypeByFileName(const String& filename)
  {
    String basename = File::basename(filename), tmp;
    // special rules for "double extensions":
    if (basename.hasSuffix(".pep.xml"))
    {
      return FileTypes::PEPXML;
    }
    if (basename.hasSuffix(".prot.xml"))
    {
      return FileTypes::PROTXML;
    }
    if (basename.hasSuffix(".xquest.xml"))
    {
      return FileTypes::XQUESTXML;
    }
    if (basename.hasSuffix(".spec.xml"))
    {
      return FileTypes::SPECXML;
    }
    try
    {
      tmp = basename.suffix('.');
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
    tmp.toUpper();
    if (tmp == "BZ2" || tmp == "GZ") // todo ZIP (not supported yet):       || tmp == "ZIP"
    {
      // do not use getTypeByContent() here, as this is deadly for output files!
      return getTypeByFileName(filename.prefix(filename.size() - tmp.size() - 1)); // check name without compression suffix (e.g. bla.mzML.gz --> bla.mzML)
    }

    return FileTypes::nameToType(tmp);
  }

  bool FileHandler::hasValidExtension(const String& filename, const FileTypes::Type type)
  {
    FileTypes::Type ft = FileHandler::getTypeByFileName(filename);
    return (ft == type || ft == FileTypes::UNKNOWN);
  }

  String FileHandler::stripExtension(const String& filename)
  {
    if (!filename.has('.'))
    {
      return filename;
    }
    // we don't just search for the last '.' and remove the suffix, because this could be wrong, e.g. bla.mzML.gz would become bla.mzML
    auto type = getTypeByFileName(filename);
    auto s_type = FileTypes::typeToName(type);
    size_t pos = String(filename).toLower().rfind(s_type.toLower()); // search backwards in entire string, because we could search for 'mzML' and have 'mzML.gz'
    if (pos == string::npos) // file type was FileTypes::UNKNOWN and we did not find '.unknown' as ending
    {
      size_t ext_pos = filename.rfind('.');
      size_t dir_sep = filename.find_last_of("/\\"); // look for '/' or '\'
      if (dir_sep != string::npos && dir_sep > ext_pos) // we found a directory separator after the last '.', e.g. '/my.dotted.dir/filename'! Ouch!
      { // do not strip anything, because there is no extension to strip
        return filename;
      }
      return filename.prefix(ext_pos);
    }
    return filename.prefix(pos - 1); // strip the '.' as well
  }

  String FileHandler::swapExtension(const String& filename, const FileTypes::Type new_type)
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

  FileTypes::Type FileHandler::getConsistentOutputfileType(const String& output_filename, const String& requested_type)
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

  FileTypes::Type FileHandler::getTypeByContent(const String& filename)
  {
    String first_line;
    String two_five;
    String all_simple;

    // only the first five lines will be set for compressed files
    // so far, compression is only supported for XML files
    vector<String> complete_file;

    // test whether the file is compressed (bzip2 or gzip)
    ifstream compressed_file(filename.c_str());
    char bz[2];
    compressed_file.read(bz, 2);
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
      String buffer_str(buffer);
      vector<String> split;
      buffer_str.split('\n', split);
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
      String buffer_str(buffer);
      vector<String> split;
      buffer_str.split('\n', split);
      split.resize(5);

      first_line = split[0];
      two_five = split[1] + ' ' + split[2] + ' ' + split[3] + ' ' + split[4];
      all_simple = first_line + ' ' + two_five;
      complete_file = split;
    }
    //else {} // TODO: ZIP
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
        two_five = two_five.chop(1);
        two_five.substitute('\t', ' ');
        all_simple = *(file.begin()) + ' ' + two_five;
        first_line = *(file.begin());
      }

      complete_file.insert(complete_file.end(), file.begin(), file.end());
    }
    //std::cerr << "\n Line1:\n" << first_line << "\nLine2-5:\n" << two_five << "\nall:\n" << all_simple << "\n\n";


    //mzXML (all lines)
    if (all_simple.hasSubstring("<mzXML"))
    {
      return FileTypes::MZXML;
    }
    //mzData (all lines)
    if (all_simple.hasSubstring("<mzData"))
    { 
      return FileTypes::MZDATA;
    }
    //mzML (all lines)
    if (all_simple.hasSubstring("<mzML"))
    {
      return FileTypes::MZML;
    }
    //"analysisXML" aka. mzid (all lines)
    if (all_simple.hasSubstring("<MzIdentML"))
    {
      return FileTypes::MZIDENTML;
    }
    //subject to change!
    if (all_simple.hasSubstring("<MzQualityMLType"))
    {
      return FileTypes::QCML;
    }
    //pepXML (all lines)
    if (all_simple.hasSubstring("xmlns=\"http://regis-web.systemsbiology.net/pepXML\""))
    {
      return FileTypes::PEPXML;
    }
    //protXML (all lines)
    if (all_simple.hasSubstring("xmlns=\"http://regis-web.systemsbiology.net/protXML\""))
    {
      return FileTypes::PROTXML;
    }
    //feature map (all lines)
    if (all_simple.hasSubstring("<featureMap"))
    {
      return FileTypes::FEATUREXML;
    }
    //idXML (all lines)
    if (all_simple.hasSubstring("<IdXML"))
    {
      return FileTypes::IDXML;
    }
    //consensusXML (all lines)
    if (all_simple.hasSubstring("<consensusXML"))
    {
      return FileTypes::CONSENSUSXML;
    }
    //TOPPAS (all lines)
    if (all_simple.hasSubstring("<PARAMETERS") && all_simple.hasSubstring("<NODE name=\"info\"") && all_simple.hasSubstring("<ITEM name=\"num_vertices\""))
    {
      return FileTypes::TOPPAS;
    }
    //INI (all lines) (must be AFTER TOPPAS) - as this is less restrictive
    if (all_simple.hasSubstring("<PARAMETERS"))
    {
      return FileTypes::INI;
    }
    //TrafoXML (all lines)
    if (all_simple.hasSubstring("<TrafoXML"))
    {
      return FileTypes::TRANSFORMATIONXML;
    }
    //GelML (all lines)
    if (all_simple.hasSubstring("<GelML"))
    {
      return FileTypes::GELML;
    }
    //traML (all lines)
    if (all_simple.hasSubstring("<TraML"))
    {
      return FileTypes::TRAML;
    }
    //OMSSAXML file
    if (all_simple.hasSubstring("<MSResponse"))
    {
      return FileTypes::OMSSAXML;
    }
    //MASCOTXML file
    if (all_simple.hasSubstring("<mascot_search_results"))
    {
      return FileTypes::MASCOTXML;
    }
    if (all_simple.hasPrefix("{"))
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
        if (complete_file[i].trim().hasPrefix(">"))
        {
          ++bigger_than;
          ++i;
        }
        else if (complete_file[i].trim().hasPrefix("#"))
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
    if (first_line.substr(1, 3) == "PNG")
    {
      return FileTypes::PNG;
    }
    //MSP (all lines)
    for (Size i = 0; i != complete_file.size(); ++i)
    {
      if (complete_file[i].hasPrefix("Name: ") && complete_file[i].hasSubstring("/"))
      {
        return FileTypes::MSP;
      }
      if (complete_file[i].hasPrefix("Num peaks: "))
      {
        return FileTypes::MSP;
      }
    }

    //tokenize lines 2-5
    vector<String> parts;
    two_five.split(' ', parts);

    //DTA
    if (parts.size() == 8)
    {
      bool conversion_error = false;
      try
      {
        for (Size i = 0; i < 8; ++i)
        {
          parts[i].toFloat();
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
          parts[i].toFloat();
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
    if (two_five.hasSubstring("BEGIN IONS"))
    {
      return FileTypes::MGF;
    }
    else
    {
      for (Size i = 0; i != complete_file.size(); ++i)
      {
        if (complete_file[i].trim() == "FORMAT=Mascot generic" || complete_file[i].trim() == "BEGIN IONS")
        {
          return FileTypes::MGF;
        }
      }
    }

    // MS2 file format
    if (all_simple.hasSubstring("CreationDate"))
    {
      if (!all_simple.empty() && all_simple[0] == 'H')
      {
        return FileTypes::MS2;
      }
    }

    // mzTab file format
    for (Size i = 0; i != complete_file.size(); ++i) {
        if (complete_file[i].hasSubstring("MTD\tmzTab-version")) {
            return FileTypes::MZTAB;
        }
    }

    // msInspect file (.tsv)
    for (Size i = 0; i != complete_file.size(); ++i)
    {
      if (complete_file[i].hasSubstring("scan\ttime\tmz\taccurateMZ\tmass\tintensity\tcharge\tchargeStates\tkl\tbackground\tmedian\tpeaks\tscanFirst\tscanLast\tscanCount\ttotalIntensity\tsumSquaresDist\tdescription"))
      {
        return FileTypes::TSV;
      }
    }

    // specArray file (.pepList)
    if (first_line.hasSubstring("       m/z\t     rt(min)\t       snr\t      charge\t   intensity"))
    {
      return FileTypes::PEPLIST;
    }

    // hardkloer file (.hardkloer)
    /**
      NOT IMPLEMENTED YET
      if (first_line.hasSubstring("File	First Scan	Last Scan	Num of Scans	Charge	Monoisotopic Mass	Base Isotope Peak	Best Intensity	Summed Intensity	First RTime	Last RTime	Best RTime	Best Correlation	Modifications"))
    {
        return FileTypes::HARDKLOER;
    }
    **/

    // kroenik file (.kroenik)
    if (first_line.hasSubstring("File\tFirst Scan\tLast Scan\tNum of Scans\tCharge\tMonoisotopic Mass\tBase Isotope Peak\tBest Intensity\tSummed Intensity\tFirst RTime\tLast RTime\tBest RTime\tBest Correlation\tModifications"))
    {
      return FileTypes::KROENIK;
    }

    // Percolator tab-delimited output (PSM level, .psms)
    if (first_line.hasPrefix("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds"))
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

  String FileHandler::computeFileHash(const String& filename)
  {
    QCryptographicHash crypto(QCryptographicHash::Sha1);
    QFile file(filename.toQString());
    file.open(QFile::ReadOnly);
    while (!file.atEnd())
    {
      crypto.addData(file.read(8192));
    }
    return String((QString)crypto.result().toHex());
  }

  void FileHandler::loadSpectrum(const String& filename, MSSpectrum& spec, const std::vector<FileTypes::Type> allowed_types)
  {
    // determine file type
    FileTypes::Type type = getType(filename);
    // If we have a restricted set of file types check that we match them
    if (allowed_types.size() != 0)
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

  void FileHandler::storeSpectrum(const String& filename, MSSpectrum& spec, const std::vector<FileTypes::Type> allowed_types)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (allowed_types.size() != 0)
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

  void FileHandler::loadExperiment(const String& filename, PeakMap& exp, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log, const bool rewrite_source_file,
                                   const bool compute_hash)
  {
    // setting the flag for hash recomputation only works if source file entries are rewritten
    OPENMS_PRECONDITION(rewrite_source_file || !compute_hash, "Can't compute hash if no SourceFile written");

    // determine file type
    FileTypes::Type type = getType(filename);
    // If we have a restricted set of file types check that we match them
    if (allowed_types.size() != 0)
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

      src_file.setNameOfFile(File::basename(filename));
      String path_to_file = File::path(File::absolutePath(filename)); // convert to absolute path and strip file name

      // make sure we end up with at most 3 forward slashes
      String uri = path_to_file.hasPrefix("/") ? String("file://") + path_to_file : String("file:///") + path_to_file;
      src_file.setPathToFile(uri);
      // this is more complicated since the data formats allowed by mzML are very verbose.
      // this is prone to changing CV's... our writer will fall back to a default if the name given here is invalid.
      src_file.setFileType(FileTypes::typeToMZML(type));

      if (compute_hash)
      {
        src_file.setChecksum(computeFileHash(filename), SourceFile::SHA1);
      }

      exp.getSourceFiles().clear();
      exp.getSourceFiles().push_back(src_file);
    }
  }

  void FileHandler::storeExperiment(const String& filename, const PeakMap& exp, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (allowed_types.size() != 0)
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

  void FileHandler::loadFeatures(const String& filename, FeatureMap& map, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    // determine file type
    FileTypes::Type type = getType(filename);
    if (allowed_types.size() != 0)
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

      default: 
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename,"type: " + FileTypes::typeToName(type) + " is not supported for loading features");
      }
    }
  }

  void FileHandler::storeFeatures(const String& filename, const FeatureMap& map, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }

    // If we have a restricted set of file types check that we match them
    if (allowed_types.size() != 0)
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

      default:
      {
          throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for storing features");
      }
    }
  }

  void FileHandler::loadConsensusFeatures(const String& filename, ConsensusMap& map, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {

    //determine file type
    FileTypes::Type type = getType(filename);

    if (allowed_types.size() != 0)
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

      default:
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not  supported for loading consensus features");
      }
    }
  }

  void FileHandler::storeConsensusFeatures(const String& filename, const ConsensusMap& map,  const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (allowed_types.size() != 0)
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
      
      default:
      {        
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for storing consensus features");
      }
    }
  }

  void FileHandler::loadIdentifications(const String& filename, std::vector<ProteinIdentification>& additional_proteins, std::vector<PeptideIdentification>& additional_peptides, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    
    //determine file type
    FileTypes::Type type = getType(filename);
    if (allowed_types.size() != 0)
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

      default:
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for loading identifications");
    }   
  }

  void FileHandler::storeIdentifications(const String& filename, const std::vector<ProteinIdentification>& additional_proteins, const std::vector<PeptideIdentification>& additional_peptides, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (allowed_types.size() != 0)
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

      default:
      {
        throw Exception::InvalidFileType(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "type: " + FileTypes::typeToName(type) + " is not supported for storing Identifications");
      }
    }   
  }

  void FileHandler::loadTransitions(const String& filename,TargetedExperiment& library, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    //determine file type
    FileTypes::Type type = getType(filename);
    if (allowed_types.size() != 0)
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

  void FileHandler::storeTransitions(const String& filename, const TargetedExperiment& library, const std::vector<FileTypes::Type> allowed_types, ProgressLogger::LogType log)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (allowed_types.size() != 0)
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

  void FileHandler::loadTransformations(const String& filename, TransformationDescription& map, bool fit_model, const std::vector<FileTypes::Type> allowed_types)
  {
    //determine file type
    FileTypes::Type type = getType(filename);
    if (allowed_types.size() != 0)
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

  void FileHandler::storeTransformations(const String& filename, const TransformationDescription& map,  const std::vector<FileTypes::Type> allowed_types)
  {
    auto type = getTypeByFileName(filename);
    if (type == FileTypes::Type::UNKNOWN && (allowed_types.size() == 1))
    { // filename is unspecific, but allowed_types is unambiguous (i.e. they do not contradict)
      type = allowed_types[0];
    }
    // If we have a restricted set of file types check that we match them
    if (allowed_types.size() != 0)
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

  void FileHandler::storeQC(const String& input_file,
               const String& filename,
               const MSExperiment& exp,
               const FeatureMap& feature_map,
               std::vector<ProteinIdentification>& prot_ids,
               std::vector<PeptideIdentification>& pep_ids,
               const ConsensusMap& consensus_map,
               const String& contact_name,
               const String& contact_address,
               const String& description,
               const String& label,
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
    if (allowed_types.size() != 0)
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
