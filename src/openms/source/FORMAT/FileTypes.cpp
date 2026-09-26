// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Andreas Bertsch, Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileTypes.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <array>
#include <cctype>
#include <list>
#include <utility>

#include <cassert>
namespace OpenMS
{

  /// connect the type to some other information
  /// We could also use paired arrays, but this way, its less likely to have mismatches if a new type is added
  struct TypeNameBinding
  {
    FileTypes::Type type;
    std::string name;
    std::string description;
    std::vector<FileTypes::FileProperties> features;
    /// additional extensions accepted for this type, besides 'name' (which stays the preferred one)
    std::vector<std::string> aliases;
    TypeNameBinding(FileTypes::Type ptype, std::string pname, std::string pdescription, std::vector<FileTypes::FileProperties> pfeatures, std::vector<std::string> paliases = {})
      : type(ptype), name(std::move(pname)), description(std::move(pdescription)), features(pfeatures), aliases(std::move(paliases))
    {
      // Check that there are no double-spaces in the description, since Qt will replace "  " with " " in filters supplied to QFileDialog::getSaveFileName.
      // And if you later ask for the selected filter, you will get a different string back.
      assert(!description.contains("  "));
    }
  };

  using PROP = FileTypes::FileProperties;   // shorten our syntax a bit
  /// Maps the FileType::Type to its preferred extension plus any additional accepted extensions.
  /// The preferred extension is what typeToName() returns and what we write; the aliases are only ever accepted on input.
  /// An alias must be unique across all types (FileTypes_test enforces this), and must not be a suffix that other formats
  /// also use (e.g. a plain 'xml' alias would make 'pepXML' and 'mzML' ambiguous).
  /// when adding new types, be sure to update the FileTypes_test typesWithProperties test to match the new files
  static const std::array<TypeNameBinding, FileTypes::SIZE_OF_TYPE> type_with_annotation__ =
  {
    TypeNameBinding(FileTypes::UNKNOWN, "unknown", "unknown file extension", {}),
    TypeNameBinding(FileTypes::DTA, "dta", "dta raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::PROVIDES_SPECTRUM, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::DTA2D, "dta2d", "dta2d raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MZDATA, "mzData", "mzData raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::MZXML, "mzXML", "mzXML raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::FEATUREXML, "featureXML", "OpenMS feature map", {PROP::PROVIDES_FEATURES, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::IDXML, "idXML", "OpenMS peptide identification file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::CONSENSUSXML, "consensusXML", "OpenMS consensus feature map", {PROP::PROVIDES_CONSENSUSFEATURES, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::MGF, "mgf", "mascot generic format file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::INI, "ini", "OpenMS parameter file", {PROP::READABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::TOPPAS, "toppas", "OpenMS TOPPAS pipeline", {PROP::READABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::TRANSFORMATIONXML, "trafoXML", "RT transformation file", {PROP::PROVIDES_TRANSFORMATIONS, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::MZML, "mzML", "mzML raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    //TODO: Add support for cachedMZML as a first class file type
    TypeNameBinding(FileTypes::CACHEDMZML, "cachedMzML", "cachedMzML raw data file", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MS2, "ms2", "ms2 file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE}),
    TypeNameBinding(FileTypes::PEPXML, "pepXML", "pepXML file", {PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}, {"pep.xml"}), //Supported for loading and storing identifications but TODO integrate this into fileHandler
    TypeNameBinding(FileTypes::PROTXML, "protXML", "protXML file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::COMPRESSED_READABLE}, {"prot.xml"}),
    TypeNameBinding(FileTypes::MZIDENTML, "mzid", "mzIdentML file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::QCML, "qcml", "quality control file", {PROP::PROVIDES_QC, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}), //TODO add load functions for QC
    TypeNameBinding(FileTypes::MZQC, "mzqc", "quality control file in json format", {PROP::PROVIDES_QC, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::GELML, "gelML", "gelML file", {}),
    TypeNameBinding(FileTypes::TRAML, "traML", "transition file", {PROP::PROVIDES_TRANSITIONS, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::MSP, "msp", "NIST spectra library file format", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::OMSSAXML, "omssaXML", "omssaXML file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::MASCOTXML, "mascotXML", "mascotXML file", {PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::PNG, "png", "portable network graphics file", {}),
    TypeNameBinding(FileTypes::XMASS, "fid", "XMass analysis file", {PROP::PROVIDES_EXPERIMENT, PROP::PROVIDES_SPECTRUM, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::TSV, "tsv", "tab-separated file", {PROP::PROVIDES_FEATURES, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MZTAB, "mzTab", "mzTab file", {}), //TODO add filehandler support for MZTAB
    TypeNameBinding(FileTypes::PEPLIST, "peplist", "SpecArray file", {PROP::PROVIDES_FEATURES, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::HARDKLOER, "hardkloer", "hardkloer file", {}),
    TypeNameBinding(FileTypes::KROENIK, "kroenik", "kroenik file", {PROP::PROVIDES_FEATURES, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::FASTA, "fasta", "FASTA file", {PROP::READABLE, PROP::WRITEABLE}, {"fa", "faa"}),
    TypeNameBinding(FileTypes::PEFF, "peff", "PEFF protein file", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::EDTA, "edta", "enhanced dta file", {PROP::PROVIDES_FEATURES, PROP::PROVIDES_CONSENSUSFEATURES, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::CSV, "csv", "comma-separated values file", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::TXT, "txt", "generic text file", {}),
    TypeNameBinding(FileTypes::OBO, "obo", "controlled vocabulary file", {}),
    TypeNameBinding(FileTypes::HTML, "html", "any HTML file", {}),
    TypeNameBinding(FileTypes::ANALYSISXML, "analysisXML", "analysisXML file", {}),
    TypeNameBinding(FileTypes::XSD, "xsd", "XSD schema format", {}),
    TypeNameBinding(FileTypes::PSQ, "psq", "NCBI binary blast db", {}),
    TypeNameBinding(FileTypes::MRM, "mrm", "SpectraST MRM list", {PROP::READABLE}),
    TypeNameBinding(FileTypes::SQMASS, "sqMass", "SQLite format for mass and chromatograms", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::PQP, "pqp", "pqp file", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::OSWPQ, "oswpq", "OpenSwath Parquet bundle", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MS, "ms", "SIRIUS file", {}),
    TypeNameBinding(FileTypes::OSW, "osw", "OpenSwath output files", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::CHROMPARQUET, "xic", "OpenSwath Parquet chromatogram output", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MOBILPARQUET, "xim", "OpenSwath Parquet mobilogram output", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::PEAKMAPPARQUET, "xipm", "OpenSwath Parquet peak-map output", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::PSMS, "psms", "Percolator tab-delimited output (PSM level)", {PROP::READABLE}),
    TypeNameBinding(FileTypes::PIN, "pin", "Percolator tab-delimited input (PSM level)", {}),
    TypeNameBinding(FileTypes::PARAMXML, "paramXML", "OpenMS internal XML file", {PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::SPLIB, "splib", "SpectraST binary spectral library file", {}),
    TypeNameBinding(FileTypes::NOVOR, "novor", "Novor custom parameter file", {}),
    TypeNameBinding(FileTypes::XQUESTXML, "xquest.xml", "xquest.xml file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::SPECXML, "spec.xml", "spec.xml file", {}),
    TypeNameBinding(FileTypes::JSON, "json", "JavaScript Object Notation file", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::RAW, "raw", "(Thermo) Raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE}),
    TypeNameBinding(FileTypes::OMS, "oms", "OpenMS SQLite file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::PROVIDES_FEATURES, PROP::PROVIDES_CONSENSUSFEATURES}),
    TypeNameBinding(FileTypes::EXE, "exe", "Windows executable", {}),
    TypeNameBinding(FileTypes::BZ2, "bz2", "bzip2 compressed file", {PROP::READABLE}),
    TypeNameBinding(FileTypes::GZ, "gz", "gzip compressed file", {PROP::READABLE}),
    TypeNameBinding(FileTypes::ZIP, "zip", "ZIP compressed file", {PROP::READABLE}),
    TypeNameBinding(FileTypes::PARQUET, "parquet", "Apache Parquet file", {PROP::READABLE, PROP::WRITEABLE}, {"pqt"}),
    TypeNameBinding(FileTypes::IDPARQUET, "idparquet", "OpenMS identification parquet bundle (directory)", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::FEATUREPARQUET, "featureparquet", "OpenMS feature map parquet bundle (directory)", {PROP::PROVIDES_FEATURES, PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::CONSENSUSPARQUET, "consensusparquet", "OpenMS consensus map parquet bundle (directory)", {PROP::PROVIDES_CONSENSUSFEATURES, PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::BRUKER_TDF, "d", "Bruker TDF", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::COMPRESSED_READABLE}), // .d.zip via ZipArchiveFile
    TypeNameBinding(FileTypes::IMZML, "imzML", "imzML mass spectrometry imaging file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE, PROP::COMPRESSED_READABLE}),
    TypeNameBinding(FileTypes::YAML, "yaml", "YAML file", {PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::XML, "xml", "any XML file", {PROP::READABLE, PROP::COMPRESSED_READABLE}),  // make sure this comes last, since the name is a suffix of other formats and should only be matched last
  };

  FileTypeList::FileTypeList(const std::vector<FileTypes::Type>& types)
    : type_list_(types)
  {
  }

  bool FileTypeList::contains(const FileTypes::Type& type) const
  {
    for (const auto& t : type_list_)
    {
      if (t == type)
      {
        return true;
      }
    }
    return false;
  }

  std::string FileTypeList::toFileDialogFilter(const FilterLayout style, bool add_all_filter) const
  {
    return ListUtils::concatenate(asFilterElements_(style, add_all_filter).items, ";;");
  }

  FileTypes::Type FileTypeList::fromFileDialogFilter(const std::string& filter, const FileTypes::Type fallback) const
  {
    auto candidates = asFilterElements_(FilterLayout::BOTH, true); // may add more filters than needed, but that's fine

    auto where = std::find(candidates.items.begin(), candidates.items.end(), filter);
    if (where == candidates.items.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filter);
    }
    const FileTypes::Type r = candidates.types[where - candidates.items.begin()];
    return r == FileTypes::Type::UNKNOWN ? fallback : r;
  }

  std::vector<FileTypes::Type> FileTypeList::typesWithProperties(std::vector<FileTypes::FileProperties> haveFeatures)
  {
    std::vector<FileTypes::Type> compatible;
    std::vector<TypeNameBinding> good_types(type_with_annotation__.begin(), type_with_annotation__.end());
    // for each feature we are looking for
    for (auto i : haveFeatures)
    {
      // Remove any types that lack the feature
      good_types.erase(std::remove_if(good_types.begin(), good_types.end(),[i](auto j) { return (std::find(j.features.begin(),j.features.end(),i) == j.features.end()); }), good_types.end());
    }
    
    for (const auto& t : good_types)
    {
      compatible.push_back(t.type);
    }

    return compatible;
  }

  
  FileTypeList::FilterElements_ FileTypeList::asFilterElements_(const FilterLayout style, bool add_all_filter) const
  {
    FilterElements_ result;

    if (style == FilterLayout::COMPACT || style == FilterLayout::BOTH)
    {
      StringList items;
      for (const auto& t : type_list_)
      {
        for (const auto& ext : FileTypes::typeToExtensions(t))
        {
          const std::string pattern = "*." + ext;
          // the same alias can only belong to one type, but a type may appear twice in type_list_
          if (!ListUtils::contains(items, pattern)) items.push_back(pattern);
        }
      }
      result.items.emplace_back("all readable files (" + ListUtils::concatenate(items, " ") + ")");
      result.types.push_back(FileTypes::Type::UNKNOWN); // cannot associate a single type to a collection
    }                                     
    if (style == FilterLayout::ONE_BY_ONE || style == FilterLayout::BOTH)
    {
      for (const auto& t : type_list_)
      {
        // e.g. 'FASTA file (*.fasta *.fa *.faa)'; the preferred extension comes first, which is the one
        // Qt appends when this filter is picked in a save dialog
        StringList patterns;
        for (const auto& ext : FileTypes::typeToExtensions(t)) patterns.push_back("*." + ext);
        result.items.push_back(FileTypes::typeToDescription(t) + " (" + ListUtils::concatenate(patterns, " ") + ")");
        result.types.push_back(t);
      }
    }
    if (add_all_filter)
    {
      result.items.emplace_back("all files (*)");
      result.types.push_back(FileTypes::Type::UNKNOWN); // cannot associate a single type to a collection
    }
    return result;
  }

  std::string FileTypes::typeToName(FileTypes::Type type)
  {
    for (const auto& t_info : type_with_annotation__)
    {
      if (t_info.type == type)
      {
        return t_info.name;
      }
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid type: Type has no name!",StringUtils::toStr(type));
  }

  std::string FileTypes::typeToDescription(Type type)
  {
    for (const auto& t_info : type_with_annotation__)
    {
      if (t_info.type == type) return t_info.description;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid type: Type has no description!",StringUtils::toStr(type));
  }


  std::vector<std::string> FileTypes::typeToExtensions(Type type)
  {
    for (const auto& t_info : type_with_annotation__)
    {
      if (t_info.type == type)
      {
        std::vector<std::string> result {t_info.name}; // preferred extension first
        result.insert(result.end(), t_info.aliases.begin(), t_info.aliases.end());
        return result;
      }
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid type: Type has no extensions!", StringUtils::toStr(type));
  }

  namespace
  {
    /// Case-insensitive equality without allocating a converted copy of either side.
    /// nameToType() runs this over the whole registry and sits in TOPPAS' per-file edge checks.
    bool equalsCaseInsensitive(const std::string& lhs, const std::string& rhs)
    {
      return lhs.size() == rhs.size()
          && std::equal(lhs.begin(), lhs.end(), rhs.begin(),
                        [](char a, char b) { return std::toupper(static_cast<unsigned char>(a)) == std::toupper(static_cast<unsigned char>(b)); });
    }
  }

  FileTypes::Type FileTypes::nameToType(const std::string& name)
  {
    // preferred extensions take precedence over aliases, so an alias can never shadow another type's canonical name
    for (const auto& t_info : type_with_annotation__)
    {
      if (equalsCaseInsensitive(t_info.name, name)) return t_info.type;
    }
    for (const auto& t_info : type_with_annotation__)
    {
      for (const auto& alias : t_info.aliases)
      {
        if (equalsCaseInsensitive(alias, name)) return t_info.type;
      }
    }

    return FileTypes::UNKNOWN;
  }


  bool FileTypes::supportsCompressedReading(Type type, Type compression)
  {
    if (compression != FileTypes::GZ && compression != FileTypes::BZ2 && compression != FileTypes::ZIP)
    {
      return false;
    }
    // Bruker's '.d.zip' is unpacked by BrukerTimsFile, not by XMLFile, and only for ZIP
    if (type == FileTypes::BRUKER_TDF)
    {
      return compression == FileTypes::ZIP;
    }
    for (const auto& t_info : type_with_annotation__)
    {
      if (t_info.type == type)
      {
        return std::find(t_info.features.begin(), t_info.features.end(), PROP::COMPRESSED_READABLE) != t_info.features.end();
      }
    }
    return false;
  }


  bool FileTypes::sameFormat(const std::string& lhs, const std::string& rhs)
  {
    const FileTypes::Type lhs_type = nameToType(lhs);
    const FileTypes::Type rhs_type = nameToType(rhs);
    if (lhs_type != FileTypes::UNKNOWN && rhs_type != FileTypes::UNKNOWN)
    {
      return lhs_type == rhs_type; // 'fasta' and 'fa' are the same format
    }
    // at least one side is a custom/unknown extension: two of those are only equal if they are spelled the same,
    // otherwise every unrecognized format would match every other one through UNKNOWN
    return equalsCaseInsensitive(lhs, rhs);
  }


  bool FileTypes::isDirectoryType(Type type)
  {
    return type == FileTypes::BRUKER_TDF
        || type == FileTypes::IDPARQUET
        || type == FileTypes::FEATUREPARQUET
        || type == FileTypes::CONSENSUSPARQUET;
  }


  std::string FileTypes::typeToMZML(FileTypes::Type type)
  {
    switch (type)
    {
      case FileTypes::DTA: return "DTA file";
      case FileTypes::DTA2D: return "DTA file"; // technically not correct, but closer than just a random CV term (currently mzData) - entry cannot be left empty
      case FileTypes::MZML: return "mzML file";
      case FileTypes::IMZML: return "mzML file"; // imzML is mzML 1.1 + IMS; reuse mzML source file term
      case FileTypes::MZDATA: return "PSI mzData file";
      case FileTypes::MZXML: return "ISB mzXML file";
      case FileTypes::MGF: return "Mascot MGF file";
      case FileTypes::XMASS: return "Bruker FID file";
      case FileTypes::BRUKER_TDF: return "Bruker TDF format";
      case FileTypes::RAW: return "Thermo RAW format";
      default: return "";
    }
  }
}
