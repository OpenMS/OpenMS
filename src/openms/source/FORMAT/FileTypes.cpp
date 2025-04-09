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
    String name;
    String description;
    std::vector<FileTypes::FileProperties> features;
    TypeNameBinding(FileTypes::Type ptype, String pname, String pdescription, std::vector<FileTypes::FileProperties> pfeatures)
      : type(ptype), name(std::move(pname)), description(std::move(pdescription)), features(pfeatures)
    {
      // Check that there are no double-spaces in the description, since Qt will replace "  " with " " in filters supplied to QFileDialog::getSaveFileName.
      // And if you later ask for the selected filter, you will get a different string back.
      assert(description.find("  ") == std::string::npos);
    }
  };

  using PROP = FileTypes::FileProperties;   // shorten our syntax a bit
  /// Maps the FileType::Type to the preferred extension.
  /// when adding new types, be sure to update the FileTypes_test typesWithProperties test to match the new files
  static const std::array<TypeNameBinding, FileTypes::SIZE_OF_TYPE> type_with_annotation__ =
  {
    TypeNameBinding(FileTypes::UNKNOWN, "unknown", "unknown file extension", {}),
    TypeNameBinding(FileTypes::DTA, "dta", "dta raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::PROVIDES_SPECTRUM, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::DTA2D, "dta2d", "dta2d raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MZDATA, "mzData", "mzData raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MZXML, "mzXML", "mzXML raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::FEATUREXML, "featureXML", "OpenMS feature map", {PROP::PROVIDES_FEATURES, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::IDXML, "idXML", "OpenMS peptide identification file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::CONSENSUSXML, "consensusXML", "OpenMS consensus feature map", {PROP::PROVIDES_CONSENSUSFEATURES, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MGF, "mgf", "mascot generic format file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::INI, "ini", "OpenMS parameter file", {PROP::READABLE}),
    TypeNameBinding(FileTypes::TOPPAS, "toppas", "OpenMS TOPPAS pipeline", {PROP::READABLE}),
    TypeNameBinding(FileTypes::TRANSFORMATIONXML, "trafoXML", "RT transformation file", {PROP::PROVIDES_TRANSFORMATIONS, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MZML, "mzML", "mzML raw data file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE}),
    //TODO: Add support for cachedMZML as a first class file type
    TypeNameBinding(FileTypes::CACHEDMZML, "cachedMzML", "cachedMzML raw data file", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MS2, "ms2", "ms2 file", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE}),
    TypeNameBinding(FileTypes::PEPXML, "pepXML", "pepXML file", {PROP::READABLE, PROP::WRITEABLE}), //Supported for loading and storing identifications but TODO integrate this into fileHandler
    TypeNameBinding(FileTypes::PROTXML, "protXML", "protXML file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE}),
    TypeNameBinding(FileTypes::MZIDENTML, "mzid", "mzIdentML file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::QCML, "qcml", "quality control file", {PROP::PROVIDES_QC, PROP::WRITEABLE}), //TODO add load functions for QC
    TypeNameBinding(FileTypes::MZQC, "mzqc", "quality control file in json format", {PROP::PROVIDES_QC, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::GELML, "gelML", "gelML file", {}),
    TypeNameBinding(FileTypes::TRAML, "traML", "transition file", {PROP::PROVIDES_TRANSITIONS, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MSP, "msp", "NIST spectra library file format", {PROP::PROVIDES_EXPERIMENT, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::OMSSAXML, "omssaXML", "omssaXML file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE}),
    TypeNameBinding(FileTypes::MASCOTXML, "mascotXML", "mascotXML file", {}),
    TypeNameBinding(FileTypes::PNG, "png", "portable network graphics file", {}),
    TypeNameBinding(FileTypes::XMASS, "fid", "XMass analysis file", {PROP::PROVIDES_EXPERIMENT, PROP::PROVIDES_SPECTRUM, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::TSV, "tsv", "tab-separated file", {PROP::PROVIDES_FEATURES, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::MZTAB, "mzTab", "mzTab file", {}), //TODO add filehandler support for MZTAB
    TypeNameBinding(FileTypes::PEPLIST, "peplist", "SpecArray file", {PROP::PROVIDES_FEATURES, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::HARDKLOER, "hardkloer", "hardkloer file", {}),
    TypeNameBinding(FileTypes::KROENIK, "kroenik", "kroenik file", {PROP::PROVIDES_FEATURES, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::FASTA, "fasta", "FASTA file", {PROP::READABLE, PROP::WRITEABLE}),
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
    TypeNameBinding(FileTypes::MS, "ms", "SIRIUS file", {}),
    TypeNameBinding(FileTypes::OSW, "osw", "OpenSwath output files", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::PSMS, "psms", "Percolator tab-delimited output (PSM level)", {PROP::READABLE}),
    TypeNameBinding(FileTypes::PIN, "pin", "Percolator tab-delimited input (PSM level)", {}),
    TypeNameBinding(FileTypes::PARAMXML, "paramXML", "OpenMS internal XML file", {}),
    TypeNameBinding(FileTypes::SPLIB, "splib", "SpectraST binary spectral library file", {}),
    TypeNameBinding(FileTypes::NOVOR, "novor", "Novor custom parameter file", {}),
    TypeNameBinding(FileTypes::XQUESTXML, "xquest.xml", "xquest.xml file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::SPECXML, "spec.xml", "spec.xml file", {}),
    TypeNameBinding(FileTypes::JSON, "json", "JavaScript Object Notation file", {PROP::READABLE, PROP::WRITEABLE}),
    TypeNameBinding(FileTypes::RAW, "raw", "(Thermo) Raw data file", {}),
    TypeNameBinding(FileTypes::OMS, "oms", "OpenMS SQLite file", {PROP::PROVIDES_IDENTIFICATIONS, PROP::PROVIDES_FEATURES, PROP::PROVIDES_CONSENSUSFEATURES}),
    TypeNameBinding(FileTypes::EXE, "exe", "Windows executable", {}),
    TypeNameBinding(FileTypes::BZ2, "bz2", "bzip2 compressed file", {PROP::READABLE}),
    TypeNameBinding(FileTypes::GZ, "gz", "gzip compressed file", {PROP::READABLE}),
    TypeNameBinding(FileTypes::XML, "xml", "any XML file", {PROP::READABLE}),  // make sure this comes last, since the name is a suffix of other formats and should only be matched last
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

  String FileTypeList::toFileDialogFilter(const FilterLayout style, bool add_all_filter) const
  {
    return ListUtils::concatenate(asFilterElements_(style, add_all_filter).items, ";;");
  }

  FileTypes::Type FileTypeList::fromFileDialogFilter(const String& filter, const FileTypes::Type fallback) const
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
    
    for (auto t : good_types)
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
        items.push_back("*." + FileTypes::typeToName(t));
      }
      result.items.emplace_back("all readable files (" + ListUtils::concatenate(items, " ") + ")");
      result.types.push_back(FileTypes::Type::UNKNOWN); // cannot associate a single type to a collection
    }                                     
    if (style == FilterLayout::ONE_BY_ONE || style == FilterLayout::BOTH)
    {
      StringList items;
      for (const auto& t : type_list_)
      {
        result.items.push_back(FileTypes::typeToDescription(t) + " (*." + FileTypes::typeToName(t) + ")");
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

  String FileTypes::typeToName(FileTypes::Type type)
  {
    for (const auto& t_info : type_with_annotation__)
    {
      if (t_info.type == type)
      {
        return t_info.name;
      }
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid type: Type has no name!", String(type));
  }

  String FileTypes::typeToDescription(Type type)
  {
    for (const auto& t_info : type_with_annotation__)
    {
      if (t_info.type == type) return t_info.description;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid type: Type has no description!", String(type));
  }


  FileTypes::Type FileTypes::nameToType(const String& name)
  {
    String name_upper = String(name).toUpper();

    for (const auto& t_info : type_with_annotation__)
    {
      if (String(t_info.name).toUpper() == name_upper)
      {
        return t_info.type;
      }
    }

    return FileTypes::UNKNOWN;
  }


  String FileTypes::typeToMZML(FileTypes::Type type)
  {
    switch (type)
    {
      case FileTypes::DTA: return "DTA file";
      case FileTypes::DTA2D: return "DTA file"; // technically not correct, but closer than just a random CV term (currently mzData) - entry cannot be left empty
      case FileTypes::MZML: return "mzML file";
      case FileTypes::MZDATA: return "PSI mzData file";
      case FileTypes::MZXML: return "ISB mzXML file";
      case FileTypes::MGF: return "Mascot MGF file";
      case FileTypes::XMASS: return "Bruker FID file";
      default: return "";
    }
  }
}
