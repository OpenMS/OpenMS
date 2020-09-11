// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Andreas Bertsch, Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileTypes.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <array>

namespace OpenMS
{

  static FileTypes::FileTypeList<1> test_type_list({ FileTypes::MZML });
  /// connect the type to some other information
  /// We could also use paired arrays, but this way, its less likely to have mismatches if a new type is added
  struct TypeNameBinding
  {
    FileTypes::Type type;
    String name;
    String description;
    TypeNameBinding(FileTypes::Type type, String name, String description)
      : type(type), name(name), description(description)
    {
    }
  };
  
  /// Maps the FileType::Type to the preferred extension.
  static const std::array<TypeNameBinding, FileTypes::SIZE_OF_TYPE> type_with_annotation__ =
  {
    TypeNameBinding(FileTypes::UNKNOWN, "unknown", "unknown file extension"),
    TypeNameBinding(FileTypes::DTA, "dta", "dta raw data file"),
    TypeNameBinding(FileTypes::DTA2D, "dta2d", "dta2d raw data file"),
    TypeNameBinding(FileTypes::MZDATA, "mzData", "mzData raw data file"),
    TypeNameBinding(FileTypes::MZXML, "mzXML", "mzXML raw data file"),
    TypeNameBinding(FileTypes::FEATUREXML, "featureXML", "OpenMS feature map"),
    TypeNameBinding(FileTypes::IDXML, "idXML", "OpenMS peptide identification  file"),
    TypeNameBinding(FileTypes::CONSENSUSXML, "consensusXML", "OpenMS consensus feature map"),
    TypeNameBinding(FileTypes::MGF, "mgf", "mascot generic format file"),
    TypeNameBinding(FileTypes::INI, "ini", "OpenMS parameter file"),
    TypeNameBinding(FileTypes::TOPPAS, "toppas", "OpenMS TOPPAS pipeline"),
    TypeNameBinding(FileTypes::TRANSFORMATIONXML, "trafoXML", "RT transformation file"),
    TypeNameBinding(FileTypes::MZML, "mzML", "mzML raw data file"),
    TypeNameBinding(FileTypes::CACHEDMZML, "cachedMzML", "cachedMzML raw data file"),
    TypeNameBinding(FileTypes::MS2, "ms2", "ms2 file"),
    TypeNameBinding(FileTypes::PEPXML, "pepXML", "pepXML file"),
    TypeNameBinding(FileTypes::PROTXML, "protXML", "protXML file"),
    TypeNameBinding(FileTypes::MZIDENTML, "mzid", "mzIdentML file"),
    TypeNameBinding(FileTypes::MZQUANTML, "mzq", "mzQuantML file"),
    TypeNameBinding(FileTypes::QCML, "qcml", "quality control file"),
    TypeNameBinding(FileTypes::GELML, "gelML", "gelML file"),
    TypeNameBinding(FileTypes::TRAML, "traML", "transition file"),
    TypeNameBinding(FileTypes::MSP, "msp", "NIST spectra library file format"),
    TypeNameBinding(FileTypes::OMSSAXML, "omssaXML", "omssaXML file"),
    TypeNameBinding(FileTypes::MASCOTXML, "mascotXML", "mascotXML file"),
    TypeNameBinding(FileTypes::PNG, "png", "portable network graphics file"),
    TypeNameBinding(FileTypes::XMASS, "fid", "XMass analysis file"),
    TypeNameBinding(FileTypes::TSV, "tsv", "tab-separated file"),
    TypeNameBinding(FileTypes::MZTAB, "mzTab", "mzTab file"),
    TypeNameBinding(FileTypes::PEPLIST, "peplist", "SpecArray file"),
    TypeNameBinding(FileTypes::HARDKLOER, "hardkloer", "hardkloer file"),
    TypeNameBinding(FileTypes::KROENIK, "kroenik", "kroenik file"),
    TypeNameBinding(FileTypes::FASTA, "fasta", "FASTA file"),
    TypeNameBinding(FileTypes::EDTA, "edta", "enhanced dta file"),
    TypeNameBinding(FileTypes::CSV, "csv", "comma-separated values file"),
    TypeNameBinding(FileTypes::TXT, "txt", "generic text file"),
    TypeNameBinding(FileTypes::OBO, "obo", "controlled vocabulary file"),
    TypeNameBinding(FileTypes::HTML, "html", "any HTML file"),
    TypeNameBinding(FileTypes::ANALYSISXML, "analysisXML", "analysisXML file"),
    TypeNameBinding(FileTypes::XSD, "xsd", "XSD schema format"),
    TypeNameBinding(FileTypes::PSQ, "psq", "NCBI binary blast db"),
    TypeNameBinding(FileTypes::MRM, "mrm", "SpectraST MRM list"),
    TypeNameBinding(FileTypes::SQMASS, "sqMass", "SqLite format for mass and chromatograms"),
    TypeNameBinding(FileTypes::PQP, "pqp", "pqp file"),
    TypeNameBinding(FileTypes::MS, "ms", "SIRIUS file"),
    TypeNameBinding(FileTypes::OSW, "osw", "OpenSwath output files"),
    TypeNameBinding(FileTypes::PSMS, "psms", "Percolator tab-delimited output (PSM level)"),
    TypeNameBinding(FileTypes::PIN, "pin", "Percolator tab-delimited input (PSM level)"),
    TypeNameBinding(FileTypes::PARAMXML, "paramXML", "OpenMS internal XML file"),
    TypeNameBinding(FileTypes::SPLIB, "splib", "SpectraST binary spectral library file"),
    TypeNameBinding(FileTypes::NOVOR, "novor", "Novor custom parameter file"),
    TypeNameBinding(FileTypes::XQUESTXML, "xquest.xml", "xquest.xml file"),
    TypeNameBinding(FileTypes::SPECXML, "spec.xml", "spec.xml file"),
    TypeNameBinding(FileTypes::JSON, "json", "JavaScript Object Notation file"),
    TypeNameBinding(FileTypes::RAW, "raw", "(Thermo) Raw data file"),
    TypeNameBinding(FileTypes::EXE, "exe", "Windows executable"),
    TypeNameBinding(FileTypes::BZ2, "bz2", "bzip2 compressed file"),
    TypeNameBinding(FileTypes::GZ, "gz", "gzip compressed file"),
    TypeNameBinding(FileTypes::XML, "xml", "any XML file")  // make sure this comes last, since the name is a suffix of other formats and should only be matched last
  };

  
  String FileTypes::typeToName(FileTypes::Type type)
  {
    for (const auto& t_info : type_with_annotation__)
    {
      if (t_info.type == type) return t_info.name;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Type has no name!", String(type));
  }

  String FileTypes::typeToDescription(Type type)
  {
    for (const auto& t_info : type_with_annotation__)
    {
      if (t_info.type == type) return t_info.description;
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Type has no description!", String(type));
  }


  FileTypes::Type FileTypes::nameToType(const String& name)
  {
    String name_upper = String(name).toUpper();
    
    for (const auto& t_info : type_with_annotation__)
    {
      if (String(t_info.name).toUpper() == name_upper) return t_info.type;
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
