// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Stephan Aiche, Andreas Bertsch, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileTypes.h>

namespace OpenMS
{
  String FileTypes::typeToName(FileTypes::Type type)
  {
    std::map<Type, String>::const_iterator it = name_of_types_.find(type);
    if (it != name_of_types_.end())
    {
      return it->second;
    }
    else
    {
      return "";
    }
  }

  String FileTypes::typeToMZML(FileTypes::Type type)
  {
    std::map<Type, String>::const_iterator it = name_of_MZMLtypes_.find(type);
    if (it != name_of_MZMLtypes_.end())
    {
      return it->second;
    }
    else
    {
      return "";
    }
  }

  FileTypes::Type FileTypes::nameToType(const String& name)
  {
    String tmp = name;
    tmp.toUpper();
    String tmp2;

    for (int i = 0; i < FileTypes::SIZE_OF_TYPE; ++i)
    {
      tmp2 = FileTypes::typeToName((FileTypes::Type)i);
      tmp2.toUpper();
      if (tmp == tmp2)
      {
        return (FileTypes::Type)i;
      }
    }

    return FileTypes::UNKNOWN;
  }

  const std::map<FileTypes::Type, String> FileTypes::name_of_types_ = FileTypes::initializeMap_();
  const std::map<FileTypes::Type, String> FileTypes::name_of_MZMLtypes_ = FileTypes::initializeMZMLMap_();

  std::map<FileTypes::Type, String> FileTypes::initializeMap_()
  {
    std::map<Type, String> targetMap;
    targetMap[FileTypes::UNKNOWN] = "unknown";
    targetMap[FileTypes::DTA] = "dta";
    targetMap[FileTypes::DTA2D] = "dta2d";
    targetMap[FileTypes::MZDATA] = "mzData";
    targetMap[FileTypes::MZXML] = "mzXML";
    targetMap[FileTypes::FEATUREXML] = "featureXML";
    targetMap[FileTypes::IDXML] = "idXML";
    targetMap[FileTypes::CONSENSUSXML] = "consensusXML";
    targetMap[FileTypes::MGF] = "mgf";
    targetMap[FileTypes::INI] = "ini";
    targetMap[FileTypes::TOPPAS] = "toppas";
    targetMap[FileTypes::TRANSFORMATIONXML] = "trafoXML";
    targetMap[FileTypes::MZML] = "mzML";
    targetMap[FileTypes::CACHEDMZML] = "cachedMzML";
    targetMap[FileTypes::MS2] = "ms2";
    targetMap[FileTypes::PEPXML] = "pepXML";
    targetMap[FileTypes::PROTXML] = "protXML";
    targetMap[FileTypes::MZIDENTML] = "mzid";
    targetMap[FileTypes::MZQUANTML] = "mzq";
    targetMap[FileTypes::QCML] = "qcml";
    targetMap[FileTypes::GELML] = "gelML";
    targetMap[FileTypes::TRAML] = "traML";
    targetMap[FileTypes::MSP] = "msp";
    targetMap[FileTypes::OMSSAXML] = "omssaXML";
    targetMap[FileTypes::MASCOTXML] = "mascotXML";
    targetMap[FileTypes::PNG] = "png";
    targetMap[FileTypes::XMASS] = "fid";
    targetMap[FileTypes::TSV] = "tsv";
    targetMap[FileTypes::PEPLIST] = "peplist";
    targetMap[FileTypes::HARDKLOER] = "hardkloer";
    targetMap[FileTypes::KROENIK] = "kroenik";
    targetMap[FileTypes::FASTA] = "fasta";
    targetMap[FileTypes::EDTA] = "edta";
    targetMap[FileTypes::CSV] = "csv";
    targetMap[FileTypes::TXT] = "txt";
    targetMap[FileTypes::OBO] = "obo";
    targetMap[FileTypes::HTML] = "html";
    targetMap[FileTypes::XML] = "xml";
    targetMap[FileTypes::ANALYSISXML] = "analysisXML";
    targetMap[FileTypes::XSD] = "xsd";
    targetMap[FileTypes::PSQ] = "psq";
    targetMap[FileTypes::MRM] = "mrm";
    targetMap[FileTypes::SQMASS] = "sqMass";
    targetMap[FileTypes::PQP] = "pqp";
    targetMap[FileTypes::OSW] = "osw";
    targetMap[FileTypes::PSMS] = "psms";
    targetMap[FileTypes::PIN] = "pin";
    targetMap[FileTypes::SPLIB] = "splib";
    targetMap[FileTypes::NOVOR] = "novor";
    targetMap[FileTypes::PARAMXML] = "paramXML";
    
    return targetMap;
  }

  std::map<FileTypes::Type, String> FileTypes::initializeMZMLMap_()
  {
    std::map<Type, String> targetMap;
    targetMap[FileTypes::DTA] = "DTA file";
    targetMap[FileTypes::DTA2D] = "DTA file"; // technically not correct, but closer than just a random CV term (currently mzData) - entry cannot be left empty
    targetMap[FileTypes::MZML] = "mzML file";
    targetMap[FileTypes::MZDATA] = "PSI mzData file";
    targetMap[FileTypes::MZXML] = "ISB mzXML file";
    targetMap[FileTypes::MGF] = "Mascot MGF file";
    targetMap[FileTypes::XMASS] = "Bruker FID file";

    return targetMap;
  }

}
