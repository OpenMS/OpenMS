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
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzQCFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

#include <nlohmann/json.hpp>

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

using namespace std;

namespace OpenMS
{ 

  void MzQCFile::store(const String & inputFileName,
                       const String & outputFileName,
                       const MSExperiment & exp,
                       const String & contactName,
                       const String & contactAddress,
                       const String & description,
                       const String & label) const
  {
    using json = nlohmann::ordered_json;
    json out;
    // required: creationDate, version
    DateTime currentTime = DateTime::now();
    out["mzQC"]["creationDate"] = currentTime.toString();
    out["mzQC"]["version"] = "1.0.0";

    // optional: contactName, contactAddress, description
    if (!contactName.empty()) 
    {
      out["mzQC"]["contactName"] = contactName;
    }
    if (!contactAddress.empty()) 
    {
      out["mzQC"]["contactAddress"] = contactAddress;
    }
    if (!description.empty()) 
    {
      out["mzQC"]["description"] = description;
    }

    out["mzQC"]["runQualities"]["metadata"]["label"] = label;
    out["mzQC"]["runQualities"]["metadata"]["inputFiles"] = {
        {
          {"location", File::absolutePath(inputFileName)},
          {"name", File::basename(inputFileName)},
          {"fileFormat",
          {
            {"accession", "MS:10000584"},
            {"name", "mzML format"}
          }
          }
        },
          "fileProperties",
          {
            {
              {"accession", "MS:1000747"},
              {"name", "completion time"},
              {"value", String(exp.getDateTime().getDate() + "T" + exp.getDateTime().getTime())}
            },
            {
              {"accession", "MS:1000569"},
              {"name", "SHA-1"},
              {"value", String(FileHandler::computeFileHash(inputFileName))}
            },
            {
              {"accession", "MS:1000031"},
              {"name", "instrument model"},
              {"value", String(exp.getInstrument().getName())}
            }
          }
    };

    VersionInfo::VersionDetails version = VersionInfo::getVersionStruct();
    out["mzQC"]["runQualities"]["metadata"]["analysisSoftware"] = { // todo
        {
          {"accession", "MS:1009001" }, // create new qc-cv for QCCalculator: MS:1009001 quality control metrics generating software
          {"name", "QCCalculator"},
          {"version", String(version.version_major)+"."+String(version.version_minor)+"."+String(version.version_patch)},
          {"uri", "https://www.openms.de"}
        }
    };

    out["mzQC"]["runQualities"]["qualityMetrics"] = {};
    
    out["mzQC"]["controlledVocabularies"] = {
        {
          {"name", "Proteomics Standards Initiative Quality Control Ontology"},
          {"uri", "https://github.com/HUPO-PSI/mzQC/blob/master/cv/qc-cv.obo"},
          {"version", "0.1.2"},
        },
        {
          {"name", "Proteomics Standards Initiative Mass Spectrometry Ontology"},
          {"uri", "https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo"},
          {"version", "4.1.49"}
        }
    };
    ControlledVocabulary cv;
    cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("QC", File::find("/CV/qc-cv.obo"));

    // create qualityMetric in json format
    auto addValue = [&cv](const String& accession, const String& name, const auto& value) -> json
      {
        json qm;
        qm["accession"] = accession;
        if (cv.exists(accession)) 
        {
          qm["name"] = cv.getTerm(accession).name;
        } 
        else
        {
          qm["name"] = name;
          cout << accession << " not found in CV." << endl;
        }
        qm["value"] = value;
        return qm;
      };

    map<Size, UInt> counts;
    for (const auto& spectrum : exp)
      {
        const Size level = spectrum.getMSLevel();
        ++counts[level];  // count MS level
      }
    // QC:4000059 Number of MS1 spectra
    out["mzQC"]["runQualities"]["qualityMetrics"] += addValue("QC:4000059", "Number of MS1 spectra", counts[1]);
    // QC:4000060 Number of MS2 spectra
    out["mzQC"]["runQualities"]["qualityMetrics"] += addValue("QC:4000060", "Number of MS2 spectra", counts[2]);
    // QC:4000135 Number of chromatograms"
    out["mzQC"]["runQualities"]["qualityMetrics"] += addValue("QC:4000135", "Number of chromatograms", exp.getChromatograms().size());
    // QC:4000053 Run time (RT duration)
    out["mzQC"]["runQualities"]["qualityMetrics"] += addValue("QC:4000053", "RT duration", UInt(exp.getMaxRT() - exp.getMinRT()));
    // QC:4000138 MZ acquisition range
    out["mzQC"]["runQualities"]["qualityMetrics"] += addValue("QC:4000138", "MZ acquisition range", tuple<UInt,UInt>{exp.getMinMZ(), exp.getMaxMZ()});
    
    const MSChromatogram& tic = exp.getTIC();
    UInt jump = 0;
    UInt fall = 0;
    vector<float> retention_times;
    vector<UInt> intensities;
    
    for (const auto& p : tic)
    {
      intensities.push_back(p.getIntensity());
      retention_times.push_back(p.getRT());
    }

    UInt tic_area = intensities[0];

    for (UInt i = 1; i < intensities.size(); ++i)
    {
      tic_area += intensities[i];
      if (intensities[i] > intensities[i-1] * 10) // detect 10x jumps between two subsequent scans
      {
        ++jump;
      }
      if (intensities[i] < intensities[i-1] / 10) // detect 10x falls between two subsequent scans
      {
        ++fall;
      }
    }
    
    json tic_values;
    tic_values["Relative intensity"] = intensities;
    tic_values["Retention time"] = retention_times;

    if (!tic.empty())
    {
      // QC:4000067 Total ion current chromatogram
      out["mzQC"]["runQualities"]["qualityMetrics"] += addValue("QC:4000067", "Total ion current chromatogram", tic_values);
      // QC:4000077 Area under TIC
      out["mzQC"]["runQualities"]["qualityMetrics"] += addValue("QC:4000077", "Area under TIC", tic_area);
      // QC:4000172 MS1 signal jump (10x) count
      out["mzQC"]["runQualities"]["qualityMetrics"] += addValue("QC:4000172", "MS1 signal jump (10x) count", jump);
      // QC:4000173 MS1 signal fall (10x) count
      out["mzQC"]["runQualities"]["qualityMetrics"] += addValue("QC:4000173", "MS1 signal fall (10x) count", fall);
    }

    // QC:4000074 Median MS1 peak FWHM for peptides
    // QC:4000257 Detected Compounds
    // MS:1000005 sample volume
    // MS:1000011 mass resolution
    // MS:1000015 scan rate
    // MS:1000026 detector type

    // HPLC:
    // QC:4000107 Column type
    // QC:4000079 Pump pressure chromatogram

    //open stream
    ofstream os(outputFileName.c_str());
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, outputFileName);
    }

    // write out the json object in proper format with indentation level of 2
    os << out.dump(2);
  }

}
