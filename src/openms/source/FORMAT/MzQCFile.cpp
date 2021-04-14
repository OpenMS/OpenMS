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

  void MzQCFile::store(const String& input_file,
                       const String& output_file,
                       const MSExperiment& exp,
                       const String& contact_name,
                       const String& contact_address,
                       const String& description,
                       const String& label) const
  {
    // ---------------------------------------------------------------
    // preparing output stream, quality metrics json object and CV
    // ---------------------------------------------------------------
    ofstream os(output_file.c_str());
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, output_file);
    }
    using json = nlohmann::ordered_json;
    json quality_metrics = {};
    ControlledVocabulary cv;
    cv.loadFromOBO("PSI-MS", File::find("/CV/psi-ms.obo"));
    cv.loadFromOBO("QC", File::find("/CV/qc-cv.obo"));

    // ---------------------------------------------------------------
    // function to add quality metrics to quality_metrics
    // ---------------------------------------------------------------
    auto addMetric = [&cv, &quality_metrics](const String& accession, const String& name, const auto& value) -> void
      {
        json qm;
        qm["accession"] = accession;
        if (cv.exists(accession)) 
        {
          qm["name"] = cv.getTerm(accession).name;
        } else
        {
          qm["name"] = name;
          cout << accession << " not found in CV." << endl;
        }
        qm["value"] = value;
        quality_metrics.push_back(qm);
      };

    // ---------------------------------------------------------------
    // collecting quality metrics
    // ---------------------------------------------------------------

    addMetric("MS:1000031", "instrument model", "LTQ Orbitrap Velos");



    // ---------------------------------------------------------------
    // writing mzQC file
    // ---------------------------------------------------------------
    json out;
    // required: creationDate, version
    DateTime currentTime = DateTime::now();
    out["mzQC"]["creationDate"] = currentTime.toString();
    out["mzQC"]["version"] = "1.0.0";

    // optional: contact_name, contact_address, description
    if (!contact_name.empty()) 
    {
      out["mzQC"]["contact_name"] = contact_name;
    }
    if (!contact_address.empty()) 
    {
      out["mzQC"]["contact_address"] = contact_address;
    }
    if (!description.empty()) 
    {
      out["mzQC"]["description"] = description;
    }
    // get OpenMS version for runQualities
    VersionInfo::VersionDetails version = VersionInfo::getVersionStruct();
    out["mzQC"]["runQualities"] =
    {
      {
        {"metadata",
          {
            {"label", label},
            {"inputFiles",
              {
                {
                  {"location", File::absolutePath(input_file)},
                  {"name", File::basename(input_file)},
                  {"fileFormat",
                    {
                      {"accession", "MS:10000584"},
                      {"name", "mzML format"}
                    }
                  },
                  {"fileProperties",
                    {
                      {
                        {"accession", "MS:1000747"},
                        {"name", "completion time"},
                        {"value", String(exp.getDateTime().getDate() + "T" + exp.getDateTime().getTime())}
                      },
                      {
                        {"accession", "MS:1000569"},
                        {"name", "SHA-1"},
                        {"value", String(FileHandler::computeFileHash(input_file))}
                      },
                      {
                        {"accession", "MS:1000031"},
                        {"name", "instrument model"},
                        {"value", String(exp.getInstrument().getName())}
                      }
                    }
                  }
                }
              }
            },
            {"analysisSoftware",
              {
                {
                  {"accession", "MS:1009001" }, // create new qc-cv for QCCalculator: MS:1009001 quality control metrics generating software
                  {"name", "QCCalculator"},
                  {"version", String(version.version_major)+"."+String(version.version_minor)+"."+String(version.version_patch)},
                  {"uri", "https://www.openms.de"}
                }
              }
            }
          }
        },
        {"qualityMetrics", quality_metrics}
      }
    };

    out["mzQC"]["controlledVocabularies"] = 
    {
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
    os << out.dump(2);
  }

}
