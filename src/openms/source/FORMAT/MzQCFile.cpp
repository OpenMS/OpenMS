// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/QC/TIC.h>
#include <OpenMS/QC/SpectrumCount.h>
#include <OpenMS/QC/FeatureSummary.h>
#include <OpenMS/QC/IdentificationSummary.h>
#include <nlohmann/json.hpp>
#include <map>
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
                       const String& label,
                       const FeatureMap& feature_map,
                       vector<ProteinIdentification>& prot_ids,
                       vector<PeptideIdentification>& pep_ids) const
  {
    // --------------------------------------------------------------------
    // preparing output stream, quality metrics json object, CV, status
    // and initialize QC metric classes
    // --------------------------------------------------------------------
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

    QCBase::Status status;
    if (!input_file.empty())
    {
      status |= QCBase::Requires::RAWMZML;
    }
    if (!feature_map.empty())
    {
      status |= QCBase::Requires::PREFDRFEAT;
    }
    if (!prot_ids.empty() && !pep_ids.empty())
    {
      status |= QCBase::Requires::ID;
    }

    TIC tic;
    SpectrumCount spectrum_count;
    FeatureSummary feature_summary;
    IdentificationSummary identification_summary;

    // ---------------------------------------------------------------
    // function to add quality metrics to quality_metrics
    // ---------------------------------------------------------------
    auto addMetric = [&cv, &quality_metrics](const String& accession, const auto& value) -> void
      {
        json qm;
        qm["accession"] = accession;
        if (cv.exists(accession)) 
        {
          qm["name"] = cv.getTerm(accession).name;
        } 
        else
        {
          cout << accession << " not found in CV." << endl;
          return;
        }
        qm["value"] = value;
        quality_metrics.push_back(qm);
      };

    // ---------------------------------------------------------------
    // collecting quality metrics
    // ---------------------------------------------------------------

    if (spectrum_count.isRunnable(status))
    {
      auto counts = spectrum_count.compute(exp);
      // Number of MS1 spectra
      addMetric("QC:4000059", counts[1]);
      // Number of MS2 spectra
      addMetric("QC:4000060", counts[2]);
    }
    // Number of chromatograms"
    addMetric("QC:4000135", exp.getChromatograms().size());
    // Run time (RT duration)
    addMetric("QC:4000053", UInt(exp.getMaxRT() - exp.getMinRT()));
    // MZ acquisition range
    addMetric("QC:4000138", tuple<UInt,UInt>{exp.getMinMZ(), exp.getMaxMZ()});
    // TICs
    if (tic.isRunnable(status))
    {
      // complete TIC (all ms levels) with area
      auto result = tic.compute(exp, 0, 0);
      if (!result.intensities.empty())
      {
        json chrom;
        chrom["Relative intensity"] = result.relative_intensities;
        chrom["Retention time"] = result.retention_times;
        // Total ion current chromatogram
        addMetric("QC:4000067", chrom);
        // Area under TIC
        addMetric("QC:4000077", result.area);
      }
      // MS1
      result = tic.compute(exp,0,1);
      if (!result.intensities.empty())
      {
        json chrom;
        chrom["Relative intensity"] = result.relative_intensities;
        chrom["Retention time"] = result.retention_times;
        // MS1 Total ion current chromatogram
        addMetric("QC:4000069", chrom);
        // MS1 signal jump (10x) count
        addMetric("QC:4000172", result.jump);
        // MS1 signal fall (10x) count
        addMetric("QC:4000173", result.fall);
      }
      // MS2
      result = tic.compute(exp, 0, 2);
      if (!result.intensities.empty())
      {
        json chrom;
        chrom["Relative intensity"] = result.relative_intensities;
        chrom["Retention time"] = result.retention_times;
        // MS2 Total ion current chromatogram
        addMetric("QC:4000070", chrom);
      }
    }
    // Metabolomics: Detected compounds from featureXML file
    if (feature_summary.isRunnable(status))
    {
      auto result = feature_summary.compute(feature_map);
      // Detected compounds
      addMetric("QC:4000257", result.feature_count);
      // Retention time mean shift (sec)
      if (result.rt_shift_mean != 0)
      {
        addMetric("QC:4000262", result.rt_shift_mean);
      }   
    }
    // peptides and proteins from idXML file
    if (identification_summary.isRunnable(status))
    {
      auto result = identification_summary.compute(prot_ids, pep_ids);
      // Total number of PSM
      addMetric("QC:4000186", result.peptide_spectrum_matches);
      // Number of identified peptides at given FDR threshold
      addMetric("QC:4000187", make_tuple(result.unique_peptides.count, result.unique_peptides.fdr_threshold));
      // Identified peptide lengths - mean
      addMetric("QC:4000214", result.peptide_length_mean);
      // Missed cleavages - mean
      addMetric("QC:4000209", result.missed_cleavages_mean);
      // Number of identified proteins at given FDR threshold
      addMetric("QC:4000185", make_tuple(result.unique_proteins.count, result.unique_proteins.fdr_threshold));
      // Identification score mean (of protein hits)
      addMetric("QC:4000204", result.protein_hit_scores_mean);
    }

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
      out["mzQC"]["contactName"] = contact_name;
    }
    if (!contact_address.empty()) 
    {
      out["mzQC"]["contactAddress"] = contact_address;
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
          {"uri", "https://raw.githubusercontent.com/HUPO-PSI/mzQC/master/cv/qc-cv.obo"},
          {"version", "0.1.2"},
        },
        {
          {"name", "Proteomics Standards Initiative Mass Spectrometry Ontology"},
          {"uri", "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"},
          {"version", "4.1.56"}
        }
    };
    os << out.dump(2);
  }

}
