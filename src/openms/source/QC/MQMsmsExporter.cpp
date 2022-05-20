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
// $Maintainer: Chris Bielow$
// $Authors: Hendrik Beschorner, Lenny Kovac, Virginia Rossow$
// --------------------------------------------------------------------------

#include <OpenMS/QC/MQMsmsExporter.h>
#include <OpenMS/QC/MQExporterHelper.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <QtCore/QDir>
#include <cmath> // isnan
#include <fstream>


using namespace OpenMS;


MQMsms::MQMsms(const String& path)
{
  if (path.empty())
  {
    return;
  }
  filename_ = path + "/msms.txt";
  try
  {
    QString evi_path = QString::fromStdString(path);
    QDir().mkpath(evi_path);
    file_ = std::fstream(filename_, std::fstream::out);
  }
  catch (...)
  {
    OPENMS_LOG_FATAL_ERROR << filename_ << " wasn’t created" << std::endl;
    throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "out_evd");
  }
  exportHeader_();
}


MQMsms::~MQMsms()
{
  file_.close();
}


void MQMsms::exportHeader_()
{

  file_ << "Raw file" << "\t"; // also in Evidence
  file_ << "Scan number" << "\t";
  file_ << "Scan index" << "\t";
  file_ << "Sequence" << "\t"; // also in Evidence
  file_ << "Length" << "\t"; // also in Evidence
  file_ << "Missed cleavages" << "\t"; // also in Evidence
  file_ << "Modifications" << "\t"; // also in Evidence
  file_ << "Modified sequence" << "\t"; // also in Evidence
  //file_ << "Oxidation (M) Probabilities" << "\t"; --> not supported by OpenMS
  //file_ << "Oxidation (M) Score diffs" << "\t"; --> not supported by OpenMS
  file_ << "Acetyl (Protein N-term)" << "\t"; // also in Evidence
  file_ << "Oxidation (M)" << "\t"; // also in Evidence
  file_ << "Proteins" << "\t"; // also in Evidence
  file_ << "Charge" << "\t"; // also in Evidence
  file_ << "Fragmentation" << "\t";
  file_ << "Mass analyzer" << "\t";
  file_ << "Type" << "\t"; // also in Evidence
  file_ << "Scan event number" << "\t";
  file_ << "Isotope index" << "\t";
  file_ << "m/z" << "\t"; // also in Evidence
  file_ << "Mass" << "\t"; // also in Evidence
  file_ << "Mass error [ppm]" << "\t"; // also in Evidence
  file_ << "Mass error [Da]" << "\t"; // also in Evidence
  file_ << "Simple mass error [ppm]" << "\t";
  file_ << "Retention time" << "\t"; // also in Evidence
  file_ << "PEP" << "\t"; // also in Evidence
  file_ << "Score" << "\t"; // also in Evidence
  file_ << "Delta score" << "\t"; // also in Evidence
  file_ << "Score diff" << "\t";
  file_ << "Localization prob" << "\t";
  //file_ << "Combinatorics" << "\t"; --> not supported by OpenMS
  // file_ << "PIF" << "\t"; --> not practical to implement
  file_ << "Fraction of total spectrum" << "\t";
  file_ << "Base peak fraction" << "\t"; // also in Evidence
  file_ << "Precursor full scan number" << "\t";
  file_ << "Precursor Intensity" << "\t";
  file_ << "Precursor apex fraction" << "\t";
  file_ << "Precursor apex offset" << "\t";
  file_ << "Precursor apex offset time" << "\t";
  file_ << "Matches Intensities" << "\t";
  file_ << "Mass deviations [Da]" << "\t";
  file_ << "Mass deviations [ppm]" << "\t";
  file_ << "Masses" << "\t";
  file_ << "Number of matches" << "\t";
  file_ << "Intensity coverage" << "\t"; // vielleicht aus intensity in evidence.txt berechenbar?
  file_ << "Peak coverage" << "\t";
  file_ << "Neutral loss level" << "\t";
  file_ << "ETD identification type" << "\t";
  file_ << "Reverse" << "\t"; // also in Evidence
  file_ << "All scores" << "\t";
  file_ << "All sequences" << "\t";
  file_ << "All modified sequences" << "\t";
  //file_ << "Reporter PIF" << "\t"; --> not supported by OpenMS
  //file_ << "Reporter fraction" << "\t"; --> not supported by OpenMS
  file_ << "id" << "\t"; // also in Evidence
  file_ << "Protein group IDs" << "\n"; // also in Evidence
  //file_ << "Peptide ID" << "\t"; --> not useful without the other MQ files
  //file_ << "Mod. peptide ID" << "\t"; --> not useful without the other MQ files
  //file_ << "Evidence ID" << "\t"; // vermutlich in evidence.txt
  //file_ << "Oxidation (M) site IDs" << "\n"; --> not useful without the other MQ files
  
}

void MQMsms::exportRowFromFeature_(
        const Feature& f,
        const ConsensusMap& cmap,
        const Size c_feature_number,
        const String& raw_file,
        const std::multimap<String, std::pair<Size, Size>>& UIDs,
        const ProteinIdentification::Mapping& mp_f,
        const OpenMS::MSExperiment& exp = {})
{

  MQExporterHelper::MQCommonOutputs common_outputs{f, cmap, c_feature_number, raw_file, UIDs, mp_f, exp, prot_mapper};

  const PeptideHit* ptr_best_hit; // the best hit referring to score
  const ConsensusFeature& cf = cmap[c_feature_number];
  Size pep_ids_size = 0;
  String type;
  if (MQExporterHelper::hasValidPepID_(f, c_feature_number, UIDs, mp_f))
  {
    for (Size i = 1; i < f.getPeptideIdentifications().size(); ++i) // for msms-count
    {
      if (!f.getPeptideIdentifications()[i].getHits().empty())
      {
        if (f.getPeptideIdentifications()[i].getHits()[0].getSequence() == f.getPeptideIdentifications()[0].getHits()[0].getSequence())
        {
          ++pep_ids_size;
        }
        else
          break;
      }
    }
    type = "MULTI-MSMS";
    ptr_best_hit = &f.getPeptideIdentifications()[0].getHits()[0];
  }
    else if (MQExporterHelper::hasPeptideIdentifications_(cf))
  {
    type = "MULTI-MATCH";
    ptr_best_hit = &cf.getPeptideIdentifications()[0].getHits()[0];
  }
  else
  {
    return; // no valid PepID; nothing to export
  }
  
  const double& max_score = ptr_best_hit->getScore();
  const AASequence& pep_seq = ptr_best_hit->getSequence();

  if (pep_seq.empty())
  {
    return;
  }

// what is written in the file in this exact order 

  file_ << raw_file << "\t"; // raw file
  file_ << "Scan number" << "\t";
  file_ << "Scan index" << "\t";
  file_ << pep_seq.toUnmodifiedString() << "\t"; // Sequence
  file_ << pep_seq.size() << "\t";               // Length
  file_ << ptr_best_hit->getMetaValue("missed_cleavages", "NA") << "\t"; // missed cleavages


  file_ << common_outputs.modifications << "\t"; // Modifications
  file_ << "_" << common_outputs.modified_sequence << "_" << "\t"; // Modified Sequence
  file_ << common_outputs.acetyl << "\t"; // Acetyl (Protein N-term)
  file_ << common_outputs.oxidation << "\t"; // Oxidation (M)
  
  const std::set<String>& accessions = ptr_best_hit->extractProteinAccessionsSet();
  file_ << ListUtils::concatenate(accessions, ";") << "\t";  // Proteins
  
  file_ << f.getCharge() << "\t"; // Charge
  
  file_ << "Fragmentation" << "\t"; // Fragmentation
  file_ << "Mass analyzer" << "\t"; // Mass analyzer
  file_ << type << "\t"; // type
  file_ << "Scan event number" << "\t"; // Scan event number
  file_ << "Isotope index" << "\t"; // Isotope index
  file_ << f.getMZ() << "\t"; // M/Z
  file_ << pep_seq.getMonoWeight() << "\t"; // Mass
  file_ << common_outputs.mass_error_ppm << "\t"; // Mass Error [ppm]
  file_ << common_outputs.mass_error_da << "\t"; // Mass error [Da]
  file_ << "Simple mass error [ppm]" << "\t"; // Simple mass error [ppm]

  f.metaValueExists("rt_raw_end") && f.metaValueExists("rt_raw_start") ?
    file_ << (double(f.getMetaValue("rt_raw_end")) - double(f.getMetaValue("rt_raw_start"))) / 60 << "\t" : file_
      << "NA" << "\t"; // Retention time
  // hier weiter...
  file_ << "PEP" << "\t"; // PEP
  file_ << "Score" << "\t"; // Score
  file_ << "Delta score" << "\t"; // Delta score
  file_ << "Score diff" << "\t"; // Score diff
  file_ << "Localization prob" << "\t"; // Localization prob
  file_ << "Fraction of total spectrum" << "\t"; // Fraction of total spectrum
  file_ << "Base peak fraction" << "\t"; // Base peak fraction
  file_ << "Precursor full scan number" << "\t"; // Precursor full scan number
  file_ << "Precursor Intensity" << "\t"; // Precursor Intensity
  file_ << "Precursor apex fraction" << "\t"; // Precursor apex fraction
  file_ << "Precursor apex offset" << "\t"; // Precursor apex offset
  file_ << "Precursor apex offset time" << "\t"; // Precursor apex offset time
  file_ << "Matches Intensities" << "\t"; // Matches Intensities
  file_ << "Mass deviations [Da]" << "\t"; // Mass deviations [Da]
  file_ << "Mass deviations [ppm]" << "\t"; // Mass deviations [ppm]
  file_ << "Masses" << "\t"; // Masses
  file_ << "Number of matches" << "\t"; // Number of matches
  file_ << "Intensity coverage" << "\t"; // Intensity coverage
  file_ << "Peak coverage" << "\t"; // Peak coverage
  file_ << "Neutral loss level" << "\t"; // Neutral loss level
  file_ << "ETD identification type" << "\t"; // ETD identification type
  ptr_best_hit->getMetaValue("target_decoy") == "decoy" ? file_ << "1"
                                                                << "\t" :
                                                          file_ << "\t"; // reverse
                                                          
  file_ << "All scores" << "\t"; // All scores
  file_ << "All sequences" << "\t"; // All sequences
  file_ << "All modified sequences" << "\t"; // All modified sequences
  
  file_ << id_ << "\t"; // ID
  ++id_;
  
  file_ << MQExporterHelper::proteinGroupID_(acessions(0)); // nicht sicher
  for (const String& p : accessions)
  {
    file_ << ";" << MQExporterHelper::proteinGroupID_(p+1); // Protein group ids
  }
  file_ << "\t";

  file_ << "Evidence ID" << "\t"; // Evidence ID


}










