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
// $Authors: Valentin Noske, Vincent Musch$
// --------------------------------------------------------------------------

#include <OpenMS/QC/MQEvidenceExporter.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
//#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>

#include <QtCore/QDir>
#include <cmath> // isnan
#include <fstream>


using namespace OpenMS;

MQEvidence::MQEvidence(const String& path)
{
  if (path.empty())
  {
    return;
  }
  filename_ = path + "/evidence.txt";
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

MQEvidence::~MQEvidence()
{
  file_.close();
}


bool MQEvidence::isValid()
{
  return File::writable(filename_);
}

void MQEvidence::exportHeader_()
{

  file_ << "Sequence" << "\t"; // ok
  file_ << "Length" << "\t"; // ok
  file_ << "Modifications" << "\t"; // ok
  file_ << "Modified sequence" << "\t"; // ok
  file_ << "Oxidation (M) Probabilities" << "\t"; // wahrscheinlich nicht in OpenM, geschaut in ModificationsDB, ResidueModification, AASequence, Residue, Unimod
  file_ << "Oxidation (M) Score Diffs" << "\t"; // siehe Probabilities
  file_ << "Acetyl (Protein N-term)" << "\t"; // ok (in beiden keine)
  file_ << "Oxidation (M)" << "\t"; // ok
  file_ << "Missed cleavages" << "\t"; // könnte richtig sein (aber schon anders)
  file_ << "Proteins" << "\t"; // ok
  file_ << "Leading Proteins" << "\t"; // - erledigt
  file_ << "Leading Razor Protein" << "\t"; // - erledigt
  file_ << "Gene Names" << "\t"; // noch unklar, ob in OpenMS, Light Compound, Peptide
  file_ << "Protein Names" << "\t"; // noch unklar, ob in OpenMS
  file_ << "Type" << "\t"; // ok
  file_ << "Raw file" << "\t"; // (wahrscheinlich) ok
  file_ << "Fraction" << "\t"; // noch unklar, ob in OpenMS
  file_ << "MS/MS m/z" << "\t"; // noch unklar, ob in OpenMS
  file_ << "Charge" << "\t"; // ok (unterschied recht gross aber andere Daten)
  file_ << "m/z" << "\t"; // ok
  file_ << "Mass" << "\t"; // wahrscheinlich falsch (ca. 1000 (OpenMS) vs ca. 100000000 (MQ))
  file_ << "Resolution" << "\t"; // wahrscheinlich falsch (ca. 0.16 (OpenMS) vs ca. 50000 (MQ))
  file_ << "Uncalibrated - Calibrated m/z [ppm]" << "\t"; // ok
  file_ << "Uncalibrated - Calibrated m/z [Da]" << "\t"; // ok
  file_ << "Mass Error [ppm]" << "\t"; // sehr unterschiedlich (0.8 vs 300000), aber nicht klar was falsch
  file_ << "Mass Error [Da]" << "\t"; // sehr unterschiedlich, aber nicht klar was falsch
  file_ << "Uncalibrated Mass Error [ppm]" << "\t"; // ok
  file_ << "Uncalibrated Mass Error [Da]" << "\t"; //ok
  file_ << "Max intensity m/z 0" << "\t"; // unklar, Summed up eXtracted Ion Current (XIC) of all isotopic clusters associated with the identified AA sequence. In case of a labeled experiment this is the total intensity of all the isotopic patterns in the label cluster.
  file_ << "Retention time" << "\t"; // sehr unterschiedlich, könnte aber passen
  file_ << "Retention length" << "\t"; // ok
  file_ << "Calibrated retention time" << "\t"; // ok
  file_ << "Calibrated retention time start" << "\t"; // ok
  file_ << "Calibrated retention time finish" << "\t"; // ok
  file_ << "Retention time calibration" << "\t"; // sehr unterschiedlich (0.1 vs ca. 15000)
  file_ << "Match time difference" << "\t"; // ok
  file_ << "Match m/z difference" << "\t"; // ok
  file_ << "Match q-value" << "\t"; // vermutlich nicht in OpenMS
  file_ << "Match score" << "\t"; // vermutlich nicht in OpenMS (andromeda)
  file_ << "Number of data points" << "\t"; // number of peak centroids collected for this peptide feature (vermutlich in PeakPickerCWT)
  file_ << "Number of scans" << "\t"; // number of MS scans that the 3d peaks of this peptide feature are overlapping with, (simpel falls MSMS.txt implementiert)
  file_ << "Number of isotopic peaks" << "\t"; // number of isotopic peaks contained in this peptide feature (Feature, convexHull) - erledigt
  file_ << "PIF" << "\t"; // Parent Ion Fraction; indicates the fraction the target peak makes up of the total intensity in the inclusion window (mögl. MaxQuant spezifisch)
  file_ << "Fraction of total spectrum" << "\t"; // The percentage the ion intensity makes up of the total intensity of the whole spectrum (wahrscheinlich in OpenMS)
  file_ << "Base peak fraction" << "\t"; // The percentage the parent ion intensity in comparison to the highest peak in the MS spectrum (wahrscheinlich in OpenMS)
  file_ << "PEP" << "\t"; // Posterior Error Probability of the identification. This value essentially operates as a p-value, where smaller is more significant (wahrscheinlich in OpenMS)
  file_ << "MS/MS Count" << "\t"; // (wahrscheinlich) ok
  file_ << "MS/MS Scan Number" << "\t"; // (simpel falls MSMS.txt implementiert), Metadaden vom file?
  file_ << "Score" << "\t"; // wahrscheinlich falsch (ca. 0.02 (OpenMS) vs ca. 100 (MQ)), Andromeda benötigt
  file_ << "Delta score" << "\t"; // wahrscheinlich nicht möglich (Andromeda benötigt)
  file_ << "Combinatorics" << "\t"; // Number of possible distributions of the modifications over the peptide sequence., unklar, ob in OpenMS
  file_ << "Intensity" << "\t"; // Summed up eXtracted Ion Current (XIC) of all isotopic clusters associated with the identified AA sequence.
                                // In case of a labeled experiment this is the total intensity of all the isotopic patterns in the label cluster
  file_ << "Reporter intensity 0" << "\t"; // unklar, was das ist und wo es herkommt (keine Dokumentation in MQ), Ionen?
  file_ << "Reporter intensity 1" << "\t";
  file_ << "Reporter intensity 2" << "\t";
  file_ << "Reporter intensity 3" << "\t";
  file_ << "Reporter intensity 4" << "\t";
  file_ << "Reporter intensity 5" << "\t";
  file_ << "Reporter intensity not corrected 0" << "\t"; // unklar, was das ist und wo es herkommt (keine Dokumentation in MQ)
  file_ << "Reporter intensity not corrected 1" << "\t";
  file_ << "Reporter intensity not corrected 2" << "\t";
  file_ << "Reporter intensity not corrected 3" << "\t";
  file_ << "Reporter intensity not corrected 4" << "\t";
  file_ << "Reporter intensity not corrected 5" << "\t";
  file_ << "Reporter PIF" << "\t"; // wahrscheinlich falsch (NA vs ca. 0.8)
  file_ << "Reporter fraction" << "\t";
  file_ << "Reverse" << "\t"; // ok
  file_ << "Potential contaminant" << "\t"; // ok
  file_ << "id" << "\t"; // ok
  file_ << "Protein group IDs" << "\t"; // ok, aber letztes ; weg
  file_ << "Peptide ID" << "\t"; // The identifier of the non-redundant peptide sequence, (in FeatureMap, wenn das die richtige id ist)
  file_ << "Mod. peptide ID" << "\t"; // unklar (vielleicht selbst implementierbar), in 'modificationSpecificPeptides.txt'
  file_ << "MS/MS IDs" << "\t"; // Identifier(s) of the associated MS/MS summary(s) stored in the file 'msms.txt', unklar (vielleicht selbst implementierbar)
  file_ << "Best MS/MS" << "\t"; // Identifier(s) of the best MS/MS associated spectrum stored in the file 'msms.txt', unklar (vielleicht selbst implementierbar)
  file_ << "AIF MS/MS IDs" << "\t"; //Identifier(s) of the associated All Ion Fragmentation MS/MS summary(s) stored in the file 'aifMsms.txt'., unklar (vielleicht selbst implementierbar)
  file_ << "Oxidation (M) site IDs" << "\n"; // unklar (vielleicht selbst implementierbar)
}

Size MQEvidence::proteinGroupID_(const String& protein_accession)
{
  auto it = protein_id_.find(protein_accession);
  if (it == protein_id_.end())
  {
    protein_id_.emplace(protein_accession, protein_id_.size() + 1);
    return protein_id_.size();
  }
  else
  {
    return it->second;
  }
}

std::map<Size, Size> MQEvidence::makeFeatureUIDtoConsensusMapIndex_(const ConsensusMap& cmap)
{
  std::map<Size, Size> f_to_ci;
  for (Size i = 0; i < cmap.size(); ++i)
  {
    for (const auto& fh : cmap[i].getFeatures())
    {
      auto[it, was_created_newly] = f_to_ci.emplace(fh.getUniqueId(), i);
      if (!was_created_newly)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Adding [" + String(it->first) + "," + String(it->second) +  "] failed. FeatureHandle exists twice in ConsensusMap!");
      }
      f_to_ci[fh.getUniqueId()] = i;
    }
  }
  return f_to_ci;
}

bool MQEvidence::hasValidPepID_(
        const Feature& f,
        const Size c_feature_number,
        const std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>>& UIDs,
        const ProteinIdentification::Mapping& mp_f)
{
  const std::vector<PeptideIdentification>& pep_ids_f = f.getPeptideIdentifications();
  if (pep_ids_f.empty())
  {
    return false;
  }
  const PeptideIdentification& best_pep_id = pep_ids_f[0]; // PeptideIdentifications are sorted
  String best_uid = PeptideIdentification::buildUIDFromPepID(best_pep_id, mp_f.identifier_to_msrunpath);
  const auto range = UIDs.equal_range(best_uid);
  for (std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>>::const_iterator it_pep = range.first;
       it_pep != range.second; ++it_pep)
  {
    if (c_feature_number == it_pep->second.first)
    {
      return !pep_ids_f[0].getHits().empty(); // checks if PeptideIdentification has at least one hit
    }
  }
  return false;
}

bool MQEvidence::hasPeptideIdentifications_(const ConsensusFeature& cf)
{
  const std::vector<PeptideIdentification>& pep_ids_c = cf.getPeptideIdentifications();
  if (!pep_ids_c.empty())
  {
    return !pep_ids_c[0].getHits().empty(); // checks if PeptideIdentification has at least one hit
  }
  return false;
}

void MQEvidence::exportRowFromFeature_(
        const Feature& f,
        const ConsensusMap& cmap,
        const Size c_feature_number,
        const String& raw_file,
        const std::multimap<String, std::pair<Size, Size>>& UIDs,
        const ProteinIdentification::Mapping& mp_f)
{
  const PeptideHit* ptr_best_hit; // the best hit referring to score
  const ConsensusFeature& cf = cmap[c_feature_number];
  Size pep_ids_size = 0;
  String type;
  if (hasValidPepID_(f, c_feature_number, UIDs, mp_f))
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
  else if (hasPeptideIdentifications_(cf))
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

  file_ << pep_seq.toUnmodifiedString() << "\t"; // Sequence
  file_ << pep_seq.size() << "\t";               // Length

  std::map<String, Size> modifications;
  if (pep_seq.hasNTerminalModification())
  {
    const String& n_terminal_modification = pep_seq.getNTerminalModificationName();
    modifications.emplace(std::make_pair(n_terminal_modification, 1));
  }
  if (pep_seq.hasCTerminalModification())
  {
    modifications.emplace(std::make_pair(pep_seq.getCTerminalModificationName(), 1));
  }
  for (Size i = 0; i < pep_seq.size(); ++i)
  {
    if (pep_seq.getResidue(i).isModified())
    {
      ++modifications[pep_seq.getResidue(i).getModification()->getFullId()];
    }
  }

  if (modifications.empty())
  {
    file_ << "Unmodified"
          << "\t";
  }
  else
  {
    for (const auto& m : modifications)
    {
      file_ << m.first << ";"; // Modification
    }
    file_ << "\t";
  }
  file_ << "_" << pep_seq << "_" << "\t"; // Modified Sequence

  file_ << "NA" << "\t"; // Oxidation (M) Probabilities
  file_ << "NA" << "\t"; // Oxidation (M) Score Diffs

  if (pep_seq.hasNTerminalModification())
  {
    const String& n_terminal_modification = pep_seq.getNTerminalModificationName();
    n_terminal_modification.hasSubstring("Acetyl") ? file_ << 1 << "\t" : file_ << 0 << "\t"; // Acetyl (Protein N-term)
  }
  else
  {
    file_ << 0 << "\t"; // Acetyl (Protein N-term)
  }
    modifications.find("Oxidation (M)") == modifications.end() ? file_ << "0"
                                                                       << "\t" :
                                                                 file_ << modifications.find("Oxidation (M)")->second << "\t";
    file_ << ptr_best_hit->getMetaValue("missed_cleavages", "NA") << "\t"; // missed cleavages
    const std::set<String>& accessions = ptr_best_hit->extractProteinAccessionsSet();
    for (const String& p : accessions)
    {
      file_ << p << ";"; // Protein
    }
    file_ << "\t";
    file_ << ptr_best_hit << "\t"; // Leading Proteins
    file_ << ptr_best_hit << "\t"; // Leading Razor Proteins
    file_ << "NA" << "\t"; // Gene Names
    file_ << "NA" << "\t"; // Proteins Names
    file_ << type << "\t"; // type

    file_ << raw_file << "\t"; // Raw File

    file_ << "NA" << "\t"; // Fraction
    file_ << "NA" << "\t"; // MS/MS m/z

    file_ << f.getCharge() << "\t";           // Charge
    file_ << f.getMZ() << "\t";               // MZ
    file_ << pep_seq.getMonoWeight() << "\t"; // Mass
    file_ << f.getWidth() / 60 << "\t";       // Resolution in min.

    const double& uncalibrated_mz_error_ppm = ptr_best_hit->getMetaValue("uncalibrated_mz_error_ppm", NAN);
    const double& calibrated_mz_error_ppm = ptr_best_hit->getMetaValue("calibrated_mz_error_ppm", NAN);

    if (std::isnan(uncalibrated_mz_error_ppm) && std::isnan(calibrated_mz_error_ppm))
    {
      file_ << "NA"
            << "\t"; // Uncalibrated - Calibrated m/z [ppm]
      file_ << "NA"
            << "\t"; // Uncalibrated - Calibrated m/z [mDa]
      file_ << "NA"
            << "\t"; // Mass error [ppm]
      file_ << "NA"
            << "\t"; // Mass error [Da]
      file_ << "NA"
            << "\t"; // Uncalibrated Mass error [ppm]
      file_ << "NA"
            << "\t"; // Uncalibrated Mass error [Da]
    }
    else if (std::isnan(calibrated_mz_error_ppm))
    {
      file_ << "NA"
            << "\t"; // Uncalibrated - Calibrated m/z [ppm]
      file_ << "NA"
            << "\t"; // Uncalibrated - Calibrated m/z [mDa]
      file_ << "NA"
            << "\t"; // Mass error [ppm]
      file_ << "NA"
            << "\t";                                                                  // Mass error [Da]
      file_ << uncalibrated_mz_error_ppm << "\t";                                     // Uncalibrated Mass error [ppm]
      file_ << OpenMS::Math::ppmToMass(uncalibrated_mz_error_ppm, f.getMZ()) << "\t"; // Uncalibrated Mass error [Da]
    }
    else if (std::isnan(uncalibrated_mz_error_ppm))
    {
      file_ << "NA"
            << "\t"; // Uncalibrated - Calibrated m/z [ppm]
      file_ << "NA"
            << "\t";                                                                // Uncalibrated - Calibrated m/z [mDa]
      file_ << calibrated_mz_error_ppm << "\t";                                     // Mass error [ppm]
      file_ << OpenMS::Math::ppmToMass(calibrated_mz_error_ppm, f.getMZ()) << "\t"; // Mass error [Da]
      file_ << "NA"
            << "\t"; // Uncalibrated Mass error [ppm]
      file_ << "NA"
            << "\t"; // Uncalibrated Mass error [Da]
    }
    else
    {
      file_ << uncalibrated_mz_error_ppm - calibrated_mz_error_ppm << "\t";                                       // Uncalibrated - Calibrated m/z [ppm]
      file_ << OpenMS::Math::ppmToMass((uncalibrated_mz_error_ppm - calibrated_mz_error_ppm), f.getMZ()) << "\t"; // Uncalibrated - Calibrated m/z [Da]
      file_ << calibrated_mz_error_ppm << "\t";                                                                   // Mass error [ppm]
      file_ << OpenMS::Math::ppmToMass(calibrated_mz_error_ppm, f.getMZ()) << "\t";                               // Mass error [Da]
      file_ << uncalibrated_mz_error_ppm << "\t";                                                                 // Uncalibrated Mass error [ppm]
      file_ << OpenMS::Math::ppmToMass(uncalibrated_mz_error_ppm, f.getMZ()) << "\t";                             // Uncalibrated Mass error [Da]
    }

    file_ << "NA" << "\t"; // Max intensity m/z 0
    file_ << f.getRT() / 60 << "\t"; // Retention time in min.

    // RET LENGTH
    f.metaValueExists("rt_raw_end") && f.metaValueExists("rt_raw_start") ?
      file_ << (double(f.getMetaValue("rt_raw_end")) - double(f.getMetaValue("rt_raw_start"))) / 60 << "\t" : file_
            << "NA" << "\t";

    if (f.metaValueExists("rt_align"))
    {
      file_ << double(f.getMetaValue("rt_align")) / 60 << "\t";   // Calibrated Retention Time
    }
    f.metaValueExists("rt_align_start") ? file_ << double(f.getMetaValue("rt_align_start")) / 60 << "\t" :
                                          file_ << "NA" << "\t"; //  Calibrated retention time start
    f.metaValueExists("rt_align_end") ? file_ << double(f.getMetaValue("rt_align_end")) / 60 << "\t" :
                                        file_ << "NA" << "\t"; // Calibrated retention time end
    if (f.metaValueExists("rt_align"))
    {
      file_ << (f.getRT() - double(f.getMetaValue("rt_align"))) / 60 << "\t"; // Retention time calibration
    }
    else
    {
      file_ << "NA" << "\t"; // calibrated retention time
      file_ << "NA" << "\t"; // Retention time calibration
    }

    if (type == "MULTI-MSMS")
    {
      file_ << "NA" << "\t"; // Match time diff
      file_ << "NA" << "\t"; // Match mz diff
    }
    else
    {
      f.metaValueExists("rt_align") ? file_ << double(f.getMetaValue("rt_align")) - cmap[c_feature_number].getRT() << "\t" :
                                      file_ << "NA"
                                            << "\t";               // Match time diff
      file_ << f.getMZ() - cmap[c_feature_number].getMZ() << "\t"; // Match mz diff
    }

    file_ << "NA" << "\t"; // Match q-value
    file_ << "NA" << "\t"; // Match score
    /*PeakPickerCWT::PeakArea_& area;
    PeakPickerCWT::getPeakCentroid_(PeakPickerCWT::PeakArea_& area)
    file_ << area.size() << "\t"; // Number of data points*/ // -> Problem: Sachen sind protected
    file_ << "NA" << "\t"; // Number of data points
    file_ << "NA" << "\t"; // Number of scans
    file_ << f.getConvexHulls().size() << "\t"; // Number of isotopic peaks
    file_ << "NA" << "\t"; // PIF
    file_ << "NA" << "\t"; // Fraction of total spectrum
    file_ << "NA" << "\t"; // Base peak fraction
    file_ << "NA" << "\t"; // PEP

    file_ << pep_ids_size << "\t"; // MS/MS count
    file_ << "NA" << "\t"; // MS/MS Scan Number
    file_ << max_score << "\t"; // Score

    file_ << "NA" << "\t"; // Delta score
    file_ << "NA" << "\t"; // Combinatorics

    file_ << f.getIntensity() << "\t"; // Intensity

    file_ << "NA" << "\t"; // Reporter intensity 0
    file_ << "NA" << "\t"; // Reporter intensity 1
    file_ << "NA" << "\t"; // Reporter intensity 2
    file_ << "NA" << "\t"; // Reporter intensity 3
    file_ << "NA" << "\t"; // Reporter intensity 4
    file_ << "NA" << "\t"; // Reporter intensity 5
    file_ << "NA" << "\t"; // Reporter intensity not corrected 0
    file_ << "NA" << "\t"; // Reporter intensity not corrected 1
    file_ << "NA" << "\t"; // Reporter intensity not corrected 2
    file_ << "NA" << "\t"; // Reporter intensity not corrected 3
    file_ << "NA" << "\t"; // Reporter intensity not corrected 4
    file_ << "NA" << "\t"; // Reporter intensity not corrected 5
    file_ << "NA" << "\t"; // Reporter PIF
    file_ << "NA" << "\t"; // Reporter fraction


    ptr_best_hit->getMetaValue("target_decoy") == "decoy" ? file_ << "1"
                                                                  << "\t" :
                                                            file_ << "\t"; // reverse


    String pot_containment = ptr_best_hit->getMetaValue("is_contaminant", "NA");
    if (pot_containment == "1")
    {
      file_ << "+"
            << "\t"; // Potential contaminant
    }
    else
    {
      file_ << "\t";
    }

    file_ << id_ << "\t"; // ID
    ++id_;

    for (const String& p : accessions)
    {
      file_ << proteinGroupID_(p) << ";"; // Protein group ids
    }
    file_ << "\t";

    file_ << "NA" << "\t"; // Peptide ID
    file_ << "NA" << "\t"; // Mod peptide ID
    file_ << "NA" << "\t"; // MS/MS IDs
    file_ << "NA" << "\t"; // Best MS/MS
    file_ << "NA" << "\t"; // AIF MS/MS IDs
    file_ << "NA" << "\n"; // Oxidation site IDs

}

void MQEvidence::exportFeatureMap(const FeatureMap& feature_map, const ConsensusMap& cmap)
{
    if (!isValid())
    {
      OPENMS_LOG_ERROR << "MqEvidence object is not valid." << std::endl;
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_);
    }
    const std::map<Size, Size>& fTc = makeFeatureUIDtoConsensusMapIndex_(cmap);
    StringList spectra_data;
    feature_map.getPrimaryMSRunPath(spectra_data);
    String raw_file = File::basename(spectra_data.empty() ? feature_map.getLoadedFilePath() : spectra_data[0]);

    ProteinIdentification::Mapping mp_f;
    mp_f.create(feature_map.getProteinIdentifications());

    std::multimap<String, std::pair<Size, Size>> UIDs = PeptideIdentification::buildUIDsFromAllPepIDs(cmap);

    for (const Feature& f : feature_map)
    {
      const Size& f_id = f.getUniqueId();
      const auto& c_id = fTc.find(f_id);
      if (c_id != fTc.end())
      {
        exportRowFromFeature_(f, cmap, c_id->second, raw_file, UIDs, mp_f);
      }
      else
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Feature in FeatureMap has no associated ConsensusFeature.");
      }
    }
    file_.flush();
}
