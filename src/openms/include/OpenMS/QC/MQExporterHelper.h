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
// $Authors: Virginia Rossow, Lenny Kovac, Hendrik Beschorner$
// --------------------------------------------------------------------------

#pragma once

#include <fstream>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

class OPENMS_DLLAPI MQExporterHelper
/**
@brief Helper class for common functions and values needed for exporting MaxQuant outputs

@ingroup Metadata
*/
{
public:

  struct MQCommonOutputs
  {
    std::stringstream modifications;
    std::stringstream modified_sequence;
    std::stringstream acetyl;
    std::stringstream oxidation;
    std::stringstream gene_names;
    std::stringstream protein_names;
    std::stringstream msms_mz;
    std::stringstream mass_error_ppm;
    std::stringstream mass_error_da;
    std::stringstream uncalibrated_mass_error_ppm;
    std::stringstream uncalibrated_mass_error_da;
    std::stringstream uncalibrated_calibrated_mz_ppm;
    std::stringstream uncalibrated_calibrated_mz_mda;

    //file_ << "Sequence" << "\t"; maybe, relativ trivial
    //file_ << "Length" << "\t";
    //file_ << "Modifications" << "\t"; implementieren
    // file_ << "Modified sequence" << "\t"; implementieren
    //file_ << "Acetyl (Protein N-term)" << "\t"; implementieren
    //file_ << "Oxidation (M)" << "\t"; implementieren
    //file_ << "Missed cleavages" << "\t"; trivial
    //file_ << "Proteins" << "\t"; trivial
    //file_ << "Gene Names" << "\t"; // in progress, aber implementieren
    //file_ << "Protein Names" << "\t"; // in progress, aber implementieren
    //file_ << "Type" << "\t"; TODO different type
    //file_ << "Raw file" << "\t"; trivial
    //file_ << "MS/MS m/z" << "\t"; implementieren TODO is m/z in MSMS MS/MS m/z
    //file_ << "Charge" << "\t"; trivial
    //file_ << "m/z" << "\t"; trivial TODO
    //file_ << "Mass" << "\t"; trivial
    //file_ << "Mass Error [ppm]" << "\t"; vielleicht, beim einen halt noch calibrated dabei
    //file_ << "Mass Error [Da]" << "\t"; vielleicht, beim einen halt noch calibrated dabei
    //file_ << "Retention time" << "\t"; trivial
    //file_ << "Fraction of total spectrum" << "\t"; trivial
    //file_ << "Base peak fraction" << "\t"; trvial
    //file_ << "PEP" << "\t"; trivial
    //file_ << "MS/MS Scan Number" << "\t"; trivial
    //file_ << "Score" << "\t"; trivial
    //file_ << "Delta score" << "\t"; trivial
    //file_ << "Reverse" << "\t";
    //file_ << "id" << "\t"; ?
    //file_ << "Protein group IDs" << "\n"; trivial

    MQCommonOutputs(
      const OpenMS::Feature& f,
      const OpenMS::ConsensusMap& cmap,
      const OpenMS::Size c_feature_number,
      const OpenMS::String& raw_file,
      const std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>>& UIDs,
      const OpenMS::ProteinIdentification::Mapping& mp_f,
      const OpenMS::MSExperiment& exp,
      const std::map<OpenMS::String,OpenMS::String>& prot_mapper)
    {
      const OpenMS::PeptideHit* ptr_best_hit; // the best hit referring to score
      const OpenMS::ConsensusFeature& cf = cmap[c_feature_number];
      if (hasValidPepID_(f, c_feature_number, UIDs, mp_f))
      {
        ptr_best_hit = &f.getPeptideIdentifications()[0].getHits()[0];
      }
      else if (hasPeptideIdentifications_(cf))
      {
        ptr_best_hit = &cf.getPeptideIdentifications()[0].getHits()[0];
      }
      else
      {
        return; // no valid PepID; nothing to export
      }

      const OpenMS::AASequence& pep_seq = ptr_best_hit->getSequence();
      if (pep_seq.empty())
      {
        return; // empty AASequence; nothing to export
      }

      //get all peptide-evidences for the best hit
      const auto& pep_evidences = ptr_best_hit->getPeptideEvidences();

      std::map<OpenMS::String, OpenMS::Size> modifications_temp;
      if (pep_seq.hasNTerminalModification())
      {
        const OpenMS::String& n_terminal_modification = pep_seq.getNTerminalModificationName();
        modifications_temp.emplace(std::make_pair(n_terminal_modification, 1));
      }
      if (pep_seq.hasCTerminalModification())
      {
        modifications_temp.emplace(std::make_pair(pep_seq.getCTerminalModificationName(), 1));
      }
      for (OpenMS::Size i = 0; i < pep_seq.size(); ++i)
      {
        if (pep_seq.getResidue(i).isModified())
        {
          ++modifications_temp[pep_seq.getResidue(i).getModification()->getFullId()];
        }
      }

      if (modifications_temp.empty())
      {
        modifications.str("Unmodified");
      }
      else
      {
        auto it = modifications_temp.begin();
        modifications.str(it->first);
        ++it;
        for (; it != modifications_temp.end(); ++it)
        {
          modifications << ";" << it->first;
        }
      }
      modified_sequence.clear();
      modified_sequence << "_" << pep_seq << "_"; // Modified Sequence

      if (pep_seq.hasNTerminalModification())
      {
        const OpenMS::String& n_terminal_modification = pep_seq.getNTerminalModificationName();
        n_terminal_modification.hasSubstring("Acetyl") ? acetyl.str("1") : acetyl.str("0"); // Acetyl (Protein N-term)
      }
      else
      {
        acetyl.str("0");
      }
      oxidation.clear();
      modifications_temp.find("Oxidation (M)") == modifications_temp.end() ? oxidation << "0" : oxidation << modifications_temp.find("Oxidation (M)")->second;

      const auto& prot_mapper_it = prot_mapper.find(pep_evidences[0].getProteinAccession());
      if(prot_mapper_it == prot_mapper.end())
      {
        gene_names.str("NA"); //Gene Names
        protein_names.str("NA"); //Protein Names
      }
      else
      {
        auto protein_description = prot_mapper_it->second;
        auto gene_name = protein_description.substr(protein_description.find("GN=") + 3);
        gene_name = gene_name.substr(0,gene_name.find(" "));

        gene_names.str(gene_name);           //Gene Names
        protein_names.str(protein_description); //Protein Names
      }

      msms_mz.clear();
      if (f.metaValueExists("spectrum_index") && !exp.empty() && exp.getNrSpectra() >= (OpenMS::Size)f.getMetaValue("spectrum_index") && !exp[f.getMetaValue("spectrum_index")].empty())
      {
        const OpenMS::MSSpectrum& ms2_spec = exp[f.getMetaValue("spectrum_index")];
        if (!ms2_spec.getPrecursors().empty())
        {
          msms_mz << ms2_spec.getPrecursors()[0].getMZ(); // MS/MS m/z
        }
      }
      else
      {
        msms_mz.str("");
      }

      const double& uncalibrated_mz_error_ppm = ptr_best_hit->getMetaValue("uncalibrated_mz_error_ppm", NAN);
      const double& calibrated_mz_error_ppm = ptr_best_hit->getMetaValue("calibrated_mz_error_ppm", NAN);

      if (std::isnan(uncalibrated_mz_error_ppm) && std::isnan(calibrated_mz_error_ppm))
      {
        uncalibrated_calibrated_mz_ppm.str("NA"); // Uncalibrated - Calibrated m/z [ppm]
        uncalibrated_calibrated_mz_mda.str("NA"); // Uncalibrated - Calibrated m/z [mDa]
        mass_error_ppm.str("NA"); // Mass error [ppm]
        mass_error_da.str("NA"); // Mass error [Da]
        uncalibrated_mass_error_ppm.str("NA"); // Uncalibrated Mass error [ppm]
        uncalibrated_mass_error_da.str("NA"); // Uncalibrated Mass error [Da]
      }
      else if (std::isnan(calibrated_mz_error_ppm))
      {
        uncalibrated_calibrated_mz_ppm.str("NA"); // Uncalibrated - Calibrated m/z [ppm]
        uncalibrated_calibrated_mz_mda.str("NA"); // Uncalibrated - Calibrated m/z [mDa]
        mass_error_ppm.str("NA"); // Mass error [ppm]
        mass_error_da.str("NA"); // Mass error [Da]
        uncalibrated_mass_error_ppm.clear();
        uncalibrated_mass_error_da.clear();
        uncalibrated_mass_error_ppm << uncalibrated_mz_error_ppm;                      // Uncalibrated Mass error [ppm]
        uncalibrated_mass_error_da << OpenMS::Math::ppmToMass(uncalibrated_mz_error_ppm, f.getMZ()); // Uncalibrated Mass error [Da]
      }
      else if (std::isnan(uncalibrated_mz_error_ppm))
      {
        uncalibrated_calibrated_mz_ppm.str("NA"); // Uncalibrated - Calibrated m/z [ppm]
        uncalibrated_calibrated_mz_mda.str("NA"); // Uncalibrated - Calibrated m/z [mDa]
        mass_error_ppm.clear();
        mass_error_da.clear();
        mass_error_ppm << calibrated_mz_error_ppm;                                     // Mass error [ppm]
        mass_error_da << OpenMS::Math::ppmToMass(calibrated_mz_error_ppm, f.getMZ()); // Mass error [Da]
        uncalibrated_mass_error_ppm.str("NA"); // Uncalibrated Mass error [ppm]
        uncalibrated_mass_error_da.str("NA"); // Uncalibrated Mass error [Da]
      }
      else
      {
        uncalibrated_calibrated_mz_ppm.clear(); // Uncalibrated - Calibrated m/z [ppm]
        uncalibrated_calibrated_mz_mda.clear(); // Uncalibrated - Calibrated m/z [mDa]
        mass_error_ppm.clear(); // Mass error [ppm]
        mass_error_da.clear(); // Mass error [Da]
        uncalibrated_mass_error_ppm.clear(); // Uncalibrated Mass error [ppm]
        uncalibrated_mass_error_da.clear(); // Uncalibrated Mass error [Da]

        uncalibrated_calibrated_mz_ppm << uncalibrated_mz_error_ppm - calibrated_mz_error_ppm; // Uncalibrated - Calibrated m/z [ppm]
        uncalibrated_calibrated_mz_mda << OpenMS::Math::ppmToMass((uncalibrated_mz_error_ppm - calibrated_mz_error_ppm), f.getMZ()); // Uncalibrated - Calibrated m/z [Da]
        mass_error_ppm << calibrated_mz_error_ppm; // Mass error [ppm]
        mass_error_da << OpenMS::Math::ppmToMass(calibrated_mz_error_ppm, f.getMZ()); // Mass error [Da]
        uncalibrated_mass_error_ppm << uncalibrated_mz_error_ppm; // Uncalibrated Mass error [ppm]
        uncalibrated_mass_error_da << OpenMS::Math::ppmToMass(uncalibrated_mz_error_ppm, f.getMZ()); // Uncalibrated Mass error [Da]
      }
    }
  };

  MQCommonOutputs mq_common_outputs; ///< struct for storing values, needed for multiple MaxQuant-exporters

  /**
    @brief returns the MaxQuant unique evidence number of a protein accession

      Obtains a unique, consecutive number for each distinct protein, which can
      be used as a protein ID in the MaxQuant output files (in lack of a proper
      proteingroup ID which maps to proteinGroups.txt)


    @param protein_accession The accession of the protein

    @return Returns distinct number for every Protein
  */
  OpenMS::Size proteinGroupID_(std::map<OpenMS::String, OpenMS::Size>& protein_id_,
                               const OpenMS::String& protein_accession);

  /**
    @brief Creates map that has the information which FeatureUID is mapped to which ConsensusFeature in ConsensusMap

    @throw Exception::Precondition if FeatureHandle exists twice in ConsensusMap

    @param cmap ConsensusMap that includes ConsensusFeatures

    @return Returns map, the index is a FeatureID, the value is the index of the ConsensusFeature
    in the vector of ConsensusMap
  */
  std::map<OpenMS::Size, OpenMS::Size> makeFeatureUIDtoConsensusMapIndex_(const OpenMS::ConsensusMap& cmap);

  /**
    @brief Checks if Feature has valid PeptideIdentifications

      If there are no PeptideIdentifications or the best hit of the Feature cannot be found in corresponding ConsensusFeature,
      the functions returns false to show that something went wrong.

    @param f Feature to extract PeptideIdentifications
    @param c_feature_number Index of corresponding ConsensusFeature in ConsensusMap
    @param UIDs UIDs of all PeptideIdentifications of the ConsensusMap
    @param mp_f Mapping between the FeatureMap and ProteinIdentifications for the UID

    @return Returns true if the PeptideIdentifications exist and are valid
  */
  static bool hasValidPepID_(
    const OpenMS::Feature& f,
    const OpenMS::Size c_feature_number,
    const std::multimap<OpenMS::String, std::pair<OpenMS::Size, OpenMS::Size>>& UIDs,
    const OpenMS::ProteinIdentification::Mapping& mp_f);

  /**
    @brief Checks if ConsensusFeature has valid PeptideIdentifications

    If there are no PeptideIdentifications,
    the functions returns false to show that something went wrong.

    @param cf is used to extract PeptideIdentifications

    @return Returns true if the ConsensusFeature has any PepIDs; otherwise false
  */
  static bool hasPeptideIdentifications_(const OpenMS::ConsensusFeature& cf);



  /**
    @brief Creates the helperclass and computes the (in multiple MaxQuant output files) common output values, which are not trivial

          The values are stored in the struct MQCommonOutputs

    @throw Exception::FileNotWritable if evidence.txt could not be created

    @param path that is the path where evidence.txt has to be stored

  */
  //explicit MQExporterHelper();

  /**
    @brief Destructor
  */
  ~MQExporterHelper();

  /**
      @brief Checks if evidence.txt is writable
             (i.e. the path in the ctor was not empty and could be created)

      @return Returns true if evidence.txt is writable
  */
  bool isValid(std::string filename_);
};