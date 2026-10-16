// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// ProSEResultSet.tables() (bound in bind_prose.cpp): the pooled results of a search as compact typed columns, three
// tables (PSMs, peptides, proteins) whose rows refer to each other by index, for a browser app with 100,000+ PSMs.
// ProSEResultSet.spectra(): the peaks of the spectra of the PSM rows, which the search reads from its input files.
// The columns are copied into JavaScript arrays, which JavaScript owns.

#include "bind_prose.h"
#include "binding_utils.h"

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <initializer_list>
#include <limits>
#include <set>
#include <optional>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace OpenMS::WASM
{
  namespace
  {
    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

    /// The name of the search engine score of ProSE, also the name of its feature column (Percolator's "score")
    const std::string HYPERSCORE = "ln(hyperscore)";

    /// The index of meta value @p name in the registry, for fast lookups of a value of many PSMs
    UInt metaIndex(const std::string& name)
    {
      return MetaInfoInterface::metaRegistry().registerName(name);
    }

    /// The meta value @p key of @p object as number: an int or double, or a string that spells a number (as
    /// Percolator reads its features); NaN if it is missing or something else.
    double metaNumber(const MetaInfoInterface& object, UInt key)
    {
      if (!object.metaValueExists(key)) return NaN;
      const DataValue& v = object.getMetaValue(key);
      switch (v.valueType())
      {
        case DataValue::INT_VALUE:
        case DataValue::DOUBLE_VALUE:
          return static_cast<double>(v);
        case DataValue::STRING_VALUE:
        {
          const std::string s = v.toString();
          char* end = nullptr;
          const double d = std::strtod(s.c_str(), &end);
          return (!s.empty() && end == s.c_str() + s.size()) ? d : NaN;
        }
        default:
          return NaN;
      }
    }

    /// The first of the meta values @p keys of @p object that is a number, or NaN
    double firstMetaNumber(const MetaInfoInterface& object, std::initializer_list<UInt> keys)
    {
      for (const UInt key : keys)
      {
        const double d = metaNumber(object, key);
        if (!std::isnan(d)) return d;
      }
      return NaN;
    }

    /// The meta value indices of the scores of a PSM
    struct ScoreKeys
    {
      UInt hyperscore = metaIndex(HYPERSCORE);                       // after Percolator (the score it replaced)
      UInt hyperscore_score = metaIndex(HYPERSCORE + "_score");      // after FDR:PSM (the score the q-value replaced)
      UInt percolator_score = metaIndex("percolator_score");
      UInt percolator_score_cv = metaIndex("MS:1001492");
      UInt percolator_q = metaIndex("percolator_q_value");
      UInt q_cv = metaIndex("MS:1001491");
      UInt q = metaIndex("q-value");
      UInt q_score = metaIndex("q-value_score");
      UInt percolator_pep = metaIndex("percolator_pep");
      UInt pep_cv = metaIndex("MS:1001493");
      UInt pep_score = metaIndex("Posterior Error Probability_score");
      UInt precursor_error = metaIndex(Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM);
      UInt isotope_error = metaIndex(Constants::UserParam::ISOTOPE_ERROR);
      UInt delta_mass = metaIndex("DeltaMass");
      UInt merge_index = metaIndex(Constants::UserParam::ID_MERGE_INDEX);
    };

    /// ProSE's ln(hyperscore) of a hit of @p id: the score of the identification, or the meta value that kept it
    /// when Percolator or the FDR:PSM replaced it
    double hyperscore(const PeptideIdentification& id, const PeptideHit& hit, const ScoreKeys& k)
    {
      if (id.getScoreType() == HYPERSCORE) return hit.getScore();
      return firstMetaNumber(hit, {k.hyperscore, k.hyperscore_score});
    }

    double qValue(const PeptideIdentification& id, const PeptideHit& hit, const ScoreKeys& k)
    {
      if (id.getScoreType() == "q-value") return hit.getScore();
      return firstMetaNumber(hit, {k.percolator_q, k.q_cv, k.q, k.q_score});
    }

    double pep(const PeptideIdentification& id, const PeptideHit& hit, const ScoreKeys& k)
    {
      if (id.getScoreType() == "Posterior Error Probability") return hit.getScore();
      return firstMetaNumber(hit, {k.percolator_pep, k.pep_cv, k.pep_score});
    }

    /// Percolator's (SVM) score of a hit, NaN if it was not rescored
    double percolatorScore(const PeptideIdentification& id, const PeptideHit& hit, const ScoreKeys& k)
    {
      if (id.getScoreType() == "svm") return hit.getScore();
      return firstMetaNumber(hit, {k.percolator_score, k.percolator_score_cv});
    }

    /// The best hit of @p id (by its score and the orientation of the score)
    const PeptideHit& bestHit(const PeptideIdentification& id)
    {
      const auto& hits = id.getHits();
      const bool higher_better = id.isHigherScoreBetter();
      return *std::min_element(hits.begin(), hits.end(), [higher_better](const PeptideHit& a, const PeptideHit& b) {
        return higher_better ? a.getScore() > b.getScore() : a.getScore() < b.getScore();
      });
    }

    template <typename T>
    T clamped(double v)
    {
      if (std::isnan(v)) return 0;
      return static_cast<T>(std::clamp(v, static_cast<double>(std::numeric_limits<T>::min()),
                                       static_cast<double>(std::numeric_limits<T>::max())));
    }

    /// A row of tables().psms (and of spectra().psmSpectrum): an identification with hits, and its best hit
    struct PsmRow
    {
      Size id;                ///< the position of the identification in the pooled PSMs
      const PeptideHit* hit;  ///< its best hit
    };

    /// The PSM rows of the pooled PSMs @p ids: one per identification with hits, in their order
    std::vector<PsmRow> psmRows(const PeptideIdentificationList& ids)
    {
      std::vector<PsmRow> rows;
      rows.reserve(ids.size());
      for (Size i = 0; i < ids.size(); ++i)
      {
        if (!ids[i].getHits().empty()) rows.push_back({i, &bestHit(ids[i])});
      }
      return rows;
    }

    /// The input file of identification @p id: its position in the input files ('id_merge_index'; 0 if absent)
    uint16_t psmFile(const PeptideIdentification& id, const ScoreKeys& k)
    {
      return clamped<uint16_t>(metaNumber(id, k.merge_index));
    }

    /// A JavaScript array of @p strings. One string split in JavaScript is much faster than one conversion per
    /// string (100,000s of them); strings with a line break are converted one by one.
    val toStringArray(const std::vector<std::string>& strings)
    {
      if (strings.empty()) return val::array();
      Size length = strings.size();
      bool joinable = true;
      for (const std::string& s : strings)
      {
        length += s.size();
        joinable = joinable && s.find('\n') == std::string::npos;
      }
      if (!joinable) return toArray(strings);
      std::string joined;
      joined.reserve(length);
      for (Size i = 0; i < strings.size(); ++i)
      {
        if (i > 0) joined += '\n';
        joined += strings[i];
      }
      return val(joined).call<val>("split", std::string("\n"));
    }

    /// The modification as ProForma names it: its Unimod name ("Oxidation") if it is a Unimod modification. (ProForma
    /// of OpenMS writes Unimod accessions, "UNIMOD:35"; the app shows the names.)
    void useUnimodName(ProForma::Modification& m, const ResidueModification* mod)
    {
      if (mod == nullptr || mod->getUniModRecordId() <= 0 || mod->getId().empty()) return;
      m.alternatives.clear();
      m.alternatives.emplace_back(ProForma::NamedMod{std::nullopt, mod->getId()}, std::nullopt);
    }

    /// @p peptide in ProForma, with Unimod names: "AVM[Oxidation]DDFAAFVEK", "[Acetyl]-DAHKSEVAHR", "PEPTIDE-[Amidated]".
    /// Other modifications as ProForma of OpenMS writes them (e.g. a mass delta "[+15.9949]").
    std::string toProForma(const AASequence& peptide)
    {
      ProForma::Peptidoform pf = ProForma::fromAASequence(peptide); // one element per residue
      for (Size i = 0; i < peptide.size() && i < pf.sequence.size(); ++i)
      {
        auto* element = std::get_if<ProForma::SequenceElement>(&pf.sequence[i]);
        if (element != nullptr && !element->modifications.empty())
        {
          useUnimodName(element->modifications.front(), peptide[i].getModification());
        }
      }
      if (!pf.n_term_mods.empty()) useUnimodName(pf.n_term_mods.front(), peptide.getNTerminalModification());
      if (!pf.c_term_mods.empty()) useUnimodName(pf.c_term_mods.front(), peptide.getCTerminalModification());
      return ProForma::toString(pf);
    }

    /// The columns of the peptides
    struct PeptideColumns
    {
      std::vector<std::string> openms, proforma, sequence;
      std::vector<uint8_t> decoy;
      std::vector<uint32_t> mod_offsets{0};
      std::vector<int16_t> mod_position;
      std::vector<std::string> mod_name;
      std::vector<int32_t> mod_unimod;
      std::vector<double> mod_mass;
      std::vector<uint32_t> protein_offsets{0};
      std::vector<uint32_t> protein, start, end;

      void addModification(int position, const ResidueModification* mod)
      {
        if (mod == nullptr) return;
        mod_position.push_back(static_cast<int16_t>(position));
        // a modification without a name (e.g. mass-only) by its full name: "[+42.0106]"
        mod_name.push_back(!mod->getId().empty() ? mod->getId() : !mod->getFullName().empty() ? mod->getFullName() : mod->getFullId());
        mod_unimod.push_back(mod->getUniModRecordId() > 0 ? mod->getUniModRecordId() : -1);
        mod_mass.push_back(mod->getDiffMonoMass());
      }
    };

    /// The sequence of @p e as ProSE searches it: without a stop codon ('*') at its end
    std::string searchedSequence(const FASTAFile::FASTAEntry& e)
    {
      std::string s = e.sequence;
      while (!s.empty() && s.back() == '*') s.pop_back();
      return s;
    }

    /// The decoy ProSE generates from the target sequence @p target, as ProSEAlgorithm builds its decoy database: the
    /// peptides of the enzyme reversed (whole protein for unspecific searches), an initial M kept in place when the
    /// initial M is clipped. Empty if it cannot be built.
    /// Keep in sync with ProSEAlgorithm::buildDecoyAugmentedDB_ (ProSEAlgorithm.cpp:~2357), whose rule this copies (it
    /// is private, and the search does not return its decoy database).
    std::string generatedDecoy(std::string target, const ProSEResultSet::DecoyGeneration& g)
    {
      const bool keep_met = g.keep_initial_met && target.size() > 1 && target[0] == 'M';
      if (keep_met) target.erase(0, 1);
      try
      {
        const DecoyGenerator generator;
        const AASequence protein = AASequence::fromString(target);
        std::string decoy = g.unspecific ? generator.reverseProtein(protein).toString()
                                         : generator.reversePeptides(protein, g.enzyme).toString();
        if (keep_met) decoy.insert(decoy.begin(), 'M');
        return decoy;
      }
      catch (const Exception::BaseException&) // not a sequence AASequence reads (ProSE would not have searched it)
      {
        return {};
      }
    }
  } // namespace

  JS::ResultTables ProSEResultSet::tables() const
  {
    const ScoreKeys keys;
    const ProteinIdentification* run = merged.protein_ids.empty() ? nullptr : &merged.protein_ids[0];

    //------------------------------------------------------------------------------
    // Proteins: the hits of the pooled run

    const std::vector<ProteinHit> no_hits;
    const std::vector<ProteinHit>& hits = run != nullptr ? run->getHits() : no_hits;
    const Size n_proteins = hits.size();
    std::unordered_map<std::string, uint32_t> protein_row;
    protein_row.reserve(n_proteins);
    for (Size p = 0; p < n_proteins; ++p)
    {
      protein_row.emplace(hits[p].getAccession(), static_cast<uint32_t>(p));
    }
    std::unordered_map<std::string, const FASTAFile::FASTAEntry*> fasta;
    fasta.reserve(database.size());
    for (const FASTAFile::FASTAEntry& e : database) fasta.emplace(e.identifier, &e);

    const auto hasDecoyMarker = [this](const std::string& accession) {
      if (decoy_string.empty() || accession.size() < decoy_string.size()) return false;
      return decoy_is_prefix ? accession.compare(0, decoy_string.size(), decoy_string) == 0
                             : accession.compare(accession.size() - decoy_string.size(), decoy_string.size(), decoy_string) == 0;
    };
    std::vector<std::string> accession(n_proteins), description(n_proteins), protein_sequence(n_proteins);
    std::vector<double> protein_score(n_proteins, NaN), protein_q(n_proteins, NaN);
    std::vector<uint8_t> protein_decoy(n_proteins, 0);
    std::vector<int32_t> group(n_proteins, -1);
    const bool protein_q_values = merged.protein_fdr_applied && run != nullptr && run->getScoreType() == "q-value";
    for (Size p = 0; p < n_proteins; ++p)
    {
      const ProteinHit& hit = hits[p];
      accession[p] = hit.getAccession();
      const bool decoy = hit.getTargetDecoyType() == ProteinHit::TargetDecoyType::DECOY
                         || (hit.getTargetDecoyType() == ProteinHit::TargetDecoyType::UNKNOWN && hasDecoyMarker(accession[p]));
      protein_decoy[p] = decoy ? 1 : 0;
      if (protein_q_values) protein_q[p] = hit.getScore();
      description[p] = hit.getDescription();
      // The sequence as searched: of the FASTA, or of the decoy that ProSE generated from its target (generated decoys
      // are the decoy prefix + the accession of the target, and also replace decoys of the FASTA).
      const FASTAFile::FASTAEntry* entry = nullptr;
      bool generate = false;
      if (decoy && decoys.generated && decoy_is_prefix && hasDecoyMarker(accession[p]))
      {
        const auto it = fasta.find(accession[p].substr(decoy_string.size()));
        if (it != fasta.end()) { entry = it->second; generate = true; }
      }
      else
      {
        const auto it = fasta.find(accession[p]);
        if (it != fasta.end()) entry = it->second;
      }
      if (entry != nullptr)
      {
        protein_sequence[p] = generate ? generatedDecoy(searchedSequence(*entry), decoys) : searchedSequence(*entry);
        if (description[p].empty()) description[p] = entry->description;
      }
    }
    if (run != nullptr)
    {
      const auto& groups = run->getIndistinguishableProteins();
      for (Size g = 0; g < groups.size(); ++g)
      {
        for (const std::string& a : groups[g].accessions)
        {
          const auto it = protein_row.find(a);
          if (it != protein_row.end() && group[it->second] < 0) group[it->second] = static_cast<int32_t>(g);
        }
      }
    }

    //------------------------------------------------------------------------------
    // PSMs (the best hit of each identification with hits) and their peptides

    const PeptideIdentificationList& ids = merged.peptide_ids;
    const std::vector<PsmRow> rows = psmRows(ids);

    // The features, Percolator's inputs as ProSE's rescoring (ProSEAlgorithm::rescoreWithPercolator()) gives them:
    // - ln(hyperscore), the search engine score (Percolator's "score");
    // - the standard features of the .pin format (PercolatorInfile::getStandardFeatureSet()) that Percolator trains on:
    //   CalcMass, mass, peplen, charge<N>, enzN, enzC, enzInt, dm, absdm. PercolatorInfile computes them on copies of
    //   the best hits (below), with the enzyme of the search and charges N over its precursor charge range (widened
    //   to the charges of the PSMs; Percolator itself takes the range of the PSMs of the file it rescores);
    // - the search engine's features of the runs ('extra_features') and PeptDeep's.
    // Only those that some PSM has are columns.
    std::vector<std::string> feature_names{HYPERSCORE};
    std::vector<bool> pin_feature{false}; // whether column j is a .pin feature that PercolatorInfile computes
    const auto addFeature = [&feature_names, &pin_feature](const std::string& name, bool pin = false) {
      if (name.empty() || name == "score" || name == HYPERSCORE) return;
      if (std::find(feature_names.begin(), feature_names.end(), name) != feature_names.end()) return;
      feature_names.push_back(name);
      pin_feature.push_back(pin);
    };
    int min_pin_charge = min_charge, max_pin_charge = max_charge;
    for (const PsmRow& row : rows)
    {
      min_pin_charge = std::min(min_pin_charge, row.hit->getCharge());
      max_pin_charge = std::max(max_pin_charge, row.hit->getCharge());
    }
    // the bookkeeping columns of the .pin format, which Percolator does not train on (Percolator::rescorePSMs())
    static const std::set<std::string> pin_not_features{"SpecId", "Label", "ScanNr", "ExpMass", "Peptide", "Proteins"};
    for (const std::string& name : PercolatorInfile::getStandardFeatureSet(min_pin_charge, max_pin_charge))
    {
      if (!pin_not_features.contains(name)) addFeature(name, true);
    }
    for (const std::string& name : extra_features) addFeature(name);
    if (run != nullptr && run->getSearchParameters().metaValueExists("extra_features"))
    {
      for (const std::string& name : ListUtils::create<std::string>(run->getSearchParameters().getMetaValue("extra_features").toString()))
      {
        addFeature(name);
      }
    }
    for (const std::string& name : {Constants::UserParam::MS2_COSINE, Constants::UserParam::MS2_SPECTRAL_ANGLE,
                                    Constants::UserParam::MS2_PEARSON, Constants::UserParam::MS2_FRAC_PRED_FOUND,
                                    Constants::UserParam::RT_ABS_ERROR})
    {
      addFeature(name);
    }
    std::vector<UInt> feature_keys;
    for (const std::string& name : feature_names) feature_keys.push_back(metaIndex(name));

    // the sign of the isotope errors of the run: ProSE records 'isotope_error_convention' when it writes observed
    // minus theoretical; without it, the sign is the opposite (earlier versions, report:isotope_error_convention)
    const bool observed_minus_theoretical = run != nullptr
      && run->getSearchParameters().getMetaValue("isotope_error_convention", DataValue("")).toString() == "observed_minus_theoretical";

    const Size n_features = feature_names.size();
    std::vector<uint16_t> file;
    std::vector<std::string> native_id;
    std::vector<double> rt, mz, score, q_value, pep_value, precursor_error, delta_mass, features;
    std::vector<int8_t> charge, isotope_error;
    std::vector<uint32_t> peptide;
    std::vector<uint8_t> psm_decoy;
    for (auto* column : {&rt, &mz, &score, &q_value, &pep_value, &precursor_error, &delta_mass}) column->reserve(rows.size());
    features.reserve(rows.size() * n_features);
    std::vector<bool> feature_present(n_features, false);
    feature_present[0] = true; // ln(hyperscore): always a column
    PeptideColumns pt;
    std::unordered_map<std::string, uint32_t> peptide_row;
    std::vector<double> peptide_best; // the best score of the PSMs of each peptide

    for (const PsmRow& row : rows)
    {
      const PeptideIdentification& id = ids[row.id];
      const PeptideHit& hit = *row.hit;
      file.push_back(psmFile(id, keys));
      native_id.push_back(id.getSpectrumReference());
      rt.push_back(id.hasRT() ? id.getRT() : NaN);
      mz.push_back(id.hasMZ() ? id.getMZ() : NaN);
      charge.push_back(clamped<int8_t>(hit.getCharge()));
      const double search_score = hyperscore(id, hit, keys);
      const double percolator = percolatorScore(id, hit, keys);
      const double s = std::isnan(percolator) ? search_score : percolator;
      score.push_back(s);
      q_value.push_back(qValue(id, hit, keys));
      pep_value.push_back(pep(id, hit, keys));
      precursor_error.push_back(metaNumber(hit, keys.precursor_error));
      const double isotope = metaNumber(hit, keys.isotope_error);
      isotope_error.push_back(clamped<int8_t>(observed_minus_theoretical ? isotope : -isotope));
      delta_mass.push_back(metaNumber(hit, keys.delta_mass));
      psm_decoy.push_back(hit.isDecoy() ? 1 : 0);
      features.push_back(search_score);
      for (Size j = 1; j < n_features; ++j)
      {
        const double v = pin_feature[j] ? NaN : metaNumber(hit, feature_keys[j]); // .pin features: below
        feature_present[j] = feature_present[j] || !std::isnan(v);
        features.push_back(v);
      }

      // its peptide, unique by the OpenMS string
      const AASequence& sequence = hit.getSequence();
      std::string openms = sequence.toString();
      const auto [it, added] = peptide_row.emplace(openms, static_cast<uint32_t>(pt.openms.size()));
      peptide.push_back(it->second);
      if (!added)
      {
        peptide_best[it->second] = std::fmax(peptide_best[it->second], s);
        continue;
      }
      peptide_best.push_back(s);
      pt.openms.push_back(std::move(openms));
      pt.proforma.push_back(toProForma(sequence));
      pt.sequence.push_back(sequence.toUnmodifiedString());
      pt.decoy.push_back(hit.isDecoy() ? 1 : 0);
      // modifications: N-terminal at -1, the residues, C-terminal at the length of the sequence
      const std::string& residues = pt.sequence.back();
      const int length = static_cast<int>(residues.size());
      if (sequence.hasNTerminalModification()) pt.addModification(-1, sequence.getNTerminalModification());
      for (int i = 0; i < length; ++i)
      {
        if (sequence[i].isModified()) pt.addModification(i, sequence[i].getModification());
      }
      if (sequence.hasCTerminalModification()) pt.addModification(length, sequence.getCTerminalModification());
      pt.mod_offsets.push_back(static_cast<uint32_t>(pt.mod_position.size()));
      // proteins: one reference per peptide evidence of a protein of the run, at the position the search recorded
      for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
      {
        const auto row = protein_row.find(evidence.getProteinAccession());
        if (row == protein_row.end()) continue; // a protein the protein FDR removed
        const std::string& protein = protein_sequence[row->second];
        // the position the search recorded (0-based, end inclusive; UNKNOWN_POSITION = -1)
        Int from = evidence.getStart();
        Int to = evidence.getEnd();
        bool known = from >= 0 && to >= from;
        const bool in_sequence = known && static_cast<Size>(to) < protein.size() && to - from + 1 == length
                                 && protein.compare(static_cast<Size>(from), residues.size(), residues) == 0;
        if (!in_sequence)
        {
          // not where the sequence has the peptide (or not recorded): its first occurrence, if the sequence has it
          const auto at = protein.find(residues);
          if (at != std::string::npos)
          {
            from = static_cast<Int>(at);
            to = from + length - 1;
            known = true;
          }
        }
        pt.protein.push_back(row->second);
        pt.start.push_back(known ? static_cast<uint32_t>(from) : 0);
        pt.end.push_back(known ? static_cast<uint32_t>(to + 1) : 0); // exclusive; unknown: an empty range at 0
      }
      pt.protein_offsets.push_back(static_cast<uint32_t>(pt.protein.size()));
    }
    const Size n_psms = file.size();
    const Size n_peptides = pt.openms.size();

    // The .pin features, as PercolatorInfile computes them for Percolator, on copies of the best hits (without their
    // peak annotations), a chunk of them at a time; the results stay as they are. A hit PercolatorInfile skips (no
    // protein reference, or not known as target or decoy) has none.
    if (n_psms > 0 && std::find(pin_feature.begin(), pin_feature.end(), true) != pin_feature.end())
    {
      const std::string enzyme = ProSEAlgorithm::pinEnzymeName(run != nullptr ? run->getSearchParameters().digestion_enzyme.getName() : "");
      constexpr Size chunk_size = 4096;
      for (Size from = 0; from < n_psms; from += chunk_size)
      {
        const Size to = std::min(n_psms, from + chunk_size);
        PeptideIdentificationList chunk;
        chunk.reserve(to - from);
        for (Size i = from; i < to; ++i)
        {
          const PeptideIdentification& id = ids[rows[i].id];
          const PeptideHit& hit = *rows[i].hit;
          PeptideIdentification copy;
          static_cast<MetaInfoInterface&>(copy) = static_cast<const MetaInfoInterface&>(id);
          copy.setIdentifier(id.getIdentifier());
          copy.setScoreType(id.getScoreType());
          copy.setHigherScoreBetter(id.isHigherScoreBetter());
          copy.setRT(id.getRT());
          copy.setMZ(id.getMZ());
          PeptideHit hit_copy(hit.getScore(), hit.getRank(), hit.getCharge(), hit.getSequence());
          static_cast<MetaInfoInterface&>(hit_copy) = static_cast<const MetaInfoInterface&>(hit);
          hit_copy.setPeptideEvidences(hit.getPeptideEvidences());
          copy.insertHit(std::move(hit_copy));
          chunk.push_back(std::move(copy));
        }
        PercolatorInfile::PinFeatureMetaValueMap added;
        const auto skipped = PercolatorInfile::stampPinFeaturesOnHits(chunk, enzyme, min_pin_charge, max_pin_charge, added);
        for (Size i = from; i < to; ++i)
        {
          if (skipped.contains({i - from, 0})) continue;
          const PeptideHit& stamped = chunk[i - from].getHits()[0];
          for (Size j = 1; j < n_features; ++j)
          {
            if (!pin_feature[j]) continue;
            const double v = metaNumber(stamped, feature_keys[j]);
            feature_present[j] = feature_present[j] || !std::isnan(v);
            features[i * n_features + j] = v;
          }
        }
      }
    }

    // the score of a protein: the best score of its PSMs (aggregation "best" of the protein inference, on psms.score)
    for (Size k = 0; k < n_peptides; ++k)
    {
      for (uint32_t r = pt.protein_offsets[k]; r < pt.protein_offsets[k + 1]; ++r)
      {
        protein_score[pt.protein[r]] = std::fmax(protein_score[pt.protein[r]], peptide_best[k]);
      }
    }

    // only the features that some PSM has
    std::vector<Size> kept;
    for (Size j = 0; j < n_features; ++j)
    {
      if (feature_present[j]) kept.push_back(j);
    }
    std::vector<std::string> kept_names;
    for (const Size j : kept) kept_names.push_back(feature_names[j]);
    if (kept.size() != n_features)
    {
      std::vector<double> compact;
      compact.reserve(n_psms * kept.size());
      for (Size i = 0; i < n_psms; ++i)
      {
        for (const Size j : kept) compact.push_back(features[i * n_features + j]);
      }
      features.swap(compact);
    }

    //------------------------------------------------------------------------------
    // JavaScript

    std::vector<std::string> paths;
    for (const ProSEFileResult& f : files) paths.push_back(f.path);

    val psms = val::object();
    psms.set("count", static_cast<double>(n_psms));
    psms.set("file", toTypedArray(file));
    psms.set("nativeId", toStringArray(native_id));
    psms.set("rt", toTypedArray(rt));
    psms.set("mz", toTypedArray(mz));
    psms.set("charge", toTypedArray(charge));
    psms.set("peptide", toTypedArray(peptide));
    psms.set("score", toTypedArray(score));
    psms.set("qValue", toTypedArray(q_value));
    psms.set("pep", toTypedArray(pep_value));
    psms.set("precursorErrorPpm", toTypedArray(precursor_error));
    psms.set("isotopeError", toTypedArray(isotope_error));
    psms.set("deltaMass", toTypedArray(delta_mass));
    psms.set("decoy", toTypedArray(psm_decoy));
    psms.set("featureNames", toStringArray(kept_names));
    psms.set("features", toTypedArray(features));

    val peptides = val::object();
    peptides.set("count", static_cast<double>(n_peptides));
    peptides.set("openms", toStringArray(pt.openms));
    peptides.set("proforma", toStringArray(pt.proforma));
    peptides.set("sequence", toStringArray(pt.sequence));
    peptides.set("decoy", toTypedArray(pt.decoy));
    peptides.set("modOffsets", toTypedArray(pt.mod_offsets));
    peptides.set("modPosition", toTypedArray(pt.mod_position));
    peptides.set("modName", toStringArray(pt.mod_name));
    peptides.set("modUnimod", toTypedArray(pt.mod_unimod));
    peptides.set("modMass", toTypedArray(pt.mod_mass));
    peptides.set("proteinOffsets", toTypedArray(pt.protein_offsets));
    peptides.set("protein", toTypedArray(pt.protein));
    peptides.set("start", toTypedArray(pt.start));
    peptides.set("end", toTypedArray(pt.end));

    val proteins = val::object();
    proteins.set("count", static_cast<double>(n_proteins));
    proteins.set("accession", toStringArray(accession));
    proteins.set("description", toStringArray(description));
    proteins.set("sequence", toStringArray(protein_sequence));
    proteins.set("score", toTypedArray(protein_score));
    proteins.set("qValue", toTypedArray(protein_q));
    proteins.set("decoy", toTypedArray(protein_decoy));
    proteins.set("group", toTypedArray(group));

    val o = val::object();
    o.set("files", toStringArray(paths));
    o.set("psms", psms);
    o.set("peptides", peptides);
    o.set("proteins", proteins);
    return o;
  }

  //------------------------------------------------------------------------------
  // The peaks of the identified spectra

  namespace
  {
    /// Keeps the peaks of the MS2 spectra of one input file whose native IDs are wanted, as they are read
    /// (ProSEResultSet::collectSpectra())
    class IdentifiedSpectraConsumer : public Interfaces::IMSDataConsumer
    {
    public:
      /// @p wanted: the PSM rows of each native ID to keep, whose spectrum is added to @p stored; a native ID found is
      /// removed (a later spectrum with the same native ID is not kept)
      IdentifiedSpectraConsumer(std::unordered_map<std::string, std::vector<Size>>& wanted, ProSEStoredSpectra& stored) :
        wanted_(wanted),
        stored_(stored)
      {
      }

      void consumeSpectrum(MSSpectrum& s) override
      {
        if (s.getMSLevel() != 2) return;
        const auto it = wanted_.find(s.getNativeID());
        if (it == wanted_.end()) return;
        if (!s.isSorted()) s.sortByPosition();
        const auto index = static_cast<int32_t>(stored_.precursor_mz.size());
        for (const Peak1D& peak : s)
        {
          stored_.mz.push_back(static_cast<float>(peak.getMZ()));
          stored_.intensity.push_back(peak.getIntensity());
        }
        stored_.offsets.push_back(static_cast<uint32_t>(stored_.mz.size));
        const std::vector<Precursor>& precursors = s.getPrecursors();
        stored_.precursor_mz.push_back(precursors.empty() ? NaN : precursors[0].getMZ());
        // OpenMS has 0 for an energy that is not known
        const double energy = precursors.empty() ? 0.0 : precursors[0].getActivationEnergy();
        stored_.collision_energy.push_back(energy > 0.0 ? static_cast<float>(energy) : std::numeric_limits<float>::quiet_NaN());
        for (const Size row : it->second) stored_.psm_spectrum[row] = index;
        wanted_.erase(it);
      }

      void consumeChromatogram(MSChromatogram& /*c*/) override {}
      void setExpectedSize(size_t /*spectra*/, size_t /*chromatograms*/) override {}
      void setExperimentalSettings(const ExperimentalSettings& /*settings*/) override {}

    private:
      std::unordered_map<std::string, std::vector<Size>>& wanted_;
      ProSEStoredSpectra& stored_;
    };

    /// A JavaScript typed array holding a copy of the values of @p v
    template <typename T>
    val toTypedArray(const BlockVector<T>& v)
    {
      val result = val::global(TypedArrayName<T>::value).new_(v.size);
      std::size_t at = 0;
      for (const std::vector<T>& block : v.blocks)
      {
        if (block.empty()) continue;
        result.call<void>("set", val(emscripten::typed_memory_view(block.size(), block.data())), static_cast<double>(at));
        at += block.size();
      }
      return result;
    }
  } // namespace

  void ProSEResultSet::collectSpectra()
  {
    stored_spectra = ProSEStoredSpectra{};
    ProSEStoredSpectra& stored = stored_spectra;
    const ScoreKeys keys;
    const PeptideIdentificationList& ids = merged.peptide_ids;
    const std::vector<PsmRow> rows = psmRows(ids);
    stored.psm_spectrum.assign(rows.size(), -1);
    if (!keep_spectra) return;

    // The rows to keep the spectra of, by file and native ID: those whose q-value passes FDR:PSM, decoys too; all rows
    // without q-values, and with FDR:PSM 0 (no filter)
    std::vector<std::unordered_map<std::string, std::vector<Size>>> wanted(files.size());
    for (Size r = 0; r < rows.size(); ++r)
    {
      const PeptideIdentification& id = ids[rows[r].id];
      const double q = qValue(id, *rows[r].hit, keys);
      if (!std::isnan(q) && psm_fdr > 0.0 && q > psm_fdr) continue;
      const Size f = psmFile(id, keys);
      if (f < files.size()) wanted[f][id.getSpectrumReference()].push_back(r);
    }

    // One file at a time, streamed: only the spectra kept stay in memory, and only they are sorted (by the consumer).
    // (Not setSkipChromatograms(): MzMLHandler then skips the contents of the spectra, too.)
    PeakFileOptions options;
    options.setMSLevels({2});
    options.setSortSpectraByMZ(false);
    for (Size f = 0; f < files.size(); ++f)
    {
      if (wanted[f].empty()) continue;
      MzMLFile mzml;
      mzml.setOptions(options);
      IdentifiedSpectraConsumer consumer(wanted[f], stored);
      try
      {
        mzml.transform(files[f].path, &consumer, true, true); // one pass: no counts, no metadata
      }
      catch (const Exception::BaseException& e) // the file was searched, but cannot be read again: its PSMs have none
      {
        OPENMS_LOG_WARN << "ProSE: cannot read the spectra of the PSMs of '" << files[f].path << "' again (" << e.what()
                        << "); the PSMs whose spectra were not read have none." << std::endl;
      }
      wanted[f].clear();
    }
    stored.mz.shrink_to_fit();
    stored.intensity.shrink_to_fit();
  }

  JS::StoredSpectra ProSEResultSet::spectra() const
  {
    const ProSEStoredSpectra& s = stored_spectra;
    val o = val::object();
    o.set("psmSpectrum", toTypedArray(s.psm_spectrum));
    o.set("offsets", toTypedArray(s.offsets));
    o.set("mz", toTypedArray(s.mz));
    o.set("intensity", toTypedArray(s.intensity));
    o.set("precursorMz", toTypedArray(s.precursor_mz));
    o.set("collisionEnergy", toTypedArray(s.collision_energy));
    return o;
  }

  void ProSEResultSet::debugSetNativeId(int row, const std::string& native_id)
  {
    const std::vector<PsmRow> rows = psmRows(merged.peptide_ids);
    if (row < 0) throw Exception::IndexUnderflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, row, rows.size());
    if (static_cast<Size>(row) >= rows.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, row, rows.size());
    merged.peptide_ids[rows[row].id].setSpectrumReference(native_id);
    collectSpectra();
  }
} // namespace OpenMS::WASM
