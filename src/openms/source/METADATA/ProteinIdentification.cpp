// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>

#include <numeric>
#include <unordered_set>

using namespace std;

namespace OpenMS
{

  const std::string ProteinIdentification::NamesOfPeakMassType[] = {"Monoisotopic", "Average"};

  ProteinIdentification::ProteinGroup::ProteinGroup() :
    probability(0.0), accessions()
  {
  }

  bool ProteinIdentification::ProteinGroup::operator==(const ProteinGroup& rhs) const
  {
    return std::tie(probability, accessions) == std::tie(rhs.probability, rhs.accessions);
  }

  bool ProteinIdentification::ProteinGroup::operator<(const ProteinGroup& rhs) const
  {
    // comparison of probabilities is intentionally "the wrong way around":
    if (probability > rhs.probability)
    {
      return true;
    }
    if (probability < rhs.probability)
    {
      return false;
    }
    if (accessions.size() < rhs.accessions.size())
    {
      return true;
    }
    if (accessions.size() > rhs.accessions.size())
    {
      return false;
    }
    return accessions < rhs.accessions;
  }

  const ProteinIdentification::ProteinGroup::FloatDataArrays& ProteinIdentification::ProteinGroup::getFloatDataArrays() const
  {
    return float_data_arrays_;
  }

  void ProteinIdentification::ProteinGroup::setFloatDataArrays(const ProteinIdentification::ProteinGroup::FloatDataArrays &fda)
  {
    float_data_arrays_ = fda;
  }

  const ProteinIdentification::ProteinGroup::StringDataArrays& ProteinIdentification::ProteinGroup::getStringDataArrays() const
  {
    return string_data_arrays_;
  }

  void ProteinIdentification::ProteinGroup::setStringDataArrays(const ProteinIdentification::ProteinGroup::StringDataArrays& sda)
  {
    string_data_arrays_ = sda;
  }

  ProteinIdentification::ProteinGroup::StringDataArrays& ProteinIdentification::ProteinGroup::getStringDataArrays()
  {
    return string_data_arrays_;
  }

  const ProteinIdentification::ProteinGroup::IntegerDataArrays& ProteinIdentification::ProteinGroup::getIntegerDataArrays() const
  {
    return integer_data_arrays_;
  }

  ProteinIdentification::ProteinGroup::IntegerDataArrays& ProteinIdentification::ProteinGroup::getIntegerDataArrays()
  {
    return integer_data_arrays_;
  }

  void ProteinIdentification::ProteinGroup::setIntegerDataArrays(const ProteinIdentification::ProteinGroup::IntegerDataArrays& ida)
  {
    integer_data_arrays_ = ida;
  }

  ProteinIdentification::SearchParameters::SearchParameters() :
      db(),
      db_version(),
      taxonomy(),
      charges(),
      mass_type(MONOISOTOPIC),
      fixed_modifications(),
      variable_modifications(),
      missed_cleavages(0),
      fragment_mass_tolerance(0.0),
      fragment_mass_tolerance_ppm(false),
      precursor_mass_tolerance(0.0),
      precursor_mass_tolerance_ppm(false),
      digestion_enzyme("unknown_enzyme", ""),
      enzyme_term_specificity(EnzymaticDigestion::SPEC_UNKNOWN)
  {
  }

  bool ProteinIdentification::SearchParameters::operator==(const SearchParameters& rhs) const
  {
    return
        std::tie(db, db_version, taxonomy, charges, mass_type, fixed_modifications, variable_modifications,
            missed_cleavages, fragment_mass_tolerance, fragment_mass_tolerance_ppm, precursor_mass_tolerance,
            precursor_mass_tolerance_ppm, digestion_enzyme, enzyme_term_specificity) ==
        std::tie(rhs.db, rhs.db_version, rhs.taxonomy, rhs.charges, rhs.mass_type, rhs.fixed_modifications,
            rhs.variable_modifications, rhs.missed_cleavages, rhs.fragment_mass_tolerance,
            rhs.fragment_mass_tolerance_ppm, rhs.precursor_mass_tolerance,
            rhs.precursor_mass_tolerance_ppm, rhs.digestion_enzyme, rhs.enzyme_term_specificity);
  }

  bool ProteinIdentification::SearchParameters::operator!=(const SearchParameters& rhs) const
  {
    return !(*this == rhs);
  }

  bool ProteinIdentification::SearchParameters::mergeable(const ProteinIdentification::SearchParameters& sp, const String& experiment_type) const
  {
    String spdb = sp.db;
    spdb.substitute("\\","/");
    String pdb = this->db;
    pdb.substitute("\\","/");

    if  (this->precursor_mass_tolerance != sp.precursor_mass_tolerance ||
        this->precursor_mass_tolerance_ppm != sp.precursor_mass_tolerance_ppm ||
        File::basename(pdb) != File::basename(spdb) ||
        this->db_version != sp.db_version ||
        this->fragment_mass_tolerance != sp.fragment_mass_tolerance ||
        this->fragment_mass_tolerance_ppm != sp.fragment_mass_tolerance_ppm ||
        this->charges != sp.charges ||
        this->digestion_enzyme != sp.digestion_enzyme ||
        this->taxonomy != sp.taxonomy ||
         this->enzyme_term_specificity != sp.enzyme_term_specificity)
    {
      return false;
    }

    set<String> fixed_mods(this->fixed_modifications.begin(), this->fixed_modifications.end());
    set<String> var_mods(this->variable_modifications.begin(), this->variable_modifications.end());
    set<String> curr_fixed_mods(sp.fixed_modifications.begin(), sp.fixed_modifications.end());
    set<String> curr_var_mods(sp.variable_modifications.begin(), sp.variable_modifications.end());
    if (fixed_mods != curr_fixed_mods ||
        var_mods != curr_var_mods)
    {
      if (experiment_type != "labeled_MS1")
      {
        return false;
      }
      else
      {
        //TODO actually introduce a flag for labelling modifications in the Mod datastructures?
        //OR put a unique ID for the used mod as a UserParam to the mapList entries (consensusHeaders)
        //TODO actually you would probably need an experimental design here, because
        //settings have to agree exactly in a FractionGroup but can slightly differ across runs.
        //Or just ignore labelling mods during the check
        return true;
      }
    }
    return true;
  }

  int ProteinIdentification::SearchParameters::getChargeValue_(String& charge_str) const
  {
    // We have to do this because some people/tools put the + or - AFTER the number...
    bool neg = charge_str.hasSubstring('-');
    neg ? charge_str.remove('-') : charge_str.remove('+');
    int val = charge_str.toInt();
    return neg ? -val : val;
  }

  std::pair<int,int> ProteinIdentification::SearchParameters::getChargeRange() const
  {
    std::pair<int,int> result{0,0};

    try // is there only one number (min = max)?
    {
      result.first = charges.toInt();
      result.second = result.first;
    }
    catch (Exception::ConversionError&) // nope, something else
    {
      if (charges.hasSubstring(',')) // it's probably a list
      {
        IntList chgs = ListUtils::create<Int>(charges);
        auto minmax = minmax_element(chgs.begin(), chgs.end());
        result.first = *minmax.first;
        result.second = *minmax.second;
      }
      else if (charges.hasSubstring(':')) // it's probably a range
      {
        StringList chgs;
        charges.split(':', chgs);
        if (chgs.size() > 2)
        {
          throw OpenMS::Exception::MissingInformation(
            __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Charge string in SearchParameters not parseable.");
        }
        result.first = getChargeValue_(chgs[0]);
        result.second = getChargeValue_(chgs[1]);
      }
      else
      {
        size_t pos = charges.find('-', 0);
        std::vector<size_t> minus_positions;
        while (pos != string::npos)
        {
          minus_positions.push_back(pos);
          pos = charges.find('-', pos + 1);
        }
        if (!minus_positions.empty() && minus_positions.size() <= 3) // it's probably a range with '-'
        {
          Size split_pos(0);
          if (minus_positions.size() <= 1) // split at first minus
          {
            split_pos = minus_positions[0];
          }
          else
          {
            split_pos = minus_positions[1];
          }
          String first = charges.substr(0, split_pos);
          String second = charges.substr(split_pos + 1, string::npos);
          result.first = getChargeValue_(first);
          result.second = getChargeValue_(second);
        }
      }
    }
    return result;
  }

  ProteinIdentification::ProteinIdentification() :
    MetaInfoInterface(),
    id_(),
    search_engine_(),
    search_engine_version_(),
    search_parameters_(),
    date_(),
    protein_score_type_(),
    higher_score_better_(true),
    protein_hits_(),
    protein_groups_(),
    indistinguishable_proteins_(),
    protein_significance_threshold_(0.0)
  {
  }

  ProteinIdentification::~ProteinIdentification() = default;

  void ProteinIdentification::setDateTime(const DateTime& date)
  {
    date_ = date;
  }

  const DateTime& ProteinIdentification::getDateTime() const
  {
    return date_;
  }

  const vector<ProteinHit>& ProteinIdentification::getHits() const
  {
    return protein_hits_;
  }

  vector<ProteinHit>& ProteinIdentification::getHits()
  {
    return protein_hits_;
  }

  void ProteinIdentification::setHits(const vector<ProteinHit>& protein_hits)
  {
    protein_hits_ = protein_hits;
  }

  vector<ProteinHit>::iterator ProteinIdentification::findHit(
    const String& accession)
  {
    vector<ProteinHit>::iterator pos = protein_hits_.begin();
    for (; pos != protein_hits_.end(); ++pos)
    {
      if (pos->getAccession() == accession)
      {
        break;
      }
    }
    return pos;
  }

  const vector<ProteinIdentification::ProteinGroup>& ProteinIdentification::getProteinGroups() const
  {
    return protein_groups_;
  }

  vector<ProteinIdentification::ProteinGroup>& ProteinIdentification::getProteinGroups()
  {
    return protein_groups_;
  }

  void ProteinIdentification::insertProteinGroup(const ProteinIdentification::ProteinGroup& group)
  {
    protein_groups_.push_back(group);
  }

  const vector<ProteinIdentification::ProteinGroup>&
  ProteinIdentification::getIndistinguishableProteins() const
  {
    return indistinguishable_proteins_;
  }

  vector<ProteinIdentification::ProteinGroup>&
  ProteinIdentification::getIndistinguishableProteins()
  {
    return indistinguishable_proteins_;
  }

  void ProteinIdentification::insertIndistinguishableProteins(
    const ProteinIdentification::ProteinGroup& group)
  {
    indistinguishable_proteins_.push_back(group);
  }

  void ProteinIdentification::fillIndistinguishableGroupsWithSingletons()
  {
    unordered_set<string> groupedAccessions;
    for (const ProteinGroup& proteinGroup : indistinguishable_proteins_)
    {
      for (const String& acc : proteinGroup.accessions)
      {
        groupedAccessions.insert(acc);
      }
    }

    for (const ProteinHit& protein : getHits())
    {
      const String& acc = protein.getAccession();
      if (groupedAccessions.find(acc) == groupedAccessions.end())
      {
        groupedAccessions.insert(acc);
        ProteinGroup pg;
        pg.accessions.push_back(acc);
        pg.probability = protein.getScore();
        indistinguishable_proteins_.push_back(pg);
      }
    }
  }

  // retrieval of the peptide significance threshold value
  double ProteinIdentification::getSignificanceThreshold() const
  {
    return protein_significance_threshold_;
  }

  // setting of the peptide significance threshold value
  void ProteinIdentification::setSignificanceThreshold(double value)
  {
    protein_significance_threshold_ = value;
  }

  void ProteinIdentification::setScoreType(const String& type)
  {
    protein_score_type_ = type;
  }

  const String& ProteinIdentification::getScoreType() const
  {
    return protein_score_type_;
  }

  void ProteinIdentification::insertHit(const ProteinHit& protein_hit)
  {
    protein_hits_.push_back(protein_hit);
  }

  void ProteinIdentification::insertHit(ProteinHit&& protein_hit)
  {
    protein_hits_.push_back(std::forward<ProteinHit>(protein_hit));
  }

  void ProteinIdentification::setPrimaryMSRunPath(const StringList& s, bool raw)
  {
    String meta_name = raw ? "spectra_data_raw" : "spectra_data";
    setMetaValue(meta_name, DataValue(StringList()));
    if (s.empty())
    {
      OPENMS_LOG_WARN << "Setting an empty value for primary MS runs paths." << std::endl;
    }
    else
    {
      addPrimaryMSRunPath(s, raw);
    }
  }

  void ProteinIdentification::setPrimaryMSRunPath(const StringList& s, MSExperiment& e)
  {
    StringList ms_path;
    e.getPrimaryMSRunPath(ms_path);
    if (ms_path.size() == 1)
    {
      FileTypes::Type filetype = FileHandler::getTypeByFileName(ms_path[0]);
      if ((filetype == FileTypes::MZML) && File::exists(ms_path[0]))
      {
        setMetaValue("spectra_data", DataValue(StringList({ms_path[0]})));
        return; // don't do anything else in this case
      }
      if (filetype == FileTypes::RAW)
      {
        setMetaValue("spectra_data_raw", DataValue(StringList({ms_path[0]})));
      }
    }
    setPrimaryMSRunPath(s);
  }

  /// get the file path to the first MS runs
  void ProteinIdentification::getPrimaryMSRunPath(StringList& output, bool raw) const
  {
    String meta_name = raw ? "spectra_data_raw" : "spectra_data";
    if (metaValueExists(meta_name))
    {
      output = getMetaValue(meta_name);
    }
  }

  void ProteinIdentification::addPrimaryMSRunPath(const StringList& s, bool raw)
  {
    String meta_name = raw ? "spectra_data_raw" : "spectra_data";
    if (!raw) // mzML files expected
    {
      for (const String &filename : s)
      {
        FileTypes::Type filetype = FileHandler::getTypeByFileName(filename);
        if (filetype != FileTypes::MZML)
        {
          OPENMS_LOG_WARN << "To ensure tracability of results please prefer mzML files as primary MS runs.\n"
                          << "Filename: '" << filename << "'" << std::endl;
        }
      }
    }
    StringList spectra_data = getMetaValue(meta_name, DataValue(StringList()));
    spectra_data.insert(spectra_data.end(), s.begin(), s.end());
    setMetaValue(meta_name, spectra_data);
  }

  void ProteinIdentification::addPrimaryMSRunPath(const String& s, bool raw)
  {
    addPrimaryMSRunPath(StringList({s}), raw);
  }

  Size ProteinIdentification::nrPrimaryMSRunPaths(bool raw) const
  {
    String meta_name = raw ? "spectra_data_raw" : "spectra_data";
    StringList spectra_data = getMetaValue(meta_name, DataValue(StringList()));
    return spectra_data.size();
  }

  //TODO find a more robust way to figure that out. CV Terms?
  bool ProteinIdentification::hasInferenceData() const
  {
    return !getInferenceEngine().empty();
  }

  bool ProteinIdentification::hasInferenceEngineAsSearchEngine() const
  {
    String se = getSearchEngine();
    return
        se == "Fido" || // for downwards compatibility: FidoAdapter overwrites when it merges several runs
        se == "BayesianProteinInference" || // for backwards compatibility
        se == "Epifany" ||
        (se == "Percolator" && !indistinguishable_proteins_.empty()) || // be careful, Percolator could be run with or without protein inference
        se == "ProteinInference";
  }

  bool ProteinIdentification::peptideIDsMergeable(const ProteinIdentification& id_run, const String& experiment_type) const
  {
    const String& warn = " You probably do not want to merge the results with this tool."
                         " For merging searches with different engines/settings please use ConsensusID or PercolatorAdapter"
                         " to create a comparable score.";
    const String& engine = this->getSearchEngine();
    const String& version = this->getSearchEngineVersion();

    bool ok = true;

    if (id_run.getSearchEngine() != engine || id_run.getSearchEngineVersion() != version)
    {
      ok = false;
      OPENMS_LOG_WARN << "Search engine " + id_run.getSearchEngine() + "from IDRun " + id_run.getIdentifier()
                         + " does not match with the others." + warn;
    }
    const ProteinIdentification::SearchParameters& params = this->getSearchParameters();
    const ProteinIdentification::SearchParameters& sp = id_run.getSearchParameters();
    if (!params.mergeable(sp, experiment_type))
    {
      ok = false;
      OPENMS_LOG_WARN << "Searchengine settings or modifications from IDRun " + id_run.getIdentifier() + " do not match with the others." + warn;
    }
    // TODO else merge as far as possible (mainly mods I guess)
    return ok;
  }

  vector<pair<String,String>> ProteinIdentification::getSearchEngineSettingsAsPairs(const String& se) const
  {
    vector<pair<String,String>> result;
    const auto& params = this->getSearchParameters();
    if (se.empty() || (this->getSearchEngine() == se
                        && this->getSearchEngine() != "Percolator" //meaningless settings
                        && !this->getSearchEngine().hasPrefix("ConsensusID"))) //meaningless settings
    {
      //TODO add spectra_data?
      result.emplace_back("db", params.db);
      result.emplace_back("db_version", params.db_version);
      result.emplace_back("fragment_mass_tolerance", params.fragment_mass_tolerance);
      result.emplace_back("fragment_mass_tolerance_unit", params.fragment_mass_tolerance_ppm ? "ppm" : "Da");
      result.emplace_back("precursor_mass_tolerance", params.precursor_mass_tolerance);
      result.emplace_back("precursor_mass_tolerance_unit", params.precursor_mass_tolerance_ppm ? "ppm" : "Da");
      result.emplace_back("enzyme", params.digestion_enzyme.getName());
      result.emplace_back("enzyme_term_specificity", EnzymaticDigestion::NamesOfSpecificity[params.enzyme_term_specificity]);
      result.emplace_back("charges", params.charges);
      result.emplace_back("missed_cleavages", params.missed_cleavages);
      result.emplace_back("fixed_modifications", ListUtils::concatenate(params.fixed_modifications,","));
      result.emplace_back("variable_modifications", ListUtils::concatenate(params.variable_modifications,","));
    }
    else
    {
      vector<String> mvkeys;
      params.getKeys(mvkeys);
      for (const String & mvkey : mvkeys)
      {
        if (mvkey.hasPrefix(se))
        {
          result.emplace_back(mvkey.substr(se.size()+1), params.getMetaValue(mvkey));
        }
      }
    }
    return result;
  }


  // Equality operator
  bool ProteinIdentification::operator==(const ProteinIdentification& rhs) const
  {
    return MetaInfoInterface::operator==(rhs) &&
            std::tie(id_, search_engine_, search_engine_version_,
                search_parameters_, date_, protein_hits_, protein_groups_,
                indistinguishable_proteins_, protein_score_type_,
                protein_significance_threshold_, higher_score_better_) ==
            std::tie(rhs.id_, rhs.search_engine_, rhs.search_engine_version_,
                     rhs.search_parameters_, rhs.date_, rhs.protein_hits_, rhs.protein_groups_,
                     rhs.indistinguishable_proteins_, rhs.protein_score_type_,
                     rhs.protein_significance_threshold_, rhs.higher_score_better_);
  }

  // Inequality operator
  bool ProteinIdentification::operator!=(const ProteinIdentification& rhs) const
  {
    return !operator==(rhs);
  }

  void ProteinIdentification::sort()
  {
    if (higher_score_better_)
    {
      std::stable_sort(protein_hits_.begin(), protein_hits_.end(), ProteinHit::ScoreMore());
    }
    else
    {
      std::stable_sort(protein_hits_.begin(), protein_hits_.end(), ProteinHit::ScoreLess());
    }
  }

  void ProteinIdentification::assignRanks()
  {
    if (protein_hits_.empty())
    {
      return;
    }
    UInt rank = 1;
    sort();
    vector<ProteinHit>::iterator lit = protein_hits_.begin();
    double tmpscore = lit->getScore();
    while (lit != protein_hits_.end())
    {
      lit->setRank(rank);
      ++lit;
      if (lit != protein_hits_.end() && lit->getScore() != tmpscore)
      {
        ++rank;
        tmpscore = lit->getScore();
      }
    }
  }


  void ProteinIdentification::fillEvidenceMapping_(unordered_map<String, set<PeptideEvidence> >& map_acc_2_evidence,
                                                   const std::vector<PeptideIdentification>& pep_ids) const
  {
    //TODO check matching identifiers?
    for (const auto & peptide_id : pep_ids)
    {
      // peptide hits
      const vector<PeptideHit>& peptide_hits = peptide_id.getHits();
      for (const auto & peptide_hit : peptide_hits)
      {
        const std::vector<PeptideEvidence>& ph_evidences = peptide_hit.getPeptideEvidences();
        // matched proteins for hit
        for (const auto & evidence : ph_evidences)
        {
          map_acc_2_evidence[evidence.getProteinAccession()].insert(evidence);
        }
      }
    }
  }

  void ProteinIdentification::computeCoverage(const std::vector<PeptideIdentification>& pep_ids)
  {
    // map protein accession to the corresponding peptide evidence
    unordered_map<String, set<PeptideEvidence> > map_acc_2_evidence;
    fillEvidenceMapping_(map_acc_2_evidence, pep_ids);
    computeCoverageFromEvidenceMapping_(map_acc_2_evidence);
  }

  void ProteinIdentification::computeCoverage(const ConsensusMap& cmap, bool use_unassigned_ids)
  {
    // map protein accession to the corresponding peptide evidence
    unordered_map<String, set<PeptideEvidence> > map_acc_2_evidence;
    for (const auto& feat : cmap)
    {
      fillEvidenceMapping_(map_acc_2_evidence,feat.getPeptideIdentifications());
    }
    if (use_unassigned_ids)
    {
      fillEvidenceMapping_(map_acc_2_evidence, cmap.getUnassignedPeptideIdentifications());
    }
    computeCoverageFromEvidenceMapping_(map_acc_2_evidence);
  }

  void ProteinIdentification::computeCoverageFromEvidenceMapping_(const unordered_map<String, set<PeptideEvidence>>& map_acc_2_evidence)
  {
    for (Size i = 0; i < this->protein_hits_.size(); ++i)
    {
      const Size protein_length = this->protein_hits_[i].getSequence().length();
      if (protein_length == 0)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, " ProteinHits do not contain a protein sequence. Cannot compute coverage! Use PeptideIndexer to annotate proteins with sequence information.");
      }
      vector<bool> covered_amino_acids(protein_length, false);

      const String& accession = this->protein_hits_[i].getAccession();
      double coverage = 0.0;
      if (map_acc_2_evidence.find(accession) != map_acc_2_evidence.end())
      {
        const set<PeptideEvidence>& evidences = map_acc_2_evidence.find(accession)->second;
        for (const auto & evidence : evidences)
        {
          int start = evidence.getStart();
          int stop = evidence.getEnd();

          if (start == PeptideEvidence::UNKNOWN_POSITION || stop == PeptideEvidence::UNKNOWN_POSITION)
          {
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              " PeptideEvidence does not contain start or end position. Cannot compute coverage!");
          }

          if (start < 0 || stop < start || stop > static_cast<int>(protein_length))
          {
            const String message = " PeptideEvidence (start/end) (" + String(start) + "/" + String(stop) +
                                   " ) are invalid or point outside of protein '" + accession +
                                   "' (length: " + String(protein_length) +
                                   "). Cannot compute coverage!";
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
          }

          std::fill(covered_amino_acids.begin() + start, covered_amino_acids.begin() + stop + 1, true);
        }
        coverage = 100.0 * (double) std::accumulate(covered_amino_acids.begin(),
                                                   covered_amino_acids.end(), 0) / (double) protein_length;
      }
      this->protein_hits_[i].setCoverage(coverage);
    }
  }

  void ProteinIdentification::computeModifications(
    const std::vector<PeptideIdentification>& pep_ids,
    const StringList& skip_modifications)
  {
    // map protein accession to observed position,modifications pairs
    unordered_map<String, set<pair<Size, ResidueModification>>> prot2mod;

    fillModMapping_(pep_ids, skip_modifications, prot2mod);

    for (auto & protein_hit : protein_hits_)
    {
      const String& accession = protein_hit.getAccession();
      if (prot2mod.find(accession) != prot2mod.end())
      {
        protein_hit.setModifications(prot2mod[accession]);
      }
    }
  }

  void ProteinIdentification::computeModifications(
    const ConsensusMap& cmap,
    const StringList& skip_modifications,
    bool use_unassigned_ids)
  {
    // map protein accession to observed position,modifications pairs
    unordered_map<String, set<pair<Size, ResidueModification>>> prot2mod;

    for (const auto& feat : cmap)
    {
      fillModMapping_(feat.getPeptideIdentifications(), skip_modifications, prot2mod);
    }
    if (use_unassigned_ids) fillModMapping_(cmap.getUnassignedPeptideIdentifications(), skip_modifications, prot2mod);

    for (auto & protein_hit : protein_hits_)
    {
      const String& accession = protein_hit.getAccession();
      if (prot2mod.find(accession) != prot2mod.end())
      {
        protein_hit.setModifications(prot2mod[accession]);
      }
    }
  }

  void ProteinIdentification::fillModMapping_(const vector<PeptideIdentification>& pep_ids, const StringList& skip_modifications,
                                              unordered_map<String, set<pair<Size, ResidueModification>>>& prot2mod) const
  {
    for (const auto & peptide_id : pep_ids)
    {
      // peptide hits
      const vector<PeptideHit>& peptide_hits = peptide_id.getHits();
      for (const auto & peptide_hit : peptide_hits)
      {
        const AASequence& aas = peptide_hit.getSequence();
        const vector<PeptideEvidence>& ph_evidences = peptide_hit.getPeptideEvidences();

        // skip unmodified peptides
        if (aas.isModified())
        {
          if (aas.hasNTerminalModification())
          {
            const ResidueModification * res_mod = aas.getNTerminalModification();
            // skip mod if Id, e.g. 'Carbamidomethyl' or full id e.g., 'Carbamidomethyl (C)' match.
            if (std::find(skip_modifications.begin(), skip_modifications.end(), res_mod->getId()) == skip_modifications.end()
             && std::find(skip_modifications.begin(), skip_modifications.end(), res_mod->getFullId()) == skip_modifications.end())
            {
              for (const auto & ph_evidence : ph_evidences)
              {
                const String& acc = ph_evidence.getProteinAccession();
                const Size mod_pos = ph_evidence.getStart(); // mod at N terminus
                prot2mod[acc].insert(make_pair(mod_pos, *res_mod));
              }
            }
          }

          for (Size ai = 0; ai != aas.size(); ++ai)
          {
            if (aas[ai].isModified())
            {
              const ResidueModification * res_mod = aas[ai].getModification();

              if (find(skip_modifications.begin(), skip_modifications.end(), res_mod->getId()) == skip_modifications.end()
               && find(skip_modifications.begin(), skip_modifications.end(), res_mod->getFullId()) == skip_modifications.end())
              {
                for (const auto & ph_evidence : ph_evidences)
                {
                  const String& acc = ph_evidence.getProteinAccession();
                  const Size mod_pos = ph_evidence.getStart() + ai; // start + ai
                  prot2mod[acc].insert(make_pair(mod_pos, *res_mod));
                }
              }
            }
          }

          if (aas.hasCTerminalModification())
          {
            const ResidueModification * res_mod = aas.getCTerminalModification();
            // skip mod?
            if (find(skip_modifications.begin(), skip_modifications.end(), res_mod->getId()) == skip_modifications.end()
             && find(skip_modifications.begin(), skip_modifications.end(), res_mod->getFullId()) == skip_modifications.end())
            {
              for (const auto & ph_evidence : ph_evidences)
              {
                const String& acc = ph_evidence.getProteinAccession();
                const Size mod_pos = ph_evidence.getEnd(); // mod at C terminus
                prot2mod[acc].insert(make_pair(mod_pos, *res_mod));
              }
            }
          }
        }
      }
    }
  }

  bool ProteinIdentification::isHigherScoreBetter() const
  {
    return higher_score_better_;
  }

  void ProteinIdentification::setHigherScoreBetter(bool value)
  {
    higher_score_better_ = value;
  }

  const String& ProteinIdentification::getIdentifier() const
  {
    return id_;
  }

  void ProteinIdentification::setIdentifier(const String& id)
  {
    id_ = id;
  }

  void ProteinIdentification::setSearchEngine(const String& search_engine)
  {
    search_engine_ = search_engine;
  }

  const String& ProteinIdentification::getSearchEngine() const
  {
    return search_engine_;
  }

  const String ProteinIdentification::getOriginalSearchEngineName() const
  {
    // TODO: extend to multiple search engines and merging
    String engine = search_engine_;
    if (!engine.hasSubstring("Percolator") && !engine.hasSubstring("ConsensusID"))
    {
      return engine;
    }

    String original_SE = "Unknown";
    vector<String> mvkeys;
    getSearchParameters().getKeys(mvkeys);
    for (const String& mvkey : mvkeys)
    {
      if (mvkey.hasPrefix("SE:") && !mvkey.hasSubstring("percolator"))
      {
        original_SE = mvkey.substr(3);
        break; // multiSE percolator before consensusID not allowed; we take first only
      }
    }
    return original_SE;
  }

  void ProteinIdentification::setSearchEngineVersion(const String& search_engine_version)
  {
    search_engine_version_ = search_engine_version;
  }

  const String& ProteinIdentification::getSearchEngineVersion() const
  {
    return search_engine_version_;
  }

  void ProteinIdentification::setSearchParameters(const SearchParameters& search_parameters)
  {
    search_parameters_ = search_parameters;
  }

  void ProteinIdentification::setSearchParameters(SearchParameters&& search_parameters)
  {
    search_parameters_ = std::move(search_parameters);
  }

  const ProteinIdentification::SearchParameters& ProteinIdentification::getSearchParameters() const
  {
    return search_parameters_;
  }

  ProteinIdentification::SearchParameters& ProteinIdentification::getSearchParameters()
  {
    return search_parameters_;
  }

  void ProteinIdentification::setInferenceEngine(const String& inference_engine)
  {
    this->search_parameters_.setMetaValue("InferenceEngine", inference_engine);
  }

  const String ProteinIdentification::getInferenceEngine() const
  {
    if (this->search_parameters_.metaValueExists("InferenceEngine"))
    {
      return this->search_parameters_.getMetaValue("InferenceEngine");
    }
    else if (hasInferenceEngineAsSearchEngine())
    {
      return search_engine_;
    }
    return "";
  }

  void ProteinIdentification::setInferenceEngineVersion(const String& search_engine_version)
  {
    this->search_parameters_.setMetaValue("InferenceEngineVersion", search_engine_version);
  }

  const String ProteinIdentification::getInferenceEngineVersion() const
  {
    if (this->search_parameters_.metaValueExists("InferenceEngineVersion"))
    {
      return this->search_parameters_.getMetaValue("InferenceEngineVersion");
    }
    else if (hasInferenceData())
    {
      return search_engine_version_;
    }
    return "";
  }

  void ProteinIdentification::copyMetaDataOnly(const ProteinIdentification& p)
  {
    id_ = p.id_;
    search_engine_ = p.search_engine_;
    search_engine_version_ = p.search_engine_version_;
    search_parameters_ = p.search_parameters_;
    date_ = p.date_;

    protein_score_type_ = p.protein_score_type_;
    higher_score_better_ = p.higher_score_better_;
    protein_significance_threshold_ = p.protein_significance_threshold_;
    MetaInfoInterface::operator=(static_cast<MetaInfoInterface>(p));
  }

} // namespace OpenMS
