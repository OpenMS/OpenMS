// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <map>
#include <set>
#include <exception>
#include <typeinfo>
#include <utility>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  FLASHHelperClasses::PrecalculatedAveragine::PrecalculatedAveragine(const double min_mass, const double max_mass, const double delta, CoarseIsotopePatternGenerator& generator,
                                                                     const bool use_RNA_averagine, const double decoy_iso_distance) :
      mass_interval_(delta),
      min_mass_(min_mass)
  {
    // the mass grid
    std::vector<double> masses;
    int i = 0;
    while (true)
    {
      double mass = i * mass_interval_;
      i++;
      if (mass < min_mass)
      {
        continue;
      }
      if (mass > max_mass)
      {
        break;
      }
      masses.push_back(mass);
    }

    const Size bin_count = masses.size();
    apex_index_.resize(bin_count);
    right_count_from_apex_.resize(bin_count);
    left_count_from_apex_.resize(bin_count);
    average_mono_mass_difference_.resize(bin_count);
    abundant_mono_mass_difference_.resize(bin_count);
    isotopes_.resize(bin_count);
    snr_mul_factor_.resize(bin_count);
    std::vector<std::exception_ptr> errors(bin_count);

    // The bins are independent of each other (the left and right counts are made cumulative below), so they are calculated in parallel,
    // each thread with its own copy of the generator.
    const auto calculate_bin = [&](const Size index, CoarseIsotopePatternGenerator& local_generator) {
      const double mass = masses[index];
      auto iso = use_RNA_averagine ? local_generator.estimateFromRNAMonoWeight(mass) : local_generator.estimateFromPeptideMonoWeight(mass);

      if (decoy_iso_distance > 0)
      {
        auto decoy_iso(iso);

        for (Size k = 0; k < iso.size(); k++)
        {
          decoy_iso[k].setMZ(iso[k].getMZ() * decoy_iso_distance);
        }
        decoy_iso.sortByMass();
        decoy_iso.renormalize();

        iso = decoy_iso;
      }

      const double min_pwr = .9999;
      const Size min_iso_length = 2;
      const int min_left_right_count = 2;
      double total_pwr = .0;
      size_t most_abundant_index_ = 0;
      double most_abundant_int = 0;

      /// sum of squared intensities to see the total power of isotope pattern. The range of isotope pattern is
      /// determined so those within range cover min_pwr of the total power.
      for (Size k = 0; k < iso.size(); k++)
      {
        total_pwr += iso[k].getIntensity() * iso[k].getIntensity();
        if (most_abundant_int >= iso[k].getIntensity())
        {
          continue;
        }
        most_abundant_int = iso[k].getIntensity();
        most_abundant_index_ = k;
      }

      int left_count = 0;
      int right_count = (int)iso.size() - 1;
      int trim_count = 0;
      double pwr = 0;
      while (iso.size() - trim_count > min_iso_length && left_count < right_count)
      {
        double lint = iso[left_count].getIntensity();
        double rint = iso[right_count].getIntensity();

        const bool trim_left = lint < rint;
        const double candidate_pwr = trim_left ? lint * lint : rint * rint;
        if (total_pwr - (pwr + candidate_pwr) < total_pwr * min_pwr)
        {
          break;
        }

        // Subtract power only for isotopes actually removed from the pattern.
        pwr += candidate_pwr;
        trim_count++;
        if (trim_left)
        {
          iso[left_count].setIntensity(0);
          left_count++;
        }
        else
        {
          iso[right_count].setIntensity(0); // for trimming
          right_count--;
        }
      }
      total_pwr -= pwr;

      left_count = (int)most_abundant_index_ - left_count;
      right_count = right_count - (int)most_abundant_index_;

      double intensity_sum = 0;
      for (auto& k : iso)
      {
        float ori_int = k.getIntensity();
        k.setIntensity(ori_int / (float)sqrt(total_pwr));
        intensity_sum += k.getIntensity();
      }
      left_count = left_count < min_left_right_count ? min_left_right_count : left_count;
      right_count = right_count < min_left_right_count ? min_left_right_count : right_count;

      apex_index_[index] = most_abundant_index_;
      right_count_from_apex_[index] = right_count;
      left_count_from_apex_[index] = left_count;
      average_mono_mass_difference_[index] = iso.averageMass() - iso[0].getMZ();
      abundant_mono_mass_difference_[index] = iso.getMostAbundant().getMZ() - iso[0].getMZ();
      isotopes_[index] = iso;
      snr_mul_factor_[index] = intensity_sum * intensity_sum;
    };

    // a generator of a derived class would be sliced by a copy: it runs serially on the caller's generator
    if (typeid(generator) != typeid(CoarseIsotopePatternGenerator))
    {
      for (Size index = 0; index < bin_count; ++index)
      {
        calculate_bin(index, generator);
      }
    }
    else
    {
      int worker_count = 1;
#ifdef _OPENMP
      if (! omp_in_parallel()) { worker_count = std::max(1, (int)std::min(bin_count, (Size)omp_get_max_threads())); }
#endif
      std::vector<CoarseIsotopePatternGenerator> generators(worker_count, generator); // copied here: no allocation inside the region
#pragma omp parallel for num_threads(worker_count) schedule(static) if (worker_count > 1)
      for (SignedSize index = 0; index < (SignedSize)bin_count; ++index)
      {
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        try
        {
          calculate_bin((Size)index, generators[worker]);
        }
        catch (...)
        {
          errors[index] = std::current_exception();
        }
      }
    }

    // the first error in mass order; the left and right counts are the maxima up to each mass
    int max_left_count = 0;
    int max_right_count = 0;
    for (Size index = 0; index < bin_count; ++index)
    {
      if (errors[index]) { std::rethrow_exception(errors[index]); }
      max_left_count = std::max(max_left_count, left_count_from_apex_[index]);
      max_right_count = std::max(max_right_count, right_count_from_apex_[index]);
      left_count_from_apex_[index] = max_left_count;
      right_count_from_apex_[index] = max_right_count;
    }
  }

  Size FLASHHelperClasses::PrecalculatedAveragine::massToIndex_(const double mass) const
  {
    Size i = (Size)round(std::max(.0, mass - min_mass_) / mass_interval_);
    i = std::min(i, isotopes_.size() - 1);
    return i;
  }

  IsotopeDistribution FLASHHelperClasses::PrecalculatedAveragine::get(const double mass) const
  {
    return isotopes_[massToIndex_(mass)];
  }

  size_t FLASHHelperClasses::PrecalculatedAveragine::getMaxIsotopeIndex() const
  {
    return max_isotope_index_;
  }

  Size FLASHHelperClasses::PrecalculatedAveragine::getLeftCountFromApex(const double mass) const
  {
    return (Size)left_count_from_apex_[massToIndex_(mass)];
  }

  double FLASHHelperClasses::PrecalculatedAveragine::getAverageMassDelta(const double mass) const
  {
    return average_mono_mass_difference_[massToIndex_(mass)];
  }

  double FLASHHelperClasses::PrecalculatedAveragine::getMostAbundantMassDelta(const double mass) const
  {
    return abundant_mono_mass_difference_[massToIndex_(mass)];
  }

  double FLASHHelperClasses::PrecalculatedAveragine::getSNRMultiplicationFactor(const double mass) const
  {
    return snr_mul_factor_[massToIndex_(mass)];
  }

  Size FLASHHelperClasses::PrecalculatedAveragine::getRightCountFromApex(const double mass) const
  {
    return (Size)right_count_from_apex_[massToIndex_(mass)];
  }

  Size FLASHHelperClasses::PrecalculatedAveragine::getApexIndex(const double mass) const
  {
    return apex_index_[massToIndex_(mass)];
  }

  Size FLASHHelperClasses::PrecalculatedAveragine::getLastIndex(const double mass) const
  {
    Size index = massToIndex_(mass);
    return apex_index_[index] + right_count_from_apex_[index];
  }

  void FLASHHelperClasses::PrecalculatedAveragine::setMaxIsotopeIndex(const int index)
  {
    max_isotope_index_ = index;
  }

  FLASHHelperClasses::LogMzPeak::LogMzPeak(const Peak1D& peak, const bool positive) :
      mz(peak.getMZ()), intensity(peak.getIntensity()), logMz(getLogMz(peak.getMZ(), positive)), abs_charge(0), is_positive(positive), isotopeIndex(0)
  {
  }

  void FLASHHelperClasses::loadExperiment(const std::string& file, PeakMap& map, const PeakFileOptions& options, ProgressLogger::LogType log)
  {
    if (FileHandler::getTypeByFileName(file) == FileTypes::RAW)
    {
#ifdef WITH_THERMO_RAW
      FileHandler file_handler;
      file_handler.setOptions(options);
      file_handler.loadExperiment(file, map, {FileTypes::RAW}, log);
      return;
#else
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "'" + file + "' is a Thermo .raw file, which this build cannot read (OpenMS built without WITH_THERMO_RAW): "
                                        "convert it to mzML first, e.g. with ThermoRawFileParser.");
#endif
    }
    MzMLFile mzml;
    mzml.setLogType(log);
    mzml.setOptions(options);
    mzml.load(file, map);
  }

  const ResidueModification* FLASHHelperClasses::resolveVariableModification(const std::string& spec_in, bool& protein_n_term)
  {
    protein_n_term = false;
    std::string spec = spec_in;
    StringUtils::trim(spec);
    // a name with the site "Protein N-term <residue>" that the database does not know: the same modification at any N-terminus, if the
    // database has it there (looked up without the name as given, whose failed lookup ModificationsDB would log as a warning)
    const std::string tag = "(Protein N-term ";
    const Size at = spec.find(tag);
    if (at != std::string::npos && ! ModificationsDB::getInstance()->has(spec))
    {
      const std::string any_n_term = spec.substr(0, at) + "(N-term " + spec.substr(at + tag.size());
      if (ModificationsDB::getInstance()->has(any_n_term))
      {
        const ResidueModification* mod = resolveModification(any_n_term);
        if (mod->getTermSpecificity() == ResidueModification::N_TERM)
        {
          protein_n_term = true;
          return mod;
        }
      }
    }
    return resolveModification(spec);
  }

  const ResidueModification* FLASHHelperClasses::resolveModification(const std::string& spec_in)
  {
    std::string spec = spec_in;
    StringUtils::trim(spec);
    const auto* db = ModificationsDB::getInstance();
    std::string upper = spec;
    StringUtils::toUpper(upper);
    // a Unimod accession: its number, with or without "UNIMOD:", with its site in parentheses or bare. Only digits before the site, so
    // that names starting with a digit (e.g. 4-ONE) stay names.
    std::string rest = StringUtils::hasPrefix(upper, "UNIMOD:") ? spec.substr(7) : spec;
    {
      std::string head = rest.substr(0, std::min(rest.find('('), rest.size()));
      StringUtils::trim(head);
      const bool digits = ! head.empty() && std::all_of(head.begin(), head.end(), [](char c) { return std::isdigit((unsigned char)c); });
      if (! StringUtils::hasPrefix(upper, "UNIMOD:") && ! digits) return db->getModification(spec);
    }
    std::string number = rest, site;
    const Size open = rest.find('(');
    if (open != std::string::npos)
    {
      const Size close = rest.find(')', open);
      if (close == std::string::npos) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "missing ')' in the modification", spec);
      number = rest.substr(0, open);
      site = rest.substr(open + 1, close - open - 1);
    }
    StringUtils::trim(number);
    StringUtils::trim(site);
    if (number.empty() || ! std::all_of(number.begin(), number.end(), [](char c) { return std::isdigit((unsigned char)c); }))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "a Unimod accession is its number (or UNIMOD:<number>), e.g. 35 (M)", spec);
    }
    std::string residue;
    ResidueModification::TermSpecificity term = ResidueModification::NUMBER_OF_TERM_SPECIFICITY;
    if (! site.empty())
    {
      std::string lower = site;
      StringUtils::toLower(lower);
      if (lower == "n-term") term = ResidueModification::N_TERM;
      else if (lower == "c-term") term = ResidueModification::C_TERM;
      else if (lower == "protein n-term") term = ResidueModification::PROTEIN_N_TERM;
      else if (lower == "protein c-term") term = ResidueModification::PROTEIN_C_TERM;
      else if (site.size() == 1 && std::isalpha((unsigned char)site[0]))
      {
        residue = std::string(1, (char)std::toupper((unsigned char)site[0]));
        term = ResidueModification::ANYWHERE;
      }
      else
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "the site of a modification is a residue letter, N-term, C-term, Protein N-term or Protein C-term", spec);
      }
    }
    std::set<const ResidueModification*> found;
    db->searchModifications(found, "UniMod:" + number, residue, term);
    if (found.empty()) throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, spec);
    // the sites found (a site can have several records, e.g. from Unimod and PSI-MOD): one site is the modification
    auto siteOf = [](const ResidueModification* m) {
      const auto t = m->getTermSpecificity();
      if (t == ResidueModification::ANYWHERE) return std::string(1, m->getOrigin());
      return std::string(m->getTermSpecificityName(t));
    };
    std::map<std::string, const ResidueModification*> by_site;
    for (const auto* m : found)
    {
      const std::string key = siteOf(m);
      auto it = by_site.find(key);
      // prefer the Unimod record of that number
      if (it == by_site.end() || (it->second->getUniModRecordId() != StringUtils::toInt32(number) && m->getUniModRecordId() == StringUtils::toInt32(number))) by_site[key] = m;
    }
    if (by_site.size() > 1)
    {
      std::string list;
      for (const auto& [key, m] : by_site) list += (list.empty() ? "" : ", ") + number + " (" + key + ")";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "the modification has several sites; give one: " + list, spec);
    }
    return by_site.begin()->second;
  }

  double FLASHHelperClasses::getLogMz(const double mz, const bool positive)
  {
    return std::log(mz - getChargeMass(positive));
  }

  bool FLASHHelperClasses::IsobaricQuantities::empty() const
  {
    return quantities.empty();
  }
  FLASHHelperClasses::Tag::Tag(const std::string& seq, double n_mass, double c_mass, const std::vector<double>& mzs, const std::vector<int>& scores,  int scan) :
      seq_(std::move(seq)), n_mass_(n_mass), c_mass_(c_mass), mzs_(mzs), scores_(scores), scan_(scan), length_(mzs.size() - 1)
  {
    upper_seq_ = StringUtils::toUpper(seq_);
  }

  FLASHHelperClasses::Tag::Tag(std::string seq, double n_mass, double c_mass, std::vector<double>& mzs, std::vector<int>& scores, int scan, Size length) :
      seq_(std::move(seq)), n_mass_(n_mass), c_mass_(c_mass), mzs_(mzs), scores_(scores), scan_(scan), length_(length)
  {
    upper_seq_ = StringUtils::toUpper(seq_);
  }

  const std::string& FLASHHelperClasses::Tag::getSequence() const
  {
    return seq_;
  }

  const std::string& FLASHHelperClasses::Tag::getUppercaseSequence() const
  {
    return upper_seq_;
  }

  Size FLASHHelperClasses::Tag::getLength() const
  {
    return length_;
  }


  const std::vector<double>& FLASHHelperClasses::Tag::getMzs() const
  {
    return mzs_;
  }

  double FLASHHelperClasses::Tag::getNtermMass() const
  {
    return n_mass_;
  }

  double FLASHHelperClasses::Tag::getCtermMass() const
  {
    return c_mass_;
  }

  int FLASHHelperClasses::Tag::getScore() const
  {
    return std::accumulate(scores_.begin(), scores_.end(), 0);
  }

  int FLASHHelperClasses::Tag::getScore(int pos) const
  {
    if (pos < 0 || pos >= scores_.size()) return 0;
    return scores_[pos];
  }

  int FLASHHelperClasses::Tag::getScan() const
  {
    return scan_;
  }

  bool FLASHHelperClasses::Tag::operator<(const Tag& a) const
  {
    if (this->length_ == a.length_)
    {
      if (this->seq_ == a.seq_)
      {
        if (this->c_mass_ >= 0 && a.c_mass_ >= 0) // c term tag
        {
          return this->c_mass_ < a.c_mass_;
        }
        else if (this->n_mass_ >= 0 && a.n_mass_ >= 0)
        {
          return this->n_mass_ < a.n_mass_;
        }
        else
        {
          return this->n_mass_ < a.n_mass_;
        }
      }
      return this->seq_ < a.seq_;
    }
    return this->length_ < a.length_;
  }

  bool FLASHHelperClasses::Tag::operator>(const Tag& a) const
  {
    if (this->length_ == a.length_)
    {
      if (this->seq_ == a.seq_)
      {
        if (this->c_mass_ >= 0 && a.c_mass_ >= 0) // c term tag
        {
          return this->c_mass_ > a.c_mass_;
        }
        else if (this->n_mass_ >= 0 && a.n_mass_ >= 0)
        {
          return this->n_mass_ > a.n_mass_;
        }
        else
        {
          return this->n_mass_ > a.n_mass_;
        }
      }
      return this->seq_ > a.seq_;
    }
    return this->length_ > a.length_;
  }

  bool FLASHHelperClasses::Tag::operator==(const Tag& a) const
  {
    return this->seq_ == a.seq_ && this->n_mass_ == a.n_mass_ && this->c_mass_ == a.c_mass_;
  }


  std::string FLASHHelperClasses::Tag::toString() const
  {
    std::string ret;
    if (n_mass_ >= 0)
      ret = '[' + std::to_string(n_mass_) + "]\t";
    ret += seq_;
    if (c_mass_ >= 0)
      ret += "\t[" + std::to_string(c_mass_) + "]";
    ret += "\tscore : " + std::to_string(getScore()) + "\tmzs : ";

    for (auto mz : mzs_)
    {
      ret += std::to_string(mz) + " ";
    }

    return ret;
  }

  FLASHHelperClasses::MassTag::MassTag(const std::string& full_seq, const std::vector<int>& positions, int scan)
  {
    if (positions.size() < 2) return;
    seq_ = std::move(full_seq);
    positions_ = positions;
    scan = scan_;
    //is_n_term_ = positions[0] < positions[1];


// String seq_, upper_seq_;
    //       float rt_ = -1;
    //       int scan_;
    //       std::vector<int> positions_;
    //       int index_;
    //       Size length_;

  }

  bool FLASHHelperClasses::MassTag::operator<(const MassTag& a) const
  {
    (void)a;
    return false;
  }

  bool FLASHHelperClasses::MassTag::operator>(const MassTag& a) const
  {
    (void)a;
    return false;
  }

  bool FLASHHelperClasses::MassTag::operator==(const MassTag& a) const
  {
    (void)a;
    return false;
  }

    const std::string& FLASHHelperClasses::MassTag::getSequence() const
    {
      return seq_;
    }

    Size FLASHHelperClasses::MassTag::getLength() const
    {
      return length_;
    }

    int FLASHHelperClasses::MassTag::getScan() const
    {
      return scan_;
    }

    // ProForma: read and written by OpenMS::ProForma; the masses of the tags as computed for the top-down tools
  namespace
  {
    using PFM = FLASHHelperClasses::ProFormaModification;

    /// the residues ProForma is read and written with here
    const std::string proforma_residues = "ACDEFGHIKLMNPQRSTVWYUO";

    [[noreturn]] void proformaFail(const std::string& input, const std::string& message)
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input, message);
    }

    /// the mass of a glycan composition from the monosaccharide names of ProForma; names that follow each other without a count
    /// (HexHexNAc) count once each, the count of a component applying to its last name
    double glycanMass(const ProForma::GlycanComposition& glycan, const std::string& input)
    {
      // longest names first, as several share a prefix (HexNAc, HexN, HexA, HexS, Hex)
      static const std::vector<std::pair<std::string, std::string>> monosaccharides = {
        {"HexNAcS", "C8H13NO8S"}, {"HexNAc", "C8H13NO5"}, {"NeuAc", "C11H17NO8"}, {"NeuGc", "C11H17NO9"}, {"HexS", "C6H10O8S"},
        {"HexP", "C6H11O8P"},     {"HexN", "C6H11NO4"},   {"HexA", "C6H8O6"},     {"dHex", "C6H10O4"},   {"Hex", "C6H10O5"},
        {"Fuc", "C6H10O4"},       {"Neu", "C9H15NO7"},    {"Pen", "C5H8O4"},      {"Kdn", "C9H14O8"},    {"Phosphate", "HPO3"},
        {"Sulfate", "SO3"}};
      double mass = 0;
      for (const auto& [monosaccharide, count] : glycan.components)
      {
        const auto* name = std::get_if<std::string>(&monosaccharide);
        if (name == nullptr) proformaFail(input, "a formula in a glycan composition is not supported");
        Size pos = 0;
        while (pos < name->size())
        {
          if ((*name)[pos] == ' ')
          {
            pos++;
            continue;
          }
          bool found = false;
          for (const auto& [mono_name, formula] : monosaccharides)
          {
            if (name->compare(pos, mono_name.size(), mono_name) != 0) continue;
            pos += mono_name.size();
            mass += (pos >= name->size() ? count : 1) * EmpiricalFormula(formula).getMonoWeight();
            found = true;
            break;
          }
          if (! found) proformaFail(input, "unknown monosaccharide in glycan composition");
        }
      }
      return mass;
    }

    /// a modification name or accession looked up in ModificationsDB (residue-specific first)
    bool lookUpMass(const std::string& name, char residue, double& mass)
    {
      const auto* db = ModificationsDB::getInstance();
      if (! db->has(name)) return false; // searchModifications would log every miss
      std::set<const ResidueModification*> mods;
      if (residue != 0) db->searchModifications(mods, name, std::string(1, residue));
      if (mods.empty()) db->searchModifications(mods, name);
      if (mods.empty()) return false;
      mass = (*mods.begin())->getDiffMonoMass();
      return true;
    }

    /// a number with an optional sign; false if the whole string is not a number
    bool toNumber(const std::string& s, double& value)
    {
      if (s.empty()) return false;
      try
      {
        Size used = 0;
        value = std::stod(s, &used);
        return used == s.size();
      }
      catch (const std::exception&)
      {
        return false;
      }
    }

    /// a name without a score suffix: Phospho(0.9) -> Phospho
    std::string withoutScore(std::string name)
    {
      StringUtils::trim(name);
      if (! name.empty() && name.back() == ')')
      {
        const Size open = name.rfind('(');
        double score;
        if (open != std::string::npos && toNumber(name.substr(open + 1, name.size() - open - 2), score)) name = name.substr(0, open);
      }
      StringUtils::trim(name);
      return name;
    }

    /// cross-links and branches are not supported
    void checkLabels(const ProForma::Modification& mod, const std::string& input)
    {
      for (const auto& [tag, label] : mod.alternatives)
      {
        if (label.has_value() && (label->type == ProForma::Label::Type::CROSSLINK || label->type == ProForma::Label::Type::BRANCH))
        {
          proformaFail(input, "cross-links and branches are not supported");
        }
      }
    }

    /// the mass of one alternative of a tag; false if it gives none (INFO, a localisation group label alone, a position, an accession of
    /// RESID, XL-MOD or GNOme, which are not in ModificationsDB)
    bool alternativeMass(const ProForma::ModificationTag& tag, char residue, double& mass, const std::string& input)
    {
      if (const auto* md = std::get_if<ProForma::MassDelta>(&tag))
      {
        switch (md->source)
        {
          case ProForma::MassDelta::Source::NONE:
          case ProForma::MassDelta::Source::OBS:
            mass = md->mass;
            return true;
          case ProForma::MassDelta::Source::U:
          case ProForma::MassDelta::Source::M:
            return lookUpMass(md->original_text, residue, mass); // a name lookup of what follows U: or M:
          default:
            return false;
        }
      }
      if (const auto* formula = std::get_if<ProForma::FormulaTag>(&tag))
      {
        if (formula->charge.has_value()) proformaFail(input, "invalid formula");
        std::string f = formula->formula_string;
        StringUtils::substitute(f, " ", "");
        try
        {
          mass = EmpiricalFormula(f).getMonoWeight();
        }
        catch (const Exception::BaseException&)
        {
          proformaFail(input, "invalid formula");
        }
        return true;
      }
      if (const auto* glycan = std::get_if<ProForma::GlycanComposition>(&tag))
      {
        mass = glycanMass(*glycan, input);
        return true;
      }
      if (const auto* cv = std::get_if<ProForma::CvAccession>(&tag))
      {
        if (cv->database == ProForma::CvDatabase::UNIMOD) return lookUpMass("UniMod:" + cv->accession, residue, mass);
        if (cv->database == ProForma::CvDatabase::MOD) return lookUpMass("MOD:" + cv->accession, residue, mass);
        return false;
      }
      if (const auto* named = std::get_if<ProForma::NamedMod>(&tag))
      {
        if (named->cv_hint.has_value() && named->cv_hint != ProForma::CvDatabase::UNIMOD && named->cv_hint != ProForma::CvDatabase::MOD)
        {
          return false;
        }
        const std::string name = withoutScore(named->name);
        if (name.empty()) return false;
        if (toNumber(name, mass)) return true;
        // a name without source, or one with a colon of its own (Cation:Na): Unimod, then PSI-MOD
        return lookUpMass(name, residue, mass);
      }
      return false; // INFO, a label alone, a position constraint
    }

    /// the mass of a modification: of the first alternative that gives one
    double modificationMass(const ProForma::Modification& mod, char residue, const std::string& input)
    {
      checkLabels(mod, input);
      for (const auto& [tag, label] : mod.alternatives)
      {
        double mass = 0;
        if (alternativeMass(tag, residue, mass, input)) return mass;
      }
      proformaFail(input, "no mass for the modification");
    }

    /// a ProForma string of one residue with @p mod, written by OpenMS::ProForma
    std::string writeOnResidue(const ProForma::Modification& mod)
    {
      ProForma::Peptidoform pf;
      ProForma::SequenceElement element;
      element.amino_acid = 'A';
      element.modifications.push_back(mod);
      pf.sequence.emplace_back(std::move(element));
      return ProForma::toString(pf);
    }

    /// the tag content of @p mod as written by OpenMS::ProForma (without its brackets)
    std::string tagText(const ProForma::Modification& mod)
    {
      const std::string written = writeOnResidue(mod);
      return written.substr(2, written.size() - 3);
    }

    /// the modification of a tag content read by OpenMS::ProForma; false if it cannot be read as one modification on a residue
    bool readModification(const std::string& tag, ProForma::Modification& mod)
    {
      ProForma::Peptidoform pf;
      try
      {
        pf = ProForma::parse("A[" + tag + "]");
      }
      catch (const Exception::BaseException&)
      {
        return false;
      }
      if (pf.sequence.size() != 1 || ! pf.global_mods.empty() || ! pf.unlocalised_mods.empty() || ! pf.labile_mods.empty() ||
          ! pf.n_term_mods.empty() || ! pf.c_term_mods.empty() || pf.name.has_value() || pf.charge.has_value())
      {
        return false;
      }
      const auto* element = std::get_if<ProForma::SequenceElement>(&pf.sequence[0]);
      if (element == nullptr || element->modifications.size() != 1) return false;
      mod = element->modifications[0];
      return true;
    }

    /// the ProForma::Modification written as @p m: a mass as formatProFormaMass writes it, or the tag as given (each alternative as
    /// OpenMS::ProForma reads it if it writes it back the same, else as a name with that text)
    ProForma::Modification toModification(const PFM& m)
    {
      ProForma::Modification mod;
      if (m.tag.empty())
      {
        ProForma::MassDelta delta;
        delta.mass = m.mass;
        delta.original_text = FLASHHelperClasses::formatProFormaMass(m.mass);
        mod.alternatives.emplace_back(delta, std::nullopt);
        return mod;
      }
      if (readModification(m.tag, mod) && tagText(mod) == m.tag) return mod;
      mod = ProForma::Modification();
      std::vector<std::string> parts;
      StringUtils::split(m.tag, '|', parts);
      for (const auto& part : parts)
      {
        ProForma::Modification alternative;
        if (readModification(part, alternative) && alternative.alternatives.size() == 1 && tagText(alternative) == part)
        {
          mod.alternatives.push_back(alternative.alternatives[0]);
          continue;
        }
        ProForma::NamedMod named;
        named.name = part;
        mod.alternatives.emplace_back(named, std::nullopt);
      }
      return mod;
    }

    /// a modification on residues start..end with its mass and its tag as written
    PFM residueModification(int start, int end, const ProForma::Modification& mod, double mass, bool global = false)
    {
      PFM m;
      m.site = PFM::Site::RESIDUE;
      m.start = start;
      m.end = end;
      m.mass = mass;
      m.tag = tagText(mod);
      m.global = global;
      return m;
    }

    /// a modification at a site other than residues
    PFM siteModification(PFM::Site site, const ProForma::Modification& mod, double mass)
    {
      PFM m;
      m.site = site;
      m.mass = mass;
      m.tag = tagText(mod);
      return m;
    }
  } // namespace

  std::string FLASHHelperClasses::formatProFormaMass(double mass)
  {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%+.6f", mass);
    std::string s(buffer);
    while (s.size() > 3 && s.back() == '0' && s[s.size() - 2] != '.') s.pop_back();
    return s;
  }

  double FLASHHelperClasses::getProFormaTagMass(const std::string& tag, char residue)
  {
    ProForma::Modification mod;
    if (readModification(tag, mod)) return modificationMass(mod, residue, tag);
    // a name or an accession that OpenMS::ProForma does not read, e.g. a Unimod name starting with a digit (2-methylbutyryl), with a '+'
    // (Met-loss+Acetyl), or OpenMS's spelling of a Unimod accession (UniMod:35): looked up as it is
    double mass = 0;
    if (lookUpMass(withoutScore(tag), residue, mass)) return mass;
    proformaFail(tag, "no mass for the modification");
  }

  ProForma::PeptidoformIon FLASHHelperClasses::parseProForma(const std::string& text)
  {
    std::string cleaned = text;
    for (const std::string& alt : {std::string("|-]"), std::string("|-|")})
    {
      for (size_t pos = cleaned.find(alt); pos != std::string::npos; pos = cleaned.find(alt, pos)) { cleaned.erase(pos, 2); }
    }
    return ProForma::parseIon(cleaned);
  }

  std::vector<FLASHHelperClasses::ProFormaModification> FLASHHelperClasses::getProFormaModifications(const ProForma::PeptidoformIon& ion,
                                                                                                      std::string& sequence)
  {
    const std::string input = ProForma::toString(ion);
    if (ion.is_chimeric || ion.chains.size() != 1) proformaFail(input, "chimeric spectra (+) and several chains (//) are not supported");
    if (ion.charge.has_value() && std::get_if<int>(&ion.charge.value()) == nullptr) proformaFail(input, "only an integer charge is supported after '/'");
    const ProForma::Peptidoform& pf = ion.chains[0];
    if (ion.name.has_value() || pf.name.has_value()) proformaFail(input, "names are not supported");
    if (pf.charge.has_value()) proformaFail(input, "a charge per chain is not supported");

    // the sequence
    sequence.clear();
    auto add_residue = [&](char residue) {
      if (proforma_residues.find(residue) == std::string::npos) proformaFail(input, std::string("unsupported residue '") + residue + "'");
      sequence += residue;
    };
    for (const auto& section : pf.sequence)
    {
      if (const auto* element = std::get_if<ProForma::SequenceElement>(&section)) add_residue(element->amino_acid);
      else if (const auto* range = std::get_if<ProForma::ModifiedRange>(&section))
      {
        if (range->elements.empty()) proformaFail(input, "empty range");
        if (range->modifications.empty()) proformaFail(input, "a range must be followed by its modification");
        for (const auto& element : range->elements) add_residue(element.amino_acid);
      }
      else proformaFail(input, "ambiguous residues are not supported");
    }
    if (sequence.empty()) proformaFail(input, "no sequence");

    // a residue tag without mass ([#g1], [INFO:x]) marks only another possible site
    auto has_mass = [&input](const ProForma::Modification& mod) {
      checkLabels(mod, input);
      for (const auto& [tag, label] : mod.alternatives)
      {
        if (! std::holds_alternative<ProForma::InfoTag>(tag)) return true;
      }
      return false;
    };

    std::vector<ProFormaModification> mods;
    for (const auto& entry : pf.global_mods)
    {
      const auto* global = std::get_if<ProForma::GlobalModification>(&entry);
      if (global == nullptr) proformaFail(input, "isotope labels <...> are not supported");
      std::vector<char> residues;
      for (auto location : global->locations)
      {
        StringUtils::trim(location);
        if (location.size() != 1 || proforma_residues.find(location[0]) == std::string::npos)
        {
          proformaFail(input, "global modification on '" + location + "' is not supported");
        }
        residues.push_back(location[0]);
      }
      if (residues.empty()) proformaFail(input, "global modification without residues");
      for (int i = 0; i < (int)sequence.size(); i++)
      {
        if (std::find(residues.begin(), residues.end(), sequence[i]) == residues.end()) continue;
        mods.push_back(residueModification(i, i, global->modification, modificationMass(global->modification, sequence[i], input), true));
      }
    }
    for (const auto& unlocalised : pf.unlocalised_mods)
    {
      for (int c = 0; c < unlocalised.occurrence.value_or(1); c++)
      {
        for (const auto& mod : unlocalised.modifications) mods.push_back(siteModification(PFM::Site::UNKNOWN, mod, modificationMass(mod, 0, input)));
      }
    }
    for (const auto& labile : pf.labile_mods)
    {
      mods.push_back(siteModification(PFM::Site::LABILE, labile.modification, modificationMass(labile.modification, 0, input)));
    }
    for (const auto& mod : pf.n_term_mods) mods.push_back(siteModification(PFM::Site::N_TERM, mod, modificationMass(mod, sequence[0], input)));
    int index = 0;
    for (const auto& section : pf.sequence)
    {
      if (const auto* element = std::get_if<ProForma::SequenceElement>(&section))
      {
        for (const auto& mod : element->modifications)
        {
          if (has_mass(mod)) mods.push_back(residueModification(index, index, mod, modificationMass(mod, sequence[index], input)));
        }
        index++;
        continue;
      }
      const auto& range = std::get<ProForma::ModifiedRange>(section);
      const int start = index;
      for (const auto& element : range.elements)
      {
        for (const auto& mod : element.modifications)
        {
          if (has_mass(mod)) mods.push_back(residueModification(index, index, mod, modificationMass(mod, sequence[index], input)));
        }
        index++;
      }
      for (const auto& mod : range.modifications)
      {
        if (! has_mass(mod)) continue;
        const int end = index - 1;
        mods.push_back(residueModification(start, end, mod, modificationMass(mod, start == end ? sequence[start] : 0, input)));
      }
    }
    for (const auto& mod : pf.c_term_mods) mods.push_back(siteModification(PFM::Site::C_TERM, mod, modificationMass(mod, sequence.back(), input)));
    return mods;
  }

  ProForma::Peptidoform FLASHHelperClasses::toPeptidoform(const std::string& sequence, const std::vector<ProFormaModification>& modifications)
  {
    for (char c : sequence)
    {
      if (proforma_residues.find(c) == std::string::npos) proformaFail(sequence, std::string("unsupported residue '") + c + "'");
    }
    // residue modifications: within the sequence, ranges the same or disjoint (a range may hold no other range)
    std::vector<ProFormaModification> mods;
    for (auto m : modifications)
    {
      if (m.site == ProFormaModification::Site::RESIDUE)
      {
        if (m.start > m.end) std::swap(m.start, m.end);
        if (m.start < 0 || m.end >= (int)sequence.size())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "modification outside the sequence",
                                        StringUtils::toStr(m.start) + "-" + StringUtils::toStr(m.end));
        }
        if (m.end > m.start)
        {
          for (const auto& other : mods)
          {
            if (other.site != ProFormaModification::Site::RESIDUE || other.end == other.start || (other.start == m.start && other.end == m.end)) continue;
            if (other.start <= m.end && m.start <= other.end)
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "overlapping ranges",
                                            StringUtils::toStr(m.start) + "-" + StringUtils::toStr(m.end));
            }
          }
        }
      }
      mods.push_back(m);
    }

    ProForma::Peptidoform pf;
    // global fixed modifications, one per tag, with their residues
    std::map<std::string, std::pair<ProForma::Modification, std::set<char>>> globals;
    for (const auto& m : mods)
    {
      if (! m.global) continue;
      const auto mod = toModification(m);
      auto& global = globals.emplace(writeOnResidue(mod), std::make_pair(mod, std::set<char>())).first->second;
      global.second.insert(sequence[m.start]);
    }
    for (const auto& [text, global] : globals)
    {
      ProForma::GlobalModification g;
      g.modification = global.first;
      for (char r : global.second) g.locations.emplace_back(1, r);
      pf.global_mods.emplace_back(std::move(g));
    }
    ProForma::UnlocalisedMod unknown;
    for (const auto& m : mods)
    {
      if (m.site == ProFormaModification::Site::UNKNOWN) unknown.modifications.push_back(toModification(m));
      if (m.site == ProFormaModification::Site::LABILE) pf.labile_mods.push_back(ProForma::LabileModification {toModification(m)});
      if (m.site == ProFormaModification::Site::N_TERM) pf.n_term_mods.push_back(toModification(m));
      if (m.site == ProFormaModification::Site::C_TERM) pf.c_term_mods.push_back(toModification(m));
    }
    if (! unknown.modifications.empty()) pf.unlocalised_mods.push_back(std::move(unknown));

    // the residues, a range with its residues and then its modifications
    auto element_at = [&](int i) {
      ProForma::SequenceElement element;
      element.amino_acid = sequence[i];
      for (const auto& m : mods)
      {
        if (m.site == ProFormaModification::Site::RESIDUE && ! m.global && m.start == i && m.end == i) element.modifications.push_back(toModification(m));
      }
      return element;
    };
    for (int i = 0; i < (int)sequence.size(); i++)
    {
      int range_end = -1;
      for (const auto& m : mods)
      {
        if (m.site == ProFormaModification::Site::RESIDUE && ! m.global && m.start == i && m.end > i) range_end = m.end;
      }
      if (range_end < 0)
      {
        pf.sequence.emplace_back(element_at(i));
        continue;
      }
      ProForma::ModifiedRange range;
      for (int j = i; j <= range_end; j++) range.elements.push_back(element_at(j));
      for (const auto& m : mods)
      {
        if (m.site == ProFormaModification::Site::RESIDUE && ! m.global && m.start == i && m.end == range_end) range.modifications.push_back(toModification(m));
      }
      pf.sequence.emplace_back(std::move(range));
      i = range_end;
    }
    return pf;
  }

} // namespace OpenMS
