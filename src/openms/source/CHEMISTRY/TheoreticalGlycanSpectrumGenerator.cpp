// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/MonosaccharideDB.h>
#include <OpenMS/CHEMISTRY/MzPAF.h>
#include <OpenMS/CHEMISTRY/TheoreticalGlycanSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <set>

namespace OpenMS
{
namespace
{
  using Generator = TheoreticalGlycanSpectrumGenerator;
  using Composition = Generator::Composition;
  using Fragment = Generator::Fragment;
  using IonType = Generator::IonType;

  [[noreturn]] void invalid(const std::string& message)
  {
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
  }

  bool nonnegative(const EmpiricalFormula& formula)
  {
    return std::all_of(formula.begin(), formula.end(), [](const auto& atom) { return atom.second >= 0; });
  }

  void validateFormula(const EmpiricalFormula& formula)
  {
    if (formula.getCharge() != 0 || formula.getNumberOfAtoms() <= 0 || ! nonnegative(formula))
    {
      invalid("Expected a nonempty neutral formula with nonnegative atom counts");
    }
  }

  struct Component
  {
    Composition::Monosaccharide monosaccharide;
    EmpiricalFormula formula;
    int count = 0;
  };
  using Components = std::map<std::string, Component>;

  Components normalize(const Composition& composition)
  {
    Components result;
    for (const auto& [mono, count] : composition.components)
    {
      if (count < 0) { invalid("Negative glycan residue count"); }
      if (count == 0) { continue; }
      std::string key;
      Component component;
      if (const auto* symbol = std::get_if<std::string>(&mono))
      {
        const auto& entry = MonosaccharideDB::getInstance()->getMonosaccharideOrThrow(*symbol);
        key = entry.symbol;
        component.monosaccharide = key;
        component.formula = EmpiricalFormula(entry.formula);
      }
      else
      {
        const auto& tag = std::get<ProForma::FormulaTag>(mono);
        if (tag.charge.value_or(0) != 0) { invalid("Charged glycan residue formula"); }
        component.formula = EmpiricalFormula(tag.formula_string);
        key = "Formula(" + component.formula.toString() + ")";
        component.monosaccharide = ProForma::FormulaTag {component.formula.toString(), std::nullopt};
      }
      validateFormula(component.formula);
      auto [it, inserted] = result.emplace(key, component);
      if (it->second.count > std::numeric_limits<int>::max() - count) { invalid("Glycan residue count overflow"); }
      it->second.count += count;
    }
    return result;
  }

  Composition compositionOf(const Components& components)
  {
    Composition result;
    for (const auto& [key, component] : components)
    {
      if (component.count > 0) { result.components.emplace_back(component.monosaccharide, component.count); }
    }
    return result;
  }

  EmpiricalFormula formulaOf(const Components& components)
  {
    EmpiricalFormula formula;
    for (const auto& [key, component] : components)
    {
      formula += component.formula * component.count;
    }
    return formula;
  }

  Size sizeOf(const Components& components)
  {
    Size result = 0;
    for (const auto& [key, component] : components)
    {
      result += component.count;
    }
    return result;
  }

  std::string nameOf(const Components& components)
  {
    std::string result;
    for (const auto& [key, component] : components)
    {
      if (component.count > 0) { result += key + std::to_string(component.count); }
    }
    return result.empty() ? "0" : result;
  }

  bool isSubset(const Components& subset, const Components& all)
  {
    for (const auto& [key, component] : subset)
    {
      const auto it = all.find(key);
      if (it == all.end() || it->second.count < component.count) { return false; }
    }
    return true;
  }

  void sortFragments(std::vector<Fragment>& fragments)
  {
    std::stable_sort(fragments.begin(), fragments.end(), [](const auto& a, const auto& b) { return a.getMZ() < b.getMZ(); });
  }
} // namespace

double TheoreticalGlycanSpectrumGenerator::Fragment::getMZ() const
{
  if (charge <= 0 || ! std::isfinite(neutral_mass) || neutral_mass <= 0) { invalid("Invalid fragment mass or charge"); }
  return neutral_mass / charge + Constants::PROTON_MASS_U;
}

std::string TheoreticalGlycanSpectrumGenerator::Fragment::getAnnotation() const
{
  // mzPAF named compounds are bracket-delimited, and its tokenizer discards whitespace.
  if (name.empty() || std::any_of(name.begin(), name.end(), [](unsigned char c) { return std::isspace(c) || c == '[' || c == ']' || c == '\0'; }))
  {
    invalid("Fragment names must be nonempty and contain no whitespace, square brackets or NUL characters");
  }
  MzPAFAnnotation annotation;
  annotation.ion_series = MzPAFIonSeries::NAMED;
  annotation.named_compound = name;
  annotation.charge = charge;
  return MzPAF::toString(annotation);
}

TheoreticalGlycanSpectrumGenerator::TheoreticalGlycanSpectrumGenerator() = default;

TheoreticalGlycanSpectrumGenerator::TheoreticalGlycanSpectrumGenerator(const Options& options)
{
  setOptions(options);
}

void TheoreticalGlycanSpectrumGenerator::setOptions(const Options& options)
{
  if (options.min_composition_size == 0 || options.min_composition_size > options.max_composition_size || options.max_cleavages == 0
      || options.max_fragments == 0 || options.max_states == 0 || options.min_charge <= 0 || options.max_charge < options.min_charge
      || options.min_oxonium_charge <= 0 || options.max_oxonium_charge < options.min_oxonium_charge
      || options.max_charge == std::numeric_limits<Int>::max() || options.max_oxonium_charge == std::numeric_limits<Int>::max())
  {
    invalid("Invalid glycan size, charge or resource limits");
  }
  for (const auto& loss : options.neutral_losses)
  {
    validateFormula(loss);
  }
  for (const auto& [symbol, losses] : options.specific_neutral_losses)
  {
    MonosaccharideDB::getInstance()->getMonosaccharideOrThrow(symbol);
    for (const auto& loss : losses)
    {
      validateFormula(loss);
    }
  }
  for (const auto& [series, retention] : options.peptide_retention)
  {
    if (std::string("bycz").find(series) == std::string::npos) { invalid("Retention rule must target b/y/c/z"); }
    for (const auto& stub : retention.stubs)
    {
      if (normalize(stub).empty()) { invalid("Use stripped retention instead of an empty stub"); }
    }
  }
  options_ = options;
}

const TheoreticalGlycanSpectrumGenerator::Options& TheoreticalGlycanSpectrumGenerator::getOptions() const
{
  return options_;
}

std::vector<Fragment> TheoreticalGlycanSpectrumGenerator::getFragments(const Composition& composition) const
{
  return generate_(composition, nullptr, nullptr, 0, FragmentationMethod::HCD);
}

std::vector<Fragment> TheoreticalGlycanSpectrumGenerator::getFragments(const GlycanStructure& structure) const
{
  return generate_(structure.getComposition(), &structure, nullptr, 0, FragmentationMethod::HCD);
}

std::vector<Fragment> TheoreticalGlycanSpectrumGenerator::getGlycopeptideFragments(const AASequence& peptide,
                                                                                   const Composition& composition,
                                                                                   Size attachment_position,
                                                                                   FragmentationMethod method) const
{
  return generate_(composition, nullptr, &peptide, attachment_position, method);
}

std::vector<Fragment> TheoreticalGlycanSpectrumGenerator::getGlycopeptideFragments(const AASequence& peptide,
                                                                                   const GlycanStructure& structure,
                                                                                   Size attachment_position,
                                                                                   FragmentationMethod method) const
{
  return generate_(structure.getComposition(), &structure, &peptide, attachment_position, method);
}

std::vector<Fragment> TheoreticalGlycanSpectrumGenerator::generate_(const Composition& composition,
                                                                    const GlycanStructure* structure,
                                                                    const AASequence* peptide,
                                                                    Size site,
                                                                    FragmentationMethod method) const
{
  if (method != FragmentationMethod::HCD && method != FragmentationMethod::ETD && method != FragmentationMethod::ETHCD)
  {
    invalid("Unknown glycopeptide fragmentation method");
  }
  const auto all = normalize(composition);
  if (peptide && (site >= peptide->size() || all.empty())) { invalid("Expected a glycan and a valid peptide attachment site"); }
  if (all.empty()) { return {}; }
  const auto water = EmpiricalFormula("H2O");
  const auto reducing_formula = peptide ? peptide->getFormula() : water;
  const double reducing_mass = peptide ? peptide->getMonoWeight() : water.getMonoWeight();
  std::vector<Fragment> result;
  std::set<std::pair<std::string, Int>> emitted;
  Size states = 0;
  auto visit = [&]() {
    if (++states > options_.max_states) { invalid("Glycan enumeration exceeds max_states"); }
  };
  const std::string residue = peptide ? peptide->getResidue(site).getOneLetterCode() : "";
  const std::string site_name = peptide ? ";site=" + residue + std::to_string(site + 1) : "";

  // Each loss is checked against the retained material, not the precursor glycan.
  auto emit = [&](Fragment fragment, const Components& retained, const EmpiricalFormula& formula, double mass, bool losses, bool attach = true) {
    fragment.composition = compositionOf(retained);
    if (peptide && attach)
    {
      fragment.attachment_position = site;
      fragment.attachment_residue = residue;
      fragment.name += site_name;
    }
    std::set<EmpiricalFormula> loss_formulas {EmpiricalFormula()};
    if (losses)
    {
      loss_formulas.insert(options_.neutral_losses.begin(), options_.neutral_losses.end());
      for (const auto& [symbol, specific] : options_.specific_neutral_losses)
      {
        const auto& canonical = MonosaccharideDB::getInstance()->getMonosaccharideOrThrow(symbol).symbol;
        const auto it = retained.find(canonical);
        if (it != retained.end() && it->second.count > 0) { loss_formulas.insert(specific.begin(), specific.end()); }
      }
    }
    const bool diagnostic = fragment.ion_type == IonType::DIAGNOSTIC;
    const Int min_charge = diagnostic ? options_.min_oxonium_charge : options_.min_charge;
    const Int max_charge = diagnostic ? options_.max_oxonium_charge : options_.max_charge;
    const std::string base_name = fragment.name;
    for (const auto& loss : loss_formulas)
    {
      if (! nonnegative(formula - loss)) { continue; }
      const double neutral_mass = mass - loss.getMonoWeight();
      if (neutral_mass <= 0 || ! std::isfinite(neutral_mass)) { continue; }
      const std::string name = base_name + (loss.isEmpty() ? "" : "-" + loss.toString());
      for (Int charge = min_charge; charge <= max_charge; ++charge)
      {
        visit();
        if (! emitted.emplace(name, charge).second) { continue; }
        if (result.size() >= options_.max_fragments) { invalid("Glycan spectrum exceeds max_fragments"); }
        fragment.name = name;
        fragment.neutral_mass = neutral_mass;
        fragment.charge = charge;
        result.push_back(fragment);
      }
    }
  };

  // ETD alone predominantly preserves the glycan; EThcD also has collision products.
  if (! peptide || method != FragmentationMethod::ETD)
  {
    if (options_.add_diagnostic_ions)
    {
      auto diagnostic = [&](const std::string& symbol, const std::vector<std::string>& losses) {
        const auto& canonical = MonosaccharideDB::getInstance()->getMonosaccharideOrThrow(symbol).symbol;
        const auto it = all.find(canonical);
        if (it == all.end()) { return; }
        Components single {{canonical, it->second}};
        single.begin()->second.count = 1;
        for (const auto& loss : losses)
        {
          const auto formula = it->second.formula - EmpiricalFormula(loss);
          Fragment fragment;
          fragment.name = "glycan:diagnostic:" + canonical + (loss.empty() ? "" : "-" + loss);
          emit(fragment, single, formula, formula.getMonoWeight(), false);
        }
      };
      // Neutral loss formulas reproduce the established 204/186/168/144/138/126 HexNAc series.
      diagnostic("HexNAc", {"", "H2O", "H4O2", "C2H4O2", "CH6O3", "C2H6O3"});
      diagnostic("Hex", {"", "H2O", "H4O2"});
      diagnostic("Fuc", {""});
      diagnostic("dHex", {""});
      diagnostic("Neu5Ac", {"", "H2O"});
      diagnostic("Neu5Gc", {"", "H2O"});
      if (all.count("HexNAc") && all.count("Hex"))
      {
        Components pair {{"HexNAc", all.at("HexNAc")}, {"Hex", all.at("Hex")}};
        for (auto& [key, component] : pair)
        {
          component.count = 1;
        }
        const auto formula = formulaOf(pair);
        Fragment fragment;
        fragment.name = "glycan:diagnostic:Hex1HexNAc1";
        emit(fragment, pair, formula, formula.getMonoWeight(), false);
      }
    }

    auto series
      = [&](const Components& retained, bool reducing, const std::string& label, std::optional<Size> root, const std::vector<Size>& branches) {
          const auto residues = formulaOf(retained);
          const auto base = reducing ? residues + reducing_formula : residues;
          const double mass = residues.getMonoWeight() + (reducing ? reducing_mass : 0.0);
          Fragment fragment;
          fragment.root_cleavage = root;
          fragment.branch_cleavages = branches;
          if (reducing ? options_.add_y_ions : options_.add_b_ions)
          {
            fragment.ion_type = reducing ? IonType::Y : IonType::B;
            fragment.name = std::string("glycan:") + (reducing ? "Y:" : "B:") + label;
            emit(fragment, retained, base, mass, true);
          }
          if (reducing ? options_.add_z_ions : options_.add_c_ions)
          {
            // A composition Z ion represents a single Z cleavage, without inferred branch count.
            const int water_count = reducing ? -static_cast<int>(std::max<Size>(1, branches.size())) : 1;
            const auto delta = water * water_count;
            if (! nonnegative(base + delta)) { return; }
            fragment.ion_type = reducing ? IonType::Z : IonType::C;
            fragment.name = std::string("glycan:") + (reducing ? "Z:" : "C:") + label;
            emit(fragment, retained, base + delta, mass + delta.getMonoWeight(), true);
          }
        };

    if (options_.add_b_ions || options_.add_y_ions || options_.add_c_ions || options_.add_z_ions)
    {
      const bool structural = structure && options_.allow_structural;
      series({}, true, structural ? "tree:cuts=0;retained=0" : "composition:0", std::nullopt,
             structural ? std::vector<Size> {0} : std::vector<Size> {});
      if (! structural)
      {
        // Counts, rather than permutations of residues, avoid duplicate compositions.
        std::vector<std::string> keys;
        for (const auto& [key, component] : all)
        {
          keys.push_back(key);
        }
        if (keys.size() > 256) { invalid("Composition fragmentation supports at most 256 distinct residue types"); }
        Components retained = all;
        for (auto& [key, component] : retained)
        {
          component.count = 0;
        }
        const Size total = sizeOf(all);
        const Size maximum = std::min(total, options_.max_composition_size);
        std::function<void(Size, Size)> enumerate = [&](Size index, Size count) {
          visit();
          if (index == keys.size())
          {
            if (count < options_.min_composition_size) { return; }
            const auto label = "composition:" + nameOf(retained);
            series(retained, false, label, std::nullopt, {});
            if (count < total) { series(retained, true, label, std::nullopt, {}); }
            return;
          }
          const auto& key = keys[index];
          const Size limit = std::min<Size>(all.at(key).count, maximum - count);
          for (Size n = 0; n <= limit; ++n)
          {
            retained.at(key).count = static_cast<int>(n);
            enumerate(index + 1, count + n);
          }
          retained.at(key).count = 0;
        };
        enumerate(0, 0);
      }
      else
      {
        const auto& nodes = structure->getNodes();
        // Bound recursion depth as well as the search, before descending into a tree.
        if (nodes.size() > 256) { invalid("Structural fragmentation supports at most 256 monosaccharides"); }
        std::vector<bool> included(nodes.size(), false);
        for (Size root = 0; root < nodes.size(); ++root)
        {
          std::fill(included.begin(), included.end(), false);
          included[root] = true;
          std::vector<Size> cuts;
          std::function<void(Size)> enumerate = [&](Size index) {
            visit();
            if (index == nodes.size())
            {
              Composition kept;
              for (Size i = root; i < nodes.size(); ++i)
              {
                if (included[i]) { kept.components.emplace_back(nodes[i].monosaccharide, 1); }
              }
              const auto retained = normalize(kept);
              std::string cut_names;
              for (Size cut : cuts)
              {
                cut_names += (cut_names.empty() ? "" : ",") + std::to_string(cut);
              }
              const auto label = "tree:root=" + std::to_string(root) + ";cuts=" + cut_names + ";retained=" + nameOf(retained);
              if (cuts.size() < options_.max_cleavages && (cuts.empty() || options_.add_internal_fragments))
              {
                series(retained, false, label, root, cuts);
              }
              if (root == 0 && ! cuts.empty()) { series(retained, true, label, std::nullopt, cuts); }
              return;
            }
            if (! included[*nodes[index].parent])
            {
              included[index] = false;
              enumerate(index + 1);
              return;
            }
            included[index] = true;
            enumerate(index + 1);
            included[index] = false;
            const Size allowance = root == 0 ? options_.max_cleavages : options_.max_cleavages - 1;
            if (cuts.size() < allowance && (root == 0 || options_.add_internal_fragments))
            {
              cuts.push_back(index);
              enumerate(index + 1);
              cuts.pop_back();
            }
          };
          enumerate(root + 1);
        }
      }
    }
  }

  if (peptide)
  {
    std::vector<std::pair<char, Residue::ResidueType>> series_types;
    if (method != FragmentationMethod::ETD)
    {
      series_types.emplace_back('b', Residue::BIon);
      series_types.emplace_back('y', Residue::YIon);
    }
    if (method != FragmentationMethod::HCD)
    {
      series_types.emplace_back('c', Residue::CIon);
      // Zp1Ion is the radical z+1 form, the main electron-transfer fragment.
      series_types.emplace_back('z', Residue::Zp1Ion);
    }
    for (const auto& [series, type] : series_types)
    {
      PeptideRetention retention;
      if (method != FragmentationMethod::HCD)
      {
        retention.intact = true;
        retention.stripped = false;
      }
      else if (all.count("HexNAc")) { retention.stubs.push_back(Composition {{{std::string("HexNAc"), 1}}}); }
      if (const auto it = options_.peptide_retention.find(series); it != options_.peptide_retention.end()) { retention = it->second; }
      std::vector<Components> retained_forms;
      if (retention.stripped) { retained_forms.emplace_back(); }
      if (retention.intact) { retained_forms.push_back(all); }
      for (const auto& stub : retention.stubs)
      {
        const auto normalized = normalize(stub);
        if (! isSubset(normalized, all)) { invalid("Peptide glycan stub exceeds the attached composition"); }
        retained_forms.push_back(normalized);
      }
      for (Size ordinal = 1; ordinal < peptide->size(); ++ordinal)
      {
        visit();
        const bool prefix = series == 'b' || series == 'c';
        const auto part = prefix ? peptide->getPrefix(ordinal) : peptide->getSuffix(ordinal);
        const auto part_formula = part.getFormula(type);
        const double part_mass = part.getMonoWeight(type);
        const bool contains_site = prefix ? site < ordinal : site >= peptide->size() - ordinal;
        const auto forms = contains_site ? retained_forms : std::vector<Components> {{}};
        for (const auto& retained : forms)
        {
          Fragment fragment;
          fragment.ion_type = IonType::PEPTIDE;
          fragment.name = "peptide:" + std::string(1, series) + std::to_string(ordinal);
          if (contains_site) { fragment.name += ";glycan=" + nameOf(retained); }
          const auto formula = formulaOf(retained);
          emit(fragment, retained, part_formula + formula, part_mass + formula.getMonoWeight(), false, contains_site);
        }
      }
    }
  }
  sortFragments(result);
  return result;
}

MSSpectrum TheoreticalGlycanSpectrumGenerator::toSpectrum(const std::vector<Fragment>& fragments)
{
  MSSpectrum spectrum;
  spectrum.setMSLevel(2);
  spectrum.getStringDataArrays().resize(1);
  spectrum.getIntegerDataArrays().resize(1);
  auto& names = spectrum.getStringDataArrays().front();
  auto& charges = spectrum.getIntegerDataArrays().front();
  names.setName(Constants::UserParam::IonNames);
  charges.setName("Charges");
  for (const auto& fragment : fragments)
  {
    Peak1D peak;
    peak.setMZ(fragment.getMZ());
    peak.setIntensity(1.0);
    spectrum.push_back(peak);
    names.push_back(fragment.getAnnotation());
    charges.push_back(fragment.charge);
  }
  spectrum.sortByPosition();
  return spectrum;
}
} // namespace OpenMS
