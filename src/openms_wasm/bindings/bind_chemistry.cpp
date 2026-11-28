// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

// Sequences, formulas, isotope distributions, modifications, digestion and theoretical
// spectra (CHEMISTRY/).

#include "binding_utils.h"

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <map>

using namespace emscripten;

namespace OpenMS::WASM
{
  namespace
  {
    //------------------------------------------------------------------------------
    // Residues and modifications as plain objects (they live in OpenMS' databases)

    JS::ResidueObject residueToJS(const Residue& r)
    {
      val o = val::object();
      o.set("name", r.getName());
      o.set("oneLetterCode", r.getOneLetterCode());
      o.set("threeLetterCode", r.getThreeLetterCode());
      o.set("formula", r.getFormula().toString());
      o.set("monoWeight", r.getMonoWeight());
      o.set("averageWeight", r.getAverageWeight());
      o.set("isModified", r.isModified());
      o.set("modificationName", r.isModified() ? val(r.getModificationName()) : val::null());
      return o;
    }

    JS::ModificationObject modificationToJS(const ResidueModification* m)
    {
      if (m == nullptr) return val::null();
      val o = val::object();
      o.set("id", m->getId());
      o.set("fullId", m->getFullId());
      o.set("fullName", m->getFullName());
      o.set("unimodAccession", m->getUniModAccession());
      o.set("psiModAccession", m->getPSIMODAccession());
      o.set("origin", std::string(1, m->getOrigin()));
      o.set("termSpecificity", m->getTermSpecificityName());
      o.set("diffMonoMass", m->getDiffMonoMass());
      o.set("monoMass", m->getMonoMass());
      o.set("diffFormula", m->getDiffFormula().toString());
      return o;
    }

    ResidueModification::TermSpecificity termSpecificityFromJS(const JS::TermSpecificityInput& v)
    {
      if (isNullish(v)) return ResidueModification::NUMBER_OF_TERM_SPECIFICITY;
      const std::string name = v.as<std::string>();
      if (name.empty() || name == "any") return ResidueModification::NUMBER_OF_TERM_SPECIFICITY;
      ResidueModification rm;
      rm.setTermSpecificity(name); // throws for unknown names
      return rm.getTermSpecificity();
    }

    //------------------------------------------------------------------------------
    // AASequence

    AASequence sequenceFromString(const std::string& s) { return AASequence::fromString(s); }
    AASequence sequenceFromStringPermissive(const std::string& s, bool permissive) { return AASequence::fromString(s, permissive); }

    double sequenceMonoWeight(const AASequence& s) { return s.getMonoWeight(); }
    double sequenceMonoWeightTyped(const AASequence& s, Residue::ResidueType type, int charge) { return s.getMonoWeight(type, charge); }
    double sequenceAverageWeight(const AASequence& s) { return s.getAverageWeight(); }
    double sequenceAverageWeightTyped(const AASequence& s, Residue::ResidueType type, int charge) { return s.getAverageWeight(type, charge); }
    double sequenceMZ(const AASequence& s, int charge) { return s.getMZ(charge); }
    double sequenceMZTyped(const AASequence& s, int charge, Residue::ResidueType type) { return s.getMZ(charge, type); }
    EmpiricalFormula sequenceFormula(const AASequence& s) { return s.getFormula(); }
    EmpiricalFormula sequenceFormulaTyped(const AASequence& s, Residue::ResidueType type, int charge) { return s.getFormula(type, charge); }
    AASequence sequencePrefix(const AASequence& s, Size n) { return s.getPrefix(n); }
    AASequence sequenceSuffix(const AASequence& s, Size n) { return s.getSuffix(n); }
    AASequence sequenceSubsequence(const AASequence& s, Size start, UInt n) { return s.getSubsequence(start, n); }
    JS::ResidueObject sequenceResidue(const AASequence& s, Size i)
    {
      if (i >= s.size()) throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, i, s.size());
      return residueToJS(s.getResidue(i));
    }
    JS::ResidueList sequenceResidues(const AASequence& s)
    {
      val result = val::array();
      for (Size i = 0; i < s.size(); ++i) result.call<void>("push", residueToJS(s.getResidue(i)));
      return result;
    }
    void sequenceSetModification(AASequence& s, Size i, const std::string& mod) { s.setModification(i, mod); }
    void sequenceSetNTermMod(AASequence& s, const std::string& mod) { s.setNTerminalModification(mod); }
    void sequenceSetCTermMod(AASequence& s, const std::string& mod) { s.setCTerminalModification(mod); }
    std::string sequenceNTermModName(const AASequence& s) { return s.getNTerminalModificationName(); }
    std::string sequenceCTermModName(const AASequence& s) { return s.getCTerminalModificationName(); }
    bool sequenceHasSubsequence(const AASequence& s, const AASequence& other) { return s.hasSubsequence(other); }
    bool sequenceHasPrefix(const AASequence& s, const AASequence& other) { return s.hasPrefix(other); }
    bool sequenceHasSuffix(const AASequence& s, const AASequence& other) { return s.hasSuffix(other); }
    bool sequenceEquals(const AASequence& s, const AASequence& other) { return s == other; }
    JS::ElementCounts sequenceAAFrequencies(const AASequence& s)
    {
      std::map<std::string, Size> table;
      s.getAAFrequencies(table);
      val o = val::object();
      for (const auto& [aa, n] : table) o.set(aa, static_cast<double>(n));
      return o;
    }

    //------------------------------------------------------------------------------
    // EmpiricalFormula

    EmpiricalFormula formulaFromString(const std::string& s) { return EmpiricalFormula(s); }
    JS::ElementCounts formulaToObject(const EmpiricalFormula& f)
    {
      val o = val::object();
      for (const auto& [element, count] : f.toMap()) o.set(element, count);
      return o;
    }
    EmpiricalFormula formulaPlus(const EmpiricalFormula& a, const EmpiricalFormula& b) { return a + b; }
    EmpiricalFormula formulaMinus(const EmpiricalFormula& a, const EmpiricalFormula& b) { return a - b; }
    double formulaNumberOfAtoms(const EmpiricalFormula& f) { return static_cast<double>(f.getNumberOfAtoms()); }
    bool formulaEquals(const EmpiricalFormula& a, const EmpiricalFormula& b) { return a == b; }
    IsotopeDistribution formulaIsotopes(const EmpiricalFormula& f, const IsotopePatternGenerator& g) { return f.getIsotopeDistribution(g); }

    //------------------------------------------------------------------------------
    // IsotopeDistribution

    JS::IsotopePeaks isotopesPeaks(const IsotopeDistribution& d)
    {
      val o = val::object();
      o.set("mz", toTypedArray<double>(d.getContainer(), [](const Peak1D& p) { return p.getMZ(); }));
      o.set("intensity", toTypedArray<double>(d.getContainer(), [](const Peak1D& p) { return p.getIntensity(); }));
      return o;
    }
    JS::PeakObjectList isotopesContainer(const IsotopeDistribution& d)
    {
      return toArray(d.getContainer(), [](const Peak1D& p) {
        val o = val::object();
        o.set("mz", p.getMZ());
        o.set("intensity", p.getIntensity());
        return o;
      });
    }
    JS::OptionalPeak isotopesMostAbundant(const IsotopeDistribution& d)
    {
      if (d.size() == 0) return val::null();
      const Peak1D p = d.getMostAbundant();
      val o = val::object();
      o.set("mz", p.getMZ());
      o.set("intensity", p.getIntensity());
      return o;
    }
    double isotopesSize(const IsotopeDistribution& d) { return static_cast<double>(d.size()); }

    IsotopeDistribution generatorRun(const IsotopePatternGenerator& g, const EmpiricalFormula& f) { return g.run(f); }
    void coarseSetMaxIsotope(CoarseIsotopePatternGenerator& g, Size n) { g.setMaxIsotope(n); }
    IsotopeDistribution coarseFromPeptideWeight(CoarseIsotopePatternGenerator& g, double average_weight) { return g.estimateFromPeptideWeight(average_weight); }

    //------------------------------------------------------------------------------
    // Databases (static functions)

    JS::ResidueObject residueDBGetResidue(const std::string& name)
    {
      const Residue* r = ResidueDB::getInstance()->getResidue(name); // throws for unknown names
      return residueToJS(*r);
    }
    JS::StringArray residueDBGetResidueNames()
    {
      std::vector<std::string> names;
      for (const Residue* r : ResidueDB::getInstance()->getResidues("Natural20")) names.push_back(r->getName());
      std::sort(names.begin(), names.end());
      return toArray(names);
    }

    JS::ModificationObject modDBGetModification(const std::string& name) { return modificationToJS(ModificationsDB::getInstance()->getModification(name)); }
    JS::ModificationObject modDBGetModificationAt(const std::string& name, const std::string& residue, const JS::TermSpecificityInput& term)
    {
      return modificationToJS(ModificationsDB::getInstance()->getModification(name, residue, termSpecificityFromJS(term)));
    }
    bool modDBHas(const std::string& name) { return ModificationsDB::getInstance()->has(name); }
    JS::StringArray modDBAllSearchModifications()
    {
      std::vector<std::string> mods;
      ModificationsDB::getInstance()->getAllSearchModifications(mods);
      return toArray(mods);
    }
    JS::StringArray modDBSearchByDiffMonoMass(double mass, double max_error, const std::string& residue)
    {
      std::vector<std::string> mods;
      ModificationsDB::getInstance()->searchModificationsByDiffMonoMassSorted(mods, mass, max_error, residue);
      return toArray(mods);
    }
    JS::ModificationObject modDBBestByDiffMonoMass(double mass, double max_error, const std::string& residue)
    {
      return modificationToJS(ModificationsDB::getInstance()->getBestModificationByDiffMonoMass(mass, max_error, residue));
    }

    JS::StringArray proteaseDBAllNames()
    {
      std::vector<std::string> names;
      ProteaseDB::getInstance()->getAllNames(names);
      std::sort(names.begin(), names.end());
      return toArray(names);
    }
    JS::EnzymeObject proteaseDBGetEnzyme(const std::string& name)
    {
      const DigestionEnzymeProtein* e = ProteaseDB::getInstance()->getEnzyme(name); // throws for unknown names
      val o = val::object();
      o.set("name", e->getName());
      o.set("regEx", e->getRegEx());
      o.set("regExDescription", e->getRegExDescription());
      o.set("psiId", e->getPSIID());
      return o;
    }

    //------------------------------------------------------------------------------
    // ProteaseDigestion

    void digestionSetMissedCleavages(ProteaseDigestion& d, Size n) { d.setMissedCleavages(n); }
    double digestionGetMissedCleavages(const ProteaseDigestion& d) { return static_cast<double>(d.getMissedCleavages()); }
    std::string digestionGetSpecificity(const ProteaseDigestion& d) { return EnzymaticDigestion::NamesOfSpecificity[d.getSpecificity()]; }
    void digestionSetSpecificity(ProteaseDigestion& d, const std::string& name) { d.setSpecificity(EnzymaticDigestion::getSpecificityByName(name)); }

    JS::AASequenceList digestionDigest(const ProteaseDigestion& d, const AASequence& protein, Size min_length, Size max_length)
    {
      std::vector<AASequence> peptides;
      d.digest(protein, peptides, min_length, max_length);
      return toArray(peptides);
    }

    JS::StringArray digestionDigestUnmodified(const ProteaseDigestion& d, const std::string& protein, Size min_length, Size max_length)
    {
      std::vector<std::string_view> peptides;
      d.digestUnmodified(protein, peptides, min_length, max_length);
      return toArray(peptides, [](std::string_view p) { return std::string(p); });
    }

    double digestionPeptideCount(ProteaseDigestion& d, const AASequence& protein) { return static_cast<double>(d.peptideCount(protein)); }

    //------------------------------------------------------------------------------
    // TheoreticalSpectrumGenerator

    void tsgGetSpectrum(const TheoreticalSpectrumGenerator& g, MSSpectrum& spec, const AASequence& peptide, int min_charge, int max_charge)
    {
      g.getSpectrum(spec, peptide, min_charge, max_charge);
    }
    void tsgGetSpectrumPrecursor(const TheoreticalSpectrumGenerator& g, MSSpectrum& spec, const AASequence& peptide, int min_charge, int max_charge, int precursor_charge)
    {
      g.getSpectrum(spec, peptide, min_charge, max_charge, precursor_charge);
    }
  } // namespace
} // namespace OpenMS::WASM

using namespace OpenMS;
using namespace OpenMS::WASM;

EMSCRIPTEN_BINDINGS(openms_chemistry)
{
  enum_<Residue::ResidueType>("ResidueType")
    .value("Full", Residue::Full)
    .value("Internal", Residue::Internal)
    .value("NTerminal", Residue::NTerminal)
    .value("CTerminal", Residue::CTerminal)
    .value("AIon", Residue::AIon)
    .value("BIon", Residue::BIon)
    .value("CIon", Residue::CIon)
    .value("XIon", Residue::XIon)
    .value("YIon", Residue::YIon)
    .value("ZIon", Residue::ZIon)
    .value("Zp1Ion", Residue::Zp1Ion)
    .value("Zp2Ion", Residue::Zp2Ion)
    .value("Precursor", Residue::Precursor);

  class_<AASequence>("AASequence")
    .constructor<>()
    .constructor<const AASequence&>()
    .class_function("fromString", OMS_FN(&sequenceFromString))
    .class_function("fromString", OMS_FN(&sequenceFromStringPermissive))
    .function("toString", OMS_METHOD(AASequence, &AASequence::toString))
    .function("toUnmodifiedString", OMS_METHOD(AASequence, &AASequence::toUnmodifiedString))
    .function("toUniModString", OMS_METHOD(AASequence, &AASequence::toUniModString))
    .function("size", OMS_METHOD(AASequence, &AASequence::size))
    .function("empty", OMS_METHOD(AASequence, &AASequence::empty))
    .function("equals", OMS_FN(&sequenceEquals))
    .function("getMonoWeight", OMS_FN(&sequenceMonoWeight))
    .function("getMonoWeight", OMS_FN(&sequenceMonoWeightTyped))
    .function("getAverageWeight", OMS_FN(&sequenceAverageWeight))
    .function("getAverageWeight", OMS_FN(&sequenceAverageWeightTyped))
    .function("getMZ", OMS_FN(&sequenceMZ))
    .function("getMZ", OMS_FN(&sequenceMZTyped))
    .function("getFormula", OMS_FN(&sequenceFormula))
    .function("getFormula", OMS_FN(&sequenceFormulaTyped))
    .function("getPrefix", OMS_FN(&sequencePrefix))
    .function("getSuffix", OMS_FN(&sequenceSuffix))
    .function("getSubsequence", OMS_FN(&sequenceSubsequence))
    .function("getResidue", OMS_FN(&sequenceResidue))
    .function("getResidues", OMS_FN(&sequenceResidues))
    .function("getAAFrequencies", OMS_FN(&sequenceAAFrequencies))
    .function("isModified", OMS_METHOD(AASequence, &AASequence::isModified))
    .function("hasNTerminalModification", OMS_METHOD(AASequence, &AASequence::hasNTerminalModification))
    .function("hasCTerminalModification", OMS_METHOD(AASequence, &AASequence::hasCTerminalModification))
    .function("getNTerminalModificationName", OMS_FN(&sequenceNTermModName))
    .function("getCTerminalModificationName", OMS_FN(&sequenceCTermModName))
    .function("setModification", OMS_FN(&sequenceSetModification))
    .function("setNTerminalModification", OMS_FN(&sequenceSetNTermMod))
    .function("setCTerminalModification", OMS_FN(&sequenceSetCTermMod))
    .function("hasSubsequence", OMS_FN(&sequenceHasSubsequence))
    .function("hasPrefix", OMS_FN(&sequenceHasPrefix))
    .function("hasSuffix", OMS_FN(&sequenceHasSuffix))
    ;

  class_<EmpiricalFormula>("EmpiricalFormula")
    .constructor<>()
    .constructor<const EmpiricalFormula&>()
    .class_function("fromString", OMS_FN(&formulaFromString))
    .function("toString", OMS_METHOD(EmpiricalFormula, &EmpiricalFormula::toString))
    .function("toObject", OMS_FN(&formulaToObject))
    .function("getMonoWeight", OMS_METHOD(EmpiricalFormula, &EmpiricalFormula::getMonoWeight))
    .function("getAverageWeight", OMS_METHOD(EmpiricalFormula, &EmpiricalFormula::getAverageWeight))
    .function("getCharge", OMS_METHOD(EmpiricalFormula, &EmpiricalFormula::getCharge))
    .function("setCharge", OMS_METHOD(EmpiricalFormula, &EmpiricalFormula::setCharge))
    .function("getNumberOfAtoms", OMS_FN(&formulaNumberOfAtoms))
    .function("isEmpty", OMS_METHOD(EmpiricalFormula, &EmpiricalFormula::isEmpty))
    .function("isCharged", OMS_METHOD(EmpiricalFormula, &EmpiricalFormula::isCharged))
    .function("contains", OMS_METHOD(EmpiricalFormula, &EmpiricalFormula::contains))
    .function("equals", OMS_FN(&formulaEquals))
    .function("plus", OMS_FN(&formulaPlus))
    .function("minus", OMS_FN(&formulaMinus))
    .function("getIsotopeDistribution", OMS_FN(&formulaIsotopes))
    ;

  class_<IsotopeDistribution>("IsotopeDistribution")
    .constructor<>()
    .function("size", OMS_FN(&isotopesSize))
    .function("getPeaks", OMS_FN(&isotopesPeaks))
    .function("getContainer", OMS_FN(&isotopesContainer))
    .function("getMostAbundant", OMS_FN(&isotopesMostAbundant))
    .function("getMin", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::getMin))
    .function("getMax", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::getMax))
    .function("averageMass", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::averageMass))
    .function("renormalize", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::renormalize))
    .function("merge", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::merge))
    .function("trimRight", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::trimRight))
    .function("trimLeft", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::trimLeft))
    .function("trimIntensities", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::trimIntensities))
    .function("sortByIntensity", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::sortByIntensity))
    .function("sortByMass", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::sortByMass))
    .function("clear", OMS_METHOD(IsotopeDistribution, &IsotopeDistribution::clear))
    ;

  class_<IsotopePatternGenerator>("IsotopePatternGenerator")
    .function("run", OMS_FN(&generatorRun));

  class_<CoarseIsotopePatternGenerator, base<IsotopePatternGenerator>>("CoarseIsotopePatternGenerator")
    .constructor<>()
    .constructor(OMS_CONSTRUCTOR(CoarseIsotopePatternGenerator, Size)) // (constructor<Size>() would take a BigInt in wasm64)
    .constructor(OMS_CONSTRUCTOR(CoarseIsotopePatternGenerator, Size, bool))
    .function("getMaxIsotope", OMS_METHOD(CoarseIsotopePatternGenerator, &CoarseIsotopePatternGenerator::getMaxIsotope))
    .function("setMaxIsotope", OMS_FN(&coarseSetMaxIsotope))
    .function("getRoundMasses", OMS_METHOD(CoarseIsotopePatternGenerator, &CoarseIsotopePatternGenerator::getRoundMasses))
    .function("setRoundMasses", OMS_METHOD(CoarseIsotopePatternGenerator, &CoarseIsotopePatternGenerator::setRoundMasses))
    .function("estimateFromPeptideWeight", OMS_FN(&coarseFromPeptideWeight))
    ;

  class_<FineIsotopePatternGenerator, base<IsotopePatternGenerator>>("FineIsotopePatternGenerator")
    .constructor<>()
    .constructor<double>()
    .constructor<double, bool>()
    .constructor<double, bool, bool>()
    .function("getThreshold", OMS_METHOD(FineIsotopePatternGenerator, &FineIsotopePatternGenerator::getThreshold))
    .function("setThreshold", OMS_METHOD(FineIsotopePatternGenerator, &FineIsotopePatternGenerator::setThreshold))
    .function("getAbsolute", OMS_METHOD(FineIsotopePatternGenerator, &FineIsotopePatternGenerator::getAbsolute))
    .function("setAbsolute", OMS_METHOD(FineIsotopePatternGenerator, &FineIsotopePatternGenerator::setAbsolute))
    .function("getTotalProbability", OMS_METHOD(FineIsotopePatternGenerator, &FineIsotopePatternGenerator::getTotalProbability))
    .function("setTotalProbability", OMS_METHOD(FineIsotopePatternGenerator, &FineIsotopePatternGenerator::setTotalProbability))
    ;

  // The databases are singletons; their functions are static ("ModificationsDB.getModification(...)").
  class_<ResidueDB>("ResidueDB")
    .class_function("getResidue", OMS_FN(&residueDBGetResidue))
    .class_function("getResidueNames", OMS_FN(&residueDBGetResidueNames));

  class_<ModificationsDB>("ModificationsDB")
    .class_function("getModification", OMS_FN(&modDBGetModification))
    .class_function("getModification", OMS_FN(&modDBGetModificationAt))
    .class_function("has", OMS_FN(&modDBHas))
    .class_function("getAllSearchModifications", OMS_FN(&modDBAllSearchModifications))
    .class_function("searchModificationsByDiffMonoMass", OMS_FN(&modDBSearchByDiffMonoMass))
    .class_function("getBestModificationByDiffMonoMass", OMS_FN(&modDBBestByDiffMonoMass));

  class_<ProteaseDB>("ProteaseDB")
    .class_function("getAllNames", OMS_FN(&proteaseDBAllNames))
    .class_function("getEnzyme", OMS_FN(&proteaseDBGetEnzyme));

  class_<ProteaseDigestion>("ProteaseDigestion")
    .constructor<>()
    .function("setEnzyme", OMS_METHOD(ProteaseDigestion, static_cast<void (ProteaseDigestion::*)(const std::string&)>(&ProteaseDigestion::setEnzyme)))
    .function("getEnzymeName", OMS_METHOD(ProteaseDigestion, &ProteaseDigestion::getEnzymeName))
    .function("setMissedCleavages", OMS_FN(&digestionSetMissedCleavages))
    .function("getMissedCleavages", OMS_FN(&digestionGetMissedCleavages))
    .function("setSpecificity", OMS_FN(&digestionSetSpecificity))
    .function("getSpecificity", OMS_FN(&digestionGetSpecificity))
    .function("digest", OMS_FN(&digestionDigest))
    .function("digestUnmodified", OMS_FN(&digestionDigestUnmodified))
    .function("peptideCount", OMS_FN(&digestionPeptideCount))
    ;

  class_<TheoreticalSpectrumGenerator, base<DefaultParamHandler>>("TheoreticalSpectrumGenerator")
    .constructor<>()
    .function("getSpectrum", OMS_FN(&tsgGetSpectrum))
    .function("getSpectrum", OMS_FN(&tsgGetSpectrumPrecursor))
    ;
}
