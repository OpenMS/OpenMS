// pyOpenMS nanobind bindings
// Domain: chemistry

#include "all_casters.h"
#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <OpenMS/CHEMISTRY/AdductInfo.h>
#include <OpenMS/CHEMISTRY/CrossLinksDB.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DecoyGenerator.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeRNA.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/IsoelectricPoint.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSElement.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/RealMassDecomposer.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>
#include <OpenMS/CHEMISTRY/ModificationDefinition.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/MzPAF.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/RNaseDB.h>
#include <OpenMS/CHEMISTRY/RNaseDigestion.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/CHEMISTRY/RibonucleotideDB.h>
#include <OpenMS/CHEMISTRY/SequenceCoverage.h>
#include <OpenMS/CHEMISTRY/Tagger.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <iomanip>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/set.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>
#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_pyopenms_chemistry, m) {
    m.doc() = "pyOpenMS chemistry bindings";

    // -----------------------------------------------------------------------
    // AAIndex
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AAIndex>(m, "AAIndex", "Representation of selected %AAIndex properties")
        .def_static("aliphatic", [](char aa) { return OpenMS::AAIndex::aliphatic(aa); }, "aa"_a)
        .def_static("acidic", [](char aa) { return OpenMS::AAIndex::acidic(aa); }, "aa"_a)
        .def_static("basic", [](char aa) { return OpenMS::AAIndex::basic(aa); }, "aa"_a)
        .def_static("polar", [](char aa) { return OpenMS::AAIndex::polar(aa); }, "aa"_a)
        .def_static("getKHAG800101", [](char aa) { return OpenMS::AAIndex::getKHAG800101(aa); }, "aa"_a)
        .def_static("getVASM830103", [](char aa) { return OpenMS::AAIndex::getVASM830103(aa); }, "aa"_a)
        .def_static("getNADH010106", [](char aa) { return OpenMS::AAIndex::getNADH010106(aa); }, "aa"_a)
        .def_static("getNADH010107", [](char aa) { return OpenMS::AAIndex::getNADH010107(aa); }, "aa"_a)
        .def_static("getWILM950102", [](char aa) { return OpenMS::AAIndex::getWILM950102(aa); }, "aa"_a)
        .def_static("getROBB760107", [](char aa) { return OpenMS::AAIndex::getROBB760107(aa); }, "aa"_a)
        .def_static("getOOBM850104", [](char aa) { return OpenMS::AAIndex::getOOBM850104(aa); }, "aa"_a)
        .def_static("getFAUJ880111", [](char aa) { return OpenMS::AAIndex::getFAUJ880111(aa); }, "aa"_a)
        .def_static("getFINA770101", [](char aa) { return OpenMS::AAIndex::getFINA770101(aa); }, "aa"_a)
        .def_static("getARGP820102", [](char aa) { return OpenMS::AAIndex::getARGP820102(aa); }, "aa"_a)
        .def_static("calculateGB", [](const OpenMS::AASequence& seq, double T) { return OpenMS::AAIndex::calculateGB(seq, T); }, "seq"_a, "T"_a)
        ;

    // -----------------------------------------------------------------------
    // AASequence
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AASequence>(m, "AASequence", 
        R"doc(
Representation of a peptide/protein sequence
This class represents amino acid sequences in OpenMS. An AASequence
instance primarily contains a sequence of residues.
)doc")
        .def("empty", [](const OpenMS::AASequence& self) { return self.empty(); }, "Check if sequence is empty")
        .def("toString", [](const OpenMS::AASequence& self) { return self.toString(); }, "Returns the peptide as string with modifications embedded in brackets")
        .def("toUnmodifiedString", [](const OpenMS::AASequence& self) { return self.toUnmodifiedString(); }, "Returns the peptide as string without any modifications")
        .def("toUniModString", [](const OpenMS::AASequence& self) { return self.toUniModString(); }, "Returns the peptide as string with UniMod-style modifications embedded in brackets")
        .def("setModification", [](OpenMS::AASequence& self, size_t index, const OpenMS::String& modification) { return self.setModification(index, modification); }, "index"_a, "modification"_a, "Sets the modification of the residue at position index. If an empty string is passed replaces the residue with its unmodified version")
        .def("setModification", [](OpenMS::AASequence& self, size_t index, OpenMS::Residue * modification) { return self.setModification(index, modification); }, "index"_a, "modification"_a, "Sets the modification of the residue at position index. If an empty string is passed replaces the residue with its unmodified version")
        .def("setModification", [](OpenMS::AASequence& self, size_t index, OpenMS::ResidueModification * modification) { return self.setModification(index, modification); }, "index"_a, "modification"_a, "Sets the modification of the residue at position index. If an empty string is passed replaces the residue with its unmodified version")
        .def("setModification", [](OpenMS::AASequence& self, size_t index, const OpenMS::ResidueModification& modification) { return self.setModification(index, modification); }, "index"_a, "modification"_a, "Sets the modification of the residue at position index. If an empty string is passed replaces the residue with its unmodified version")
        .def("setModificationByDiffMonoMass", [](OpenMS::AASequence& self, size_t index, double diffMonoMass) { return self.setModificationByDiffMonoMass(index, diffMonoMass); }, "index"_a, "diffMonoMass"_a, "Modifies the residue at index in the sequence and potentially in the ResidueDB")
        .def("setNTerminalModification", [](OpenMS::AASequence& self, const OpenMS::String& modification) { return self.setNTerminalModification(modification); }, "modification"_a, "Sets the N-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)")
        .def("setNTerminalModification", [](OpenMS::AASequence& self, OpenMS::ResidueModification * modification) { return self.setNTerminalModification(modification); }, "modification"_a, "Sets the N-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)")
        .def("setNTerminalModification", [](OpenMS::AASequence& self, const OpenMS::ResidueModification& mod) { return self.setNTerminalModification(mod); }, "mod"_a, "Sets the N-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)")
        .def("setNTerminalModificationByDiffMonoMass", [](OpenMS::AASequence& self, double diffMonoMass, bool protein_term) { return self.setNTerminalModificationByDiffMonoMass(diffMonoMass, protein_term); }, "diffMonoMass"_a, "protein_term"_a, 
            R"doc(
Sets the N-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)
)doc")
        .def("getNTerminalModificationName", [](const OpenMS::AASequence& self) { return self.getNTerminalModificationName(); }, "Returns the name (ID) of the N-terminal modification, or an empty string if none is set")
        .def("getNTerminalModification", [](const OpenMS::AASequence& self) { return self.getNTerminalModification(); }, nb::rv_policy::reference_internal, "Returns a copy of the name N-terminal modification object, or None")
        .def("setCTerminalModification", [](OpenMS::AASequence& self, const OpenMS::String& modification) { return self.setCTerminalModification(modification); }, "modification"_a, "Sets the C-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)")
        .def("setCTerminalModification", [](OpenMS::AASequence& self, OpenMS::ResidueModification * modification) { return self.setCTerminalModification(modification); }, "modification"_a, "Sets the C-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)")
        .def("setCTerminalModification", [](OpenMS::AASequence& self, const OpenMS::ResidueModification& mod) { return self.setCTerminalModification(mod); }, "mod"_a, "Sets the C-terminal modification (by lookup in the mod names of the ModificationsDB). Throws if nothing is found (since the name is not enough information to create a new mod)")
        .def("setCTerminalModificationByDiffMonoMass", [](OpenMS::AASequence& self, double diffMonoMass, bool protein_term) { return self.setCTerminalModificationByDiffMonoMass(diffMonoMass, protein_term); }, "diffMonoMass"_a, "protein_term"_a, 
            R"doc(
Sets the C-terminal modification by the monoisotopic mass difference it introduces (creates a "user-defined" mod if not present)
)doc")
        .def("getCTerminalModificationName", [](const OpenMS::AASequence& self) { return self.getCTerminalModificationName(); }, "Returns the name (ID) of the C-terminal modification, or an empty string if none is set")
        .def("getCTerminalModification", [](const OpenMS::AASequence& self) { return self.getCTerminalModification(); }, nb::rv_policy::reference_internal, "Returns a copy of the name C-terminal modification object, or None")
        .def("getResidue", [](const OpenMS::AASequence& self, size_t index) -> const OpenMS::Residue & { return self.getResidue(index); }, "index"_a, nb::rv_policy::reference_internal, "Returns the residue at position index")
        .def(nb::self + nb::self)
        .def("size", [](const OpenMS::AASequence& self) { return self.size(); }, "Returns the number of residues")
        .def("getPrefix", [](const OpenMS::AASequence& self, size_t index) { return self.getPrefix(index); }, "index"_a, "Returns a peptide sequence of the first index residues")
        .def("getSuffix", [](const OpenMS::AASequence& self, size_t index) { return self.getSuffix(index); }, "index"_a, "Returns a peptide sequence of the last index residues")
        .def("getSubsequence", [](const OpenMS::AASequence& self, size_t index, unsigned int number) { return self.getSubsequence(index, number); }, "index"_a, "number"_a, "Returns a peptide sequence of number residues, beginning at position index")
        .def("has", [](const OpenMS::AASequence& self, const OpenMS::Residue& residue) { return self.has(residue); }, "residue"_a, "Returns true if the peptide contains the given residue")
        .def("hasSubsequence", [](const OpenMS::AASequence& self, const OpenMS::AASequence& peptide) { return self.hasSubsequence(peptide); }, "peptide"_a, "Returns true if the peptide contains the given peptide")
        .def("hasPrefix", [](const OpenMS::AASequence& self, const OpenMS::AASequence& peptide) { return self.hasPrefix(peptide); }, "peptide"_a, "Returns true if the peptide has the given prefix")
        .def("hasSuffix", [](const OpenMS::AASequence& self, const OpenMS::AASequence& peptide) { return self.hasSuffix(peptide); }, "peptide"_a, "Returns true if the peptide has the given suffix")
        .def("hasNTerminalModification", [](const OpenMS::AASequence& self) { return self.hasNTerminalModification(); }, "Predicate which is true if the peptide is N-term modified")
        .def("hasCTerminalModification", [](const OpenMS::AASequence& self) { return self.hasCTerminalModification(); }, "Predicate which is true if the peptide is C-term modified")
        .def("isModified", [](const OpenMS::AASequence& self) { return self.isModified(); }, "Returns true if any of the residues or termini are modified")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("__hash__", [](const OpenMS::AASequence& self) { return std::hash<OpenMS::AASequence>{}(self); })
        .def("__iter__", [](OpenMS::AASequence& self) { return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::AASequence>(), "AASequence_iter", self.begin(), self.end()); })
        .def("__len__", [](OpenMS::AASequence& self) { return self.size(); })
        .def("__getitem__", [](OpenMS::AASequence& self, size_t i) -> const OpenMS::Residue & { 
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)

        .def(nb::init<>(), "Default constructor - creates empty sequence")
        .def(nb::init<const OpenMS::AASequence&>(), "Copy constructor")
        .def("__copy__", [](const OpenMS::AASequence& self) { return OpenMS::AASequence(self); })
        .def("__deepcopy__", [](const OpenMS::AASequence& self, nb::dict) { return OpenMS::AASequence(self); }, "memo"_a)
        .def("__init__", [](OpenMS::AASequence* self, const std::string& s) {
            new (self) OpenMS::AASequence(OpenMS::String(s));
        }, "sequence"_a, "Create AASequence from string (e.g., 'PEPTIDE')")
        .def(nb::init<OpenMS::String>())
        .def(nb::init<OpenMS::String, bool>())

        .def_static("fromString", [](const std::string& s) {
            return OpenMS::AASequence::fromString(OpenMS::String(s));
        }, "s"_a, "Create AASequence from string (deprecated - use AASequence(s) constructor)")

        .def("getMZ", [](const OpenMS::AASequence& self, int charge) {
            return self.getMZ(charge, OpenMS::Residue::ResidueType::Full);
        }, "charge"_a, "Returns the mass-to-charge ratio of the peptide")
        .def("getMZ", [](const OpenMS::AASequence& self, int charge, OpenMS::Residue::ResidueType type) {
            return self.getMZ(charge, type);
        }, "charge"_a, "type"_a, "Returns the mass-to-charge ratio of the peptide with specified residue type")

        .def("getMonoWeight", [](const OpenMS::AASequence& self) {
            return self.getMonoWeight(OpenMS::Residue::ResidueType::Full, 0);
        }, "Returns the monoisotopic weight of the peptide")
        .def("getMonoWeight", [](const OpenMS::AASequence& self, OpenMS::Residue::ResidueType type, int charge) {
            return self.getMonoWeight(type, charge);
        }, "type"_a, "charge"_a, "Returns the monoisotopic weight of the peptide with specified residue type and charge")

        .def("getAverageWeight", [](const OpenMS::AASequence& self) {
            return self.getAverageWeight(OpenMS::Residue::ResidueType::Full, 0);
        }, "Returns the average weight of the peptide")
        .def("getAverageWeight", [](const OpenMS::AASequence& self, OpenMS::Residue::ResidueType type, int charge) {
            return self.getAverageWeight(type, charge);
        }, "type"_a, "charge"_a, "Returns the average weight of the peptide with specified residue type and charge")

        .def_static("fromStringPermissive", [](const std::string& s, bool permissive) {
            return OpenMS::AASequence::fromString(OpenMS::String(s), permissive);
        }, "s"_a, "permissive"_a, "Create AASequence from string with permissive mode")

        .def("getFormula", [](const OpenMS::AASequence& self) {
            return self.getFormula(OpenMS::Residue::ResidueType::Full, 0);
        }, "Returns the empirical formula of the peptide")
        .def("getFormula", [](const OpenMS::AASequence& self, OpenMS::Residue::ResidueType type, int charge) {
            return self.getFormula(type, charge);
        }, "type"_a, "charge"_a, "Returns the empirical formula of the peptide with specified residue type and charge")

        .def("toBracketString", &OpenMS::AASequence::toBracketString,
            "integer_mass"_a = true, "mass_delta"_a = false, "fixed_modifications"_a = std::vector<OpenMS::String>(),
            "Returns the bracket string representation of the peptide")
        .def("getAAFrequencies", [](const OpenMS::AASequence& self) {
            std::map<OpenMS::String, OpenMS::Size> freq;
            self.getAAFrequencies(freq);
            return freq;
        }, "Returns the amino acid frequencies of the peptide")
        .def("__iadd__", [](OpenMS::AASequence& self, const OpenMS::AASequence& rhs) -> OpenMS::AASequence& { return self += rhs; }, "rhs"_a, "In-place concatenation of sequences", nb::rv_policy::reference_internal)
        .def("__repr__", [](const OpenMS::AASequence& self) {
            std::ostringstream oss;
            oss << "AASequence(sequence=" << '\'' << std::string(self.toString()) << '\''
                << ", length=" << self.size()
                << ", mono_mass=" << self.getMonoWeight();
            if (self.isModified()) oss << ", modified=True";
            oss << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::AASequence& self) { return std::string(self.toString()); })
        ;

    // -----------------------------------------------------------------------
    // CoarseIsotopePatternGenerator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CoarseIsotopePatternGenerator>(m, "CoarseIsotopePatternGenerator", "OpenMS class CoarseIsotopePatternGenerator")
        .def(nb::init<size_t, bool>())
        .def("setMaxIsotope", [](OpenMS::CoarseIsotopePatternGenerator& self, const size_t& max_isotope) { return self.setMaxIsotope(max_isotope); }, "max_isotope"_a, "Sets the maximal isotope with 'max_isotope'")
        .def("setRoundMasses", [](OpenMS::CoarseIsotopePatternGenerator& self, bool round_masses) { return self.setRoundMasses(round_masses); }, "round_masses"_a, "Sets the round_masses flag to round masses to integer values (true) or return accurate masses (false)")
        .def("getMaxIsotope", [](const OpenMS::CoarseIsotopePatternGenerator& self) { return self.getMaxIsotope(); }, "Returns the currently set maximum isotope")
        .def("getRoundMasses", [](const OpenMS::CoarseIsotopePatternGenerator& self) { return self.getRoundMasses(); }, "Returns the current value of the flag to round masses to integer values (true) or return accurate masses (false)")
        .def("run", [](const OpenMS::CoarseIsotopePatternGenerator& self, const OpenMS::EmpiricalFormula& ef) { return self.run(ef); }, "ef"_a)
        .def("estimateFromPeptideWeight", [](OpenMS::CoarseIsotopePatternGenerator& self, double average_weight) { return self.estimateFromPeptideWeight(average_weight); }, "average_weight"_a, "Estimate Peptide Isotopedistribution from weight and number of isotopes that should be reported")
        .def("estimateFromPeptideWeightAndS", [](OpenMS::CoarseIsotopePatternGenerator& self, double average_weight, unsigned int S) { return self.estimateFromPeptideWeightAndS(average_weight, S); }, "average_weight"_a, "S"_a, "Estimate peptide IsotopeDistribution from average weight and exact number of sulfurs")
        .def_static("approximateFromPeptideWeight", [](double mass, unsigned int num_peaks, unsigned int charge) { return OpenMS::CoarseIsotopePatternGenerator::approximateFromPeptideWeight(mass, num_peaks, charge); }, "mass"_a, "num_peaks"_a, "charge"_a, 
            R"doc(
Estimate peptide fragment IsotopeDistribution from the precursor's average weight,
number of sulfurs in the precursor, fragment's average weight, number of sulfurs in the fragment,
and a set of isolated precursor isotopes.
)doc")
        .def_static("approximateIntensities", [](double mass, unsigned int num_peaks) { return OpenMS::CoarseIsotopePatternGenerator::approximateIntensities(mass, num_peaks); }, "mass"_a, "num_peaks"_a, 
            R"doc(
Roughly approximate peptide IsotopeDistribution from monoisotopic weight using Poisson distribution.
m/z values approximated by adding one neutron mass (divided by charge) for every peak, starting at
the given monoisotopic weight. Foundation from: Bellew et al, https://dx.doi.org/10.1093/bioinformatics/btl276
This method is around 50 times faster than estimateFromPeptideWeight, but only an approximation.
The following are the intensities of the first 6 peaks generated for a monoisotopic mass of 1000:
estimateFromPeptideWeight:    0.571133000;0.306181000;0.095811100;0.022036900;0.004092170;0.000644568
approximateFromPeptideWeight: 0.573753000;0.318752000;0.088542200;0.016396700;0.002277320;0.000253036
KL divergences of the first 20 intensities of estimateFromPeptideWeight and this approximation range from 4.97E-5 for a
monoisotopic mass of 20 to 0.0144 for a mass of 2500. For comparison, when comparing an observed pattern with a
theoretical ground truth, the observed pattern is said to be an isotopic pattern if the KL between the two is below 0.05
for 2 peaks and below 0.6 for >=6 peaks by Guo Ci Teo et al.
)doc")
        .def("estimateFromRNAWeight", [](OpenMS::CoarseIsotopePatternGenerator& self, double average_weight) { return self.estimateFromRNAWeight(average_weight); }, "average_weight"_a, "Estimate Nucleotide Isotopedistribution from weight")
        .def("estimateFromDNAWeight", [](OpenMS::CoarseIsotopePatternGenerator& self, double average_weight) { return self.estimateFromDNAWeight(average_weight); }, "average_weight"_a, "Estimate Nucleotide Isotopedistribution from weight")
        .def("estimateFromWeightAndComp", [](OpenMS::CoarseIsotopePatternGenerator& self, double average_weight, double C, double H, double N, double O, double S, double P) { return self.estimateFromWeightAndComp(average_weight, C, H, N, O, S, P); }, "average_weight"_a, "C"_a, "H"_a, "N"_a, "O"_a, "S"_a, "P"_a)
        .def("estimateFromWeightAndCompAndS", [](OpenMS::CoarseIsotopePatternGenerator& self, double average_weight, unsigned int S, double C, double H, double N, double O, double P) { return self.estimateFromWeightAndCompAndS(average_weight, S, C, H, N, O, P); }, "average_weight"_a, "S"_a, "C"_a, "H"_a, "N"_a, "O"_a, "P"_a, "Estimate IsotopeDistribution from weight, exact number of sulfurs, and average remaining composition")

        .def(nb::init<>())
        .def(nb::init<OpenMS::Size>(), "max_isotope"_a)
        .def("calcFragmentIsotopeDist", [](const OpenMS::CoarseIsotopePatternGenerator& self, const OpenMS::IsotopeDistribution& fragment_isotope_dist, const OpenMS::IsotopeDistribution& comp_fragment_isotope_dist, const std::set<unsigned int>& precursor_isotopes, double fragment_mono_mass) { return self.calcFragmentIsotopeDist(fragment_isotope_dist, comp_fragment_isotope_dist, precursor_isotopes, fragment_mono_mass); }, "fragment_isotope_dist"_a, "comp_fragment_isotope_dist"_a, "precursor_isotopes"_a, "fragment_mono_mass"_a, "Calculates fragment isotope distribution")
        .def("estimateForFragmentFromPeptideWeight", [](OpenMS::CoarseIsotopePatternGenerator& self, double average_weight_precursor, double average_weight_fragment, const std::set<unsigned int>& precursor_isotopes) { return self.estimateForFragmentFromPeptideWeight(average_weight_precursor, average_weight_fragment, precursor_isotopes); }, "average_weight_precursor"_a, "average_weight_fragment"_a, "precursor_isotopes"_a, "Estimate fragment isotope distribution from peptide weights")
        .def("estimateForFragmentFromPeptideWeightAndS", [](const OpenMS::CoarseIsotopePatternGenerator& self, double average_weight_precursor, unsigned int S_precursor, double average_weight_fragment, unsigned int S_fragment, const std::set<unsigned int>& precursor_isotopes) { return self.estimateForFragmentFromPeptideWeightAndS(average_weight_precursor, S_precursor, average_weight_fragment, S_fragment, precursor_isotopes); }, "average_weight_precursor"_a, "S_precursor"_a, "average_weight_fragment"_a, "S_fragment"_a, "precursor_isotopes"_a, "Estimate fragment isotope distribution from peptide weight and sulfur count")
        .def("estimateForFragmentFromDNAWeight", [](OpenMS::CoarseIsotopePatternGenerator& self, double average_weight_precursor, double average_weight_fragment, const std::set<unsigned int>& precursor_isotopes) { return self.estimateForFragmentFromDNAWeight(average_weight_precursor, average_weight_fragment, precursor_isotopes); }, "average_weight_precursor"_a, "average_weight_fragment"_a, "precursor_isotopes"_a, "Estimate fragment isotope distribution from DNA weights")
        .def("estimateForFragmentFromRNAWeight", [](OpenMS::CoarseIsotopePatternGenerator& self, double average_weight_precursor, double average_weight_fragment, const std::set<unsigned int>& precursor_isotopes) { return self.estimateForFragmentFromRNAWeight(average_weight_precursor, average_weight_fragment, precursor_isotopes); }, "average_weight_precursor"_a, "average_weight_fragment"_a, "precursor_isotopes"_a, "Estimate fragment isotope distribution from RNA weights")
        .def("estimateForFragmentFromWeightAndComp", [](const OpenMS::CoarseIsotopePatternGenerator& self, double average_weight_precursor, double average_weight_fragment, const std::set<unsigned int>& precursor_isotopes, double C, double H, double N, double O, double S, double P) { return self.estimateForFragmentFromWeightAndComp(average_weight_precursor, average_weight_fragment, precursor_isotopes, C, H, N, O, S, P); }, "average_weight_precursor"_a, "average_weight_fragment"_a, "precursor_isotopes"_a, "C"_a, "H"_a, "N"_a, "O"_a, "S"_a, "P"_a, "Estimate fragment isotope distribution from weight and composition")
        ;

    // -----------------------------------------------------------------------
    // ConversionIssue
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::ConversionIssue>(m, "ConversionIssue", "Description of a conversion issue from Peptidoform to AASequence")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::ConversionIssue &>())
        .def("__copy__", [](const OpenMS::ProForma::ConversionIssue& self) { return OpenMS::ProForma::ConversionIssue(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::ConversionIssue& self, nb::dict) { return OpenMS::ProForma::ConversionIssue(self); }, "memo"_a)
        .def_rw("type", &OpenMS::ProForma::ConversionIssue::type)
        .def_rw("description", &OpenMS::ProForma::ConversionIssue::description)
        .def_rw("position", &OpenMS::ProForma::ConversionIssue::position)
        ;

    // -----------------------------------------------------------------------
    // CvAccession
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::CvAccession>(m, "CvAccession", "Controlled vocabulary accession for a modification (e.g., UNIMOD:35)")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::CvAccession &>())
        .def("__copy__", [](const OpenMS::ProForma::CvAccession& self) { return OpenMS::ProForma::CvAccession(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::CvAccession& self, nb::dict) { return OpenMS::ProForma::CvAccession(self); }, "memo"_a)
        .def_rw("database", &OpenMS::ProForma::CvAccession::database)
        .def_rw("accession", &OpenMS::ProForma::CvAccession::accession)
        ;

    // -----------------------------------------------------------------------
    // DecoyGenerator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DecoyGenerator>(m, "DecoyGenerator", 
        R"doc(
Methods to generate isobaric decoy sequences for DDA target-decoy
searches
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::DecoyGenerator& self) { return OpenMS::DecoyGenerator(self); })
        .def("__deepcopy__", [](const OpenMS::DecoyGenerator& self, nb::dict) { return OpenMS::DecoyGenerator(self); }, "memo"_a)
        .def("setSeed", [](OpenMS::DecoyGenerator& self, size_t seed) { return self.setSeed(seed); }, "seed"_a)
        .def("reverseProtein", [](const OpenMS::DecoyGenerator& self, const OpenMS::AASequence& protein) { return self.reverseProtein(protein); }, "protein"_a, "Reverses the protein sequence")
        .def("reversePeptides", [](const OpenMS::DecoyGenerator& self, const OpenMS::AASequence& protein, const OpenMS::String& protease) { return self.reversePeptides(protein, protease); }, "protein"_a, "protease"_a, "Reverses the protein's peptide sequences between enzymatic cutting positions")
        .def("shuffle", [](OpenMS::DecoyGenerator& self, const OpenMS::AASequence& protein, const OpenMS::String& protease, int decoy_factor) { return self.shuffle(protein, protease, decoy_factor); }, "protein"_a, "protease"_a, "decoy_factor"_a = 1, 
            R"doc(
Generate decoy protein sequences using shuffle algorithm. Digests protein using specified protease and shuffles each peptide. For top-down proteomics use "no cleavage". decoy_factor is the number of complete decoy proteins to generate. Returns vector of AASequence
)doc")
        .def("shufflePeptides", [](OpenMS::DecoyGenerator& self, const OpenMS::AASequence& aas, const OpenMS::String& protease, int max_attempts) { return self.shufflePeptides(aas, protease, max_attempts); }, "aas"_a, "protease"_a, "max_attempts"_a = 100, "Shuffle the protein's peptide sequences between enzymatic cutting positions. Each peptide is shuffled `max_attempts` times to minimize sequence identity")
        ;

    // -----------------------------------------------------------------------
    // DigestionEnzyme
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DigestionEnzyme>(m, "DigestionEnzyme", "Base class for digestion enzymes")
        .def(nb::init<const OpenMS::DigestionEnzyme &>())
        .def("__copy__", [](const OpenMS::DigestionEnzyme& self) { return OpenMS::DigestionEnzyme(self); })
        .def("__deepcopy__", [](const OpenMS::DigestionEnzyme& self, nb::dict) { return OpenMS::DigestionEnzyme(self); }, "memo"_a)
        .def(nb::init<OpenMS::String, OpenMS::String, std::set<OpenMS::String>, OpenMS::String>())
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::String, OpenMS::String, std::set<OpenMS::String>, OpenMS::String>())
        .def("setName", [](OpenMS::DigestionEnzyme& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the enzyme")
        .def("getName", [](const OpenMS::DigestionEnzyme& self) { return self.getName(); }, "Returns the name of the enzyme")
        .def("setSynonyms", [](OpenMS::DigestionEnzyme& self, const std::set<OpenMS::String>& synonyms) { return self.setSynonyms(synonyms); }, "synonyms"_a, "Sets the synonyms")
        .def("addSynonym", [](OpenMS::DigestionEnzyme& self, const OpenMS::String& synonym) { return self.addSynonym(synonym); }, "synonym"_a, "Adds a synonym")
        .def("getSynonyms", [](const OpenMS::DigestionEnzyme& self) -> const std::set<OpenMS::String> & { return self.getSynonyms(); }, nb::rv_policy::reference_internal, "Returns the synonyms")
        .def("setRegEx", [](OpenMS::DigestionEnzyme& self, const OpenMS::String& cleavage_regex) { return self.setRegEx(cleavage_regex); }, "cleavage_regex"_a, "Sets the cleavage regex")
        .def("getRegEx", [](const OpenMS::DigestionEnzyme& self) { return self.getRegEx(); }, "Returns the cleavage regex")
        .def("setRegExDescription", [](OpenMS::DigestionEnzyme& self, const OpenMS::String& value) { return self.setRegExDescription(value); }, "value"_a, "Sets the regex description")
        .def("getRegExDescription", [](const OpenMS::DigestionEnzyme& self) { return self.getRegExDescription(); }, "Returns the regex description")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self < nb::self)
        .def("setValueFromFile", [](OpenMS::DigestionEnzyme& self, const OpenMS::String& key, const OpenMS::String& value) { return self.setValueFromFile(key, value); }, "key"_a, "value"_a, "Sets the value of a member variable based on an entry from an input file")
        .def("__hash__", [](const OpenMS::DigestionEnzyme& self) { return std::hash<OpenMS::DigestionEnzyme>{}(self); })
        ;

    // -----------------------------------------------------------------------
    // DigestionEnzymeProtein
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DigestionEnzymeProtein, OpenMS::DigestionEnzyme>(m, "DigestionEnzymeProtein", 
        R"doc(
DigestionEnzyme

Representation of a digestion enzyme for proteins (protease)
)doc")
        .def(nb::init<>())
        .def(nb::init<OpenMS::DigestionEnzyme>())
        .def(nb::init<const OpenMS::DigestionEnzymeProtein &>())
        .def("__copy__", [](const OpenMS::DigestionEnzymeProtein& self) { return OpenMS::DigestionEnzymeProtein(self); })
        .def("__deepcopy__", [](const OpenMS::DigestionEnzymeProtein& self, nb::dict) { return OpenMS::DigestionEnzymeProtein(self); }, "memo"_a)
        .def(nb::init<OpenMS::String, OpenMS::String, std::set<OpenMS::String>, OpenMS::String, OpenMS::EmpiricalFormula, OpenMS::EmpiricalFormula, OpenMS::String, OpenMS::String, int, int, int>())
        .def("setNTermGain", [](OpenMS::DigestionEnzymeProtein& self, const OpenMS::EmpiricalFormula& value) { return self.setNTermGain(value); }, "value"_a, "Sets the N-term gain")
        .def("getNTermGain", [](const OpenMS::DigestionEnzymeProtein& self) { return self.getNTermGain(); }, "Returns the N-term gain")
        .def("setCTermGain", [](OpenMS::DigestionEnzymeProtein& self, const OpenMS::EmpiricalFormula& value) { return self.setCTermGain(value); }, "value"_a, "Sets the C-term gain")
        .def("getCTermGain", [](const OpenMS::DigestionEnzymeProtein& self) { return self.getCTermGain(); }, "Returns the C-term gain")
        .def("setPSIID", [](OpenMS::DigestionEnzymeProtein& self, const OpenMS::String& value) { return self.setPSIID(value); }, "value"_a, "Sets the PSI ID")
        .def("getPSIID", [](const OpenMS::DigestionEnzymeProtein& self) { return self.getPSIID(); }, "Returns the PSI ID")
        .def("setXTandemID", [](OpenMS::DigestionEnzymeProtein& self, const OpenMS::String& value) { return self.setXTandemID(value); }, "value"_a, "Sets the X! Tandem enzyme ID")
        .def("getXTandemID", [](const OpenMS::DigestionEnzymeProtein& self) { return self.getXTandemID(); }, "Returns the X! Tandem enzyme ID")
        .def("setCometID", [](OpenMS::DigestionEnzymeProtein& self, int value) { return self.setCometID(value); }, "value"_a, "Sets the Comet enzyme ID")
        .def("getCometID", [](const OpenMS::DigestionEnzymeProtein& self) { return self.getCometID(); }, "Returns the Comet enzyme ID")
        .def("setMSGFID", [](OpenMS::DigestionEnzymeProtein& self, int value) { return self.setMSGFID(value); }, "value"_a, "Sets the MSGFPlus enzyme id")
        .def("getMSGFID", [](const OpenMS::DigestionEnzymeProtein& self) { return self.getMSGFID(); }, "Returns the MSGFPlus enzyme id")
        .def("setOMSSAID", [](OpenMS::DigestionEnzymeProtein& self, int value) { return self.setOMSSAID(value); }, "value"_a, "Sets the OMSSA enzyme ID")
        .def("getOMSSAID", [](const OpenMS::DigestionEnzymeProtein& self) { return self.getOMSSAID(); }, "Returns the OMSSA enzyme ID")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self < nb::self)
        .def("setValueFromFile", [](OpenMS::DigestionEnzymeProtein& self, const OpenMS::String& key, const OpenMS::String& value) { return self.setValueFromFile(key, value); }, "key"_a, "value"_a, "Sets the value of a member variable based on an entry from an input file")
        .def("setName", [](OpenMS::DigestionEnzymeProtein& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the enzyme")
        .def("getName", [](const OpenMS::DigestionEnzymeProtein& self) { return self.getName(); }, "Returns the name of the enzyme")
        .def("setSynonyms", [](OpenMS::DigestionEnzymeProtein& self, const std::set<OpenMS::String>& synonyms) { return self.setSynonyms(synonyms); }, "synonyms"_a, "Sets the synonyms")
        .def("addSynonym", [](OpenMS::DigestionEnzymeProtein& self, const OpenMS::String& synonym) { return self.addSynonym(synonym); }, "synonym"_a, "Adds a synonym")
        .def("getSynonyms", [](const OpenMS::DigestionEnzymeProtein& self) -> const std::set<OpenMS::String> & { return self.getSynonyms(); }, nb::rv_policy::reference_internal, "Returns the synonyms")
        .def("setRegEx", [](OpenMS::DigestionEnzymeProtein& self, const OpenMS::String& cleavage_regex) { return self.setRegEx(cleavage_regex); }, "cleavage_regex"_a, "Sets the cleavage regex")
        .def("getRegEx", [](const OpenMS::DigestionEnzymeProtein& self) { return self.getRegEx(); }, "Returns the cleavage regex")
        .def("setRegExDescription", [](OpenMS::DigestionEnzymeProtein& self, const OpenMS::String& value) { return self.setRegExDescription(value); }, "value"_a, "Sets the regex description")
        .def("getRegExDescription", [](const OpenMS::DigestionEnzymeProtein& self) { return self.getRegExDescription(); }, "Returns the regex description")
        ;

    // -----------------------------------------------------------------------
    // DigestionEnzymeRNA
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DigestionEnzymeRNA, OpenMS::DigestionEnzyme>(m, "DigestionEnzymeRNA", 
        R"doc(
DigestionEnzyme

Representation of a digestion enzyme for RNA (RNase)
The cutting sites of these enzymes are defined using two different mechanisms:
First, a single regular expression that is applied to strings of unmodified RNA sequence and defines cutting sites via zero-length matches (using lookahead/lookbehind assertions).
This is the same mechanism that is used for proteases (see ProteaseDigestion).
However, due to the complex notation involved, this approach is not practical for modification-aware digestion.
Thus, the second mechanism uses two regular expressions ("cuts after"/"cuts before"), which are applied to the short codes (e.g. "m6A") of sequential ribonucleotides.
If both expressions match, then there is a cutting site between the two ribonucleotides.
There is support for terminal (5'/3') modifications that may be generated on fragments as a result of RNase cleavage.
A typical example is 3'-phosphate, resulting from cleavage of the phosphate backbone.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DigestionEnzymeRNA &>())
        .def("__copy__", [](const OpenMS::DigestionEnzymeRNA& self) { return OpenMS::DigestionEnzymeRNA(self); })
        .def("__deepcopy__", [](const OpenMS::DigestionEnzymeRNA& self, nb::dict) { return OpenMS::DigestionEnzymeRNA(self); }, "memo"_a)
        .def("setCutsAfterRegEx", [](OpenMS::DigestionEnzymeRNA& self, const OpenMS::String& value) { return self.setCutsAfterRegEx(value); }, "value"_a, 
            R"doc(
Sets the "cuts after ..." regular expression
)doc")
        .def("getCutsAfterRegEx", [](const OpenMS::DigestionEnzymeRNA& self) { return self.getCutsAfterRegEx(); }, 
            R"doc(
Returns the "cuts after ..." regular expression
)doc")
        .def("setCutsBeforeRegEx", [](OpenMS::DigestionEnzymeRNA& self, const OpenMS::String& value) { return self.setCutsBeforeRegEx(value); }, "value"_a, 
            R"doc(
Sets the "cuts before ..." regular expression
)doc")
        .def("getCutsBeforeRegEx", [](const OpenMS::DigestionEnzymeRNA& self) { return self.getCutsBeforeRegEx(); }, 
            R"doc(
Returns the "cuts before ..." regular expression
)doc")
        .def("setThreePrimeGain", [](OpenMS::DigestionEnzymeRNA& self, const OpenMS::String& value) { return self.setThreePrimeGain(value); }, "value"_a, "Sets the 3' gain")
        .def("getThreePrimeGain", [](const OpenMS::DigestionEnzymeRNA& self) { return self.getThreePrimeGain(); }, "Returns the 3' gain")
        .def("setFivePrimeGain", [](OpenMS::DigestionEnzymeRNA& self, const OpenMS::String& value) { return self.setFivePrimeGain(value); }, "value"_a, "Sets the 5' gain")
        .def("getFivePrimeGain", [](const OpenMS::DigestionEnzymeRNA& self) { return self.getFivePrimeGain(); }, "Returns the 5' gain")
        .def("setValueFromFile", [](OpenMS::DigestionEnzymeRNA& self, const OpenMS::String& key, const OpenMS::String& value) { return self.setValueFromFile(key, value); }, "key"_a, "value"_a, "Sets the value of a member variable based on an entry from an input file")
        .def("setName", [](OpenMS::DigestionEnzymeRNA& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the enzyme")
        .def("getName", [](const OpenMS::DigestionEnzymeRNA& self) { return self.getName(); }, "Returns the name of the enzyme")
        .def("setSynonyms", [](OpenMS::DigestionEnzymeRNA& self, const std::set<OpenMS::String>& synonyms) { return self.setSynonyms(synonyms); }, "synonyms"_a, "Sets the synonyms")
        .def("addSynonym", [](OpenMS::DigestionEnzymeRNA& self, const OpenMS::String& synonym) { return self.addSynonym(synonym); }, "synonym"_a, "Adds a synonym")
        .def("getSynonyms", [](const OpenMS::DigestionEnzymeRNA& self) -> const std::set<OpenMS::String> & { return self.getSynonyms(); }, nb::rv_policy::reference_internal, "Returns the synonyms")
        .def("setRegEx", [](OpenMS::DigestionEnzymeRNA& self, const OpenMS::String& cleavage_regex) { return self.setRegEx(cleavage_regex); }, "cleavage_regex"_a, "Sets the cleavage regex")
        .def("getRegEx", [](const OpenMS::DigestionEnzymeRNA& self) { return self.getRegEx(); }, "Returns the cleavage regex")
        .def("setRegExDescription", [](OpenMS::DigestionEnzymeRNA& self, const OpenMS::String& value) { return self.setRegExDescription(value); }, "value"_a, "Sets the regex description")
        .def("getRegExDescription", [](const OpenMS::DigestionEnzymeRNA& self) { return self.getRegExDescription(); }, "Returns the regex description")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self < nb::self)
        ;

    // -----------------------------------------------------------------------
    // Element
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Element>(m, "Element", "Representation of an element")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Element &>())
        .def("__copy__", [](const OpenMS::Element& self) { return OpenMS::Element(self); })
        .def("__deepcopy__", [](const OpenMS::Element& self, nb::dict) { return OpenMS::Element(self); }, "memo"_a)
        .def(nb::init<OpenMS::String, OpenMS::String, unsigned int, double, double, OpenMS::IsotopeDistribution>())
        .def("setAtomicNumber", [](OpenMS::Element& self, unsigned int atomic_number) { return self.setAtomicNumber(atomic_number); }, "atomic_number"_a, "Sets unique atomic number")
        .def("getAtomicNumber", [](const OpenMS::Element& self) { return self.getAtomicNumber(); }, "Returns the unique atomic number")
        .def("setAverageWeight", [](OpenMS::Element& self, double weight) { return self.setAverageWeight(weight); }, "weight"_a, "Sets the average weight of the element")
        .def("getAverageWeight", [](const OpenMS::Element& self) { return self.getAverageWeight(); }, "Returns the average weight of the element")
        .def("setMonoWeight", [](OpenMS::Element& self, double weight) { return self.setMonoWeight(weight); }, "weight"_a, "Sets the mono isotopic weight of the element")
        .def("getMonoWeight", [](const OpenMS::Element& self) { return self.getMonoWeight(); }, "Returns the mono isotopic weight of the element")
        .def("setIsotopeDistribution", [](OpenMS::Element& self, const OpenMS::IsotopeDistribution& isotopes) { return self.setIsotopeDistribution(isotopes); }, "isotopes"_a, "Sets the isotope distribution of the element")
        .def("getIsotopeDistribution", [](const OpenMS::Element& self) -> const OpenMS::IsotopeDistribution & { return self.getIsotopeDistribution(); }, nb::rv_policy::reference_internal, "Returns the isotope distribution of the element")
        .def("setName", [](OpenMS::Element& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the element")
        .def("getName", [](const OpenMS::Element& self) { return self.getName(); }, "Returns the name of the element")
        .def("setSymbol", [](OpenMS::Element& self, const OpenMS::String& symbol) { return self.setSymbol(symbol); }, "symbol"_a, "Sets symbol of the element")
        .def("getSymbol", [](const OpenMS::Element& self) { return self.getSymbol(); }, "Returns symbol of the element")
        .def("__hash__", [](const OpenMS::Element& self) { return std::hash<OpenMS::Element>{}(self); })
        ;

    // -----------------------------------------------------------------------
    // ElementDB
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ElementDB>(m, "ElementDB", 
        R"doc(
Database of chemical elements and their isotopes. This is a singleton class.
The elements are initialized with data from IUPAC tables.
)doc")

        .def("getElement", [](const OpenMS::ElementDB& self, const std::string& name) -> const OpenMS::Element* {
            return self.getElement(name);
        }, "name"_a, nb::rv_policy::reference, "Get element by name or symbol")

        .def("getElement", [](const OpenMS::ElementDB& self, unsigned int atomic_number) -> const OpenMS::Element* {
            return self.getElement(atomic_number);
        }, "atomic_number"_a, nb::rv_policy::reference, "Get element by atomic number")

        .def("hasElement", [](const OpenMS::ElementDB& self, const std::string& name) {
            return self.hasElement(name);
        }, "name"_a, "Check if element exists by name or symbol")

        .def("hasElement", [](const OpenMS::ElementDB& self, unsigned int atomic_number) {
            return self.hasElement(atomic_number);
        }, "atomic_number"_a, "Check if element exists by atomic number")

        .def("addElement", [](OpenMS::ElementDB& self, const std::string& name, const std::string& symbol,
                unsigned int an, const std::map<unsigned int, double>& abundance,
                const std::map<unsigned int, double>& mass, bool replace_existing) {
            self.addElement(name, symbol, an, abundance, mass, replace_existing);
        }, "name"_a, "symbol"_a, "atomic_number"_a, "abundance"_a, "mass"_a, "replace_existing"_a,
        "Add a new element to the database")
        .def_static("getInstance", []() -> OpenMS::ElementDB* { return OpenMS::ElementDB::getInstance(); }, nb::rv_policy::reference, "Returns the singleton instance")
        ;

    // -----------------------------------------------------------------------
    // EmpiricalFormula
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::EmpiricalFormula>(m, "EmpiricalFormula", "Representation of an empirical formula")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::EmpiricalFormula &>())
        .def("__copy__", [](const OpenMS::EmpiricalFormula& self) { return OpenMS::EmpiricalFormula(self); })
        .def("__deepcopy__", [](const OpenMS::EmpiricalFormula& self, nb::dict) { return OpenMS::EmpiricalFormula(self); }, "memo"_a)
        .def(nb::init<OpenMS::String>())
        .def("getMonoWeight", [](const OpenMS::EmpiricalFormula& self) { return self.getMonoWeight(); }, "Returns the mono isotopic weight of the formula (includes proton charges)")
        .def("getAverageWeight", [](const OpenMS::EmpiricalFormula& self) { return self.getAverageWeight(); }, "Returns the average weight of the formula (includes proton charges)")
        .def("getLightestIsotopeWeight", [](const OpenMS::EmpiricalFormula& self) { return self.getLightestIsotopeWeight(); }, "Returns the sum of the lightest isotope weight of all elements in the formula (includes proton charges)")
        .def("calculateTheoreticalIsotopesNumber", [](const OpenMS::EmpiricalFormula& self) { return self.calculateTheoreticalIsotopesNumber(); })
        .def("estimateFromWeightAndComp", [](OpenMS::EmpiricalFormula& self, double average_weight, double C, double H, double N, double O, double S, double P) { return self.estimateFromWeightAndComp(average_weight, C, H, N, O, S, P); }, "average_weight"_a, "C"_a, "H"_a, "N"_a, "O"_a, "S"_a, "P"_a, "Fills this EmpiricalFormula with an approximate elemental composition for a given average weight and approximate elemental stoichiometry")
        .def("estimateFromWeightAndCompAndS", [](OpenMS::EmpiricalFormula& self, double average_weight, unsigned int S, double C, double H, double N, double O, double P) { return self.estimateFromWeightAndCompAndS(average_weight, S, C, H, N, O, P); }, "average_weight"_a, "S"_a, "C"_a, "H"_a, "N"_a, "O"_a, "P"_a, "Fills this EmpiricalFormula with an approximate elemental composition for a given average weight, exact number of sulfurs, and approximate elemental stoichiometry")
        .def("getNumberOfAtoms", [](const OpenMS::EmpiricalFormula& self) { return self.getNumberOfAtoms(); }, "Returns the total number of atoms")
        .def("getCharge", [](const OpenMS::EmpiricalFormula& self) { return self.getCharge(); }, "Returns the total charge")
        .def("setCharge", [](OpenMS::EmpiricalFormula& self, int charge) { return self.setCharge(charge); }, "charge"_a, "Sets the charge")
        .def("toString", [](const OpenMS::EmpiricalFormula& self) { return self.toString(); }, "Returns the formula as a string (charges are not included)")
        .def("getElementalComposition", [](const OpenMS::EmpiricalFormula& self) { return self.toMap(); }, "Get elemental composition as a hash {'Symbol' -> NrAtoms}")
        .def(nb::self + nb::self)
        .def(nb::self - nb::self)
        .def("isEmpty", [](const OpenMS::EmpiricalFormula& self) { return self.isEmpty(); }, "Returns true if the formula does not contain a element")
        .def("isCharged", [](const OpenMS::EmpiricalFormula& self) { return self.isCharged(); }, "Returns true if charge is not equal to zero")
        .def("hasElement", [](const OpenMS::EmpiricalFormula& self, OpenMS::Element * element) { return self.hasElement(element); }, "element"_a, "Returns true if the formula contains the element")
        .def("getNumberOf", [](const OpenMS::EmpiricalFormula& self, OpenMS::Element * element) { return self.getNumberOf(element); }, "element"_a, "Returns the number of atoms of the given element (can be negative)")
        .def("contains", [](const OpenMS::EmpiricalFormula& self, const OpenMS::EmpiricalFormula& ef) { return self.contains(ef); }, "ef"_a, "Returns true if all elements from `ef` ( empirical formula ) are LESS abundant (negative allowed) than the corresponding elements of this EmpiricalFormula")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("__hash__", [](const OpenMS::EmpiricalFormula& self) { return std::hash<OpenMS::EmpiricalFormula>{}(self); })

        .def("getElementalComposition", [](const OpenMS::EmpiricalFormula& self) {
            return self.toMap();
        }, "Get elemental composition as a dict {'Symbol': NrAtoms}")

        .def("getIsotopeDistribution", [](const OpenMS::EmpiricalFormula& self, const OpenMS::CoarseIsotopePatternGenerator& method) {
            return self.getIsotopeDistribution(method);
        }, "method"_a, "Returns isotope distribution using CoarseIsotopePatternGenerator")
        .def("getIsotopeDistribution", [](const OpenMS::EmpiricalFormula& self, const OpenMS::FineIsotopePatternGenerator& method) {
            return self.getIsotopeDistribution(method);
        }, "method"_a, "Returns isotope distribution using FineIsotopePatternGenerator")

        .def("__repr__", [](const OpenMS::EmpiricalFormula& self) {
            std::ostringstream oss;
            oss << "EmpiricalFormula(formula='" << std::string(self.toString())
                << "', mono_mass=" << self.getMonoWeight() << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::EmpiricalFormula& self) {
            return std::string(self.toString());
        }, "Returns the formula as a string")
        .def("__iadd__", [](OpenMS::EmpiricalFormula& self, const OpenMS::EmpiricalFormula& rhs) -> OpenMS::EmpiricalFormula& { return self += rhs; }, "rhs"_a, "In-place addition of empirical formulas", nb::rv_policy::reference_internal)
        .def("getConditionalFragmentIsotopeDist", [](const OpenMS::EmpiricalFormula& self, const OpenMS::EmpiricalFormula& precursor, const std::set<unsigned int>& precursor_isotopes, const OpenMS::CoarseIsotopePatternGenerator& method) { return self.getConditionalFragmentIsotopeDist(precursor, precursor_isotopes, method); }, "precursor"_a, "precursor_isotopes"_a, "method"_a, "Returns the conditional fragment isotope distribution")
        .def("__isub__", [](OpenMS::EmpiricalFormula& self, const OpenMS::EmpiricalFormula& rhs) -> OpenMS::EmpiricalFormula& { return self -= rhs; }, "rhs"_a, "In-place subtraction of empirical formulas", nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // EnzymaticDigestion
    // -----------------------------------------------------------------------
    auto enzymaticdigestion_class = nb::class_<OpenMS::EnzymaticDigestion>(m, "EnzymaticDigestion", 
        R"doc(
Class for the enzymatic digestion of proteins
Digestion can be performed using simple regular expressions, e.g. [KR] | [^P] for trypsin.
Also missed cleavages can be modeled, i.e. adjacent peptides are not cleaved
due to enzyme malfunction/access restrictions. If n missed cleavages are allowed, all possible resulting
peptides (cleaved and uncleaved) with up to n missed cleavages are returned.
Thus no random selection of just n specific missed cleavage sites is performed.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::EnzymaticDigestion &>())
        .def("__copy__", [](const OpenMS::EnzymaticDigestion& self) { return OpenMS::EnzymaticDigestion(self); })
        .def("__deepcopy__", [](const OpenMS::EnzymaticDigestion& self, nb::dict) { return OpenMS::EnzymaticDigestion(self); }, "memo"_a)
        .def("getMissedCleavages", [](const OpenMS::EnzymaticDigestion& self) { return self.getMissedCleavages(); }, "Returns the max. number of allowed missed cleavages for the digestion")
        .def("setMissedCleavages", [](OpenMS::EnzymaticDigestion& self, size_t missed_cleavages) { return self.setMissedCleavages(missed_cleavages); }, "missed_cleavages"_a, "Sets the max. number of allowed missed cleavages for the digestion (default is 0). This setting is ignored when log model is used")
        .def("getEnzymeName", [](const OpenMS::EnzymaticDigestion& self) { return self.getEnzymeName(); }, "Returns the enzyme for the digestion")
        .def("setEnzyme", [](OpenMS::EnzymaticDigestion& self, OpenMS::DigestionEnzyme * enzyme) { return self.setEnzyme(enzyme); }, "enzyme"_a, "Sets the enzyme for the digestion")
        .def("getSpecificity", [](const OpenMS::EnzymaticDigestion& self) { return self.getSpecificity(); }, "Returns the specificity for the digestion")
        .def("setSpecificity", [](OpenMS::EnzymaticDigestion& self, OpenMS::EnzymaticDigestion::Specificity spec) { return self.setSpecificity(spec); }, "spec"_a, "Sets the specificity for the digestion (default is SPEC_FULL)")
        .def_static("getSpecificityByName", [](const OpenMS::String& name) { return OpenMS::EnzymaticDigestion::getSpecificityByName(name); }, "name"_a, "Returns the specificity by name. Returns SPEC_UNKNOWN if name is not valid")
        .def("isValidProduct", [](const OpenMS::EnzymaticDigestion& self, const OpenMS::String& protein, int pep_pos, int pep_length, bool ignore_missed_cleavages) { return self.isValidProduct(protein, pep_pos, pep_length, ignore_missed_cleavages); }, "protein"_a, "pep_pos"_a, "pep_length"_a, "ignore_missed_cleavages"_a = true, 
            R"doc(
Performs the enzymatic digestion of an unmodified sequence\n
By returning only references into the original string this is very fast
:param sequence: Sequence to digest
:param output: Digestion products
:param min_length: Minimal length of reported products
:param max_length: Maximal length of reported products (0 = no restriction)
:return: Number of discarded digestion products (which are not matching length restrictions)
)doc")
        .def("countInternalCleavageSites", [](const OpenMS::EnzymaticDigestion& self, const OpenMS::String& sequence) { return self.countInternalCleavageSites(sequence); }, "sequence"_a, "Returns the number of internal cleavage sites for this sequence.")
        .def("digestUnmodified", [](const OpenMS::EnzymaticDigestion& self, const std::string& sequence_str, size_t min_length, size_t max_length) {
            OpenMS::StringView sequence(sequence_str);
            std::vector<OpenMS::StringView> output;
            OpenMS::Size discarded = self.digestUnmodified(sequence, output, min_length, max_length);
            // Convert StringViews to strings for Python
            std::vector<std::string> result;
            result.reserve(output.size());
            for (const auto& sv : output) result.push_back(std::string(sv.getString()));
            return nb::make_tuple(result, discarded);
        }, "sequence"_a, "min_length"_a = 1, "max_length"_a = 0, "Digest unmodified sequence, returns (products, num_discarded)")
        ;
    // Specificity enum nested under EnzymaticDigestion
    nb::enum_<OpenMS::EnzymaticDigestion::Specificity>(enzymaticdigestion_class, "Specificity", nb::is_arithmetic())
        .value("SPEC_NONE", OpenMS::EnzymaticDigestion::Specificity::SPEC_NONE)
        .value("SPEC_SEMI", OpenMS::EnzymaticDigestion::Specificity::SPEC_SEMI)
        .value("SPEC_FULL", OpenMS::EnzymaticDigestion::Specificity::SPEC_FULL)
        .value("SPEC_UNKNOWN", OpenMS::EnzymaticDigestion::Specificity::SPEC_UNKNOWN)
        .value("SPEC_NOCTERM", OpenMS::EnzymaticDigestion::Specificity::SPEC_NOCTERM)
        .value("SPEC_NONTERM", OpenMS::EnzymaticDigestion::Specificity::SPEC_NONTERM)
        .value("SIZE_OF_SPECIFICITY", OpenMS::EnzymaticDigestion::Specificity::SIZE_OF_SPECIFICITY)
        .export_values();

    // -----------------------------------------------------------------------
    // FineIsotopePatternGenerator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FineIsotopePatternGenerator>(m, "FineIsotopePatternGenerator", 
        R"doc(
Isotope pattern generator for fine isotope distributions.
Generates isotopes until a stop condition (threshold) is reached,
the lower the threshold the more isotopes are generated. The
parameter use_total_prob defines whether the stop condition is
interpreted as the total probability that the distribution should
cover (default) or as a threshold for individual peaks. Finally,
the absolute parameter specifies for individual peak thresholding
if the threshold is absolute or relative.
)doc")
        .def(nb::init<>())
        .def(nb::init<double, bool, bool>())
        .def("run", [](const OpenMS::FineIsotopePatternGenerator& self, const OpenMS::EmpiricalFormula& ef) { return self.run(ef); }, "ef"_a)
        .def("setThreshold", [](OpenMS::FineIsotopePatternGenerator& self, double stop_condition) { return self.setThreshold(stop_condition); }, "stop_condition"_a)
        .def("getThreshold", [](const OpenMS::FineIsotopePatternGenerator& self) { return self.getThreshold(); })
        .def("setAbsolute", [](OpenMS::FineIsotopePatternGenerator& self, bool absolute) { return self.setAbsolute(absolute); }, "absolute"_a)
        .def("getAbsolute", [](const OpenMS::FineIsotopePatternGenerator& self) { return self.getAbsolute(); })
        .def("setTotalProbability", [](OpenMS::FineIsotopePatternGenerator& self, bool total) { return self.setTotalProbability(total); }, "total"_a)
        .def("getTotalProbability", [](const OpenMS::FineIsotopePatternGenerator& self) { return self.getTotalProbability(); })
        ;

    // -----------------------------------------------------------------------
    // FormulaTag
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::FormulaTag>(m, "FormulaTag", "Chemical formula modification tag")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::FormulaTag &>())
        .def("__copy__", [](const OpenMS::ProForma::FormulaTag& self) { return OpenMS::ProForma::FormulaTag(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::FormulaTag& self, nb::dict) { return OpenMS::ProForma::FormulaTag(self); }, "memo"_a)
        .def_rw("formula_string", &OpenMS::ProForma::FormulaTag::formula_string)
        ;

    // -----------------------------------------------------------------------
    // IMSElement
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ims::IMSElement>(m, "IMSElement", "OpenMS class IMSElement")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ims::IMSElement &>())
        .def("__copy__", [](const OpenMS::ims::IMSElement& self) { return OpenMS::ims::IMSElement(self); })
        .def("__deepcopy__", [](const OpenMS::ims::IMSElement& self, nb::dict) { return OpenMS::ims::IMSElement(self); }, "memo"_a)
        .def(nb::init<OpenMS::String, OpenMS::ims::IMSIsotopeDistribution>())
        .def(nb::init<OpenMS::String, double>())
        .def(nb::init<OpenMS::String, unsigned int>())
        .def("getName", [](const OpenMS::ims::IMSElement& self) { return self.getName(); }, "Gets element's name")
        .def("setName", [](OpenMS::ims::IMSElement& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets element's name")
        .def("getSequence", [](const OpenMS::ims::IMSElement& self) { return self.getSequence(); }, "Gets element's sequence")
        .def("setSequence", [](OpenMS::ims::IMSElement& self, const OpenMS::String& sequence) { return self.setSequence(sequence); }, "sequence"_a, "Sets element's sequence")
        .def("getNominalMass", [](const OpenMS::ims::IMSElement& self) { return self.getNominalMass(); }, "Gets element's nominal mass")
        .def("getMass", [](const OpenMS::ims::IMSElement& self, size_t index) { return self.getMass(index); }, "index"_a = 0, "Gets mass of element's isotope 'index'")
        .def("getAverageMass", [](const OpenMS::ims::IMSElement& self) { return self.getAverageMass(); }, "Gets element's average mass")
        .def("getIonMass", [](const OpenMS::ims::IMSElement& self, int electrons_number) { return self.getIonMass(electrons_number); }, "electrons_number"_a = 1, "Gets ion mass of element. By default ion lacks 1 electron, but this can be changed by setting other 'electrons_number'")
        .def("getIsotopeDistribution", [](const OpenMS::ims::IMSElement& self) -> const OpenMS::ims::IMSIsotopeDistribution & { return self.getIsotopeDistribution(); }, nb::rv_policy::reference_internal, "Gets element's isotope distribution")
        .def("setIsotopeDistribution", [](OpenMS::ims::IMSElement& self, const OpenMS::ims::IMSIsotopeDistribution& isotopes) { return self.setIsotopeDistribution(isotopes); }, "isotopes"_a, "Sets element's isotope distribution")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        ;

    // -----------------------------------------------------------------------
    // IMSIsotopeDistribution
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ims::IMSIsotopeDistribution>(m, "IMSIsotopeDistribution", 
        R"doc(
Represents a distribution of isotopes restricted to the first K elements
Represents a distribution of isotopes of chemical elements as a list
of peaks each as a pair of mass and abundance. 'IsotopeDistribution'
unlike 'IsotopeSpecies' has one abundance per a nominal mass.
Here is an example in the format (mass; abundance %)
for molecule H2O (values are taken randomly):
- IsotopeDistribution
(18.00221; 99.03 %)
(19.00334; 0.8 %)
(20.00476; 0.17 %)
- IsotopeSpecies
(18.00197; 98.012 %)
(18.00989; 1.018 %)
(19.00312; 0.683 %)
(19.00531; 0.117 %)
(20.00413; 0.134 %)
(20.00831; 0.036 %)
To the sake of faster computations distribution is restricted
to the first K elements, where K can be set by adjusting size
'SIZE' of distribution. @note For the elements most abundant in
living beings (CHNOPS) this restriction is negligible, since abundances
decrease dramatically in isotopes order and are usually of no interest
starting from +10 isotope.
'IsotopeDistribution' implements folding with other distribution using an
algorithm described in details in paper:
Boecker et al. "Decomposing metabolic isotope patterns" WABI 2006. doi: 10.1007/11851561_2
Folding with itself is done using Russian Multiplication Scheme
)doc")
        .def(nb::init<unsigned int>())
        .def(nb::init<double>())
        .def(nb::init<std::vector<OpenMS::ims::IMSIsotopeDistribution::Peak>, unsigned int>())
        .def(nb::init<const OpenMS::ims::IMSIsotopeDistribution &>())
        .def("__copy__", [](const OpenMS::ims::IMSIsotopeDistribution& self) { return OpenMS::ims::IMSIsotopeDistribution(self); })
        .def("__deepcopy__", [](const OpenMS::ims::IMSIsotopeDistribution& self, nb::dict) { return OpenMS::ims::IMSIsotopeDistribution(self); }, "memo"_a)
        .def("size", [](const OpenMS::ims::IMSIsotopeDistribution& self) { return self.size(); })
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getMass", [](const OpenMS::ims::IMSIsotopeDistribution& self, size_t i) { return self.getMass(i); }, "i"_a)
        .def("getAbundance", [](const OpenMS::ims::IMSIsotopeDistribution& self, size_t i) { return self.getAbundance(i); }, "i"_a)
        .def("getAverageMass", [](const OpenMS::ims::IMSIsotopeDistribution& self) { return self.getAverageMass(); })
        .def("getNominalMass", [](const OpenMS::ims::IMSIsotopeDistribution& self) { return self.getNominalMass(); })
        .def("setNominalMass", [](OpenMS::ims::IMSIsotopeDistribution& self, unsigned int nominalMass) { return self.setNominalMass(nominalMass); }, "nominalMass"_a)
        .def("getMasses", [](const OpenMS::ims::IMSIsotopeDistribution& self) { return self.getMasses(); }, "Gets a mass of isotope 'i'")
        .def("getAbundances", [](const OpenMS::ims::IMSIsotopeDistribution& self) { return self.getAbundances(); }, "Gets an abundance of isotope 'i'")
        .def("normalize", [](OpenMS::ims::IMSIsotopeDistribution& self) { return self.normalize(); }, "Normalizes distribution, i.e. scaling abundances to be summed up to 1 with an error")
        .def("empty", [](const OpenMS::ims::IMSIsotopeDistribution& self) { return self.empty(); }, "Returns true if the distribution has no peaks, false - otherwise")
        .def("__len__", [](OpenMS::ims::IMSIsotopeDistribution& self) { return self.size(); })
        .def_ro_static("ABUNDANCES_SUM_ERROR", &OpenMS::ims::IMSIsotopeDistribution::ABUNDANCES_SUM_ERROR)
        .def_ro_static("SIZE", &OpenMS::ims::IMSIsotopeDistribution::SIZE)
        ;

    // -----------------------------------------------------------------------
    // IMSIsotopeDistribution_Peak
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ims::IMSIsotopeDistribution::Peak>(m, "IMSIsotopeDistribution_Peak", "OpenMS class IMSIsotopeDistribution_Peak")
        .def(nb::init<double, double>())
        .def(nb::self == nb::self)
        .def_rw("mass", &OpenMS::ims::IMSIsotopeDistribution::Peak::mass)
        .def_rw("abundance", &OpenMS::ims::IMSIsotopeDistribution::Peak::abundance)
        ;

    // -----------------------------------------------------------------------
    // InfoTag
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::InfoTag>(m, "InfoTag", "Info tag for arbitrary text annotations")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::InfoTag &>())
        .def("__copy__", [](const OpenMS::ProForma::InfoTag& self) { return OpenMS::ProForma::InfoTag(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::InfoTag& self, nb::dict) { return OpenMS::ProForma::InfoTag(self); }, "memo"_a)
        .def_rw("text", &OpenMS::ProForma::InfoTag::text)
        ;

    // -----------------------------------------------------------------------
    // IsotopeDistribution
    // -----------------------------------------------------------------------
    auto isotopedistribution_class = nb::class_<OpenMS::IsotopeDistribution>(m, "IsotopeDistribution", 
        R"doc(
Isotope distribution class
A container that holds an isotope distribution. It consists of mass values
and their correspondent probabilities (stored in the intensity slot)
Isotope distributions can be calculated using either the
CoarseIsotopePatternGenerator for quantized atomic masses which group
isotopes with the same atomic number. Alternatively, the
FineIsotopePatternGenerator can be used that calculates hyperfine isotopic
distributions
This class only describes the container that holds the isotopic
distribution, calculations are done using classes derived from
IsotopePatternGenerator
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IsotopeDistribution &>())
        .def("__copy__", [](const OpenMS::IsotopeDistribution& self) { return OpenMS::IsotopeDistribution(self); })
        .def("__deepcopy__", [](const OpenMS::IsotopeDistribution& self, nb::dict) { return OpenMS::IsotopeDistribution(self); }, "memo"_a)
        .def("set", [](OpenMS::IsotopeDistribution& self, const std::vector<OpenMS::Peak1D>& distribution) { return self.set(distribution); }, "distribution"_a, "Overwrites the container which holds the distribution using 'distribution'")
        .def("set", [](OpenMS::IsotopeDistribution& self, std::vector<OpenMS::Peak1D>& distribution) { return self.set(distribution); }, "distribution"_a, "Overwrites the container which holds the distribution using 'distribution'")
        .def("getContainer", [](const OpenMS::IsotopeDistribution& self) -> const std::vector<OpenMS::Peak1D> & { return self.getContainer(); }, nb::rv_policy::reference_internal, "Returns the container which holds the distribution")
        .def("getMax", [](const OpenMS::IsotopeDistribution& self) { return self.getMax(); }, "Returns the maximal weight isotope which is stored in the distribution")
        .def("getMin", [](const OpenMS::IsotopeDistribution& self) { return self.getMin(); }, "Returns the minimal weight isotope which is stored in the distribution")
        .def("getMostAbundant", [](const OpenMS::IsotopeDistribution& self) { return self.getMostAbundant(); }, "Returns the most abundant isotope which is stored in the distribution")
        .def("size", [](const OpenMS::IsotopeDistribution& self) { return self.size(); }, "Returns the size of the distribution which is the number of isotopes in the distribution")
        .def("clear", [](OpenMS::IsotopeDistribution& self) { return self.clear(); }, "Clears the distribution and resets max isotope to 0")
        .def("resize", [](OpenMS::IsotopeDistribution& self, unsigned int size) { return self.resize(size); }, "size"_a, "Resizes distribution container")
        .def("trimIntensities", [](OpenMS::IsotopeDistribution& self, double cutoff) { return self.trimIntensities(cutoff); }, "cutoff"_a, "Remove intensities below the cutoff")
        .def("sortByIntensity", [](OpenMS::IsotopeDistribution& self) { return self.sortByIntensity(); }, "Sort isotope distribution by intensity")
        .def("sortByMass", [](OpenMS::IsotopeDistribution& self) { return self.sortByMass(); }, "Sort isotope distribution by mass")
        .def("renormalize", [](OpenMS::IsotopeDistribution& self) { return self.renormalize(); }, "Renormalizes the sum of the probabilities of the isotopes to 1")
        .def("merge", [](OpenMS::IsotopeDistribution& self, double resolution, double min_prob) { return self.merge(resolution, min_prob); }, "resolution"_a, "min_prob"_a, "Merges distributions of arbitrary data points with constant defined resolution")
        .def("trimRight", [](OpenMS::IsotopeDistribution& self, double cutoff) { return self.trimRight(cutoff); }, "cutoff"_a, "Trims the right side of the isotope distribution to isotopes with a significant contribution")
        .def("trimLeft", [](OpenMS::IsotopeDistribution& self, double cutoff) { return self.trimLeft(cutoff); }, "cutoff"_a, "Trims the left side of the isotope distribution to isotopes with a significant contribution")
        .def("averageMass", [](const OpenMS::IsotopeDistribution& self) { return self.averageMass(); }, "Compute average mass of isotope distribution (weighted average of all isotopes)")
        .def("begin", [](const OpenMS::IsotopeDistribution& self) { return self.begin(); })
        .def("end", [](const OpenMS::IsotopeDistribution& self) { return self.end(); })
        .def("insert", [](OpenMS::IsotopeDistribution& self, const double& mass, const float& intensity) { return self.insert(mass, intensity); }, "mass"_a, "intensity"_a)
        .def("__hash__", [](const OpenMS::IsotopeDistribution& self) { return std::hash<OpenMS::IsotopeDistribution>{}(self); })
        .def("__iter__", [](OpenMS::IsotopeDistribution& self) { return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::IsotopeDistribution>(), "IsotopeDistribution_iter", self.begin(), self.end()); })
        .def("__len__", [](OpenMS::IsotopeDistribution& self) { return self.size(); })
        .def("__getitem__", [](OpenMS::IsotopeDistribution& self, size_t i) -> OpenMS::Peak1D & {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)
        .def("__repr__", [](const OpenMS::IsotopeDistribution& self) {
            std::ostringstream oss;
            oss << "IsotopeDistribution(num_isotopes=" << self.size()
                << ", mass_range=(" << self.getMin() << ", " << self.getMax() << "))";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::IsotopeDistribution& self) { return nb::cast(self).attr("__repr__")(); })
        ;
    // Sorted enum nested under IsotopeDistribution
    nb::enum_<OpenMS::IsotopeDistribution::Sorted>(isotopedistribution_class, "Sorted", nb::is_arithmetic())
        .value("INTENSITY", OpenMS::IsotopeDistribution::Sorted::INTENSITY)
        .value("MASS", OpenMS::IsotopeDistribution::Sorted::MASS)
        .value("UNDEFINED", OpenMS::IsotopeDistribution::Sorted::UNDEFINED)
        .export_values();

    // -----------------------------------------------------------------------
    // IsotopeReplacement
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::IsotopeReplacement>(m, "IsotopeReplacement", "Isotope replacement for stable isotope labeling")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::IsotopeReplacement &>())
        .def("__copy__", [](const OpenMS::ProForma::IsotopeReplacement& self) { return OpenMS::ProForma::IsotopeReplacement(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::IsotopeReplacement& self, nb::dict) { return OpenMS::ProForma::IsotopeReplacement(self); }, "memo"_a)
        .def_rw("isotope", &OpenMS::ProForma::IsotopeReplacement::isotope)
        ;

    // -----------------------------------------------------------------------
    // Label
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::Label>(m, "Label", "Label for cross-links, branches, or ambiguous grouping")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::Label &>())
        .def("__copy__", [](const OpenMS::ProForma::Label& self) { return OpenMS::ProForma::Label(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::Label& self, nb::dict) { return OpenMS::ProForma::Label(self); }, "memo"_a)
        .def_rw("identifier", &OpenMS::ProForma::Label::identifier)
        ;

    // -----------------------------------------------------------------------
    // LabileModification
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::LabileModification>(m, "LabileModification", "Labile modification that may be lost during fragmentation")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::LabileModification &>())
        .def("__copy__", [](const OpenMS::ProForma::LabileModification& self) { return OpenMS::ProForma::LabileModification(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::LabileModification& self, nb::dict) { return OpenMS::ProForma::LabileModification(self); }, "memo"_a)
        .def_rw("modification", &OpenMS::ProForma::LabileModification::modification)
        ;

    // -----------------------------------------------------------------------
    // MassDecomposition
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MassDecomposition>(m, "MassDecomposition", 
        R"doc(
Class represents a decomposition of a mass into amino acids
This class represents a mass decomposition into amino acids. A
decomposition are amino acids given with frequencies which add
up to a specific mass.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MassDecomposition &>())
        .def("__copy__", [](const OpenMS::MassDecomposition& self) { return OpenMS::MassDecomposition(self); })
        .def("__deepcopy__", [](const OpenMS::MassDecomposition& self, nb::dict) { return OpenMS::MassDecomposition(self); }, "memo"_a)
        .def(nb::init<OpenMS::String>())
        .def("toString", [](const OpenMS::MassDecomposition& self) { return self.toString(); }, "Returns the decomposition as a string")
        .def("toExpandedString", [](const OpenMS::MassDecomposition& self) { return self.toExpandedString(); }, "Returns the decomposition as a string; instead of frequencies the amino acids are repeated")
        .def("getNumberOfMaxAA", [](const OpenMS::MassDecomposition& self) { return self.getNumberOfMaxAA(); }, "Returns the max frequency of this composition")
        .def("containsTag", [](const OpenMS::MassDecomposition& self, const OpenMS::String& tag) { return self.containsTag(tag); }, "tag"_a, "Returns true if tag is contained in the mass decomposition")
        .def("compatible", [](const OpenMS::MassDecomposition& self, const OpenMS::MassDecomposition& deco) { return self.compatible(deco); }, "deco"_a, "Returns true if the mass decomposition if contained in this instance")
        ;

    // -----------------------------------------------------------------------
    // MassDelta
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::MassDelta>(m, "MassDelta", "Mass delta modification with optional source hint")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::MassDelta &>())
        .def("__copy__", [](const OpenMS::ProForma::MassDelta& self) { return OpenMS::ProForma::MassDelta(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::MassDelta& self, nb::dict) { return OpenMS::ProForma::MassDelta(self); }, "memo"_a)
        .def_rw("mass", &OpenMS::ProForma::MassDelta::mass)
        .def_rw("original_text", &OpenMS::ProForma::MassDelta::original_text)
        ;

    // -----------------------------------------------------------------------
    // Modification
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::Modification>(m, "Modification", "A modification with one or more alternative tags")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::Modification &>())
        .def("__copy__", [](const OpenMS::ProForma::Modification& self) { return OpenMS::ProForma::Modification(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::Modification& self, nb::dict) { return OpenMS::ProForma::Modification(self); }, "memo"_a)
        ;

    // -----------------------------------------------------------------------
    // ModificationDefinition
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ModificationDefinition>(m, "ModificationDefinition", "Representation of modification definition")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ModificationDefinition &>())
        .def("__copy__", [](const OpenMS::ModificationDefinition& self) { return OpenMS::ModificationDefinition(self); })
        .def("__deepcopy__", [](const OpenMS::ModificationDefinition& self, nb::dict) { return OpenMS::ModificationDefinition(self); }, "memo"_a)
        .def(nb::init<OpenMS::String, bool, unsigned int>())
        .def(nb::init<OpenMS::ResidueModification, bool, unsigned int>())
        .def("setFixedModification", [](OpenMS::ModificationDefinition& self, bool fixed) { return self.setFixedModification(fixed); }, "fixed"_a, "Sets whether this modification definition is fixed or variable (modification must occur vs. can occur)")
        .def("isFixedModification", [](const OpenMS::ModificationDefinition& self) { return self.isFixedModification(); }, "Returns if the modification if fixed true, else false")
        .def("setMaxOccurrences", [](OpenMS::ModificationDefinition& self, unsigned int num) { return self.setMaxOccurrences(num); }, "num"_a, "Sets the maximal number of occurrences per peptide (unbounded if 0)")
        .def("getMaxOccurrences", [](const OpenMS::ModificationDefinition& self) { return self.getMaxOccurrences(); }, "Returns the maximal number of occurrences per peptide")
        .def("getModificationName", [](const OpenMS::ModificationDefinition& self) { return self.getModificationName(); }, "Returns the name of the modification")
        .def("setModification", [](OpenMS::ModificationDefinition& self, const OpenMS::String& modification) { return self.setModification(modification); }, "modification"_a, "Sets the modification, allowed are unique names provided by ModificationsDB")
        .def("getModification", [](const OpenMS::ModificationDefinition& self) -> const OpenMS::ResidueModification & { return self.getModification(); }, nb::rv_policy::reference_internal)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self < nb::self)
        .def("__hash__", [](const OpenMS::ModificationDefinition& self) { return std::hash<OpenMS::ModificationDefinition>{}(self); })
        ;

    // -----------------------------------------------------------------------
    // ModificationDefinitionsSet
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ModificationDefinitionsSet>(m, "ModificationDefinitionsSet", 
        R"doc(
Representation of a set of modification definitions
This class enhances the modification definitions as defined in the
class ModificationDefinition into a set of definitions. This is also
e.g. used as input parameters in search engines.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ModificationDefinitionsSet &>())
        .def("__copy__", [](const OpenMS::ModificationDefinitionsSet& self) { return OpenMS::ModificationDefinitionsSet(self); })
        .def("__deepcopy__", [](const OpenMS::ModificationDefinitionsSet& self, nb::dict) { return OpenMS::ModificationDefinitionsSet(self); }, "memo"_a)
        .def(nb::init<std::vector<OpenMS::String>, std::vector<OpenMS::String>>())
        .def("setMaxModifications", [](OpenMS::ModificationDefinitionsSet& self, size_t max_mod) { return self.setMaxModifications(max_mod); }, "max_mod"_a, "Sets the maximal number of modifications allowed per peptide")
        .def("getMaxModifications", [](const OpenMS::ModificationDefinitionsSet& self) { return self.getMaxModifications(); }, "Return the maximal number of modifications allowed per peptide")
        .def("getNumberOfModifications", [](const OpenMS::ModificationDefinitionsSet& self) { return self.getNumberOfModifications(); }, "Returns the number of modifications stored in this set")
        .def("getNumberOfFixedModifications", [](const OpenMS::ModificationDefinitionsSet& self) { return self.getNumberOfFixedModifications(); }, "Returns the number of fixed modifications stored in this set")
        .def("getNumberOfVariableModifications", [](const OpenMS::ModificationDefinitionsSet& self) { return self.getNumberOfVariableModifications(); }, "Returns the number of variable modifications stored in this set")
        .def("addModification", [](OpenMS::ModificationDefinitionsSet& self, const OpenMS::ModificationDefinition& mod_def) { return self.addModification(mod_def); }, "mod_def"_a, "Adds a modification definition to the set")
        .def("setModifications", [](OpenMS::ModificationDefinitionsSet& self, const std::set<OpenMS::ModificationDefinition>& mod_defs) { return self.setModifications(mod_defs); }, "mod_defs"_a, "Sets the modification definitions")
        .def("setModifications", [](OpenMS::ModificationDefinitionsSet& self, const OpenMS::String& fixed_modifications, const OpenMS::String& variable_modifications) { return self.setModifications(fixed_modifications, variable_modifications); }, "fixed_modifications"_a, "variable_modifications"_a)
        .def("setModifications", [](OpenMS::ModificationDefinitionsSet& self, const std::vector<OpenMS::String>& fixed_modifications, const std::vector<OpenMS::String>& variable_modifications) { return self.setModifications(fixed_modifications, variable_modifications); }, "fixed_modifications"_a, "variable_modifications"_a)
        .def("getModifications", [](const OpenMS::ModificationDefinitionsSet& self) { return self.getModifications(); }, "Returns the stored modification definitions")
        .def("getFixedModifications", [](const OpenMS::ModificationDefinitionsSet& self) -> const std::set<OpenMS::ModificationDefinition> & { return self.getFixedModifications(); }, nb::rv_policy::reference_internal, "Returns the stored fixed modification definitions")
        .def("getVariableModifications", [](const OpenMS::ModificationDefinitionsSet& self) -> const std::set<OpenMS::ModificationDefinition> & { return self.getVariableModifications(); }, nb::rv_policy::reference_internal, "Returns the stored variable modification definitions")
        .def("getModificationNames", [](const OpenMS::ModificationDefinitionsSet& self) { return self.getModificationNames(); }, "Returns only the names of the modifications stored in the set")
        .def("getFixedAndVariableModificationNames", [](const OpenMS::ModificationDefinitionsSet& self) { std::vector<OpenMS::String> fixed_modifications, variable_modifications; self.getModificationNames(fixed_modifications, variable_modifications); return nb::make_tuple(fixed_modifications, variable_modifications); }, "Returns a tuple of (fixed_modification_names, variable_modification_names)")
        .def("getFixedModificationNames", [](const OpenMS::ModificationDefinitionsSet& self) { return self.getFixedModificationNames(); }, "Returns only the names of the fixed modifications")
        .def("getVariableModificationNames", [](const OpenMS::ModificationDefinitionsSet& self) { return self.getVariableModificationNames(); }, "Returns only the names of the variable modifications")
        .def("isCompatible", [](const OpenMS::ModificationDefinitionsSet& self, const OpenMS::AASequence& peptide) { return self.isCompatible(peptide); }, "peptide"_a, "Returns true if the peptide is compatible with the definitions, e.g. does not contain other modifications")
        .def("inferFromPeptides", [](OpenMS::ModificationDefinitionsSet& self, const OpenMS::PeptideIdentificationList& peptides) { return self.inferFromPeptides(peptides); }, "peptides"_a, "Infers the sets of defined modifications from the modifications present on peptide identifications")
        ;

    // -----------------------------------------------------------------------
    // ModificationsDB
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ModificationsDB>(m, "ModificationsDB", 
        R"doc(
Database of modifications (amino acid modifications). This is a singleton class.
The modifications are read from the unimod.xml file on construction.
)doc")
        .def_static("isInstantiated", []() { return OpenMS::ModificationsDB::isInstantiated(); }, "Check whether ModificationsDB was instantiated before")
        .def("has", [](const OpenMS::ModificationsDB& self, const OpenMS::String& modification) { return self.has(modification); }, "modification"_a, "Returns True if the modification exists")
        .def("findModificationIndex", [](const OpenMS::ModificationsDB& self, const OpenMS::String& mod_name) { return self.findModificationIndex(mod_name); }, "mod_name"_a, "Returns the index of the modification in the mods_ vector; a unique name must be given")

        .def("searchModifications", [](const OpenMS::ModificationsDB& self, const OpenMS::String& mod_name, const OpenMS::String& residue, nb::object term_spec_obj) {
            int term_spec = nb::cast<int>(nb::int_(term_spec_obj));
            std::set<const OpenMS::ResidueModification*> mods;
            self.searchModifications(mods, mod_name, residue, static_cast<OpenMS::ResidueModification::TermSpecificity>(term_spec));
            nb::list result;
            for (auto* m : mods) {
                result.append(nb::cast(m, nb::rv_policy::reference));
            }
            return result;
        }, "mod_name"_a, "residue"_a = "", "term_spec"_a = nb::int_(static_cast<int>(OpenMS::ResidueModification::TermSpecificity::NUMBER_OF_TERM_SPECIFICITY)), "Search for modifications by name")

        .def("getNumberOfModifications", [](const OpenMS::ModificationsDB& self) {
            return self.getNumberOfModifications();
        }, "Get the number of modifications")
        .def("getModification", [](const OpenMS::ModificationsDB& self, OpenMS::Size index) {
            return self.getModification(index);
        }, "index"_a, nb::rv_policy::reference, "Returns the modification with the given index")
        .def("getModification", [](const OpenMS::ModificationsDB& self, const OpenMS::String& mod_name, const OpenMS::String& residue, nb::object term_spec_obj) {
            int term_spec = nb::cast<int>(nb::int_(term_spec_obj));
            return self.getModification(mod_name, residue, static_cast<OpenMS::ResidueModification::TermSpecificity>(term_spec));
        }, "mod_name"_a, "residue"_a = "", "term_spec"_a = nb::int_(static_cast<int>(OpenMS::ResidueModification::TermSpecificity::NUMBER_OF_TERM_SPECIFICITY)), nb::rv_policy::reference, "Returns the modification with the given name, residue, and term specificity")
        .def("getAllSearchModifications", [](const OpenMS::ModificationsDB& self) {
            std::vector<OpenMS::String> mods;
            self.getAllSearchModifications(mods);
            return mods;
        }, "Returns all modifications that can be used for identification searches")
        .def("getBestModificationByDiffMonoMass", [](OpenMS::ModificationsDB& self, double mass, double max_error, const OpenMS::String& residue, nb::object term_spec_obj) {
            int term_spec = nb::cast<int>(nb::int_(term_spec_obj));
            return self.getBestModificationByDiffMonoMass(mass, max_error, residue, static_cast<OpenMS::ResidueModification::TermSpecificity>(term_spec));
        }, "mass"_a, "max_error"_a, "residue"_a = "", "term_spec"_a = nb::int_(static_cast<int>(OpenMS::ResidueModification::TermSpecificity::NUMBER_OF_TERM_SPECIFICITY)), nb::rv_policy::reference, "Returns the best modification by diff mono mass")
        .def("searchModificationsByDiffMonoMass", [](OpenMS::ModificationsDB& self, double mass, double max_error, const OpenMS::String& residue, nb::object term_spec_obj) {
            int term_spec = nb::cast<int>(nb::int_(term_spec_obj));
            std::vector<OpenMS::String> mods;
            self.searchModificationsByDiffMonoMass(mods, mass, max_error, residue, static_cast<OpenMS::ResidueModification::TermSpecificity>(term_spec));
            return mods;
        }, "mass"_a, "max_error"_a, "residue"_a = "", "term_spec"_a = nb::int_(static_cast<int>(OpenMS::ResidueModification::TermSpecificity::NUMBER_OF_TERM_SPECIFICITY)), "Search modifications by difference in mono mass (returns list of modification names)")
        .def("addModification", [](OpenMS::ModificationsDB& self, const OpenMS::ResidueModification& new_mod) {
            return self.addModification(new_mod);
        }, "new_mod"_a, nb::rv_policy::reference, "Add a new modification to the database (makes a copy)")
        .def_static("getInstance", []() -> OpenMS::ModificationsDB* { return OpenMS::ModificationsDB::getInstance(); }, nb::rv_policy::reference, "Returns the singleton instance")
        ;

    // -----------------------------------------------------------------------
    // CrossLinksDB
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CrossLinksDB, OpenMS::ModificationsDB>(m, "CrossLinksDB",
        R"doc(
Database of cross-linker modifications. This is a singleton class that inherits from ModificationsDB.
The cross-linker modifications are read from an OBO file.
)doc")
        .def_static("getInstance", []() -> OpenMS::CrossLinksDB* { return OpenMS::CrossLinksDB::getInstance(); }, nb::rv_policy::reference, "Returns the singleton instance")
        .def_static("isInstantiated", []() {
            // Access getInstance and check — CrossLinksDB doesn't override isInstantiated,
            // but we need to ensure the DB is loaded
            return OpenMS::CrossLinksDB::getInstance() != nullptr;
        }, "Check whether CrossLinksDB was instantiated")
        .def("getAllSearchModifications", [](const OpenMS::CrossLinksDB& self) {
            std::vector<OpenMS::String> mods;
            self.getAllSearchModifications(mods);
            return mods;
        }, "Returns all cross-link modifications that can be used for identification searches")
        ;

    // -----------------------------------------------------------------------
    // ModifiedPeptideGenerator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ModifiedPeptideGenerator>(m, "ModifiedPeptideGenerator", "Generates modified peptides/proteins.")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ModifiedPeptideGenerator &>())
        .def("__copy__", [](const OpenMS::ModifiedPeptideGenerator& self) { return OpenMS::ModifiedPeptideGenerator(self); })
        .def("__deepcopy__", [](const OpenMS::ModifiedPeptideGenerator& self, nb::dict) { return OpenMS::ModifiedPeptideGenerator(self); }, "memo"_a)
        .def_static("getModifications", [](const std::vector<OpenMS::String>& modNames) { return OpenMS::ModifiedPeptideGenerator::getModifications(modNames); }, "modNames"_a)
        .def_static("applyFixedModifications", [](const OpenMS::ModifiedPeptideGenerator::MapToResidueType& fixed_mods, OpenMS::AASequence& peptide) { return OpenMS::ModifiedPeptideGenerator::applyFixedModifications(fixed_mods, peptide); }, "fixed_mods"_a, "peptide"_a)
        .def_static("applyVariableModifications", [](const OpenMS::ModifiedPeptideGenerator::MapToResidueType& var_mods, const OpenMS::AASequence& peptide, size_t max_variable_mods_per_peptide, bool keep_original) { std::vector<OpenMS::AASequence> all_modified_peptides; OpenMS::ModifiedPeptideGenerator::applyVariableModifications(var_mods, peptide, max_variable_mods_per_peptide, all_modified_peptides, keep_original); return all_modified_peptides; }, "var_mods"_a, "peptide"_a, "max_variable_mods_per_peptide"_a, "keep_original"_a)
        ;

    // -----------------------------------------------------------------------
    // ModifiedPeptideGenerator_MapToResidueType
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ModifiedPeptideGenerator::MapToResidueType>(m, "ModifiedPeptideGenerator_MapToResidueType", "OpenMS class ModifiedPeptideGenerator_MapToResidueType")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ModifiedPeptideGenerator::MapToResidueType &>())
        .def("__copy__", [](const OpenMS::ModifiedPeptideGenerator::MapToResidueType& self) { return OpenMS::ModifiedPeptideGenerator::MapToResidueType(self); })
        .def("__deepcopy__", [](const OpenMS::ModifiedPeptideGenerator::MapToResidueType& self, nb::dict) { return OpenMS::ModifiedPeptideGenerator::MapToResidueType(self); }, "memo"_a)
        .def_ro("val", &OpenMS::ModifiedPeptideGenerator::MapToResidueType::val)
        ;

    // -----------------------------------------------------------------------
    // MzPAFAnnotation
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzPAFAnnotation>(m, "MzPAFAnnotation", 
        R"doc(
A single mzPAF peak annotation.
Represents one annotation for a peak in mzPAF (Peak Annotation Format),
the HUPO-PSI standard for fragment ion annotations.
Examples:
- y4 - Simple y-ion at position 4
- b2-H2O - b-ion with neutral loss
- y4^2 - Doubly charged y-ion
- y4/0.001*0.75 - With mass delta and confidence
- IY - Immonium ion (tyrosine)
- m3:6 - Internal fragment
- r[TMT127N] - Reporter ion
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MzPAFAnnotation &>())
        .def("__copy__", [](const OpenMS::MzPAFAnnotation& self) { return OpenMS::MzPAFAnnotation(self); })
        .def("__deepcopy__", [](const OpenMS::MzPAFAnnotation& self, nb::dict) { return OpenMS::MzPAFAnnotation(self); }, "memo"_a)
        .def("isValid", [](const OpenMS::MzPAFAnnotation& self) { return self.isValid(); })
        .def(nb::self == nb::self)
        .def_rw("analyte_index", &OpenMS::MzPAFAnnotation::analyte_index)
        .def_rw("ion_series", &OpenMS::MzPAFAnnotation::ion_series)
        .def_rw("ordinal", &OpenMS::MzPAFAnnotation::ordinal)
        .def_rw("immonium_residue", &OpenMS::MzPAFAnnotation::immonium_residue)
        .def_rw("internal_range", &OpenMS::MzPAFAnnotation::internal_range)
        .def_rw("reporter_name", &OpenMS::MzPAFAnnotation::reporter_name)
        .def_rw("formula", &OpenMS::MzPAFAnnotation::formula)
        .def_rw("named_compound", &OpenMS::MzPAFAnnotation::named_compound)
        .def_rw("neutral_losses", &OpenMS::MzPAFAnnotation::neutral_losses)
        .def_rw("isotope_offset", &OpenMS::MzPAFAnnotation::isotope_offset)
        .def_rw("adduct", &OpenMS::MzPAFAnnotation::adduct)
        .def_rw("charge", &OpenMS::MzPAFAnnotation::charge)
        .def_rw("mass_delta", &OpenMS::MzPAFAnnotation::mass_delta)
        .def_rw("confidence", &OpenMS::MzPAFAnnotation::confidence)
        .def_rw("embedded_sequence", &OpenMS::MzPAFAnnotation::embedded_sequence)
        ;

    // -----------------------------------------------------------------------
    // MzPAFMassDelta
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzPAFMassDelta>(m, "MzPAFMassDelta", "Mass delta in an mzPAF annotation (e.g., /0.001, /-1.4ppm)")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MzPAFMassDelta &>())
        .def("__copy__", [](const OpenMS::MzPAFMassDelta& self) { return OpenMS::MzPAFMassDelta(self); })
        .def("__deepcopy__", [](const OpenMS::MzPAFMassDelta& self, nb::dict) { return OpenMS::MzPAFMassDelta(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def_rw("value", &OpenMS::MzPAFMassDelta::value)
        .def_rw("unit", &OpenMS::MzPAFMassDelta::unit)
        ;

    // -----------------------------------------------------------------------
    // MzPAFNeutralLoss
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzPAFNeutralLoss>(m, "MzPAFNeutralLoss", "Neutral loss in an mzPAF annotation (e.g., -H2O, -NH3)")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MzPAFNeutralLoss &>())
        .def("__copy__", [](const OpenMS::MzPAFNeutralLoss& self) { return OpenMS::MzPAFNeutralLoss(self); })
        .def("__deepcopy__", [](const OpenMS::MzPAFNeutralLoss& self, nb::dict) { return OpenMS::MzPAFNeutralLoss(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def_rw("formula", &OpenMS::MzPAFNeutralLoss::formula)
        ;

    // -----------------------------------------------------------------------
    // MzPAFPeakAnnotations
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzPAFPeakAnnotations>(m, "MzPAFPeakAnnotations", "Multiple mzPAF annotations for a single peak (comma-separated alternatives)")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MzPAFPeakAnnotations &>())
        .def("__copy__", [](const OpenMS::MzPAFPeakAnnotations& self) { return OpenMS::MzPAFPeakAnnotations(self); })
        .def("__deepcopy__", [](const OpenMS::MzPAFPeakAnnotations& self, nb::dict) { return OpenMS::MzPAFPeakAnnotations(self); }, "memo"_a)
        .def("empty", [](const OpenMS::MzPAFPeakAnnotations& self) { return self.empty(); })
        .def("size", [](const OpenMS::MzPAFPeakAnnotations& self) { return self.size(); })
        .def(nb::self == nb::self)
        .def_rw("annotations", &OpenMS::MzPAFPeakAnnotations::annotations)
        .def("__len__", [](OpenMS::MzPAFPeakAnnotations& self) { return self.size(); })
        ;

    // -----------------------------------------------------------------------
    // MzPAFIonSeries
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::MzPAFIonSeries>(m, "MzPAFIonSeries", "Ion series types for mzPAF peak annotations", nb::is_arithmetic())
        .value("A", OpenMS::MzPAFIonSeries::A)
        .value("B", OpenMS::MzPAFIonSeries::B)
        .value("C", OpenMS::MzPAFIonSeries::C)
        .value("X", OpenMS::MzPAFIonSeries::X)
        .value("Y", OpenMS::MzPAFIonSeries::Y)
        .value("Z", OpenMS::MzPAFIonSeries::Z)
        .value("PRECURSOR", OpenMS::MzPAFIonSeries::PRECURSOR)
        .value("IMMONIUM", OpenMS::MzPAFIonSeries::IMMONIUM)
        .value("INTERNAL", OpenMS::MzPAFIonSeries::INTERNAL)
        .value("REPORTER", OpenMS::MzPAFIonSeries::REPORTER)
        .value("FORMULA", OpenMS::MzPAFIonSeries::FORMULA)
        .value("NAMED", OpenMS::MzPAFIonSeries::NAMED)
        .value("UNKNOWN", OpenMS::MzPAFIonSeries::UNKNOWN)

        ;

    // -----------------------------------------------------------------------
    // MzPAFDeltaUnit
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::MzPAFDeltaUnit>(m, "MzPAFDeltaUnit", "Unit for mass delta values in mzPAF annotations", nb::is_arithmetic())
        .value("DALTON", OpenMS::MzPAFDeltaUnit::DALTON)
        .value("PPM", OpenMS::MzPAFDeltaUnit::PPM)

        ;

    // -----------------------------------------------------------------------
    // MzPAFErrorCode
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::MzPAFErrorCode>(m, "MzPAFErrorCode", "Error codes for mzPAF parsing errors", nb::is_arithmetic())
        .value("UNEXPECTED_CHARACTER", OpenMS::MzPAFErrorCode::UNEXPECTED_CHARACTER)
        .value("UNCLOSED_BRACKET", OpenMS::MzPAFErrorCode::UNCLOSED_BRACKET)
        .value("INVALID_ION_SERIES", OpenMS::MzPAFErrorCode::INVALID_ION_SERIES)
        .value("INVALID_NUMBER", OpenMS::MzPAFErrorCode::INVALID_NUMBER)
        .value("INVALID_FORMULA", OpenMS::MzPAFErrorCode::INVALID_FORMULA)
        .value("INVALID_CHARGE", OpenMS::MzPAFErrorCode::INVALID_CHARGE)
        .value("INVALID_DELTA", OpenMS::MzPAFErrorCode::INVALID_DELTA)
        .value("INVALID_CONFIDENCE", OpenMS::MzPAFErrorCode::INVALID_CONFIDENCE)
        .value("EMPTY_INPUT", OpenMS::MzPAFErrorCode::EMPTY_INPUT)
        .value("UNEXPECTED_END_OF_INPUT", OpenMS::MzPAFErrorCode::UNEXPECTED_END_OF_INPUT)
        .value("INTERNAL_ERROR", OpenMS::MzPAFErrorCode::INTERNAL_ERROR)
        ;

    // -----------------------------------------------------------------------
    // MzPAF (static utility class)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzPAF>(m, "MzPAF", "Parser and writer for mzPAF (Peak Annotation Format) notation")
        .def_static("parse", [](const OpenMS::String& input) { return OpenMS::MzPAF::parse(input); }, "input"_a, "Parse an mzPAF string into a single annotation")
        .def_static("parseMultiple", [](const OpenMS::String& input) { return OpenMS::MzPAF::parseMultiple(input); }, "input"_a, "Parse an mzPAF string with potentially multiple annotations")
        .def_static("tryParse", [](const OpenMS::String& input) { return OpenMS::MzPAF::tryParse(input); }, "input"_a, "Try to parse an mzPAF string (returns None on failure)")
        .def_static("tryParseMultiple", [](const OpenMS::String& input) { return OpenMS::MzPAF::tryParseMultiple(input); }, "input"_a, "Try to parse multiple annotations (returns None on failure)")
        .def_static("toString", [](const OpenMS::MzPAFAnnotation& ann) { return OpenMS::MzPAF::toString(ann); }, "ann"_a, "Convert an annotation to mzPAF string")
        .def_static("toStringMultiple", [](const OpenMS::MzPAFPeakAnnotations& anns) { return OpenMS::MzPAF::toString(anns); }, "anns"_a, "Convert multiple annotations to mzPAF string")
        .def_static("toPeakAnnotation", [](const OpenMS::MzPAFAnnotation& mzpaf, double mz, double intensity) { return OpenMS::MzPAF::toPeakAnnotation(mzpaf, mz, intensity); }, "mzpaf"_a, "mz"_a, "intensity"_a, "Create a PeakAnnotation from mzPAF data")
        .def_static("fromPeakAnnotation", [](const OpenMS::PeptideHit::PeakAnnotation& peak_annotation) { return OpenMS::MzPAF::fromPeakAnnotation(peak_annotation); }, "peak_annotation"_a, "Parse mzPAF annotations from a PeakAnnotation")
        .def_static("isMzPAFFormat", [](const OpenMS::String& annotation) { return OpenMS::MzPAF::isMzPAFFormat(annotation); }, "annotation"_a, "Check if a string appears to be in mzPAF format")
        .def_static("isStandardFragmentIon", [](OpenMS::MzPAFIonSeries series) { return OpenMS::MzPAF::isStandardFragmentIon(series); }, "series"_a, "Check if ion series is a standard fragment ion (a, b, c, x, y, z)")
        .def_static("ionSeriesToChar", [](OpenMS::MzPAFIonSeries series) { return OpenMS::MzPAF::ionSeriesToChar(series); }, "series"_a, "Get the ion series character for an annotation")
        .def_static("charToIonSeries", [](char c) -> std::optional<OpenMS::MzPAFIonSeries> {
            OpenMS::MzPAFIonSeries series;
            if (OpenMS::MzPAF::charToIonSeries(c, series)) return series;
            return std::nullopt;
        }, "c"_a, "Parse ion series from character (returns None if invalid)")
        ;

    // -----------------------------------------------------------------------
    // NamedMod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::NamedMod>(m, "NamedMod", "Named modification with optional CV prefix hint")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::NamedMod &>())
        .def("__copy__", [](const OpenMS::ProForma::NamedMod& self) { return OpenMS::ProForma::NamedMod(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::NamedMod& self, nb::dict) { return OpenMS::ProForma::NamedMod(self); }, "memo"_a)
        .def_rw("name", &OpenMS::ProForma::NamedMod::name)
        ;

    // -----------------------------------------------------------------------
    // Peptidoform
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::Peptidoform>(m, "Peptidoform", 
        R"doc(
A single peptidoform (one peptide chain) with modifications
Represents a complete peptide chain including global modifications,
unlocalised modifications, labile modifications, terminal modifications,
and the amino acid sequence with modifications.
Note: The following fields are intentionally not exposed and should be
accessed via the ProForma parser/writer API functions:
- name: optional peptidoform name (v2.1 extension)
- sequence: the amino acid sequence with modifications (complex type)
- global_mods: global modifications like isotope labels (complex type)
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::Peptidoform &>())
        .def("__copy__", [](const OpenMS::ProForma::Peptidoform& self) { return OpenMS::ProForma::Peptidoform(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::Peptidoform& self, nb::dict) { return OpenMS::ProForma::Peptidoform(self); }, "memo"_a)
        .def_rw("unlocalised_mods", &OpenMS::ProForma::Peptidoform::unlocalised_mods)
        .def_rw("labile_mods", &OpenMS::ProForma::Peptidoform::labile_mods)
        .def_rw("n_term_mods", &OpenMS::ProForma::Peptidoform::n_term_mods)
        .def_rw("c_term_mods", &OpenMS::ProForma::Peptidoform::c_term_mods)
        .def("__repr__", [](const OpenMS::ProForma::Peptidoform& self) {
            std::ostringstream oss;
            oss << "Peptidoform(length=" << self.sequence.size() << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::ProForma::Peptidoform& self) { return nb::cast(self).attr("__repr__")(); })
        ;

    // -----------------------------------------------------------------------
    // PeptidoformIon
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::PeptidoformIon>(m, "PeptidoformIon", 
        R"doc(
A peptidoform ion (one or more chains with optional charge)
Note: The following fields are intentionally not exposed and should be
accessed via the ProForma parser/writer API functions:
- name: optional ion name (v2.1 extension)
- charge: charge state specification (complex variant type)
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::PeptidoformIon &>())
        .def("__copy__", [](const OpenMS::ProForma::PeptidoformIon& self) { return OpenMS::ProForma::PeptidoformIon(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::PeptidoformIon& self, nb::dict) { return OpenMS::ProForma::PeptidoformIon(self); }, "memo"_a)
        .def_rw("chains", &OpenMS::ProForma::PeptidoformIon::chains)
        .def("__repr__", [](const OpenMS::ProForma::PeptidoformIon& self) {
            std::ostringstream oss;
            oss << "PeptidoformIon(num_chains=" << self.chains.size() << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::ProForma::PeptidoformIon& self) { return nb::cast(self).attr("__repr__")(); })
        ;

    // -----------------------------------------------------------------------
    // ProForma
    // -----------------------------------------------------------------------
    auto proforma_class = nb::class_<OpenMS::ProForma>(m, "ProForma", 
        R"doc(
ProForma v2 peptidoform notation parser and data structures.
This class provides parsing, serialization, and conversion functionality for
the ProForma v2 peptidoform notation standard. It contains nested types that
form the Abstract Syntax Tree (AST) representation of parsed ProForma strings.
All methods are static. Use ProForma.parse() to parse a ProForma string.
Usage example:
.. code-block:: python
pf = ProForma.parse("EM[UNIMOD:35]K")
# pf now contains the parsed Peptidoform AST
s = ProForma.toString(pf, ProForma.WriteMode.LOSSLESS)
# s is "EM[UNIMOD:35]K"
)doc")
        .def_static("parse", [](const OpenMS::String& input) { return OpenMS::ProForma::parse(input); }, "input"_a, "Parse a ProForma string into a Peptidoform AST")
        .def_static("parseIon", [](const OpenMS::String& input) { return OpenMS::ProForma::parseIon(input); }, "input"_a, "Parse a ProForma string into a PeptidoformIon AST (with charge state)")
        .def_static("toString", [](const OpenMS::ProForma::Peptidoform& pf, OpenMS::ProForma::WriteMode mode) { return OpenMS::ProForma::toString(pf, mode); }, "pf"_a, "mode"_a, "Convert a Peptidoform AST back to ProForma string notation with specified mode")
        .def_static("toString", [](const OpenMS::ProForma::PeptidoformIon& pfi, OpenMS::ProForma::WriteMode mode) { return OpenMS::ProForma::toString(pfi, mode); }, "pfi"_a, "mode"_a, "Convert a Peptidoform AST back to ProForma string notation with specified mode")
        .def_static("peptidoformFromJSON", [](const OpenMS::String& json_str) { return OpenMS::ProForma::peptidoformFromJSON(json_str); }, "json_str"_a, "Construct Peptidoform from JSON string")
        .def_static("peptidoformIonFromJSON", [](const OpenMS::String& json_str) { return OpenMS::ProForma::peptidoformIonFromJSON(json_str); }, "json_str"_a, "Construct PeptidoformIon from JSON string")
        .def_static("resolveModifications", [](OpenMS::ProForma::Peptidoform& pf) { return OpenMS::ProForma::resolveModifications(pf); }, "pf"_a, "Resolve all modifications in a Peptidoform using ModificationsDB")
        .def_static("toAASequence", [](const OpenMS::ProForma::Peptidoform& pf, OpenMS::ProForma::ConversionPolicy policy) { return OpenMS::ProForma::toAASequence(pf, policy); }, "pf"_a, "policy"_a, "Convert a Peptidoform to an OpenMS AASequence")
        .def_static("fromAASequence", [](const OpenMS::AASequence& seq) { return OpenMS::ProForma::fromAASequence(seq); }, "seq"_a, "Create a Peptidoform from an OpenMS AASequence")
        .def_static("isRepresentableAsAASequence", [](const OpenMS::ProForma::Peptidoform& pf) { return OpenMS::ProForma::isRepresentableAsAASequence(pf); }, "pf"_a, "Check if a Peptidoform can be fully represented as an AASequence")
        .def_static("getAASequenceConversionIssues", [](const OpenMS::ProForma::Peptidoform& pf) { return OpenMS::ProForma::getAASequenceConversionIssues(pf); }, "pf"_a, "Get a list of all issues that would arise during AASequence conversion")
        .def_static("canCalculateMass", [](const OpenMS::ProForma::Peptidoform& pf) { return OpenMS::ProForma::canCalculateMass(pf); }, "pf"_a, "Check if mass can be calculated for a Peptidoform")
        .def_static("canCalculateMass", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::canCalculateMass(pfi); }, "pfi"_a, "Check if mass can be calculated for a Peptidoform")
        .def_static("getMassCalculationIssues", [](const OpenMS::ProForma::Peptidoform& pf) { return OpenMS::ProForma::getMassCalculationIssues(pf); }, "pf"_a, "Get issues preventing mass calculation for a Peptidoform")
        .def_static("getMassCalculationIssues", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::getMassCalculationIssues(pfi); }, "pfi"_a, "Get issues preventing mass calculation for a Peptidoform")
        .def_static("getMonoWeight", [](const OpenMS::ProForma::Peptidoform& pf) { return OpenMS::ProForma::getMonoWeight(pf); }, "pf"_a, "Calculate monoisotopic mass of a Peptidoform in Daltons")
        .def_static("getMonoWeight", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::getMonoWeight(pfi); }, "pfi"_a, "Calculate monoisotopic mass of a Peptidoform in Daltons")
        .def_static("getMZ", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::getMZ(pfi); }, "pfi"_a, "Calculate m/z for a PeptidoformIon using its charge state")
        .def_static("getMZ", [](const OpenMS::ProForma::Peptidoform& pf, int charge) { return OpenMS::ProForma::getMZ(pf, charge); }, "pf"_a, "charge"_a, "Calculate m/z for a PeptidoformIon using its charge state")
        .def_static("canGenerateSpectrum", [](const OpenMS::ProForma::Peptidoform& pf) { return OpenMS::ProForma::canGenerateSpectrum(pf); }, "pf"_a, "Check if a theoretical spectrum can be generated for a Peptidoform")
        .def_static("canGenerateSpectrum", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::canGenerateSpectrum(pfi); }, "pfi"_a, "Check if a theoretical spectrum can be generated for a Peptidoform")
        .def_static("getSpectrumGenerationIssues", [](const OpenMS::ProForma::Peptidoform& pf) { return OpenMS::ProForma::getSpectrumGenerationIssues(pf); }, "pf"_a, "Get issues preventing spectrum generation for a Peptidoform")
        .def_static("getSpectrumGenerationIssues", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::getSpectrumGenerationIssues(pfi); }, "pfi"_a, "Get issues preventing spectrum generation for a Peptidoform")
        .def_static("generateSpectrum", [](const OpenMS::ProForma::Peptidoform& pf, int min_charge, int max_charge, const OpenMS::String& ion_types, bool add_losses, bool add_metainfo) { return OpenMS::ProForma::generateSpectrum(pf, min_charge, max_charge, ion_types, add_losses, add_metainfo); }, "pf"_a, "min_charge"_a, "max_charge"_a, "ion_types"_a, "add_losses"_a, "add_metainfo"_a, 
            R"doc(
Generate theoretical MS/MS spectrum for a Peptidoform. ion_types uses chars a,b,c,x,y,z for ion series, M for precursor, I for immonium (e.g. "by" or "abyM")
)doc")
        .def_static("generateSpectrum", [](const OpenMS::ProForma::PeptidoformIon& pfi, int min_charge, int max_charge, const OpenMS::String& ion_types, bool add_losses, bool add_metainfo) { return OpenMS::ProForma::generateSpectrum(pfi, min_charge, max_charge, ion_types, add_losses, add_metainfo); }, "pfi"_a, "min_charge"_a, "max_charge"_a, "ion_types"_a, "add_losses"_a, "add_metainfo"_a, 
            R"doc(
Generate theoretical MS/MS spectrum for a Peptidoform. ion_types uses chars a,b,c,x,y,z for ion series, M for precursor, I for immonium (e.g. "by" or "abyM")
)doc")
        .def_static("toStringIon", [](const OpenMS::ProForma::PeptidoformIon& pfi, OpenMS::ProForma::WriteMode mode) { return OpenMS::ProForma::toString(pfi, mode); }, "pfi"_a, "mode"_a, "Convert a PeptidoformIon AST back to ProForma string notation with specified mode")
        .def_static("canCalculateMassIon", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::canCalculateMass(pfi); }, "pfi"_a, "Check if mass can be calculated for a PeptidoformIon")
        .def_static("getMassCalculationIssuesIon", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::getMassCalculationIssues(pfi); }, "pfi"_a, "Get issues preventing mass calculation for a PeptidoformIon")
        .def_static("getMonoWeightIon", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::getMonoWeight(pfi); }, "pfi"_a, "Calculate monoisotopic mass of a PeptidoformIon in Daltons")
        .def_static("getMZCharge", [](const OpenMS::ProForma::Peptidoform& pf, int charge) { return OpenMS::ProForma::getMZ(pf, charge); }, "pf"_a, "charge"_a, "Calculate m/z for a Peptidoform at a given charge state")
        .def_static("canGenerateSpectrumIon", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::canGenerateSpectrum(pfi); }, "pfi"_a, "Check if a theoretical spectrum can be generated for a PeptidoformIon")
        .def_static("getSpectrumGenerationIssuesIon", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::getSpectrumGenerationIssues(pfi); }, "pfi"_a, "Get issues preventing spectrum generation for a PeptidoformIon")
        .def_static("generateSpectrumIon", [](const OpenMS::ProForma::PeptidoformIon& pfi, int min_charge, int max_charge, const std::string& ion_types, bool add_losses, bool add_metainfo) { return OpenMS::ProForma::generateSpectrum(pfi, min_charge, max_charge, ion_types, add_losses, add_metainfo); }, "pfi"_a, "min_charge"_a, "max_charge"_a, "ion_types"_a, "add_losses"_a, "add_metainfo"_a, "Generate theoretical MS/MS spectrum for a PeptidoformIon (supports cross-linked peptides). ion_types uses chars a,b,c,x,y,z for ion series, M for precursor, I for immonium")
        .def_static("peptidoformToJSON", [](const OpenMS::ProForma::Peptidoform& pf) { return OpenMS::ProForma::toJSON(pf); }, "pf"_a, "Convert Peptidoform to JSON string representation")
        .def_static("peptidoformIonToJSON", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::toJSON(pfi); }, "pfi"_a, "Convert PeptidoformIon to JSON string representation")
        .def_static("toJSON", [](const OpenMS::ProForma::Peptidoform& pf) { return OpenMS::ProForma::toJSON(pf); }, "pf"_a, "Convert Peptidoform to JSON string")
        .def_static("toJSONIon", [](const OpenMS::ProForma::PeptidoformIon& pfi) { return OpenMS::ProForma::toJSON(pfi); }, "pfi"_a, "Convert PeptidoformIon to JSON string")
        ;
    // WriteMode enum nested under ProForma
    nb::enum_<OpenMS::ProForma::WriteMode>(proforma_class, "WriteMode", nb::is_arithmetic())
        .value("LOSSLESS", OpenMS::ProForma::WriteMode::LOSSLESS)
        .value("CANONICAL", OpenMS::ProForma::WriteMode::CANONICAL)
        ;
    // ConversionPolicy enum nested under ProForma
    nb::enum_<OpenMS::ProForma::ConversionPolicy>(proforma_class, "ConversionPolicy", nb::is_arithmetic())
        .value("FAIL_ON_LOSS", OpenMS::ProForma::ConversionPolicy::FAIL_ON_LOSS)
        .value("DROP_UNLOCALISED", OpenMS::ProForma::ConversionPolicy::DROP_UNLOCALISED)
        .value("BEST_EFFORT", OpenMS::ProForma::ConversionPolicy::BEST_EFFORT)
        ;
    // ConversionIssueType enum nested under ProForma
    nb::enum_<OpenMS::ProForma::ConversionIssueType>(proforma_class, "ConversionIssueType", nb::is_arithmetic())
        .value("UNRESOLVED_MOD", OpenMS::ProForma::ConversionIssueType::UNRESOLVED_MOD)
        .value("UNLOCALISED_MOD", OpenMS::ProForma::ConversionIssueType::UNLOCALISED_MOD)
        .value("LABILE_MOD", OpenMS::ProForma::ConversionIssueType::LABILE_MOD)
        .value("GLOBAL_MOD", OpenMS::ProForma::ConversionIssueType::GLOBAL_MOD)
        .value("AMBIGUOUS_MOD", OpenMS::ProForma::ConversionIssueType::AMBIGUOUS_MOD)
        .value("AMBIGUOUS_REGION", OpenMS::ProForma::ConversionIssueType::AMBIGUOUS_REGION)
        .value("MODIFIED_RANGE", OpenMS::ProForma::ConversionIssueType::MODIFIED_RANGE)
        .value("CROSS_LINK", OpenMS::ProForma::ConversionIssueType::CROSS_LINK)
        .value("MULTIPLE_CHAINS", OpenMS::ProForma::ConversionIssueType::MULTIPLE_CHAINS)
        .value("ALTERNATIVE_MODS", OpenMS::ProForma::ConversionIssueType::ALTERNATIVE_MODS)
        .value("UNSUPPORTED_FEATURE", OpenMS::ProForma::ConversionIssueType::UNSUPPORTED_FEATURE)
        ;
    // CvDatabase enum nested under ProForma
    nb::enum_<OpenMS::ProForma::CvDatabase>(proforma_class, "CvDatabase", nb::is_arithmetic())
        .value("UNIMOD", OpenMS::ProForma::CvDatabase::UNIMOD)
        .value("MOD", OpenMS::ProForma::CvDatabase::MOD)
        .value("RESID", OpenMS::ProForma::CvDatabase::RESID)
        .value("XLMOD", OpenMS::ProForma::CvDatabase::XLMOD)
        .value("GNO", OpenMS::ProForma::CvDatabase::GNO)
        ;

    // -----------------------------------------------------------------------
    // ProteaseDB
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProteaseDB>(m, "ProteaseDB", 
        R"doc(
Database of enzymes that digest proteins (proteases). This is a singleton class.
The enzymes are read from share/CHEMISTRY/Enzymes.xml.
)doc")
        .def("getAllXTandemNames", [](const OpenMS::ProteaseDB& self) { std::vector<OpenMS::String> all_names; self.getAllXTandemNames(all_names); return all_names; }, "Returns all the enzyme names available for XTandem")
        .def("getAllCometNames", [](const OpenMS::ProteaseDB& self) { std::vector<OpenMS::String> all_names; self.getAllCometNames(all_names); return all_names; }, "Returns all the enzyme names available for Comet")
        .def("getAllOMSSANames", [](const OpenMS::ProteaseDB& self) { std::vector<OpenMS::String> all_names; self.getAllOMSSANames(all_names); return all_names; }, "Returns all the enzyme names available for OMSSA")
        .def("getAllMSGFNames", [](const OpenMS::ProteaseDB& self) { std::vector<OpenMS::String> all_names; self.getAllMSGFNames(all_names); return all_names; }, "Returns all the enzyme names available for MSGFPlus")

        .def_static("getInstance", []() -> OpenMS::ProteaseDB* { return OpenMS::ProteaseDB::getInstance(); }, nb::rv_policy::reference, "Returns the singleton instance")

        .def("getAllNames", [](const OpenMS::ProteaseDB& self, nb::list output) {
            std::vector<OpenMS::String> all_names;
            self.getAllNames(all_names);
            for (const auto& name : all_names) {
                output.append(nb::str(name.c_str()));
            }
        }, "output"_a, "Returns all enzyme names (appends to output list)")
        .def("getAllNames", [](const OpenMS::ProteaseDB& self) {
            std::vector<OpenMS::String> all_names;
            self.getAllNames(all_names);
            return all_names;
        }, "Returns all enzyme names")

        .def("getEnzyme", [](const OpenMS::ProteaseDB& self, const OpenMS::String& name) -> const OpenMS::DigestionEnzymeProtein* {
            return self.getEnzyme(name);
        }, "name"_a, nb::rv_policy::reference, "Returns the enzyme with the given name (supports synonym names)")

        .def("hasEnzyme", [](const OpenMS::ProteaseDB& self, const OpenMS::String& name) {
            return self.hasEnzyme(name);
        }, "name"_a, "Check if an enzyme with the given name exists")
        .def("getEnzymeByRegEx", [](const OpenMS::ProteaseDB& self, const OpenMS::String& cleavage_regex) -> const OpenMS::DigestionEnzymeProtein* {
            return self.getEnzymeByRegEx(cleavage_regex);
        }, "cleavage_regex"_a, nb::rv_policy::reference, "Returns the enzyme with the given cleavage regex")
        .def("hasRegEx", [](const OpenMS::ProteaseDB& self, const OpenMS::String& cleavage_regex) {
            return self.hasRegEx(cleavage_regex);
        }, "cleavage_regex"_a, "Check if an enzyme with the given regex exists")
        ;

    // -----------------------------------------------------------------------
    // ProteaseDigestion
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProteaseDigestion, OpenMS::EnzymaticDigestion>(m, "ProteaseDigestion", 
        R"doc(
EnzymaticDigestion

Class for the enzymatic digestion of proteins
Digestion can be performed using simple regular expressions, e.g. [KR] | [^P] for trypsin.
Also missed cleavages can be modeled, i.e. adjacent peptides are not cleaved
due to enzyme malfunction/access restrictions. If n missed cleavages are allowed, all possible resulting
peptides (cleaved and uncleaved) with up to n missed cleavages are returned.
Thus no random selection of just n specific missed cleavage sites is performed.
If specificity is set to semi-specific, digestion also returns semi-specific products,
i.e. with only one end at actual cleavage sites.
Usage:
.. code-block:: python
from pyopenms import *
from urllib.request import urlretrieve
#
urlretrieve ("http://www.uniprot.org/uniprot/P02769.fasta", "bsa.fasta")
#
dig = ProteaseDigestion()
dig.setEnzyme('Lys-C')
bsa_string = "".join([l.strip() for l in open("bsa.fasta").readlines()[1:]])
bsa_oms_string = String(bsa_string) # convert python string to OpenMS::String for further processing
#
minlen = 6
maxlen = 30
#
# Using AASequence and digest
result_digest = []
result_digest_min_max = []
bsa_aaseq = AASequence.fromString(bsa_oms_string)
dig.digest(bsa_aaseq, result_digest)
dig.digest(bsa_aaseq, result_digest_min_max, minlen, maxlen)
print(result_digest[4].toString()) # GLVLIAFSQYLQQCPFDEHVK
print(len(result_digest)) # 57 peptides
print(result_digest_min_max[4].toString()) # LVNELTEFAK
print(len(result_digest_min_max)) # 42 peptides
#
# Semi-specific digestion
result_semispecific = []
dig.setSpecificity(EnzymaticDigestion.SPEC_SEMI)
dig.digest(bsa_aaseq, result_semispecific)
#
# Using digestUnmodified without the need for AASequence from the EnzymaticDigestion base class
result_digest_unmodified = []
dig.digestUnmodified(StringView(bsa_oms_string), result_digest_unmodified, minlen, maxlen)
print(result_digest_unmodified[4].getString()) # LVNELTEFAK
print(len(result_digest_unmodified)) # 42 peptides
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProteaseDigestion &>())
        .def("__copy__", [](const OpenMS::ProteaseDigestion& self) { return OpenMS::ProteaseDigestion(self); })
        .def("__deepcopy__", [](const OpenMS::ProteaseDigestion& self, nb::dict) { return OpenMS::ProteaseDigestion(self); }, "memo"_a)
        .def("setEnzyme", [](OpenMS::ProteaseDigestion& self, const OpenMS::String& name) { return self.setEnzyme(name); }, "name"_a, "Sets the enzyme for the digestion (by name)")
        .def("peptideCount", [](OpenMS::ProteaseDigestion& self, const OpenMS::AASequence& protein) { return self.peptideCount(protein); }, "protein"_a, "Returns the number of peptides a digestion of protein would yield under the current enzyme and missed cleavage settings")
        .def("getMissedCleavages", [](const OpenMS::ProteaseDigestion& self) { return self.getMissedCleavages(); }, "Returns the max. number of allowed missed cleavages for the digestion")
        .def("setMissedCleavages", [](OpenMS::ProteaseDigestion& self, size_t missed_cleavages) { return self.setMissedCleavages(missed_cleavages); }, "missed_cleavages"_a, "Sets the max. number of allowed missed cleavages for the digestion (default is 0). This setting is ignored when log model is used")
        .def("getEnzymeName", [](const OpenMS::ProteaseDigestion& self) { return self.getEnzymeName(); }, "Returns the enzyme for the digestion")
        .def("getSpecificity", [](const OpenMS::ProteaseDigestion& self) { return self.getSpecificity(); }, "Returns the specificity for the digestion")
        .def("setSpecificity", [](OpenMS::ProteaseDigestion& self, OpenMS::EnzymaticDigestion::Specificity spec) { return self.setSpecificity(spec); }, "spec"_a, "Sets the specificity for the digestion (default is SPEC_FULL)")
        .def_static("getSpecificityByName", [](const OpenMS::String& name) { return OpenMS::ProteaseDigestion::getSpecificityByName(name); }, "name"_a, "Returns the specificity by name. Returns SPEC_UNKNOWN if name is not valid")
        .def("countInternalCleavageSites", [](const OpenMS::ProteaseDigestion& self, const OpenMS::String& sequence) { return self.countInternalCleavageSites(sequence); }, "sequence"_a, "Returns the number of internal cleavage sites for this sequence.")

        .def("digest", [](OpenMS::ProteaseDigestion& self, const OpenMS::AASequence& protein, nb::list output, size_t min_length, size_t max_length) {
            // Backward-compatible: appends results to the provided Python list (with length filters)
            std::vector<OpenMS::AASequence> result;
            self.digest(protein, result, min_length, max_length);
            for (const auto& seq : result) {
                output.append(nb::cast(seq));
            }
        }, "protein"_a, "output"_a, "min_length"_a, "max_length"_a, "Performs the enzymatic digestion of a protein sequence (appends to output list, with length filters)")
        .def("digest", [](OpenMS::ProteaseDigestion& self, const OpenMS::AASequence& protein, nb::list output) {
            // Backward-compatible: appends results to the provided Python list
            std::vector<OpenMS::AASequence> result;
            self.digest(protein, result);
            for (const auto& seq : result) {
                output.append(nb::cast(seq));
            }
        }, "protein"_a, "output"_a, "Performs the enzymatic digestion of a protein sequence (appends to output list)")
        .def("digest", [](OpenMS::ProteaseDigestion& self, const OpenMS::AASequence& protein, size_t min_length, size_t max_length) {
            std::vector<OpenMS::AASequence> output;
            self.digest(protein, output, min_length, max_length);
            return output;
        }, "protein"_a, "min_length"_a, "max_length"_a, "Performs the enzymatic digestion of a protein sequence with length filters, returns list of peptides")
        .def("digest", [](OpenMS::ProteaseDigestion& self, const OpenMS::AASequence& protein) {
            std::vector<OpenMS::AASequence> output;
            self.digest(protein, output);
            return output;
        }, "protein"_a, "Performs the enzymatic digestion of a protein sequence, returns list of peptides")

        .def("isValidProduct", nb::overload_cast<const OpenMS::AASequence&, int, int, bool, bool, bool>(&OpenMS::ProteaseDigestion::isValidProduct, nb::const_),
            "protein"_a, "pep_pos"_a, "pep_length"_a, "ignore_missed_cleavages"_a = true, "allow_nterm_protein_cleavage"_a = false, "allow_random_asp_pro_cleavage"_a = false,
            "Check if peptide is a valid digestion product of protein (AASequence version)")
        .def("isValidProduct", nb::overload_cast<const OpenMS::String&, int, int, bool, bool, bool>(&OpenMS::ProteaseDigestion::isValidProduct, nb::const_),
            "protein"_a, "pep_pos"_a, "pep_length"_a, "ignore_missed_cleavages"_a = true, "allow_nterm_protein_cleavage"_a = false, "allow_random_asp_pro_cleavage"_a = false,
            "Check if peptide is a valid digestion product of protein (String version)")
        .def("digestUnmodified", [](const OpenMS::ProteaseDigestion& self, const std::string& sequence_str, size_t min_length, size_t max_length) {
            OpenMS::StringView sequence(sequence_str);
            std::vector<OpenMS::StringView> output;
            OpenMS::Size discarded = self.digestUnmodified(sequence, output, min_length, max_length);
            std::vector<std::string> result;
            result.reserve(output.size());
            for (const auto& sv : output) result.push_back(std::string(sv.getString()));
            return nb::make_tuple(result, discarded);
        }, "sequence"_a, "min_length"_a = 1, "max_length"_a = 0, "Digest unmodified sequence, returns (products, num_discarded)")
        ;

    // -----------------------------------------------------------------------
    // RNaseDB
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::RNaseDB>(m, "RNaseDB", 
        R"doc(
Database of enzymes that digest RNA (RNases). This is a singleton class.
The enzymes are read from share/CHEMISTRY/Enzymes_RNA.xml.
)doc")

        .def("hasEnzyme", [](const OpenMS::RNaseDB& self, const OpenMS::String& name) {
            return self.hasEnzyme(name);
        }, "name"_a, "Check if an enzyme with the given name exists")

        .def("getEnzyme", [](const OpenMS::RNaseDB& self, const OpenMS::String& name) -> const OpenMS::DigestionEnzymeRNA* {
            return self.getEnzyme(name);
        }, "name"_a, nb::rv_policy::reference, "Get the enzyme with the given name")

        .def("getEnzymeByRegEx", [](const OpenMS::RNaseDB& self, const OpenMS::String& cleavage_regex) -> const OpenMS::DigestionEnzymeRNA* {
            return self.getEnzymeByRegEx(cleavage_regex);
        }, "cleavage_regex"_a, nb::rv_policy::reference, "Get the enzyme with the given cleavage regex")

        .def("getAllNames", [](const OpenMS::RNaseDB& self) {
            std::vector<OpenMS::String> names;
            self.getAllNames(names);
            nb::list result;
            for (const auto& n : names) result.append(nb::str(n.c_str()));
            return result;
        }, "Get all enzyme names")
        .def_static("getInstance", []() -> OpenMS::RNaseDB* { return OpenMS::RNaseDB::getInstance(); }, nb::rv_policy::reference, "Returns the singleton instance")
        ;

    // -----------------------------------------------------------------------
    // RNaseDigestion
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::RNaseDigestion, OpenMS::EnzymaticDigestion>(m, "RNaseDigestion", 
        R"doc(
EnzymaticDigestion

Class for the enzymatic digestion of RNA
Usage:
.. code-block:: python
from pyopenms import *
oligo = NASequence.fromString("pAUGUCGCAG");
dig = RNaseDigestion()
dig.setEnzyme("RNase_T1")
result = []
dig.digest(oligo, result)
for fragment in result:
print (fragment)
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::RNaseDigestion &>())
        .def("__copy__", [](const OpenMS::RNaseDigestion& self) { return OpenMS::RNaseDigestion(self); })
        .def("__deepcopy__", [](const OpenMS::RNaseDigestion& self, nb::dict) { return OpenMS::RNaseDigestion(self); }, "memo"_a)
        .def("setEnzyme", [](OpenMS::RNaseDigestion& self, OpenMS::DigestionEnzyme * enzyme) { return self.setEnzyme(enzyme); }, "enzyme"_a, "Sets the enzyme for the digestion (by name)")
        .def("setEnzyme", [](OpenMS::RNaseDigestion& self, const OpenMS::String& name) { return self.setEnzyme(name); }, "name"_a, "Sets the enzyme for the digestion (by name)")
        .def("getMissedCleavages", [](const OpenMS::RNaseDigestion& self) { return self.getMissedCleavages(); }, "Returns the max. number of allowed missed cleavages for the digestion")
        .def("setMissedCleavages", [](OpenMS::RNaseDigestion& self, size_t missed_cleavages) { return self.setMissedCleavages(missed_cleavages); }, "missed_cleavages"_a, "Sets the max. number of allowed missed cleavages for the digestion (default is 0). This setting is ignored when log model is used")
        .def("getEnzymeName", [](const OpenMS::RNaseDigestion& self) { return self.getEnzymeName(); }, "Returns the enzyme for the digestion")
        .def("getSpecificity", [](const OpenMS::RNaseDigestion& self) { return self.getSpecificity(); }, "Returns the specificity for the digestion")
        .def("setSpecificity", [](OpenMS::RNaseDigestion& self, OpenMS::EnzymaticDigestion::Specificity spec) { return self.setSpecificity(spec); }, "spec"_a, "Sets the specificity for the digestion (default is SPEC_FULL)")
        .def_static("getSpecificityByName", [](const OpenMS::String& name) { return OpenMS::RNaseDigestion::getSpecificityByName(name); }, "name"_a, "Returns the specificity by name. Returns SPEC_UNKNOWN if name is not valid")
        .def("isValidProduct", [](const OpenMS::RNaseDigestion& self, const OpenMS::String& protein, int pep_pos, int pep_length, bool ignore_missed_cleavages) { return self.isValidProduct(protein, pep_pos, pep_length, ignore_missed_cleavages); }, "protein"_a, "pep_pos"_a, "pep_length"_a, "ignore_missed_cleavages"_a, 
            R"doc(
Performs the enzymatic digestion of an unmodified sequence\n
By returning only references into the original string this is very fast
:param sequence: Sequence to digest
:param output: Digestion products
:param min_length: Minimal length of reported products
:param max_length: Maximal length of reported products (0 = no restriction)
:return: Number of discarded digestion products (which are not matching length restrictions)
)doc")
        .def("countInternalCleavageSites", [](const OpenMS::RNaseDigestion& self, const OpenMS::String& sequence) { return self.countInternalCleavageSites(sequence); }, "sequence"_a, "Returns the number of internal cleavage sites for this sequence.")
        .def("digest", [](const OpenMS::RNaseDigestion& self, const OpenMS::NASequence& rna, OpenMS::Size min_length, OpenMS::Size max_length) {
            std::vector<OpenMS::NASequence> output;
            self.digest(rna, output, min_length, max_length);
            return output;
        }, "rna"_a, "min_length"_a = 0, "max_length"_a = 0, "Digest an RNA sequence and return the fragments")
        .def("digestUnmodified", [](const OpenMS::RNaseDigestion& self, const std::string& sequence_str, size_t min_length, size_t max_length) {
            OpenMS::StringView sequence(sequence_str);
            std::vector<OpenMS::StringView> output;
            OpenMS::Size discarded = self.digestUnmodified(sequence, output, min_length, max_length);
            std::vector<std::string> result;
            result.reserve(output.size());
            for (const auto& sv : output) result.push_back(std::string(sv.getString()));
            return nb::make_tuple(result, discarded);
        }, "sequence"_a, "min_length"_a = 1, "max_length"_a = 0, "Digest unmodified sequence, returns (products, num_discarded)")
        ;

    // -----------------------------------------------------------------------
    // RealMassDecomposer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ims::RealMassDecomposer>(m, "RealMassDecomposer", 
        R"doc(
Handles decomposing of non-integer values/masses over a set of
non-integer weights with an error allowed
)doc")
        .def(nb::init<OpenMS::ims::Weights>())
        .def("getNumberOfDecompositions", [](OpenMS::ims::RealMassDecomposer& self, double mass, double error) { return self.getNumberOfDecompositions(mass, error); }, "mass"_a, "error"_a)
        ;

    // -----------------------------------------------------------------------
    // Residue
    // -----------------------------------------------------------------------
    auto residue_class = nb::class_<OpenMS::Residue>(m, "Residue", "Representation of an amino acid residue")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Residue &>())
        .def("__copy__", [](const OpenMS::Residue& self) { return OpenMS::Residue(self); })
        .def("__deepcopy__", [](const OpenMS::Residue& self, nb::dict) { return OpenMS::Residue(self); }, "memo"_a)
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::String, OpenMS::EmpiricalFormula, double, double, double, double, double, double, std::set<OpenMS::String>>())
        .def_static("getInternalToFull", []() { return OpenMS::Residue::getInternalToFull(); })
        .def_static("getInternalToNTerm", []() { return OpenMS::Residue::getInternalToNTerm(); })
        .def_static("getInternalToCTerm", []() { return OpenMS::Residue::getInternalToCTerm(); })
        .def_static("getInternalToAIon", []() { return OpenMS::Residue::getInternalToAIon(); })
        .def_static("getInternalToBIon", []() { return OpenMS::Residue::getInternalToBIon(); })
        .def_static("getInternalToCIon", []() { return OpenMS::Residue::getInternalToCIon(); })
        .def_static("getInternalToXIon", []() { return OpenMS::Residue::getInternalToXIon(); })
        .def_static("getInternalToYIon", []() { return OpenMS::Residue::getInternalToYIon(); })
        .def_static("getInternalToZIon", []() { return OpenMS::Residue::getInternalToZIon(); })
        .def_static("getResidueTypeName", [](OpenMS::Residue::ResidueType res_type) { return OpenMS::Residue::getResidueTypeName(res_type); }, "res_type"_a, "Returns the ion name given as a residue type")
        .def("setName", [](OpenMS::Residue& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the residue")
        .def("getName", [](const OpenMS::Residue& self) { return self.getName(); }, "Returns the name of the residue")
        .def("setSynonyms", [](OpenMS::Residue& self, const std::set<OpenMS::String>& synonyms) { return self.setSynonyms(synonyms); }, "synonyms"_a, "Sets the synonyms")
        .def("addSynonym", [](OpenMS::Residue& self, const OpenMS::String& synonym) { return self.addSynonym(synonym); }, "synonym"_a, "Adds a synonym")
        .def("getSynonyms", [](const OpenMS::Residue& self) -> const std::set<OpenMS::String> & { return self.getSynonyms(); }, nb::rv_policy::reference_internal, "Returns the sysnonyms")
        .def("setThreeLetterCode", [](OpenMS::Residue& self, const OpenMS::String& three_letter_code) { return self.setThreeLetterCode(three_letter_code); }, "three_letter_code"_a, "Sets the name of the residue as three letter code")
        .def("getThreeLetterCode", [](const OpenMS::Residue& self) { return self.getThreeLetterCode(); }, "Returns the name of the residue as three letter code")
        .def("setOneLetterCode", [](OpenMS::Residue& self, const OpenMS::String& one_letter_code) { return self.setOneLetterCode(one_letter_code); }, "one_letter_code"_a, "Sets the name as one letter code")
        .def("getOneLetterCode", [](const OpenMS::Residue& self) { return self.getOneLetterCode(); }, "Returns the name as one letter code")
        .def("addLossFormula", [](OpenMS::Residue& self, const OpenMS::EmpiricalFormula& p0) { return self.addLossFormula(p0); }, "Adds a neutral loss formula")
        .def("setLossFormulas", [](OpenMS::Residue& self, const std::vector<OpenMS::EmpiricalFormula>& p0) { return self.setLossFormulas(p0); }, "Sets the neutral loss formulas")
        .def("addNTermLossFormula", [](OpenMS::Residue& self, const OpenMS::EmpiricalFormula& p0) { return self.addNTermLossFormula(p0); }, "Adds N-terminal losses")
        .def("setNTermLossFormulas", [](OpenMS::Residue& self, const std::vector<OpenMS::EmpiricalFormula>& p0) { return self.setNTermLossFormulas(p0); }, "Sets the N-terminal losses")
        .def("getLossFormulas", [](const OpenMS::Residue& self) -> const std::vector<OpenMS::EmpiricalFormula> & { return self.getLossFormulas(); }, nb::rv_policy::reference_internal, "Returns the neutral loss formulas")
        .def("getNTermLossFormulas", [](const OpenMS::Residue& self) -> const std::vector<OpenMS::EmpiricalFormula> & { return self.getNTermLossFormulas(); }, nb::rv_policy::reference_internal, "Returns N-terminal loss formulas")
        .def("setLossNames", [](OpenMS::Residue& self, const std::vector<OpenMS::String>& name) { return self.setLossNames(name); }, "name"_a, "Sets the neutral loss molecule name")
        .def("setNTermLossNames", [](OpenMS::Residue& self, const std::vector<OpenMS::String>& name) { return self.setNTermLossNames(name); }, "name"_a, "Sets the N-terminal loss names")
        .def("addLossName", [](OpenMS::Residue& self, const OpenMS::String& name) { return self.addLossName(name); }, "name"_a, "Adds neutral loss molecule name")
        .def("addNTermLossName", [](OpenMS::Residue& self, const OpenMS::String& name) { return self.addNTermLossName(name); }, "name"_a, "Adds a N-terminal loss name")
        .def("getLossNames", [](const OpenMS::Residue& self) -> const std::vector<OpenMS::String> & { return self.getLossNames(); }, nb::rv_policy::reference_internal, "Gets neutral loss name (if there is one, else returns an empty string)")
        .def("getNTermLossNames", [](const OpenMS::Residue& self) -> const std::vector<OpenMS::String> & { return self.getNTermLossNames(); }, nb::rv_policy::reference_internal, "Returns the N-terminal loss names")
        .def("setFormula", [](OpenMS::Residue& self, const OpenMS::EmpiricalFormula& formula) { return self.setFormula(formula); }, "formula"_a, "Sets empirical formula of the residue (must be full, with N and C-terminus)")
        .def("getFormula", [](const OpenMS::Residue& self, OpenMS::Residue::ResidueType res_type) { return self.getFormula(res_type); }, "res_type"_a)
        .def("setAverageWeight", [](OpenMS::Residue& self, double weight) { return self.setAverageWeight(weight); }, "weight"_a, "Sets average weight of the residue (must be full, with N and C-terminus)")
        .def("getAverageWeight", [](const OpenMS::Residue& self, OpenMS::Residue::ResidueType res_type) { return self.getAverageWeight(res_type); }, "res_type"_a)
        .def("setMonoWeight", [](OpenMS::Residue& self, double weight) { return self.setMonoWeight(weight); }, "weight"_a, "Sets monoisotopic weight of the residue (must be full, with N and C-terminus)")
        .def("getMonoWeight", [](const OpenMS::Residue& self, OpenMS::Residue::ResidueType res_type) { return self.getMonoWeight(res_type); }, "res_type"_a)
        .def("getModification", [](const OpenMS::Residue& self) { return self.getModification(); }, nb::rv_policy::reference_internal)
        .def("setModification", [](OpenMS::Residue& self, const OpenMS::String& name) { return self.setModification(name); }, "name"_a, "Sets the modification by name; the mod should be present in ModificationsDB")
        .def("setModification", [](OpenMS::Residue& self, OpenMS::ResidueModification * mod) { return self.setModification(mod); }, "mod"_a, "Sets the modification by name; the mod should be present in ModificationsDB")
        .def("setModification", [](OpenMS::Residue& self, const OpenMS::ResidueModification& mod) { return self.setModification(mod); }, "mod"_a, "Sets the modification by name; the mod should be present in ModificationsDB")
        .def("setModificationByDiffMonoMass", [](OpenMS::Residue& self, double diffMonoMass) { return self.setModificationByDiffMonoMass(diffMonoMass); }, "diffMonoMass"_a, 
            R"doc(
Sets the modification by monoisotopic mass difference in Da; checks if present in ModificationsDB with tolerance and adds a "user-defined" modification if not (for later lookups).
)doc")
        .def("getModificationName", [](const OpenMS::Residue& self) { return self.getModificationName(); }, "Returns the name of the modification to the modification")
        .def("setLowMassIons", [](OpenMS::Residue& self, const std::vector<OpenMS::EmpiricalFormula>& low_mass_ions) { return self.setLowMassIons(low_mass_ions); }, "low_mass_ions"_a, "Sets the low mass marker ions as a vector of formulas")
        .def("getLowMassIons", [](const OpenMS::Residue& self) -> const std::vector<OpenMS::EmpiricalFormula> & { return self.getLowMassIons(); }, nb::rv_policy::reference_internal, "Returns a vector of formulas with the low mass markers of the residue")
        .def("setResidueSets", [](OpenMS::Residue& self, const std::set<OpenMS::String>& residues_sets) { return self.setResidueSets(residues_sets); }, "residues_sets"_a, "Sets the residue sets the amino acid is contained in")
        .def("addResidueSet", [](OpenMS::Residue& self, const OpenMS::String& residue_sets) { return self.addResidueSet(residue_sets); }, "residue_sets"_a, "Adds a residue set to the residue sets")
        .def("getResidueSets", [](const OpenMS::Residue& self) -> const std::set<OpenMS::String> & { return self.getResidueSets(); }, nb::rv_policy::reference_internal, "Returns the residue sets this residue is contained in")
        .def("getPka", [](const OpenMS::Residue& self) { return self.getPka(); }, "Returns the pka of the residue")
        .def("getPkb", [](const OpenMS::Residue& self) { return self.getPkb(); }, "Returns the pkb of the residue")
        .def("getPkc", [](const OpenMS::Residue& self) { return self.getPkc(); }, "Returns the pkc of the residue if it exists otherwise -1")
        .def("getPiValue", [](const OpenMS::Residue& self) { return self.getPiValue(); }, "Calculates the isoelectric point using the pk values")
        .def("setPka", [](OpenMS::Residue& self, double value) { return self.setPka(value); }, "value"_a, "Sets the pka of the residue")
        .def("setPkb", [](OpenMS::Residue& self, double value) { return self.setPkb(value); }, "value"_a, "Sets the pkb of the residue")
        .def("setPkc", [](OpenMS::Residue& self, double value) { return self.setPkc(value); }, "value"_a, "Sets the pkc of the residue")
        .def("getSideChainBasicity", [](const OpenMS::Residue& self) { return self.getSideChainBasicity(); }, "Returns the side chain basicity")
        .def("setSideChainBasicity", [](OpenMS::Residue& self, double gb_sc) { return self.setSideChainBasicity(gb_sc); }, "gb_sc"_a, "Sets the side chain basicity")
        .def("getBackboneBasicityLeft", [](const OpenMS::Residue& self) { return self.getBackboneBasicityLeft(); }, "Returns the backbone basicitiy if located in N-terminal direction")
        .def("setBackboneBasicityLeft", [](OpenMS::Residue& self, double gb_bb_l) { return self.setBackboneBasicityLeft(gb_bb_l); }, "gb_bb_l"_a, "Sets the N-terminal direction backbone basicitiy")
        .def("getBackboneBasicityRight", [](const OpenMS::Residue& self) { return self.getBackboneBasicityRight(); }, "Returns the C-terminal direction backbone basicitiy")
        .def("setBackboneBasicityRight", [](OpenMS::Residue& self, double gb_bb_r) { return self.setBackboneBasicityRight(gb_bb_r); }, "gb_bb_r"_a, "Sets the C-terminal direction backbone basicity")
        .def("hasNeutralLoss", [](const OpenMS::Residue& self) { return self.hasNeutralLoss(); }, "True if the residue has neutral loss")
        .def("hasNTermNeutralLosses", [](const OpenMS::Residue& self) { return self.hasNTermNeutralLosses(); }, "True if N-terminal neutral losses are set")
        .def("getHydrophobicity", [](const OpenMS::Residue& self, OpenMS::HydrophobicityScaleMethod scale) { return self.getHydrophobicity(scale); }, "scale"_a, "Returns the hydrophobicity value of the residue for the given scale (throws for non-standard residues)")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("isModified", [](const OpenMS::Residue& self) { return self.isModified(); }, "True if the residue is a modified one")
        .def("isInResidueSet", [](OpenMS::Residue& self, const OpenMS::String& residue_set) { return self.isInResidueSet(residue_set); }, "residue_set"_a, "True if the residue is contained in the set")
        .def_static("residueTypeToIonLetter", [](const OpenMS::Residue::ResidueType& res_type) { return OpenMS::Residue::residueTypeToIonLetter(res_type); }, "res_type"_a, "Helper for mapping residue types to letters for Text annotations and labels")
        .def("toString", [](const OpenMS::Residue& self) { return self.toString(); }, "Returns the residue as string (one letter code with optional modification)")
        .def("__hash__", [](const OpenMS::Residue& self) { return std::hash<OpenMS::Residue>{}(self); })
        .def("__repr__", [](const OpenMS::Residue& self) {
            std::ostringstream oss;
            oss << "Residue(name='" << std::string(self.getName())
                << "', one_letter='" << std::string(self.getOneLetterCode())
                << "', three_letter='" << std::string(self.getThreeLetterCode())
                << "', mono_mass=" << self.getMonoWeight(OpenMS::Residue::ResidueType::Full) << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::Residue& self) { return std::string(self.getOneLetterCode()); })
        ;
    // ResidueType enum nested under Residue
    nb::enum_<OpenMS::Residue::ResidueType>(residue_class, "ResidueType", nb::is_arithmetic())
        .value("Full", OpenMS::Residue::ResidueType::Full)
        .value("Internal", OpenMS::Residue::ResidueType::Internal)
        .value("NTerminal", OpenMS::Residue::ResidueType::NTerminal)
        .value("CTerminal", OpenMS::Residue::ResidueType::CTerminal)
        .value("AIon", OpenMS::Residue::ResidueType::AIon)
        .value("BIon", OpenMS::Residue::ResidueType::BIon)
        .value("CIon", OpenMS::Residue::ResidueType::CIon)
        .value("XIon", OpenMS::Residue::ResidueType::XIon)
        .value("YIon", OpenMS::Residue::ResidueType::YIon)
        .value("ZIon", OpenMS::Residue::ResidueType::ZIon)
        .value("Zp1Ion", OpenMS::Residue::ResidueType::Zp1Ion)
        .value("Zp2Ion", OpenMS::Residue::ResidueType::Zp2Ion)
        .value("Precursor", OpenMS::Residue::ResidueType::Precursor)
        .value("BIonMinusH20", OpenMS::Residue::ResidueType::BIonMinusH20)
        .value("YIonMinusH20", OpenMS::Residue::ResidueType::YIonMinusH20)
        .value("BIonMinusNH3", OpenMS::Residue::ResidueType::BIonMinusNH3)
        .value("YIonMinusNH3", OpenMS::Residue::ResidueType::YIonMinusNH3)
        .value("NonIdentified", OpenMS::Residue::ResidueType::NonIdentified)
        .value("Unannotated", OpenMS::Residue::ResidueType::Unannotated)
        .value("SizeOfResidueType", OpenMS::Residue::ResidueType::SizeOfResidueType)
        .export_values();

    // HydrophobicityScaleMethod enum (namespace-scoped, used by Residue::getHydrophobicity)
    nb::enum_<OpenMS::HydrophobicityScaleMethod>(m, "HydrophobicityScaleMethod", nb::is_arithmetic())
        .value("KYTE_DOOLITTLE", OpenMS::HydrophobicityScaleMethod::KYTE_DOOLITTLE)
        .value("EISENBERG", OpenMS::HydrophobicityScaleMethod::EISENBERG)
        .value("HOPP_WOODS", OpenMS::HydrophobicityScaleMethod::HOPP_WOODS)
        .value("BULL_BREESE", OpenMS::HydrophobicityScaleMethod::BULL_BREESE)
        .value("BLACK_MOULD", OpenMS::HydrophobicityScaleMethod::BLACK_MOULD)
        .value("GUY", OpenMS::HydrophobicityScaleMethod::GUY)
        .value("EISENBERG_CONSENSUS", OpenMS::HydrophobicityScaleMethod::EISENBERG_CONSENSUS)
        .export_values();

    // ProteomicsPkaScale enum (namespace-scoped, used by IsoelectricPoint)
    nb::enum_<OpenMS::ProteomicsPkaScale>(m, "ProteomicsPkaScale", nb::is_arithmetic())
        .value("LEHNINGER", OpenMS::ProteomicsPkaScale::LEHNINGER)
        .value("EMBOSS", OpenMS::ProteomicsPkaScale::EMBOSS)
        .value("SILLERO", OpenMS::ProteomicsPkaScale::SILLERO)
        .value("BJELLQVIST", OpenMS::ProteomicsPkaScale::BJELLQVIST)
        .export_values();

    // IsoelectricPoint
    auto isoelectricpoint_class = nb::class_<OpenMS::IsoelectricPoint>(m, "IsoelectricPoint",
        "Utility class for computing isoelectric point (pI) and net charge of peptides")
        .def_static("computeCharge",
            [](const OpenMS::AASequence& seq, double pH, OpenMS::ProteomicsPkaScale scale) {
                return OpenMS::IsoelectricPoint::computeCharge(seq, pH, scale);
            },
            "seq"_a, "pH"_a, "scale"_a = OpenMS::ProteomicsPkaScale::LEHNINGER,
            "Computes the net charge of an amino acid sequence at a given pH")
        .def_static("computePI",
            [](const OpenMS::AASequence& seq, OpenMS::ProteomicsPkaScale scale, double tolerance) {
                return OpenMS::IsoelectricPoint::computePI(seq, scale, tolerance);
            },
            "seq"_a, "scale"_a = OpenMS::ProteomicsPkaScale::LEHNINGER, "tolerance"_a = 1e-4,
            "Computes the isoelectric point (pI) of an amino acid sequence via bisection")
        ;

    // -----------------------------------------------------------------------
    // ResidueDB
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ResidueDB>(m, "ResidueDB", 
        R"doc(
Database of residues (amino acids). This is a singleton class.
All unmodified residues are added to the database on construction.
Modified residues get created and added if getModifiedResidue is called.
)doc")
        .def("getNumberOfResidues", [](const OpenMS::ResidueDB& self) { return self.getNumberOfResidues(); }, "Returns the number of residues stored")
        .def("getNumberOfModifiedResidues", [](const OpenMS::ResidueDB& self) { return self.getNumberOfModifiedResidues(); }, "Returns the number of modified residues stored")
        .def("getResidue", [](const OpenMS::ResidueDB& self, const OpenMS::String& name) { return self.getResidue(name); }, "name"_a, nb::rv_policy::reference, "Returns a pointer to the residue with name, 3 letter code or 1 letter code name")
        .def("getResidue", [](const OpenMS::ResidueDB& self, const unsigned char& one_letter_code) { return self.getResidue(one_letter_code); }, "one_letter_code"_a, nb::rv_policy::reference, "Returns a pointer to the residue with name, 3 letter code or 1 letter code name")
        .def("getModifiedResidue", [](OpenMS::ResidueDB& self, const OpenMS::String& name) { return self.getModifiedResidue(name); }, "name"_a, nb::rv_policy::reference, "Returns a pointer to a modified residue given a modification name")
        .def("getModifiedResidue", [](OpenMS::ResidueDB& self, OpenMS::Residue * residue, const OpenMS::String& name) { return self.getModifiedResidue(residue, name); }, "residue"_a, "name"_a, nb::rv_policy::reference, "Returns a pointer to a modified residue given a residue and a modification name")
        .def("getModifiedResidue", [](OpenMS::ResidueDB& self, OpenMS::Residue * residue, OpenMS::ResidueModification * mod) { return self.getModifiedResidue(residue, mod); }, "residue"_a, "mod"_a, nb::rv_policy::reference, "Returns a pointer to a modified residue given a residue and a modification name")
        .def("getResidues", [](const OpenMS::ResidueDB& self, const OpenMS::String& residue_set) { return self.getResidues(residue_set); }, "residue_set"_a = "All", nb::rv_policy::reference, "Returns a set of all residues stored in this residue db")
        .def("getResidueSets", [](const OpenMS::ResidueDB& self) { return self.getResidueSets(); }, "Returns all residue sets that are registered which this instance")
        .def("hasResidue", [](const OpenMS::ResidueDB& self, const OpenMS::String& name) { return self.hasResidue(name); }, "name"_a, "Returns True if the db contains a residue with the given name")
        .def("hasResidue", [](const OpenMS::ResidueDB& self, OpenMS::Residue * residue) { return self.hasResidue(residue); }, "residue"_a, "Returns True if the db contains a residue with the given name")
        .def_static("getInstance", []() -> OpenMS::ResidueDB* { return OpenMS::ResidueDB::getInstance(); }, nb::rv_policy::reference, "Returns the singleton instance")
        ;

    // -----------------------------------------------------------------------
    // ResidueModification
    // -----------------------------------------------------------------------
    auto residuemodification_class = nb::class_<OpenMS::ResidueModification>(m, "ResidueModification", "Representation of a modification on an amino acid residue")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ResidueModification &>())
        .def("__copy__", [](const OpenMS::ResidueModification& self) { return OpenMS::ResidueModification(self); })
        .def("__deepcopy__", [](const OpenMS::ResidueModification& self, nb::dict) { return OpenMS::ResidueModification(self); }, "memo"_a)
        .def("setId", [](OpenMS::ResidueModification& self, const OpenMS::String& id) { return self.setId(id); }, "id"_a, "Sets the identifier of the modification")
        .def("getId", [](const OpenMS::ResidueModification& self) { return self.getId(); }, "Returns the identifier of the modification")
        .def("setFullId", [](OpenMS::ResidueModification& self, const OpenMS::String& full_id) { return self.setFullId(full_id); }, "full_id"_a = "", "Sets the full identifier (Unimod Accession + origin, if available)")
        .def("getFullId", [](const OpenMS::ResidueModification& self) { return self.getFullId(); })
        .def("setUniModRecordId", [](OpenMS::ResidueModification& self, const int& id) { return self.setUniModRecordId(id); }, "id"_a, "Sets the unimod record id")
        .def("getUniModRecordId", [](const OpenMS::ResidueModification& self) { return self.getUniModRecordId(); }, "Gets the unimod record id")
        .def("getUniModAccession", [](const OpenMS::ResidueModification& self) { return self.getUniModAccession(); }, "Returns the unimod accession if available")
        .def("setPSIMODAccession", [](OpenMS::ResidueModification& self, const OpenMS::String& id) { return self.setPSIMODAccession(id); }, "id"_a, "Sets the MOD-XXXXX accession of PSI-MOD")
        .def("getPSIMODAccession", [](const OpenMS::ResidueModification& self) { return self.getPSIMODAccession(); }, "Returns the PSI-MOD accession if available")
        .def("setFullName", [](OpenMS::ResidueModification& self, const OpenMS::String& full_name) { return self.setFullName(full_name); }, "full_name"_a, "Sets the full name of the modification; must NOT contain the origin (or . for terminals!)")
        .def("getFullName", [](const OpenMS::ResidueModification& self) { return self.getFullName(); }, "Returns the full name of the modification")
        .def("setName", [](OpenMS::ResidueModification& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of modification")
        .def("getName", [](const OpenMS::ResidueModification& self) { return self.getName(); }, "Returns the PSI-MS-label if available; e.g. Mascot uses this name")
        .def("setTermSpecificity", [](OpenMS::ResidueModification& self, OpenMS::ResidueModification::TermSpecificity term_spec) { return self.setTermSpecificity(term_spec); }, "term_spec"_a, "Sets the term specificity")
        .def("setTermSpecificity", [](OpenMS::ResidueModification& self, const OpenMS::String& name) { return self.setTermSpecificity(name); }, "name"_a, "Sets the term specificity")
        .def("getTermSpecificity", [](const OpenMS::ResidueModification& self) { return self.getTermSpecificity(); }, "Returns terminal specificity")
        .def("getTermSpecificityName", [](const OpenMS::ResidueModification& self, OpenMS::ResidueModification::TermSpecificity term_spec) { return self.getTermSpecificityName(term_spec); }, "term_spec"_a, "Returns the name of the terminal specificity")
        .def("setOrigin", [](OpenMS::ResidueModification& self, char origin) { return self.setOrigin(origin); }, "origin"_a, "Sets the origin (i.e. modified amino acid)")
        .def("getOrigin", [](const OpenMS::ResidueModification& self) { return self.getOrigin(); }, "Returns the origin (i.e. modified amino acid)")
        .def("setSourceClassification", [](OpenMS::ResidueModification& self, const OpenMS::String& classification) { return self.setSourceClassification(classification); }, "classification"_a, "Classification as defined by the PSI-MOD")
        .def("setSourceClassification", [](OpenMS::ResidueModification& self, OpenMS::ResidueModification::SourceClassification classification) { return self.setSourceClassification(classification); }, "classification"_a, "Classification as defined by the PSI-MOD")
        .def("getSourceClassification", [](const OpenMS::ResidueModification& self) { return self.getSourceClassification(); }, "Returns the source classification, if none was set, it is unspecific")
        .def("getSourceClassificationName", [](const OpenMS::ResidueModification& self, OpenMS::ResidueModification::SourceClassification classification) { return self.getSourceClassificationName(classification); }, "classification"_a, "Returns the classification")
        .def("setAverageMass", [](OpenMS::ResidueModification& self, double mass) { return self.setAverageMass(mass); }, "mass"_a, "Sets the average mass")
        .def("getAverageMass", [](const OpenMS::ResidueModification& self) { return self.getAverageMass(); }, "Returns the average mass if set")
        .def("setMonoMass", [](OpenMS::ResidueModification& self, double mass) { return self.setMonoMass(mass); }, "mass"_a, "Sets the monoisotopic mass (this must include the weight of the residue itself!)")
        .def("getMonoMass", [](const OpenMS::ResidueModification& self) { return self.getMonoMass(); }, "Return the monoisotopic mass, or 0.0 if not set")
        .def("setDiffAverageMass", [](OpenMS::ResidueModification& self, double mass) { return self.setDiffAverageMass(mass); }, "mass"_a, "Sets the difference average mass")
        .def("getDiffAverageMass", [](const OpenMS::ResidueModification& self) { return self.getDiffAverageMass(); }, "Returns the difference average mass, or 0.0 if not set")
        .def("setDiffMonoMass", [](OpenMS::ResidueModification& self, double mass) { return self.setDiffMonoMass(mass); }, "mass"_a, "Sets the difference monoisotopic mass")
        .def("getDiffMonoMass", [](const OpenMS::ResidueModification& self) { return self.getDiffMonoMass(); }, "Returns the diff monoisotopic mass, or 0.0 if not set")
        .def("setFormula", [](OpenMS::ResidueModification& self, const OpenMS::String& composition) { return self.setFormula(composition); }, "composition"_a, "Sets the formula (no masses will be changed)")
        .def("getFormula", [](const OpenMS::ResidueModification& self) { return self.getFormula(); }, "Returns the chemical formula if set")
        .def("setDiffFormula", [](OpenMS::ResidueModification& self, const OpenMS::EmpiricalFormula& diff_formula) { return self.setDiffFormula(diff_formula); }, "diff_formula"_a, "Sets diff formula (no masses will be changed)")
        .def("getDiffFormula", [](const OpenMS::ResidueModification& self) -> const OpenMS::EmpiricalFormula& { return self.getDiffFormula(); }, nb::rv_policy::reference_internal, "Returns the diff formula if one was set")
        .def("setSynonyms", [](OpenMS::ResidueModification& self, const std::set<OpenMS::String>& synonyms) { return self.setSynonyms(synonyms); }, "synonyms"_a, "Sets the synonyms of that modification")
        .def("addSynonym", [](OpenMS::ResidueModification& self, const OpenMS::String& synonym) { return self.addSynonym(synonym); }, "synonym"_a, "Adds a synonym to the unique list")
        .def("getSynonyms", [](const OpenMS::ResidueModification& self) -> const std::set<OpenMS::String> & { return self.getSynonyms(); }, nb::rv_policy::reference_internal, "Returns the set of synonyms")
        .def("setNeutralLossDiffFormulas", [](OpenMS::ResidueModification& self, const std::vector<OpenMS::EmpiricalFormula>& diff_formulas) { return self.setNeutralLossDiffFormulas(diff_formulas); }, "diff_formulas"_a, "Sets the neutral loss formula")
        .def("getNeutralLossDiffFormulas", [](const OpenMS::ResidueModification& self) -> const std::vector<OpenMS::EmpiricalFormula> & { return self.getNeutralLossDiffFormulas(); }, nb::rv_policy::reference_internal, "Returns the neutral loss diff formula (if available)")
        .def("setNeutralLossMonoMasses", [](OpenMS::ResidueModification& self, std::vector<double> mono_masses) { return self.setNeutralLossMonoMasses(mono_masses); }, "mono_masses"_a, "Sets the neutral loss mono weight")
        .def("getNeutralLossMonoMasses", [](const OpenMS::ResidueModification& self) { return self.getNeutralLossMonoMasses(); }, "Returns the neutral loss mono weight")
        .def("setNeutralLossAverageMasses", [](OpenMS::ResidueModification& self, std::vector<double> average_masses) { return self.setNeutralLossAverageMasses(average_masses); }, "average_masses"_a, "Sets the neutral loss average weight")
        .def("getNeutralLossAverageMasses", [](const OpenMS::ResidueModification& self) { return self.getNeutralLossAverageMasses(); }, "Returns the neutral loss average weight")
        .def("hasNeutralLoss", [](const OpenMS::ResidueModification& self) { return self.hasNeutralLoss(); }, "Returns true if a neutral loss formula is set")
        .def("isUserDefined", [](const OpenMS::ResidueModification& self) { return self.isUserDefined(); }, "Returns true if it is a user-defined modification (empty id)")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("__hash__", [](const OpenMS::ResidueModification& self) { return std::hash<OpenMS::ResidueModification>{}(self); })
        .def("__repr__", [](const OpenMS::ResidueModification& self) {
            std::ostringstream oss;
            oss << "ResidueModification(id='" << std::string(self.getId())
                << "', full_name='" << std::string(self.getFullName())
                << "', mono_mass=" << self.getDiffMonoMass() << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::ResidueModification& self) { return nb::cast(self).attr("__repr__")(); })
        ;
    // TermSpecificity enum nested under ResidueModification
    nb::enum_<OpenMS::ResidueModification::TermSpecificity>(residuemodification_class, "TermSpecificity", nb::is_arithmetic())
        .value("ANYWHERE", OpenMS::ResidueModification::TermSpecificity::ANYWHERE)
        .value("C_TERM", OpenMS::ResidueModification::TermSpecificity::C_TERM)
        .value("N_TERM", OpenMS::ResidueModification::TermSpecificity::N_TERM)
        .value("PROTEIN_C_TERM", OpenMS::ResidueModification::TermSpecificity::PROTEIN_C_TERM)
        .value("PROTEIN_N_TERM", OpenMS::ResidueModification::TermSpecificity::PROTEIN_N_TERM)
        .value("NUMBER_OF_TERM_SPECIFICITY", OpenMS::ResidueModification::TermSpecificity::NUMBER_OF_TERM_SPECIFICITY)
        .export_values();
    // SourceClassification enum nested under ResidueModification
    nb::enum_<OpenMS::ResidueModification::SourceClassification>(residuemodification_class, "SourceClassification", nb::is_arithmetic())
        .value("ARTIFACT", OpenMS::ResidueModification::SourceClassification::ARTIFACT)
        .value("HYPOTHETICAL", OpenMS::ResidueModification::SourceClassification::HYPOTHETICAL)
        .value("NATURAL", OpenMS::ResidueModification::SourceClassification::NATURAL)
        .value("POSTTRANSLATIONAL", OpenMS::ResidueModification::SourceClassification::POSTTRANSLATIONAL)
        .value("MULTIPLE", OpenMS::ResidueModification::SourceClassification::MULTIPLE)
        .value("CHEMICAL_DERIVATIVE", OpenMS::ResidueModification::SourceClassification::CHEMICAL_DERIVATIVE)
        .value("ISOTOPIC_LABEL", OpenMS::ResidueModification::SourceClassification::ISOTOPIC_LABEL)
        .value("PRETRANSLATIONAL", OpenMS::ResidueModification::SourceClassification::PRETRANSLATIONAL)
        .value("OTHER_GLYCOSYLATION", OpenMS::ResidueModification::SourceClassification::OTHER_GLYCOSYLATION)
        .value("NLINKED_GLYCOSYLATION", OpenMS::ResidueModification::SourceClassification::NLINKED_GLYCOSYLATION)
        .value("AA_SUBSTITUTION", OpenMS::ResidueModification::SourceClassification::AA_SUBSTITUTION)
        .value("OTHER", OpenMS::ResidueModification::SourceClassification::OTHER)
        .value("NONSTANDARD_RESIDUE", OpenMS::ResidueModification::SourceClassification::NONSTANDARD_RESIDUE)
        .value("COTRANSLATIONAL", OpenMS::ResidueModification::SourceClassification::COTRANSLATIONAL)
        .value("OLINKED_GLYCOSYLATION", OpenMS::ResidueModification::SourceClassification::OLINKED_GLYCOSYLATION)
        .value("UNKNOWN", OpenMS::ResidueModification::SourceClassification::UNKNOWN)
        .value("NUMBER_OF_SOURCE_CLASSIFICATIONS", OpenMS::ResidueModification::SourceClassification::NUMBER_OF_SOURCE_CLASSIFICATIONS)
        .export_values();

    // -----------------------------------------------------------------------
    // Ribonucleotide
    // -----------------------------------------------------------------------
    auto ribonucleotide_class = nb::class_<OpenMS::Ribonucleotide>(m, "Ribonucleotide", "Ribonucleotide")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Ribonucleotide &>())
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::String, OpenMS::String, OpenMS::EmpiricalFormula, char, double, double, OpenMS::Ribonucleotide::TermSpecificityNuc, OpenMS::EmpiricalFormula>())
        .def("__copy__", [](const OpenMS::Ribonucleotide& self) { return OpenMS::Ribonucleotide(self); })
        .def("__deepcopy__", [](const OpenMS::Ribonucleotide& self, nb::dict) { return OpenMS::Ribonucleotide(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("getCode", [](const OpenMS::Ribonucleotide& self) { return self.getCode(); }, "Returns the short name")
        .def("setCode", [](OpenMS::Ribonucleotide& self, const OpenMS::String& code) { return self.setCode(code); }, "code"_a, "Sets the short name")
        .def("getName", [](const OpenMS::Ribonucleotide& self) { return self.getName(); }, "Returns the name of the ribonucleotide")
        .def("setName", [](OpenMS::Ribonucleotide& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the ribonucleotide")
        .def("getFormula", [](const OpenMS::Ribonucleotide& self) { return self.getFormula(); }, "Returns the empirical formula of the residue")
        .def("setFormula", [](OpenMS::Ribonucleotide& self, const OpenMS::EmpiricalFormula& formula) { return self.setFormula(formula); }, "formula"_a, "Sets empirical formula of the ribonucleotide (must be full, with N and C-terminus)")
        .def("getMonoMass", [](const OpenMS::Ribonucleotide& self) { return self.getMonoMass(); }, "Returns monoisotopic mass of the ribonucleotide")
        .def("setMonoMass", [](OpenMS::Ribonucleotide& self, double mono_mass) { return self.setMonoMass(mono_mass); }, "mono_mass"_a, "Sets monoisotopic mass of the ribonucleotide")
        .def("getAvgMass", [](const OpenMS::Ribonucleotide& self) { return self.getAvgMass(); }, "Returns average mass of the ribonucleotide")
        .def("setAvgMass", [](OpenMS::Ribonucleotide& self, double avg_mass) { return self.setAvgMass(avg_mass); }, "avg_mass"_a, "Sets average mass of the ribonucleotide")
        .def("getNewCode", [](const OpenMS::Ribonucleotide& self) { return self.getNewCode(); }, "Returns the new code")
        .def("setNewCode", [](OpenMS::Ribonucleotide& self, const OpenMS::String& new_code) { return self.setNewCode(new_code); }, "new_code"_a, "Sets the new code")
        .def("getOrigin", [](const OpenMS::Ribonucleotide& self) { return self.getOrigin(); }, 
            R"doc(
Returns the code of the unmodified base (e.g., "A", "C", ...)
)doc")
        .def("setOrigin", [](OpenMS::Ribonucleotide& self, char origin) { return self.setOrigin(origin); }, "origin"_a, 
            R"doc(
Sets the code of the unmodified base (e.g., "A", "C", ...)
)doc")
        .def("getHTMLCode", [](const OpenMS::Ribonucleotide& self) { return self.getHTMLCode(); }, "Returns the HTML (RNAMods) code")
        .def("setHTMLCode", [](OpenMS::Ribonucleotide& self, const OpenMS::String& html_code) { return self.setHTMLCode(html_code); }, "html_code"_a, "Sets the HTML (RNAMods) code")
        .def("getTermSpecificity", [](const OpenMS::Ribonucleotide& self) { return self.getTermSpecificity(); }, "Returns the terminal specificity")
        .def("setTermSpecificity", [](OpenMS::Ribonucleotide& self, OpenMS::Ribonucleotide::TermSpecificityNuc term_spec) { return self.setTermSpecificity(term_spec); }, "term_spec"_a, "Sets the terminal specificity")
        .def("getBaselossFormula", [](const OpenMS::Ribonucleotide& self) { return self.getBaselossFormula(); }, "Returns sum formula after loss of the nucleobase")
        .def("setBaselossFormula", [](OpenMS::Ribonucleotide& self, const OpenMS::EmpiricalFormula& formula) { return self.setBaselossFormula(formula); }, "formula"_a, "Sets sum formula after loss of the nucleobase")
        .def("isModified", [](const OpenMS::Ribonucleotide& self) { return self.isModified(); }, "True if the ribonucleotide is a modified one")
        .def("__hash__", [](const OpenMS::Ribonucleotide& self) { return std::hash<OpenMS::Ribonucleotide>{}(self); })
        ;
    // TermSpecificityNuc enum nested under Ribonucleotide
    nb::enum_<OpenMS::Ribonucleotide::TermSpecificityNuc>(ribonucleotide_class, "TermSpecificityNuc", nb::is_arithmetic())
        .value("ANYWHERE", OpenMS::Ribonucleotide::TermSpecificityNuc::ANYWHERE)
        .value("FIVE_PRIME", OpenMS::Ribonucleotide::TermSpecificityNuc::FIVE_PRIME)
        .value("THREE_PRIME", OpenMS::Ribonucleotide::TermSpecificityNuc::THREE_PRIME)
        .value("NUMBER_OF_TERM_SPECIFICITY", OpenMS::Ribonucleotide::TermSpecificityNuc::NUMBER_OF_TERM_SPECIFICITY)
        .export_values();

    // -----------------------------------------------------------------------
    // RibonucleotideDB
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::RibonucleotideDB>(m, "RibonucleotideDB", 
        R"doc(
Database of ribonucleotides (modified and unmodified). This is a singleton class.
The ribonucleotides are read from data/CHEMISTRY/Modomics.tsv and Custom_RNA_modifications.tsv.
)doc")

        .def_static("getInstance", []() -> OpenMS::RibonucleotideDB* { return OpenMS::RibonucleotideDB::getInstance(); }, nb::rv_policy::reference, "Returns the singleton instance")

        .def("getRibonucleotide", [](OpenMS::RibonucleotideDB& self, const OpenMS::String& code) -> const OpenMS::Ribonucleotide* {
            return self.getRibonucleotide(code);
        }, "code"_a, nb::rv_policy::reference, "Returns the ribonucleotide with the given code")

        .def("getRibonucleotidePrefix", [](OpenMS::RibonucleotideDB& self, const OpenMS::String& seq) -> const OpenMS::Ribonucleotide* {
            return self.getRibonucleotidePrefix(seq);
        }, "seq"_a, nb::rv_policy::reference, "Returns the ribonucleotide matching the longest prefix")
        .def("getRibonucleotideAlternatives", [](OpenMS::RibonucleotideDB& self, const std::string& code) {
            auto result = self.getRibonucleotideAlternatives(code);
            return nb::make_tuple(
                nb::cast(result.first, nb::rv_policy::reference),
                nb::cast(result.second, nb::rv_policy::reference)
            );
        }, "code"_a, "Returns the two alternatives for an ambiguous modification code as a tuple")
        ;

    // -----------------------------------------------------------------------
    // SequenceCoverage
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SequenceCoverage>(m, "SequenceCoverage", "Compute sequence coverage of a protein by peptide sequences")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SequenceCoverage &>())
        .def("__copy__", [](const OpenMS::SequenceCoverage& self) { return OpenMS::SequenceCoverage(self); })
        .def("__deepcopy__", [](const OpenMS::SequenceCoverage& self, nb::dict) { return OpenMS::SequenceCoverage(self); }, "memo"_a)

        .def_static("getCoverage", [](const OpenMS::AASequence& protein, const std::vector<OpenMS::AASequence>& peptides) {
            return OpenMS::SequenceCoverage::getCoverage(protein, peptides);
        }, "protein"_a, "peptides"_a, "Compute sequence coverage percentage (0-100)")
        ;

    // -----------------------------------------------------------------------
    // SequenceElement
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::SequenceElement>(m, "SequenceElement", "A single amino acid with its modifications")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::SequenceElement &>())
        .def("__copy__", [](const OpenMS::ProForma::SequenceElement& self) { return OpenMS::ProForma::SequenceElement(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::SequenceElement& self, nb::dict) { return OpenMS::ProForma::SequenceElement(self); }, "memo"_a)
        .def_rw("amino_acid", &OpenMS::ProForma::SequenceElement::amino_acid)
        .def_rw("modifications", &OpenMS::ProForma::SequenceElement::modifications)
        ;

    // -----------------------------------------------------------------------
    // Tagger
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Tagger>(m, "Tagger", 
        R"doc(
Constructor for Tagger
The parameter `max_charge_` should be >= `min_charge_`
Also `max_tag_length` should be >= `min_tag_length`
:param min_tag_length: The minimal sequence tag length
:param ppm: The tolerance for matching residue masses to peak delta masses
:param max_tag_length: The maximal sequence tag length
:param min_charge: Minimal fragment charge considered for each sequence tag
:param max_charge: Maximal fragment charge considered for each sequence tag
:param fixed_mods: A list of modification names. The modified residues replace the unmodified versions
:param var_mods: A list of modification names. The modified residues are added as additional entries to the list of residues
)doc")
        .def("__copy__", [](const OpenMS::Tagger& self) { return OpenMS::Tagger(self); })
        .def("__deepcopy__", [](const OpenMS::Tagger& self, nb::dict) { return OpenMS::Tagger(self); }, "memo"_a)
        .def(nb::init<size_t, double, size_t, size_t, size_t, std::vector<OpenMS::String>, std::vector<OpenMS::String>, bool>())
        .def("getTag", [](const OpenMS::Tagger& self, const std::vector<double>& mzs) { std::vector<std::string> tags; self.getTag(mzs, tags); return tags; }, "mzs"_a)
        .def("getTag", [](const OpenMS::Tagger& self, const OpenMS::MSSpectrum& spec) { std::vector<std::string> tags; self.getTag(spec, tags); return tags; }, "spec"_a)
        .def("setMaxCharge", [](OpenMS::Tagger& self, size_t max_charge) { return self.setMaxCharge(max_charge); }, "max_charge"_a, 
            R"doc(
Generate tags from an MSSpectrum
The parameter `tags` is filled with one string per sequence tag
It uses the standard residues from ResidueDB including
the fixed and variable modifications given to the constructor
:param spec: A centroided fragment spectrum
:param tags: The vector of tags, that is filled with this function
)doc")
        ;

    // -----------------------------------------------------------------------
    // UnlocalisedMod
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProForma::UnlocalisedMod>(m, "UnlocalisedMod", "Unlocalised modification with optional occurrence count")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProForma::UnlocalisedMod &>())
        .def("__copy__", [](const OpenMS::ProForma::UnlocalisedMod& self) { return OpenMS::ProForma::UnlocalisedMod(self); })
        .def("__deepcopy__", [](const OpenMS::ProForma::UnlocalisedMod& self, nb::dict) { return OpenMS::ProForma::UnlocalisedMod(self); }, "memo"_a)
        .def_rw("modifications", &OpenMS::ProForma::UnlocalisedMod::modifications)
        ;


    // -----------------------------------------------------------------------
    // NASequence
    // -----------------------------------------------------------------------
    auto nasequence_class = nb::class_<OpenMS::NASequence>(m, "NASequence", "OpenMS class NASequence")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::NASequence&>())
        .def("__copy__", [](const OpenMS::NASequence& self) { return OpenMS::NASequence(self); })
        .def("__deepcopy__", [](const OpenMS::NASequence& self, nb::dict) { return OpenMS::NASequence(self); }, "memo"_a)
        .def("toString", &OpenMS::NASequence::toString, "Get string representation")
        .def("__str__", [](const OpenMS::NASequence& self) { return self.toString(); })
        .def("size", &OpenMS::NASequence::size, "Get number of residues")
        .def("empty", &OpenMS::NASequence::empty, "Check if empty")
        .def("getFormula", [](const OpenMS::NASequence& self) { return self.getFormula(); }, "Get empirical formula")
        .def("getFormula", [](const OpenMS::NASequence& self, OpenMS::NASequence::NASFragmentType type, int charge) { return self.getFormula(type, charge); }, "type"_a, "charge"_a = 0, "Get empirical formula for a fragment type")
        .def("getMonoWeight", [](const OpenMS::NASequence& self) { return self.getMonoWeight(); }, "Get monoisotopic weight")
        .def("getMonoWeight", [](const OpenMS::NASequence& self, OpenMS::NASequence::NASFragmentType type, int charge) { return self.getMonoWeight(type, charge); }, "type"_a, "charge"_a = 0, "Get monoisotopic weight for a fragment type")
        .def("getAverageWeight", [](const OpenMS::NASequence& self) { return self.getAverageWeight(); }, "Get average weight")
        .def("getAverageWeight", [](const OpenMS::NASequence& self, OpenMS::NASequence::NASFragmentType type, int charge) { return self.getAverageWeight(type, charge); }, "type"_a, "charge"_a = 0, "Get average weight for a fragment type")
        .def_static("fromString", [](const std::string& s) { return OpenMS::NASequence::fromString(s); }, "s"_a, "Create NASequence from string")
        .def("__eq__", [](const OpenMS::NASequence& self, const OpenMS::NASequence& other) { return self == other; }, "other"_a)
        .def("__ne__", [](const OpenMS::NASequence& self, const OpenMS::NASequence& other) { return self != other; }, "other"_a)
        .def("__len__", [](const OpenMS::NASequence& self) { return self.size(); })
        .def("__getitem__", [](const OpenMS::NASequence& self, size_t i) -> const OpenMS::Ribonucleotide* {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference, "i"_a)
        .def("__iter__", [](const OpenMS::NASequence& self) {
            return nb::make_iterator<nb::rv_policy::reference>(nb::type<OpenMS::NASequence>(), "NASequence_iter",
                self.begin(), self.end());
        })
        .def("get", [](OpenMS::NASequence& self, size_t index) -> const OpenMS::Ribonucleotide* {
            return self.get(index);
        }, "index"_a, nb::rv_policy::reference, "Returns the ribonucleotide at the given index")
        .def("getPrefix", [](const OpenMS::NASequence& self, size_t length) { return self.getPrefix(length); }, "length"_a, "Returns the prefix of the given length")
        .def("getSuffix", [](const OpenMS::NASequence& self, size_t length) { return self.getSuffix(length); }, "length"_a, "Returns the suffix of the given length")
        .def("getSubsequence", [](const OpenMS::NASequence& self, size_t start, size_t length) { return self.getSubsequence(start, length); }, "start"_a, "length"_a, "Returns a subsequence starting at start with the given length")
        .def("set", [](OpenMS::NASequence& self, size_t index, const OpenMS::Ribonucleotide* r) { self.set(index, r); }, "index"_a, "ribonucleotide"_a, "Sets the ribonucleotide at the given index")
        .def("getFivePrimeMod", [](const OpenMS::NASequence& self) -> const OpenMS::Ribonucleotide* { return self.getFivePrimeMod(); }, nb::rv_policy::reference, "Returns the 5' modification, or None if not set")
        .def("getThreePrimeMod", [](const OpenMS::NASequence& self) -> const OpenMS::Ribonucleotide* { return self.getThreePrimeMod(); }, nb::rv_policy::reference, "Returns the 3' modification, or None if not set")
        .def("setFivePrimeMod", [](OpenMS::NASequence& self, const OpenMS::Ribonucleotide* mod) { self.setFivePrimeMod(mod); }, "mod"_a, "Sets the 5' modification")
        .def("setThreePrimeMod", [](OpenMS::NASequence& self, const OpenMS::Ribonucleotide* mod) { self.setThreePrimeMod(mod); }, "mod"_a, "Sets the 3' modification")
        .def("hasFivePrimeMod", &OpenMS::NASequence::hasFivePrimeMod, "Returns true if the sequence has a 5' modification")
        .def("hasThreePrimeMod", &OpenMS::NASequence::hasThreePrimeMod, "Returns true if the sequence has a 3' modification")
        .def("getSequence", [](const OpenMS::NASequence& self) { return self.getSequence(); }, nb::rv_policy::reference_internal, "Returns the sequence of ribonucleotides")
        .def("setSequence", [](OpenMS::NASequence& self, const std::vector<const OpenMS::Ribonucleotide*>& seq) { self.setSequence(seq); }, "seq"_a, "Sets the sequence of ribonucleotides")
        .def("__repr__", [](const OpenMS::NASequence& self) {
            std::ostringstream oss;
            oss << "NASequence(sequence='" << std::string(self.toString())
                << "', length=" << self.size()
                << ", mono_mass=" << self.getMonoWeight() << ")";
            return oss.str();
        })
        ;

    // NASFragmentType enum nested under NASequence
    nb::enum_<OpenMS::NASequence::NASFragmentType>(nasequence_class, "NASFragmentType", nb::is_arithmetic())
        .value("Full", OpenMS::NASequence::NASFragmentType::Full)
        .value("Internal", OpenMS::NASequence::NASFragmentType::Internal)
        .value("FivePrime", OpenMS::NASequence::NASFragmentType::FivePrime)
        .value("ThreePrime", OpenMS::NASequence::NASFragmentType::ThreePrime)
        .value("AIon", OpenMS::NASequence::NASFragmentType::AIon)
        .value("BIon", OpenMS::NASequence::NASFragmentType::BIon)
        .value("CIon", OpenMS::NASequence::NASFragmentType::CIon)
        .value("XIon", OpenMS::NASequence::NASFragmentType::XIon)
        .value("YIon", OpenMS::NASequence::NASFragmentType::YIon)
        .value("ZIon", OpenMS::NASequence::NASFragmentType::ZIon)
        .value("Precursor", OpenMS::NASequence::NASFragmentType::Precursor)
        .value("WIon", OpenMS::NASequence::NASFragmentType::WIon)
        .value("DIon", OpenMS::NASequence::NASFragmentType::DIon)
        ;

    // -----------------------------------------------------------------------
    // IMSAlphabet
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ims::IMSAlphabet>(m, "IMSAlphabet",
        "Indexed container of bio-chemical elements for mass decomposition")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ims::IMSAlphabet&>())
        .def("__copy__", [](const OpenMS::ims::IMSAlphabet& self) { return OpenMS::ims::IMSAlphabet(self); })
        .def("__deepcopy__", [](const OpenMS::ims::IMSAlphabet& self, nb::dict) { return OpenMS::ims::IMSAlphabet(self); }, "memo"_a)
        .def("size", [](const OpenMS::ims::IMSAlphabet& self) { return self.size(); })
        .def("getElement", [](const OpenMS::ims::IMSAlphabet& self, size_t index) -> const OpenMS::ims::IMSElement& { return self.getElement(index); }, "index"_a, nb::rv_policy::reference_internal)
        .def("getName", [](const OpenMS::ims::IMSAlphabet& self, size_t index) { return self.getName(index); }, "index"_a)
        .def("getMass", [](const OpenMS::ims::IMSAlphabet& self, size_t index) { return self.getMass(index); }, "index"_a)
        .def("hasName", [](const OpenMS::ims::IMSAlphabet& self, const std::string& name) { return self.hasName(name); }, "name"_a)
        .def("push_back", [](OpenMS::ims::IMSAlphabet& self, const std::string& name, double mass) { self.push_back(name, mass); }, "name"_a, "mass"_a)
        .def("clear", [](OpenMS::ims::IMSAlphabet& self) { self.clear(); })
        .def("sortByNames", [](OpenMS::ims::IMSAlphabet& self) { self.sortByNames(); })
        .def("sortByValues", [](OpenMS::ims::IMSAlphabet& self) { self.sortByValues(); })
        .def("__len__", [](const OpenMS::ims::IMSAlphabet& self) { return self.size(); })
        .def("erase", [](OpenMS::ims::IMSAlphabet& self, const std::string& name) { return self.erase(name); }, "name"_a, "Removes the element with the given name")
        .def("getMasses", [](const OpenMS::ims::IMSAlphabet& self, size_t isotope_index) { return self.getMasses(isotope_index); }, "isotope_index"_a = 0, "Gets masses of elements isotopes given by isotope_index")
        .def("getAverageMasses", [](const OpenMS::ims::IMSAlphabet& self) { return self.getAverageMasses(); }, "Gets average masses of elements")
        .def("load", [](OpenMS::ims::IMSAlphabet& self, const std::string& fname) { self.load(fname); }, "fname"_a, "Loads the alphabet data from the file")
        .def("setElement", [](OpenMS::ims::IMSAlphabet& self, const std::string& name, double mass, bool forced) { self.setElement(name, mass, forced); }, "name"_a, "mass"_a, "forced"_a = false, "Overwrites an element in the alphabet")
        ;

    // -----------------------------------------------------------------------
    // SimplePeak
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SimpleTSGXLMS::SimplePeak>(m, "SimplePeak",
        "Simple peak struct with m/z and charge")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::SimpleTSGXLMS::SimplePeak& self) { return OpenMS::SimpleTSGXLMS::SimplePeak(self); })
        .def("__deepcopy__", [](const OpenMS::SimpleTSGXLMS::SimplePeak& self, nb::dict) { return OpenMS::SimpleTSGXLMS::SimplePeak(self); }, "memo"_a)
        .def(nb::init<double, int>(), "mz"_a, "charge"_a)
        .def_rw("mz", &OpenMS::SimpleTSGXLMS::SimplePeak::mz)
        .def_rw("charge", &OpenMS::SimpleTSGXLMS::SimplePeak::charge)
        ;

    // -----------------------------------------------------------------------
    // __static_* module-level wrappers for ProForma
    // -----------------------------------------------------------------------
    m.def("__static_ProForma_parse", [](const OpenMS::String& input) -> OpenMS::ProForma::Peptidoform { return OpenMS::ProForma::parse(input); }, "input"_a);
    m.def("__static_ProForma_parseIon", [](const OpenMS::String& input) -> OpenMS::ProForma::PeptidoformIon { return OpenMS::ProForma::parseIon(input); }, "input"_a);
    m.def("__static_ProForma_toString", [](const OpenMS::ProForma::Peptidoform& pf, OpenMS::ProForma::WriteMode mode) -> OpenMS::String { return OpenMS::ProForma::toString(pf, mode); }, "pf"_a, "mode"_a);
    m.def("__static_ProForma_toStringIon", [](const OpenMS::ProForma::PeptidoformIon& pfi, OpenMS::ProForma::WriteMode mode) -> OpenMS::String { return OpenMS::ProForma::toString(pfi, mode); }, "pfi"_a, "mode"_a);
    m.def("__static_ProForma_peptidoformFromJSON", [](const OpenMS::String& json_str) -> OpenMS::ProForma::Peptidoform { return OpenMS::ProForma::peptidoformFromJSON(json_str); }, "json_str"_a);
    m.def("__static_ProForma_peptidoformIonFromJSON", [](const OpenMS::String& json_str) -> OpenMS::ProForma::PeptidoformIon { return OpenMS::ProForma::peptidoformIonFromJSON(json_str); }, "json_str"_a);
    m.def("__static_ProForma_resolveModifications", [](OpenMS::ProForma::Peptidoform& pf) -> void { OpenMS::ProForma::resolveModifications(pf); }, "pf"_a);
    m.def("__static_ProForma_toAASequence", [](const OpenMS::ProForma::Peptidoform& pf, OpenMS::ProForma::ConversionPolicy policy) -> OpenMS::AASequence { return OpenMS::ProForma::toAASequence(pf, policy); }, "pf"_a, "policy"_a);
    m.def("__static_ProForma_fromAASequence", [](const OpenMS::AASequence& seq) -> OpenMS::ProForma::Peptidoform { return OpenMS::ProForma::fromAASequence(seq); }, "seq"_a);
    m.def("__static_ProForma_isRepresentableAsAASequence", [](const OpenMS::ProForma::Peptidoform& pf) -> bool { return OpenMS::ProForma::isRepresentableAsAASequence(pf); }, "pf"_a);
    m.def("__static_ProForma_getAASequenceConversionIssues", [](const OpenMS::ProForma::Peptidoform& pf) -> std::vector<OpenMS::ProForma::ConversionIssue> { return OpenMS::ProForma::getAASequenceConversionIssues(pf); }, "pf"_a);
    m.def("__static_ProForma_canCalculateMass", [](const OpenMS::ProForma::Peptidoform& pf) -> bool { return OpenMS::ProForma::canCalculateMass(pf); }, "pf"_a);
    m.def("__static_ProForma_getMassCalculationIssues", [](const OpenMS::ProForma::Peptidoform& pf) -> std::vector<OpenMS::ProForma::ConversionIssue> { return OpenMS::ProForma::getMassCalculationIssues(pf); }, "pf"_a);
    m.def("__static_ProForma_getMonoWeight", [](const OpenMS::ProForma::Peptidoform& pf) -> double { return OpenMS::ProForma::getMonoWeight(pf); }, "pf"_a);
    m.def("__static_ProForma_getMZ", [](const OpenMS::ProForma::PeptidoformIon& pfi) -> double { return OpenMS::ProForma::getMZ(pfi); }, "pfi"_a);
    m.def("__static_ProForma_getMZCharge", [](const OpenMS::ProForma::Peptidoform& pf, int charge) -> double { return OpenMS::ProForma::getMZ(pf, charge); }, "pf"_a, "charge"_a);
    m.def("__static_ProForma_canGenerateSpectrum", [](const OpenMS::ProForma::Peptidoform& pf) -> bool { return OpenMS::ProForma::canGenerateSpectrum(pf); }, "pf"_a);
    m.def("__static_ProForma_getSpectrumGenerationIssues", [](const OpenMS::ProForma::Peptidoform& pf) -> std::vector<OpenMS::ProForma::ConversionIssue> { return OpenMS::ProForma::getSpectrumGenerationIssues(pf); }, "pf"_a);
    m.def("__static_ProForma_generateSpectrum", [](const OpenMS::ProForma::Peptidoform& pf, int min_charge, int max_charge, const OpenMS::String& ion_types, bool add_losses, bool add_metainfo) -> OpenMS::MSSpectrum { return OpenMS::ProForma::generateSpectrum(pf, min_charge, max_charge, ion_types, add_losses, add_metainfo); }, "pf"_a, "min_charge"_a, "max_charge"_a, "ion_types"_a, "add_losses"_a, "add_metainfo"_a);
    m.def("__static_ProForma_peptidoformToJSON", [](const OpenMS::ProForma::Peptidoform& pf) -> OpenMS::String { return OpenMS::ProForma::toJSON(pf); }, "pf"_a);
    m.def("__static_ProForma_peptidoformIonToJSON", [](const OpenMS::ProForma::PeptidoformIon& pfi) -> OpenMS::String { return OpenMS::ProForma::toJSON(pfi); }, "pfi"_a);
    m.def("__static_ProForma_toJSON", [](const OpenMS::ProForma::Peptidoform& pf) -> OpenMS::String { return OpenMS::ProForma::toJSON(pf); }, "pf"_a);

    // -----------------------------------------------------------------------
    // __static_* module-level wrappers for MzPAF
    // -----------------------------------------------------------------------
    m.def("__static_MzPAF_parse", [](const OpenMS::String& annotation) -> OpenMS::MzPAFAnnotation { return OpenMS::MzPAF::parse(annotation); }, "annotation"_a);
    m.def("__static_MzPAF_parseMultiple", [](const OpenMS::String& annotation) -> OpenMS::MzPAFPeakAnnotations { return OpenMS::MzPAF::parseMultiple(annotation); }, "annotation"_a);
    m.def("__static_MzPAF_tryParse", [](const OpenMS::String& annotation) -> std::optional<OpenMS::MzPAFAnnotation> { return OpenMS::MzPAF::tryParse(annotation); }, "annotation"_a);
    m.def("__static_MzPAF_tryParseMultiple", [](const OpenMS::String& annotation) -> std::optional<OpenMS::MzPAFPeakAnnotations> { return OpenMS::MzPAF::tryParseMultiple(annotation); }, "annotation"_a);
    m.def("__static_MzPAF_toString", [](const OpenMS::MzPAFAnnotation& annotation) -> OpenMS::String { return OpenMS::MzPAF::toString(annotation); }, "annotation"_a);
    m.def("__static_MzPAF_toStringMultiple", [](const OpenMS::MzPAFPeakAnnotations& annotations) -> OpenMS::String { return OpenMS::MzPAF::toString(annotations); }, "annotations"_a);
    m.def("__static_MzPAF_toPeakAnnotation", [](const OpenMS::MzPAFAnnotation& annotation, double mz, double intensity) -> OpenMS::PeptideHit::PeakAnnotation { return OpenMS::MzPAF::toPeakAnnotation(annotation, mz, intensity); }, "annotation"_a, "mz"_a, "intensity"_a);
    m.def("__static_MzPAF_fromPeakAnnotation", [](const OpenMS::PeptideHit::PeakAnnotation& pa) -> OpenMS::MzPAFPeakAnnotations { return OpenMS::MzPAF::fromPeakAnnotation(pa); }, "pa"_a);
    m.def("__static_MzPAF_isMzPAFFormat", [](const OpenMS::String& annotation) -> bool { return OpenMS::MzPAF::isMzPAFFormat(annotation); }, "annotation"_a);

    // -----------------------------------------------------------------------
    // AdductInfo (AMSE_AdductInfo)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AdductInfo>(m, "AMSE_AdductInfo",
        "Representation of an adduct for metabolite identification")
        .def(nb::init<const OpenMS::String&, const OpenMS::EmpiricalFormula&, int, OpenMS::UInt>(),
            "name"_a, "adduct"_a, "charge"_a, "mol_multiplier"_a = 1)
        .def("getNeutralMass", &OpenMS::AdductInfo::getNeutralMass, "observed_mz"_a,
            "Returns the neutral mass of the small molecule without adduct")
        .def("getMZ", &OpenMS::AdductInfo::getMZ, "neutral_mass"_a,
            "Returns the m/z of the small molecule with the adduct added")
        .def("getMassShift", &OpenMS::AdductInfo::getMassShift, "use_avg_mass"_a = false,
            "Returns the mass shift caused by this adduct")
        .def("isCompatible", &OpenMS::AdductInfo::isCompatible, "db_entry"_a,
            "Checks if adduct is compatible with the given formula")
        .def("getCharge", &OpenMS::AdductInfo::getCharge, "Returns the charge of the adduct")
        .def("getName", &OpenMS::AdductInfo::getName, "Returns the original name string")
        .def("getEmpiricalFormula", &OpenMS::AdductInfo::getEmpiricalFormula, "Returns the sum formula of the adduct")
        .def("getMolMultiplier", &OpenMS::AdductInfo::getMolMultiplier, "Returns the molecular multiplier")
        .def_static("parseAdductString", &OpenMS::AdductInfo::parseAdductString, "adduct"_a,
            "Parse an adduct string containing a formula and charge")
        .def("__eq__", &OpenMS::AdductInfo::operator==)
        ;
    m.def("__static_AdductInfo_parseAdductString", [](const OpenMS::String& adduct) -> OpenMS::AdductInfo { return OpenMS::AdductInfo::parseAdductString(adduct); }, "adduct"_a);

}