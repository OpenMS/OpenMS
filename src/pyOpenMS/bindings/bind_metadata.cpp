// pyOpenMS nanobind bindings
// Domain: metadata

#include "all_casters.h"
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/AbsoluteQuantitationStandards.h>
#include <OpenMS/METADATA/AnnotatedMSRun.h>
#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/Gradient.h>
#include <OpenMS/METADATA/HPLC.h>
#include <OpenMS/METADATA/IdentifierMSRunMapper.h>
#include <OpenMS/METADATA/MetaInfo.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>
#include <OpenMS/METADATA/USI.h>
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

NB_MODULE(_pyopenms_metadata, m) {
    m.doc() = "pyOpenMS metadata bindings";

    // -----------------------------------------------------------------------
    // AQS_runConcentration (AbsoluteQuantitationStandards::runConcentration)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AbsoluteQuantitationStandards::runConcentration>(m, "AQS_runConcentration",
        "Structure to hold a single run with its known concentration")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::AbsoluteQuantitationStandards::runConcentration& self) { return OpenMS::AbsoluteQuantitationStandards::runConcentration(self); })
        .def("__deepcopy__", [](const OpenMS::AbsoluteQuantitationStandards::runConcentration& self, nb::dict) { return OpenMS::AbsoluteQuantitationStandards::runConcentration(self); }, "memo"_a)
        .def_rw("sample_name", &OpenMS::AbsoluteQuantitationStandards::runConcentration::sample_name)
        .def_rw("component_name", &OpenMS::AbsoluteQuantitationStandards::runConcentration::component_name)
        .def_rw("IS_component_name", &OpenMS::AbsoluteQuantitationStandards::runConcentration::IS_component_name)
        .def_rw("actual_concentration", &OpenMS::AbsoluteQuantitationStandards::runConcentration::actual_concentration)
        .def_rw("IS_actual_concentration", &OpenMS::AbsoluteQuantitationStandards::runConcentration::IS_actual_concentration)
        .def_rw("concentration_units", &OpenMS::AbsoluteQuantitationStandards::runConcentration::concentration_units)
        .def_rw("dilution_factor", &OpenMS::AbsoluteQuantitationStandards::runConcentration::dilution_factor)
        ;

    // -----------------------------------------------------------------------
    // AQS_featureConcentration (AbsoluteQuantitationStandards::featureConcentration)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AbsoluteQuantitationStandards::featureConcentration>(m, "AQS_featureConcentration",
        "Structure to hold a single component with its corresponding known concentration")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::AbsoluteQuantitationStandards::featureConcentration& self) { return OpenMS::AbsoluteQuantitationStandards::featureConcentration(self); })
        .def("__deepcopy__", [](const OpenMS::AbsoluteQuantitationStandards::featureConcentration& self, nb::dict) { return OpenMS::AbsoluteQuantitationStandards::featureConcentration(self); }, "memo"_a)
        .def_rw("feature", &OpenMS::AbsoluteQuantitationStandards::featureConcentration::feature)
        .def_rw("IS_feature", &OpenMS::AbsoluteQuantitationStandards::featureConcentration::IS_feature)
        .def_rw("actual_concentration", &OpenMS::AbsoluteQuantitationStandards::featureConcentration::actual_concentration)
        .def_rw("IS_actual_concentration", &OpenMS::AbsoluteQuantitationStandards::featureConcentration::IS_actual_concentration)
        .def_rw("concentration_units", &OpenMS::AbsoluteQuantitationStandards::featureConcentration::concentration_units)
        .def_rw("dilution_factor", &OpenMS::AbsoluteQuantitationStandards::featureConcentration::dilution_factor)
        ;

    // -----------------------------------------------------------------------
    // AbsoluteQuantitationStandards
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AbsoluteQuantitationStandards>(m, "AbsoluteQuantitationStandards", 
        R"doc(
AbsoluteQuantitationStandards is a class to handle the relationship
between runs, components, and their actual concentrations
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::AbsoluteQuantitationStandards& self) { return OpenMS::AbsoluteQuantitationStandards(self); })
        .def("__deepcopy__", [](const OpenMS::AbsoluteQuantitationStandards& self, nb::dict) { return OpenMS::AbsoluteQuantitationStandards(self); }, "memo"_a)
        .def("getComponentFeatureConcentrations", [](const OpenMS::AbsoluteQuantitationStandards& self,
            const std::vector<OpenMS::AbsoluteQuantitationStandards::runConcentration>& run_concentrations,
            const std::vector<OpenMS::FeatureMap>& feature_maps,
            const OpenMS::String& component_name) {
            std::vector<OpenMS::AbsoluteQuantitationStandards::featureConcentration> feature_concentrations;
            self.getComponentFeatureConcentrations(run_concentrations, feature_maps, component_name, feature_concentrations);
            return feature_concentrations;
        }, "run_concentrations"_a, "feature_maps"_a, "component_name"_a, "Gets the feature concentrations from run concentrations and feature maps")
        ;

    // -----------------------------------------------------------------------
    // AnnotatedMSRun
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AnnotatedMSRun>(m, "AnnotatedMSRun", 
        R"doc(
Class for storing MS run data with peptide and protein identifications
This class stores an MSExperiment (containing spectra) along with peptide and protein
identifications. Each spectrum in the MSExperiment is associated with a single
PeptideIdentification object. Object gets typically not manually created but generated
by the IDMapper class.
Usage:
.. code-block:: python
run = AnnotatedMSRun()
exp = MSExperiment()
MzMLFile().load(path_to_file, exp)
run.setMSExperiment(exp)
run.setPeptideIdentifications(my_peptide_ids)
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::AnnotatedMSRun &>())
        .def("__copy__", [](const OpenMS::AnnotatedMSRun& self) { return OpenMS::AnnotatedMSRun(self); })
        .def("__deepcopy__", [](const OpenMS::AnnotatedMSRun& self, nb::dict) { return OpenMS::AnnotatedMSRun(self); }, "memo"_a)
        .def("getProteinIdentifications", [](OpenMS::AnnotatedMSRun& self) -> std::vector<OpenMS::ProteinIdentification> & { return self.getProteinIdentifications(); }, nb::rv_policy::reference_internal)
        .def("setProteinIdentifications", [](OpenMS::AnnotatedMSRun& self, const std::vector<OpenMS::ProteinIdentification>& ids) { return self.setProteinIdentifications(ids); }, "ids"_a)
        .def("setProteinIdentifications", [](OpenMS::AnnotatedMSRun& self, std::vector<OpenMS::ProteinIdentification>& ids) { return self.setProteinIdentifications(ids); }, "ids"_a)
        .def("getPeptideIdentifications", [](OpenMS::AnnotatedMSRun& self) -> OpenMS::PeptideIdentificationList & { return self.getPeptideIdentifications(); }, nb::rv_policy::reference_internal)
        .def("setPeptideIdentifications", [](OpenMS::AnnotatedMSRun& self, OpenMS::PeptideIdentificationList& ids) { return self.setPeptideIdentifications(ids); }, "ids"_a)
        .def("setPeptideIdentifications", [](OpenMS::AnnotatedMSRun& self, const OpenMS::PeptideIdentificationList& ids) { return self.setPeptideIdentifications(ids); }, "ids"_a)
        .def("getMSExperiment", [](OpenMS::AnnotatedMSRun& self) -> OpenMS::MSExperiment & { return self.getMSExperiment(); }, nb::rv_policy::reference_internal)
        .def("setMSExperiment", [](OpenMS::AnnotatedMSRun& self, OpenMS::MSExperiment& experiment) { return self.setMSExperiment(experiment); }, "experiment"_a)
        .def("setMSExperiment", [](OpenMS::AnnotatedMSRun& self, const OpenMS::MSExperiment& experiment) { return self.setMSExperiment(experiment); }, "experiment"_a)
        .def("__hash__", [](const OpenMS::AnnotatedMSRun& self) { return std::hash<OpenMS::AnnotatedMSRun>{}(self); })
        ;

    // -----------------------------------------------------------------------
    // Unit (CVTerm::Unit)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CVTerm::Unit>(m, "Unit",
        "Unit for a controlled vocabulary term")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::String&, const OpenMS::String&, const OpenMS::String&>(), "accession"_a, "name"_a, "cv_ref"_a)
        .def(nb::init<const OpenMS::CVTerm::Unit &>())
        .def("__copy__", [](const OpenMS::CVTerm::Unit& self) { return OpenMS::CVTerm::Unit(self); })
        .def("__deepcopy__", [](const OpenMS::CVTerm::Unit& self, nb::dict) { return OpenMS::CVTerm::Unit(self); }, "memo"_a)
        .def_rw("accession", &OpenMS::CVTerm::Unit::accession)
        .def_rw("name", &OpenMS::CVTerm::Unit::name)
        .def_rw("cv_ref", &OpenMS::CVTerm::Unit::cv_ref)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        ;

    // -----------------------------------------------------------------------
    // CVTerm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CVTerm>(m, "CVTerm", "Representation of controlled vocabulary term")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::CVTerm &>())
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::String, OpenMS::String, OpenMS::CVTerm::Unit>())
        .def("__copy__", [](const OpenMS::CVTerm& self) { return OpenMS::CVTerm(self); })
        .def("__deepcopy__", [](const OpenMS::CVTerm& self, nb::dict) { return OpenMS::CVTerm(self); }, "memo"_a)
        .def("setAccession", [](OpenMS::CVTerm& self, const OpenMS::String& accession) { return self.setAccession(accession); }, "accession"_a, "Sets the accession string of the term")
        .def("getAccession", [](const OpenMS::CVTerm& self) { return self.getAccession(); }, "Returns the accession string of the term")
        .def("setName", [](OpenMS::CVTerm& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the term")
        .def("getName", [](const OpenMS::CVTerm& self) { return self.getName(); }, "Returns the name of the term")
        .def("setCVIdentifierRef", [](OpenMS::CVTerm& self, const OpenMS::String& cv_identifier_ref) { return self.setCVIdentifierRef(cv_identifier_ref); }, "cv_identifier_ref"_a, "Sets the CV identifier reference string, e.g. UO for unit obo")
        .def("getCVIdentifierRef", [](const OpenMS::CVTerm& self) { return self.getCVIdentifierRef(); }, "Returns the CV identifier reference string")
        .def("setValue", [](OpenMS::CVTerm& self, const OpenMS::DataValue& value) { return self.setValue(value); }, "value"_a, "Sets the value of the term")
        .def("getValue", [](const OpenMS::CVTerm& self) { return self.getValue(); }, "Returns the value of the term")
        .def("setUnit", [](OpenMS::CVTerm& self, const OpenMS::CVTerm::Unit& unit) { return self.setUnit(unit); }, "unit"_a, "Sets the unit of the term")
        .def("getUnit", [](const OpenMS::CVTerm& self) -> const OpenMS::CVTerm::Unit & { return self.getUnit(); }, nb::rv_policy::reference_internal, "Returns the unit")
        .def(nb::self == nb::self)
        .def("hasValue", [](const OpenMS::CVTerm& self) { return self.hasValue(); }, "Checks whether the term has a value")
        .def("hasUnit", [](const OpenMS::CVTerm& self) { return self.hasUnit(); }, "Checks whether the term has a unit")
        .def("__hash__", [](const OpenMS::CVTerm& self) { return std::hash<OpenMS::CVTerm>{}(self); })
        ;

    // -----------------------------------------------------------------------
    // ExperimentalDesign
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ExperimentalDesign>(m, "ExperimentalDesign", "Representation of an experimental design in OpenMS. Instances can be loaded with the ExperimentalDesignFile class")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::ExperimentalDesign& self) { return OpenMS::ExperimentalDesign(self); })
        .def("__deepcopy__", [](const OpenMS::ExperimentalDesign& self, nb::dict) { return OpenMS::ExperimentalDesign(self); }, "memo"_a)
        .def(nb::init<std::vector<OpenMS::ExperimentalDesign::MSFileSectionEntry>, OpenMS::ExperimentalDesign::SampleSection>())
        .def("getNumberOfSamples", [](const OpenMS::ExperimentalDesign& self) { return self.getNumberOfSamples(); }, "Returns the number of samples measured (= highest sample index)")
        .def("getNumberOfFractions", [](const OpenMS::ExperimentalDesign& self) { return self.getNumberOfFractions(); }, "Returns the number of fractions (= highest fraction index)")
        .def("getNumberOfLabels", [](const OpenMS::ExperimentalDesign& self) { return self.getNumberOfLabels(); }, "Returns the number of labels per file")
        .def("getNumberOfMSFiles", [](const OpenMS::ExperimentalDesign& self) { return self.getNumberOfMSFiles(); }, "Returns the number of MS files (= fractions * fraction_groups)")
        .def("getNumberOfFractionGroups", [](const OpenMS::ExperimentalDesign& self) { return self.getNumberOfFractionGroups(); }, "Allows to group fraction ids and source files. Return the number of fraction_groups")
        .def("getSample", [](OpenMS::ExperimentalDesign& self, unsigned int fraction_group, unsigned int label) { return self.getSample(fraction_group, label); }, "fraction_group"_a, "label"_a = 1, "Returns sample index (depends on fraction_group and label)")
        .def("isFractionated", [](const OpenMS::ExperimentalDesign& self) { return self.isFractionated(); }, "Returns whether at least one fraction_group in this experimental design is fractionated")
        .def("sameNrOfMSFilesPerFraction", [](const OpenMS::ExperimentalDesign& self) { return self.sameNrOfMSFilesPerFraction(); }, "Returns if each fraction number is associated with the same number of fraction_group")

        .def_static("fromConsensusMap", [](const OpenMS::ConsensusMap& c) {
            return OpenMS::ExperimentalDesign::fromConsensusMap(c);
        }, "c"_a, "Extract experimental design from consensus map")

        .def_static("fromFeatureMap", [](const OpenMS::FeatureMap& f) {
            return OpenMS::ExperimentalDesign::fromFeatureMap(f);
        }, "f"_a, "Extract experimental design from feature map")

        .def_static("fromIdentifications", [](const std::vector<OpenMS::ProteinIdentification>& proteins) {
            return OpenMS::ExperimentalDesign::fromIdentifications(proteins);
        }, "proteins"_a, "Extract experimental design from identifications")
        .def("getMSFileSection", [](const OpenMS::ExperimentalDesign& self) { return self.getMSFileSection(); }, "Returns the MS file section")
        .def("setMSFileSection", [](OpenMS::ExperimentalDesign& self, const OpenMS::ExperimentalDesign::MSFileSection& msfile_section) { self.setMSFileSection(msfile_section); }, "msfile_section"_a, "Sets the MS file section")
        .def("getSampleSection", [](const OpenMS::ExperimentalDesign& self) { return self.getSampleSection(); }, "Returns the sample section")
        .def("setSampleSection", [](OpenMS::ExperimentalDesign& self, const OpenMS::ExperimentalDesign::SampleSection& sample_section) { self.setSampleSection(sample_section); }, "sample_section"_a, "Sets the sample section")
        ;

    // -----------------------------------------------------------------------
    // ExperimentalDesign_MSFileSectionEntry
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ExperimentalDesign::MSFileSectionEntry>(m, "ExperimentalDesign_MSFileSectionEntry", "OpenMS class ExperimentalDesign_MSFileSectionEntry")
        .def(nb::init<>())
        .def_rw("fraction_group", &OpenMS::ExperimentalDesign::MSFileSectionEntry::fraction_group)
        .def_rw("fraction", &OpenMS::ExperimentalDesign::MSFileSectionEntry::fraction)
        .def_rw("path", &OpenMS::ExperimentalDesign::MSFileSectionEntry::path)
        .def_rw("label", &OpenMS::ExperimentalDesign::MSFileSectionEntry::label)
        .def_rw("sample", &OpenMS::ExperimentalDesign::MSFileSectionEntry::sample)
        .def_rw("sample_name", &OpenMS::ExperimentalDesign::MSFileSectionEntry::sample_name)
        ;

    // -----------------------------------------------------------------------
    // ExperimentalDesign_SampleSection
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ExperimentalDesign::SampleSection>(m, "ExperimentalDesign_SampleSection", "OpenMS class ExperimentalDesign_SampleSection")
        .def(nb::init<>())
        .def(nb::init<std::vector<std::vector<OpenMS::String>>, std::map<OpenMS::String, size_t>, std::map<OpenMS::String, size_t>>())
        .def("getSamples", [](const OpenMS::ExperimentalDesign::SampleSection& self) { return self.getSamples(); }, "Returns a set of all samples that are present in the sample section")
        .def("getFactors", [](const OpenMS::ExperimentalDesign::SampleSection& self) { return self.getFactors(); }, "Returns a set of all factors (column names) that were defined for the sample section")
        .def("hasSample", [](const OpenMS::ExperimentalDesign::SampleSection& self, const OpenMS::String& sample) { return self.hasSample(sample); }, "sample"_a, "Checks whether sample section has row for a sample number")
        .def("hasFactor", [](const OpenMS::ExperimentalDesign::SampleSection& self, const OpenMS::String& factor) { return self.hasFactor(factor); }, "factor"_a, "Checks whether Sample Section has a specific factor (i.e. column name)")
        .def("getFactorValue", [](const OpenMS::ExperimentalDesign::SampleSection& self, const OpenMS::String& sample_name, const OpenMS::String& factor) { return self.getFactorValue(sample_name, factor); }, "sample_name"_a, "factor"_a, "Returns value of factor for given sample and factor name")
        .def("getFactorValue", [](const OpenMS::ExperimentalDesign::SampleSection& self, unsigned int sample_idx, const OpenMS::String& factor) { return self.getFactorValue(sample_idx, factor); }, "sample_idx"_a, "factor"_a, "Returns value of factor for given sample and factor name")
        ;

    // -----------------------------------------------------------------------
    // Gradient
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Gradient>(m, "Gradient", "Representation of a HPLC gradient")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Gradient &>())
        .def("__copy__", [](const OpenMS::Gradient& self) { return OpenMS::Gradient(self); })
        .def("__deepcopy__", [](const OpenMS::Gradient& self, nb::dict) { return OpenMS::Gradient(self); }, "memo"_a)
        .def("addEluent", [](OpenMS::Gradient& self, const OpenMS::String& eluent) { return self.addEluent(eluent); }, "eluent"_a, "Adds an eluent at the end of the eluent array")
        .def("clearEluents", [](OpenMS::Gradient& self) { return self.clearEluents(); }, "Removes all eluents")
        .def("getEluents", [](const OpenMS::Gradient& self) -> const std::vector<OpenMS::String> & { return self.getEluents(); }, nb::rv_policy::reference_internal, "Returns a reference to the list of eluents")
        .def("addTimepoint", [](OpenMS::Gradient& self, int timepoint) { return self.addTimepoint(timepoint); }, "timepoint"_a, "Adds a timepoint at the end of the timepoint array")
        .def("clearTimepoints", [](OpenMS::Gradient& self) { return self.clearTimepoints(); }, "Removes all timepoints")
        .def("getTimepoints", [](const OpenMS::Gradient& self) -> const std::vector<int> & { return self.getTimepoints(); }, nb::rv_policy::reference_internal, "Returns a reference to the list of timepoints")
        .def("setPercentage", [](OpenMS::Gradient& self, const OpenMS::String& eluent, int timepoint, unsigned int percentage) { return self.setPercentage(eluent, timepoint, percentage); }, "eluent"_a, "timepoint"_a, "percentage"_a, "Sets the percentage of 'eluent' at 'timepoint'")
        .def("getPercentage", [](const OpenMS::Gradient& self, const OpenMS::String& eluent, int timepoint) { return self.getPercentage(eluent, timepoint); }, "eluent"_a, "timepoint"_a, "Returns a const reference to the percentages")
        .def("clearPercentages", [](OpenMS::Gradient& self) { return self.clearPercentages(); }, "Sets all percentage values to 0")
        .def("isValid", [](const OpenMS::Gradient& self) { return self.isValid(); }, "Checks if the percentages of all timepoints add up to 100%")
        ;

    // -----------------------------------------------------------------------
    // HPLC
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::HPLC>(m, "HPLC", "Representation of a HPLC experiment")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::HPLC &>())
        .def("__copy__", [](const OpenMS::HPLC& self) { return OpenMS::HPLC(self); })
        .def("__deepcopy__", [](const OpenMS::HPLC& self, nb::dict) { return OpenMS::HPLC(self); }, "memo"_a)
        .def("getInstrument", [](const OpenMS::HPLC& self) { return self.getInstrument(); }, "Returns a reference to the instument name")
        .def("setInstrument", [](OpenMS::HPLC& self, const OpenMS::String& instrument) { return self.setInstrument(instrument); }, "instrument"_a, "Sets the instument name")
        .def("getColumn", [](const OpenMS::HPLC& self) { return self.getColumn(); }, "Returns a reference to the column description")
        .def("setColumn", [](OpenMS::HPLC& self, const OpenMS::String& column) { return self.setColumn(column); }, "column"_a, "Sets the column description")
        .def("getTemperature", [](const OpenMS::HPLC& self) { return self.getTemperature(); }, "Returns the temperature (in degree C)")
        .def("setTemperature", [](OpenMS::HPLC& self, int temperature) { return self.setTemperature(temperature); }, "temperature"_a, "Sets the temperature (in degree C)")
        .def("getPressure", [](const OpenMS::HPLC& self) { return self.getPressure(); }, "Returns the pressure (in bar)")
        .def("setPressure", [](OpenMS::HPLC& self, unsigned int pressure) { return self.setPressure(pressure); }, "pressure"_a, "Sets the pressure (in bar)")
        .def("getFlux", [](const OpenMS::HPLC& self) { return self.getFlux(); }, "Returns the flux (in microliter/sec)")
        .def("setFlux", [](OpenMS::HPLC& self, unsigned int flux) { return self.setFlux(flux); }, "flux"_a, "Sets the flux (in microliter/sec)")
        .def("getComment", [](const OpenMS::HPLC& self) { return self.getComment(); }, "Returns the comments")
        .def("setComment", [](OpenMS::HPLC& self, OpenMS::String comment) { return self.setComment(comment); }, "comment"_a, "Sets the comments")
        .def("getGradient", [](OpenMS::HPLC& self) -> OpenMS::Gradient & { return self.getGradient(); }, nb::rv_policy::reference_internal, "Returns a mutable reference to the used gradient")
        .def("setGradient", [](OpenMS::HPLC& self, const OpenMS::Gradient& gradient) { return self.setGradient(gradient); }, "gradient"_a, "Sets the used gradient")
        ;

    // -----------------------------------------------------------------------
    // IdentifierMSRunMapper
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IdentifierMSRunMapper>(m, "IdentifierMSRunMapper", 
        R"doc(
Maps protein identification identifiers to MS run file paths.
This class encapsulates the mapping between ProteinIdentification identifiers
and their associated MS run paths. It is useful for resolving the correct
source file for peptide identifications, especially in merged identification results.
Example usage:
.. code-block:: python
# Create mapping from protein identifications
mapper = oms.IdentifierMSRunMapper(protein_ids)
# Get MS run paths for a peptide's identifier
paths = mapper.getMSRunPaths(pep_id.getIdentifier())
# Build a USI using the mapping
usi = pep_id.buildUSI(mapper, "PXD000561", False)
)doc")
        .def(nb::init<>())
        .def(nb::init<std::vector<OpenMS::ProteinIdentification>>())
        .def("create", [](OpenMS::IdentifierMSRunMapper& self, const std::vector<OpenMS::ProteinIdentification>& prot_ids) { return self.create(prot_ids); }, "prot_ids"_a, 
            R"doc(
Construct mapper from a list of ProteinIdentifications.
:param prot_ids: List of ProteinIdentification objects
)doc")
        .def("hasIdentifier", [](const OpenMS::IdentifierMSRunMapper& self, const OpenMS::String& identifier) { return self.hasIdentifier(identifier); }, "identifier"_a,
            R"doc(
Check if the mapping contains an entry for the given identifier.
:param identifier: ProteinIdentification identifier
:return: True if identifier exists in mapping
)doc")
        .def("empty", [](const OpenMS::IdentifierMSRunMapper& self) { return self.empty(); },
            R"doc(
Check if the mapping is empty.
:return: True if no mappings exist
)doc")
        .def("size", [](const OpenMS::IdentifierMSRunMapper& self) { return self.size(); },
            R"doc(
Get the number of identifier mappings.
:return: Number of identifiers in the mapping
)doc")
        .def("getMSRunPaths", [](const OpenMS::IdentifierMSRunMapper& self, const OpenMS::String& identifier) -> const std::vector<OpenMS::String> & { return self.getMSRunPaths(identifier); }, "identifier"_a, nb::rv_policy::reference_internal,
            R"doc(
Get the MS run paths associated with the given identifier.
:param identifier: ProteinIdentification identifier
:return: List of MS run file paths associated with this identifier
)doc")
        .def("__len__", [](OpenMS::IdentifierMSRunMapper& self) { return self.size(); })
        ;

    // -----------------------------------------------------------------------
    // MetaInfo
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MetaInfo>(m, "MetaInfo", 
        R"doc(
A Type-Name-Value tuple class
MetaInfo maps an index (an integer corresponding to a string) to
DataValue objects.  The mapping of strings to the index is performed by
the MetaInfoRegistry, which can be accessed by the method registry()
There are two versions of nearly all members. One which operates with a
string name and another one which operates on an index. The index version
is always faster, as it does not need to look up the index corresponding
to the string in the MetaInfoRegistry
If you wish to add a MetaInfo member to a class, consider deriving that
class from MetaInfoInterface, instead of simply adding MetaInfo as
member. MetaInfoInterface implements a full interface to a MetaInfo
member and is more memory efficient if no meta info gets added
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MetaInfo &>())
        .def("__copy__", [](const OpenMS::MetaInfo& self) { return OpenMS::MetaInfo(self); })
        .def("__deepcopy__", [](const OpenMS::MetaInfo& self, nb::dict) { return OpenMS::MetaInfo(self); }, "memo"_a)
        .def("getValue", [](const OpenMS::MetaInfo& self, const OpenMS::String& name, const OpenMS::DataValue& default_value) { return self.getValue(name, default_value); }, "name"_a, "default_value"_a, "Returns the value corresponding to a string")
        .def("getValue", [](const OpenMS::MetaInfo& self, unsigned int index, const OpenMS::DataValue& default_value) { return self.getValue(index, default_value); }, "index"_a, "default_value"_a, "Returns the value corresponding to a string")
        .def("exists", [](const OpenMS::MetaInfo& self, const OpenMS::String& name) { return self.exists(name); }, "name"_a, "Returns if this MetaInfo is set")
        .def("exists", [](const OpenMS::MetaInfo& self, unsigned int index) { return self.exists(index); }, "index"_a, "Returns if this MetaInfo is set")
        .def("setValue", [](OpenMS::MetaInfo& self, const OpenMS::String& name, const OpenMS::DataValue& value) { return self.setValue(name, value); }, "name"_a, "value"_a, "Sets the DataValue corresponding to a name")
        .def("setValue", [](OpenMS::MetaInfo& self, unsigned int index, const OpenMS::DataValue& value) { return self.setValue(index, value); }, "index"_a, "value"_a, "Sets the DataValue corresponding to a name")
        .def("removeValue", [](OpenMS::MetaInfo& self, const OpenMS::String& name) { return self.removeValue(name); }, "name"_a, "Removes the DataValue corresponding to `name` if it exists")
        .def("removeValue", [](OpenMS::MetaInfo& self, unsigned int index) { return self.removeValue(index); }, "index"_a, "Removes the DataValue corresponding to `name` if it exists")
        .def_static("registry", []() { return OpenMS::MetaInfo::registry(); })
        
        .def("getKeys", [](const OpenMS::MetaInfo& self, nb::list py_keys) {
            std::vector<OpenMS::String> keys;
            self.getKeys(keys);
            py_keys.attr("clear")();
            for (const auto& k : keys) {
                py_keys.append(nb::str(k.c_str()));
            }
        }, "keys"_a, "Fills the given list with all meta value keys")
        .def("getKeys", [](const OpenMS::MetaInfo& self) {
            std::vector<OpenMS::String> keys;
            self.getKeys(keys);
            nb::list result;
            for (const auto& k : keys) {
                result.append(nb::str(k.c_str()));
            }
            return result;
        }, "Returns all meta value keys as a list")
        
        .def("getKeys", [](const OpenMS::MetaInfo& self, nb::list py_keys) {
            std::vector<OpenMS::String> keys;
            self.getKeys(keys);
            py_keys.attr("clear")();
            for (const auto& k : keys) {
                py_keys.append(nb::str(k.c_str()));
            }
        }, "keys"_a, "Fills the given list with all meta value keys")
        .def("getKeys", [](const OpenMS::MetaInfo& self) {
            std::vector<OpenMS::String> keys;
            self.getKeys(keys);
            nb::list result;
            for (const auto& k : keys) {
                result.append(nb::str(k.c_str()));
            }
            return result;
        }, "Returns all meta value keys as a list")
        .def("empty", [](const OpenMS::MetaInfo& self) { return self.empty(); }, "Returns if the MetaInfo is empty")
        .def("clear", [](OpenMS::MetaInfo& self) { return self.clear(); }, "Removes all meta values")
        .def("__hash__", [](const OpenMS::MetaInfo& self) { return std::hash<OpenMS::MetaInfo>{}(self); })
        .def("getKeysAsIntegers", [](const OpenMS::MetaInfo& self) {
            std::vector<OpenMS::UInt> keys;
            self.getKeys(keys);
            return keys;
        }, "Returns a list of all integer keys for which a value is set")
        ;

    // -----------------------------------------------------------------------
    // MetaInfoRegistry
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MetaInfoRegistry>(m, "MetaInfoRegistry", 
        R"doc(
Registry which assigns unique integer indices to strings
When registering a new name an index >= 1024 is assigned.
Indices from 1 to 1023 are reserved for fast access and will never change:
1 - isotopic_range
2 - cluster_id
3 - label
4 - icon
5 - color
6 - RT
7 - MZ
8 - predicted_RT
9 - predicted_RT_p_value
10 - spectrum_reference
11 - ID
12 - low_quality
13 - charge
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MetaInfoRegistry &>())
        .def("__copy__", [](const OpenMS::MetaInfoRegistry& self) { return OpenMS::MetaInfoRegistry(self); })
        .def("__deepcopy__", [](const OpenMS::MetaInfoRegistry& self, nb::dict) { return OpenMS::MetaInfoRegistry(self); }, "memo"_a)
        .def("registerName", [](OpenMS::MetaInfoRegistry& self, const OpenMS::String& name, const OpenMS::String& description, const OpenMS::String& unit) { return self.registerName(name, description, unit); }, "name"_a, "description"_a = "", "unit"_a = "", "Registers a string, stores its description and unit, and returns the corresponding index. If the string is already registered, it returns the index of the string")
        .def("setDescription", [](OpenMS::MetaInfoRegistry& self, unsigned int index, const OpenMS::String& description) { return self.setDescription(index, description); }, "index"_a, "description"_a, "Sets the description (String), corresponding to an index")
        .def("setDescription", [](OpenMS::MetaInfoRegistry& self, const OpenMS::String& name, const OpenMS::String& description) { return self.setDescription(name, description); }, "name"_a, "description"_a, "Sets the description (String), corresponding to an index")
        .def("setUnit", [](OpenMS::MetaInfoRegistry& self, unsigned int index, const OpenMS::String& unit) { return self.setUnit(index, unit); }, "index"_a, "unit"_a, "Sets the unit (String), corresponding to an index")
        .def("setUnit", [](OpenMS::MetaInfoRegistry& self, const OpenMS::String& name, const OpenMS::String& unit) { return self.setUnit(name, unit); }, "name"_a, "unit"_a, "Sets the unit (String), corresponding to an index")
        .def("getIndex", [](const OpenMS::MetaInfoRegistry& self, const OpenMS::String& name) { return self.getIndex(name); }, "name"_a, "Returns the integer index corresponding to a string. If the string is not registered, returns UInt(-1) (= UINT_MAX)")
        .def("getName", [](const OpenMS::MetaInfoRegistry& self, unsigned int index) { return self.getName(index); }, "index"_a, "Returns the corresponding name to an index")
        .def("getDescription", [](const OpenMS::MetaInfoRegistry& self, unsigned int index) { return self.getDescription(index); }, "index"_a, "Returns the description of an index")
        .def("getDescription", [](const OpenMS::MetaInfoRegistry& self, const OpenMS::String& name) { return self.getDescription(name); }, "name"_a, "Returns the description of an index")
        .def("getUnit", [](const OpenMS::MetaInfoRegistry& self, unsigned int index) { return self.getUnit(index); }, "index"_a, "Returns the unit of an index")
        .def("getUnit", [](const OpenMS::MetaInfoRegistry& self, const OpenMS::String& name) { return self.getUnit(name); }, "name"_a, "Returns the unit of an index")
        ;

    // -----------------------------------------------------------------------
    // PeptideEvidence
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeptideEvidence>(m, "PeptideEvidence", "Representation of a peptide evidence")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PeptideEvidence &>())
        .def(nb::init<OpenMS::String, int, int, char, char>())
        .def("__copy__", [](const OpenMS::PeptideEvidence& self) { return OpenMS::PeptideEvidence(self); })
        .def("__deepcopy__", [](const OpenMS::PeptideEvidence& self, nb::dict) { return OpenMS::PeptideEvidence(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("hasValidLimits", [](const OpenMS::PeptideEvidence& self) { return self.hasValidLimits(); }, "Start and end numbers in evidence represent actual numeric indices")
        .def("getProteinAccession", [](const OpenMS::PeptideEvidence& self) { return self.getProteinAccession(); }, "Returns the protein accession the peptide matches to. If not available the empty string is returned")
        .def("setProteinAccession", [](OpenMS::PeptideEvidence& self, const OpenMS::String& s) { return self.setProteinAccession(s); }, "s"_a, "Sets the protein accession the peptide matches to. If not available set to empty string")
        .def("setStart", [](OpenMS::PeptideEvidence& self, int a) { return self.setStart(a); }, "a"_a, "Sets the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set to UNKNOWN_POSITION. N-terminal positions must be marked with `N_TERMINAL_AA`")
        .def("getStart", [](const OpenMS::PeptideEvidence& self) { return self.getStart(); }, "Returns the position in the protein (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned")
        .def("setEnd", [](OpenMS::PeptideEvidence& self, int a) { return self.setEnd(a); }, "a"_a, "Sets the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available, set UNKNOWN_POSITION. C-terminal positions must be marked with C_TERMINAL_AA")
        .def("getEnd", [](const OpenMS::PeptideEvidence& self) { return self.getEnd(); }, "Returns the position of the last AA of the peptide in protein coordinates (starting at 0 for the N-terminus). If not available UNKNOWN_POSITION constant is returned")
        .def("setAABefore", [](OpenMS::PeptideEvidence& self, char acid) { return self.setAABefore(acid); }, "acid"_a, "Sets the amino acid single letter code before the sequence (preceding amino acid in the protein). If not available, set to UNKNOWN_AA. If N-terminal set to N_TERMINAL_AA")
        .def("getAABefore", [](const OpenMS::PeptideEvidence& self) { return self.getAABefore(); }, "Returns the amino acid single letter code before the sequence (preceding amino acid in the protein). If not available, UNKNOWN_AA is returned. If N-terminal, N_TERMINAL_AA is returned")
        .def("setAAAfter", [](OpenMS::PeptideEvidence& self, char acid) { return self.setAAAfter(acid); }, "acid"_a, "Sets the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, set to UNKNOWN_AA. If C-terminal set to C_TERMINAL_AA")
        .def("getAAAfter", [](const OpenMS::PeptideEvidence& self) { return self.getAAAfter(); }, "Returns the amino acid single letter code after the sequence (subsequent amino acid in the protein). If not available, UNKNOWN_AA is returned. If C-terminal, C_TERMINAL_AA is returned")
        .def("__hash__", [](const OpenMS::PeptideEvidence& self) { return std::hash<OpenMS::PeptideEvidence>{}(self); })
        .def("__repr__", [](const OpenMS::PeptideEvidence& self) {
            std::ostringstream oss;
            oss << "PeptideEvidence(protein='" << std::string(self.getProteinAccession())
                << "', start=" << self.getStart()
                << ", end=" << self.getEnd()
                << ", aa_before='" << self.getAABefore()
                << "', aa_after='" << self.getAAAfter() << "')";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::PeptideEvidence& self) { return nb::cast(self).attr("__repr__")(); })
        ;

    // -----------------------------------------------------------------------
    // PeptideHit_AnalysisResult (PeptideHit::PepXMLAnalysisResult)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeptideHit::PepXMLAnalysisResult>(m, "PeptideHit_AnalysisResult",
        "Analysis result from pepXML post-processing tools (e.g. PeptideProphet, iProphet)")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::PeptideHit::PepXMLAnalysisResult& self) { return OpenMS::PeptideHit::PepXMLAnalysisResult(self); })
        .def("__deepcopy__", [](const OpenMS::PeptideHit::PepXMLAnalysisResult& self, nb::dict) { return OpenMS::PeptideHit::PepXMLAnalysisResult(self); }, "memo"_a)
        .def_rw("score_type", &OpenMS::PeptideHit::PepXMLAnalysisResult::score_type)
        .def_rw("higher_is_better", &OpenMS::PeptideHit::PepXMLAnalysisResult::higher_is_better)
        .def_rw("main_score", &OpenMS::PeptideHit::PepXMLAnalysisResult::main_score)
        .def_rw("sub_scores", &OpenMS::PeptideHit::PepXMLAnalysisResult::sub_scores)
        .def(nb::self == nb::self)
        ;

    // -----------------------------------------------------------------------
    // PeptideHit_PeakAnnotation
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeptideHit::PeakAnnotation>(m, "PeptideHit_PeakAnnotation", "OpenMS class PeptideHit_PeakAnnotation")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PeptideHit::PeakAnnotation &>())
        .def("__copy__", [](const OpenMS::PeptideHit::PeakAnnotation& self) { return OpenMS::PeptideHit::PeakAnnotation(self); })
        .def("__deepcopy__", [](const OpenMS::PeptideHit::PeakAnnotation& self, nb::dict) { return OpenMS::PeptideHit::PeakAnnotation(self); }, "memo"_a)
        .def(nb::self < nb::self)
        .def(nb::self == nb::self)
        .def_static("writePeakAnnotationsString_", [](OpenMS::String& annotation_string, std::vector<OpenMS::PeptideHit::PeakAnnotation> annotations) { return OpenMS::PeptideHit::PeakAnnotation::writePeakAnnotationsString_(annotation_string, annotations); }, "annotation_string"_a, "annotations"_a)
        .def_rw("annotation", &OpenMS::PeptideHit::PeakAnnotation::annotation)
        .def_rw("charge", &OpenMS::PeptideHit::PeakAnnotation::charge)
        .def_rw("mz", &OpenMS::PeptideHit::PeakAnnotation::mz)
        .def_rw("intensity", &OpenMS::PeptideHit::PeakAnnotation::intensity)
        ;

    // -----------------------------------------------------------------------
    // PeptideIdentificationList
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeptideIdentificationList>(m, "PeptideIdentificationList", 
        R"doc(
A container for peptide identifications from multiple spectra.
This class provides a vector-like interface for managing collections
of peptide identifications, typically obtained from tandem mass
spectrometry experiments or peptide database searches. Each PeptideIdentification represents
the identification results from a single spectrum.
This class supports direct iteration in Python.
)doc")
        .def(nb::init<std::vector<OpenMS::PeptideIdentification>>())
        .def(nb::init<std::initializer_list<OpenMS::PeptideIdentification>>())

        .def(nb::init<>(), "Default constructor - creates an empty list")

        .def("size", [](const OpenMS::PeptideIdentificationList& self) {
            return self.size();
        }, "Returns the number of peptide identifications")

        .def("__getitem__", [](OpenMS::PeptideIdentificationList& self, size_t i) -> OpenMS::PeptideIdentification& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)

        .def("append", [](OpenMS::PeptideIdentificationList& self, const OpenMS::PeptideIdentification& id) {
            self.push_back(id);
        }, "id"_a, "Add a peptide identification (alias for push_back)")

        .def("extend", [](OpenMS::PeptideIdentificationList& self, const OpenMS::PeptideIdentificationList& other) {
            for (const auto& id : other) self.push_back(id);
        }, "other"_a, "Extend with items from another PeptideIdentificationList")
        .def("extend", [](OpenMS::PeptideIdentificationList& self, const nb::list& items) {
            for (auto item : items) {
                self.push_back(nb::cast<OpenMS::PeptideIdentification>(item));
            }
        }, "items"_a, "Extend with items from a Python list")

        .def("push_back", [](OpenMS::PeptideIdentificationList& self, const OpenMS::PeptideIdentification& id) {
            self.push_back(id);
        }, "id"_a, "Add a peptide identification")

        .def("clear", [](OpenMS::PeptideIdentificationList& self) {
            self.clear();
        }, "Clear all identifications")

        .def("empty", [](const OpenMS::PeptideIdentificationList& self) {
            return self.empty();
        }, "Returns True if the list is empty")

        .def("__len__", [](const OpenMS::PeptideIdentificationList& self) {
            return self.size();
        })
        .def("__iter__", [](OpenMS::PeptideIdentificationList& self) {
            return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::PeptideIdentificationList>(), "PeptideIdentificationList_iter", self.begin(), self.end());
        })
        .def("__setitem__", [](OpenMS::PeptideIdentificationList& self, size_t i, const OpenMS::PeptideIdentification& val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        }, "i"_a, "val"_a)
        .def("at", [](OpenMS::PeptideIdentificationList& self, size_t i) -> OpenMS::PeptideIdentification& {
            return self.at(i);
        }, "i"_a, nb::rv_policy::reference_internal, "Returns reference to element at index (with bounds checking)")
        .def("front", [](OpenMS::PeptideIdentificationList& self) -> OpenMS::PeptideIdentification& {
            return self.front();
        }, nb::rv_policy::reference_internal, "Returns reference to first element")
        .def("back", [](OpenMS::PeptideIdentificationList& self) -> OpenMS::PeptideIdentification& {
            return self.back();
        }, nb::rv_policy::reference_internal, "Returns reference to last element")
        .def("__repr__", [](const OpenMS::PeptideIdentificationList& self) {
            return "PeptideIdentificationList(size=" + std::to_string(self.size()) + ")";
        })
        .def("__str__", [](const OpenMS::PeptideIdentificationList& self) { return nb::cast(self).attr("__repr__")(); })
        ;

    // -----------------------------------------------------------------------
    // ProteinGroup
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProteinIdentification::ProteinGroup>(m, "ProteinGroup", "ProteinIdentification")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProteinIdentification::ProteinGroup &>())
        .def("__copy__", [](const OpenMS::ProteinIdentification::ProteinGroup& self) { return OpenMS::ProteinIdentification::ProteinGroup(self); })
        .def("__deepcopy__", [](const OpenMS::ProteinIdentification::ProteinGroup& self, nb::dict) { return OpenMS::ProteinIdentification::ProteinGroup(self); }, "memo"_a)
        .def_rw("probability", &OpenMS::ProteinIdentification::ProteinGroup::probability)
        .def_rw("accessions", &OpenMS::ProteinIdentification::ProteinGroup::accessions)
        ;

    // -----------------------------------------------------------------------
    // SpectrumLookup
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumLookup>(m, "SpectrumLookup", "Helper class for looking up spectra based on different attributes")
        .def(nb::init<>())
        .def("empty", [](const OpenMS::SpectrumLookup& self) { return self.empty(); }, "Check if any spectra were set")
        .def("findByRT", [](const OpenMS::SpectrumLookup& self, double rt) { return self.findByRT(rt); }, "rt"_a, 
            R"doc(
Read and index spectra for later look-up
:param spectra: Container of spectra
:param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>". For example, "scan=(?<SCAN>\\d+)").
)doc")
        .def("findByNativeID", [](const OpenMS::SpectrumLookup& self, const OpenMS::String& native_id) { return self.findByNativeID(native_id); }, "native_id"_a, 
            R"doc(
Look up spectrum by retention time (RT)
:param rt: Retention time to look up
:returns: Index of the spectrum that matched
)doc")
        .def("findByIndex", [](const OpenMS::SpectrumLookup& self, size_t index, bool count_from_one) { return self.findByIndex(index, count_from_one); }, "index"_a, "count_from_one"_a = false, 
            R"doc(
Look up spectrum by native ID
:param native_id: Native ID to look up
:returns: Index of the spectrum that matched
Size findByIndex(Size index)
)doc")
        .def("findByScanNumber", [](const OpenMS::SpectrumLookup& self, size_t scan_number) { return self.findByScanNumber(scan_number); }, "scan_number"_a, 
            R"doc(
Look up spectrum by index (position in the vector of spectra)
:param index: Index to look up
:param count_from_one: Do indexes start counting at one (default zero)?
:returns: Index of the spectrum that matched
)doc")
        .def("findByReference", [](const OpenMS::SpectrumLookup& self, const OpenMS::String& spectrum_ref) { return self.findByReference(spectrum_ref); }, "spectrum_ref"_a, 
            R"doc(
Look up spectrum by scan number (extracted from the native ID)
:param scan_number: Scan number to look up
:returns: Index of the spectrum that matched
)doc")
        .def("addReferenceFormat", [](OpenMS::SpectrumLookup& self, const OpenMS::String& regexp) { return self.addReferenceFormat(regexp); }, "regexp"_a, 
            R"doc(
Look up spectrum by reference
:param spectrum_ref: Spectrum reference to parse
:returns: Index of the spectrum that matched
)doc")
        .def_static("extractScanNumber", [](const OpenMS::String& native_id, const boost::basic_regex<char>& scan_regexp, bool no_error) { return OpenMS::SpectrumLookup::extractScanNumber(native_id, scan_regexp, no_error); }, "native_id"_a, "scan_regexp"_a, "no_error"_a, 
            R"doc(
Extract scan number from a native ID using a regular expression
:param native_id: The native spectrum ID string
:param scan_regexp: Regular expression for extracting the scan number
:param no_error: If true, do not throw on failure
)doc")
        .def_static("extractScanNumber", [](const OpenMS::String& native_id, const OpenMS::String& native_id_type_accession) { return OpenMS::SpectrumLookup::extractScanNumber(native_id, native_id_type_accession); }, "native_id"_a, "native_id_type_accession"_a,
            R"doc(
Extract scan number from a native ID using the accession type
:param native_id: The native spectrum ID string
:param native_id_type_accession: The native ID type accession
)doc")
        .def_rw("rt_tolerance", &OpenMS::SpectrumLookup::rt_tolerance)
        .def("readSpectra", [](OpenMS::SpectrumLookup& self, const OpenMS::MSExperiment& spectra, const OpenMS::String& scan_regexp) { self.readSpectra(spectra, scan_regexp); }, "spectra"_a, "scan_regexp"_a = OpenMS::SpectrumLookup::default_scan_regexp, "Read and index spectra for later look-up")
        ;

    // -----------------------------------------------------------------------
    // SpectrumNativeIDParser
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumNativeIDParser>(m, "SpectrumNativeIDParser", "OpenMS class SpectrumNativeIDParser")
        .def_static("extractScanNumber", [](const OpenMS::String& native_id, const boost::basic_regex<char>& scan_regexp, bool no_error) { return OpenMS::SpectrumNativeIDParser::extractScanNumber(native_id, scan_regexp, no_error); }, "native_id"_a, "scan_regexp"_a, "no_error"_a, 
            R"doc(
wrap-attach:
SpectrumNativeIDParser
)doc")
        .def_static("extractScanNumber", [](const OpenMS::String& native_id, const OpenMS::String& native_id_type_accession) { return OpenMS::SpectrumNativeIDParser::extractScanNumber(native_id, native_id_type_accession); }, "native_id"_a, "native_id_type_accession"_a, 
            R"doc(
wrap-attach:
SpectrumNativeIDParser
)doc")
        .def_static("getRegExFromNativeID", [](const OpenMS::String& native_id) { return OpenMS::SpectrumNativeIDParser::getRegExFromNativeID(native_id); }, "native_id"_a, 
            R"doc(
wrap-attach:
SpectrumNativeIDParser
)doc")
        .def_static("isNativeID", [](const OpenMS::String& id) { return OpenMS::SpectrumNativeIDParser::isNativeID(id); }, "id"_a, 
            R"doc(
wrap-attach:
SpectrumNativeIDParser
)doc")
        ;

    // -----------------------------------------------------------------------
    // USI
    // -----------------------------------------------------------------------
    auto usi_class = nb::class_<OpenMS::USI>(m, "USI", 
        R"doc(
Utility class for handling Universal Spectrum Identifiers (USI).
USI format (PSI-MS MS:1003063):
mzspec:<collection>:<ms_run>:<index_type>:<index>[:interpretation]
The optional interpretation part uses ProForma proteoform-ion notation.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::USI &>())
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::USI::IndexType, OpenMS::String, OpenMS::String>())
        .def(nb::init<OpenMS::String>())
        .def("__copy__", [](const OpenMS::USI& self) { return OpenMS::USI(self); })
        .def("__deepcopy__", [](const OpenMS::USI& self, nb::dict) { return OpenMS::USI(self); }, "memo"_a)
        .def("isValid", [](const OpenMS::USI& self) { return self.isValid(); }, "Return True if required fields are set")
        .def_static("isValidUSI", [](const OpenMS::String& usi_string) { return OpenMS::USI::isValidUSI(usi_string); }, "usi_string"_a, "Validate a USI string format")
        .def("getCollection", [](const OpenMS::USI& self) { return self.getCollection(); }, "Get the dataset/library identifier")
        .def("setCollection", [](OpenMS::USI& self, const OpenMS::String& collection) { return self.setCollection(collection); }, "collection"_a, "Set the dataset/library identifier")
        .def("getMSRun", [](const OpenMS::USI& self) { return self.getMSRun(); }, "Get the MS run file name")
        .def("setMSRun", [](OpenMS::USI& self, const OpenMS::String& ms_run) { return self.setMSRun(ms_run); }, "ms_run"_a, "Set the MS run file name")
        .def("getIndexType", [](const OpenMS::USI& self) { return self.getIndexType(); }, "Get the index type (scan/index/nativeId)")
        .def("setIndexType", [](OpenMS::USI& self, OpenMS::USI::IndexType index_type) { return self.setIndexType(index_type); }, "index_type"_a, "Set the index type")
        .def("getIndex", [](const OpenMS::USI& self) { return self.getIndex(); }, "Get the spectrum index value")
        .def("setIndex", [](OpenMS::USI& self, const OpenMS::String& index) { return self.setIndex(index); }, "index"_a, "Set the spectrum index value")
        .def("getInterpretation", [](const OpenMS::USI& self) { return self.getInterpretation(); }, "Get the optional ProForma interpretation")
        .def("setInterpretation", [](OpenMS::USI& self, const OpenMS::String& interpretation) { return self.setInterpretation(interpretation); }, "interpretation"_a, "Set the optional ProForma interpretation")
        .def("hasInterpretation", [](const OpenMS::USI& self) { return self.hasInterpretation(); }, "Return True if interpretation is present")
        .def("toString", [](const OpenMS::USI& self) { return self.toString(); }, "Convert this USI to its string representation (empty if invalid)")
        .def("fromString", [](OpenMS::USI& self, const OpenMS::String& usi_string) { return self.fromString(usi_string); }, "usi_string"_a, "Parse a USI string into this object")
        .def_static("indexTypeToString", [](OpenMS::USI::IndexType index_type) { return OpenMS::USI::indexTypeToString(index_type); }, "index_type"_a, "Convert index type enum to string")
        .def_static("indexTypeFromString", [](const OpenMS::String& type_string) { return OpenMS::USI::indexTypeFromString(type_string); }, "type_string"_a, "Parse index type from string")
        .def_static("extractBasename", [](const OpenMS::String& filepath) { return OpenMS::USI::extractBasename(filepath); }, "filepath"_a, "Extract basename from file path/URI for use as ms_run")
        .def_static("getCVAccession", []() { return OpenMS::USI::getCVAccession(); }, "Get PSI-MS CV accession for USI (MS:1003063)")
        .def_static("getCVName", []() { return OpenMS::USI::getCVName(); }, "Get PSI-MS CV name for USI")
        ;
    // IndexType enum nested under USI
    nb::enum_<OpenMS::USI::IndexType>(usi_class, "IndexType", nb::is_arithmetic())
        .value("SCAN", OpenMS::USI::IndexType::SCAN)
        .value("INDEX", OpenMS::USI::IndexType::INDEX)
        .value("NATIVEID", OpenMS::USI::IndexType::NATIVEID)
        ;


    // -----------------------------------------------------------------------
    // SpectrumMetaData (nested struct of SpectrumMetaDataLookup)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumMetaDataLookup::SpectrumMetaData>(m, "SpectrumMetaData", "Spectrum metadata (RT, precursor info, MS level, etc.)")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::SpectrumMetaDataLookup::SpectrumMetaData& self) { return OpenMS::SpectrumMetaDataLookup::SpectrumMetaData(self); })
        .def("__deepcopy__", [](const OpenMS::SpectrumMetaDataLookup::SpectrumMetaData& self, nb::dict) { return OpenMS::SpectrumMetaDataLookup::SpectrumMetaData(self); }, "memo"_a)
        .def_rw("rt", &OpenMS::SpectrumMetaDataLookup::SpectrumMetaData::rt)
        .def_rw("precursor_rt", &OpenMS::SpectrumMetaDataLookup::SpectrumMetaData::precursor_rt)
        .def_rw("precursor_mz", &OpenMS::SpectrumMetaDataLookup::SpectrumMetaData::precursor_mz)
        .def_rw("precursor_charge", &OpenMS::SpectrumMetaDataLookup::SpectrumMetaData::precursor_charge)
        .def_rw("ms_level", &OpenMS::SpectrumMetaDataLookup::SpectrumMetaData::ms_level)
        .def_rw("scan_number", &OpenMS::SpectrumMetaDataLookup::SpectrumMetaData::scan_number)
        .def_rw("native_id", &OpenMS::SpectrumMetaDataLookup::SpectrumMetaData::native_id)
        ;

    // -----------------------------------------------------------------------
    // SpectrumMetaDataLookup
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumMetaDataLookup, OpenMS::SpectrumLookup>(m, "SpectrumMetaDataLookup",
        R"doc(Helper class for looking up spectrum meta data.

Provides functions for extracting and looking up meta data of spectra,
including retention time, precursor m/z, MS level, scan number, and native ID.
Inherits lookup-by-RT, lookup-by-native-ID, and lookup-by-index from SpectrumLookup.
)doc")
        .def(nb::init<>())
        .def("readSpectra", [](OpenMS::SpectrumMetaDataLookup& self, const OpenMS::MSExperiment& spectra, const OpenMS::String& scan_regexp, bool get_precursor_rt) {
            self.readSpectra(spectra, scan_regexp, get_precursor_rt);
        }, "spectra"_a, "scan_regexp"_a = OpenMS::SpectrumLookup::default_scan_regexp, "get_precursor_rt"_a = false,
            R"doc(Read spectra and store their meta data for later look-up.

:param spectra: MSExperiment containing the spectra
:param scan_regexp: Regular expression for matching scan numbers in spectrum native IDs
:param get_precursor_rt: Assign precursor retention times?
)doc")
        .def("findByNativeID", [](const OpenMS::SpectrumMetaDataLookup& self, const OpenMS::String& native_id) {
            return self.findByNativeID(native_id);
        }, "native_id"_a,
            R"doc(Look up spectrum by native ID.

:param native_id: Native ID to look up
:returns: Index of the spectrum that matched
)doc")
        .def("findByRT", [](const OpenMS::SpectrumMetaDataLookup& self, double rt) {
            return self.findByRT(rt);
        }, "rt"_a,
            R"doc(Look up spectrum by retention time (RT).

:param rt: Retention time to look up
:returns: Index of the spectrum that matched
)doc")
        .def("findByIndex", [](const OpenMS::SpectrumMetaDataLookup& self, size_t index, bool count_from_one) {
            return self.findByIndex(index, count_from_one);
        }, "index"_a, "count_from_one"_a = false,
            R"doc(Look up spectrum by index (position in the vector of spectra).

:param index: Index to look up
:param count_from_one: Do indexes start counting at one (default: zero)?
:returns: Index of the spectrum that matched
)doc")
        .def("findByScanNumber", [](const OpenMS::SpectrumMetaDataLookup& self, size_t scan_number) {
            return self.findByScanNumber(scan_number);
        }, "scan_number"_a,
            R"doc(Look up spectrum by scan number (extracted from the native ID).

:param scan_number: Scan number to look up
:returns: Index of the spectrum that matched
)doc")
        .def("findByReference", [](const OpenMS::SpectrumMetaDataLookup& self, const OpenMS::String& spectrum_ref) {
            return self.findByReference(spectrum_ref);
        }, "spectrum_ref"_a,
            R"doc(Look up spectrum by reference string.

:param spectrum_ref: Spectrum reference to parse
:returns: Index of the spectrum that matched
)doc")
        .def("addReferenceFormat", [](OpenMS::SpectrumMetaDataLookup& self, const OpenMS::String& regexp) {
            self.addReferenceFormat(regexp);
        }, "regexp"_a,
            R"doc(Register a possible format for a spectrum reference.

:param regexp: Regular expression defining the format
)doc")
        .def("empty", [](const OpenMS::SpectrumMetaDataLookup& self) { return self.empty(); }, "Check if any spectra were set")
        .def("setSpectraDataRef", [](OpenMS::SpectrumMetaDataLookup& self, const OpenMS::String& spectra_data) {
            self.setSpectraDataRef(spectra_data);
        }, "spectra_data"_a, "Set spectra data reference (filename)")
        .def_static("addMissingRTsToPeptideIDs", [](OpenMS::PeptideIdentificationList& peptides, const OpenMS::MSExperiment& exp) -> bool {
            return OpenMS::SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(peptides, exp);
        }, "peptides"_a, "exp"_a, "Add missing RTs to peptide IDs")
        .def_static("addMissingSpectrumReferences", [](OpenMS::PeptideIdentificationList& peptides,
                const std::string& filename, bool stop_on_error, bool override_spectra_data, bool override_spectra_references) -> bool {
            std::vector<OpenMS::ProteinIdentification> proteins;
            return OpenMS::SpectrumMetaDataLookup::addMissingSpectrumReferences(peptides, filename, stop_on_error, override_spectra_data, override_spectra_references, proteins);
        }, "peptides"_a, "filename"_a, "stop_on_error"_a = false, "override_spectra_data"_a = false, "override_spectra_references"_a = false, "Add missing spectrum references")
        .def_static("getSpectrumMetaData", [](const OpenMS::MSSpectrum& spectrum) {
            OpenMS::SpectrumMetaDataLookup::SpectrumMetaData result;
            OpenMS::SpectrumMetaDataLookup::getSpectrumMetaData(spectrum, result);
            return result;
        }, "spectrum"_a, "Extract meta data from a spectrum")
        .def_static("addMissingIMToPeptideIDs", [](OpenMS::PeptideIdentificationList& pep_ids, const OpenMS::MSExperiment& exp) {
            return OpenMS::SpectrumMetaDataLookup::addMissingIMToPeptideIDs(pep_ids, exp);
        }, "peptide_ids"_a, "exp"_a, "Add missing ion mobility information to peptide identifications")
        ;

    // Free function aliases for backward compatibility
    m.def("extractScanNumber", [](const OpenMS::String& native_id, const OpenMS::String& native_id_type_accession) {
        return OpenMS::SpectrumNativeIDParser::extractScanNumber(native_id, native_id_type_accession);
    }, "native_id"_a, "native_id_type_accession"_a, "Extract scan number from native ID string");
    m.def("getRegExFromNativeID", [](const OpenMS::String& native_id) {
        return OpenMS::SpectrumNativeIDParser::getRegExFromNativeID(native_id);
    }, "native_id"_a, "Get regular expression from native ID string");
    m.def("isNativeID", [](const OpenMS::String& id) {
        return OpenMS::SpectrumNativeIDParser::isNativeID(id);
    }, "id"_a, "Check if string is a native ID");

}