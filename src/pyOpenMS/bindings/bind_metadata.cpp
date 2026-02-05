// pyOpenMS nanobind bindings
// Domain: metadata

#include "all_casters.h"
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
    // AbsoluteQuantitationStandards
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AbsoluteQuantitationStandards>(m, "AbsoluteQuantitationStandards", 
        R"doc(
AbsoluteQuantitationStandards is a class to handle the relationship
between runs, components, and their actual concentrations
)doc")
        .def(nb::init<>())
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
    // CVTerm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CVTerm>(m, "CVTerm", "Representation of controlled vocabulary term")
        .def(nb::init<>())
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::String, OpenMS::String, OpenMS::CVTerm::Unit>())
        .def(nb::init<const OpenMS::CVTerm &>())
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
        .def(nb::init<std::vector<OpenMS::ExperimentalDesign::MSFileSectionEntry>, OpenMS::ExperimentalDesign::SampleSection>())
        .def("getNumberOfSamples", [](const OpenMS::ExperimentalDesign& self) { return self.getNumberOfSamples(); }, "Returns the number of samples measured (= highest sample index)")
        .def("getNumberOfFractions", [](const OpenMS::ExperimentalDesign& self) { return self.getNumberOfFractions(); }, "Returns the number of fractions (= highest fraction index)")
        .def("getNumberOfLabels", [](const OpenMS::ExperimentalDesign& self) { return self.getNumberOfLabels(); }, "Returns the number of labels per file")
        .def("getNumberOfMSFiles", [](const OpenMS::ExperimentalDesign& self) { return self.getNumberOfMSFiles(); }, "Returns the number of MS files (= fractions * fraction_groups)")
        .def("getNumberOfFractionGroups", [](const OpenMS::ExperimentalDesign& self) { return self.getNumberOfFractionGroups(); }, "Allows to group fraction ids and source files. Return the number of fraction_groups")
        .def("getSample", [](OpenMS::ExperimentalDesign& self, unsigned int fraction_group, unsigned int label) { return self.getSample(fraction_group, label); }, "fraction_group"_a, "label"_a, "Returns sample index (depends on fraction_group and label)")
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
        ;

    // -----------------------------------------------------------------------
    // ExperimentalDesign_MSFileSectionEntry
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ExperimentalDesign::MSFileSectionEntry>(m, "ExperimentalDesign_MSFileSectionEntry")
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
    nb::class_<OpenMS::ExperimentalDesign::SampleSection>(m, "ExperimentalDesign_SampleSection")
        .def(nb::init<>())
        .def(nb::init<std::vector<std::vector<OpenMS::String>>, std::map<OpenMS::String, unsigned long>, std::map<OpenMS::String, unsigned long>>())
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
Create/update mapping from a list of ProteinIdentifications.
:param prot_ids: List of ProteinIdentification objects
)doc")
        .def("empty", [](const OpenMS::IdentifierMSRunMapper& self) { return self.empty(); }, 
            R"doc(
Check if the mapping contains an entry for the given identifier.
:param identifier: ProteinIdentification identifier
:return: True if identifier exists in mapping
)doc")
        .def("size", [](const OpenMS::IdentifierMSRunMapper& self) { return self.size(); }, 
            R"doc(
Check if the mapping is empty.
:return: True if no mappings exist
)doc")
        .def("getMSRunPaths", [](const OpenMS::IdentifierMSRunMapper& self, const OpenMS::String& identifier) -> const std::vector<OpenMS::String> & { return self.getMSRunPaths(identifier); }, "identifier"_a, nb::rv_policy::reference_internal, 
            R"doc(
Get the number of identifier mappings.
:return: Number of identifiers in the mapping
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
        
        .def("getKeys", [](const OpenMS::MetaInfo& self, nb::list py_keys) {
            std::vector<OpenMS::String> keys;
            self.getKeys(keys);
            py_keys.attr("clear")();
            for (const auto& k : keys) {
                py_keys.append(nb::str(k.c_str()));
            }
        }, "keys"_a, "Fills the given list with all meta value keys")
        .def("empty", [](const OpenMS::MetaInfo& self) { return self.empty(); }, "Returns if the MetaInfo is empty")
        .def("clear", [](OpenMS::MetaInfo& self) { return self.clear(); }, "Removes all meta values")
        .def("__hash__", [](const OpenMS::MetaInfo& self) { return std::hash<OpenMS::MetaInfo>{}(self); })
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
        .def("registerName", [](OpenMS::MetaInfoRegistry& self, const OpenMS::String& name, const OpenMS::String& description, const OpenMS::String& unit) { return self.registerName(name, description, unit); }, "name"_a, "description"_a, "unit"_a, "Registers a string, stores its description and unit, and returns the corresponding index. If the string is already registered, it returns the index of the string")
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
        .def(nb::init<OpenMS::String, int, int, char, char>())
        .def(nb::init<const OpenMS::PeptideEvidence &>())
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
        ;

    // -----------------------------------------------------------------------
    // PeptideHit_PeakAnnotation
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeptideHit::PeakAnnotation>(m, "PeptideHit_PeakAnnotation")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PeptideHit::PeakAnnotation &>())
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
        }, "other"_a, "Extend with items from another list")

        .def("push_back", [](OpenMS::PeptideIdentificationList& self, const OpenMS::PeptideIdentification& id) {
            self.push_back(id);
        }, "id"_a, "Add a peptide identification")

        .def("clear", [](OpenMS::PeptideIdentificationList& self) {
            self.clear();
        }, "Clear all identifications")
        ;

    // -----------------------------------------------------------------------
    // ProteinGroup
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProteinIdentification::ProteinGroup>(m, "ProteinGroup", "ProteinIdentification")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProteinIdentification::ProteinGroup &>())
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
        .def("findByIndex", [](const OpenMS::SpectrumLookup& self, unsigned long index, bool count_from_one) { return self.findByIndex(index, count_from_one); }, "index"_a, "count_from_one"_a, 
            R"doc(
Look up spectrum by native ID
:param native_id: Native ID to look up
:returns: Index of the spectrum that matched
Size findByIndex(Size index)
)doc")
        .def("findByScanNumber", [](const OpenMS::SpectrumLookup& self, unsigned long scan_number) { return self.findByScanNumber(scan_number); }, "scan_number"_a, 
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
Register a possible format for a spectrum reference
:param regexp: Regular expression defining the format
# NAMESPACE # Int extractScanNumber(const String & native_id, boost::regex & scan_regexp, bool no_error) except + nogil
)doc")
        .def_static("extractScanNumber", [](const OpenMS::String& native_id, const OpenMS::String& native_id_type_accession) { return OpenMS::SpectrumLookup::extractScanNumber(native_id, native_id_type_accession); }, "native_id"_a, "native_id_type_accession"_a, 
            R"doc(
Register a possible format for a spectrum reference
:param regexp: Regular expression defining the format
# NAMESPACE # Int extractScanNumber(const String & native_id, boost::regex & scan_regexp, bool no_error) except + nogil
)doc")
        .def_rw("rt_tolerance", &OpenMS::SpectrumLookup::rt_tolerance)
        ;

    // -----------------------------------------------------------------------
    // SpectrumNativeIDParser
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumNativeIDParser>(m, "SpectrumNativeIDParser")
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
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::USI::IndexType, OpenMS::String, OpenMS::String>())
        .def(nb::init<OpenMS::String>())
        .def(nb::init<const OpenMS::USI &>())
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


    struct SpectrumMetaDataLookup_Dummy {};

    auto smld_addMissingRTs = [](OpenMS::PeptideIdentificationList& peptides, const OpenMS::MSExperiment& exp) -> bool {
        return OpenMS::SpectrumMetaDataLookup::addMissingRTsToPeptideIDs(peptides, exp);
    };

    auto smld_addMissingRefs = [](OpenMS::PeptideIdentificationList& peptides,
            const std::string& filename, bool stop_on_error, bool override_spectra_data, bool override_spectra_references) -> bool {
        std::vector<OpenMS::ProteinIdentification> proteins;
        return OpenMS::SpectrumMetaDataLookup::addMissingSpectrumReferences(peptides, filename, stop_on_error, override_spectra_data, override_spectra_references, proteins);
    };

    // -----------------------------------------------------------------------
    // SpectrumMetaDataLookup
    // -----------------------------------------------------------------------
    nb::class_<SpectrumMetaDataLookup_Dummy>(m, "SpectrumMetaDataLookup")
        .def_static("addMissingRTsToPeptideIDs", smld_addMissingRTs, "peptides"_a, "exp"_a, "Add missing RTs to peptide IDs")
        .def_static("addMissingSpectrumReferences", smld_addMissingRefs, "peptides"_a, "filename"_a, "stop_on_error"_a = false, "override_spectra_data"_a = false, "override_spectra_references"_a = false, "Add missing spectrum references")
        ;

}