// Template helpers for binding common interface methods
// Used by Category B classes (no nanobind base declared due to multiple inheritance)
// to avoid duplicating MetaInfoInterface, UniqueIdInterface, DocumentIdentifier,
// DefaultParamHandler, and CVTermList method bindings.

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FileTypes.h>

namespace nb = nanobind;
using namespace nb::literals;

/// Bind MetaInfoInterface methods (~9 methods) onto a nanobind class.
/// Use for classes that inherit MetaInfoInterface but can't declare it as
/// a nanobind base (e.g. due to multiple inheritance).
template <typename T, typename Class>
void def_MetaInfoInterface(Class& cls)
{
    cls
        .def("getMetaValue", [](const T& self, const OpenMS::String& name) {
            return self.getMetaValue(name);
        }, "name"_a, "Returns the value corresponding to a string, or DataValue::EMPTY if not found")
        .def("metaValueExists", [](const T& self, const OpenMS::String& name) {
            return self.metaValueExists(name);
        }, "name"_a, "Returns whether an entry with the given name exists")
        .def("setMetaValue", [](T& self, const OpenMS::String& name, const OpenMS::DataValue& value) {
            return self.setMetaValue(name, value);
        }, "name"_a, "value"_a, "Sets the DataValue corresponding to a name")
        .def("removeMetaValue", [](T& self, const OpenMS::String& name) {
            return self.removeMetaValue(name);
        }, "name"_a, "Removes the DataValue corresponding to `name` if it exists")
        .def_static("metaRegistry", []() {
            return T::metaRegistry();
        }, "Returns a reference to the MetaInfoRegistry")
        .def("getKeys", [](const T& self, nb::list py_keys) {
            std::vector<OpenMS::String> keys;
            self.getKeys(keys);
            py_keys.attr("clear")();
            for (const auto& k : keys) {
                py_keys.append(nb::str(k.c_str()));
            }
        }, "keys"_a, "Fills the given list with all meta value keys")
        .def("getKeys", [](const T& self) {
            std::vector<OpenMS::String> keys;
            self.getKeys(keys);
            nb::list result;
            for (const auto& k : keys) {
                result.append(nb::str(k.c_str()));
            }
            return result;
        }, "Returns all meta value keys as a list")
        .def("isMetaEmpty", [](const T& self) {
            return self.isMetaEmpty();
        }, "Returns if the MetaInfo is empty")
        .def("clearMetaInfo", [](T& self) {
            return self.clearMetaInfo();
        }, "Removes all meta values")
        ;
}

/// Bind UniqueIdInterface methods (~7 methods) onto a nanobind class.
template <typename T, typename Class>
void def_UniqueIdInterface(Class& cls)
{
    cls
        .def_static("isValid", [](size_t unique_id) {
            return T::isValid(unique_id);
        }, "unique_id"_a, "Returns true if the unique_id is valid, false otherwise")
        .def("getUniqueId", [](const T& self) {
            return self.getUniqueId();
        }, "Returns the unique id")
        .def("clearUniqueId", [](T& self) {
            return self.clearUniqueId();
        }, "Clear the unique id. The new unique id will be invalid.  Returns 1 if the unique id was changed, 0 otherwise")
        .def("hasValidUniqueId", [](const T& self) {
            return self.hasValidUniqueId();
        }, "Returns whether the unique id is valid.  Returns 1 if the unique id is valid, 0 otherwise")
        .def("hasInvalidUniqueId", [](const T& self) {
            return self.hasInvalidUniqueId();
        }, "Returns whether the unique id is invalid.  Returns 1 if the unique id is invalid, 0 otherwise")
        .def("setUniqueId", [](T& self) {
            return self.setUniqueId();
        }, "Assigns a new, valid unique id.  Always returns 1")
        .def("setUniqueId", [](T& self, size_t rhs) {
            self.setUniqueId(rhs);
        }, "rhs"_a, "Assigns the given unique id")
        .def("setUniqueId", [](T& self, const OpenMS::String& rhs) {
            self.setUniqueId(rhs);
        }, "rhs"_a, "Assigns a unique id from a string (extracts UInt64 after last underscore)")
        .def("ensureUniqueId", [](T& self) {
            return self.ensureUniqueId();
        }, "Assigns a valid unique id, but only if the present one is invalid.  Returns 1 if the unique id was changed, 0 otherwise")
        ;
}

/// Bind DocumentIdentifier methods (6 methods) onto a nanobind class.
template <typename T, typename Class>
void def_DocumentIdentifier(Class& cls)
{
    cls
        .def("setIdentifier", [](T& self, const OpenMS::String& id) {
            return self.setIdentifier(id);
        }, "id"_a, "Sets document identifier (e.g. an LSID)")
        .def("getIdentifier", [](const T& self) {
            return self.getIdentifier();
        }, "Retrieve document identifier (e.g. an LSID)")
        .def("setLoadedFilePath", [](T& self, const OpenMS::String& file_name) {
            return self.setLoadedFilePath(file_name);
        }, "file_name"_a, "Sets the file_name according to absolute path of the file loaded, preferably done whilst loading")
        .def("getLoadedFilePath", [](const T& self) {
            return self.getLoadedFilePath();
        }, "Returns the file_name which is the absolute path to the file loaded")
        .def("setLoadedFileType", [](T& self, const OpenMS::String& file_name) {
            return self.setLoadedFileType(file_name);
        }, "file_name"_a, "Sets the file_type according to the type of the file loaded from, preferably done whilst loading")
        .def("getLoadedFileType", [](const T& self) -> const OpenMS::FileTypes::Type & {
            return self.getLoadedFileType();
        }, nb::rv_policy::reference_internal, "Returns the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded")
        ;
}

/// Bind DefaultParamHandler methods (6 methods) onto a nanobind class.
template <typename T, typename Class>
void def_DefaultParamHandler(Class& cls)
{
    cls
        .def("setParameters", [](T& self, const OpenMS::Param& param) {
            return self.setParameters(param);
        }, "param"_a, "Sets the parameters")
        .def("getParameters", [](const T& self) -> const OpenMS::Param & {
            return self.getParameters();
        }, nb::rv_policy::reference_internal, "Returns the parameters")
        .def("getDefaults", [](const T& self) -> const OpenMS::Param & {
            return self.getDefaults();
        }, nb::rv_policy::reference_internal, "Returns the default parameters")
        .def("getName", [](const T& self) {
            return self.getName();
        }, "Returns the name")
        .def("setName", [](T& self, const OpenMS::String& name) {
            return self.setName(name);
        }, "name"_a, "Sets the name")
        .def("getSubsections", [](const T& self) -> const std::vector<OpenMS::String> & {
            return self.getSubsections();
        }, nb::rv_policy::reference_internal)
        ;
}

/// Bind ProgressLogger methods (6 methods) onto a nanobind class.
/// Use for classes that inherit ProgressLogger but can't declare it as
/// a nanobind base (e.g. due to multiple inheritance with DefaultParamHandler).
template <typename T, typename Class>
void def_ProgressLogger(Class& cls)
{
    cls
        .def("setLogType", [](const T& self, OpenMS::ProgressLogger::LogType type) {
            return self.setLogType(type);
        }, "type"_a, "Sets the progress log that should be used. The default type is NONE!")
        .def("getLogType", [](const T& self) {
            return self.getLogType();
        }, "Returns the type of progress log being used")
        .def("startProgress", [](const T& self, long begin, long end, const OpenMS::String& label) {
            return self.startProgress(begin, end, label);
        }, "begin"_a, "end"_a, "label"_a)
        .def("setProgress", [](const T& self, long value) {
            return self.setProgress(value);
        }, "value"_a, "Sets the current progress")
        .def("endProgress", [](const T& self, size_t bytes_processed) {
            return self.endProgress(bytes_processed);
        }, "bytes_processed"_a = 0, "Ends the progress display")
        .def("nextProgress", [](const T& self) {
            return self.nextProgress();
        }, "Increment progress by 1 (according to range begin-end)")
        ;
}

/// Bind CVTermList methods (~8 methods) onto a nanobind class.
/// Use for classes that inherit CVTermList but can't declare it as
/// a nanobind base (e.g. Compound, Peptide via PeptideCompound).
template <typename T, typename Class>
void def_CVTermList(Class& cls)
{
    cls
        .def("replaceCVTerm", [](T& self, const OpenMS::CVTerm& cv_term) {
            return self.replaceCVTerm(cv_term);
        }, "cv_term"_a, "Replaces the specified CV term")
        .def("replaceCVTerms", [](T& self, const std::map<OpenMS::String, std::vector<OpenMS::CVTerm>>& cv_term_map) {
            self.replaceCVTerms(cv_term_map);
        }, "cv_term_map"_a, "Replaces all CV terms with the given map")
        .def("consumeCVTerms", [](T& self, const std::map<OpenMS::String, std::vector<OpenMS::CVTerm>>& cv_term_map) {
            return self.consumeCVTerms(cv_term_map);
        }, "cv_term_map"_a, "Merges the given map into the member map, no duplicate checking")
        .def("getCVTerms", [](const T& self) -> const std::map<OpenMS::String, std::vector<OpenMS::CVTerm>> & {
            return self.getCVTerms();
        }, nb::rv_policy::reference_internal, "Returns the accession string of the term")
        .def("addCVTerm", [](T& self, const OpenMS::CVTerm& term) {
            return self.addCVTerm(term);
        }, "term"_a, "Adds a CV term")
        .def("setCVTerms", [](T& self, const std::vector<OpenMS::CVTerm>& terms) {
            return self.setCVTerms(terms);
        }, "terms"_a, "Sets the CV terms from a vector")
        .def("hasCVTerm", [](const T& self, const OpenMS::String& accession) {
            return self.hasCVTerm(accession);
        }, "accession"_a)
        .def("empty", [](const T& self) {
            return self.empty();
        })
        ;
}
