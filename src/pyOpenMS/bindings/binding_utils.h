// Template helpers for binding common interface methods
// Used by Category B classes (no nanobind base declared due to multiple inheritance)
// to avoid duplicating MetaInfoInterface, UniqueIdInterface, DocumentIdentifier,
// DefaultParamHandler, and CVTermList method bindings.

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <algorithm>
#include <optional>
#include <stdexcept>
#include <string>
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
        .def("getMetaValue", [](const T& self, const std::string& name) {
            return self.getMetaValue(name);
        }, "name"_a, "Returns the value corresponding to a string, or DataValue::EMPTY if not found")
        .def("metaValueExists", [](const T& self, const std::string& name) {
            return self.metaValueExists(name);
        }, "name"_a, "Returns whether an entry with the given name exists")
        .def("setMetaValue", [](T& self, const std::string& name, const OpenMS::DataValue& value) {
            return self.setMetaValue(name, value);
        }, "name"_a, "value"_a, "Sets the DataValue corresponding to a name")
        .def("removeMetaValue", [](T& self, const std::string& name) {
            return self.removeMetaValue(name);
        }, "name"_a, "Removes the DataValue corresponding to `name` if it exists")
        .def_static("metaRegistry", []() {
            return T::metaRegistry();
        }, "Returns a reference to the MetaInfoRegistry")
        .def("getKeys", [](const T& self, nb::list py_keys) {
            std::vector<std::string> keys;
            self.getKeys(keys);
            py_keys.attr("clear")();
            for (const auto& k : keys) {
                py_keys.append(nb::str(k.c_str()));
            }
        }, "keys"_a, "Fills the given list with all meta value keys")
        .def("getKeys", [](const T& self) {
            std::vector<std::string> keys;
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
        .def("setUniqueId", [](T& self, const std::string& rhs) {
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
        .def("setIdentifier", [](T& self, const std::string& id) {
            return self.setIdentifier(id);
        }, "id"_a, "Sets document identifier (e.g. an LSID)")
        .def("getIdentifier", [](const T& self) {
            return self.getIdentifier();
        }, "Retrieve document identifier (e.g. an LSID)")
        .def("setLoadedFilePath", [](T& self, const std::string& file_name) {
            return self.setLoadedFilePath(file_name);
        }, "file_name"_a, "Sets the file_name according to absolute path of the file loaded, preferably done whilst loading")
        .def("getLoadedFilePath", [](const T& self) {
            return self.getLoadedFilePath();
        }, "Returns the file_name which is the absolute path to the file loaded")
        .def("setLoadedFileType", [](T& self, const std::string& file_name) {
            return self.setLoadedFileType(file_name);
        }, "file_name"_a, "Sets the file_type according to the type of the file loaded from, preferably done whilst loading")
        .def("getLoadedFileType", [](const T& self) -> OpenMS::FileTypes::Type {
            return self.getLoadedFileType();
        }, "Returns the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded")
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
        .def("getParameters", [](const T& self) -> OpenMS::Param {
            return self.getParameters();
        }, "Returns a copy of the parameters")
        .def("getDefaults", [](const T& self) -> OpenMS::Param {
            return self.getDefaults();
        }, "Returns a copy of the default parameters")
        .def("getName", [](const T& self) {
            return self.getName();
        }, "Returns the name")
        .def("setName", [](T& self, const std::string& name) {
            return self.setName(name);
        }, "name"_a, "Sets the name")
        .def("getSubsections", [](const T& self) -> std::vector<std::string> {
            return self.getSubsections();
        })
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
        .def("startProgress", [](const T& self, long begin, long end, const std::string& label) {
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
        .def("replaceCVTerms", [](T& self, const std::map<std::string, std::vector<OpenMS::CVTerm>>& cv_term_map) {
            self.replaceCVTerms(cv_term_map);
        }, "cv_term_map"_a, "Replaces all CV terms with the given map")
        .def("consumeCVTerms", [](T& self, const std::map<std::string, std::vector<OpenMS::CVTerm>>& cv_term_map) {
            return self.consumeCVTerms(cv_term_map);
        }, "cv_term_map"_a, "Merges the given map into the member map, no duplicate checking")
        .def("getCVTerms", [](const T& self) -> std::map<std::string, std::vector<OpenMS::CVTerm>> {
            return self.getCVTerms();
        }, "Returns the accession string of the term")
        .def("addCVTerm", [](T& self, const OpenMS::CVTerm& term) {
            return self.addCVTerm(term);
        }, "term"_a, "Adds a CV term")
        .def("setCVTerms", [](T& self, const std::vector<OpenMS::CVTerm>& terms) {
            return self.setCVTerms(terms);
        }, "terms"_a, "Sets the CV terms from a vector")
        .def("hasCVTerm", [](const T& self, const std::string& accession) {
            return self.hasCVTerm(accession);
        }, "accession"_a)
        .def("empty", [](const T& self) {
            return self.empty();
        })
        ;
}

/// Trailing paragraph shared by the six set_peaks() docstrings, appended by literal concatenation.
#define PYOPENMS_SET_PEAKS_METADATA_DOC                                                                        \
    "\n\nThe binary data arrays are not resized. ``metadata`` decides what happens when the new peak count "    \
    "would leave a non-empty float, string or integer data array mis-sized: ``error`` (the default) raises "    \
    "and leaves the container untouched, ``clear`` drops the arrays that no longer fit, and ``keep`` keeps "    \
    "them, leaving the container internally inconsistent. Arrays that are already the right length, and empty " \
    "ones, are kept in every case, so replacing peaks with an equally long set never touches them."

/// Trailing paragraph for the two MSSpectrum set_peaks() docstrings, which also accept ion mobility.
#define PYOPENMS_SET_PEAKS_IM_DOC                                                                                \
    "\n\n``ion_mobility`` sets the per-peak ion mobility array in the same call, replacing any existing one; "    \
    "it must have one entry per peak. ``ion_mobility_unit`` names the array through the PSI-MS controlled "       \
    "vocabulary and may be omitted only when the spectrum already has an ion mobility array to inherit the "      \
    "unit from. Only ``DriftTimeUnit.MILLISECOND`` (drift time) and ``DriftTimeUnit.VSSC`` (1/K0) have a CV "     \
    "term for this, so other units are rejected. Because the ion mobility array is replaced rather than "         \
    "stranded, it is exempt from the ``metadata`` check above."

/// What set_peaks() should do with binary data arrays that the new peak count leaves mis-sized.
enum class PeakMetadataPolicy
{
    Error, ///< refuse the call and leave the container untouched
    Clear, ///< drop all float/string/integer data arrays
    Keep   ///< keep them as they are, leaving the container internally inconsistent
};

/// Translate the Python-facing `metadata=` string into a PeakMetadataPolicy.
/// @throws std::invalid_argument (ValueError in Python) for any other string.
inline PeakMetadataPolicy parsePeakMetadataPolicy(const std::string& policy)
{
    if (policy == "error") return PeakMetadataPolicy::Error;
    if (policy == "clear") return PeakMetadataPolicy::Clear;
    if (policy == "keep") return PeakMetadataPolicy::Keep;
    throw std::invalid_argument("Invalid metadata policy '" + policy + "'. Must be 'error', 'clear' or 'keep'.");
}

/// Is @p array still usable as a per-peak annotation for a peak count of @p n?
///
/// An empty array carries no per-peak association and so can never be mis-associated; it is
/// exempt, matching the `!arrays[i].empty()` guard in MSChromatogram::checkDataArraySizes_().
template <typename ArrayT>
inline bool peakMetadataArrayIsAligned(const ArrayT& array, size_t n)
{
    return array.empty() || array.size() == n;
}

/// Reject a set_peaks() that would strand a data array. Only PeakMetadataPolicy::Error throws.
///
/// Call this *before* the peak storage is resized, so that a rejected call leaves the container
/// exactly as it was. @p n is therefore the new, incoming count -- checking against self.size()
/// would compare the arrays against the count they are already aligned with and accept everything.
///
/// The predicate is the one MSChromatogram::checkDataArraySizes_() applies, so this moves the
/// core's own rejection to the call that causes it rather than adding a stricter rule. That helper
/// cannot be reused directly: it is protected, MSSpectrum does not have it, and it can only check
/// against the current count.
///
/// Works for any container exposing get{Float,String,Integer}DataArrays(). @p container is the
/// lowercase noun used in the message ("spectrum", "chromatogram", "mobilogram").
///
/// @p replaced_float_index, when set, names a float array the caller is about to overwrite --
/// set_peaks(..., ion_mobility=...) replaces the ion mobility array wholesale. Its current length
/// is therefore irrelevant: you cannot strand an array you are overwriting, and rejecting the call
/// would leave no way to swap peaks and their ion mobilities in one step.
template <typename Container>
void checkPeakMetadataAlignment(const Container& self, size_t n, PeakMetadataPolicy policy, const char* container,
                                std::optional<size_t> replaced_float_index = std::nullopt)
{
    if (policy != PeakMetadataPolicy::Error) return;

    auto check = [n, container, replaced_float_index](const auto& arrays, const char* what,
                                                      bool honour_replacement = false) {
        for (size_t i = 0; i < arrays.size(); ++i)
        {
            if (honour_replacement && replaced_float_index && *replaced_float_index == i) continue;
            if (!peakMetadataArrayIsAligned(arrays[i], n))
            {
                throw std::runtime_error(std::string(what) + "[" + std::to_string(i) + "] size (" +
                                         std::to_string(arrays[i].size()) + ") does not match the new " + container +
                                         " size (" + std::to_string(n) +
                                         "). set_peaks() does not resize binary data arrays; pass metadata='clear' to "
                                         "drop them or metadata='keep' to keep them and realign them yourself.");
            }
        }
    };
    check(self.getFloatDataArrays(), "FloatDataArray", /* honour_replacement = */ true);
    check(self.getStringDataArrays(), "StringDataArray");
    check(self.getIntegerDataArrays(), "IntegerDataArray");
}

/// Drop the data arrays that a peak count of @p n stranded. Only PeakMetadataPolicy::Clear acts.
///
/// Call this *after* the peaks have been written, never before. The incoming intensity/position
/// arrays may alias the very storage this releases -- FloatDataArray.get_data_view() hands out a
/// zero-copy view into a data array, and getFloatDataArrays() exposes the container's own arrays
/// by reference -- so releasing it up front would leave the fill loop reading freed memory.
/// Deferring also means a throwing resize() cannot leave the arrays already dropped.
///
/// Arrays that still fit, and empty ones, are kept: they remain valid annotations, and the
/// docstring promises as much. Dropping rather than resizing to @p n is deliberate, since padding
/// would fabricate annotations that nothing marks as made up.
template <typename Container>
void dropStrandedPeakMetadata(Container& self, size_t n, PeakMetadataPolicy policy)
{
    if (policy != PeakMetadataPolicy::Clear) return;

    auto drop = [n](auto& arrays) {
        arrays.erase(std::remove_if(arrays.begin(), arrays.end(),
                                    [n](const auto& array) { return !peakMetadataArrayIsAligned(array, n); }),
                     arrays.end());
    };
    drop(self.getFloatDataArrays());
    drop(self.getStringDataArrays());
    drop(self.getIntegerDataArrays());
}

/// Replace the peaks of @p self with @p n new ones, honouring @p policy for the data arrays.
///
/// This exists so the ordering below is written down once instead of once per set_peaks()
/// overload. Every step is load-bearing and none of them commute:
///
///  - the alignment check runs against the *incoming* count and before resize(), so a rejected
///    call is a strict no-op rather than leaving the peaks already replaced;
///  - @p fill runs after resize() and is the only container-specific part (m/z, RT or mobility);
///  - the stranded-array drop runs after @p fill, because the caller's arrays may alias the very
///    storage it releases -- see dropStrandedPeakMetadata().
///
/// Callers are responsible for validating their own inputs first, and for anything that has to
/// happen after the peaks exist (installing a replacement ion mobility array, for instance).
template <typename Container, typename Fill>
void writePeaksWithPolicy(Container& self, size_t n, PeakMetadataPolicy policy, const char* container,
                          Fill fill, std::optional<size_t> replaced_float_index = std::nullopt)
{
    checkPeakMetadataAlignment(self, n, policy, container, replaced_float_index);
    self.resize(n);
    fill(self, n);
    dropStrandedPeakMetadata(self, n, policy);
}
