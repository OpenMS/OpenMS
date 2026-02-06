// pyOpenMS nanobind bindings
// Domain: processing

#include "all_casters.h"
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/PROCESSING/CALIBRATION/InternalCalibration.h>
#include <OpenMS/PROCESSING/CALIBRATION/MZTrafoModel.h>
#include <OpenMS/PROCESSING/CALIBRATION/PrecursorCorrection.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/PROCESSING/FEATURE/FeatureOverlapFilter.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/PROCESSING/MISC/DataFilters.h>
#include <OpenMS/PROCESSING/MISC/SplineInterpolatedPeaks.h>
#include <OpenMS/PROCESSING/MISC/SplinePackage.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>
#include <limits>
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

NB_MODULE(_pyopenms_processing, m) {
    m.doc() = "pyOpenMS processing bindings";

    // -----------------------------------------------------------------------
    // DataFilters
    // -----------------------------------------------------------------------
    auto datafilters_class = nb::class_<OpenMS::DataFilters>(m, "DataFilters", "DataFilter array providing some convenience functions")
        .def(nb::init<>())
        .def("size", [](const OpenMS::DataFilters& self) { return self.size(); })
        .def("__getitem__", [](OpenMS::DataFilters& self, size_t i) { if (i >= self.size()) throw nb::index_error(); return self[i]; })
        .def("add", [](OpenMS::DataFilters& self, const OpenMS::DataFilters::DataFilter& filter) { return self.add(filter); }, "filter"_a)
        .def("remove", [](OpenMS::DataFilters& self, unsigned long index) { return self.remove(index); }, "index"_a)
        .def("replace", [](OpenMS::DataFilters& self, unsigned long index, const OpenMS::DataFilters::DataFilter& filter) { return self.replace(index, filter); }, "index"_a, "filter"_a)
        .def("clear", [](OpenMS::DataFilters& self) { return self.clear(); })
        .def("setActive", [](OpenMS::DataFilters& self, bool is_active) { return self.setActive(is_active); }, "is_active"_a)
        .def("isActive", [](const OpenMS::DataFilters& self) { return self.isActive(); })
        .def("__len__", [](OpenMS::DataFilters& self) { return self.size(); })
        .def("__getitem__", [](OpenMS::DataFilters& self, size_t i) -> const OpenMS::DataFilters::DataFilter & { 
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)
        ;
    // FilterType enum nested under DataFilters
    nb::enum_<OpenMS::DataFilters::FilterType>(datafilters_class, "FilterType", nb::is_arithmetic())
        .value("INTENSITY", OpenMS::DataFilters::FilterType::INTENSITY)
        .value("QUALITY", OpenMS::DataFilters::FilterType::QUALITY)
        .value("CHARGE", OpenMS::DataFilters::FilterType::CHARGE)
        .value("SIZE", OpenMS::DataFilters::FilterType::SIZE)
        .value("META_DATA", OpenMS::DataFilters::FilterType::META_DATA)
        .export_values();
    // FilterOperation enum nested under DataFilters
    nb::enum_<OpenMS::DataFilters::FilterOperation>(datafilters_class, "FilterOperation", nb::is_arithmetic())
        .value("GREATER_EQUAL", OpenMS::DataFilters::FilterOperation::GREATER_EQUAL)
        .value("EQUAL", OpenMS::DataFilters::FilterOperation::EQUAL)
        .value("LESS_EQUAL", OpenMS::DataFilters::FilterOperation::LESS_EQUAL)
        .value("EXISTS", OpenMS::DataFilters::FilterOperation::EXISTS)
        .export_values();

    // -----------------------------------------------------------------------
    // Deisotoper
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Deisotoper>(m, "Deisotoper", "OpenMS class Deisotoper")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Deisotoper &>())

        .def_static("deisotopeAndSingleCharge", [](OpenMS::MSSpectrum& spectrum,
             double fragment_tolerance, bool fragment_unit_ppm,
             int min_charge, int max_charge, bool keep_only_deisotoped,
             unsigned int min_isopeaks, unsigned int max_isopeaks,
             bool make_single_charged, bool annotate_charge,
             bool annotate_iso_peak_count, bool use_decreasing_model,
             unsigned int start_intensity_check, bool add_up_intensity,
             bool annotate_features) {
            OpenMS::Deisotoper::deisotopeAndSingleCharge(spectrum,
                fragment_tolerance, fragment_unit_ppm, min_charge, max_charge,
                keep_only_deisotoped, min_isopeaks, max_isopeaks,
                make_single_charged, annotate_charge, annotate_iso_peak_count,
                use_decreasing_model, start_intensity_check, add_up_intensity,
                annotate_features);
        }, "spectrum"_a, "fragment_tolerance"_a, "fragment_unit_ppm"_a,
           "min_charge"_a = 1, "max_charge"_a = 3, "keep_only_deisotoped"_a = false,
           "min_isopeaks"_a = 3, "max_isopeaks"_a = 10,
           "make_single_charged"_a = true, "annotate_charge"_a = false,
           "annotate_iso_peak_count"_a = false, "use_decreasing_model"_a = true,
           "start_intensity_check"_a = 2, "add_up_intensity"_a = false,
           "annotate_features"_a = false,
           "Deisotope and single charge a spectrum")
        ;

    // -----------------------------------------------------------------------
    // FeatureOverlapFilter
    // -----------------------------------------------------------------------
    auto featureoverlapfilter_class = nb::class_<OpenMS::FeatureOverlapFilter>(m, "FeatureOverlapFilter", 
        R"doc(
Filters and merges overlapping features in a FeatureMap.
This class provides functionality for detecting and handling overlapping features
based on various criteria including convex hull overlap, trace-level overlap,
and centroid-based distance thresholds.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FeatureOverlapFilter &>())
        ;
    // MergeIntensityMode enum nested under FeatureOverlapFilter
    nb::enum_<OpenMS::MergeIntensityMode>(featureoverlapfilter_class, "MergeIntensityMode", nb::is_arithmetic())
        .value("SUM", OpenMS::MergeIntensityMode::SUM)
        .value("MAX", OpenMS::MergeIntensityMode::MAX)
        ;

    // -----------------------------------------------------------------------
    // IDFilter
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IDFilter>(m, "IDFilter", 
        R"doc(
Finds the best-scoring hit in a vector of peptide or protein identifications\n
This class provides functions for filtering collections of peptide or protein identifications according to various criteria.
It also contains helper functions and classes (functors that implement predicates) that are used in this context.\n
The filter functions modify their inputs, rather than creating filtered copies.\n
Most filters work on the hit level, i.e. they remove peptide or protein hits from peptide or protein identifications (IDs).
A few filters work on the ID level instead, i.e. they remove peptide or protein IDs from vectors thereof.
Independent of this, the inputs for all filter functions are vectors of IDs, because the data most often comes in this form.
This design also allows many helper objects to be set up only once per vector, rather than once per ID.\n
The filter functions for vectors of peptide/protein IDs do not include clean-up steps (e.g. removal of IDs without hits, reassignment of hit ranks, ...).
They only carry out their specific filtering operations.
This is so filters can be chained without having to repeat clean-up operations.
The group of clean-up functions provides helpers that are useful to ensure data integrity after filters have been applied, but it is up to the individual developer to use them when necessary.\n
The filter functions for MS/MS experiments do include clean-up steps, because they filter peptide and protein IDs in conjunction and potential contradictions between the two must be eliminated.
)doc")
        .def(nb::init<>())
        .def_static("filterHitsByRank", [](OpenMS::PeptideIdentificationList& ids, unsigned long min_rank, unsigned long max_rank) { return OpenMS::IDFilter::filterHitsByRank(ids, min_rank, max_rank); }, "ids"_a, "min_rank"_a, "max_rank"_a)
        .def_static("removeHitsMatchingProteins", [](OpenMS::PeptideIdentificationList& ids, const std::set<OpenMS::String>& accessions) { return OpenMS::IDFilter::removeHitsMatchingProteins(ids, accessions); }, "ids"_a, "accessions"_a, "Filters peptide or protein identifications according to the given proteins (negative)")
        .def_static("keepHitsMatchingProteins", [](OpenMS::PeptideIdentificationList& ids, const std::set<OpenMS::String>& accessions) { return OpenMS::IDFilter::keepHitsMatchingProteins(ids, accessions); }, "ids"_a, "accessions"_a, "Filters peptide or protein identifications according to the given proteins (positive)")
        .def_static("getBestHit", [](OpenMS::PeptideIdentificationList& ids, bool assume_sorted, OpenMS::PeptideHit& best_hit) { return OpenMS::IDFilter::getBestHit(ids, assume_sorted, best_hit); }, "ids"_a, "assume_sorted"_a, "best_hit"_a)
        .def_static("removeEmptyIdentifications", [](OpenMS::PeptideIdentificationList& ids) { return OpenMS::IDFilter::removeEmptyIdentifications(ids); }, "ids"_a, "Removes peptide or protein identifications that have no hits in them")
        .def_static("extractPeptideSequences", [](const OpenMS::PeptideIdentificationList& peptides, bool ignore_mods) { std::set<OpenMS::String> sequences; OpenMS::IDFilter::extractPeptideSequences(peptides, sequences, ignore_mods); return sequences; }, "peptides"_a, "ignore_mods"_a, 
            R"doc(
Finds the best-scoring hit in a vector of peptide or protein identifications
If there are several hits with the best score, the first one is taken
:param identifications: Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
:param assume_sorted: Are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
:param best_hit: Contains the best hit if successful in a vector of protein identifications
:return: true if a hit was present, false otherwise
)doc")
        .def_static("removeUnreferencedProteins", [](OpenMS::ConsensusMap& cmap, bool include_unassigned) { return OpenMS::IDFilter::removeUnreferencedProteins(cmap, include_unassigned); }, "cmap"_a, "include_unassigned"_a, "Removes protein hits from the protein IDs in a 'cmap' that are not referenced by a peptide in the features or if requested in the unassigned peptide list")
        .def_static("removeUnreferencedProteins", [](std::vector<OpenMS::ProteinIdentification> proteins, const OpenMS::PeptideIdentificationList& peptides) { OpenMS::IDFilter::removeUnreferencedProteins(proteins, peptides); return proteins; }, "proteins"_a, "peptides"_a, "Removes protein hits from the protein IDs in a 'cmap' that are not referenced by a peptide in the features or if requested in the unassigned peptide list")
        .def_static("removeUnreferencedProteins", [](OpenMS::ProteinIdentification& proteins, const OpenMS::PeptideIdentificationList& peptides) { return OpenMS::IDFilter::removeUnreferencedProteins(proteins, peptides); }, "proteins"_a, "peptides"_a, "Removes protein hits from the protein IDs in a 'cmap' that are not referenced by a peptide in the features or if requested in the unassigned peptide list")
        .def_static("removeDanglingProteinReferences", [](OpenMS::PeptideIdentificationList& peptides, const std::vector<OpenMS::ProteinIdentification>& proteins, bool remove_peptides_without_reference) { return OpenMS::IDFilter::removeDanglingProteinReferences(peptides, proteins, remove_peptides_without_reference); }, "peptides"_a, "proteins"_a, "remove_peptides_without_reference"_a)
        .def_static("removeDanglingProteinReferences", [](OpenMS::ConsensusMap& cmap, bool remove_peptides_without_reference) { return OpenMS::IDFilter::removeDanglingProteinReferences(cmap, remove_peptides_without_reference); }, "cmap"_a, "remove_peptides_without_reference"_a)
        .def_static("removeDanglingProteinReferences", [](OpenMS::ConsensusMap& cmap, const OpenMS::ProteinIdentification& ref_run, bool remove_peptides_without_reference) { return OpenMS::IDFilter::removeDanglingProteinReferences(cmap, ref_run, remove_peptides_without_reference); }, "cmap"_a, "ref_run"_a, "remove_peptides_without_reference"_a)
        .def_static("updateProteinGroups", [](std::vector<OpenMS::ProteinIdentification::ProteinGroup> groups, const std::vector<OpenMS::ProteinHit>& hits) { auto result = OpenMS::IDFilter::updateProteinGroups(groups, hits); return nb::make_tuple(result, groups); }, "groups"_a, "hits"_a,
            R"doc(
Removes dangling protein references from peptide hits
Cleans up PeptideEvidence entries by removing references to proteins that no longer exist
in the provided protein identifications. This is typically called after filtering protein hits
to maintain consistency between peptide-to-protein mappings.
:param peptides: The peptide identifications to process (in/out)
:param proteins: The protein identifications containing valid protein hits (in)
:param remove_peptides_without_reference: If true, peptide hits that have no remaining
protein references after cleanup are also removed (default: false) (in)
)doc")
        .def_static("keepNBestHits", [](OpenMS::PeptideIdentificationList& pep_ids, unsigned long n) { return OpenMS::IDFilter::keepNBestHits(pep_ids, n); }, "pep_ids"_a, "n"_a)
        .def_static("keepBestPeptideHits", [](OpenMS::PeptideIdentificationList& peptides, bool strict) { return OpenMS::IDFilter::keepBestPeptideHits(peptides, strict); }, "peptides"_a, "strict"_a)
        .def_static("filterPeptidesByLength", [](OpenMS::PeptideIdentificationList& peptides, unsigned long min_length, unsigned long max_length) { return OpenMS::IDFilter::filterPeptidesByLength(peptides, min_length, max_length); }, "peptides"_a, "min_length"_a, "max_length"_a, "Filters peptide identifications according to peptide sequence length")
        .def_static("filterPeptidesByCharge", [](OpenMS::PeptideIdentificationList& peptides, int min_charge, int max_charge) { return OpenMS::IDFilter::filterPeptidesByCharge(peptides, min_charge, max_charge); }, "peptides"_a, "min_charge"_a, "max_charge"_a, "Filters peptide identifications according to charge state")
        .def_static("filterPeptidesByRT", [](OpenMS::PeptideIdentificationList& peptides, double min_rt, double max_rt) { return OpenMS::IDFilter::filterPeptidesByRT(peptides, min_rt, max_rt); }, "peptides"_a, "min_rt"_a, "max_rt"_a, "Filters peptide identifications by precursor RT, keeping only IDs in the given range")
        .def_static("filterPeptidesByMZ", [](OpenMS::PeptideIdentificationList& peptides, double min_mz, double max_mz) { return OpenMS::IDFilter::filterPeptidesByMZ(peptides, min_mz, max_mz); }, "peptides"_a, "min_mz"_a, "max_mz"_a, "Filters peptide identifications by precursor m/z, keeping only IDs in the given range")
        .def_static("filterPeptidesByMZError", [](OpenMS::PeptideIdentificationList& peptides, double mass_error, bool unit_ppm) { return OpenMS::IDFilter::filterPeptidesByMZError(peptides, mass_error, unit_ppm); }, "peptides"_a, "mass_error"_a, "unit_ppm"_a, "Filter peptide identifications according to mass deviation")
        .def_static("filterPeptidesByRTPredictPValue", [](OpenMS::PeptideIdentificationList& peptides, const OpenMS::String& metavalue_key, double threshold) { return OpenMS::IDFilter::filterPeptidesByRTPredictPValue(peptides, metavalue_key, threshold); }, "peptides"_a, "metavalue_key"_a, "threshold"_a)
        .def_static("removePeptidesWithMatchingModifications", [](OpenMS::PeptideIdentificationList& peptides, const std::set<OpenMS::String>& modifications) { return OpenMS::IDFilter::removePeptidesWithMatchingModifications(peptides, modifications); }, "peptides"_a, "modifications"_a, "Removes all peptide hits that have at least one of the given modifications")
        .def_static("keepPeptidesWithMatchingModifications", [](OpenMS::PeptideIdentificationList& peptides, const std::set<OpenMS::String>& modifications) { return OpenMS::IDFilter::keepPeptidesWithMatchingModifications(peptides, modifications); }, "peptides"_a, "modifications"_a, "Keeps only peptide hits that have at least one of the given modifications")
        .def_static("removePeptidesWithMatchingSequences", [](OpenMS::PeptideIdentificationList& peptides, const OpenMS::PeptideIdentificationList& bad_peptides, bool ignore_mods) { return OpenMS::IDFilter::removePeptidesWithMatchingSequences(peptides, bad_peptides, ignore_mods); }, "peptides"_a, "bad_peptides"_a, "ignore_mods"_a, "Removes all peptide hits with a sequence that matches one in 'bad_peptides'")
        .def_static("keepPeptidesWithMatchingSequences", [](OpenMS::PeptideIdentificationList& peptides, const OpenMS::PeptideIdentificationList& good_peptides, bool ignore_mods) { return OpenMS::IDFilter::keepPeptidesWithMatchingSequences(peptides, good_peptides, ignore_mods); }, "peptides"_a, "good_peptides"_a, "ignore_mods"_a, "Removes all peptide hits with a sequence that does not match one in 'good_peptides'")
        .def_static("keepUniquePeptidesPerProtein", [](OpenMS::PeptideIdentificationList& peptides) { return OpenMS::IDFilter::keepUniquePeptidesPerProtein(peptides); }, "peptides"_a, "Removes all peptides that are not annotated as unique for a protein (by PeptideIndexer)")
        .def_static("removeDuplicatePeptideHits", [](OpenMS::PeptideIdentificationList& peptides, bool seq_only) { return OpenMS::IDFilter::removeDuplicatePeptideHits(peptides, seq_only); }, "peptides"_a, "seq_only"_a, "Removes duplicate peptide hits from each peptide identification, keeping only unique hits (per ID)")
        .def_static("filterHitsByScore", [](OpenMS::AnnotatedMSRun& annotated_data, double peptide_threshold_score, double protein_threshold_score) { return OpenMS::IDFilter::filterHitsByScore(annotated_data, peptide_threshold_score, protein_threshold_score); }, "annotated_data"_a, "peptide_threshold_score"_a, "protein_threshold_score"_a, "Filters an MS/MS experiment according to score thresholds")
        .def_static("keepNBestHits", [](OpenMS::AnnotatedMSRun& annotated_data, unsigned long n) { return OpenMS::IDFilter::keepNBestHits(annotated_data, n); }, "annotated_data"_a, "n"_a)
        .def_static("keepNBestSpectra", [](OpenMS::PeptideIdentificationList& peptides, unsigned long n) { return OpenMS::IDFilter::keepNBestSpectra(peptides, n); }, "peptides"_a, "n"_a, 
            R"doc(
Filter identifications by "N best" PeptideIdentification objects (better PeptideIdentification means better [best] PeptideHit than other)
)doc")
        .def_static("keepBestPerPeptide", [](OpenMS::PeptideIdentificationList& pep_ids, bool ignore_mods, bool ignore_charges, unsigned long nr_best_spectrum) { return OpenMS::IDFilter::keepBestPerPeptide(pep_ids, ignore_mods, ignore_charges, nr_best_spectrum); }, "pep_ids"_a, "ignore_mods"_a, "ignore_charges"_a, "nr_best_spectrum"_a, "Filters PeptideHits from PeptideIdentification by keeping only the best peptide hits for every peptide sequence")
        .def_static("keepBestPerPeptidePerRun", [](std::vector<OpenMS::ProteinIdentification> prot_ids, OpenMS::PeptideIdentificationList& pep_ids, bool ignore_mods, bool ignore_charges, unsigned long nr_best_spectrum) { OpenMS::IDFilter::keepBestPerPeptidePerRun(prot_ids, pep_ids, ignore_mods, ignore_charges, nr_best_spectrum); return prot_ids; }, "prot_ids"_a, "pep_ids"_a, "ignore_mods"_a, "ignore_charges"_a, "nr_best_spectrum"_a, "Filters PeptideHits from PeptideIdentification by keeping only the best peptide hits for every peptide sequence on a per run basis")
        .def_static("removeDecoyHits", [](OpenMS::PeptideIdentificationList& ids) { return OpenMS::IDFilter::removeDecoyHits(ids); }, "ids"_a, 
            R"doc(
Removes hits annotated as decoys from peptide or protein identifications. Checks for meta values named "target_decoy" and "isDecoy", and removes protein/peptide hits if the values are "decoy" and "true", respectively
)doc")
        .def_static("filterHitsByScore", [](OpenMS::PeptideIdentificationList& ids, double threshold_score) { return OpenMS::IDFilter::filterHitsByScore(ids, threshold_score); }, "ids"_a, "threshold_score"_a, "Filters peptide or protein identifications according to the score of the hits. The score orientation has to be set to higherscorebetter in each PeptideIdentification. Only peptide/protein hits with a score at least as good as 'threshold_score' are kept")
        .def_static("removeUnreferencedProteins", [](std::vector<OpenMS::ProteinIdentification> proteins, OpenMS::PeptideIdentificationList& ids) { OpenMS::IDFilter::removeUnreferencedProteins(proteins, ids); return proteins; }, "proteins"_a, "ids"_a, "Removes protein hits from the protein IDs in a 'cmap' that are not referenced by a peptide in the features or if requested in the unassigned peptide list")
        ;

    // -----------------------------------------------------------------------
    // InternalCalibration_LockMass
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::InternalCalibration::LockMass>(m, "InternalCalibration_LockMass", "OpenMS class InternalCalibration_LockMass")
        .def(nb::init<double, int, int>())
        .def_rw("mz", &OpenMS::InternalCalibration::LockMass::mz)
        .def_rw("ms_level", &OpenMS::InternalCalibration::LockMass::ms_level)
        .def_rw("charge", &OpenMS::InternalCalibration::LockMass::charge)
        ;

    // -----------------------------------------------------------------------
    // MZTrafoModel
    // -----------------------------------------------------------------------
    auto mztrafomodel_class = nb::class_<OpenMS::MZTrafoModel>(m, "MZTrafoModel", 
        R"doc(
Create and apply models of a mass recalibration function
The input is a list of calibration points (ideally spanning a wide m/z range to prevent extrapolation when applying to model)
Models (LINEAR, LINEAR_WEIGHTED, QUADRATIC, QUADRATIC_WEIGHTED) can be trained using CalData points (or a subset of them)
Calibration points can have different retention time points, and a model should be build such that it captures
the local (in time) decalibration of the instrument, i.e. choose appropriate time windows along RT to calibrate the
spectra in this RT region
From the available calibrant data, a model is build. Later, any uncalibrated m/z value can be fed to the model, to obtain
a calibrated m/z
The input domain can either be absolute mass differences in [Th], or relative differences in [ppm]
The models are build based on this input
Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models
)doc")
        .def(nb::init<>())
        .def(nb::init<bool>())
        .def_static("nameToEnum", [](const OpenMS::String& name) { return OpenMS::MZTrafoModel::nameToEnum(name); }, "name"_a)
        .def_static("enumToName", [](OpenMS::MZTrafoModel::MODELTYPE mt) { return OpenMS::MZTrafoModel::enumToName(mt); }, "mt"_a)
        .def_static("setRANSACParams", [](const OpenMS::Math::RANSACParam& p) { return OpenMS::MZTrafoModel::setRANSACParams(p); }, "p"_a)
        .def_static("setCoefficientLimits", [](double offset, double scale, double power) { return OpenMS::MZTrafoModel::setCoefficientLimits(offset, scale, power); }, "offset"_a, "scale"_a, "power"_a)
        .def_static("isValidModel", [](const OpenMS::MZTrafoModel& trafo) { return OpenMS::MZTrafoModel::isValidModel(trafo); }, "trafo"_a)
        .def("isTrained", [](const OpenMS::MZTrafoModel& self) { return self.isTrained(); }, "Returns true if the model have coefficients (i.e. was trained successfully)")
        .def("getRT", [](const OpenMS::MZTrafoModel& self) { return self.getRT(); }, "Get RT associated with the model (training region)")
        .def("predict", [](const OpenMS::MZTrafoModel& self, double mz) { return self.predict(mz); }, "mz"_a)
        .def("train", [](OpenMS::MZTrafoModel& self, const OpenMS::CalibrationData& cd, OpenMS::MZTrafoModel::MODELTYPE md, bool use_RANSAC, double rt_left, double rt_right) { return self.train(cd, md, use_RANSAC, rt_left, rt_right); }, "cd"_a, "md"_a, "use_RANSAC"_a, "rt_left"_a = -std::numeric_limits<double>::max(), "rt_right"_a = std::numeric_limits<double>::max(), 
            R"doc(
Apply the model to an uncalibrated m/z value
Make sure the model was trained (train()) and is valid (isValidModel()) before calling this function!
Applies the function y = intercept + slope*mz + power*mz^2
and returns y
:param mz: The uncalibrated m/z value
:return: The calibrated m/z value
)doc")
        .def("train", [](OpenMS::MZTrafoModel& self, std::vector<double> error_mz, std::vector<double> theo_mz, std::vector<double> weights, OpenMS::MZTrafoModel::MODELTYPE md, bool use_RANSAC) { return self.train(error_mz, theo_mz, weights, md, use_RANSAC); }, "error_mz"_a, "theo_mz"_a, "weights"_a, "md"_a, "use_RANSAC"_a, 
            R"doc(
Apply the model to an uncalibrated m/z value
Make sure the model was trained (train()) and is valid (isValidModel()) before calling this function!
Applies the function y = intercept + slope*mz + power*mz^2
and returns y
:param mz: The uncalibrated m/z value
:return: The calibrated m/z value
)doc")
        .def("getCoefficients", [](OpenMS::MZTrafoModel& self, double& intercept, double& slope, double& power) { return self.getCoefficients(intercept, slope, power); }, "intercept"_a, "slope"_a, "power"_a, 
            R"doc(
Train a model using calibrant data
Given theoretical and observed mass values (and corresponding weights),
a model (linear, quadratic, ...) is build
Outlier removal is applied before
The 'obs_mz' can be either given as absolute masses in [Th] or relative deviations in [ppm]
The MZTrafoModel must be constructed accordingly (see constructor). This has no influence on the model building itself, but
rather on how 'predict()' works internally
Outlier detection before model building via the RANSAC algorithm is supported for LINEAR and QUADRATIC models
Internally, these steps take place:
- [apply RANSAC] (depending on 'use_RANSAC')
- build model and store its parameters internally
:param error_mz: Observed Mass error (in ppm or Th)
:param theo_mz: Theoretical m/z values, corresponding to 'error_mz'
:param weights: For weighted models only: weight of calibrants; ignored otherwise
:param md: Type of model (linear, quadratic, ...)
:param use_RANSAC: Remove outliers before computing the model?
:return: True if model was build, false otherwise
)doc")
        .def("setCoefficients", [](OpenMS::MZTrafoModel& self, const OpenMS::MZTrafoModel& rhs) { return self.setCoefficients(rhs); }, "rhs"_a, "Copy model coefficients from another model")
        .def("setCoefficients", [](OpenMS::MZTrafoModel& self, double intercept, double slope, double power) { return self.setCoefficients(intercept, slope, power); }, "intercept"_a, "slope"_a, "power"_a)
        .def("toString", [](const OpenMS::MZTrafoModel& self) { return self.toString(); }, 
            R"doc(
Manually set model coefficients
Can be used instead of train(), so manually set coefficients
It must be exactly three values. If you want a linear model, set 'power' to zero
If you want a constant model, set slope to zero in addition
:param intercept: The offset
:param slope: The slope
:param power: The x*x coefficient (for quadratic models)
)doc")
        ;
    // MODELTYPE enum nested under MZTrafoModel
    nb::enum_<OpenMS::MZTrafoModel::MODELTYPE>(mztrafomodel_class, "MODELTYPE", nb::is_arithmetic())
        .value("LINEAR", OpenMS::MZTrafoModel::MODELTYPE::LINEAR)
        .value("LINEAR_WEIGHTED", OpenMS::MZTrafoModel::MODELTYPE::LINEAR_WEIGHTED)
        .value("QUADRATIC", OpenMS::MZTrafoModel::MODELTYPE::QUADRATIC)
        .value("QUADRATIC_WEIGHTED", OpenMS::MZTrafoModel::MODELTYPE::QUADRATIC_WEIGHTED)
        .value("SIZE_OF_MODELTYPE", OpenMS::MZTrafoModel::MODELTYPE::SIZE_OF_MODELTYPE)
        ;

    // -----------------------------------------------------------------------
    // PrecursorCorrection
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PrecursorCorrection>(m, "PrecursorCorrection", "This class provides methods for precursor correction")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PrecursorCorrection &>())
        .def_static("getPrecursors", [](const OpenMS::MSExperiment& exp) { std::vector<OpenMS::Precursor> precursors; std::vector<double> precursors_rt; std::vector<unsigned long> precursor_scan_index; OpenMS::PrecursorCorrection::getPrecursors(exp, precursors, precursors_rt, precursor_scan_index); return nb::make_tuple(precursors, precursors_rt, precursor_scan_index); }, "exp"_a)
        .def_static("writeHist", [](const OpenMS::String& out_csv, const std::vector<double>& delta_mzs, const std::vector<double>& mzs, const std::vector<double>& rts) { return OpenMS::PrecursorCorrection::writeHist(out_csv, delta_mzs, mzs, rts); }, "out_csv"_a, "delta_mzs"_a, "mzs"_a, "rts"_a)
        .def_static("correctToNearestMS1Peak", [](OpenMS::MSExperiment& exp, double mz_tolerance, bool ppm) { std::vector<double> delta_mzs, mzs, rts; auto result = OpenMS::PrecursorCorrection::correctToNearestMS1Peak(exp, mz_tolerance, ppm, delta_mzs, mzs, rts); return nb::make_tuple(result, delta_mzs, mzs, rts); }, "exp"_a, "mz_tolerance"_a, "ppm"_a)
        .def_static("correctToHighestIntensityMS1Peak", [](OpenMS::MSExperiment& exp, double mz_tolerance, bool ppm) { std::vector<double> delta_mzs, mzs, rts; auto result = OpenMS::PrecursorCorrection::correctToHighestIntensityMS1Peak(exp, mz_tolerance, ppm, delta_mzs, mzs, rts); return nb::make_tuple(result, delta_mzs, mzs, rts); }, "exp"_a, "mz_tolerance"_a, "ppm"_a)
        ;

    // -----------------------------------------------------------------------
    // SignalToNoiseEstimatorMeanIterative
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SignalToNoiseEstimatorMeanIterative<OpenMS::MSSpectrum>>(m, "SignalToNoiseEstimatorMeanIterative", "OpenMS class SignalToNoiseEstimatorMeanIterative")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SignalToNoiseEstimatorMeanIterative<OpenMS::MSSpectrum>&>())
        ;


    // -----------------------------------------------------------------------
    // SignalToNoiseEstimatorMedian
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SignalToNoiseEstimatorMedian<OpenMS::MSSpectrum>>(m, "SignalToNoiseEstimatorMedian", "OpenMS class SignalToNoiseEstimatorMedian")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SignalToNoiseEstimatorMedian<OpenMS::MSSpectrum>&>())
        ;


    // -----------------------------------------------------------------------
    // SignalToNoiseEstimatorMedianRapid
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SignalToNoiseEstimatorMedianRapid>(m, "SignalToNoiseEstimatorMedianRapid", 
        R"doc(
Estimates the signal/noise (S/N) ratio of each data point in a scan
by using the median (window based)
)doc")
        .def(nb::init<double>())
        .def("estimateNoise", [](OpenMS::SignalToNoiseEstimatorMedianRapid& self, std::shared_ptr<OpenMS::Interfaces::Spectrum> spectrum) { return self.estimateNoise(spectrum); }, "spectrum"_a)
        .def("estimateNoise", [](OpenMS::SignalToNoiseEstimatorMedianRapid& self, std::shared_ptr<OpenMS::Interfaces::Chromatogram> chrom) { return self.estimateNoise(chrom); }, "chrom"_a)
        .def("estimateNoise", [](OpenMS::SignalToNoiseEstimatorMedianRapid& self, const std::vector<double>& mz_array, const std::vector<double>& int_array) { return self.estimateNoise(mz_array, int_array); }, "mz_array"_a, "int_array"_a)
        ;

    // -----------------------------------------------------------------------
    // SplineInterpolatedPeaks
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SplineInterpolatedPeaks>(m, "SplineInterpolatedPeaks", 
        R"doc(
Data structure for spline interpolation of MS1 spectra and
chromatograms * * The data structure consists of a set of splines,
each interpolating the MS1 spectrum (or chromatogram) in a * certain
m/z (or RT) range. Between these splines no raw data points exist and
the intensity is identical to zero. * * A spline on non-equi-distant
input data is not well supported in regions without data points.
Hence, a spline tends to * swing wildly in these regions and cannot be
used for reliable interpolation. We assume that in m/z (or RT) regions
* without data points, the spectrum (or chromatogram) is identical to
zero. * * @see SplinePackage * @see MSSpectrum * @see MSChromatogram
)doc")
        .def(nb::init<std::vector<double>, std::vector<double>>())
        .def(nb::init<OpenMS::MSSpectrum>())
        .def(nb::init<OpenMS::MSChromatogram>())
        .def("getPosMin", [](const OpenMS::SplineInterpolatedPeaks& self) { return self.getPosMin(); })
        .def("getPosMax", [](const OpenMS::SplineInterpolatedPeaks& self) { return self.getPosMax(); })
        .def("size", [](const OpenMS::SplineInterpolatedPeaks& self) { return self.size(); })
        .def("getNavigator", [](OpenMS::SplineInterpolatedPeaks& self, double scaling) { return self.getNavigator(scaling); }, "scaling"_a = 0.7)
        .def("__len__", [](OpenMS::SplineInterpolatedPeaks& self) { return self.size(); })
        ;

    // -----------------------------------------------------------------------
    // SplinePackage
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SplinePackage>(m, "SplinePackage", 
        R"doc(
fundamental data structure for SplineInterpolatedPeaks * * In many
cases, data points in MS spectra (or chromatograms) are not
equidistant in m/z (or RT) * but consist of packages of data points
separated by wide m/z (or RT) ranges with zero intensity. *
SplinePackage contains the spline fit of a single set of such data
points. * * @see SplineInterpolatedPeaks
)doc")
        .def(nb::init<std::vector<double>, std::vector<double>>())
        .def("getPosMin", [](const OpenMS::SplinePackage& self) { return self.getPosMin(); }, "Returns the minimum position for which the spline fit is valid")
        .def("getPosMax", [](const OpenMS::SplinePackage& self) { return self.getPosMax(); }, "Returns the maximum position for which the spline fit is valid")
        .def("getPosStepWidth", [](const OpenMS::SplinePackage& self) { return self.getPosStepWidth(); }, "Returns a sensible position step width for the package")
        .def("isInPackage", [](const OpenMS::SplinePackage& self, double pos) { return self.isInPackage(pos); }, "pos"_a, "Returns true if position in [posMin:posMax] interval else false")
        .def("eval", [](const OpenMS::SplinePackage& self, double pos) { return self.eval(pos); }, "pos"_a, "Returns interpolated intensity position `pos`")
        ;

}