// pyOpenMS nanobind bindings
// Domain: qc

#include "all_casters.h"
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/QC/Contaminants.h>
#include <OpenMS/QC/DBSuitability.h>
#include <OpenMS/QC/FWHM.h>
#include <OpenMS/QC/FeatureSummary.h>
#include <OpenMS/QC/FragmentMassError.h>
#include <OpenMS/QC/IdentificationSummary.h>
#include <OpenMS/QC/MissedCleavages.h>
#include <OpenMS/QC/Ms2IdentificationRate.h>
#include <OpenMS/QC/Ms2SpectrumStats.h>
#include <OpenMS/QC/MzCalibration.h>
#include <OpenMS/QC/PSMExplainedIonCurrent.h>
#include <OpenMS/QC/PeptideMass.h>
#include <OpenMS/QC/RTAlignment.h>
#include <OpenMS/QC/SpectrumCount.h>
#include <OpenMS/QC/TIC.h>

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/shared_ptr.h>
#include "binding_utils.h"

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_pyopenms_qc, m) {
    m.doc() = "pyOpenMS Quality Control (QC) bindings";

    // -------------------------------------------------------------------------
    // QCBase Enums & Internal Classes
    // -------------------------------------------------------------------------
    nb::enum_<OpenMS::QCBase::Requires>(m, "QCBaseRequires")
        .value("NOTHING", OpenMS::QCBase::Requires::NOTHING)
        .value("RAWMZML", OpenMS::QCBase::Requires::RAWMZML)
        .value("POSTFDRFEAT", OpenMS::QCBase::Requires::POSTFDRFEAT)
        .value("PREFDRFEAT", OpenMS::QCBase::Requires::PREFDRFEAT)
        .value("CONTAMINANTS", OpenMS::QCBase::Requires::CONTAMINANTS)
        .value("TRAFOALIGN", OpenMS::QCBase::Requires::TRAFOALIGN)
        .value("ID", OpenMS::QCBase::Requires::ID)
        ;

    nb::enum_<OpenMS::QCBase::ToleranceUnit>(m, "QCBaseToleranceUnit")
        .value("AUTO", OpenMS::QCBase::ToleranceUnit::AUTO)
        .value("PPM", OpenMS::QCBase::ToleranceUnit::PPM)
        .value("DA", OpenMS::QCBase::ToleranceUnit::DA)
        ;

    // FlagSet<Requires> / Status
    nb::class_<OpenMS::FlagSet<OpenMS::QCBase::Requires>>(m, "QCBaseStatus")
        .def(nb::init<>())
        .def(nb::init<OpenMS::QCBase::Requires>())
        .def(nb::self == nb::self)
        .def(nb::self & OpenMS::QCBase::Requires())
        .def(nb::self & nb::self)
        .def(nb::self | OpenMS::QCBase::Requires())
        .def(nb::self | nb::self)
        .def(nb::self - OpenMS::QCBase::Requires())
        .def(nb::self - nb::self)
        .def("__ior__", [](OpenMS::FlagSet<OpenMS::QCBase::Requires>& self, const OpenMS::QCBase::Requires& rhs) { self |= rhs; return self; })
        .def("__ior__", [](OpenMS::FlagSet<OpenMS::QCBase::Requires>& self, const OpenMS::FlagSet<OpenMS::QCBase::Requires>& rhs) { self |= rhs; return self; })
        .def("__iand__", [](OpenMS::FlagSet<OpenMS::QCBase::Requires>& self, const OpenMS::QCBase::Requires& rhs) { self &= rhs; return self; })
        .def("__iand__", [](OpenMS::FlagSet<OpenMS::QCBase::Requires>& self, const OpenMS::FlagSet<OpenMS::QCBase::Requires>& rhs) { self &= rhs; return self; })
        .def("isSuperSetOf", [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& self, const OpenMS::FlagSet<OpenMS::QCBase::Requires>& required) { return self.isSuperSetOf(required); })
        .def("isSuperSetOf", [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& self, const OpenMS::QCBase::Requires& required) { return self.isSuperSetOf(required); })
        .def("empty", &OpenMS::FlagSet<OpenMS::QCBase::Requires>::empty)
        .def("value", &OpenMS::FlagSet<OpenMS::QCBase::Requires>::value)
        ;

    nb::class_<OpenMS::QCBase::SpectraMap>(m, "QCSpectraMap")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MSExperiment&>())
        .def("calculateMap", &OpenMS::QCBase::SpectraMap::calculateMap)
        .def("at", &OpenMS::QCBase::SpectraMap::at)
        .def("clear", &OpenMS::QCBase::SpectraMap::clear)
        .def("empty", &OpenMS::QCBase::SpectraMap::empty)
        .def("size", &OpenMS::QCBase::SpectraMap::size)
        ;

    // QCBase
    nb::class_<OpenMS::QCBase>(m, "QCBase")
        .def("getName", &OpenMS::QCBase::getName, nb::rv_policy::reference_internal)
        .def("requirements", &OpenMS::QCBase::requirements)
        .def("isRunnable", &OpenMS::QCBase::isRunnable)
        .def_static("isLabeledExperiment", &OpenMS::QCBase::isLabeledExperiment)
        ;

    // -------------------------------------------------------------------------
    // FWHM
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::FWHM, OpenMS::QCBase>(m, "FWHM")
        .def(nb::init<>())
        .def("compute", &OpenMS::FWHM::compute, "features"_a)
        ;

    // -------------------------------------------------------------------------
    // TIC
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::TIC::Result>(m, "TICResult")
        .def(nb::init<>())
        .def_rw("intensities", &OpenMS::TIC::Result::intensities)
        .def_rw("relative_intensities", &OpenMS::TIC::Result::relative_intensities)
        .def_rw("retention_times", &OpenMS::TIC::Result::retention_times)
        .def_rw("area", &OpenMS::TIC::Result::area)
        .def_rw("fall", &OpenMS::TIC::Result::fall)
        .def_rw("jump", &OpenMS::TIC::Result::jump)
        .def(nb::self == nb::self)
        ;

    nb::class_<OpenMS::TIC, OpenMS::QCBase>(m, "TIC")
        .def(nb::init<>())
        .def("compute", &OpenMS::TIC::compute, "exp"_a, "bin_size"_a = 0.0f, "ms_level"_a = 1)
        .def("getResults", &OpenMS::TIC::getResults, nb::rv_policy::reference_internal)
        .def("addMetaDataMetricsToMzTab", &OpenMS::TIC::addMetaDataMetricsToMzTab, "meta"_a, "tics"_a)
        ;

    // -------------------------------------------------------------------------
    // MissedCleavages
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::MissedCleavages, OpenMS::QCBase>(m, "MissedCleavages")
        .def(nb::init<>())
        .def("compute", nb::overload_cast<OpenMS::FeatureMap&>(&OpenMS::MissedCleavages::compute), "fmap"_a)
        .def("compute", nb::overload_cast<std::vector<OpenMS::ProteinIdentification>&, OpenMS::PeptideIdentificationList&>(&OpenMS::MissedCleavages::compute), "prot_ids"_a, "pep_ids"_a)
        .def("getResults", &OpenMS::MissedCleavages::getResults, nb::rv_policy::reference_internal)
        ;

    // -------------------------------------------------------------------------
    // FragmentMassError
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::FragmentMassError::Statistics>(m, "FragmentMassErrorStatistics")
        .def(nb::init<>())
        .def_rw("average_ppm", &OpenMS::FragmentMassError::Statistics::average_ppm)
        .def_rw("variance_ppm", &OpenMS::FragmentMassError::Statistics::variance_ppm)
        ;

    nb::class_<OpenMS::FragmentMassError, OpenMS::QCBase>(m, "FragmentMassError")
        .def(nb::init<>())
        .def("compute", nb::overload_cast<OpenMS::FeatureMap&, const OpenMS::MSExperiment&, const OpenMS::QCBase::SpectraMap&, OpenMS::QCBase::ToleranceUnit, double>(&OpenMS::FragmentMassError::compute),
             "fmap"_a, "exp"_a, "map_to_spectrum"_a, "tolerance_unit"_a = OpenMS::QCBase::ToleranceUnit::AUTO, "tolerance"_a = 20.0)
        .def("compute", nb::overload_cast<OpenMS::PeptideIdentificationList&, const OpenMS::ProteinIdentification::SearchParameters&, const OpenMS::MSExperiment&, const OpenMS::QCBase::SpectraMap&, OpenMS::QCBase::ToleranceUnit, double>(&OpenMS::FragmentMassError::compute),
             "pep_ids"_a, "search_params"_a, "exp"_a, "map_to_spectrum"_a, "tolerance_unit"_a = OpenMS::QCBase::ToleranceUnit::AUTO, "tolerance"_a = 20.0)
        .def("getResults", &OpenMS::FragmentMassError::getResults, nb::rv_policy::reference_internal)
        ;

    // -------------------------------------------------------------------------
    // SpectrumCount
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumCount, OpenMS::QCBase>(m, "SpectrumCount")
        .def(nb::init<>())
        .def("compute", &OpenMS::SpectrumCount::compute, "exp"_a)
        ;

    // -------------------------------------------------------------------------
    // PeptideMass
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::PeptideMass, OpenMS::QCBase>(m, "PeptideMass")
        .def(nb::init<>())
        .def("compute", &OpenMS::PeptideMass::compute, "features"_a)
        ;

    // -------------------------------------------------------------------------
    // RTAlignment
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::RTAlignment, OpenMS::QCBase>(m, "RTAlignment")
        .def(nb::init<>())
        .def("compute", nb::overload_cast<OpenMS::FeatureMap&, const OpenMS::TransformationDescription&>(&OpenMS::RTAlignment::compute, nb::const_), "fm"_a, "trafo"_a)
        .def("compute", nb::overload_cast<OpenMS::PeptideIdentificationList&, const OpenMS::TransformationDescription&>(&OpenMS::RTAlignment::compute, nb::const_), "ids"_a, "trafo"_a)
        ;

    // -------------------------------------------------------------------------
    // MzCalibration
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::MzCalibration, OpenMS::QCBase>(m, "MzCalibration")
        .def(nb::init<>())
        .def("compute", &OpenMS::MzCalibration::compute, "features"_a, "exp"_a, "map_to_spectrum"_a)
        ;

    // -------------------------------------------------------------------------
    // Contaminants
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::Contaminants::ContaminantsSummary>(m, "ContaminantsSummary")
        .def(nb::init<>())
        .def_rw("assigned_contaminants_ratio", &OpenMS::Contaminants::ContaminantsSummary::assigned_contaminants_ratio)
        .def_rw("unassigned_contaminants_ratio", &OpenMS::Contaminants::ContaminantsSummary::unassigned_contaminants_ratio)
        .def_rw("all_contaminants_ratio", &OpenMS::Contaminants::ContaminantsSummary::all_contaminants_ratio)
        .def_rw("assigned_contaminants_intensity_ratio", &OpenMS::Contaminants::ContaminantsSummary::assigned_contaminants_intensity_ratio)
        .def_rw("empty_features", &OpenMS::Contaminants::ContaminantsSummary::empty_features)
        ;

    nb::class_<OpenMS::Contaminants, OpenMS::QCBase>(m, "Contaminants")
        .def(nb::init<>())
        .def("compute", &OpenMS::Contaminants::compute, "features"_a, "contaminants"_a)
        .def("getResults", &OpenMS::Contaminants::getResults, nb::rv_policy::reference_internal)
        ;

    // -------------------------------------------------------------------------
    // DBSuitability
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::DBSuitability::SuitabilityData>(m, "DBSuitabilityData")
        .def(nb::init<>())
        .def_rw("num_top_novo", &OpenMS::DBSuitability::SuitabilityData::num_top_novo)
        .def_rw("num_top_db", &OpenMS::DBSuitability::SuitabilityData::num_top_db)
        .def_rw("num_interest", &OpenMS::DBSuitability::SuitabilityData::num_interest)
        .def_rw("num_re_ranked", &OpenMS::DBSuitability::SuitabilityData::num_re_ranked)
        .def_rw("cut_off", &OpenMS::DBSuitability::SuitabilityData::cut_off)
        .def_rw("suitability", &OpenMS::DBSuitability::SuitabilityData::suitability)
        .def_rw("suitability_no_rerank", &OpenMS::DBSuitability::SuitabilityData::suitability_no_rerank)
        .def_rw("suitability_corr_no_rerank", &OpenMS::DBSuitability::SuitabilityData::suitability_corr_no_rerank)
        .def("clear", &OpenMS::DBSuitability::SuitabilityData::clear)
        .def("setCorrectionFactor", &OpenMS::DBSuitability::SuitabilityData::setCorrectionFactor, "factor"_a)
        .def("getCorrectionFactor", &OpenMS::DBSuitability::SuitabilityData::getCorrectionFactor)
        .def("getCorrectedNovoHits", &OpenMS::DBSuitability::SuitabilityData::getCorrectedNovoHits)
        .def("getCorrectedSuitability", &OpenMS::DBSuitability::SuitabilityData::getCorrectedSuitability)
        .def("simulateNoReRanking", &OpenMS::DBSuitability::SuitabilityData::simulateNoReRanking)
        ;

    auto dbsuitability_class = nb::class_<OpenMS::DBSuitability>(m, "DBSuitability")
        .def(nb::init<>())
        .def("compute", [](OpenMS::DBSuitability& self, const OpenMS::PeptideIdentificationList& pep_ids, const OpenMS::MSExperiment& exp, const std::vector<OpenMS::FASTAFile::FASTAEntry>& original_fasta, const std::vector<OpenMS::FASTAFile::FASTAEntry>& novo_fasta, const OpenMS::ProteinIdentification::SearchParameters& search_params) {
            OpenMS::PeptideIdentificationList p = pep_ids;
            self.compute(std::move(p), exp, original_fasta, novo_fasta, search_params);
        }, "pep_ids"_a, "exp"_a, "original_fasta"_a, "novo_fasta"_a, "search_params"_a)
        .def("getResults", &OpenMS::DBSuitability::getResults, nb::rv_policy::reference_internal)
        ;
    def_DefaultParamHandler<OpenMS::DBSuitability>(dbsuitability_class);

    // -------------------------------------------------------------------------
    // FeatureSummary
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::FeatureSummary::Result>(m, "FeatureSummaryResult")
        .def(nb::init<>())
        .def_rw("feature_count", &OpenMS::FeatureSummary::Result::feature_count)
        .def_rw("rt_shift_mean", &OpenMS::FeatureSummary::Result::rt_shift_mean)
        .def(nb::self == nb::self)
        ;

    nb::class_<OpenMS::FeatureSummary, OpenMS::QCBase>(m, "FeatureSummary")
        .def(nb::init<>())
        .def("compute", &OpenMS::FeatureSummary::compute, "feature_map"_a)
        ;

    // -------------------------------------------------------------------------
    // IdentificationSummary
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::IdentificationSummary::UniqueID>(m, "UniqueID")
        .def(nb::init<>())
        .def_rw("count", &OpenMS::IdentificationSummary::UniqueID::count)
        .def_rw("fdr_threshold", &OpenMS::IdentificationSummary::UniqueID::fdr_threshold)
        ;

    nb::class_<OpenMS::IdentificationSummary::Result>(m, "IdentificationSummaryResult")
        .def(nb::init<>())
        .def_rw("peptide_spectrum_matches", &OpenMS::IdentificationSummary::Result::peptide_spectrum_matches)
        .def_rw("unique_peptides", &OpenMS::IdentificationSummary::Result::unique_peptides)
        .def_rw("unique_proteins", &OpenMS::IdentificationSummary::Result::unique_proteins)
        .def_rw("missed_cleavages_mean", &OpenMS::IdentificationSummary::Result::missed_cleavages_mean)
        .def_rw("protein_hit_scores_mean", &OpenMS::IdentificationSummary::Result::protein_hit_scores_mean)
        .def_rw("peptide_length_mean", &OpenMS::IdentificationSummary::Result::peptide_length_mean)
        .def(nb::self == nb::self)
        ;

    nb::class_<OpenMS::IdentificationSummary, OpenMS::QCBase>(m, "IdentificationSummary")
        .def(nb::init<>())
        .def("compute", &OpenMS::IdentificationSummary::compute, "prot_ids"_a, "pep_ids"_a)
        ;

    // -------------------------------------------------------------------------
    // Ms2IdentificationRate
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::Ms2IdentificationRate::IdentificationRateData>(m, "IdentificationRateData")
        .def(nb::init<>())
        .def_rw("num_peptide_identification", &OpenMS::Ms2IdentificationRate::IdentificationRateData::num_peptide_identification)
        .def_rw("num_ms2_spectra", &OpenMS::Ms2IdentificationRate::IdentificationRateData::num_ms2_spectra)
        .def_rw("identification_rate", &OpenMS::Ms2IdentificationRate::IdentificationRateData::identification_rate)
        ;

    nb::class_<OpenMS::Ms2IdentificationRate, OpenMS::QCBase>(m, "Ms2IdentificationRate")
        .def(nb::init<>())
        .def("compute", nb::overload_cast<const OpenMS::FeatureMap&, const OpenMS::MSExperiment&, bool>(&OpenMS::Ms2IdentificationRate::compute), "feature_map"_a, "exp"_a, "assume_all_target"_a = false)
        .def("compute", nb::overload_cast<const OpenMS::PeptideIdentificationList&, const OpenMS::MSExperiment&, bool>(&OpenMS::Ms2IdentificationRate::compute), "pep_ids"_a, "exp"_a, "assume_all_target"_a = false)
        .def("getResults", &OpenMS::Ms2IdentificationRate::getResults, nb::rv_policy::reference_internal)
        .def("addMetaDataMetricsToMzTab", &OpenMS::Ms2IdentificationRate::addMetaDataMetricsToMzTab, "meta"_a)
        ;

    // -------------------------------------------------------------------------
    // Ms2SpectrumStats
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::Ms2SpectrumStats::ScanEvent>(m, "Ms2SpectrumStatsScanEvent")
        .def(nb::init<UInt32, bool>())
        .def_rw("scan_event_number", &OpenMS::Ms2SpectrumStats::ScanEvent::scan_event_number)
        .def_rw("ms2_presence", &OpenMS::Ms2SpectrumStats::ScanEvent::ms2_presence)
        ;

    nb::class_<OpenMS::Ms2SpectrumStats, OpenMS::QCBase>(m, "Ms2SpectrumStats")
        .def(nb::init<>())
        .def("compute", &OpenMS::Ms2SpectrumStats::compute, "exp"_a, "features"_a, "map_to_spectrum"_a)
        ;

    // -------------------------------------------------------------------------
    // PSMExplainedIonCurrent
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::PSMExplainedIonCurrent::Statistics>(m, "PSMExplainedIonCurrentStatistics")
        .def(nb::init<>())
        .def_rw("average_correctness", &OpenMS::PSMExplainedIonCurrent::Statistics::average_correctness)
        .def_rw("variance_correctness", &OpenMS::PSMExplainedIonCurrent::Statistics::variance_correctness)
        ;

    nb::class_<OpenMS::PSMExplainedIonCurrent, OpenMS::QCBase>(m, "PSMExplainedIonCurrent")
        .def(nb::init<>())
        .def("compute", nb::overload_cast<OpenMS::FeatureMap&, const OpenMS::MSExperiment&, const OpenMS::QCBase::SpectraMap&, OpenMS::QCBase::ToleranceUnit, double>(&OpenMS::PSMExplainedIonCurrent::compute),
             "fmap"_a, "exp"_a, "map_to_spectrum"_a, "tolerance_unit"_a = OpenMS::QCBase::ToleranceUnit::AUTO, "tolerance"_a = 20.0)
        .def("compute", nb::overload_cast<OpenMS::PeptideIdentificationList&, const OpenMS::ProteinIdentification::SearchParameters&, const OpenMS::MSExperiment&, const OpenMS::QCBase::SpectraMap&, OpenMS::QCBase::ToleranceUnit, double>(&OpenMS::PSMExplainedIonCurrent::compute),
             "pep_ids"_a, "search_params"_a, "exp"_a, "map_to_spectrum"_a, "tolerance_unit"_a = OpenMS::QCBase::ToleranceUnit::AUTO, "tolerance"_a = 20.0)
        .def("getResults", &OpenMS::PSMExplainedIonCurrent::getResults, nb::rv_policy::reference_internal)
        ;
}
