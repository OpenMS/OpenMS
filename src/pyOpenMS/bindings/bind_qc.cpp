// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Krrish Verma $
// $Authors: Krrish Verma $
// --------------------------------------------------------------------------

// pyOpenMS nanobind bindings
// Domain: qc

#include "all_casters.h"

// QC headers (forward-declare many types; pull in full definitions below)
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

// Full definitions required by nanobind for forward-declared types in QC headers
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include <nanobind/operators.h>

#include "binding_utils.h"

using namespace nb::literals;

NB_MODULE(_pyopenms_qc, m) {
    m.doc() = "pyOpenMS Quality Control (QC) bindings";

    // -------------------------------------------------------------------------
    // QCBase Enums
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

    // -------------------------------------------------------------------------
    // FlagSet<QCBase::Requires> a.k.a. QCBase::Status
    // Use explicit lambdas instead of nb::self for non-const operators
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::FlagSet<OpenMS::QCBase::Requires>>(m, "QCBaseStatus")
        .def(nb::init<>())
        .def(nb::init<OpenMS::QCBase::Requires>())
        .def(nb::init<const OpenMS::FlagSet<OpenMS::QCBase::Requires>&>())
        .def("__copy__",  [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& self, nb::dict&) { return self; }, "memo"_a)
        .def("__eq__",  [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& a, const OpenMS::FlagSet<OpenMS::QCBase::Requires>& b) { return a == b; })
        // bitwise OR (returns new object)
        .def("__or__",  [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& a, const OpenMS::QCBase::Requires& b) { return a | b; })
        .def("__or__",  [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& a, const OpenMS::FlagSet<OpenMS::QCBase::Requires>& b) { return a | b; })
        // bitwise AND (returns new object)
        .def("__and__", [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& a, const OpenMS::QCBase::Requires& b) { return a & b; })
        .def("__and__", [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& a, const OpenMS::FlagSet<OpenMS::QCBase::Requires>& b) { return a & b; })
        // subtraction (returns new object; note: operator- is non-const in FlagSet)
        .def("__sub__", [](OpenMS::FlagSet<OpenMS::QCBase::Requires> a, const OpenMS::QCBase::Requires& b) { return a - b; })
        .def("__sub__", [](OpenMS::FlagSet<OpenMS::QCBase::Requires> a, const OpenMS::FlagSet<OpenMS::QCBase::Requires>& b) { return a - b; })
        // in-place operators
        .def("__ior__",  [](OpenMS::FlagSet<OpenMS::QCBase::Requires>& a, const OpenMS::QCBase::Requires& b) -> OpenMS::FlagSet<OpenMS::QCBase::Requires>& { a |= b; return a; })
        .def("__ior__",  [](OpenMS::FlagSet<OpenMS::QCBase::Requires>& a, const OpenMS::FlagSet<OpenMS::QCBase::Requires>& b) -> OpenMS::FlagSet<OpenMS::QCBase::Requires>& { a |= b; return a; })
        .def("__iand__", [](OpenMS::FlagSet<OpenMS::QCBase::Requires>& a, const OpenMS::QCBase::Requires& b) -> OpenMS::FlagSet<OpenMS::QCBase::Requires>& { a &= b; return a; })
        .def("__iand__", [](OpenMS::FlagSet<OpenMS::QCBase::Requires>& a, const OpenMS::FlagSet<OpenMS::QCBase::Requires>& b) -> OpenMS::FlagSet<OpenMS::QCBase::Requires>& { a &= b; return a; })
        .def("isSuperSetOf", [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& self, const OpenMS::FlagSet<OpenMS::QCBase::Requires>& required) { return self.isSuperSetOf(required); })
        .def("isSuperSetOf", [](const OpenMS::FlagSet<OpenMS::QCBase::Requires>& self, const OpenMS::QCBase::Requires& required) { return self.isSuperSetOf(required); })
        .def("empty", &OpenMS::FlagSet<OpenMS::QCBase::Requires>::empty)
        .def("value", &OpenMS::FlagSet<OpenMS::QCBase::Requires>::value)
        ;

    // -------------------------------------------------------------------------
    // QCBase::SpectraMap
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::QCBase::SpectraMap>(m, "QCSpectraMap")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MSExperiment&>())
        .def(nb::init<const OpenMS::QCBase::SpectraMap&>())
        .def("__copy__",  [](const OpenMS::QCBase::SpectraMap& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::QCBase::SpectraMap& self, nb::dict&) { return self; }, "memo"_a)
        .def("calculateMap", &OpenMS::QCBase::SpectraMap::calculateMap, "exp"_a)
        .def("at",   &OpenMS::QCBase::SpectraMap::at,   "identifier"_a)
        .def("clear",&OpenMS::QCBase::SpectraMap::clear)
        .def("empty",&OpenMS::QCBase::SpectraMap::empty)
        .def("size", &OpenMS::QCBase::SpectraMap::size)
        ;

    // -------------------------------------------------------------------------
    // QCBase (abstract base)
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::QCBase>(m, "QCBase")
        .def("getName",        &OpenMS::QCBase::getName, nb::rv_policy::reference_internal)
        .def("requirements",   &OpenMS::QCBase::requirements)
        .def("isRunnable",     &OpenMS::QCBase::isRunnable, "s"_a)
        .def_static("isLabeledExperiment", &OpenMS::QCBase::isLabeledExperiment, "cm"_a)
        ;

    // -------------------------------------------------------------------------
    // FWHM
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::FWHM, OpenMS::QCBase>(m, "FWHM")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FWHM&>())
        .def("__copy__",     [](const OpenMS::FWHM& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::FWHM& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute",    &OpenMS::FWHM::compute, "features"_a)
        ;

    // -------------------------------------------------------------------------
    // TIC::Result and TIC
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::TIC::Result>(m, "TICResult")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TIC::Result&>())
        .def("__copy__",     [](const OpenMS::TIC::Result& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::TIC::Result& self, nb::dict&) { return self; }, "memo"_a)
        .def_rw("intensities",          &OpenMS::TIC::Result::intensities)
        .def_rw("relative_intensities", &OpenMS::TIC::Result::relative_intensities)
        .def_rw("retention_times",      &OpenMS::TIC::Result::retention_times)
        .def_rw("area",  &OpenMS::TIC::Result::area)
        .def_rw("fall",  &OpenMS::TIC::Result::fall)
        .def_rw("jump",  &OpenMS::TIC::Result::jump)
        .def("__eq__",   [](const OpenMS::TIC::Result& a, const OpenMS::TIC::Result& b) { return a == b; })
        ;

    nb::class_<OpenMS::TIC, OpenMS::QCBase>(m, "TIC")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TIC&>())
        .def("__copy__",     [](const OpenMS::TIC& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::TIC& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute",    &OpenMS::TIC::compute, "exp"_a, "bin_size"_a = 0.0f, "ms_level"_a = 1)
        .def("addMetaDataMetricsToMzTab", &OpenMS::TIC::addMetaDataMetricsToMzTab, "meta"_a, "tics"_a)
        ;

    // -------------------------------------------------------------------------
    // MissedCleavages
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::MissedCleavages, OpenMS::QCBase>(m, "MissedCleavages")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MissedCleavages&>())
        .def("__copy__",     [](const OpenMS::MissedCleavages& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::MissedCleavages& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute", nb::overload_cast<OpenMS::FeatureMap&>(&OpenMS::MissedCleavages::compute), "fmap"_a)
        .def("compute", nb::overload_cast<std::vector<OpenMS::ProteinIdentification>&, OpenMS::PeptideIdentificationList&>(&OpenMS::MissedCleavages::compute), "prot_ids"_a, "pep_ids"_a)
        .def("getResults", &OpenMS::MissedCleavages::getResults, nb::rv_policy::reference_internal)
        ;

    // -------------------------------------------------------------------------
    // FragmentMassError::Statistics and FragmentMassError
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::FragmentMassError::Statistics>(m, "FragmentMassErrorStatistics")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FragmentMassError::Statistics&>())
        .def("__copy__",     [](const OpenMS::FragmentMassError::Statistics& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::FragmentMassError::Statistics& self, nb::dict&) { return self; }, "memo"_a)
        .def_rw("average_ppm",  &OpenMS::FragmentMassError::Statistics::average_ppm)
        .def_rw("variance_ppm", &OpenMS::FragmentMassError::Statistics::variance_ppm)
        ;

    nb::class_<OpenMS::FragmentMassError, OpenMS::QCBase>(m, "FragmentMassError")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FragmentMassError&>())
        .def("__copy__",     [](const OpenMS::FragmentMassError& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::FragmentMassError& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute",
             nb::overload_cast<OpenMS::FeatureMap&, const OpenMS::MSExperiment&, const OpenMS::QCBase::SpectraMap&, OpenMS::QCBase::ToleranceUnit, double>(&OpenMS::FragmentMassError::compute),
             "fmap"_a, "exp"_a, "map_to_spectrum"_a,
             "tolerance_unit"_a = OpenMS::QCBase::ToleranceUnit::AUTO, "tolerance"_a = 20.0)
        .def("compute",
             nb::overload_cast<OpenMS::PeptideIdentificationList&, const OpenMS::ProteinIdentification::SearchParameters&, const OpenMS::MSExperiment&, const OpenMS::QCBase::SpectraMap&, OpenMS::QCBase::ToleranceUnit, double>(&OpenMS::FragmentMassError::compute),
             "pep_ids"_a, "search_params"_a, "exp"_a, "map_to_spectrum"_a,
             "tolerance_unit"_a = OpenMS::QCBase::ToleranceUnit::AUTO, "tolerance"_a = 20.0)
        .def("getResults", &OpenMS::FragmentMassError::getResults, nb::rv_policy::reference_internal)
        ;

    // -------------------------------------------------------------------------
    // SpectrumCount
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumCount, OpenMS::QCBase>(m, "SpectrumCount")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SpectrumCount&>())
        .def("__copy__",     [](const OpenMS::SpectrumCount& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::SpectrumCount& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute", &OpenMS::SpectrumCount::compute, "exp"_a)
        ;

    // -------------------------------------------------------------------------
    // PeptideMass
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::PeptideMass, OpenMS::QCBase>(m, "PeptideMass")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PeptideMass&>())
        .def("__copy__",     [](const OpenMS::PeptideMass& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::PeptideMass& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute", &OpenMS::PeptideMass::compute, "features"_a)
        ;

    // -------------------------------------------------------------------------
    // RTAlignment
    // RTAlignment::compute is const-overloaded; use explicit const pointer cast
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::RTAlignment, OpenMS::QCBase>(m, "RTAlignment")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::RTAlignment&>())
        .def("__copy__",     [](const OpenMS::RTAlignment& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::RTAlignment& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute",
             [](const OpenMS::RTAlignment& self, OpenMS::FeatureMap& fm, const OpenMS::TransformationDescription& trafo) {
                 self.compute(fm, trafo);
             }, "fm"_a, "trafo"_a)
        .def("compute",
             [](const OpenMS::RTAlignment& self, OpenMS::PeptideIdentificationList& ids, const OpenMS::TransformationDescription& trafo) {
                 self.compute(ids, trafo);
             }, "ids"_a, "trafo"_a)
        ;

    // -------------------------------------------------------------------------
    // MzCalibration
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::MzCalibration, OpenMS::QCBase>(m, "MzCalibration")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MzCalibration&>())
        .def("__copy__",     [](const OpenMS::MzCalibration& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::MzCalibration& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute", &OpenMS::MzCalibration::compute, "features"_a, "exp"_a, "map_to_spectrum"_a)
        ;

    // -------------------------------------------------------------------------
    // Contaminants::ContaminantsSummary and Contaminants
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::Contaminants::ContaminantsSummary>(m, "ContaminantsSummary")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Contaminants::ContaminantsSummary&>())
        .def("__copy__",     [](const OpenMS::Contaminants::ContaminantsSummary& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::Contaminants::ContaminantsSummary& self, nb::dict&) { return self; }, "memo"_a)
        .def_rw("assigned_contaminants_ratio",           &OpenMS::Contaminants::ContaminantsSummary::assigned_contaminants_ratio)
        .def_rw("unassigned_contaminants_ratio",         &OpenMS::Contaminants::ContaminantsSummary::unassigned_contaminants_ratio)
        .def_rw("all_contaminants_ratio",                &OpenMS::Contaminants::ContaminantsSummary::all_contaminants_ratio)
        .def_rw("assigned_contaminants_intensity_ratio", &OpenMS::Contaminants::ContaminantsSummary::assigned_contaminants_intensity_ratio)
        .def_rw("empty_features",                        &OpenMS::Contaminants::ContaminantsSummary::empty_features)
        ;

    nb::class_<OpenMS::Contaminants, OpenMS::QCBase>(m, "Contaminants")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Contaminants&>())
        .def("__copy__",     [](const OpenMS::Contaminants& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::Contaminants& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute",    &OpenMS::Contaminants::compute, "features"_a, "contaminants"_a)
        .def("getResults", &OpenMS::Contaminants::getResults, nb::rv_policy::reference_internal)
        ;

    // -------------------------------------------------------------------------
    // DBSuitability::SuitabilityData and DBSuitability
    // Note: DBSuitability inherits from DefaultParamHandler, NOT QCBase
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::DBSuitability::SuitabilityData>(m, "DBSuitabilityData")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DBSuitability::SuitabilityData&>())
        .def("__copy__",     [](const OpenMS::DBSuitability::SuitabilityData& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::DBSuitability::SuitabilityData& self, nb::dict&) { return self; }, "memo"_a)
        .def_rw("num_top_novo",              &OpenMS::DBSuitability::SuitabilityData::num_top_novo)
        .def_rw("num_top_db",                &OpenMS::DBSuitability::SuitabilityData::num_top_db)
        .def_rw("num_interest",              &OpenMS::DBSuitability::SuitabilityData::num_interest)
        .def_rw("num_re_ranked",             &OpenMS::DBSuitability::SuitabilityData::num_re_ranked)
        .def_rw("cut_off",                   &OpenMS::DBSuitability::SuitabilityData::cut_off)
        .def_rw("suitability",               &OpenMS::DBSuitability::SuitabilityData::suitability)
        .def_rw("suitability_no_rerank",     &OpenMS::DBSuitability::SuitabilityData::suitability_no_rerank)
        .def_rw("suitability_corr_no_rerank",&OpenMS::DBSuitability::SuitabilityData::suitability_corr_no_rerank)
        .def("clear",                 &OpenMS::DBSuitability::SuitabilityData::clear)
        .def("setCorrectionFactor",   &OpenMS::DBSuitability::SuitabilityData::setCorrectionFactor,   "factor"_a)
        .def("getCorrectionFactor",   &OpenMS::DBSuitability::SuitabilityData::getCorrectionFactor)
        .def("getCorrectedNovoHits",  &OpenMS::DBSuitability::SuitabilityData::getCorrectedNovoHits)
        .def("getCorrectedSuitability",&OpenMS::DBSuitability::SuitabilityData::getCorrectedSuitability)
        .def("simulateNoReRanking",   &OpenMS::DBSuitability::SuitabilityData::simulateNoReRanking)
        ;

    auto dbs_cls = nb::class_<OpenMS::DBSuitability>(m, "DBSuitability")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DBSuitability&>())
        .def("__copy__",     [](const OpenMS::DBSuitability& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::DBSuitability& self, nb::dict&) { return self; }, "memo"_a)
        // compute takes PeptideIdentificationList&&; copy in Python, move in C++
        .def("compute",
             [](OpenMS::DBSuitability& self,
                OpenMS::PeptideIdentificationList pep_ids,
                const OpenMS::MSExperiment& exp,
                const std::vector<OpenMS::FASTAFile::FASTAEntry>& original_fasta,
                const std::vector<OpenMS::FASTAFile::FASTAEntry>& novo_fasta,
                const OpenMS::ProteinIdentification::SearchParameters& search_params) {
                 self.compute(std::move(pep_ids), exp, original_fasta, novo_fasta, search_params);
             }, "pep_ids"_a, "exp"_a, "original_fasta"_a, "novo_fasta"_a, "search_params"_a)
        .def("getResults", &OpenMS::DBSuitability::getResults, nb::rv_policy::reference_internal)
        ;
    def_DefaultParamHandler<OpenMS::DBSuitability>(dbs_cls);

    // -------------------------------------------------------------------------
    // FeatureSummary::Result and FeatureSummary
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::FeatureSummary::Result>(m, "FeatureSummaryResult")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FeatureSummary::Result&>())
        .def("__copy__",     [](const OpenMS::FeatureSummary::Result& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::FeatureSummary::Result& self, nb::dict&) { return self; }, "memo"_a)
        .def_rw("feature_count", &OpenMS::FeatureSummary::Result::feature_count)
        .def_rw("rt_shift_mean", &OpenMS::FeatureSummary::Result::rt_shift_mean)
        .def("__eq__", [](const OpenMS::FeatureSummary::Result& a, const OpenMS::FeatureSummary::Result& b) { return a == b; })
        ;

    nb::class_<OpenMS::FeatureSummary, OpenMS::QCBase>(m, "FeatureSummary")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FeatureSummary&>())
        .def("__copy__",     [](const OpenMS::FeatureSummary& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::FeatureSummary& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute", &OpenMS::FeatureSummary::compute, "feature_map"_a)
        ;

    // -------------------------------------------------------------------------
    // IdentificationSummary::UniqueID, ::Result, and IdentificationSummary
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::IdentificationSummary::UniqueID>(m, "UniqueID")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IdentificationSummary::UniqueID&>())
        .def("__copy__",     [](const OpenMS::IdentificationSummary::UniqueID& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::IdentificationSummary::UniqueID& self, nb::dict&) { return self; }, "memo"_a)
        .def_rw("count",         &OpenMS::IdentificationSummary::UniqueID::count)
        .def_rw("fdr_threshold", &OpenMS::IdentificationSummary::UniqueID::fdr_threshold)
        ;

    nb::class_<OpenMS::IdentificationSummary::Result>(m, "IdentificationSummaryResult")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IdentificationSummary::Result&>())
        .def("__copy__",     [](const OpenMS::IdentificationSummary::Result& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::IdentificationSummary::Result& self, nb::dict&) { return self; }, "memo"_a)
        .def_rw("peptide_spectrum_matches", &OpenMS::IdentificationSummary::Result::peptide_spectrum_matches)
        .def_rw("unique_peptides",          &OpenMS::IdentificationSummary::Result::unique_peptides)
        .def_rw("unique_proteins",          &OpenMS::IdentificationSummary::Result::unique_proteins)
        .def_rw("missed_cleavages_mean",    &OpenMS::IdentificationSummary::Result::missed_cleavages_mean)
        .def_rw("protein_hit_scores_mean",  &OpenMS::IdentificationSummary::Result::protein_hit_scores_mean)
        .def_rw("peptide_length_mean",      &OpenMS::IdentificationSummary::Result::peptide_length_mean)
        .def("__eq__", [](const OpenMS::IdentificationSummary::Result& a, const OpenMS::IdentificationSummary::Result& b) { return a == b; })
        ;

    nb::class_<OpenMS::IdentificationSummary, OpenMS::QCBase>(m, "IdentificationSummary")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IdentificationSummary&>())
        .def("__copy__",     [](const OpenMS::IdentificationSummary& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::IdentificationSummary& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute", &OpenMS::IdentificationSummary::compute, "prot_ids"_a, "pep_ids"_a)
        ;

    // -------------------------------------------------------------------------
    // Ms2IdentificationRate::IdentificationRateData and Ms2IdentificationRate
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::Ms2IdentificationRate::IdentificationRateData>(m, "IdentificationRateData")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Ms2IdentificationRate::IdentificationRateData&>())
        .def("__copy__",     [](const OpenMS::Ms2IdentificationRate::IdentificationRateData& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::Ms2IdentificationRate::IdentificationRateData& self, nb::dict&) { return self; }, "memo"_a)
        .def_rw("num_peptide_identification", &OpenMS::Ms2IdentificationRate::IdentificationRateData::num_peptide_identification)
        .def_rw("num_ms2_spectra",            &OpenMS::Ms2IdentificationRate::IdentificationRateData::num_ms2_spectra)
        .def_rw("identification_rate",        &OpenMS::Ms2IdentificationRate::IdentificationRateData::identification_rate)
        ;

    nb::class_<OpenMS::Ms2IdentificationRate, OpenMS::QCBase>(m, "Ms2IdentificationRate")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Ms2IdentificationRate&>())
        .def("__copy__",     [](const OpenMS::Ms2IdentificationRate& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::Ms2IdentificationRate& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute",
             nb::overload_cast<const OpenMS::FeatureMap&, const OpenMS::MSExperiment&, bool>(&OpenMS::Ms2IdentificationRate::compute),
             "feature_map"_a, "exp"_a, "assume_all_target"_a = false)
        .def("compute",
             nb::overload_cast<const OpenMS::PeptideIdentificationList&, const OpenMS::MSExperiment&, bool>(&OpenMS::Ms2IdentificationRate::compute),
             "pep_ids"_a, "exp"_a, "assume_all_target"_a = false)
        .def("getResults",               &OpenMS::Ms2IdentificationRate::getResults, nb::rv_policy::reference_internal)
        .def("addMetaDataMetricsToMzTab",&OpenMS::Ms2IdentificationRate::addMetaDataMetricsToMzTab, "meta"_a)
        ;

    // -------------------------------------------------------------------------
    // Ms2SpectrumStats::ScanEvent and Ms2SpectrumStats
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::Ms2SpectrumStats::ScanEvent>(m, "Ms2SpectrumStatsScanEvent")
        .def(nb::init<OpenMS::UInt32, bool>())
        .def(nb::init<const OpenMS::Ms2SpectrumStats::ScanEvent&>())
        .def("__copy__",     [](const OpenMS::Ms2SpectrumStats::ScanEvent& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::Ms2SpectrumStats::ScanEvent& self, nb::dict&) { return self; }, "memo"_a)
        .def_rw("scan_event_number", &OpenMS::Ms2SpectrumStats::ScanEvent::scan_event_number)
        .def_rw("ms2_presence",      &OpenMS::Ms2SpectrumStats::ScanEvent::ms2_presence)
        ;

    nb::class_<OpenMS::Ms2SpectrumStats, OpenMS::QCBase>(m, "Ms2SpectrumStats")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Ms2SpectrumStats&>())
        .def("__copy__",     [](const OpenMS::Ms2SpectrumStats& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::Ms2SpectrumStats& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute", &OpenMS::Ms2SpectrumStats::compute, "exp"_a, "features"_a, "map_to_spectrum"_a)
        ;

    // -------------------------------------------------------------------------
    // PSMExplainedIonCurrent::Statistics and PSMExplainedIonCurrent
    // -------------------------------------------------------------------------
    nb::class_<OpenMS::PSMExplainedIonCurrent::Statistics>(m, "PSMExplainedIonCurrentStatistics")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PSMExplainedIonCurrent::Statistics&>())
        .def("__copy__",     [](const OpenMS::PSMExplainedIonCurrent::Statistics& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::PSMExplainedIonCurrent::Statistics& self, nb::dict&) { return self; }, "memo"_a)
        .def_rw("average_correctness",  &OpenMS::PSMExplainedIonCurrent::Statistics::average_correctness)
        .def_rw("variance_correctness", &OpenMS::PSMExplainedIonCurrent::Statistics::variance_correctness)
        ;

    nb::class_<OpenMS::PSMExplainedIonCurrent, OpenMS::QCBase>(m, "PSMExplainedIonCurrent")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PSMExplainedIonCurrent&>())
        .def("__copy__",     [](const OpenMS::PSMExplainedIonCurrent& self) { return self; })
        .def("__deepcopy__", [](const OpenMS::PSMExplainedIonCurrent& self, nb::dict&) { return self; }, "memo"_a)
        .def("compute",
             nb::overload_cast<OpenMS::FeatureMap&, const OpenMS::MSExperiment&, const OpenMS::QCBase::SpectraMap&, OpenMS::QCBase::ToleranceUnit, double>(&OpenMS::PSMExplainedIonCurrent::compute),
             "fmap"_a, "exp"_a, "map_to_spectrum"_a,
             "tolerance_unit"_a = OpenMS::QCBase::ToleranceUnit::AUTO, "tolerance"_a = 20.0)
        .def("compute",
             nb::overload_cast<OpenMS::PeptideIdentificationList&, const OpenMS::ProteinIdentification::SearchParameters&, const OpenMS::MSExperiment&, const OpenMS::QCBase::SpectraMap&, OpenMS::QCBase::ToleranceUnit, double>(&OpenMS::PSMExplainedIonCurrent::compute),
             "pep_ids"_a, "search_params"_a, "exp"_a, "map_to_spectrum"_a,
             "tolerance_unit"_a = OpenMS::QCBase::ToleranceUnit::AUTO, "tolerance"_a = 20.0)
        .def("getResults", &OpenMS::PSMExplainedIonCurrent::getResults, nb::rv_policy::reference_internal)
        ;
}
