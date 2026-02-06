// pyOpenMS nanobind bindings
// Domain: featurefinder

#include "all_casters.h"
#include <OpenMS/FEATUREFINDER/BiGaussFitter1D.h>
#include <OpenMS/FEATUREFINDER/BiGaussModel.h>
#include <OpenMS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/FEATUREFINDER/EmgScoring.h>
#include <OpenMS/FEATUREFINDER/GaussTraceFitter.h>
#include <OpenMS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/FEATUREFINDER/IsotopeFitter1D.h>
#include <OpenMS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>
#include <OpenMS/FEATUREFINDER/PeakWidthEstimator.h>
#include <OpenMS/FEATUREFINDER/SeedListGenerator.h>
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

NB_MODULE(_pyopenms_featurefinder, m) {
    m.doc() = "pyOpenMS featurefinder bindings";

    // -----------------------------------------------------------------------
    // BiGaussFitter1D
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::BiGaussFitter1D>(m, "BiGaussFitter1D", 
        R"doc(
BiGaussian distribution fitter (1-dim.) approximated using linear
interpolation
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::BiGaussFitter1D &>())
        ;

    // -----------------------------------------------------------------------
    // EmgScoring
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::EmgScoring>(m, "EmgScoring", "OpenMS class EmgScoring")
        .def(nb::init<>())
        .def("setFitterParam", [](OpenMS::EmgScoring& self, const OpenMS::Param& param) { return self.setFitterParam(param); }, "param"_a)
        .def("getDefaults", [](OpenMS::EmgScoring& self) { return self.getDefaults(); })
        .def("elutionModelFit", [](const OpenMS::EmgScoring& self, const std::vector<OpenMS::DPosition<2>>& current_section, bool smooth_data) { return self.elutionModelFit(current_section, smooth_data); }, "current_section"_a, "smooth_data"_a)
        ;

    // -----------------------------------------------------------------------
    // GaussTraceFitter
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::GaussTraceFitter>(m, "GaussTraceFitter", 
        R"doc(
Fitter for RT profiles using a Gaussian background model * *
@htmlinclude OpenMS_GaussTraceFitter.parameters * * @todo More docu
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::GaussTraceFitter &>())
        .def("fit", [](OpenMS::GaussTraceFitter& self, OpenMS::FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces) { return self.fit(traces); }, "traces"_a, "Override important methods")
        .def("getLowerRTBound", [](const OpenMS::GaussTraceFitter& self) { return self.getLowerRTBound(); }, "Returns the lower RT bound")
        .def("getUpperRTBound", [](const OpenMS::GaussTraceFitter& self) { return self.getUpperRTBound(); }, "Returns the upper RT bound")
        .def("getHeight", [](const OpenMS::GaussTraceFitter& self) { return self.getHeight(); }, "Returns height of the fitted gaussian model")
        .def("getCenter", [](const OpenMS::GaussTraceFitter& self) { return self.getCenter(); }, "Returns center of the fitted gaussian model")
        .def("getFWHM", [](const OpenMS::GaussTraceFitter& self) { return self.getFWHM(); }, "Returns FWHM of the fitted gaussian model")
        .def("getSigma", [](const OpenMS::GaussTraceFitter& self) { return self.getSigma(); }, "Returns Sigma of the fitted gaussian model")
        .def("checkMaximalRTSpan", [](OpenMS::GaussTraceFitter& self, double max_rt_span) { return self.checkMaximalRTSpan(max_rt_span); }, "max_rt_span"_a)
        .def("checkMinimalRTSpan", [](OpenMS::GaussTraceFitter& self, const std::pair<double, double>& rt_bounds, double min_rt_span) { return self.checkMinimalRTSpan(rt_bounds, min_rt_span); }, "rt_bounds"_a, "min_rt_span"_a)
        .def("getValue", [](const OpenMS::GaussTraceFitter& self, double rt) { return self.getValue(rt); }, "rt"_a, "Returns value of the fitted gaussian model")
        .def("getArea", [](OpenMS::GaussTraceFitter& self) { return self.getArea(); }, "Returns area of the fitted gaussian model")
        .def("getGnuplotFormula", [](OpenMS::GaussTraceFitter& self, const OpenMS::FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, char function_name, double baseline, double rt_shift) { return self.getGnuplotFormula(trace, function_name, baseline, rt_shift); }, "trace"_a, "function_name"_a, "baseline"_a, "rt_shift"_a)
        .def("computeTheoretical", [](const OpenMS::GaussTraceFitter& self, const OpenMS::FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, unsigned long k) { return self.computeTheoretical(trace, k); }, "trace"_a, "k"_a)
        ;

    // -----------------------------------------------------------------------
    // InterpolationModel
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::InterpolationModel>(m, "InterpolationModel", 
        R"doc(
Abstract class for 1D-models that are approximated using linear
interpolation
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::InterpolationModel &>())
        .def("getIntensity", [](const OpenMS::InterpolationModel& self, const OpenMS::DPosition<1>& pos) { return self.getIntensity(pos); }, "pos"_a, "Access model predicted intensity at position 'pos'")
        .def("getIntensity", [](const OpenMS::InterpolationModel& self, double coord) { return self.getIntensity(coord); }, "coord"_a, "Access model predicted intensity at position 'pos'")
        .def("getInterpolation", [](const OpenMS::InterpolationModel& self) -> const OpenMS::Math::LinearInterpolation<> & { return self.getInterpolation(); }, nb::rv_policy::reference_internal, "Returns the interpolation class")
        .def("getScalingFactor", [](const OpenMS::InterpolationModel& self) { return self.getScalingFactor(); }, "Returns the interpolation class")
        .def("setOffset", [](OpenMS::InterpolationModel& self, double offset) { return self.setOffset(offset); }, "offset"_a, "Sets the offset of the model")
        .def("getCenter", [](const OpenMS::InterpolationModel& self) { return self.getCenter(); }, 
            R"doc(
Returns the "center" of the model, particular definition (depends on the derived model)
)doc")
        .def("setSamples", [](OpenMS::InterpolationModel& self) { return self.setSamples(); }, "Sets sample/supporting points of interpolation wrt params")
        .def("setInterpolationStep", [](OpenMS::InterpolationModel& self, double interpolation_step) { return self.setInterpolationStep(interpolation_step); }, "interpolation_step"_a, "Sets the interpolation step for the linear interpolation of the model")
        .def("setScalingFactor", [](OpenMS::InterpolationModel& self, double scaling) { return self.setScalingFactor(scaling); }, "scaling"_a, "Sets the scaling factor of the model")
        ;

    // -----------------------------------------------------------------------
    // BiGaussModel
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::BiGaussModel, OpenMS::InterpolationModel>(m, "BiGaussModel", "BiGaussian distribution approximated using linear interpolation")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::BiGaussModel &>())
        .def("setOffset", [](OpenMS::BiGaussModel& self, double offset) { return self.setOffset(offset); }, "offset"_a)
        .def("setSamples", [](OpenMS::BiGaussModel& self) { return self.setSamples(); })
        .def("getCenter", [](const OpenMS::BiGaussModel& self) { return self.getCenter(); })
        ;

    // -----------------------------------------------------------------------
    // EmgModel
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::EmgModel, OpenMS::InterpolationModel>(m, "EmgModel", 
        R"doc(
Exponentially modified gaussian distribution model for elution
profiles
InterpolationModel
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::EmgModel &>())
        .def("setOffset", [](OpenMS::EmgModel& self, double offset) { return self.setOffset(offset); }, "offset"_a, "Sets the offset of the model")
        .def("setSamples", [](OpenMS::EmgModel& self) { return self.setSamples(); }, "Sets sample/supporting points of interpolation wrt params")
        .def("getCenter", [](const OpenMS::EmgModel& self) { return self.getCenter(); }, 
            R"doc(
Returns the "center" of the model, particular definition (depends on the derived model)
)doc")
        .def("getIntensity", [](const OpenMS::EmgModel& self, const OpenMS::DPosition<1>& pos) { return self.getIntensity(pos); }, "pos"_a, "Access model predicted intensity at position 'pos'")
        .def("getInterpolation", [](const OpenMS::EmgModel& self) -> const OpenMS::Math::LinearInterpolation<> & { return self.getInterpolation(); }, nb::rv_policy::reference_internal, "Returns the interpolation class")
        .def("getScalingFactor", [](const OpenMS::EmgModel& self) { return self.getScalingFactor(); }, "Returns the interpolation class")
        .def("setInterpolationStep", [](OpenMS::EmgModel& self, double interpolation_step) { return self.setInterpolationStep(interpolation_step); }, "interpolation_step"_a, "Sets the interpolation step for the linear interpolation of the model")
        .def("setScalingFactor", [](OpenMS::EmgModel& self, double scaling) { return self.setScalingFactor(scaling); }, "scaling"_a, "Sets the scaling factor of the model")
        ;

    // -----------------------------------------------------------------------
    // IsotopeFitter1D
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IsotopeFitter1D>(m, "IsotopeFitter1D", 
        R"doc(
Isotope distribution fitter (1-dim.) approximated using linear
interpolation
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IsotopeFitter1D &>())
        ;

    // -----------------------------------------------------------------------
    // IsotopeModel
    // -----------------------------------------------------------------------
    auto isotopemodel_class = nb::class_<OpenMS::IsotopeModel, OpenMS::InterpolationModel>(m, "IsotopeModel", 
        R"doc(
Isotope distribution approximated using linear interpolation
This models a smoothed (widened) distribution, i.e. can be used to sample actual raw peaks (depending on the points you query)
If you only want the distribution (no widening), use either
EmpiricalFormula::getIsotopeDistribution() // for a certain sum formula
or
IsotopeDistribution::estimateFromPeptideWeight (double average_weight)  // for averagine
Peak widening is achieved by either a Gaussian or Lorentzian shape
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IsotopeModel &>())
        .def("getCharge", [](const OpenMS::IsotopeModel& self) { return self.getCharge(); })
        .def("setOffset", [](OpenMS::IsotopeModel& self, double offset) { return self.setOffset(offset); }, "offset"_a)
        .def("getOffset", [](OpenMS::IsotopeModel& self) { return self.getOffset(); }, "Get the offset of the model")
        .def("getFormula", [](OpenMS::IsotopeModel& self) { return self.getFormula(); }, "Return the Averagine peptide formula (mass calculated from mean mass and charge -- use .setParameters() to set them)")
        .def("setSamples", [](OpenMS::IsotopeModel& self, const OpenMS::EmpiricalFormula& formula) { return self.setSamples(formula); }, "formula"_a, "Set sample/supporting points of interpolation")
        .def("getCenter", [](const OpenMS::IsotopeModel& self) { return self.getCenter(); })
        .def("getIsotopeDistribution", [](const OpenMS::IsotopeModel& self) -> const OpenMS::IsotopeDistribution & { return self.getIsotopeDistribution(); }, nb::rv_policy::reference_internal, 
            R"doc(
Get the center of the Isotope model
This is a m/z-value not necessarily the monoisotopic mass
)doc")
        ;
    // Averagines enum nested under IsotopeModel
    nb::enum_<OpenMS::IsotopeModel::Averagines>(isotopemodel_class, "Averagines", nb::is_arithmetic())
        .value("C", OpenMS::IsotopeModel::Averagines::C)
        .value("H", OpenMS::IsotopeModel::Averagines::H)
        .value("N", OpenMS::IsotopeModel::Averagines::N)
        .value("O", OpenMS::IsotopeModel::Averagines::O)
        .value("S", OpenMS::IsotopeModel::Averagines::S)
        .value("AVERAGINE_NUM", OpenMS::IsotopeModel::Averagines::AVERAGINE_NUM)
        .export_values();

    // -----------------------------------------------------------------------
    // MultiplexDeltaMasses
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MultiplexDeltaMasses>(m, "MultiplexDeltaMasses", 
        R"doc(
Data structure for mass shift pattern
Groups of labelled peptides appear with characteristic mass shifts
For example, for an Arg6 labeled SILAC peptide pair we expect to see
mass shifts of 0 and 6 Da. Or as second example, for a
peptide pair of a dimethyl labelled sample with a single lysine
we will see mass shifts of 56 Da and 64 Da.
28 Da (N-term) + 28 Da (K) and 34 Da (N-term) + 34 Da (K)
for light and heavy partners respectively
The data structure stores the mass shifts and corresponding labels
for a group of matching peptide features
)doc")
        .def(nb::init<>())
        .def(nb::init<std::vector<OpenMS::MultiplexDeltaMasses::DeltaMass>>())
        .def("getDeltaMasses", [](OpenMS::MultiplexDeltaMasses& self) -> std::vector<OpenMS::MultiplexDeltaMasses::DeltaMass> & { return self.getDeltaMasses(); }, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // MultiplexDeltaMassesGenerator_Label
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MultiplexDeltaMassesGenerator::Label>(m, "MultiplexDeltaMassesGenerator_Label", "OpenMS class MultiplexDeltaMassesGenerator_Label")
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::String, double>())
        .def_rw("short_name", &OpenMS::MultiplexDeltaMassesGenerator::Label::short_name)
        .def_rw("long_name", &OpenMS::MultiplexDeltaMassesGenerator::Label::long_name)
        .def_rw("description", &OpenMS::MultiplexDeltaMassesGenerator::Label::description)
        .def_rw("delta_mass", &OpenMS::MultiplexDeltaMassesGenerator::Label::delta_mass)
        ;

    // -----------------------------------------------------------------------
    // MultiplexDeltaMasses_DeltaMass
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MultiplexDeltaMasses::DeltaMass>(m, "MultiplexDeltaMasses_DeltaMass", "OpenMS class MultiplexDeltaMasses_DeltaMass")
        .def(nb::init<double, std::multiset<OpenMS::String>>())
        .def(nb::init<double, OpenMS::String>())
        .def_rw("delta_mass", &OpenMS::MultiplexDeltaMasses::DeltaMass::delta_mass)
        .def_rw("label_set", &OpenMS::MultiplexDeltaMasses::DeltaMass::label_set)
        ;

    // -----------------------------------------------------------------------
    // MultiplexIsotopicPeakPattern
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MultiplexIsotopicPeakPattern>(m, "MultiplexIsotopicPeakPattern", 
        R"doc(
data structure for pattern of isotopic peaks * * Groups of peptides
appear as characteristic patterns of isotopic peaks * in MS1 spectra.
For example, for an Arg6 labeled SILAC peptide pair * of charge 2+
with three isotopic peaks we expect peaks * at relative m/z shifts of
0, 0.5, 1, 3, 3.5 and 4 Th
)doc")
        .def(nb::init<int, int, OpenMS::MultiplexDeltaMasses, int>())
        .def("getCharge", [](const OpenMS::MultiplexIsotopicPeakPattern& self) { return self.getCharge(); }, "Returns charge")
        .def("getPeaksPerPeptide", [](const OpenMS::MultiplexIsotopicPeakPattern& self) { return self.getPeaksPerPeptide(); }, "Returns peaks per peptide")
        .def("getMassShifts", [](const OpenMS::MultiplexIsotopicPeakPattern& self) { return self.getMassShifts(); }, "Returns mass shifts")
        .def("getMassShiftIndex", [](const OpenMS::MultiplexIsotopicPeakPattern& self) { return self.getMassShiftIndex(); }, "Returns mass shift index")
        .def("getMassShiftCount", [](const OpenMS::MultiplexIsotopicPeakPattern& self) { return self.getMassShiftCount(); }, "Returns number of mass shifts i.e. the number of peptides in the multiplet")
        .def("getMassShiftAt", [](const OpenMS::MultiplexIsotopicPeakPattern& self, unsigned long i) { return self.getMassShiftAt(i); }, "i"_a, "Returns mass shift at position i")
        .def("getMZShiftAt", [](const OpenMS::MultiplexIsotopicPeakPattern& self, unsigned long i) { return self.getMZShiftAt(i); }, "i"_a, "Returns m/z shift at position i")
        .def("getMZShiftCount", [](const OpenMS::MultiplexIsotopicPeakPattern& self) { return self.getMZShiftCount(); }, "Returns number of m/z shifts")
        ;

    // -----------------------------------------------------------------------
    // PeakWidthEstimator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakWidthEstimator>(m, "PeakWidthEstimator", 
        R"doc(
Rough estimation of the peak width at m/z
Based on the peaks of the dataset (peak position & width) and the peak
boundaries as reported by the PeakPickerHiRes, the typical peak width is
estimated for arbitrary m/z using a spline interpolation.
private
PeakWidthEstimator() except + nogil
)doc")
        .def(nb::init<OpenMS::MSExperiment, std::vector<std::vector<OpenMS::PeakPickerHiRes::PeakBoundary>>>())
        .def("getPeakWidth", [](OpenMS::PeakWidthEstimator& self, double mz) { return self.getPeakWidth(mz); }, "mz"_a, "Returns the estimated peak width at m/z")
        ;

    // -----------------------------------------------------------------------
    // SeedListGenerator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SeedListGenerator>(m, "SeedListGenerator", "Generate seed lists for feature detection")
        .def(nb::init<>())
        .def("generateSeedList", [](OpenMS::SeedListGenerator& self, const OpenMS::PeakMap& experiment) {
            OpenMS::SeedListGenerator::SeedList seeds;
            self.generateSeedList(experiment, seeds);
            return seeds;
        }, "experiment"_a, "Generate a seed list from an experiment using precursor positions")
        .def("generateSeedList", [](OpenMS::SeedListGenerator& self, OpenMS::PeptideIdentificationList& peptides, bool use_peptide_mass) {
            OpenMS::SeedListGenerator::SeedList seeds;
            self.generateSeedList(peptides, seeds, use_peptide_mass);
            return seeds;
        }, "peptides"_a, "use_peptide_mass"_a = false, "Generate a seed list from peptide identifications")
        .def("generateSeedLists", [](OpenMS::SeedListGenerator& self, const OpenMS::ConsensusMap& consensus) {
            std::map<OpenMS::UInt64, OpenMS::SeedListGenerator::SeedList> seed_lists;
            self.generateSeedLists(consensus, seed_lists);
            return seed_lists;
        }, "consensus"_a, "Generate one seed list per constituent map in a consensus map")
        .def("convertSeedList", [](OpenMS::SeedListGenerator& self, const OpenMS::SeedListGenerator::SeedList& seeds) {
            OpenMS::FeatureMap features;
            self.convertSeedList(seeds, features);
            return features;
        }, "seeds"_a, "Convert seed positions to a FeatureMap")
        .def("convertSeedList", [](OpenMS::SeedListGenerator& self, const OpenMS::FeatureMap& features) {
            OpenMS::SeedListGenerator::SeedList seeds;
            self.convertSeedList(features, seeds);
            return seeds;
        }, "features"_a, "Convert a FeatureMap of seeds back to seed positions")
        ;

}
