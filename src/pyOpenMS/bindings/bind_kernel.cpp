// pyOpenMS nanobind bindings
// Domain: kernel

#include "all_casters.h"
#include <type_traits>
#include <limits>
#include <memory>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScores.h>
#include <OpenMS/ANALYSIS/TARGETED/IncludeExcludeTarget.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/BinnedSpectrum.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/KERNEL/ChromatogramRangeManager.h>
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConversionHelper.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/MobilityPeak1D.h>
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/PeakIndex.h>
#include <OpenMS/KERNEL/RangeManager.h>
#include <OpenMS/KERNEL/RichPeak2D.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/KERNEL/SpectrumRangeManager.h>
#include <OpenMS/METADATA/Acquisition.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/METADATA/CVTermListInterface.h>
#include <OpenMS/METADATA/ChromatogramSettings.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/HPLC.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Product.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/ScanWindow.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
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
#include "binding_utils.h"

namespace nb = nanobind;
using namespace nb::literals;

// Helper: ensure an object is a C-contiguous numpy array of type T.
// Accepts numpy arrays (zero-copy if already correct dtype and contiguous) and
// Python lists/tuples (converts via numpy.ascontiguousarray). This maintains
// backward compatibility with code that passes plain Python lists to set_peaks() etc.
template <typename T>
nb::ndarray<nb::numpy, T, nb::ndim<1>, nb::c_contig> as_numpy_array(nb::object obj) {
    nb::ndarray<nb::numpy, T, nb::ndim<1>, nb::c_contig> arr;
    if (nb::try_cast(obj, arr)) return arr;  // already the right type and contiguous — zero-copy
    // Convert lists/sequences/non-contiguous arrays to contiguous numpy array
    auto np = nb::module_::import_("numpy");
    nb::object dtype;
    if constexpr (std::is_same_v<T, float>) dtype = np.attr("float32");
    else dtype = np.attr("float64");
    nb::object np_arr = np.attr("ascontiguousarray")(obj, dtype);
    if (!nb::try_cast(np_arr, arr)) {
        throw std::runtime_error("Failed to convert input to contiguous numpy array");
    }
    return arr;
}

NB_MODULE(_pyopenms_kernel, m) {
    // ABI guards for zero-copy structured array access (get_peaks_struct dtype depends on these)
    static_assert(std::is_standard_layout_v<OpenMS::MobilityPeak1D>,
                  "MobilityPeak1D must be standard-layout for zero-copy struct views (guarantees member order matches dtype)");
    static_assert(sizeof(OpenMS::MobilityPeak1D) == 16,
                  "MobilityPeak1D must be 16 bytes for zero-copy structured array access");
    static_assert(std::is_same_v<OpenMS::MobilityPeak1D::CoordinateType, double>,
                  "MobilityPeak1D::CoordinateType must be double (dtype assumes float64 for position)");
    static_assert(std::is_same_v<OpenMS::MobilityPeak1D::IntensityType, float>,
                  "MobilityPeak1D::IntensityType must be float (dtype assumes float32 for intensity)");

    m.doc() = "pyOpenMS kernel bindings";


    // DriftTimeUnit enum
    // -----------------------------------------------------------------------
    // DriftTimeUnit
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::DriftTimeUnit>(m, "DriftTimeUnit", nb::is_arithmetic(), "Unit of drift time for ion mobility")
        .value("NONE", OpenMS::DriftTimeUnit::NONE)
        .value("MILLISECOND", OpenMS::DriftTimeUnit::MILLISECOND)
        .value("VSSC", OpenMS::DriftTimeUnit::VSSC)
        .value("FAIMS_COMPENSATION_VOLTAGE", OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE)
        .value("CCS", OpenMS::DriftTimeUnit::CCS)
        .value("SIZE_OF_DRIFTTIMEUNIT", OpenMS::DriftTimeUnit::SIZE_OF_DRIFTTIMEUNIT)

        .export_values();


    // IMFormat enum
    // -----------------------------------------------------------------------
    // IMFormat
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::IMFormat>(m, "IMFormat", "Ion mobility data format in an experiment")
        .value("NONE", OpenMS::IMFormat::NONE)
        .value("IM_PEAK", OpenMS::IMFormat::IM_PEAK)
        .value("IM_SPECTRUM", OpenMS::IMFormat::IM_SPECTRUM)
        .value("UNKNOWN", OpenMS::IMFormat::UNKNOWN)
        .value("SIZE_OF_IMFORMAT", OpenMS::IMFormat::SIZE_OF_IMFORMAT)

        .export_values();

    // IMPeakType enum
    // -----------------------------------------------------------------------
    // IMPeakType
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::IMPeakType>(m, "IMPeakType",
        "Processing state of ion mobility data (profile vs centroided in IM dimension)")
        .value("IM_PROFILE", OpenMS::IMPeakType::IM_PROFILE,
               "Raw/unprocessed IM data (e.g. full TIMS frame)")
        .value("IM_CENTROIDED", OpenMS::IMPeakType::IM_CENTROIDED,
               "IM data after centroiding in the IM dimension")
        .value("UNKNOWN", OpenMS::IMPeakType::UNKNOWN,
               "IM peak type not yet determined")
        .value("SIZE_OF_IMPEAKTYPE", OpenMS::IMPeakType::SIZE_OF_IMPEAKTYPE);

    // -----------------------------------------------------------------------
    // AnnotationStatistics
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AnnotationStatistics>(m, "AnnotationStatistics", "OpenMS class AnnotationStatistics")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::AnnotationStatistics &>())
        .def("__copy__", [](const OpenMS::AnnotationStatistics& self) { return OpenMS::AnnotationStatistics(self); })
        .def("__deepcopy__", [](const OpenMS::AnnotationStatistics& self, nb::dict) { return OpenMS::AnnotationStatistics(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def_rw("states", &OpenMS::AnnotationStatistics::states)
        ;

    // -----------------------------------------------------------------------
    // BinnedSpectrum
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::BinnedSpectrum>(m, "BinnedSpectrum", "This is a binned representation of a PeakSpectrum")
        .def(nb::init<>())
        .def(nb::init<OpenMS::MSSpectrum, float, bool, unsigned int, float>())
        .def(nb::init<const OpenMS::BinnedSpectrum &>())
        .def("__copy__", [](const OpenMS::BinnedSpectrum& self) { return OpenMS::BinnedSpectrum(self); })
        .def("__deepcopy__", [](const OpenMS::BinnedSpectrum& self, nb::dict) { return OpenMS::BinnedSpectrum(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getBinIntensity", [](OpenMS::BinnedSpectrum& self, double mz) { return self.getBinIntensity(mz); }, "mz"_a, "Returns the bin intensity at a given m/z position")
        .def("getBinIndex", [](const OpenMS::BinnedSpectrum& self, float mz) { return self.getBinIndex(mz); }, "mz"_a, "Returns the bin index of a given m/z position")
        .def("getBinLowerMZ", [](const OpenMS::BinnedSpectrum& self, size_t i) { return self.getBinLowerMZ(i); }, "i"_a, "Returns the lower m/z of a bin given its index")
        .def("getBinSize", [](const OpenMS::BinnedSpectrum& self) { return self.getBinSize(); }, "Returns the bin size")
        .def("getBinSpread", [](const OpenMS::BinnedSpectrum& self) { return self.getBinSpread(); }, "Returns the bin spread")
        .def("getOffset", [](const OpenMS::BinnedSpectrum& self) { return self.getOffset(); }, "Returns offset")
        .def("getPrecursors", [](OpenMS::BinnedSpectrum& self) -> std::vector<OpenMS::Precursor> & { return self.getPrecursors(); }, nb::rv_policy::reference_internal)
        .def_static("isCompatible", [](const OpenMS::BinnedSpectrum& a, const OpenMS::BinnedSpectrum& b) { return OpenMS::BinnedSpectrum::isCompatible(a, b); }, "a"_a, "b"_a)
        ;

    // -----------------------------------------------------------------------
    // ChromatogramPeak
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ChromatogramPeak>(m, "ChromatogramPeak", "A 1-dimensional raw data point or peak for chromatograms")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ChromatogramPeak &>())
        .def("__copy__", [](const OpenMS::ChromatogramPeak& self) { return OpenMS::ChromatogramPeak(self); })
        .def("__deepcopy__", [](const OpenMS::ChromatogramPeak& self, nb::dict) { return OpenMS::ChromatogramPeak(self); }, "memo"_a)
        .def(nb::init<OpenMS::DPosition<1>, float>())
        .def("getIntensity", [](const OpenMS::ChromatogramPeak& self) { return self.getIntensity(); }, "Returns the intensity")
        .def("setIntensity", [](OpenMS::ChromatogramPeak& self, float intensity) { return self.setIntensity(intensity); }, "intensity"_a, "Sets the intensity")
        .def("getRT", [](const OpenMS::ChromatogramPeak& self) { return self.getRT(); }, "Returns the retention time")
        .def("setRT", [](OpenMS::ChromatogramPeak& self, double rt) { return self.setRT(rt); }, "rt"_a, "Sets retention time")
        .def("getPos", [](const OpenMS::ChromatogramPeak& self) { return self.getPos(); }, "Alias for getRT()")
        .def("setPos", [](OpenMS::ChromatogramPeak& self, double pos) { return self.setPos(pos); }, "pos"_a, "Alias for setRT()")
        .def("getPosition", [](OpenMS::ChromatogramPeak& self) -> OpenMS::DPosition<1> & { return self.getPosition(); }, nb::rv_policy::reference_internal, "Returns the position (RT)")
        .def("setPosition", [](OpenMS::ChromatogramPeak& self, const OpenMS::DPosition<1>& position) { return self.setPosition(position); }, "position"_a, "Sets the position (RT)")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("__hash__", [](const OpenMS::ChromatogramPeak& self) {
            // Content-based hash using rt and intensity
            size_t h1 = std::hash<double>{}(self.getRT());
            size_t h2 = std::hash<float>{}(self.getIntensity());
            return h1 ^ (h2 << 1);
        })
        .def("__repr__", [](const OpenMS::ChromatogramPeak& self) {
            std::ostringstream os;
            os << "ChromatogramPeak(rt=" << self.getRT() << ", intensity=" << self.getIntensity() << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::ChromatogramPeak& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;

    // -----------------------------------------------------------------------
    // ChromatogramRangeManager
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ChromatogramRangeManager>(m, "ChromatogramRangeManager", 
        R"doc(
Range manager for chromatograms
This class manages retention time, m/z, and intensity ranges for multiple chromatograms.
It extends the basic RangeManager to provide specialized functionality for chromatogram data.
The template parameters for the base RangeManager are ordered differently than in SpectrumRangeManager:
- RangeRT (retention time) is the first parameter, as it's the primary dimension for chromatograms
- RangeIntensity is the second parameter
- RangeMZ is the third parameter
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ChromatogramRangeManager &>())
        .def("__copy__", [](const OpenMS::ChromatogramRangeManager& self) { return OpenMS::ChromatogramRangeManager(self); })
        .def("__deepcopy__", [](const OpenMS::ChromatogramRangeManager& self, nb::dict) { return OpenMS::ChromatogramRangeManager(self); }, "memo"_a)
        .def("clearRanges", [](OpenMS::ChromatogramRangeManager& self) { return self.clearRanges(); }, "Clear all ranges")
        .def("getMinRT", [](const OpenMS::ChromatogramRangeManager& self) { return self.getMinRT(); }, "Get the minimum RT value")
        .def("getMaxRT", [](const OpenMS::ChromatogramRangeManager& self) { return self.getMaxRT(); }, "Get the maximum RT value")
        .def("getMinMZ", [](const OpenMS::ChromatogramRangeManager& self) { return self.getMinMZ(); }, "Get the minimum m/z value")
        .def("getMaxMZ", [](const OpenMS::ChromatogramRangeManager& self) { return self.getMaxMZ(); }, "Get the maximum m/z value")
        .def("getMinIntensity", [](const OpenMS::ChromatogramRangeManager& self) { return self.getMinIntensity(); }, "Get the minimum intensity value")
        .def("getMaxIntensity", [](const OpenMS::ChromatogramRangeManager& self) { return self.getMaxIntensity(); }, "Get the maximum intensity value")
        ;

    // -----------------------------------------------------------------------
    // ChromatogramTools
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ChromatogramTools>(m, "ChromatogramTools", "Conversion class to convert chromatograms")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ChromatogramTools &>())
        .def("__copy__", [](const OpenMS::ChromatogramTools& self) { return OpenMS::ChromatogramTools(self); })
        .def("__deepcopy__", [](const OpenMS::ChromatogramTools& self, nb::dict) { return OpenMS::ChromatogramTools(self); }, "memo"_a)
        .def("convertChromatogramsToSpectra", [](OpenMS::ChromatogramTools& self, OpenMS::MSExperiment& exp) { self.convertChromatogramsToSpectra(exp); }, "exp"_a, "Converts the chromatograms in the experiment to a list of spectra with instrument settings")
        .def("convertSpectraToChromatograms", [](OpenMS::ChromatogramTools& self, OpenMS::MSExperiment& exp, bool remove_spectra, bool force_conversion) { self.convertSpectraToChromatograms(exp, remove_spectra, force_conversion); }, "exp"_a, "remove_spectra"_a = false, "force_conversion"_a = false, "Converts e.g. SRM spectra to chromatograms")
        ;

    // -----------------------------------------------------------------------
    // DocumentIdentifier
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DocumentIdentifier>(m, "DocumentIdentifier", "Manage source document information")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DocumentIdentifier &>())
        .def("__copy__", [](const OpenMS::DocumentIdentifier& self) { return OpenMS::DocumentIdentifier(self); })
        .def("__deepcopy__", [](const OpenMS::DocumentIdentifier& self, nb::dict) { return OpenMS::DocumentIdentifier(self); }, "memo"_a)
        .def("setIdentifier", [](OpenMS::DocumentIdentifier& self, const OpenMS::String& id) { return self.setIdentifier(id); }, "id"_a, "Sets document identifier (e.g. an LSID)")
        .def("getIdentifier", [](const OpenMS::DocumentIdentifier& self) { return self.getIdentifier(); }, "Retrieve document identifier (e.g. an LSID)")
        .def("setLoadedFilePath", [](OpenMS::DocumentIdentifier& self, const OpenMS::String& file_name) { return self.setLoadedFilePath(file_name); }, "file_name"_a, "Sets the file_name according to absolute path of the file loaded, preferably done whilst loading")
        .def("getLoadedFilePath", [](const OpenMS::DocumentIdentifier& self) { return self.getLoadedFilePath(); }, "Returns the file_name which is the absolute path to the file loaded")
        .def("setLoadedFileType", [](OpenMS::DocumentIdentifier& self, const OpenMS::String& file_name) { return self.setLoadedFileType(file_name); }, "file_name"_a, "Sets the file_type according to the type of the file loaded from, preferably done whilst loading")
        .def("getLoadedFileType", [](const OpenMS::DocumentIdentifier& self) -> const OpenMS::FileTypes::Type & { return self.getLoadedFileType(); }, nb::rv_policy::reference_internal, "Returns the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded")
        ;

    // -----------------------------------------------------------------------
    // MRMTransitionGroupCP
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>>(m, "MRMTransitionGroupCP", "OpenMS class MRMTransitionGroupCP")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>&>())
        .def("__copy__", [](const OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) { return OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>(self); })
        .def("__deepcopy__", [](const OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, nb::dict) { return OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>(self); }, "memo"_a)
        .def("size", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) { return self.size(); })
        .def("getTransitionGroupID", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) { return self.getTransitionGroupID(); })
        .def("setTransitionGroupID", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::String tr_gr_id) { self.setTransitionGroupID(tr_gr_id); }, "tr_gr_id"_a)
        .def("getTransitions", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) { return self.getTransitions(); })
        .def("getTransitionsMuteable", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) -> std::vector<OpenMS::ReactionMonitoringTransition>& { return self.getTransitionsMuteable(); }, nb::rv_policy::reference_internal)
        .def("addTransition", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::ReactionMonitoringTransition transition, OpenMS::String key) { self.addTransition(transition, key); }, "transition"_a, "key"_a)
        .def("getTransition", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::String key) { return self.getTransition(key); }, "key"_a)
        .def("hasTransition", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::String key) { return self.hasTransition(key); }, "key"_a)
        .def("getChromatograms", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) -> std::vector<OpenMS::MSChromatogram>& { return self.getChromatograms(); }, nb::rv_policy::reference_internal)
        .def("addChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::MSChromatogram chromatogram, OpenMS::String key) { self.addChromatogram(chromatogram, key); }, "chromatogram"_a, "key"_a)
        .def("getChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::String key) { return self.getChromatogram(key); }, "key"_a)
        .def("hasChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::String key) { return self.hasChromatogram(key); }, "key"_a)
        .def("getPrecursorChromatograms", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) -> std::vector<OpenMS::MSChromatogram>& { return self.getPrecursorChromatograms(); }, nb::rv_policy::reference_internal)
        .def("addPrecursorChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::MSChromatogram chromatogram, OpenMS::String key) { self.addPrecursorChromatogram(chromatogram, key); }, "chromatogram"_a, "key"_a)
        .def("getPrecursorChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::String key) { return self.getPrecursorChromatogram(key); }, "key"_a)
        .def("hasPrecursorChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::String key) { return self.hasPrecursorChromatogram(key); }, "key"_a)
        .def("getFeatures", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) { return self.getFeatures(); })
        .def("getFeaturesMuteable", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) -> std::vector<OpenMS::MRMFeature>& { return self.getFeaturesMuteable(); }, nb::rv_policy::reference_internal)
        .def("addFeature", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, OpenMS::MRMFeature feature) { self.addFeature(feature); }, "feature"_a)
        .def("getBestFeature", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) { return self.getBestFeature(); })
        .def("getLibraryIntensity", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, std::vector<double> result) { self.getLibraryIntensity(result); }, "result"_a)
        .def("subset", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self, std::vector<std::string> tr_ids) { return self.subset(tr_ids); }, "tr_ids"_a)
        .def("isInternallyConsistent", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) { return self.isInternallyConsistent(); })
        .def("chromatogramIdsMatch", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) { return self.chromatogramIdsMatch(); })
        .def("__len__", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) { return self.size(); })
        .def("__repr__", [](const OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) {
            return "MRMTransitionGroupCP(id='" + self.getTransitionGroupID() + "', size=" + std::to_string(self.size()) + ")";
        })
        .def("__str__", [](const OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenMS::ReactionMonitoringTransition>& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;

    // -----------------------------------------------------------------------
    // LightMRMTransitionGroupCP
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>>(m, "LightMRMTransitionGroupCP", "OpenMS class LightMRMTransitionGroupCP")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>&>())
        .def("__copy__", [](const OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) { return OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>(self); })
        .def("__deepcopy__", [](const OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, nb::dict) { return OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>(self); }, "memo"_a)
        .def("size", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) { return self.size(); })
        .def("getTransitionGroupID", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) { return self.getTransitionGroupID(); })
        .def("setTransitionGroupID", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenMS::String tr_gr_id) { self.setTransitionGroupID(tr_gr_id); }, "tr_gr_id"_a)
        .def("getTransitions", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) { return self.getTransitions(); })
        .def("getTransitionsMuteable", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) -> std::vector<OpenSwath::LightTransition>& { return self.getTransitionsMuteable(); }, nb::rv_policy::reference_internal)
        .def("addTransition", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenSwath::LightTransition transition, OpenMS::String key) { self.addTransition(transition, key); }, "transition"_a, "key"_a)
        .def("getTransition", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenMS::String key) { return self.getTransition(key); }, "key"_a)
        .def("hasTransition", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenMS::String key) { return self.hasTransition(key); }, "key"_a)
        .def("getChromatograms", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) -> std::vector<OpenMS::MSChromatogram>& { return self.getChromatograms(); }, nb::rv_policy::reference_internal)
        .def("addChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenMS::MSChromatogram chromatogram, OpenMS::String key) { self.addChromatogram(chromatogram, key); }, "chromatogram"_a, "key"_a)
        .def("getChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenMS::String key) { return self.getChromatogram(key); }, "key"_a)
        .def("hasChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenMS::String key) { return self.hasChromatogram(key); }, "key"_a)
        .def("getPrecursorChromatograms", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) -> std::vector<OpenMS::MSChromatogram>& { return self.getPrecursorChromatograms(); }, nb::rv_policy::reference_internal)
        .def("addPrecursorChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenMS::MSChromatogram chromatogram, OpenMS::String key) { self.addPrecursorChromatogram(chromatogram, key); }, "chromatogram"_a, "key"_a)
        .def("getPrecursorChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenMS::String key) { return self.getPrecursorChromatogram(key); }, "key"_a)
        .def("hasPrecursorChromatogram", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenMS::String key) { return self.hasPrecursorChromatogram(key); }, "key"_a)
        .def("getFeatures", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) { return self.getFeatures(); })
        .def("getFeaturesMuteable", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) -> std::vector<OpenMS::MRMFeature>& { return self.getFeaturesMuteable(); }, nb::rv_policy::reference_internal)
        .def("addFeature", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, OpenMS::MRMFeature feature) { self.addFeature(feature); }, "feature"_a)
        .def("getBestFeature", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) { return self.getBestFeature(); })
        .def("getLibraryIntensity", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, std::vector<double> result) { self.getLibraryIntensity(result); }, "result"_a)
        .def("subset", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self, std::vector<std::string> tr_ids) { return self.subset(tr_ids); }, "tr_ids"_a)
        .def("isInternallyConsistent", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) { return self.isInternallyConsistent(); })
        .def("chromatogramIdsMatch", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) { return self.chromatogramIdsMatch(); })
        .def("__len__", [](OpenMS::MRMTransitionGroup<OpenMS::MSChromatogram, OpenSwath::LightTransition>& self) { return self.size(); })
        ;


    // -----------------------------------------------------------------------
    // MapConversion
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MapConversion>(m, "MapConversion", "OpenMS class MapConversion")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MapConversion &>())
        .def("__copy__", [](const OpenMS::MapConversion& self) { return OpenMS::MapConversion(self); })
        .def("__deepcopy__", [](const OpenMS::MapConversion& self, nb::dict) { return OpenMS::MapConversion(self); }, "memo"_a)
        .def_static("convert", [](size_t input_map_index, OpenMS::MSExperiment& input_map, OpenMS::ConsensusMap& output_map, size_t n) { return OpenMS::MapConversion::convert(input_map_index, input_map, output_map, n); }, "input_map_index"_a, "input_map"_a, "output_map"_a, "n"_a)
        .def_static("convert", [](const OpenMS::ConsensusMap& input_map, bool keep_uids, OpenMS::FeatureMap& output_map) { return OpenMS::MapConversion::convert(input_map, keep_uids, output_map); }, "input_map"_a, "keep_uids"_a, "output_map"_a)
        .def_static("convert", [](size_t input_map_index, const OpenMS::FeatureMap& input_map, OpenMS::ConsensusMap& output_map, size_t n) { return OpenMS::MapConversion::convert(input_map_index, input_map, output_map, n); }, "input_map_index"_a, "input_map"_a, "output_map"_a, "n"_a)
        ;

    // -----------------------------------------------------------------------
    // MassTrace
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MassTrace>(m, "MassTrace", "OpenMS class MassTrace")
        .def(nb::init<>())
        .def(nb::init<std::list<OpenMS::Peak2D>>())
        .def(nb::init<std::vector<OpenMS::Peak2D>>())
        .def(nb::init<const OpenMS::MassTrace &>())
        .def("__copy__", [](const OpenMS::MassTrace& self) { return OpenMS::MassTrace(self); })
        .def("__deepcopy__", [](const OpenMS::MassTrace& self, nb::dict) { return OpenMS::MassTrace(self); }, "memo"_a)
        .def("getSize", [](const OpenMS::MassTrace& self) { return self.getSize(); })
        .def("getLabel", [](const OpenMS::MassTrace& self) { return self.getLabel(); })
        .def("setLabel", [](OpenMS::MassTrace& self, const OpenMS::String& label) { self.setLabel(label); }, "label"_a)
        .def("getCentroidMZ", [](const OpenMS::MassTrace& self) { return self.getCentroidMZ(); })
        .def("getCentroidRT", [](const OpenMS::MassTrace& self) { return self.getCentroidRT(); })
        .def("getCentroidSD", [](const OpenMS::MassTrace& self) { return self.getCentroidSD(); })
        .def("setCentroidSD", [](OpenMS::MassTrace& self, double sd) { self.setCentroidSD(sd); }, "sd"_a)
        .def("getFWHM", [](const OpenMS::MassTrace& self) { return self.getFWHM(); })
        .def("getTraceLength", [](const OpenMS::MassTrace& self) { return self.getTraceLength(); })
        .def("getFWHMborders", [](const OpenMS::MassTrace& self) { return self.getFWHMborders(); })
        .def("getSmoothedIntensities", [](const OpenMS::MassTrace& self) { return self.getSmoothedIntensities(); })
        .def("setSmoothedIntensities", [](OpenMS::MassTrace& self, const std::vector<double>& ints) { self.setSmoothedIntensities(ints); }, "intensities"_a)
        .def("getAverageMS1CycleTime", [](const OpenMS::MassTrace& self) { return self.getAverageMS1CycleTime(); })
        .def("computeSmoothedPeakArea", [](const OpenMS::MassTrace& self) { return self.computeSmoothedPeakArea(); })
        .def("computePeakArea", [](const OpenMS::MassTrace& self) { return self.computePeakArea(); })
        .def("findMaxByIntPeak", [](const OpenMS::MassTrace& self, bool use_smoothed) { return self.findMaxByIntPeak(use_smoothed); }, "use_smoothed_ints"_a = false)
        .def("estimateFWHM", [](OpenMS::MassTrace& self, bool use_smoothed) { return self.estimateFWHM(use_smoothed); }, "use_smoothed_ints"_a = false)
        .def("computeFwhmArea", [](const OpenMS::MassTrace& self) { return self.computeFwhmArea(); })
        .def("computeFwhmAreaSmooth", [](const OpenMS::MassTrace& self) { return self.computeFwhmAreaSmooth(); })
        .def("getIntensity", [](const OpenMS::MassTrace& self, bool smoothed) { return self.getIntensity(smoothed); }, "smoothed"_a)
        .def("getMaxIntensity", [](const OpenMS::MassTrace& self, bool smoothed) { return self.getMaxIntensity(smoothed); }, "smoothed"_a)
        .def("getConvexhull", [](const OpenMS::MassTrace& self) { return self.getConvexhull(); })
        .def("updateSmoothedMaxRT", [](OpenMS::MassTrace& self) { self.updateSmoothedMaxRT(); })
        .def("updateWeightedMeanRT", [](OpenMS::MassTrace& self) { self.updateWeightedMeanRT(); })
        .def("updateSmoothedWeightedMeanRT", [](OpenMS::MassTrace& self) { self.updateSmoothedWeightedMeanRT(); })
        .def("updateMedianRT", [](OpenMS::MassTrace& self) { self.updateMedianRT(); })
        .def("updateMedianMZ", [](OpenMS::MassTrace& self) { self.updateMedianMZ(); })
        .def("updateMeanMZ", [](OpenMS::MassTrace& self) { self.updateMeanMZ(); })
        .def("updateWeightedMeanMZ", [](OpenMS::MassTrace& self) { self.updateWeightedMeanMZ(); })
        .def("updateWeightedMZsd", [](OpenMS::MassTrace& self) { self.updateWeightedMZsd(); })
        .def_rw("fwhm_mz_avg", &OpenMS::MassTrace::fwhm_mz_avg)
        .def_rw("fwhm_im_avg", &OpenMS::MassTrace::fwhm_im_avg)
        .def("computeIntensitySum", [](const OpenMS::MassTrace& self) { return self.computeIntensitySum(); }, "Sum all peak intensities in the mass trace")
        .def("getQuantMethod", [](const OpenMS::MassTrace& self) { return self.getQuantMethod(); }, "Returns the quantification method")
        .def("setQuantMethod", [](OpenMS::MassTrace& self, OpenMS::MassTrace::MT_QUANTMETHOD method) { self.setQuantMethod(method); }, "method"_a, "Sets the quantification method")
        .def("getAvgMZ", [](const OpenMS::MassTrace& self) { return self.getCentroidMZ(); }, "Returns the centroid m/z (alias for getCentroidMZ)")
        ;

    // -----------------------------------------------------------------------
    // MT_QUANTMETHOD (MassTrace::MT_QUANTMETHOD)
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::MassTrace::MT_QUANTMETHOD>(m, "MT_QUANTMETHOD",
        "Quantification method for mass traces")
        .value("MT_QUANT_AREA", OpenMS::MassTrace::MT_QUANT_AREA)
        .value("MT_QUANT_MEDIAN", OpenMS::MassTrace::MT_QUANT_MEDIAN)
        .value("MT_QUANT_HEIGHT", OpenMS::MassTrace::MT_QUANT_HEIGHT)
        .value("SIZE_OF_MT_QUANTMETHOD", OpenMS::MassTrace::SIZE_OF_MT_QUANTMETHOD)

        ;

    // -----------------------------------------------------------------------
    // MetaInfoInterface
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MetaInfoInterface>(m, "MetaInfoInterface", 
        R"doc(
Interface for classes that can store arbitrary meta information
(Type-Name-Value tuples).
MetaInfoInterface is a base class for all classes that use one MetaInfo
object as member.  If you want to add meta information to a class, let it
publicly inherit the MetaInfoInterface.  Meta information is an array of
Type-Name-Value tuples.
Usage:
.. code-block:: python
k = []
exp.getKeys(k) # explore available key-value pairs
exp.getMetaValue("someMetaName")
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MetaInfoInterface &>())
        .def("__copy__", [](const OpenMS::MetaInfoInterface& self) { return OpenMS::MetaInfoInterface(self); })
        .def("__deepcopy__", [](const OpenMS::MetaInfoInterface& self, nb::dict) { return OpenMS::MetaInfoInterface(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getMetaValue", [](const OpenMS::MetaInfoInterface& self, const OpenMS::String& name) { return self.getMetaValue(name); }, "name"_a, "Returns the value corresponding to a string, or DataValue::EMPTY if not found")
        .def("getMetaValue", [](const OpenMS::MetaInfoInterface& self, const OpenMS::String& name, const OpenMS::DataValue& default_value) { return self.getMetaValue(name, default_value); }, "name"_a, "default_value"_a, "Returns the value corresponding to a string, or DataValue::EMPTY if not found")
        .def("getMetaValue", [](const OpenMS::MetaInfoInterface& self, unsigned int index) { return self.getMetaValue(index); }, "index"_a, "Returns the value corresponding to a string, or DataValue::EMPTY if not found")
        .def("getMetaValue", [](const OpenMS::MetaInfoInterface& self, unsigned int index, const OpenMS::DataValue& default_value) { return self.getMetaValue(index, default_value); }, "index"_a, "default_value"_a, "Returns the value corresponding to a string, or DataValue::EMPTY if not found")
        .def("metaValueExists", [](const OpenMS::MetaInfoInterface& self, const OpenMS::String& name) { return self.metaValueExists(name); }, "name"_a, "Returns whether an entry with the given name exists")
        .def("metaValueExists", [](const OpenMS::MetaInfoInterface& self, unsigned int index) { return self.metaValueExists(index); }, "index"_a, "Returns whether an entry with the given name exists")
        .def("setMetaValue", [](OpenMS::MetaInfoInterface& self, const OpenMS::String& name, const OpenMS::DataValue& value) { return self.setMetaValue(name, value); }, "name"_a, "value"_a, "Sets the DataValue corresponding to a name")
        .def("setMetaValue", [](OpenMS::MetaInfoInterface& self, unsigned int index, const OpenMS::DataValue& value) { return self.setMetaValue(index, value); }, "index"_a, "value"_a, "Sets the DataValue corresponding to a name")
        .def("removeMetaValue", [](OpenMS::MetaInfoInterface& self, const OpenMS::String& name) { return self.removeMetaValue(name); }, "name"_a, "Removes the DataValue corresponding to `name` if it exists")
        .def("removeMetaValue", [](OpenMS::MetaInfoInterface& self, unsigned int index) { return self.removeMetaValue(index); }, "index"_a, "Removes the DataValue corresponding to `name` if it exists")
        .def_static("metaRegistry", []() { return OpenMS::MetaInfoInterface::metaRegistry(); }, "Returns a reference to the MetaInfoRegistry")
        
        .def("getKeys", [](const OpenMS::MetaInfoInterface& self, nb::list py_keys) {
            std::vector<OpenMS::String> keys;
            self.getKeys(keys);
            py_keys.attr("clear")();
            for (const auto& k : keys) {
                py_keys.append(nb::str(k.c_str()));
            }
        }, "keys"_a, "Fills the given list with all meta value keys")
        .def("isMetaEmpty", [](const OpenMS::MetaInfoInterface& self) { return self.isMetaEmpty(); }, "Returns if the MetaInfo is empty")
        .def("clearMetaInfo", [](OpenMS::MetaInfoInterface& self) { return self.clearMetaInfo(); }, "Removes all meta values")
        .def("getMetaValues", [](const OpenMS::MetaInfoInterface& self) {
            nb::dict result;
            std::vector<OpenMS::String> keys;
            self.getKeys(keys);
            for (const auto& key : keys)
            {
                result[nb::str(key.c_str())] = nb::cast(self.getMetaValue(key));
            }
            return result;
        }, "Returns all meta values as a dictionary of {name: value}")
        .def("setMetaValues", [](OpenMS::MetaInfoInterface& self, nb::dict values) {
            for (auto item : values)
            {
                std::string key = nb::cast<std::string>(item.first);
                OpenMS::DataValue val = nb::cast<OpenMS::DataValue>(item.second);
                self.setMetaValue(key, val);
            }
        }, "values"_a, "Sets multiple meta values from a dictionary of {name: value}")
        .def("__hash__", [](const OpenMS::MetaInfoInterface& self) { return std::hash<OpenMS::MetaInfoInterface>{}(self); })
        ;

    // -----------------------------------------------------------------------
    // Acquisition
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Acquisition, OpenMS::MetaInfoInterface>(m, "Acquisition", 
        R"doc(
Information about one raw data spectrum that was combined with
several other raw data spectra
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Acquisition &>())
        .def("__copy__", [](const OpenMS::Acquisition& self) { return OpenMS::Acquisition(self); })
        .def("__deepcopy__", [](const OpenMS::Acquisition& self, nb::dict) { return OpenMS::Acquisition(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getIdentifier", [](const OpenMS::Acquisition& self) { return self.getIdentifier(); })
        .def("setIdentifier", [](OpenMS::Acquisition& self, const OpenMS::String& identifier) { return self.setIdentifier(identifier); }, "identifier"_a)
        
        ;

    // -----------------------------------------------------------------------
    // AcquisitionInfo
    // -----------------------------------------------------------------------
    auto acquisitioninfo_class = nb::class_<OpenMS::AcquisitionInfo>(m, "AcquisitionInfo", 
        R"doc(
Description of the combination of raw data to a single spectrum
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::AcquisitionInfo &>())
        .def("__copy__", [](const OpenMS::AcquisitionInfo& self) { return OpenMS::AcquisitionInfo(self); })
        .def("__deepcopy__", [](const OpenMS::AcquisitionInfo& self, nb::dict) { return OpenMS::AcquisitionInfo(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getMethodOfCombination", [](const OpenMS::AcquisitionInfo& self) { return self.getMethodOfCombination(); }, "Returns the method of combination")
        .def("setMethodOfCombination", [](OpenMS::AcquisitionInfo& self, const OpenMS::String& method_of_combination) { return self.setMethodOfCombination(method_of_combination); }, "method_of_combination"_a, "Sets the method of combination")
        

        .def("size", [](const OpenMS::AcquisitionInfo& self) -> size_t {
            return self.size();
        }, "Returns the number of Acquisition objects")

        .def("push_back", [](OpenMS::AcquisitionInfo& self, const OpenMS::Acquisition& acq) {
            self.push_back(acq);
        }, "acq"_a, "Append an Acquisition object")

        .def("resize", [](OpenMS::AcquisitionInfo& self, size_t n) {
            self.resize(n);
        }, "n"_a, "Resize the AcquisitionInfo")

        .def("__getitem__", [](OpenMS::AcquisitionInfo& self, size_t i) -> OpenMS::Acquisition& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, "i"_a, nb::rv_policy::reference_internal)

        .def("__setitem__", [](OpenMS::AcquisitionInfo& self, size_t i, const OpenMS::Acquisition& acq) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = acq;
        }, "i"_a, "acq"_a)
        ;
    def_MetaInfoInterface<OpenMS::AcquisitionInfo>(acquisitioninfo_class);

    // -----------------------------------------------------------------------
    // CVTermList
    // -----------------------------------------------------------------------
    auto cvtermlist_class = nb::class_<OpenMS::CVTermList>(m, "CVTermList", "OpenMS class CVTermList")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::CVTermList &>())
        .def("__copy__", [](const OpenMS::CVTermList& self) { return OpenMS::CVTermList(self); })
        .def("__deepcopy__", [](const OpenMS::CVTermList& self, nb::dict) { return OpenMS::CVTermList(self); }, "memo"_a)
        .def("replaceCVTerm", [](OpenMS::CVTermList& self, const OpenMS::CVTerm& cv_term) { return self.replaceCVTerm(cv_term); }, "cv_term"_a, "Replaces the specified CV term")
        .def("replaceCVTerms", [](OpenMS::CVTermList& self, const std::map<OpenMS::String, std::vector<OpenMS::CVTerm>>& cv_term_map) { return self.replaceCVTerms(cv_term_map); }, "cv_term_map"_a)
        .def("consumeCVTerms", [](OpenMS::CVTermList& self, const std::map<OpenMS::String, std::vector<OpenMS::CVTerm>>& cv_term_map) { return self.consumeCVTerms(cv_term_map); }, "cv_term_map"_a, "Merges the given map into the member map, no duplicate checking")
        .def("getCVTerms", [](const OpenMS::CVTermList& self) -> const std::map<OpenMS::String, std::vector<OpenMS::CVTerm>> & { return self.getCVTerms(); }, nb::rv_policy::reference_internal, "Returns the accession string of the term")
        .def("addCVTerm", [](OpenMS::CVTermList& self, const OpenMS::CVTerm& term) { return self.addCVTerm(term); }, "term"_a, "Adds a CV term")
        .def("setCVTerms", [](OpenMS::CVTermList& self, const std::vector<OpenMS::CVTerm>& terms) { return self.setCVTerms(terms); }, "terms"_a, "Sets the CV terms from a vector")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("hasCVTerm", [](const OpenMS::CVTermList& self, const OpenMS::String& accession) { return self.hasCVTerm(accession); }, "accession"_a)
        .def("empty", [](const OpenMS::CVTermList& self) { return self.empty(); })
        
        .def("__hash__", [](const OpenMS::CVTermList& self) { return std::hash<OpenMS::CVTermList>{}(self); })
        ;
    def_MetaInfoInterface<OpenMS::CVTermList>(cvtermlist_class);

    // -----------------------------------------------------------------------
    // CVTermListInterface
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CVTermListInterface, OpenMS::MetaInfoInterface>(m, "CVTermListInterface", 
        R"doc(
Interface to the controlled vocabulary term list
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::CVTermListInterface &>())
        .def("__copy__", [](const OpenMS::CVTermListInterface& self) { return OpenMS::CVTermListInterface(self); })
        .def("__deepcopy__", [](const OpenMS::CVTermListInterface& self, nb::dict) { return OpenMS::CVTermListInterface(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("replaceCVTerms", [](OpenMS::CVTermListInterface& self, std::map<OpenMS::String, std::vector<OpenMS::CVTerm>>& cv_terms) { return self.replaceCVTerms(cv_terms); }, "cv_terms"_a)
        .def("replaceCVTerm", [](OpenMS::CVTermListInterface& self, const OpenMS::CVTerm& cv_term) { return self.replaceCVTerm(cv_term); }, "cv_term"_a)
        .def("replaceCVTerms", [](OpenMS::CVTermListInterface& self, const std::map<OpenMS::String, std::vector<OpenMS::CVTerm>>& cv_term_map) { return self.replaceCVTerms(cv_term_map); }, "cv_term_map"_a)
        .def("consumeCVTerms", [](OpenMS::CVTermListInterface& self, const std::map<OpenMS::String, std::vector<OpenMS::CVTerm>>& cv_term_map) { return self.consumeCVTerms(cv_term_map); }, "cv_term_map"_a, "Merges the given map into the member map, no duplicate checking")
        .def("getCVTerms", [](const OpenMS::CVTermListInterface& self) -> const std::map<OpenMS::String, std::vector<OpenMS::CVTerm>> & { return self.getCVTerms(); }, nb::rv_policy::reference_internal)
        .def("addCVTerm", [](OpenMS::CVTermListInterface& self, const OpenMS::CVTerm& term) { return self.addCVTerm(term); }, "term"_a, "Adds a CV term")
        .def("setCVTerms", [](OpenMS::CVTermListInterface& self, const std::vector<OpenMS::CVTerm>& terms) { return self.setCVTerms(terms); }, "terms"_a, "Sets the CV terms from a vector")
        .def("hasCVTerm", [](const OpenMS::CVTermListInterface& self, const OpenMS::String& accession) { return self.hasCVTerm(accession); }, "accession"_a, "Checks whether the term has a value")
        .def("empty", [](const OpenMS::CVTermListInterface& self) { return self.empty(); })
        
        ;

    // -----------------------------------------------------------------------
    // ChromatogramSettings
    // -----------------------------------------------------------------------
    auto chromatogramsettings_class = nb::class_<OpenMS::ChromatogramSettings>(m, "ChromatogramSettings", 
        R"doc(
MetaInfoInterface

Description of the chromatogram settings, provides meta-information
about a single chromatogram.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ChromatogramSettings &>())
        .def("__copy__", [](const OpenMS::ChromatogramSettings& self) { return OpenMS::ChromatogramSettings(self); })
        .def("__deepcopy__", [](const OpenMS::ChromatogramSettings& self, nb::dict) { return OpenMS::ChromatogramSettings(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getNativeID", [](const OpenMS::ChromatogramSettings& self) { return self.getNativeID(); }, "Returns the native identifier for the spectrum, used by the acquisition software.")
        .def("setNativeID", [](OpenMS::ChromatogramSettings& self, const OpenMS::String& native_id) { return self.setNativeID(native_id); }, "native_id"_a, "Sets the native identifier for the spectrum, used by the acquisition software.")
        .def("getComment", [](const OpenMS::ChromatogramSettings& self) { return self.getComment(); }, "Returns the free-text comment")
        .def("setComment", [](OpenMS::ChromatogramSettings& self, const OpenMS::String& comment) { return self.setComment(comment); }, "comment"_a, "Sets the free-text comment")
        .def("getInstrumentSettings", [](OpenMS::ChromatogramSettings& self) -> OpenMS::InstrumentSettings & { return self.getInstrumentSettings(); }, nb::rv_policy::reference_internal, "Returns the instrument settings of the current spectrum")
        .def("setInstrumentSettings", [](OpenMS::ChromatogramSettings& self, const OpenMS::InstrumentSettings& instrument_settings) { return self.setInstrumentSettings(instrument_settings); }, "instrument_settings"_a, "Sets the instrument settings of the current spectrum")
        .def("getAcquisitionInfo", [](OpenMS::ChromatogramSettings& self) -> OpenMS::AcquisitionInfo & { return self.getAcquisitionInfo(); }, nb::rv_policy::reference_internal, "Returns the acquisition info")
        .def("setAcquisitionInfo", [](OpenMS::ChromatogramSettings& self, const OpenMS::AcquisitionInfo& acquisition_info) { return self.setAcquisitionInfo(acquisition_info); }, "acquisition_info"_a, "Sets the acquisition info")
        .def("getSourceFile", [](OpenMS::ChromatogramSettings& self) -> OpenMS::SourceFile & { return self.getSourceFile(); }, nb::rv_policy::reference_internal, "Returns the source file")
        .def("setSourceFile", [](OpenMS::ChromatogramSettings& self, const OpenMS::SourceFile& source_file) { return self.setSourceFile(source_file); }, "source_file"_a, "Sets the source file")
        .def("getPrecursor", [](OpenMS::ChromatogramSettings& self) -> OpenMS::Precursor & { return self.getPrecursor(); }, nb::rv_policy::reference_internal, "Returns the precursors")
        .def("setPrecursor", [](OpenMS::ChromatogramSettings& self, const OpenMS::Precursor& precursor) { return self.setPrecursor(precursor); }, "precursor"_a, "Sets the precursors")
        .def("getProduct", [](OpenMS::ChromatogramSettings& self) -> OpenMS::Product & { return self.getProduct(); }, nb::rv_policy::reference_internal, "Returns the product ion")
        .def("setProduct", [](OpenMS::ChromatogramSettings& self, const OpenMS::Product& product) { return self.setProduct(product); }, "product"_a, "Sets the product ion")
        .def("getChromatogramType", [](const OpenMS::ChromatogramSettings& self) { return self.getChromatogramType(); }, "Get the chromatogram type")
        .def("setChromatogramType", [](OpenMS::ChromatogramSettings& self, OpenMS::ChromatogramSettings::ChromatogramType type) { return self.setChromatogramType(type); }, "type"_a, "Sets the chromatogram type")
        .def("setDataProcessing", [](OpenMS::ChromatogramSettings& self, const std::vector<std::shared_ptr<OpenMS::DataProcessing>>& data_processing) { return self.setDataProcessing(data_processing); }, "data_processing"_a, "Sets the description of the applied processing")
        .def("getDataProcessing", [](OpenMS::ChromatogramSettings& self) -> std::vector<std::shared_ptr<OpenMS::DataProcessing>> & { return self.getDataProcessing(); }, nb::rv_policy::reference_internal, "Returns the description of the applied processing")
        
        .def("__hash__", [](const OpenMS::ChromatogramSettings& self) { return std::hash<OpenMS::ChromatogramSettings>{}(self); })
        ;
    // ChromatogramType enum nested under ChromatogramSettings
    nb::enum_<OpenMS::ChromatogramSettings::ChromatogramType>(chromatogramsettings_class, "ChromatogramType", nb::is_arithmetic())
        .value("MASS_CHROMATOGRAM", OpenMS::ChromatogramSettings::ChromatogramType::MASS_CHROMATOGRAM)
        .value("TOTAL_ION_CURRENT_CHROMATOGRAM", OpenMS::ChromatogramSettings::ChromatogramType::TOTAL_ION_CURRENT_CHROMATOGRAM)
        .value("SELECTED_ION_CURRENT_CHROMATOGRAM", OpenMS::ChromatogramSettings::ChromatogramType::SELECTED_ION_CURRENT_CHROMATOGRAM)
        .value("BASEPEAK_CHROMATOGRAM", OpenMS::ChromatogramSettings::ChromatogramType::BASEPEAK_CHROMATOGRAM)
        .value("SELECTED_ION_MONITORING_CHROMATOGRAM", OpenMS::ChromatogramSettings::ChromatogramType::SELECTED_ION_MONITORING_CHROMATOGRAM)
        .value("SELECTED_REACTION_MONITORING_CHROMATOGRAM", OpenMS::ChromatogramSettings::ChromatogramType::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
        .value("ELECTROMAGNETIC_RADIATION_CHROMATOGRAM", OpenMS::ChromatogramSettings::ChromatogramType::ELECTROMAGNETIC_RADIATION_CHROMATOGRAM)
        .value("ABSORPTION_CHROMATOGRAM", OpenMS::ChromatogramSettings::ChromatogramType::ABSORPTION_CHROMATOGRAM)
        .value("EMISSION_CHROMATOGRAM", OpenMS::ChromatogramSettings::ChromatogramType::EMISSION_CHROMATOGRAM)
        .value("SIZE_OF_CHROMATOGRAM_TYPE", OpenMS::ChromatogramSettings::ChromatogramType::SIZE_OF_CHROMATOGRAM_TYPE)
        ;
    def_MetaInfoInterface<OpenMS::ChromatogramSettings>(chromatogramsettings_class);

    // -----------------------------------------------------------------------
    // ColumnHeader
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusMap::ColumnHeader, OpenMS::MetaInfoInterface>(m, "ColumnHeader", "OpenMS class ColumnHeader")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ConsensusMap::ColumnHeader &>())
        .def("__copy__", [](const OpenMS::ConsensusMap::ColumnHeader& self) { return OpenMS::ConsensusMap::ColumnHeader(self); })
        .def("__deepcopy__", [](const OpenMS::ConsensusMap::ColumnHeader& self, nb::dict) { return OpenMS::ConsensusMap::ColumnHeader(self); }, "memo"_a)
        .def_rw("filename", &OpenMS::ConsensusMap::ColumnHeader::filename)
        .def_rw("label", &OpenMS::ConsensusMap::ColumnHeader::label)
        .def_rw("size", &OpenMS::ConsensusMap::ColumnHeader::size)
        .def_rw("unique_id", &OpenMS::ConsensusMap::ColumnHeader::unique_id)
        ;

    // -----------------------------------------------------------------------
    // Configuration
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedExperimentHelper::Configuration, OpenMS::CVTermList>(m, "Configuration", "OpenMS class Configuration")
        .def(nb::init<>())
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        
        .def_rw("contact_ref", &OpenMS::TargetedExperimentHelper::Configuration::contact_ref)
        .def_rw("instrument_ref", &OpenMS::TargetedExperimentHelper::Configuration::instrument_ref)
        .def_rw("validations", &OpenMS::TargetedExperimentHelper::Configuration::validations)
        ;

    // -----------------------------------------------------------------------
    // Contact
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedExperimentHelper::Contact, OpenMS::CVTermList>(m, "Contact", "OpenMS class Contact")
        .def(nb::init<>())
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        
        .def_rw("id", &OpenMS::TargetedExperimentHelper::Contact::id)
        ;

    // -----------------------------------------------------------------------
    // ContactPerson
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ContactPerson, OpenMS::MetaInfoInterface>(m, "ContactPerson", 
        R"doc(
Contact person information
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ContactPerson &>())
        .def("__copy__", [](const OpenMS::ContactPerson& self) { return OpenMS::ContactPerson(self); })
        .def("__deepcopy__", [](const OpenMS::ContactPerson& self, nb::dict) { return OpenMS::ContactPerson(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getFirstName", [](const OpenMS::ContactPerson& self) { return self.getFirstName(); }, "Returns the first name of the person")
        .def("setFirstName", [](OpenMS::ContactPerson& self, const OpenMS::String& name) { return self.setFirstName(name); }, "name"_a, "Sets the first name of the person")
        .def("getLastName", [](const OpenMS::ContactPerson& self) { return self.getLastName(); }, "Returns the last name of the person")
        .def("setLastName", [](OpenMS::ContactPerson& self, const OpenMS::String& name) { return self.setLastName(name); }, "name"_a, "Sets the last name of the person")
        .def("setName", [](OpenMS::ContactPerson& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the full name of the person (gets split into first and last name internally)")
        .def("getInstitution", [](const OpenMS::ContactPerson& self) { return self.getInstitution(); }, "Returns the affiliation")
        .def("setInstitution", [](OpenMS::ContactPerson& self, const OpenMS::String& institution) { return self.setInstitution(institution); }, "institution"_a, "Sets the affiliation")
        .def("getEmail", [](const OpenMS::ContactPerson& self) { return self.getEmail(); }, "Returns the email address")
        .def("setEmail", [](OpenMS::ContactPerson& self, const OpenMS::String& email) { return self.setEmail(email); }, "email"_a, "Sets the email address")
        .def("getURL", [](const OpenMS::ContactPerson& self) { return self.getURL(); }, 
            R"doc(
Returns the URL associated with the contact person (e.g., the institute webpage "https://www.higglesworth.edu/")
)doc")
        .def("setURL", [](OpenMS::ContactPerson& self, const OpenMS::String& email) { return self.setURL(email); }, "email"_a, 
            R"doc(
Sets the URL associated with the contact person (e.g., the institute webpage "https://www.higglesworth.edu/")
)doc")
        .def("getAddress", [](const OpenMS::ContactPerson& self) { return self.getAddress(); }, "Returns the address")
        .def("setAddress", [](OpenMS::ContactPerson& self, const OpenMS::String& email) { return self.setAddress(email); }, "email"_a, "Sets the address")
        .def("getContactInfo", [](const OpenMS::ContactPerson& self) { return self.getContactInfo(); }, "Returns miscellaneous info about the contact person")
        .def("setContactInfo", [](OpenMS::ContactPerson& self, const OpenMS::String& contact_info) { return self.setContactInfo(contact_info); }, "contact_info"_a, "Sets miscellaneous info about the contact person")
        
        ;

    // -----------------------------------------------------------------------
    // DataProcessing
    // -----------------------------------------------------------------------
    auto dataprocessing_class = nb::class_<OpenMS::DataProcessing, OpenMS::MetaInfoInterface>(m, "DataProcessing", 
        R"doc(
Description of the applied preprocessing steps
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DataProcessing &>())
        .def("__copy__", [](const OpenMS::DataProcessing& self) { return OpenMS::DataProcessing(self); })
        .def("__deepcopy__", [](const OpenMS::DataProcessing& self, nb::dict) { return OpenMS::DataProcessing(self); }, "memo"_a)
        .def_static("getAllNamesOfProcessingAction", []() { return OpenMS::DataProcessing::getAllNamesOfProcessingAction(); }, "Returns all processing action names known to OpenMS")
        .def_static("processingActionToString", [](OpenMS::DataProcessing::ProcessingAction action) { return OpenMS::DataProcessing::processingActionToString(action); }, "action"_a, "Convert a ProcessingAction enum to String. Throws Exception::InvalidValue if action is SIZE_OF_PROCESSINGACTION")
        .def_static("toProcessingAction", [](const OpenMS::String& name) { return OpenMS::DataProcessing::toProcessingAction(name); }, "name"_a, "Convert a string to ProcessingAction enum. Throws Exception::InvalidValue if name is not found")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getSoftware", [](OpenMS::DataProcessing& self) -> OpenMS::Software & { return self.getSoftware(); }, nb::rv_policy::reference_internal)
        .def("setSoftware", [](OpenMS::DataProcessing& self, const OpenMS::Software& software) { return self.setSoftware(software); }, "software"_a)
        .def("getProcessingActions", [](OpenMS::DataProcessing& self) -> std::set<OpenMS::DataProcessing::ProcessingAction> & { return self.getProcessingActions(); }, nb::rv_policy::reference_internal)
        .def("setProcessingActions", [](OpenMS::DataProcessing& self, const std::set<OpenMS::DataProcessing::ProcessingAction>& actions) { return self.setProcessingActions(actions); }, "actions"_a)
        .def("getCompletionTime", [](const OpenMS::DataProcessing& self) -> const OpenMS::DateTime & { return self.getCompletionTime(); }, nb::rv_policy::reference_internal)
        .def("setCompletionTime", [](OpenMS::DataProcessing& self, const OpenMS::DateTime& completion_time) { return self.setCompletionTime(completion_time); }, "completion_time"_a)
        
        ;
    // ProcessingAction enum nested under DataProcessing
    nb::enum_<OpenMS::DataProcessing::ProcessingAction>(dataprocessing_class, "ProcessingAction")
        .value("DATA_PROCESSING", OpenMS::DataProcessing::ProcessingAction::DATA_PROCESSING)
        .value("CHARGE_DECONVOLUTION", OpenMS::DataProcessing::ProcessingAction::CHARGE_DECONVOLUTION)
        .value("DEISOTOPING", OpenMS::DataProcessing::ProcessingAction::DEISOTOPING)
        .value("SMOOTHING", OpenMS::DataProcessing::ProcessingAction::SMOOTHING)
        .value("CHARGE_CALCULATION", OpenMS::DataProcessing::ProcessingAction::CHARGE_CALCULATION)
        .value("PRECURSOR_RECALCULATION", OpenMS::DataProcessing::ProcessingAction::PRECURSOR_RECALCULATION)
        .value("BASELINE_REDUCTION", OpenMS::DataProcessing::ProcessingAction::BASELINE_REDUCTION)
        .value("PEAK_PICKING", OpenMS::DataProcessing::ProcessingAction::PEAK_PICKING)
        .value("ALIGNMENT", OpenMS::DataProcessing::ProcessingAction::ALIGNMENT)
        .value("CALIBRATION", OpenMS::DataProcessing::ProcessingAction::CALIBRATION)
        .value("NORMALIZATION", OpenMS::DataProcessing::ProcessingAction::NORMALIZATION)
        .value("FILTERING", OpenMS::DataProcessing::ProcessingAction::FILTERING)
        .value("QUANTITATION", OpenMS::DataProcessing::ProcessingAction::QUANTITATION)
        .value("FEATURE_GROUPING", OpenMS::DataProcessing::ProcessingAction::FEATURE_GROUPING)
        .value("IDENTIFICATION_MAPPING", OpenMS::DataProcessing::ProcessingAction::IDENTIFICATION_MAPPING)
        .value("FORMAT_CONVERSION", OpenMS::DataProcessing::ProcessingAction::FORMAT_CONVERSION)
        .value("CONVERSION_MZDATA", OpenMS::DataProcessing::ProcessingAction::CONVERSION_MZDATA)
        .value("CONVERSION_MZML", OpenMS::DataProcessing::ProcessingAction::CONVERSION_MZML)
        .value("CONVERSION_MZXML", OpenMS::DataProcessing::ProcessingAction::CONVERSION_MZXML)
        .value("CONVERSION_DTA", OpenMS::DataProcessing::ProcessingAction::CONVERSION_DTA)
        .value("IDENTIFICATION", OpenMS::DataProcessing::ProcessingAction::IDENTIFICATION)
        .value("ION_MOBILITY_BINNING", OpenMS::DataProcessing::ProcessingAction::ION_MOBILITY_BINNING)
        .value("SIZE_OF_PROCESSINGACTION", OpenMS::DataProcessing::ProcessingAction::SIZE_OF_PROCESSINGACTION)
        .export_values();

    // -----------------------------------------------------------------------
    // ExperimentalSettings
    // -----------------------------------------------------------------------
    auto experimentalsettings_class = nb::class_<OpenMS::ExperimentalSettings>(m, "ExperimentalSettings", 
        R"doc(
DocumentIdentifier
MetaInfoInterface

Description of the experimental settings, provides meta-information
about an LC-MS/MS injection.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ExperimentalSettings &>())
        .def("__copy__", [](const OpenMS::ExperimentalSettings& self) { return OpenMS::ExperimentalSettings(self); })
        .def("__deepcopy__", [](const OpenMS::ExperimentalSettings& self, nb::dict) { return OpenMS::ExperimentalSettings(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getSample", [](OpenMS::ExperimentalSettings& self) -> OpenMS::Sample & { return self.getSample(); }, nb::rv_policy::reference_internal, "Returns a reference to the sample description")
        .def("setSample", [](OpenMS::ExperimentalSettings& self, const OpenMS::Sample& sample) { return self.setSample(sample); }, "sample"_a, "Sets the sample description")
        .def("getSourceFiles", [](OpenMS::ExperimentalSettings& self) -> std::vector<OpenMS::SourceFile> & { return self.getSourceFiles(); }, nb::rv_policy::reference_internal, "Returns a reference to the source data file")
        .def("setSourceFiles", [](OpenMS::ExperimentalSettings& self, const std::vector<OpenMS::SourceFile>& source_files) { return self.setSourceFiles(source_files); }, "source_files"_a, "Sets the source data file")
        .def("getContacts", [](OpenMS::ExperimentalSettings& self) -> std::vector<OpenMS::ContactPerson> & { return self.getContacts(); }, nb::rv_policy::reference_internal, "Returns a reference to the list of contact persons")
        .def("setContacts", [](OpenMS::ExperimentalSettings& self, const std::vector<OpenMS::ContactPerson>& contacts) { return self.setContacts(contacts); }, "contacts"_a, "Sets the list of contact persons")
        .def("getInstrument", [](OpenMS::ExperimentalSettings& self) -> OpenMS::Instrument & { return self.getInstrument(); }, nb::rv_policy::reference_internal, "Returns a reference to the MS instrument description")
        .def("setInstrument", [](OpenMS::ExperimentalSettings& self, const OpenMS::Instrument& instrument) { return self.setInstrument(instrument); }, "instrument"_a, "Sets the MS instrument description")
        .def("getHPLC", [](OpenMS::ExperimentalSettings& self) -> OpenMS::HPLC & { return self.getHPLC(); }, nb::rv_policy::reference_internal, "Returns a reference to the description of the HPLC run")
        .def("setHPLC", [](OpenMS::ExperimentalSettings& self, const OpenMS::HPLC& hplc) { return self.setHPLC(hplc); }, "hplc"_a, "Sets the description of the HPLC run")
        .def("getDateTime", [](const OpenMS::ExperimentalSettings& self) -> const OpenMS::DateTime & { return self.getDateTime(); }, nb::rv_policy::reference_internal, "Returns the date the experiment was performed")
        .def("setDateTime", [](OpenMS::ExperimentalSettings& self, const OpenMS::DateTime& date) { return self.setDateTime(date); }, "date"_a, "Sets the date the experiment was performed")
        .def("getComment", [](const OpenMS::ExperimentalSettings& self) { return self.getComment(); }, "Returns the free-text comment")
        .def("setComment", [](OpenMS::ExperimentalSettings& self, const OpenMS::String& comment) { return self.setComment(comment); }, "comment"_a, "Sets the free-text comment")
        .def("getFractionIdentifier", [](const OpenMS::ExperimentalSettings& self) { return self.getFractionIdentifier(); }, "Returns fraction identifier")
        .def("setFractionIdentifier", [](OpenMS::ExperimentalSettings& self, const OpenMS::String& fraction_identifier) { return self.setFractionIdentifier(fraction_identifier); }, "fraction_identifier"_a, "Sets the fraction identifier")
        
        ;
    def_MetaInfoInterface<OpenMS::ExperimentalSettings>(experimentalsettings_class);
    def_DocumentIdentifier<OpenMS::ExperimentalSettings>(experimentalsettings_class);

    // -----------------------------------------------------------------------
    // IncludeExcludeTarget
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IncludeExcludeTarget, OpenMS::CVTermList>(m, "IncludeExcludeTarget", "This class stores a SRM/MRM transition")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IncludeExcludeTarget &>())
        .def("__copy__", [](const OpenMS::IncludeExcludeTarget& self) { return OpenMS::IncludeExcludeTarget(self); })
        .def("__deepcopy__", [](const OpenMS::IncludeExcludeTarget& self, nb::dict) { return OpenMS::IncludeExcludeTarget(self); }, "memo"_a)
        .def("setName", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a)
        .def("getName", [](const OpenMS::IncludeExcludeTarget& self) { return self.getName(); })
        .def("setPeptideRef", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::String& peptide_ref) { return self.setPeptideRef(peptide_ref); }, "peptide_ref"_a)
        .def("getPeptideRef", [](const OpenMS::IncludeExcludeTarget& self) { return self.getPeptideRef(); })
        .def("setCompoundRef", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::String& compound_ref) { return self.setCompoundRef(compound_ref); }, "compound_ref"_a)
        .def("getCompoundRef", [](const OpenMS::IncludeExcludeTarget& self) { return self.getCompoundRef(); })
        .def("setPrecursorMZ", [](OpenMS::IncludeExcludeTarget& self, double mz) { return self.setPrecursorMZ(mz); }, "mz"_a)
        .def("getPrecursorMZ", [](const OpenMS::IncludeExcludeTarget& self) { return self.getPrecursorMZ(); })
        .def("setPrecursorCVTermList", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::CVTermList& list) { return self.setPrecursorCVTermList(list); }, "list"_a)
        .def("addPrecursorCVTerm", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::CVTerm& cv_term) { return self.addPrecursorCVTerm(cv_term); }, "cv_term"_a)
        .def("getPrecursorCVTermList", [](const OpenMS::IncludeExcludeTarget& self) -> const OpenMS::CVTermList & { return self.getPrecursorCVTermList(); }, nb::rv_policy::reference_internal)
        .def("setProductMZ", [](OpenMS::IncludeExcludeTarget& self, double mz) { return self.setProductMZ(mz); }, "mz"_a)
        .def("getProductMZ", [](const OpenMS::IncludeExcludeTarget& self) { return self.getProductMZ(); })
        .def("setProductCVTermList", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::CVTermList& list) { return self.setProductCVTermList(list); }, "list"_a)
        .def("addProductCVTerm", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::CVTerm& cv_term) { return self.addProductCVTerm(cv_term); }, "cv_term"_a)
        .def("getProductCVTermList", [](const OpenMS::IncludeExcludeTarget& self) -> const OpenMS::CVTermList & { return self.getProductCVTermList(); }, nb::rv_policy::reference_internal)
        .def("setInterpretations", [](OpenMS::IncludeExcludeTarget& self, const std::vector<OpenMS::CVTermList>& interpretations) { return self.setInterpretations(interpretations); }, "interpretations"_a)
        .def("getInterpretations", [](const OpenMS::IncludeExcludeTarget& self) -> const std::vector<OpenMS::CVTermList> & { return self.getInterpretations(); }, nb::rv_policy::reference_internal)
        .def("addInterpretation", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::CVTermList& interpretation) { return self.addInterpretation(interpretation); }, "interpretation"_a)
        .def("setConfigurations", [](OpenMS::IncludeExcludeTarget& self, const std::vector<OpenMS::TargetedExperimentHelper::Configuration>& configuration) { return self.setConfigurations(configuration); }, "configuration"_a)
        .def("getConfigurations", [](const OpenMS::IncludeExcludeTarget& self) -> const std::vector<OpenMS::TargetedExperimentHelper::Configuration> & { return self.getConfigurations(); }, nb::rv_policy::reference_internal)
        .def("addConfiguration", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::TargetedExperimentHelper::Configuration& configuration) { return self.addConfiguration(configuration); }, "configuration"_a)
        .def("setPrediction", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::CVTermList& prediction) { return self.setPrediction(prediction); }, "prediction"_a)
        .def("addPredictionTerm", [](OpenMS::IncludeExcludeTarget& self, const OpenMS::CVTerm& prediction) { return self.addPredictionTerm(prediction); }, "prediction"_a)
        .def("getPrediction", [](const OpenMS::IncludeExcludeTarget& self) -> const OpenMS::CVTermList & { return self.getPrediction(); }, nb::rv_policy::reference_internal)
        .def("setRetentionTime", [](OpenMS::IncludeExcludeTarget& self, OpenMS::TargetedExperimentHelper::RetentionTime rt) { return self.setRetentionTime(rt); }, "rt"_a)
        .def("getRetentionTime", [](const OpenMS::IncludeExcludeTarget& self) -> const OpenMS::TargetedExperimentHelper::RetentionTime & { return self.getRetentionTime(); }, nb::rv_policy::reference_internal)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("__hash__", [](const OpenMS::IncludeExcludeTarget& self) { return std::hash<OpenMS::IncludeExcludeTarget>{}(self); })
        ;

    // -----------------------------------------------------------------------
    // Instrument
    // -----------------------------------------------------------------------
    // --- TargetedExperiment_Instrument (TargetedExperimentHelper::Instrument) ---
    nb::class_<OpenMS::TargetedExperimentHelper::Instrument, OpenMS::CVTermList>(m, "TargetedExperiment_Instrument",
        R"doc(
CVTermList

Instrument description used in targeted experiments (TraML).
This is a lightweight instrument reference with just an id field.
For the full MS instrument description, use the Instrument class instead.
)doc")
        .def(nb::init<>())
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)

        .def_rw("id", &OpenMS::TargetedExperimentHelper::Instrument::id)
        ;

    // --- Instrument (OpenMS::Instrument from METADATA) ---
    auto instrument_class = nb::class_<OpenMS::Instrument, OpenMS::MetaInfoInterface>(m, "Instrument",
        R"doc(
MetaInfoInterface

Description of a MS instrument.
Contains information about ion sources, mass analyzers, ion detectors,
software, vendor, model, and ion optics configuration.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Instrument &>())
        .def("__copy__", [](const OpenMS::Instrument& self) { return OpenMS::Instrument(self); })
        .def("__deepcopy__", [](const OpenMS::Instrument& self, nb::dict) { return OpenMS::Instrument(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getName", [](const OpenMS::Instrument& self) { return self.getName(); }, "Returns the name of the instrument")
        .def("setName", [](OpenMS::Instrument& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the instrument")
        .def("getVendor", [](const OpenMS::Instrument& self) { return self.getVendor(); }, "Returns the instrument vendor")
        .def("setVendor", [](OpenMS::Instrument& self, const OpenMS::String& vendor) { return self.setVendor(vendor); }, "vendor"_a, "Sets the instrument vendor")
        .def("getModel", [](const OpenMS::Instrument& self) { return self.getModel(); }, "Returns the instrument model")
        .def("setModel", [](OpenMS::Instrument& self, const OpenMS::String& model) { return self.setModel(model); }, "model"_a, "Sets the instrument model")
        .def("getCustomizations", [](const OpenMS::Instrument& self) { return self.getCustomizations(); }, "Returns a description of customizations")
        .def("setCustomizations", [](OpenMS::Instrument& self, const OpenMS::String& customizations) { return self.setCustomizations(customizations); }, "customizations"_a, "Sets a description of customizations")
        .def("getIonSources", [](const OpenMS::Instrument& self) -> const std::vector<OpenMS::IonSource>& { return self.getIonSources(); }, nb::rv_policy::reference_internal, "Returns a reference to the list of ion sources")
        .def("setIonSources", [](OpenMS::Instrument& self, const std::vector<OpenMS::IonSource>& ion_sources) { return self.setIonSources(ion_sources); }, "ion_sources"_a, "Sets the list of ion sources")
        .def("getIonDetectors", [](const OpenMS::Instrument& self) -> const std::vector<OpenMS::IonDetector>& { return self.getIonDetectors(); }, nb::rv_policy::reference_internal, "Returns a reference to the list of ion detectors")
        .def("setIonDetectors", [](OpenMS::Instrument& self, const std::vector<OpenMS::IonDetector>& ion_detectors) { return self.setIonDetectors(ion_detectors); }, "ion_detectors"_a, "Sets the list of ion detectors")
        .def("getMassAnalyzers", [](const OpenMS::Instrument& self) -> const std::vector<OpenMS::MassAnalyzer>& { return self.getMassAnalyzers(); }, nb::rv_policy::reference_internal, "Returns a reference to the list of mass analyzers")
        .def("setMassAnalyzers", [](OpenMS::Instrument& self, const std::vector<OpenMS::MassAnalyzer>& mass_analyzers) { return self.setMassAnalyzers(mass_analyzers); }, "mass_analyzers"_a, "Sets the list of mass analyzers")
        .def("getSoftware", [](const OpenMS::Instrument& self) -> const OpenMS::Software& { return self.getSoftware(); }, nb::rv_policy::reference_internal, "Returns a reference to the instrument software")
        .def("setSoftware", [](OpenMS::Instrument& self, const OpenMS::Software& software) { return self.setSoftware(software); }, "software"_a, "Sets the instrument software")
        .def("getIonOptics", [](const OpenMS::Instrument& self) { return self.getIonOptics(); }, "Returns the ion optics type")
        .def("setIonOptics", [](OpenMS::Instrument& self, OpenMS::Instrument::IonOpticsType ion_optics) { return self.setIonOptics(ion_optics); }, "ion_optics"_a, "Sets the ion optics type")
        .def_static("getAllNamesOfIonOpticsType", []() { return OpenMS::Instrument::getAllNamesOfIonOpticsType(); }, "Returns all ion optics type names known to OpenMS")
        ;
    // IonOpticsType enum nested under Instrument
    nb::enum_<OpenMS::Instrument::IonOpticsType>(instrument_class, "IonOpticsType")
        .value("UNKNOWN", OpenMS::Instrument::IonOpticsType::UNKNOWN)
        .value("MAGNETIC_DEFLECTION", OpenMS::Instrument::IonOpticsType::MAGNETIC_DEFLECTION)
        .value("DELAYED_EXTRACTION", OpenMS::Instrument::IonOpticsType::DELAYED_EXTRACTION)
        .value("COLLISION_QUADRUPOLE", OpenMS::Instrument::IonOpticsType::COLLISION_QUADRUPOLE)
        .value("SELECTED_ION_FLOW_TUBE", OpenMS::Instrument::IonOpticsType::SELECTED_ION_FLOW_TUBE)
        .value("TIME_LAG_FOCUSING", OpenMS::Instrument::IonOpticsType::TIME_LAG_FOCUSING)
        .value("REFLECTRON", OpenMS::Instrument::IonOpticsType::REFLECTRON)
        .value("EINZEL_LENS", OpenMS::Instrument::IonOpticsType::EINZEL_LENS)
        .value("FIRST_STABILITY_REGION", OpenMS::Instrument::IonOpticsType::FIRST_STABILITY_REGION)
        .value("FRINGING_FIELD", OpenMS::Instrument::IonOpticsType::FRINGING_FIELD)
        .value("KINETIC_ENERGY_ANALYZER", OpenMS::Instrument::IonOpticsType::KINETIC_ENERGY_ANALYZER)
        .value("STATIC_FIELD", OpenMS::Instrument::IonOpticsType::STATIC_FIELD)
        .value("SIZE_OF_IONOPTICSTYPE", OpenMS::Instrument::IonOpticsType::SIZE_OF_IONOPTICSTYPE)
        ;

    // -----------------------------------------------------------------------
    // InstrumentSettings
    // -----------------------------------------------------------------------
    auto instrumentsettings_class = nb::class_<OpenMS::InstrumentSettings, OpenMS::MetaInfoInterface>(m, "InstrumentSettings", 
        R"doc(
Description of the settings a MS Instrument was run with
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::InstrumentSettings &>())
        .def("__copy__", [](const OpenMS::InstrumentSettings& self) { return OpenMS::InstrumentSettings(self); })
        .def("__deepcopy__", [](const OpenMS::InstrumentSettings& self, nb::dict) { return OpenMS::InstrumentSettings(self); }, "memo"_a)
        .def_static("getAllNamesOfScanMode", []() { return OpenMS::InstrumentSettings::getAllNamesOfScanMode(); }, "Returns all scan mode names known to OpenMS")
        .def_static("scanModeToString", [](OpenMS::InstrumentSettings::ScanMode mode) { return OpenMS::InstrumentSettings::scanModeToString(mode); }, "mode"_a, "Convert a ScanMode enum to String. Throws Exception::InvalidValue if value is SIZE_OF_SCANMODE")
        .def_static("toScanMode", [](const OpenMS::String& name) { return OpenMS::InstrumentSettings::toScanMode(name); }, "name"_a, "Convert a string to ScanMode enum. Throws Exception::InvalidValue if name is not a valid scan mode")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getScanMode", [](const OpenMS::InstrumentSettings& self) { return self.getScanMode(); }, "Returns the scan mode")
        .def("setScanMode", [](OpenMS::InstrumentSettings& self, OpenMS::InstrumentSettings::ScanMode scan_mode) { return self.setScanMode(scan_mode); }, "scan_mode"_a, "Sets the scan mode")
        .def("getZoomScan", [](const OpenMS::InstrumentSettings& self) { return self.getZoomScan(); }, "Returns if this scan is a zoom (enhanced resolution) scan")
        .def("setZoomScan", [](OpenMS::InstrumentSettings& self, bool zoom_scan) { return self.setZoomScan(zoom_scan); }, "zoom_scan"_a, "Sets if this scan is a zoom (enhanced resolution) scan")
        .def("getPolarity", [](const OpenMS::InstrumentSettings& self) { return self.getPolarity(); }, "Returns the polarity")
        .def("setPolarity", [](OpenMS::InstrumentSettings& self, OpenMS::IonSource::Polarity polarity) { return self.setPolarity(polarity); }, "polarity"_a, "Sets the polarity")
        .def("getScanWindows", [](const OpenMS::InstrumentSettings& self) -> const std::vector<OpenMS::ScanWindow>& { return self.getScanWindows(); }, nb::rv_policy::reference_internal, "Returns a reference to the list of scan windows")
        .def("setScanWindows", [](OpenMS::InstrumentSettings& self, std::vector<OpenMS::ScanWindow> scan_windows) { return self.setScanWindows(scan_windows); }, "scan_windows"_a, "Sets the scan windows")

        ;
    // ScanMode enum nested under InstrumentSettings
    nb::enum_<OpenMS::InstrumentSettings::ScanMode>(instrumentsettings_class, "ScanMode")
        .value("UNKNOWN", OpenMS::InstrumentSettings::ScanMode::UNKNOWN)
        .value("MASSSPECTRUM", OpenMS::InstrumentSettings::ScanMode::MASSSPECTRUM)
        .value("MS1SPECTRUM", OpenMS::InstrumentSettings::ScanMode::MS1SPECTRUM)
        .value("MSNSPECTRUM", OpenMS::InstrumentSettings::ScanMode::MSNSPECTRUM)
        .value("SIM", OpenMS::InstrumentSettings::ScanMode::SIM)
        .value("SRM", OpenMS::InstrumentSettings::ScanMode::SRM)
        .value("CRM", OpenMS::InstrumentSettings::ScanMode::CRM)
        .value("CNG", OpenMS::InstrumentSettings::ScanMode::CNG)
        .value("CNL", OpenMS::InstrumentSettings::ScanMode::CNL)
        .value("PRECURSOR", OpenMS::InstrumentSettings::ScanMode::PRECURSOR)
        .value("EMC", OpenMS::InstrumentSettings::ScanMode::EMC)
        .value("TDF", OpenMS::InstrumentSettings::ScanMode::TDF)
        .value("EMR", OpenMS::InstrumentSettings::ScanMode::EMR)
        .value("EMISSION", OpenMS::InstrumentSettings::ScanMode::EMISSION)
        .value("ABSORPTION", OpenMS::InstrumentSettings::ScanMode::ABSORPTION)
        .value("SIZE_OF_SCANMODE", OpenMS::InstrumentSettings::ScanMode::SIZE_OF_SCANMODE)
        ;

    // -----------------------------------------------------------------------
    // IonDetector
    // -----------------------------------------------------------------------
    auto iondetector_class = nb::class_<OpenMS::IonDetector, OpenMS::MetaInfoInterface>(m, "IonDetector", 
        R"doc(
Description of a ion detector (part of a MS Instrument)
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IonDetector &>())
        .def("__copy__", [](const OpenMS::IonDetector& self) { return OpenMS::IonDetector(self); })
        .def("__deepcopy__", [](const OpenMS::IonDetector& self, nb::dict) { return OpenMS::IonDetector(self); }, "memo"_a)
        .def_static("getAllNamesOfType", []() { return OpenMS::IonDetector::getAllNamesOfType(); }, "Returns all detector type names known to OpenMS")
        .def_static("getAllNamesOfAcquisitionMode", []() { return OpenMS::IonDetector::getAllNamesOfAcquisitionMode(); }, "Returns all acquisition mode names known to OpenMS")
        .def_static("typeToString", [](OpenMS::IonDetector::Type type) { return OpenMS::IonDetector::typeToString(type); }, "type"_a, "Convert a Type enum to its string representation. Throws Exception::InvalidValue if type is SIZE_OF_TYPE")
        .def_static("toType", [](const OpenMS::String& name) { return OpenMS::IonDetector::toType(name); }, "name"_a, "Convert a string to a Type enum. Throws Exception::InvalidValue if name is not found")
        .def_static("acquisitionModeToString", [](OpenMS::IonDetector::AcquisitionMode mode) { return OpenMS::IonDetector::acquisitionModeToString(mode); }, "mode"_a, "Convert an AcquisitionMode enum to its string representation. Throws Exception::InvalidValue if mode is SIZE_OF_ACQUISITIONMODE")
        .def_static("toAcquisitionMode", [](const OpenMS::String& name) { return OpenMS::IonDetector::toAcquisitionMode(name); }, "name"_a, "Convert a string to an AcquisitionMode enum. Throws Exception::InvalidValue if name is not found")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getType", [](const OpenMS::IonDetector& self) { return self.getType(); }, "Returns the detector type")
        .def("setType", [](OpenMS::IonDetector& self, OpenMS::IonDetector::Type type) { return self.setType(type); }, "type"_a, "Sets the detector type")
        .def("getAcquisitionMode", [](const OpenMS::IonDetector& self) { return self.getAcquisitionMode(); }, "Returns the acquisition mode")
        .def("setAcquisitionMode", [](OpenMS::IonDetector& self, OpenMS::IonDetector::AcquisitionMode acquisition_mode) { return self.setAcquisitionMode(acquisition_mode); }, "acquisition_mode"_a, "Sets the acquisition mode")
        .def("getResolution", [](const OpenMS::IonDetector& self) { return self.getResolution(); }, "Returns the resolution (in ns)")
        .def("setResolution", [](OpenMS::IonDetector& self, double resolution) { return self.setResolution(resolution); }, "resolution"_a, "Sets the resolution (in ns)")
        .def("getADCSamplingFrequency", [](const OpenMS::IonDetector& self) { return self.getADCSamplingFrequency(); }, "Returns the analog-to-digital converter sampling frequency (in Hz)")
        .def("setADCSamplingFrequency", [](OpenMS::IonDetector& self, double ADC_sampling_frequency) { return self.setADCSamplingFrequency(ADC_sampling_frequency); }, "ADC_sampling_frequency"_a, "Sets the analog-to-digital converter sampling frequency (in Hz)")
        .def("getOrder", [](const OpenMS::IonDetector& self) { return self.getOrder(); }, "Returns the order")
        .def("setOrder", [](OpenMS::IonDetector& self, int order) { return self.setOrder(order); }, "order"_a, "Sets the order")
        
        .def("__hash__", [](const OpenMS::IonDetector& self) { return std::hash<OpenMS::IonDetector>{}(self); })
        ;
    // Type enum nested under IonDetector
    nb::enum_<OpenMS::IonDetector::Type>(iondetector_class, "Type")
        .value("TYPENULL", OpenMS::IonDetector::Type::TYPENULL)
        .value("ELECTRONMULTIPLIER", OpenMS::IonDetector::Type::ELECTRONMULTIPLIER)
        .value("PHOTOMULTIPLIER", OpenMS::IonDetector::Type::PHOTOMULTIPLIER)
        .value("FOCALPLANEARRAY", OpenMS::IonDetector::Type::FOCALPLANEARRAY)
        .value("FARADAYCUP", OpenMS::IonDetector::Type::FARADAYCUP)
        .value("CONVERSIONDYNODEELECTRONMULTIPLIER", OpenMS::IonDetector::Type::CONVERSIONDYNODEELECTRONMULTIPLIER)
        .value("CONVERSIONDYNODEPHOTOMULTIPLIER", OpenMS::IonDetector::Type::CONVERSIONDYNODEPHOTOMULTIPLIER)
        .value("MULTICOLLECTOR", OpenMS::IonDetector::Type::MULTICOLLECTOR)
        .value("CHANNELELECTRONMULTIPLIER", OpenMS::IonDetector::Type::CHANNELELECTRONMULTIPLIER)
        .value("CHANNELTRON", OpenMS::IonDetector::Type::CHANNELTRON)
        .value("DALYDETECTOR", OpenMS::IonDetector::Type::DALYDETECTOR)
        .value("MICROCHANNELPLATEDETECTOR", OpenMS::IonDetector::Type::MICROCHANNELPLATEDETECTOR)
        .value("ARRAYDETECTOR", OpenMS::IonDetector::Type::ARRAYDETECTOR)
        .value("CONVERSIONDYNODE", OpenMS::IonDetector::Type::CONVERSIONDYNODE)
        .value("DYNODE", OpenMS::IonDetector::Type::DYNODE)
        .value("FOCALPLANECOLLECTOR", OpenMS::IonDetector::Type::FOCALPLANECOLLECTOR)
        .value("IONTOPHOTONDETECTOR", OpenMS::IonDetector::Type::IONTOPHOTONDETECTOR)
        .value("POINTCOLLECTOR", OpenMS::IonDetector::Type::POINTCOLLECTOR)
        .value("POSTACCELERATIONDETECTOR", OpenMS::IonDetector::Type::POSTACCELERATIONDETECTOR)
        .value("PHOTODIODEARRAYDETECTOR", OpenMS::IonDetector::Type::PHOTODIODEARRAYDETECTOR)
        .value("INDUCTIVEDETECTOR", OpenMS::IonDetector::Type::INDUCTIVEDETECTOR)
        .value("ELECTRONMULTIPLIERTUBE", OpenMS::IonDetector::Type::ELECTRONMULTIPLIERTUBE)
        .value("SIZE_OF_TYPE", OpenMS::IonDetector::Type::SIZE_OF_TYPE)
        ;
    // AcquisitionMode enum nested under IonDetector
    nb::enum_<OpenMS::IonDetector::AcquisitionMode>(iondetector_class, "AcquisitionMode")
        .value("ACQMODENULL", OpenMS::IonDetector::AcquisitionMode::ACQMODENULL)
        .value("PULSECOUNTING", OpenMS::IonDetector::AcquisitionMode::PULSECOUNTING)
        .value("ADC", OpenMS::IonDetector::AcquisitionMode::ADC)
        .value("TDC", OpenMS::IonDetector::AcquisitionMode::TDC)
        .value("TRANSIENTRECORDER", OpenMS::IonDetector::AcquisitionMode::TRANSIENTRECORDER)
        .value("SIZE_OF_ACQUISITIONMODE", OpenMS::IonDetector::AcquisitionMode::SIZE_OF_ACQUISITIONMODE)
        ;

    // -----------------------------------------------------------------------
    // IonSource
    // -----------------------------------------------------------------------
    auto ionsource_class = nb::class_<OpenMS::IonSource, OpenMS::MetaInfoInterface>(m, "IonSource", 
        R"doc(
Description of an ion source (part of a MS Instrument)
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IonSource &>())
        .def("__copy__", [](const OpenMS::IonSource& self) { return OpenMS::IonSource(self); })
        .def("__deepcopy__", [](const OpenMS::IonSource& self, nb::dict) { return OpenMS::IonSource(self); }, "memo"_a)
        .def_static("getAllNamesOfInletType", []() { return OpenMS::IonSource::getAllNamesOfInletType(); }, "Returns all inlet type names known to OpenMS")
        .def_static("getAllNamesOfIonizationMethod", []() { return OpenMS::IonSource::getAllNamesOfIonizationMethod(); }, "Returns all ionization method names known to OpenMS")
        .def_static("getAllNamesOfPolarity", []() { return OpenMS::IonSource::getAllNamesOfPolarity(); }, "Returns all polarity names known to OpenMS")
        .def_static("inletTypeToString", [](OpenMS::IonSource::InletType type) { return OpenMS::IonSource::inletTypeToString(type); }, "type"_a, "Convert an InletType enum to its string representation")
        .def_static("toInletType", [](const OpenMS::String& name) { return OpenMS::IonSource::toInletType(name); }, "name"_a, "Convert a string to an InletType enum")
        .def_static("ionizationMethodToString", [](OpenMS::IonSource::IonizationMethod method) { return OpenMS::IonSource::ionizationMethodToString(method); }, "method"_a, "Convert an IonizationMethod enum to its string representation")
        .def_static("toIonizationMethod", [](const OpenMS::String& name) { return OpenMS::IonSource::toIonizationMethod(name); }, "name"_a, "Convert a string to an IonizationMethod enum")
        .def_static("polarityToString", [](OpenMS::IonSource::Polarity polarity) { return OpenMS::IonSource::polarityToString(polarity); }, "polarity"_a, "Convert a Polarity enum to its string representation")
        .def_static("toPolarity", [](const OpenMS::String& name) { return OpenMS::IonSource::toPolarity(name); }, "name"_a, "Convert a string to a Polarity enum")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getInletType", [](const OpenMS::IonSource& self) { return self.getInletType(); }, "Returns the inlet type")
        .def("setInletType", [](OpenMS::IonSource& self, OpenMS::IonSource::InletType inlet_type) { return self.setInletType(inlet_type); }, "inlet_type"_a, "Sets the  inlet type")
        .def("getIonizationMethod", [](const OpenMS::IonSource& self) { return self.getIonizationMethod(); }, "Returns the ionization method")
        .def("setIonizationMethod", [](OpenMS::IonSource& self, OpenMS::IonSource::IonizationMethod ionization_type) { return self.setIonizationMethod(ionization_type); }, "ionization_type"_a, "Sets the ionization method")
        .def("getPolarity", [](const OpenMS::IonSource& self) { return self.getPolarity(); }, "Returns the ionization mode")
        .def("setPolarity", [](OpenMS::IonSource& self, OpenMS::IonSource::Polarity polarity) { return self.setPolarity(polarity); }, "polarity"_a, "Sets the ionization mode")
        .def("getOrder", [](const OpenMS::IonSource& self) { return self.getOrder(); })
        .def("setOrder", [](OpenMS::IonSource& self, int order) { return self.setOrder(order); }, "order"_a, "Sets the order")
        
        .def("__hash__", [](const OpenMS::IonSource& self) { return std::hash<OpenMS::IonSource>{}(self); })
        ;
    // Polarity enum nested under IonSource
    nb::enum_<OpenMS::IonSource::Polarity>(ionsource_class, "Polarity")
        .value("POLNULL", OpenMS::IonSource::Polarity::POLNULL)
        .value("POSITIVE", OpenMS::IonSource::Polarity::POSITIVE)
        .value("NEGATIVE", OpenMS::IonSource::Polarity::NEGATIVE)
        .value("SIZE_OF_POLARITY", OpenMS::IonSource::Polarity::SIZE_OF_POLARITY)
        ;
    // InletType enum nested under IonSource
    nb::enum_<OpenMS::IonSource::InletType>(ionsource_class, "InletType")
        .value("INLETNULL", OpenMS::IonSource::InletType::INLETNULL)
        .value("DIRECT", OpenMS::IonSource::InletType::DIRECT)
        .value("BATCH", OpenMS::IonSource::InletType::BATCH)
        .value("CHROMATOGRAPHY", OpenMS::IonSource::InletType::CHROMATOGRAPHY)
        .value("PARTICLEBEAM", OpenMS::IonSource::InletType::PARTICLEBEAM)
        .value("MEMBRANESEPARATOR", OpenMS::IonSource::InletType::MEMBRANESEPARATOR)
        .value("OPENSPLIT", OpenMS::IonSource::InletType::OPENSPLIT)
        .value("JETSEPARATOR", OpenMS::IonSource::InletType::JETSEPARATOR)
        .value("SEPTUM", OpenMS::IonSource::InletType::SEPTUM)
        .value("RESERVOIR", OpenMS::IonSource::InletType::RESERVOIR)
        .value("MOVINGBELT", OpenMS::IonSource::InletType::MOVINGBELT)
        .value("MOVINGWIRE", OpenMS::IonSource::InletType::MOVINGWIRE)
        .value("FLOWINJECTIONANALYSIS", OpenMS::IonSource::InletType::FLOWINJECTIONANALYSIS)
        .value("ELECTROSPRAYINLET", OpenMS::IonSource::InletType::ELECTROSPRAYINLET)
        .value("THERMOSPRAYINLET", OpenMS::IonSource::InletType::THERMOSPRAYINLET)
        .value("INFUSION", OpenMS::IonSource::InletType::INFUSION)
        .value("CONTINUOUSFLOWFASTATOMBOMBARDMENT", OpenMS::IonSource::InletType::CONTINUOUSFLOWFASTATOMBOMBARDMENT)
        .value("INDUCTIVELYCOUPLEDPLASMA", OpenMS::IonSource::InletType::INDUCTIVELYCOUPLEDPLASMA)
        .value("MEMBRANE", OpenMS::IonSource::InletType::MEMBRANE)
        .value("NANOSPRAY", OpenMS::IonSource::InletType::NANOSPRAY)
        .value("SIZE_OF_INLETTYPE", OpenMS::IonSource::InletType::SIZE_OF_INLETTYPE)
        ;
    // IonizationMethod enum nested under IonSource
    nb::enum_<OpenMS::IonSource::IonizationMethod>(ionsource_class, "IonizationMethod")
        .value("IONMETHODNULL", OpenMS::IonSource::IonizationMethod::IONMETHODNULL)
        .value("ESI", OpenMS::IonSource::IonizationMethod::ESI)
        .value("EI", OpenMS::IonSource::IonizationMethod::EI)
        .value("CI", OpenMS::IonSource::IonizationMethod::CI)
        .value("FAB", OpenMS::IonSource::IonizationMethod::FAB)
        .value("TSP", OpenMS::IonSource::IonizationMethod::TSP)
        .value("LD", OpenMS::IonSource::IonizationMethod::LD)
        .value("FD", OpenMS::IonSource::IonizationMethod::FD)
        .value("FI", OpenMS::IonSource::IonizationMethod::FI)
        .value("PD", OpenMS::IonSource::IonizationMethod::PD)
        .value("SI", OpenMS::IonSource::IonizationMethod::SI)
        .value("TI", OpenMS::IonSource::IonizationMethod::TI)
        .value("API", OpenMS::IonSource::IonizationMethod::API)
        .value("ISI", OpenMS::IonSource::IonizationMethod::ISI)
        .value("CID", OpenMS::IonSource::IonizationMethod::CID)
        .value("CAD", OpenMS::IonSource::IonizationMethod::CAD)
        .value("HN", OpenMS::IonSource::IonizationMethod::HN)
        .value("APCI", OpenMS::IonSource::IonizationMethod::APCI)
        .value("APPI", OpenMS::IonSource::IonizationMethod::APPI)
        .value("ICP", OpenMS::IonSource::IonizationMethod::ICP)
        .value("NESI", OpenMS::IonSource::IonizationMethod::NESI)
        .value("MESI", OpenMS::IonSource::IonizationMethod::MESI)
        .value("SELDI", OpenMS::IonSource::IonizationMethod::SELDI)
        .value("SEND", OpenMS::IonSource::IonizationMethod::SEND)
        .value("FIB", OpenMS::IonSource::IonizationMethod::FIB)
        .value("MALDI", OpenMS::IonSource::IonizationMethod::MALDI)
        .value("MPI", OpenMS::IonSource::IonizationMethod::MPI)
        .value("DI", OpenMS::IonSource::IonizationMethod::DI)
        .value("FA", OpenMS::IonSource::IonizationMethod::FA)
        .value("FII", OpenMS::IonSource::IonizationMethod::FII)
        .value("GD_MS", OpenMS::IonSource::IonizationMethod::GD_MS)
        .value("NICI", OpenMS::IonSource::IonizationMethod::NICI)
        .value("NRMS", OpenMS::IonSource::IonizationMethod::NRMS)
        .value("PI", OpenMS::IonSource::IonizationMethod::PI)
        .value("PYMS", OpenMS::IonSource::IonizationMethod::PYMS)
        .value("REMPI", OpenMS::IonSource::IonizationMethod::REMPI)
        .value("AI", OpenMS::IonSource::IonizationMethod::AI)
        .value("ASI", OpenMS::IonSource::IonizationMethod::ASI)
        .value("AD", OpenMS::IonSource::IonizationMethod::AD)
        .value("AUI", OpenMS::IonSource::IonizationMethod::AUI)
        .value("CEI", OpenMS::IonSource::IonizationMethod::CEI)
        .value("CHEMI", OpenMS::IonSource::IonizationMethod::CHEMI)
        .value("DISSI", OpenMS::IonSource::IonizationMethod::DISSI)
        .value("LSI", OpenMS::IonSource::IonizationMethod::LSI)
        .value("PEI", OpenMS::IonSource::IonizationMethod::PEI)
        .value("SOI", OpenMS::IonSource::IonizationMethod::SOI)
        .value("SPI", OpenMS::IonSource::IonizationMethod::SPI)
        .value("SUI", OpenMS::IonSource::IonizationMethod::SUI)
        .value("VI", OpenMS::IonSource::IonizationMethod::VI)
        .value("AP_MALDI", OpenMS::IonSource::IonizationMethod::AP_MALDI)
        .value("SILI", OpenMS::IonSource::IonizationMethod::SILI)
        .value("SALDI", OpenMS::IonSource::IonizationMethod::SALDI)
        .value("SIZE_OF_IONIZATIONMETHOD", OpenMS::IonSource::IonizationMethod::SIZE_OF_IONIZATIONMETHOD)
        ;

    // -----------------------------------------------------------------------
    // MassAnalyzer
    // -----------------------------------------------------------------------
    auto massanalyzer_class = nb::class_<OpenMS::MassAnalyzer, OpenMS::MetaInfoInterface>(m, "MassAnalyzer", 
        R"doc(
Description of a mass analyzer (part of a MS Instrument)
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MassAnalyzer &>())
        .def("__copy__", [](const OpenMS::MassAnalyzer& self) { return OpenMS::MassAnalyzer(self); })
        .def("__deepcopy__", [](const OpenMS::MassAnalyzer& self, nb::dict) { return OpenMS::MassAnalyzer(self); }, "memo"_a)
        .def_static("getAllNamesOfAnalyzerType", []() { return OpenMS::MassAnalyzer::getAllNamesOfAnalyzerType(); }, "Returns all analyzer type names known to OpenMS")
        .def_static("getAllNamesOfResolutionMethod", []() { return OpenMS::MassAnalyzer::getAllNamesOfResolutionMethod(); }, "Returns all resolution method names known to OpenMS")
        .def_static("getAllNamesOfResolutionType", []() { return OpenMS::MassAnalyzer::getAllNamesOfResolutionType(); }, "Returns all resolution type names known to OpenMS")
        .def_static("getAllNamesOfScanDirection", []() { return OpenMS::MassAnalyzer::getAllNamesOfScanDirection(); }, "Returns all scan direction names known to OpenMS")
        .def_static("getAllNamesOfScanLaw", []() { return OpenMS::MassAnalyzer::getAllNamesOfScanLaw(); }, "Returns all scan law names known to OpenMS")
        .def_static("getAllNamesOfReflectronState", []() { return OpenMS::MassAnalyzer::getAllNamesOfReflectronState(); }, "Returns all reflectron state names known to OpenMS")
        .def_static("analyzerTypeToString", [](OpenMS::MassAnalyzer::AnalyzerType type) { return OpenMS::MassAnalyzer::analyzerTypeToString(type); }, "type"_a, "Convert AnalyzerType enum to string")
        .def_static("toAnalyzerType", [](const OpenMS::String& name) { return OpenMS::MassAnalyzer::toAnalyzerType(name); }, "name"_a, "Convert string to AnalyzerType enum")
        .def_static("resolutionMethodToString", [](OpenMS::MassAnalyzer::ResolutionMethod method) { return OpenMS::MassAnalyzer::resolutionMethodToString(method); }, "method"_a, "Convert ResolutionMethod enum to string")
        .def_static("toResolutionMethod", [](const OpenMS::String& name) { return OpenMS::MassAnalyzer::toResolutionMethod(name); }, "name"_a, "Convert string to ResolutionMethod enum")
        .def_static("resolutionTypeToString", [](OpenMS::MassAnalyzer::ResolutionType type) { return OpenMS::MassAnalyzer::resolutionTypeToString(type); }, "type"_a, "Convert ResolutionType enum to string")
        .def_static("toResolutionType", [](const OpenMS::String& name) { return OpenMS::MassAnalyzer::toResolutionType(name); }, "name"_a, "Convert string to ResolutionType enum")
        .def_static("scanDirectionToString", [](OpenMS::MassAnalyzer::ScanDirection direction) { return OpenMS::MassAnalyzer::scanDirectionToString(direction); }, "direction"_a, "Convert ScanDirection enum to string")
        .def_static("toScanDirection", [](const OpenMS::String& name) { return OpenMS::MassAnalyzer::toScanDirection(name); }, "name"_a, "Convert string to ScanDirection enum")
        .def_static("scanLawToString", [](OpenMS::MassAnalyzer::ScanLaw law) { return OpenMS::MassAnalyzer::scanLawToString(law); }, "law"_a, "Convert ScanLaw enum to string")
        .def_static("toScanLaw", [](const OpenMS::String& name) { return OpenMS::MassAnalyzer::toScanLaw(name); }, "name"_a, "Convert string to ScanLaw enum")
        .def_static("reflectronStateToString", [](OpenMS::MassAnalyzer::ReflectronState state) { return OpenMS::MassAnalyzer::reflectronStateToString(state); }, "state"_a, "Convert ReflectronState enum to string")
        .def_static("toReflectronState", [](const OpenMS::String& name) { return OpenMS::MassAnalyzer::toReflectronState(name); }, "name"_a, "Convert string to ReflectronState enum")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getType", [](const OpenMS::MassAnalyzer& self) { return self.getType(); }, "Returns the analyzer type")
        .def("setType", [](OpenMS::MassAnalyzer& self, OpenMS::MassAnalyzer::AnalyzerType type) { return self.setType(type); }, "type"_a, "Sets the analyzer type")
        .def("getResolutionMethod", [](const OpenMS::MassAnalyzer& self) { return self.getResolutionMethod(); }, "Returns the method used for determination of the resolution")
        .def("setResolutionMethod", [](OpenMS::MassAnalyzer& self, OpenMS::MassAnalyzer::ResolutionMethod resolution_method) { return self.setResolutionMethod(resolution_method); }, "resolution_method"_a, "Sets the method used for determination of the resolution")
        .def("getResolutionType", [](const OpenMS::MassAnalyzer& self) { return self.getResolutionType(); }, "Returns the resolution type")
        .def("setResolutionType", [](OpenMS::MassAnalyzer& self, OpenMS::MassAnalyzer::ResolutionType resolution_type) { return self.setResolutionType(resolution_type); }, "resolution_type"_a, "Sets the resolution type")
        .def("getScanDirection", [](const OpenMS::MassAnalyzer& self) { return self.getScanDirection(); }, "Returns the direction of scanning")
        .def("setScanDirection", [](OpenMS::MassAnalyzer& self, OpenMS::MassAnalyzer::ScanDirection scan_direction) { return self.setScanDirection(scan_direction); }, "scan_direction"_a, "Sets the direction of scanning")
        .def("getScanLaw", [](const OpenMS::MassAnalyzer& self) { return self.getScanLaw(); }, "Returns the scan law")
        .def("setScanLaw", [](OpenMS::MassAnalyzer& self, OpenMS::MassAnalyzer::ScanLaw scan_law) { return self.setScanLaw(scan_law); }, "scan_law"_a, "Sets the scan law")
        .def("getReflectronState", [](const OpenMS::MassAnalyzer& self) { return self.getReflectronState(); }, "Returns the reflectron state (for TOF)")
        .def("setReflectronState", [](OpenMS::MassAnalyzer& self, OpenMS::MassAnalyzer::ReflectronState reflecton_state) { return self.setReflectronState(reflecton_state); }, "reflecton_state"_a, "Sets the reflectron state (for TOF)")
        .def("getResolution", [](const OpenMS::MassAnalyzer& self) { return self.getResolution(); }, "Returns the resolution. The maximum m/z value at which two peaks can be resolved, according to one of the standard measures")
        .def("setResolution", [](OpenMS::MassAnalyzer& self, double resolution) { return self.setResolution(resolution); }, "resolution"_a, "Sets the resolution")
        .def("getAccuracy", [](const OpenMS::MassAnalyzer& self) { return self.getAccuracy(); }, "Returns the mass accuracy i.e. how much the theoretical mass may differ from the measured mass (in ppm)")
        .def("setAccuracy", [](OpenMS::MassAnalyzer& self, double accuracy) { return self.setAccuracy(accuracy); }, "accuracy"_a, "Sets the accuracy i.e. how much the theoretical mass may differ from the measured mass (in ppm)")
        .def("getScanRate", [](const OpenMS::MassAnalyzer& self) { return self.getScanRate(); }, "Returns the scan rate (in s)")
        .def("setScanRate", [](OpenMS::MassAnalyzer& self, double scan_rate) { return self.setScanRate(scan_rate); }, "scan_rate"_a, "Sets the scan rate (in s)")
        .def("getScanTime", [](const OpenMS::MassAnalyzer& self) { return self.getScanTime(); }, "Returns the scan time for a single scan (in s)")
        .def("setScanTime", [](OpenMS::MassAnalyzer& self, double scan_time) { return self.setScanTime(scan_time); }, "scan_time"_a, "Sets the scan time for a single scan (in s)")
        .def("getTOFTotalPathLength", [](const OpenMS::MassAnalyzer& self) { return self.getTOFTotalPathLength(); }, "Returns the path length for a TOF mass analyzer (in meter)")
        .def("setTOFTotalPathLength", [](OpenMS::MassAnalyzer& self, double TOF_total_path_length) { return self.setTOFTotalPathLength(TOF_total_path_length); }, "TOF_total_path_length"_a, "Sets the path length for a TOF mass analyzer (in meter)")
        .def("getIsolationWidth", [](const OpenMS::MassAnalyzer& self) { return self.getIsolationWidth(); }, "Returns the isolation width i.e. in which m/z range the precursor ion is selected for MS to the n (in m/z)")
        .def("setIsolationWidth", [](OpenMS::MassAnalyzer& self, double isolation_width) { return self.setIsolationWidth(isolation_width); }, "isolation_width"_a, "Sets the isolation width i.e. in which m/z range the precursor ion is selected for MS to the n (in m/z)")
        .def("getFinalMSExponent", [](const OpenMS::MassAnalyzer& self) { return self.getFinalMSExponent(); }, "Returns the final MS exponent")
        .def("setFinalMSExponent", [](OpenMS::MassAnalyzer& self, int final_MS_exponent) { return self.setFinalMSExponent(final_MS_exponent); }, "final_MS_exponent"_a, "Sets the final MS exponent")
        .def("getMagneticFieldStrength", [](const OpenMS::MassAnalyzer& self) { return self.getMagneticFieldStrength(); }, "Returns the strength of the magnetic field (in T)")
        .def("setMagneticFieldStrength", [](OpenMS::MassAnalyzer& self, double magnetic_field_strength) { return self.setMagneticFieldStrength(magnetic_field_strength); }, "magnetic_field_strength"_a, "Sets the strength of the magnetic field (in T)")
        .def("getOrder", [](const OpenMS::MassAnalyzer& self) { return self.getOrder(); }, "Returns the position of this part in the whole Instrument")
        .def("setOrder", [](OpenMS::MassAnalyzer& self, int order) { return self.setOrder(order); }, "order"_a, "Sets the order")
        
        .def("__hash__", [](const OpenMS::MassAnalyzer& self) { return std::hash<OpenMS::MassAnalyzer>{}(self); })
        ;
    // AnalyzerType enum nested under MassAnalyzer
    nb::enum_<OpenMS::MassAnalyzer::AnalyzerType>(massanalyzer_class, "AnalyzerType")
        .value("ANALYZERNULL", OpenMS::MassAnalyzer::AnalyzerType::ANALYZERNULL)
        .value("QUADRUPOLE", OpenMS::MassAnalyzer::AnalyzerType::QUADRUPOLE)
        .value("PAULIONTRAP", OpenMS::MassAnalyzer::AnalyzerType::PAULIONTRAP)
        .value("RADIALEJECTIONLINEARIONTRAP", OpenMS::MassAnalyzer::AnalyzerType::RADIALEJECTIONLINEARIONTRAP)
        .value("AXIALEJECTIONLINEARIONTRAP", OpenMS::MassAnalyzer::AnalyzerType::AXIALEJECTIONLINEARIONTRAP)
        .value("TOF", OpenMS::MassAnalyzer::AnalyzerType::TOF)
        .value("SECTOR", OpenMS::MassAnalyzer::AnalyzerType::SECTOR)
        .value("FOURIERTRANSFORM", OpenMS::MassAnalyzer::AnalyzerType::FOURIERTRANSFORM)
        .value("IONSTORAGE", OpenMS::MassAnalyzer::AnalyzerType::IONSTORAGE)
        .value("ESA", OpenMS::MassAnalyzer::AnalyzerType::ESA)
        .value("IT", OpenMS::MassAnalyzer::AnalyzerType::IT)
        .value("SWIFT", OpenMS::MassAnalyzer::AnalyzerType::SWIFT)
        .value("CYCLOTRON", OpenMS::MassAnalyzer::AnalyzerType::CYCLOTRON)
        .value("ORBITRAP", OpenMS::MassAnalyzer::AnalyzerType::ORBITRAP)
        .value("LIT", OpenMS::MassAnalyzer::AnalyzerType::LIT)
        .value("SIZE_OF_ANALYZERTYPE", OpenMS::MassAnalyzer::AnalyzerType::SIZE_OF_ANALYZERTYPE)
        ;
    // ResolutionMethod enum nested under MassAnalyzer
    nb::enum_<OpenMS::MassAnalyzer::ResolutionMethod>(massanalyzer_class, "ResolutionMethod")
        .value("RESMETHNULL", OpenMS::MassAnalyzer::ResolutionMethod::RESMETHNULL)
        .value("FWHM", OpenMS::MassAnalyzer::ResolutionMethod::FWHM)
        .value("TENPERCENTVALLEY", OpenMS::MassAnalyzer::ResolutionMethod::TENPERCENTVALLEY)
        .value("BASELINE", OpenMS::MassAnalyzer::ResolutionMethod::BASELINE)
        .value("SIZE_OF_RESOLUTIONMETHOD", OpenMS::MassAnalyzer::ResolutionMethod::SIZE_OF_RESOLUTIONMETHOD)
        ;
    // ResolutionType enum nested under MassAnalyzer
    nb::enum_<OpenMS::MassAnalyzer::ResolutionType>(massanalyzer_class, "ResolutionType")
        .value("RESTYPENULL", OpenMS::MassAnalyzer::ResolutionType::RESTYPENULL)
        .value("CONSTANT", OpenMS::MassAnalyzer::ResolutionType::CONSTANT)
        .value("PROPORTIONAL", OpenMS::MassAnalyzer::ResolutionType::PROPORTIONAL)
        .value("SIZE_OF_RESOLUTIONTYPE", OpenMS::MassAnalyzer::ResolutionType::SIZE_OF_RESOLUTIONTYPE)
        ;
    // ScanDirection enum nested under MassAnalyzer
    nb::enum_<OpenMS::MassAnalyzer::ScanDirection>(massanalyzer_class, "ScanDirection")
        .value("SCANDIRNULL", OpenMS::MassAnalyzer::ScanDirection::SCANDIRNULL)
        .value("UP", OpenMS::MassAnalyzer::ScanDirection::UP)
        .value("DOWN", OpenMS::MassAnalyzer::ScanDirection::DOWN)
        .value("SIZE_OF_SCANDIRECTION", OpenMS::MassAnalyzer::ScanDirection::SIZE_OF_SCANDIRECTION)
        ;
    // ScanLaw enum nested under MassAnalyzer
    nb::enum_<OpenMS::MassAnalyzer::ScanLaw>(massanalyzer_class, "ScanLaw")
        .value("SCANLAWNULL", OpenMS::MassAnalyzer::ScanLaw::SCANLAWNULL)
        .value("EXPONENTIAL", OpenMS::MassAnalyzer::ScanLaw::EXPONENTIAL)
        .value("LINEAR", OpenMS::MassAnalyzer::ScanLaw::LINEAR)
        .value("QUADRATIC", OpenMS::MassAnalyzer::ScanLaw::QUADRATIC)
        .value("SIZE_OF_SCANLAW", OpenMS::MassAnalyzer::ScanLaw::SIZE_OF_SCANLAW)
        ;
    // ReflectronState enum nested under MassAnalyzer
    nb::enum_<OpenMS::MassAnalyzer::ReflectronState>(massanalyzer_class, "ReflectronState")
        .value("REFLSTATENULL", OpenMS::MassAnalyzer::ReflectronState::REFLSTATENULL)
        .value("ON", OpenMS::MassAnalyzer::ReflectronState::ON)
        .value("OFF", OpenMS::MassAnalyzer::ReflectronState::OFF)
        .value("NONE", OpenMS::MassAnalyzer::ReflectronState::NONE)
        .value("SIZE_OF_REFLECTRONSTATE", OpenMS::MassAnalyzer::ReflectronState::SIZE_OF_REFLECTRONSTATE)
        ;

    // -----------------------------------------------------------------------
    // MetaInfoDescription
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MetaInfoDescription, OpenMS::MetaInfoInterface>(m, "MetaInfoDescription", 
        R"doc(
Description of the meta data arrays of MSSpectrum
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MetaInfoDescription &>())
        .def("__copy__", [](const OpenMS::MetaInfoDescription& self) { return OpenMS::MetaInfoDescription(self); })
        .def("__deepcopy__", [](const OpenMS::MetaInfoDescription& self, nb::dict) { return OpenMS::MetaInfoDescription(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getName", [](const OpenMS::MetaInfoDescription& self) { return self.getName(); }, "Returns the name of the peak annotations")
        .def("setName", [](OpenMS::MetaInfoDescription& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the peak annotations")
        .def("getDataProcessing", [](OpenMS::MetaInfoDescription& self) -> std::vector<std::shared_ptr<OpenMS::DataProcessing>> & { return self.getDataProcessing(); }, nb::rv_policy::reference_internal, "Returns a reference to the description of the applied processing")
        .def("setDataProcessing", [](OpenMS::MetaInfoDescription& self, const std::vector<std::shared_ptr<OpenMS::DataProcessing>>& data_processing) { return self.setDataProcessing(data_processing); }, "data_processing"_a, "Sets the description of the applied processing")
        
        ;

    // -----------------------------------------------------------------------
    // FloatDataArray
    // -----------------------------------------------------------------------
    auto floatdataarray_class = nb::class_<OpenMS::DataArrays::FloatDataArray>(m, "FloatDataArray", 
        R"doc(
MetaInfoDescription

The representation of extra float data attached to a spectrum or chromatogram.
Raw data access is provided by `get_peaks` and `set_peaks`, which yields numpy arrays.
Commonly used for storing ion mobility values or other per-peak float annotations.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DataArrays::FloatDataArray &>())
        .def("__copy__", [](const OpenMS::DataArrays::FloatDataArray& self) { return OpenMS::DataArrays::FloatDataArray(self); })
        .def("__deepcopy__", [](const OpenMS::DataArrays::FloatDataArray& self, nb::dict) { return OpenMS::DataArrays::FloatDataArray(self); }, "memo"_a)
        .def(nb::self != nb::self)
        .def(nb::self == nb::self)
        .def("setDataProcessing", [](OpenMS::DataArrays::FloatDataArray& self, const std::vector<std::shared_ptr<OpenMS::DataProcessing>>& data_processing) { return self.setDataProcessing(data_processing); }, "data_processing"_a, "Sets the description of the applied processing")
        
        .def("__len__", [](OpenMS::DataArrays::FloatDataArray& self) { return self.size(); })
        .def("size", [](const OpenMS::DataArrays::FloatDataArray& self) { return self.size(); }, "Returns the number of elements")
        .def("__getitem__", [](OpenMS::DataArrays::FloatDataArray& self, size_t i) -> float& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)
        .def("__setitem__", [](OpenMS::DataArrays::FloatDataArray& self, size_t i, float val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        })

        .def(nb::init<>(), "Default constructor")
        .def(nb::init<const OpenMS::DataArrays::FloatDataArray&>(), "Copy constructor")
        .def("__copy__", [](const OpenMS::DataArrays::FloatDataArray& self) { return OpenMS::DataArrays::FloatDataArray(self); })
        .def("__deepcopy__", [](const OpenMS::DataArrays::FloatDataArray& self, nb::dict) { return OpenMS::DataArrays::FloatDataArray(self); }, "memo"_a)

        .def("setName", [](OpenMS::DataArrays::FloatDataArray& self, const OpenMS::String& name) {
            self.setName(name);
        }, "name"_a, "Set the name")

        .def("getName", [](const OpenMS::DataArrays::FloatDataArray& self) {
            return self.getName();
        }, "Get the name")

        .def("get_data", [](const OpenMS::DataArrays::FloatDataArray& self) {
            // Return as numpy array - copy data into new array
            size_t n = self.size();
            float* data = new float[n];
            std::copy(self.begin(), self.end(), data);
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(data, {n}, owner);
        }, "Returns a copy of the data as numpy array")

        .def("get_data_view", [](nb::object self_obj) {
            // Return a writable view (zero-copy), empty array if empty
            auto& self = nb::cast<OpenMS::DataArrays::FloatDataArray&>(self_obj);
            float* data_ptr = self.empty() ? nullptr : self.data();
            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                data_ptr, {self.size()}, self_obj
            );
        }, "Returns a zero-copy writable view of the data as numpy array (empty array if empty)")

        .def("set_data", [](OpenMS::DataArrays::FloatDataArray& self, nb::object data_obj) {
            // Fast path: float32 numpy array — direct memcpy
            nb::ndarray<nb::numpy, float, nb::ndim<1>> f32_arr;
            if (nb::try_cast(data_obj, f32_arr)) {
                const size_t n = f32_arr.shape(0);
                self.resize(n);
                std::copy(static_cast<const float*>(f32_arr.data()),
                          static_cast<const float*>(f32_arr.data()) + n, self.data());
                return;
            }
            // Fast path: float64 numpy array — single-pass narrowing copy
            nb::ndarray<nb::numpy, double, nb::ndim<1>> f64_arr;
            if (nb::try_cast(data_obj, f64_arr)) {
                const size_t n = f64_arr.shape(0);
                self.resize(n);
                const double* ptr = static_cast<const double*>(f64_arr.data());
                for (size_t i = 0; i < n; ++i) self[i] = static_cast<float>(ptr[i]);
                return;
            }
            // Fallback: any iterable (list, etc.)
            std::vector<float> vec = nb::cast<std::vector<float>>(data_obj);
            self.assign(vec.begin(), vec.end());
        }, "data"_a, "Set data from a numpy array or list")

        .def("push_back", [](OpenMS::DataArrays::FloatDataArray& self, float val) {
            self.push_back(val);
        }, "val"_a, "Add a value to the array")

        .def("resize", [](OpenMS::DataArrays::FloatDataArray& self, size_t n) {
            self.resize(n);
        }, "n"_a, "Resize the array")

        .def("clear", [](OpenMS::DataArrays::FloatDataArray& self) {
            self.clear();
        }, "Clear the array")
        .def("reserve", [](OpenMS::DataArrays::FloatDataArray& self, size_t n) { self.reserve(n); }, "n"_a, "Reserve memory for n elements")
        .def("getDataProcessing", [](const OpenMS::DataArrays::FloatDataArray& self) { return self.getDataProcessing(); }, "Returns the data processing objects")
        .def("__repr__", [](const OpenMS::DataArrays::FloatDataArray& self) {
            return "FloatDataArray(name='" + self.getName() + "', size=" + std::to_string(self.size()) + ")";
        })
        .def("__str__", [](const OpenMS::DataArrays::FloatDataArray& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;
    def_MetaInfoInterface<OpenMS::DataArrays::FloatDataArray>(floatdataarray_class);

    // -----------------------------------------------------------------------
    // IntegerDataArray
    // -----------------------------------------------------------------------
    auto integerdataarray_class = nb::class_<OpenMS::DataArrays::IntegerDataArray>(m, "IntegerDataArray", 
        R"doc(
MetaInfoDescription

The representation of extra integer data attached to a spectrum or chromatogram.
Raw data access is provided by `get_peaks` and `set_peaks`, which yields numpy arrays.
Used for storing per-peak integer annotations.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DataArrays::IntegerDataArray &>())
        .def("__copy__", [](const OpenMS::DataArrays::IntegerDataArray& self) { return OpenMS::DataArrays::IntegerDataArray(self); })
        .def("__deepcopy__", [](const OpenMS::DataArrays::IntegerDataArray& self, nb::dict) { return OpenMS::DataArrays::IntegerDataArray(self); }, "memo"_a)
        .def(nb::self != nb::self)
        .def(nb::self == nb::self)
        .def("getName", [](const OpenMS::DataArrays::IntegerDataArray& self) { return self.getName(); }, "Returns the name of the peak annotations")
        .def("setName", [](OpenMS::DataArrays::IntegerDataArray& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the peak annotations")
        .def("setDataProcessing", [](OpenMS::DataArrays::IntegerDataArray& self, const std::vector<std::shared_ptr<OpenMS::DataProcessing>>& data_processing) { return self.setDataProcessing(data_processing); }, "data_processing"_a, "Sets the description of the applied processing")
        
        .def("__len__", [](OpenMS::DataArrays::IntegerDataArray& self) { return self.size(); })
        .def("size", [](const OpenMS::DataArrays::IntegerDataArray& self) { return self.size(); }, "Returns the number of elements")
        .def("__getitem__", [](OpenMS::DataArrays::IntegerDataArray& self, size_t i) -> int& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)
        .def("__setitem__", [](OpenMS::DataArrays::IntegerDataArray& self, size_t i, int val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        })

        .def(nb::init<>(), "Default constructor")
        .def(nb::init<const OpenMS::DataArrays::IntegerDataArray&>(), "Copy constructor")
        .def("__copy__", [](const OpenMS::DataArrays::IntegerDataArray& self) { return OpenMS::DataArrays::IntegerDataArray(self); })
        .def("__deepcopy__", [](const OpenMS::DataArrays::IntegerDataArray& self, nb::dict) { return OpenMS::DataArrays::IntegerDataArray(self); }, "memo"_a)

        .def("get_data", [](const OpenMS::DataArrays::IntegerDataArray& self) {
            // Return as numpy array — bulk memcpy instead of element-by-element conversion
            size_t n = self.size();
            int32_t* data = new int32_t[n];
            std::copy(self.begin(), self.end(), data);
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<int32_t*>(p); });
            return nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(data, {n}, owner);
        }, "Returns a copy of the data as numpy array")

        .def("get_data_view", [](nb::object self_obj) {
            // Return a writable view (zero-copy), empty array if empty
            auto& self = nb::cast<OpenMS::DataArrays::IntegerDataArray&>(self_obj);
            int32_t* data_ptr = self.empty() ? nullptr : self.data();
            return nb::ndarray<nb::numpy, int32_t, nb::ndim<1>>(
                data_ptr, {self.size()}, self_obj
            );
        }, "Returns a zero-copy writable view of the data as numpy array (empty array if empty)")

        .def("set_data", [](OpenMS::DataArrays::IntegerDataArray& self, nb::object data_obj) {
            // Fast path: int32 numpy array — direct memcpy
            nb::ndarray<nb::numpy, int32_t, nb::ndim<1>> i32_arr;
            if (nb::try_cast(data_obj, i32_arr)) {
                const size_t n = i32_arr.shape(0);
                self.resize(n);
                std::copy(static_cast<const int32_t*>(i32_arr.data()),
                          static_cast<const int32_t*>(i32_arr.data()) + n, self.data());
                return;
            }
            // Fast path: int64 numpy array — single-pass narrowing copy
            nb::ndarray<nb::numpy, int64_t, nb::ndim<1>> i64_arr;
            if (nb::try_cast(data_obj, i64_arr)) {
                const size_t n = i64_arr.shape(0);
                self.resize(n);
                const int64_t* ptr = static_cast<const int64_t*>(i64_arr.data());
                for (size_t i = 0; i < n; ++i) {
                    if (ptr[i] > std::numeric_limits<OpenMS::Int>::max() || ptr[i] < std::numeric_limits<OpenMS::Int>::min()) {
                        throw std::overflow_error("Integer value at index " + std::to_string(i) + " overflows Int32");
                    }
                    self[i] = static_cast<OpenMS::Int>(ptr[i]);
                }
                return;
            }
            // Fallback: any iterable (list, etc.)
            std::vector<OpenMS::Int> vec = nb::cast<std::vector<OpenMS::Int>>(data_obj);
            self.assign(vec.begin(), vec.end());
        }, "data"_a, "Set data from a numpy array or list")

        .def("push_back", [](OpenMS::DataArrays::IntegerDataArray& self, OpenMS::Int val) {
            self.push_back(val);
        }, "val"_a, "Add a value to the array")

        .def("resize", [](OpenMS::DataArrays::IntegerDataArray& self, size_t n) {
            self.resize(n);
        }, "n"_a, "Resize the array")

        .def("clear", [](OpenMS::DataArrays::IntegerDataArray& self) {
            self.clear();
        }, "Clear the array")
        .def("reserve", [](OpenMS::DataArrays::IntegerDataArray& self, size_t n) { self.reserve(n); }, "n"_a, "Reserve memory for n elements")
        .def("getDataProcessing", [](const OpenMS::DataArrays::IntegerDataArray& self) { return self.getDataProcessing(); }, "Returns the data processing objects")
        .def("__repr__", [](const OpenMS::DataArrays::IntegerDataArray& self) {
            return "IntegerDataArray(name='" + self.getName() + "', size=" + std::to_string(self.size()) + ")";
        })
        .def("__str__", [](const OpenMS::DataArrays::IntegerDataArray& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;
    def_MetaInfoInterface<OpenMS::DataArrays::IntegerDataArray>(integerdataarray_class);

    // -----------------------------------------------------------------------
    // MobilityPeak1D
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MobilityPeak1D>(m, "MobilityPeak1D", 
        R"doc(
A 1-dimensional raw data mobility point or peak. The unit (ms, 1/K_0,
etc) is implicit
)doc")
        .def(nb::init<>())
        .def(nb::init<OpenMS::DPosition<1>, float>())
        .def(nb::init<const OpenMS::MobilityPeak1D &>())
        .def("__copy__", [](const OpenMS::MobilityPeak1D& self) { return OpenMS::MobilityPeak1D(self); })
        .def("__deepcopy__", [](const OpenMS::MobilityPeak1D& self, nb::dict) { return OpenMS::MobilityPeak1D(self); }, "memo"_a)
        .def("getIntensity", [](const OpenMS::MobilityPeak1D& self) { return self.getIntensity(); })
        .def("setIntensity", [](OpenMS::MobilityPeak1D& self, float intensity) { return self.setIntensity(intensity); }, "intensity"_a)
        .def("getMobility", [](const OpenMS::MobilityPeak1D& self) { return self.getMobility(); })
        .def("setMobility", [](OpenMS::MobilityPeak1D& self, double mobility) { return self.setMobility(mobility); }, "mobility"_a)
        .def("getPos", [](const OpenMS::MobilityPeak1D& self) { return self.getPos(); })
        .def("setPos", [](OpenMS::MobilityPeak1D& self, double pos) { return self.setPos(pos); }, "pos"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("__hash__", [](const OpenMS::MobilityPeak1D& self) { return std::hash<OpenMS::MobilityPeak1D>{}(self); })
        .def("__repr__", [](const OpenMS::MobilityPeak1D& self) {
            std::ostringstream os;
            os << "MobilityPeak1D(mobility=" << self.getMobility() << ", intensity=" << self.getIntensity() << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::MobilityPeak1D& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;

    // -----------------------------------------------------------------------
    // Mobilogram
    // -----------------------------------------------------------------------
    
    // ABI guard for zero-copy access (sizeof + standard-layout at top of module)
    nb::class_<OpenMS::Mobilogram>(m, "Mobilogram",
                                   
        R"doc(
RangeManagerMobInt

The representation of a 1D ion mobilogram.
Raw data access is provided by `get_peaks`, `get_peaks_struct`, and `set_peaks`.
Iterations yields access to underlying peak objects but is slower
Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
Usage:
.. code-block:: python
rt = mobilogram.getRT()
drift_time_unit = mobilogram.getDriftTimeUnit()
mobility, intensities = mobilogram.get_peaks()
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Mobilogram &>())
        .def("__copy__", [](const OpenMS::Mobilogram& self) { return OpenMS::Mobilogram(self); })
        .def("__deepcopy__", [](const OpenMS::Mobilogram& self, nb::dict) { return OpenMS::Mobilogram(self); }, "memo"_a)
        .def("resize", [](OpenMS::Mobilogram& self, size_t new_size) { return self.resize(new_size); }, "new_size"_a, "Resize the peak array")
        .def("reserve", [](OpenMS::Mobilogram& self, size_t new_size) { return self.reserve(new_size); }, "new_size"_a)
        .def("updateRanges", [](OpenMS::Mobilogram& self) { return self.updateRanges(); })
        .def("clearRanges", [](OpenMS::Mobilogram& self) { return self.clearRanges(); }, "Clear all ranges")
        .def("getRT", [](const OpenMS::Mobilogram& self) { return self.getRT(); }, "Returns the retention time (in seconds)")
        .def("setRT", [](OpenMS::Mobilogram& self, double rt) { return self.setRT(rt); }, "rt"_a, "Sets the retention time (in seconds)")
        .def("getDriftTimeUnit", [](const OpenMS::Mobilogram& self) { return self.getDriftTimeUnit(); }, "Returns the ion mobility drift time unit")
        .def("getDriftTimeUnitAsString", [](const OpenMS::Mobilogram& self) { return self.getDriftTimeUnitAsString(); }, "Returns the ion mobility drift time unit as string")
        .def("setDriftTimeUnit", [](OpenMS::Mobilogram& self, OpenMS::DriftTimeUnit dt) { return self.setDriftTimeUnit(dt); }, "dt"_a, "Sets the ion mobility drift time unit")
        .def("sortByIntensity", [](OpenMS::Mobilogram& self, bool reverse) { return self.sortByIntensity(reverse); }, "reverse"_a = false)
        .def("sortByPosition", [](OpenMS::Mobilogram& self) { return self.sortByPosition(); }, 
            R"doc(
Lexicographically sorts the peaks by their intensity
Sorts the peaks according to ascending intensity. Meta data arrays will be sorted accordingly
)doc")
        .def("isSorted", [](const OpenMS::Mobilogram& self) { return self.isSorted(); }, "Checks if all peaks are sorted with respect to ascending mobility")
        .def("calculateTIC", [](const OpenMS::Mobilogram& self) { return self.calculateTIC(); }, "Compute the total ion count (sum of all peak intensities)")
        .def("__iter__", [](OpenMS::Mobilogram& self) { return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::Mobilogram>(), "Mobilogram_iter", self.begin(), self.end()); })
        .def("__len__", [](OpenMS::Mobilogram& self) { return self.size(); })
        .def("__getitem__", [](OpenMS::Mobilogram& self, size_t i) -> OpenMS::MobilityPeak1D& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)
        .def("__setitem__", [](OpenMS::Mobilogram& self, size_t i, const OpenMS::MobilityPeak1D& val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        }, "i"_a, "val"_a, "Sets peak at index i")
        .def("findNearest", [](const OpenMS::Mobilogram& self, double mb) { return self.findNearest(mb); }, "mb"_a, "Returns the index of the closest peak in mobility")
        .def("findNearest", [](const OpenMS::Mobilogram& self, double mb, double tolerance) { return self.findNearest(mb, tolerance); }, "mb"_a, "tolerance"_a, "Returns the index of the closest peak in the provided tolerance window (-1 if none match)")
        .def("findNearest", [](const OpenMS::Mobilogram& self, double mb, double tolerance_left, double tolerance_right) { return self.findNearest(mb, tolerance_left, tolerance_right); }, "mb"_a, "tolerance_left"_a, "tolerance_right"_a, "Returns the index of the closest peak in the provided tolerance window (-1 if none match)")
        .def("getStringDataArrays", [](OpenMS::Mobilogram& self) -> OpenMS::Mobilogram::StringDataArrays& { return self.getStringDataArrays(); }, nb::rv_policy::reference_internal, "Returns a reference to the string data arrays")
        .def("setStringDataArrays", [](OpenMS::Mobilogram& self, const OpenMS::Mobilogram::StringDataArrays& sda) { return self.setStringDataArrays(sda); }, "sda"_a, "Sets the string data arrays")
        .def("getIntegerDataArrays", [](OpenMS::Mobilogram& self) -> OpenMS::Mobilogram::IntegerDataArrays& { return self.getIntegerDataArrays(); }, nb::rv_policy::reference_internal, "Returns a reference to the integer data arrays")
        .def("setIntegerDataArrays", [](OpenMS::Mobilogram& self, const OpenMS::Mobilogram::IntegerDataArrays& ida) { return self.setIntegerDataArrays(ida); }, "ida"_a, "Sets the integer data arrays")

        .def("size", [](const OpenMS::Mobilogram& self) {
            return self.size();
        }, "Returns the number of peaks")

        .def("push_back", [](OpenMS::Mobilogram& self, const OpenMS::MobilityPeak1D& p) {
            self.push_back(p);
        }, "peak"_a, "Add a peak")

        .def("clear", [](OpenMS::Mobilogram& self) {
            self.clear();
        }, "Clear all peaks")

        .def("get_peaks", [](const OpenMS::Mobilogram& self) {
            // Single allocation + single capsule to reduce overhead for small arrays
            const size_t n = self.size();
            const size_t mob_bytes = n * sizeof(double);
            const size_t int_bytes = n * sizeof(float);
            char* buf = new char[mob_bytes + int_bytes];
            double* mob_data = reinterpret_cast<double*>(buf);
            float* int_data = reinterpret_cast<float*>(buf + mob_bytes);
            for (size_t i = 0; i < n; ++i) {
                mob_data[i] = self[i].getMobility();
                int_data[i] = self[i].getIntensity();
            }
            nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<char*>(p); });
            auto mob_arr = nb::ndarray<nb::numpy, double, nb::ndim<1>>(mob_data, {n}, owner);
            auto int_arr = nb::ndarray<nb::numpy, float, nb::ndim<1>>(int_data, {n}, owner);
            return nb::make_tuple(mob_arr, int_arr);
        }, "Get mobility and intensity arrays as numpy arrays")

        .def("_get_peaks_view", [](nb::object self_obj) {
            auto& self = nb::cast<OpenMS::Mobilogram&>(self_obj);
            uint8_t* data_ptr = self.empty() ? nullptr : reinterpret_cast<uint8_t*>(&self[0]);
            size_t shape[1] = { self.size() * sizeof(OpenMS::MobilityPeak1D) };
            return nb::ndarray<nb::numpy, uint8_t, nb::c_contig>(
                data_ptr,
                1,
                shape,
                self_obj
            );
        },
        nb::rv_policy::reference_internal,
        "Returns a raw byte view of the underlying MobilityPeak1D array (AoS layout).")

        .def("get_peaks_struct",
            [](nb::object self_obj) -> nb::object {
                auto& self = nb::cast<OpenMS::Mobilogram&>(self_obj);
                size_t n = self.size();
                auto np = nb::module_::import_("numpy");
                nb::dict dtype_dict;
                // Derive dtype from C++ layout (validated by static_asserts at module init)
                constexpr size_t pos_offset = 0; // standard-layout: first member at offset 0
                constexpr size_t int_offset = sizeof(OpenMS::MobilityPeak1D::PositionType);
                dtype_dict["names"] = nb::make_tuple("mobility", "intensity");
                dtype_dict["formats"] = nb::make_tuple(np.attr("float64"), np.attr("float32"));
                dtype_dict["offsets"] = nb::make_tuple(pos_offset, int_offset);
                dtype_dict["itemsize"] = sizeof(OpenMS::MobilityPeak1D);
                auto py_dtype = np.attr("dtype")(dtype_dict);
                if (n == 0) {
                    return np.attr("empty")(0, py_dtype);
                }
                uint8_t* data_ptr = reinterpret_cast<uint8_t*>(&self[0]);
                size_t byte_shape[1] = { n * sizeof(OpenMS::MobilityPeak1D) };
                auto raw = nb::ndarray<nb::numpy, uint8_t, nb::c_contig>(
                    data_ptr, 1, byte_shape, self_obj
                );
                return np.attr("frombuffer")(raw, py_dtype);
            },
            nb::rv_policy::reference_internal,
            "Returns zero-copy structured array with fields 'mobility' (float64) and 'intensity' (float32)."
        )

        .def("set_peaks", [](OpenMS::Mobilogram& self, nb::object mob_obj, nb::object int_obj) {
            auto mob_arr = as_numpy_array<double>(mob_obj);
            auto int_arr = as_numpy_array<float>(int_obj);
            size_t n = mob_arr.shape(0);
            if (int_arr.shape(0) != n) throw std::runtime_error("Mobility and intensity arrays must have the same length");
            self.resize(n);
            const double* mob_ptr = static_cast<const double*>(mob_arr.data());
            const float* int_ptr = static_cast<const float*>(int_arr.data());
            for (size_t i = 0; i < n; ++i) {
                self[i].setMobility(mob_ptr[i]);
                self[i].setIntensity(int_ptr[i]);
            }
        }, "mob"_a, "intensity"_a, "Set mobility and intensity from numpy arrays")
        .def("set_peaks", [](OpenMS::Mobilogram& self, nb::object peaks_seq) {
            nb::object item0, item1;
            if (nb::isinstance<nb::tuple>(peaks_seq)) {
                auto tup = nb::cast<nb::tuple>(peaks_seq);
                if (nb::len(tup) != 2) throw std::runtime_error("set_peaks sequence must contain exactly 2 arrays");
                item0 = tup[0]; item1 = tup[1];
            } else {
                auto lst = nb::cast<nb::list>(peaks_seq);
                if (nb::len(lst) != 2) throw std::runtime_error("set_peaks sequence must contain exactly 2 arrays");
                item0 = lst[0]; item1 = lst[1];
            }
            auto mob_arr = as_numpy_array<double>(item0);
            auto int_arr = as_numpy_array<float>(item1);
            size_t n = mob_arr.shape(0);
            if (int_arr.shape(0) != n) throw std::runtime_error("Mobility and intensity arrays must have the same length");
            self.resize(n);
            const double* mob_ptr = static_cast<const double*>(mob_arr.data());
            const float* int_ptr = static_cast<const float*>(int_arr.data());
            for (size_t i = 0; i < n; ++i) {
                self[i].setMobility(mob_ptr[i]);
                self[i].setIntensity(int_ptr[i]);
            }
        }, "peaks"_a, "Set peaks from [mobility_array, intensity_array]")

        .def("getFloatDataArrays", [](const OpenMS::Mobilogram& self) {
            return self.getFloatDataArrays();
        }, "Get float data arrays")

        .def("setFloatDataArrays", [](OpenMS::Mobilogram& self, std::vector<OpenMS::DataArrays::FloatDataArray> fda) {
            self.setFloatDataArrays(fda);
        }, "fda"_a, "Set float data arrays")

        .def("getMinMobility", [](const OpenMS::Mobilogram& self) -> double {
            if (self.empty()) return 0.0;
            double min_mob = self[0].getMobility();
            for (size_t i = 1; i < self.size(); ++i) {
                if (self[i].getMobility() < min_mob) min_mob = self[i].getMobility();
            }
            return min_mob;
        }, "Get minimum mobility value")

        .def("getMaxMobility", [](const OpenMS::Mobilogram& self) -> double {
            if (self.empty()) return 0.0;
            double max_mob = self[0].getMobility();
            for (size_t i = 1; i < self.size(); ++i) {
                if (self[i].getMobility() > max_mob) max_mob = self[i].getMobility();
            }
            return max_mob;
        }, "Get maximum mobility value")

        .def("getMinIntensity", [](const OpenMS::Mobilogram& self) -> double {
            if (self.empty()) return 0.0;
            double min_int = self[0].getIntensity();
            for (size_t i = 1; i < self.size(); ++i) {
                if (self[i].getIntensity() < min_int) min_int = self[i].getIntensity();
            }
            return min_int;
        }, "Get minimum intensity value")

        .def("getMaxIntensity", [](const OpenMS::Mobilogram& self) -> double {
            if (self.empty()) return 0.0;
            double max_int = self[0].getIntensity();
            for (size_t i = 1; i < self.size(); ++i) {
                if (self[i].getIntensity() > max_int) max_int = self[i].getIntensity();
            }
            return max_int;
        }, "Get maximum intensity value")
        .def("__repr__", [](const OpenMS::Mobilogram& self) {
            std::ostringstream os;
            os << "Mobilogram(num_peaks=" << self.size() << ", rt=" << self.getRT() << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::Mobilogram& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;

    // -----------------------------------------------------------------------
    // OnDiscMSExperiment
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OnDiscMSExperiment>(m, "OnDiscMSExperiment", "Representation of a mass spectrometry experiment on disk.")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::OnDiscMSExperiment &>())
        .def("__copy__", [](const OpenMS::OnDiscMSExperiment& self) { return OpenMS::OnDiscMSExperiment(self); })
        .def("__deepcopy__", [](const OpenMS::OnDiscMSExperiment& self, nb::dict) { return OpenMS::OnDiscMSExperiment(self); }, "memo"_a)
        .def("openFile", [](OpenMS::OnDiscMSExperiment& self, const OpenMS::String& filename, bool skipMetaData) { return self.openFile(filename, skipMetaData); }, "filename"_a, "skipMetaData"_a = false)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getNrSpectra", [](const OpenMS::OnDiscMSExperiment& self) { return self.getNrSpectra(); }, "Returns the total number of spectra available")
        .def("getNrChromatograms", [](const OpenMS::OnDiscMSExperiment& self) { return self.getNrChromatograms(); }, "Returns the total number of chromatograms available")
        .def("getSpectrum", [](OpenMS::OnDiscMSExperiment& self, size_t id) { return self.getSpectrum(id); }, "id"_a)
        .def("getChromatogram", [](OpenMS::OnDiscMSExperiment& self, size_t id) { return self.getChromatogram(id); }, "id"_a, 
            R"doc(
Returns a single spectrum
:param id: The native identifier of the spectrum
)doc")
        .def("getChromatogramByNativeId", [](OpenMS::OnDiscMSExperiment& self, const OpenMS::String& id) { return self.getChromatogramByNativeId(id); }, "id"_a, 
            R"doc(
Returns a single chromatogram
:param id: The index of the chromatogram
)doc")
        .def("getSpectrumByNativeId", [](OpenMS::OnDiscMSExperiment& self, const OpenMS::String& id) { return self.getSpectrumByNativeId(id); }, "id"_a, 
            R"doc(
Returns a single spectrum
:param id: The index of the spectrum
)doc")
        .def("setSkipXMLChecks", [](OpenMS::OnDiscMSExperiment& self, bool skip) { return self.setSkipXMLChecks(skip); }, "skip"_a, "Sets whether to skip some XML checks and be fast instead")
        .def("getOptions", [](OpenMS::OnDiscMSExperiment& self) -> OpenMS::PeakFileOptions & { return self.getOptions(); }, nb::rv_policy::reference_internal, "Returns the options for loading/storing")
        .def("setOptions", [](OpenMS::OnDiscMSExperiment& self, const OpenMS::PeakFileOptions& options) { return self.setOptions(options); }, "options"_a, "Sets the options for loading/storing")
        .def("getExperimentalSettings", [](const OpenMS::OnDiscMSExperiment& self) { return self.getExperimentalSettings(); },
            R"doc(Returns the meta information of this experiment (as shared pointer to ExperimentalSettings).
)doc")
        .def("getMetaData", [](const OpenMS::OnDiscMSExperiment& self) { return self.getMetaData(); },
            R"doc(Returns the meta data as a shared pointer to a PeakMap (MSExperiment).

The returned MSExperiment contains the meta data (spectra headers, chromatogram headers,
experimental settings) but no actual peak data.

:return: Shared pointer to the meta data MSExperiment
)doc")
        .def("size", [](const OpenMS::OnDiscMSExperiment& self) { return self.size(); }, "Returns the number of spectra available")
        .def("empty", [](const OpenMS::OnDiscMSExperiment& self) { return self.empty(); }, "Returns whether the experiment contains no spectra")
        .def("isSortedByRT", [](const OpenMS::OnDiscMSExperiment& self) { return self.isSortedByRT(); }, "Checks if all spectra are sorted with respect to ascending RT")
        // Aliases for Cython API compatibility: getSpectrumById/getChromatogramById take a native ID string
        .def("getSpectrumById", [](OpenMS::OnDiscMSExperiment& self, const OpenMS::String& id) { return self.getSpectrumByNativeId(id); }, "id"_a,
            R"doc(
Returns a single spectrum by its native ID string.
This is an alias for getSpectrumByNativeId().
:param id: The native identifier of the spectrum
)doc")
        .def("getChromatogramById", [](OpenMS::OnDiscMSExperiment& self, const OpenMS::String& id) { return self.getChromatogramByNativeId(id); }, "id"_a,
            R"doc(
Returns a single chromatogram by its native ID string.
This is an alias for getChromatogramByNativeId().
:param id: The native identifier of the chromatogram
)doc")
        ;

    // -----------------------------------------------------------------------
    // Peak1D
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Peak1D>(m, "Peak1D", 
        R"doc(
A 1-dimensional raw data point or peak.
This data structure is intended for continuous data or peak data.
If you want to annotate single peaks with meta data, use RichPeak1D instead.
)doc")
        .def(nb::init<>())
        .def(nb::init<OpenMS::DPosition<1>, float>())
        .def(nb::init<const OpenMS::Peak1D &>())
        .def("__copy__", [](const OpenMS::Peak1D& self) { return OpenMS::Peak1D(self); })
        .def("__deepcopy__", [](const OpenMS::Peak1D& self, nb::dict) { return OpenMS::Peak1D(self); }, "memo"_a)
        .def("getIntensity", [](const OpenMS::Peak1D& self) { return self.getIntensity(); }, "Returns the intensity (height) of the peak")
        .def("setIntensity", [](OpenMS::Peak1D& self, float intensity) { return self.setIntensity(intensity); }, "intensity"_a, "Sets the intensity (height) of the peak")
        .def("getMZ", [](const OpenMS::Peak1D& self) { return self.getMZ(); }, "Returns the m/z (mass-to-charge) value of the peak")
        .def("setMZ", [](OpenMS::Peak1D& self, double mz) { return self.setMZ(mz); }, "mz"_a, "Sets the m/z (mass-to-charge) value of the peak")
        .def("getPos", [](const OpenMS::Peak1D& self) { return self.getPos(); }, "Returns the position (alias for getMZ)")
        .def("setPos", [](OpenMS::Peak1D& self, double pos) { return self.setPos(pos); }, "pos"_a, "Sets the position (alias for setMZ)")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("__hash__", [](const OpenMS::Peak1D& self) {
            // Content-based hash using mz and intensity
            size_t h1 = std::hash<double>{}(self.getMZ());
            size_t h2 = std::hash<float>{}(self.getIntensity());
            return h1 ^ (h2 << 1);
        })
        .def("__repr__", [](const OpenMS::Peak1D& self) {
            std::ostringstream os;
            os << std::fixed << std::setprecision(4) << "Peak1D(mz=" << self.getMZ() << ", intensity=" << std::setprecision(2) << self.getIntensity() << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::Peak1D& self) {
            std::ostringstream os;
            os << std::fixed << std::setprecision(4) << "(" << self.getMZ() << ", " << std::setprecision(2) << self.getIntensity() << ")";
            return os.str();
        })
        ;

    // -----------------------------------------------------------------------
    // Peak2D
    // -----------------------------------------------------------------------
    auto peak2d_class = nb::class_<OpenMS::Peak2D>(m, "Peak2D", 
        R"doc(
A 2-dimensional raw data point or peak.
This data structure is intended for continuous data or peak data.
If you want to annotated single peaks with meta data, use RichPeak2D instead
)doc")
        .def(nb::init<>())
        .def(nb::init<OpenMS::DPosition<2>, float>())
        .def(nb::init<const OpenMS::Peak2D &>())
        .def("__copy__", [](const OpenMS::Peak2D& self) { return OpenMS::Peak2D(self); })
        .def("__deepcopy__", [](const OpenMS::Peak2D& self, nb::dict) { return OpenMS::Peak2D(self); }, "memo"_a)
        .def("getIntensity", [](const OpenMS::Peak2D& self) { return self.getIntensity(); }, "Returns the data point intensity (height)")
        .def("setIntensity", [](OpenMS::Peak2D& self, float intensity) { return self.setIntensity(intensity); }, "intensity"_a, "Sets the data point intensity (height)")
        .def("getMZ", [](const OpenMS::Peak2D& self) { return self.getMZ(); }, "Returns the m/z coordinate (index 1)")
        .def("setMZ", [](OpenMS::Peak2D& self, double coordinate) { return self.setMZ(coordinate); }, "coordinate"_a, "Sets the m/z coordinate (index 1)")
        .def("getRT", [](const OpenMS::Peak2D& self) { return self.getRT(); }, "Returns the RT coordinate (index 0)")
        .def("setRT", [](OpenMS::Peak2D& self, double coordinate) { return self.setRT(coordinate); }, "coordinate"_a, "Sets the RT coordinate (index 0)")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("__hash__", [](const OpenMS::Peak2D& self) {
            // Content-based hash using mz, rt, and intensity
            size_t h1 = std::hash<double>{}(self.getMZ());
            size_t h2 = std::hash<double>{}(self.getRT());
            size_t h3 = std::hash<float>{}(self.getIntensity());
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        })
        .def("__repr__", [](const OpenMS::Peak2D& self) {
            std::ostringstream os;
            os << "Peak2D(rt=" << self.getRT() << ", mz=" << self.getMZ() << ", intensity=" << self.getIntensity() << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::Peak2D& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;
    // DimensionDescription enum nested under Peak2D
    nb::enum_<OpenMS::Peak2D::DimensionDescription>(peak2d_class, "DimensionDescription", "Names of the dimensions of a 2D peak")
        .value("RT", OpenMS::Peak2D::DimensionDescription::RT)
        .value("MZ", OpenMS::Peak2D::DimensionDescription::MZ)
        .value("DIMENSION", OpenMS::Peak2D::DimensionDescription::DIMENSION)

        .export_values();

    // -----------------------------------------------------------------------
    // PeakIndex
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakIndex>(m, "PeakIndex", 
        R"doc(
Index of a peak or feature
This struct can be used to store both peak or feature indices
)doc")
        .def(nb::init<>())
        .def(nb::init<size_t>())
        .def(nb::init<size_t, size_t>())
        .def("isValid", [](const OpenMS::PeakIndex& self) { return self.isValid(); }, "Returns if the current peak ref is valid")
        .def("clear", [](OpenMS::PeakIndex& self) { return self.clear(); }, "Invalidates the current index")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def_rw("peak", &OpenMS::PeakIndex::peak)
        .def_rw("spectrum", &OpenMS::PeakIndex::spectrum)
        .def("getFeature", [](const OpenMS::PeakIndex& self, const OpenMS::FeatureMap& map) -> const OpenMS::Feature& { return self.getFeature(map); }, "map"_a, nb::rv_policy::reference_internal, "Returns the feature in the given map")
        .def("getPeak", [](const OpenMS::PeakIndex& self, const OpenMS::MSExperiment& map) -> const OpenMS::Peak1D& { return self.getPeak(map); }, "map"_a, nb::rv_policy::reference_internal, "Returns the peak in the given map")
        .def("getSpectrum", [](const OpenMS::PeakIndex& self, const OpenMS::MSExperiment& map) -> const OpenMS::MSSpectrum& { return self.getSpectrum(map); }, "map"_a, nb::rv_policy::reference_internal, "Returns the spectrum in the given map")
        ;

    // -----------------------------------------------------------------------
    // PeptideHit
    // -----------------------------------------------------------------------
    auto peptidehit_class = nb::class_<OpenMS::PeptideHit>(m, "PeptideHit",
        R"doc(
MetaInfoInterface

Represents a single peptide identification hit from a database search
A PeptideHit stores information about a candidate peptide sequence that was
matched to a spectrum. Each hit contains:
- The peptide sequence (as AASequence)
- A score from the search engine
- The rank among all candidates
- The charge state
- Protein mappings (PeptideEvidence objects)
Multiple PeptideHit objects are typically stored in a PeptideIdentification,
sorted by score to show the most likely candidates first.
Example usage:
.. code-block:: python
hit = oms.PeptideHit()
hit.setSequence(oms.AASequence.fromString("PEPTIDER"))
hit.setScore(95.5)
hit.setRank(1)
hit.setCharge(2)
# Access information
print(f"Sequence: {hit.getSequence().toString()}")
print(f"Score: {hit.getScore()}, Rank: {hit.getRank()}")
print(f"Charge: {hit.getCharge()}")
)doc")
        .def(nb::init<>())
        .def(nb::init<double, unsigned int, int, OpenMS::AASequence>())
        .def(nb::init<double, unsigned int, int, OpenMS::AASequence>())
        .def(nb::init<const OpenMS::PeptideHit &>())
        .def("__copy__", [](const OpenMS::PeptideHit& self) { return OpenMS::PeptideHit(self); })
        .def("__deepcopy__", [](const OpenMS::PeptideHit& self, nb::dict) { return OpenMS::PeptideHit(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getSequence", [](OpenMS::PeptideHit& self) -> OpenMS::AASequence & { return self.getSequence(); }, nb::rv_policy::reference_internal, 
            R"doc(
Returns the rank of this hit among all candidates
:return: Rank (1 = best hit, 2 = second best, etc.)
)doc")
        .def("setSequence", [](OpenMS::PeptideHit& self, const OpenMS::AASequence& sequence) { return self.setSequence(sequence); }, "sequence"_a, 
            R"doc(
Sets the rank of this hit
:param rank: Rank among all candidates (1 = best)
)doc")
        .def("setSequence", [](OpenMS::PeptideHit& self, OpenMS::AASequence& sequence) { return self.setSequence(sequence); }, "sequence"_a, 
            R"doc(
Sets the rank of this hit
:param rank: Rank among all candidates (1 = best)
)doc")
        .def("getCharge", [](const OpenMS::PeptideHit& self) { return self.getCharge(); }, 
            R"doc(
Returns the peptide sequence
:return: The peptide amino acid sequence with modifications
)doc")
        .def("setCharge", [](OpenMS::PeptideHit& self, int charge) { return self.setCharge(charge); }, "charge"_a, 
            R"doc(
Sets the peptide sequence
:param sequence: The peptide amino acid sequence
)doc")
        .def("getPeptideEvidences", [](const OpenMS::PeptideHit& self) -> const std::vector<OpenMS::PeptideEvidence> & { return self.getPeptideEvidences(); }, nb::rv_policy::reference_internal, 
            R"doc(
Returns the charge state of the peptide
:return: Charge state (e.g., 2 for doubly charged)
)doc")
        .def("setPeptideEvidences", [](OpenMS::PeptideHit& self, const std::vector<OpenMS::PeptideEvidence>& peptide_evidences) { return self.setPeptideEvidences(peptide_evidences); }, "peptide_evidences"_a, 
            R"doc(
Returns protein mapping information for this peptide
:return: List of proteins where this peptide was found
Each evidence contains protein accession, start/end positions, and if it's a decoy
)doc")
        .def("setPeptideEvidences", [](OpenMS::PeptideHit& self, std::vector<OpenMS::PeptideEvidence>& peptide_evidences) { return self.setPeptideEvidences(peptide_evidences); }, "peptide_evidences"_a, 
            R"doc(
Returns protein mapping information for this peptide
:return: List of proteins where this peptide was found
Each evidence contains protein accession, start/end positions, and if it's a decoy
)doc")
        .def("addPeptideEvidence", [](OpenMS::PeptideHit& self, const OpenMS::PeptideEvidence& peptide_evidence) { return self.addPeptideEvidence(peptide_evidence); }, "peptide_evidence"_a, 
            R"doc(
Sets the protein mapping information
:param evidences: Protein locations for this peptide
)doc")
        .def("getScore", [](const OpenMS::PeptideHit& self) { return self.getScore(); })
        .def("setScore", [](OpenMS::PeptideHit& self, double score) { return self.setScore(score); }, "score"_a, 
            R"doc(
Returns fragment ion annotations
:return: Annotated fragment peaks
)doc")
        .def("setAnalysisResults", [](OpenMS::PeptideHit& self, const std::vector<OpenMS::PeptideHit::PepXMLAnalysisResult>& aresults) { self.setAnalysisResults(aresults); }, "aresults"_a,
            R"doc(
Set information on (search engine) sub scores associated with this PSM.
:param aresults: List of PepXMLAnalysisResult objects
)doc")
        .def("addAnalysisResults", [](OpenMS::PeptideHit& self, const OpenMS::PeptideHit::PepXMLAnalysisResult& aresult) { return self.addAnalysisResults(aresult); }, "aresult"_a,
            R"doc(
Sets search engine sub-scores
:param aresult: Sub-score information from search engine
)doc")
        .def("getAnalysisResults", [](const OpenMS::PeptideHit& self) { return self.getAnalysisResults(); },
            R"doc(
Adds a search engine sub-score
:param aresult: Sub-score to add
)doc")
        .def("getRank", [](const OpenMS::PeptideHit& self) { return self.getRank(); }, 
            R"doc(
Returns the score of this peptide-spectrum match (PSM)
:return: The search engine score
Interpretation depends on the score type (check isHigherScoreBetter)
)doc")
        .def("setRank", [](OpenMS::PeptideHit& self, unsigned int newrank) { return self.setRank(newrank); }, "newrank"_a, 
            R"doc(
Sets the PSM score
:param score: The search engine score to set
)doc")
        .def("getPeakAnnotations", [](OpenMS::PeptideHit& self) -> std::vector<OpenMS::PeptideHit::PeakAnnotation> & { return self.getPeakAnnotations(); }, nb::rv_policy::reference_internal, 
            R"doc(
Sets fragment ion annotations
:param annotations: Fragment peak annotations
)doc")
        .def("setPeakAnnotations", [](OpenMS::PeptideHit& self, std::vector<OpenMS::PeptideHit::PeakAnnotation> frag_annotations) { return self.setPeakAnnotations(frag_annotations); }, "frag_annotations"_a, 
            R"doc(
Returns all search engine sub-scores
:return: Sub-score information
)doc")
        .def("isDecoy", [](const OpenMS::PeptideHit& self) { return self.isDecoy(); }, 
            R"doc(
Extracts all unique protein accessions
:return: Set of unique protein accession strings
Empty accessions are excluded from the result
)doc")
        .def("extractProteinAccessionsSet", [](const OpenMS::PeptideHit& self) { return self.extractProteinAccessionsSet(); }, 
            R"doc(
Adds a single protein mapping
:param evidence: Protein location information to add
)doc")
        
        .def("__hash__", [](const OpenMS::PeptideHit& self) { return std::hash<OpenMS::PeptideHit>{}(self); })
        .def("__repr__", [](const OpenMS::PeptideHit& self) {
            std::ostringstream os;
            os << "PeptideHit(score=" << self.getScore()
               << ", sequence='" << self.getSequence().toString() << "'"
               << ", charge=" << self.getCharge();
            const auto& evs = self.getPeptideEvidences();
            if (!evs.empty()) {
                os << ", evidences=[";
                for (size_t i = 0; i < evs.size(); ++i) {
                    if (i > 0) os << ", ";
                    os << "PeptideEvidence(accession='" << evs[i].getProteinAccession() << "')";
                }
                os << "]";
            }
            os << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::PeptideHit& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;
    def_MetaInfoInterface<OpenMS::PeptideHit>(peptidehit_class);

    // -----------------------------------------------------------------------
    // PeptideIdentification
    // -----------------------------------------------------------------------
    auto peptideidentification_class = nb::class_<OpenMS::PeptideIdentification>(m, "PeptideIdentification",
        R"doc(
MetaInfoInterface

Represents peptide identification results for a single spectrum or feature
PeptideIdentification stores the results of peptide identification from database
search engines (e.g., Mascot, X!Tandem, MSGF+). Each PeptideIdentification contains:
- A list of peptide hits (candidate sequences) ranked by score
- The precursor m/z and retention time
- Score type and significance threshold
- Link to the ProteinIdentification (via identifier)
Multiple PeptideIdentifications can belong to one ProteinIdentification, which
stores the search parameters and protein-level results.
Example usage:
.. code-block:: python
pep_id = oms.PeptideIdentification()
pep_id.setRT(1234.5)  # Set retention time
pep_id.setMZ(445.678)  # Set precursor m/z
pep_id.setScoreType("XTandem")
# Add a peptide hit
hit = oms.PeptideHit()
hit.setScore(50.5)
hit.setRank(1)
hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
hit.setCharge(2)
pep_id.insertHit(hit)
# Access hits
for hit in pep_id.getHits():
print(f"Sequence: {hit.getSequence().toString()}, Score: {hit.getScore()}")
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PeptideIdentification &>())
        .def("__copy__", [](const OpenMS::PeptideIdentification& self) { return OpenMS::PeptideIdentification(self); })
        .def("__deepcopy__", [](const OpenMS::PeptideIdentification& self, nb::dict) { return OpenMS::PeptideIdentification(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getRT", [](const OpenMS::PeptideIdentification& self) { return self.getRT(); }, 
            R"doc(
Checks if retention time is set
:return: True if RT is available, False otherwise
)doc")
        .def("setRT", [](OpenMS::PeptideIdentification& self, double rt) { return self.setRT(rt); }, "rt"_a, 
            R"doc(
Returns the retention time of the precursor
:return: Retention time in seconds
)doc")
        .def("hasRT", [](const OpenMS::PeptideIdentification& self) { return self.hasRT(); }, 
            R"doc(
Sets the precursor m/z value
:param mz: Mass-to-charge ratio of the precursor
)doc")
        .def("getMZ", [](const OpenMS::PeptideIdentification& self) { return self.getMZ(); }, 
            R"doc(
Checks if m/z value is set
:return: True if m/z is available, False otherwise
)doc")
        .def("setMZ", [](OpenMS::PeptideIdentification& self, double mz) { return self.setMZ(mz); }, "mz"_a, 
            R"doc(
Returns the precursor m/z value
:return: Mass-to-charge ratio of the precursor ion
)doc")
        .def("hasMZ", [](const OpenMS::PeptideIdentification& self) { return self.hasMZ(); }, 
            R"doc(
Set the spectrum reference (native ID) for this identification.
:param ref: Spectrum reference string (native ID)
)doc")
        .def("getHits", [](OpenMS::PeptideIdentification& self) -> std::vector<OpenMS::PeptideHit> & { return self.getHits(); }, nb::rv_policy::reference_internal)
        .def("insertHit", [](OpenMS::PeptideIdentification& self, const OpenMS::PeptideHit& hit) { return self.insertHit(hit); }, "hit"_a, 
            R"doc(
Returns all peptide hits (candidate sequences)
:return: List of peptide candidates ranked by score
Hits are typically sorted by score, with the best hit at index 0
)doc")
        .def("insertHit", [](OpenMS::PeptideIdentification& self, OpenMS::PeptideHit& hit) { return self.insertHit(hit); }, "hit"_a, 
            R"doc(
Returns all peptide hits (candidate sequences)
:return: List of peptide candidates ranked by score
Hits are typically sorted by score, with the best hit at index 0
)doc")
        .def("setHits", [](OpenMS::PeptideIdentification& self, const std::vector<OpenMS::PeptideHit>& hits) { return self.setHits(hits); }, "hits"_a, 
            R"doc(
Appends a peptide hit to the list
:param hit: The peptide hit to add
)doc")
        .def("setHits", [](OpenMS::PeptideIdentification& self, std::vector<OpenMS::PeptideHit>& hits) { return self.setHits(hits); }, "hits"_a, 
            R"doc(
Appends a peptide hit to the list
:param hit: The peptide hit to add
)doc")
        .def("getSignificanceThreshold", [](const OpenMS::PeptideIdentification& self) { return self.getSignificanceThreshold(); }, 
            R"doc(
Sets all peptide hits at once
:param hits: List of peptide hits to store
)doc")
        .def("setSignificanceThreshold", [](OpenMS::PeptideIdentification& self, double value) { return self.setSignificanceThreshold(value); }, "value"_a, 
            R"doc(
Returns the significance threshold value
:return: The threshold value (interpretation depends on score type)
Hits with scores below/above this threshold (depending on score direction) may be considered insignificant
)doc")
        .def("getScoreType", [](const OpenMS::PeptideIdentification& self) { return self.getScoreType(); }, 
            R"doc(
Sets the significance threshold value
:param value: The threshold value to set
)doc")
        .def("setScoreType", [](OpenMS::PeptideIdentification& self, const OpenMS::String& type) { return self.setScoreType(type); }, "type"_a, 
            R"doc(
Returns the type of score (e.g., "Mascot", "XTandem", "q-value")
:return: Name of the score type
)doc")
        .def("isHigherScoreBetter", [](const OpenMS::PeptideIdentification& self) { return self.isHigherScoreBetter(); }, 
            R"doc(
Sets the score type
:param score_type: Name of the score type (e.g., "Mascot", "XTandem")
)doc")
        .def("setHigherScoreBetter", [](OpenMS::PeptideIdentification& self, bool value) { return self.setHigherScoreBetter(value); }, "value"_a, 
            R"doc(
Returns whether higher scores are better
:return: True if higher scores indicate better matches, False if lower is better
)doc")
        .def("getIdentifier", [](const OpenMS::PeptideIdentification& self) { return self.getIdentifier(); }, 
            R"doc(
Sets whether higher scores are better
:param higher_better: True if higher scores are better, False otherwise
)doc")
        .def("setIdentifier", [](OpenMS::PeptideIdentification& self, const OpenMS::String& id) { return self.setIdentifier(id); }, "id"_a, 
            R"doc(
Returns the identifier linking to the parent ProteinIdentification
:return: Unique identifier string
Use this to find the corresponding ProteinIdentification with search parameters
)doc")
        .def("getBaseName", [](const OpenMS::PeptideIdentification& self) { return self.getBaseName(); }, 
            R"doc(
Sets the retention time of the precursor
:param rt: Retention time in seconds
)doc")
        .def("setBaseName", [](OpenMS::PeptideIdentification& self, const OpenMS::String& base_name) { return self.setBaseName(base_name); }, "base_name"_a)
        .def("getExperimentLabel", [](const OpenMS::PeptideIdentification& self) { return self.getExperimentLabel(); })
        .def("setExperimentLabel", [](OpenMS::PeptideIdentification& self, const OpenMS::String& type) { return self.setExperimentLabel(type); }, "type"_a)
        .def("getSpectrumReference", [](const OpenMS::PeptideIdentification& self) { return self.getSpectrumReference(); }, 
            R"doc(
Sets the identifier linking to a ProteinIdentification
:param identifier: Unique identifier string
)doc")
        .def("setSpectrumReference", [](OpenMS::PeptideIdentification& self, const OpenMS::String& ref) { return self.setSpectrumReference(ref); }, "ref"_a, 
            R"doc(
Get the spectrum reference (native ID) for this identification.
:return: Spectrum reference string (native ID)
)doc")
        .def("sort", [](OpenMS::PeptideIdentification& self) { return self.sort(); })
        .def("empty", [](const OpenMS::PeptideIdentification& self) { return self.empty(); })
        .def_static("getReferencingHits", [](const std::vector<OpenMS::PeptideHit>& p0, const std::set<OpenMS::String>& accession) { return OpenMS::PeptideIdentification::getReferencingHits(p0, accession); }, "Returns all peptide hits which reference to a given protein accession (i.e. filter by protein accession)")
        .def("buildUSI", [](const OpenMS::PeptideIdentification& self, const OpenMS::String& ms_run_name, const OpenMS::String& dataset_id, bool include_interpretation) { return self.buildUSI(ms_run_name, dataset_id, include_interpretation); }, "ms_run_name"_a, "dataset_id"_a = "local", "include_interpretation"_a = false)
        .def("buildUSI", [](const OpenMS::PeptideIdentification& self, const OpenMS::IdentifierMSRunMapper& mapping, const OpenMS::String& dataset_id, bool include_interpretation) { return self.buildUSI(mapping, dataset_id, include_interpretation); }, "mapping"_a, "dataset_id"_a = "local", "include_interpretation"_a = false)
        
        .def("__hash__", [](const OpenMS::PeptideIdentification& self) { return std::hash<OpenMS::PeptideIdentification>{}(self); })
        .def("__repr__", [](const OpenMS::PeptideIdentification& self) {
            std::ostringstream os;
            os << "PeptideIdentification(rt=" << self.getRT()
               << ", mz=" << self.getMZ()
               << ", score_type='" << self.getScoreType() << "'"
               << ", num_hits=" << self.getHits().size();
            const auto& hits = self.getHits();
            if (!hits.empty()) {
                os << ", top_hit='" << hits[0].getSequence().toString() << "'";
            }
            os << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::PeptideIdentification& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;
    def_MetaInfoInterface<OpenMS::PeptideIdentification>(peptideidentification_class);

    // -----------------------------------------------------------------------
    // Precursor
    // -----------------------------------------------------------------------
    auto precursor_class = nb::class_<OpenMS::Precursor, OpenMS::CVTermList>(m, "Precursor", 
        R"doc(
Peak1D
CVTermList

Precursor meta information
This class contains precursor information:
- isolation window
- activation
- selected ion (m/z, intensity, charge, possible charge states)
- ion mobility drift time
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Precursor &>())
        .def("__copy__", [](const OpenMS::Precursor& self) { return OpenMS::Precursor(self); })
        .def("__deepcopy__", [](const OpenMS::Precursor& self, nb::dict) { return OpenMS::Precursor(self); }, "memo"_a)
        .def_static("getAllNamesOfActivationMethods", []() { return OpenMS::Precursor::getAllNamesOfActivationMethods(); }, 
            R"doc(
Returns the full names (e.g., "Collision-induced dissociation") of ALL possible activation methods, not just those set on this instance
)doc")
        .def_static("getAllShortNamesOfActivationMethods", []() { return OpenMS::Precursor::getAllShortNamesOfActivationMethods(); }, 
            R"doc(
Returns the abbreviations (e.g., "CID") of ALL possible activation methods, not just those set on this instance
)doc")
        .def_static("activationMethodToString", [](OpenMS::Precursor::ActivationMethod m) { return OpenMS::Precursor::activationMethodToString(m); }, "m"_a, "Convert an ActivationMethod enum to its full name string")
        .def_static("activationMethodToShortString", [](OpenMS::Precursor::ActivationMethod m) { return OpenMS::Precursor::activationMethodToShortString(m); }, "m"_a, "Convert an ActivationMethod enum to its short (abbreviated) name string")
        .def_static("toActivationMethod", [](const OpenMS::String& name) { return OpenMS::Precursor::toActivationMethod(name); }, "name"_a, "Convert a string (full name or short name) to an ActivationMethod enum")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getActivationMethods", [](OpenMS::Precursor& self) -> std::set<OpenMS::Precursor::ActivationMethod> & { return self.getActivationMethods(); }, nb::rv_policy::reference_internal, "Returns the activation methods")
        .def("getActivationMethodsAsString", [](const OpenMS::Precursor& self) { return self.getActivationMethodsAsString(); }, 
            R"doc(
Returns the full names (e.g., "Collision-induced dissociation") of the activation methods set on this instance
)doc")
        .def("getActivationMethodsAsShortString", [](const OpenMS::Precursor& self) { return self.getActivationMethodsAsShortString(); }, 
            R"doc(
Returns the abbreviations (e.g., "CID") of the activation methods set on this instance
)doc")
        .def("setActivationMethods", [](OpenMS::Precursor& self, const std::set<OpenMS::Precursor::ActivationMethod>& activation_methods) { return self.setActivationMethods(activation_methods); }, "activation_methods"_a, "Sets the activation methods")
        .def("getActivationEnergy", [](const OpenMS::Precursor& self) { return self.getActivationEnergy(); }, "Returns the activation energy (in electronvolt)")
        .def("setActivationEnergy", [](OpenMS::Precursor& self, double activation_energy) { return self.setActivationEnergy(activation_energy); }, "activation_energy"_a, "Sets the activation energy (in electronvolt)")
        .def("getIsolationWindowLowerOffset", [](const OpenMS::Precursor& self) { return self.getIsolationWindowLowerOffset(); }, "Returns the lower offset from the target m/z")
        .def("setIsolationWindowLowerOffset", [](OpenMS::Precursor& self, double bound) { return self.setIsolationWindowLowerOffset(bound); }, "bound"_a, "Sets the lower offset from the target m/z")
        .def("getIsolationWindowUpperOffset", [](const OpenMS::Precursor& self) { return self.getIsolationWindowUpperOffset(); }, "Returns the upper offset from the target m/z")
        .def("setIsolationWindowUpperOffset", [](OpenMS::Precursor& self, double bound) { return self.setIsolationWindowUpperOffset(bound); }, "bound"_a, "Sets the upper offset from the target m/z")
        .def("getDriftTime", [](const OpenMS::Precursor& self) { return self.getDriftTime(); }, "Returns the ion mobility drift time in milliseconds (-1 means it is not set)")
        .def("setDriftTime", [](OpenMS::Precursor& self, double drift_time) { return self.setDriftTime(drift_time); }, "drift_time"_a, "Sets the ion mobility drift time in milliseconds")
        .def("getDriftTimeWindowLowerOffset", [](const OpenMS::Precursor& self) { return self.getDriftTimeWindowLowerOffset(); }, "Returns the lower offset from the target ion mobility in milliseconds")
        .def("setDriftTimeWindowLowerOffset", [](OpenMS::Precursor& self, double drift_time) { return self.setDriftTimeWindowLowerOffset(drift_time); }, "drift_time"_a, "Sets the lower offset from the target ion mobility")
        .def("getDriftTimeWindowUpperOffset", [](const OpenMS::Precursor& self) { return self.getDriftTimeWindowUpperOffset(); }, "Returns the upper offset from the target ion mobility in milliseconds")
        .def("setDriftTimeWindowUpperOffset", [](OpenMS::Precursor& self, double drift_time) { return self.setDriftTimeWindowUpperOffset(drift_time); }, "drift_time"_a, "Sets the upper offset from the target ion mobility")
        .def("getCharge", [](const OpenMS::Precursor& self) { return self.getCharge(); }, "Returns the charge")
        .def("setCharge", [](OpenMS::Precursor& self, int charge) { return self.setCharge(charge); }, "charge"_a, "Sets the charge")
        .def("getPossibleChargeStates", [](OpenMS::Precursor& self) -> std::vector<int> & { return self.getPossibleChargeStates(); }, nb::rv_policy::reference_internal, "Returns the possible charge states")
        .def("setPossibleChargeStates", [](OpenMS::Precursor& self, const std::vector<int>& possible_charge_states) { return self.setPossibleChargeStates(possible_charge_states); }, "possible_charge_states"_a, "Sets the possible charge states")
        .def("getUnchargedMass", [](const OpenMS::Precursor& self) { return self.getUnchargedMass(); }, "Returns the uncharged mass of the precursor, if charge is unknown, i.e. 0 best guess is its doubly charged")
        
        .def("getIntensity", [](const OpenMS::Precursor& self) { return self.getIntensity(); }, "Returns the intensity (height) of the peak")
        .def("setIntensity", [](OpenMS::Precursor& self, float intensity) { return self.setIntensity(intensity); }, "intensity"_a, "Sets the intensity (height) of the peak")
        .def("getMZ", [](const OpenMS::Precursor& self) { return self.getMZ(); }, "Returns the m/z (mass-to-charge) value of the peak")
        .def("setMZ", [](OpenMS::Precursor& self, double mz) { return self.setMZ(mz); }, "mz"_a, "Sets the m/z (mass-to-charge) value of the peak")
        .def("getPos", [](const OpenMS::Precursor& self) { return self.getPos(); }, "Returns the position (alias for getMZ)")
        .def("setPos", [](OpenMS::Precursor& self, double pos) { return self.setPos(pos); }, "pos"_a, "Sets the position (alias for setMZ)")
        .def("getMZ", [](OpenMS::Precursor& self) { return self.getMZ(); }, "Returns the m/z")
        .def("setMZ", [](OpenMS::Precursor& self, double mz) { self.setMZ(mz); }, "mz"_a, "Sets the m/z")
        .def("getIntensity", [](OpenMS::Precursor& self) { return self.getIntensity(); }, "Returns the intensity")
        .def("setIntensity", [](OpenMS::Precursor& self, float intensity) { self.setIntensity(intensity); }, "intensity"_a, "Sets the intensity")
        .def("getPos", [](OpenMS::Precursor& self) { return self.getPos(); }, "Returns the position (m/z)")
        .def("setPos", [](OpenMS::Precursor& self, double pos) { self.setPos(pos); }, "pos"_a, "Sets the position (m/z)")
        .def("__hash__", [](const OpenMS::Precursor& self) { return std::hash<OpenMS::Precursor>{}(self); })
        ;
    // ActivationMethod enum nested under Precursor
    nb::enum_<OpenMS::Precursor::ActivationMethod>(precursor_class, "ActivationMethod")
        .value("CID", OpenMS::Precursor::ActivationMethod::CID)
        .value("PSD", OpenMS::Precursor::ActivationMethod::PSD)
        .value("PD", OpenMS::Precursor::ActivationMethod::PD)
        .value("SID", OpenMS::Precursor::ActivationMethod::SID)
        .value("BIRD", OpenMS::Precursor::ActivationMethod::BIRD)
        .value("ECD", OpenMS::Precursor::ActivationMethod::ECD)
        .value("IMD", OpenMS::Precursor::ActivationMethod::IMD)
        .value("SORI", OpenMS::Precursor::ActivationMethod::SORI)
        .value("HCID", OpenMS::Precursor::ActivationMethod::HCID)
        .value("LCID", OpenMS::Precursor::ActivationMethod::LCID)
        .value("PHD", OpenMS::Precursor::ActivationMethod::PHD)
        .value("ETD", OpenMS::Precursor::ActivationMethod::ETD)
        .value("ETciD", OpenMS::Precursor::ActivationMethod::ETciD)
        .value("EThcD", OpenMS::Precursor::ActivationMethod::EThcD)
        .value("PQD", OpenMS::Precursor::ActivationMethod::PQD)
        .value("TRAP", OpenMS::Precursor::ActivationMethod::TRAP)
        .value("HCD", OpenMS::Precursor::ActivationMethod::HCD)
        .value("INSOURCE", OpenMS::Precursor::ActivationMethod::INSOURCE)
        .value("LIFT", OpenMS::Precursor::ActivationMethod::LIFT)
        .value("SIZE_OF_ACTIVATIONMETHOD", OpenMS::Precursor::ActivationMethod::SIZE_OF_ACTIVATIONMETHOD)
        ;

    // -----------------------------------------------------------------------
    // Prediction
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedExperimentHelper::Prediction, OpenMS::CVTermList>(m, "Prediction", "OpenMS class Prediction")
        .def(nb::init<>())
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        
        .def_rw("software_ref", &OpenMS::TargetedExperimentHelper::Prediction::software_ref)
        .def_rw("contact_ref", &OpenMS::TargetedExperimentHelper::Prediction::contact_ref)
        ;

    // -----------------------------------------------------------------------
    // Product
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Product, OpenMS::CVTermList>(m, "Product", "This class describes the product isolation window for special scan types, such as MRM")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Product &>())
        .def("__copy__", [](const OpenMS::Product& self) { return OpenMS::Product(self); })
        .def("__deepcopy__", [](const OpenMS::Product& self, nb::dict) { return OpenMS::Product(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getMZ", [](const OpenMS::Product& self) { return self.getMZ(); }, "Returns the target m/z")
        .def("setMZ", [](OpenMS::Product& self, double mz) { return self.setMZ(mz); }, "mz"_a, "Sets the target m/z")
        .def("getIsolationWindowLowerOffset", [](const OpenMS::Product& self) { return self.getIsolationWindowLowerOffset(); }, "Returns the lower offset from the target m/z")
        .def("setIsolationWindowLowerOffset", [](OpenMS::Product& self, double bound) { return self.setIsolationWindowLowerOffset(bound); }, "bound"_a, "Sets the lower offset from the target m/z")
        .def("getIsolationWindowUpperOffset", [](const OpenMS::Product& self) { return self.getIsolationWindowUpperOffset(); }, "Returns the upper offset from the target m/z")
        .def("setIsolationWindowUpperOffset", [](OpenMS::Product& self, double bound) { return self.setIsolationWindowUpperOffset(bound); }, "bound"_a, "Sets the upper offset from the target m/z")
        .def("__hash__", [](const OpenMS::Product& self) { return std::hash<OpenMS::Product>{}(self); })
        ;

    // -----------------------------------------------------------------------
    // Protein
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedExperimentHelper::Protein, OpenMS::CVTermList>(m, "Protein", "OpenMS class Protein")
        .def(nb::init<>())
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        
        .def_rw("id", &OpenMS::TargetedExperimentHelper::Protein::id)
        .def_rw("sequence", &OpenMS::TargetedExperimentHelper::Protein::sequence)
        ;

    // -----------------------------------------------------------------------
    // ProteinHit
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProteinHit, OpenMS::MetaInfoInterface>(m, "ProteinHit", 
        R"doc(
MetaInfoInterface

Represents a single protein identification hit from a database search
A ProteinHit stores information about a protein that was identified based on
peptide evidence. Each hit contains:
- Protein accession (database identifier)
- Score from protein inference
- Rank among protein candidates
- Protein sequence (optional)
- Sequence coverage percentage
Multiple ProteinHit objects are stored in a ProteinIdentification, typically
sorted by score to show the most confident identifications first.
Example usage:
.. code-block:: python
protein_hit = oms.ProteinHit()
protein_hit.setAccession("P12345")
protein_hit.setScore(150.5)
protein_hit.setRank(1)
protein_hit.setCoverage(45.2)  # 45.2% coverage
protein_hit.setDescription("Example protein")
# Access information
print(f"Accession: {protein_hit.getAccession()}")
print(f"Score: {protein_hit.getScore()}, Coverage: {protein_hit.getCoverage()}%")
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProteinHit &>())
        .def(nb::init<double, unsigned int, OpenMS::String, OpenMS::String>())
        .def("__copy__", [](const OpenMS::ProteinHit& self) { return OpenMS::ProteinHit(self); })
        .def("__deepcopy__", [](const OpenMS::ProteinHit& self, nb::dict) { return OpenMS::ProteinHit(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getScore", [](const OpenMS::ProteinHit& self) { return self.getScore(); })
        .def("getRank", [](const OpenMS::ProteinHit& self) { return self.getRank(); }, 
            R"doc(
Returns the protein inference score
:return: Score from protein inference algorithm
)doc")
        .def("getSequence", [](const OpenMS::ProteinHit& self) { return self.getSequence(); }, 
            R"doc(
Returns the rank of this protein hit
:return: Rank (1 = best hit, 2 = second best, etc.)
)doc")
        .def("getAccession", [](const OpenMS::ProteinHit& self) { return self.getAccession(); }, 
            R"doc(
Returns the protein sequence
:return: Full amino acid sequence of the protein (if available)
)doc")
        .def("getDescription", [](const OpenMS::ProteinHit& self) { return self.getDescription(); }, 
            R"doc(
Returns the protein accession
:return: Database accession/identifier (e.g., "P12345" for UniProt)
)doc")
        .def("getCoverage", [](const OpenMS::ProteinHit& self) { return self.getCoverage(); }, 
            R"doc(
Returns the protein description
:return: Human-readable protein name/description from database
)doc")
        .def("setScore", [](OpenMS::ProteinHit& self, double score) { return self.setScore(score); }, "score"_a, 
            R"doc(
Returns the sequence coverage percentage
:return: Percentage of protein sequence covered by identified peptides
Value is in range 0-100 (e.g., 45.2 means 45.2% coverage)
)doc")
        .def("setRank", [](OpenMS::ProteinHit& self, unsigned int newrank) { return self.setRank(newrank); }, "newrank"_a, 
            R"doc(
Sets the protein inference score
:param score: Score to set
)doc")
        .def("setSequence", [](OpenMS::ProteinHit& self, const OpenMS::String& sequence) { return self.setSequence(sequence); }, "sequence"_a, 
            R"doc(
Sets the rank
:param rank: Rank among all protein candidates
)doc")
        .def("setSequence", [](OpenMS::ProteinHit& self, OpenMS::String& sequence) { return self.setSequence(sequence); }, "sequence"_a, 
            R"doc(
Sets the rank
:param rank: Rank among all protein candidates
)doc")
        .def("setAccession", [](OpenMS::ProteinHit& self, const OpenMS::String& accession) { return self.setAccession(accession); }, "accession"_a, 
            R"doc(
Sets the protein sequence
:param sequence: Full amino acid sequence
)doc")
        .def("setDescription", [](OpenMS::ProteinHit& self, const OpenMS::String& description) { return self.setDescription(description); }, "description"_a, 
            R"doc(
Sets the protein accession
:param accession: Database accession/identifier
)doc")
        .def("setCoverage", [](OpenMS::ProteinHit& self, double coverage) { return self.setCoverage(coverage); }, "coverage"_a, 
            R"doc(
Sets the protein description
:param description: Human-readable protein name/description
)doc")
        .def("isDecoy", [](const OpenMS::ProteinHit& self) { return self.isDecoy(); }, 
            R"doc(
Sets the sequence coverage percentage
:param coverage: Percentage (0-100) of sequence covered by peptides
)doc")
        
        .def("__hash__", [](const OpenMS::ProteinHit& self) { return std::hash<OpenMS::ProteinHit>{}(self); })
        .def("__repr__", [](const OpenMS::ProteinHit& self) {
            std::ostringstream os;
            os << "ProteinHit(accession='" << self.getAccession() << "'"
               << ", score=" << self.getScore()
               << ", coverage=" << self.getCoverage() << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::ProteinHit& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;

    // -----------------------------------------------------------------------
    // ProteinIdentification
    // -----------------------------------------------------------------------
    auto proteinidentification_class = nb::class_<OpenMS::ProteinIdentification>(m, "ProteinIdentification",
        R"doc(
Representation of a protein identification run
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProteinIdentification &>())
        .def("__copy__", [](const OpenMS::ProteinIdentification& self) { return OpenMS::ProteinIdentification(self); })
        .def("__deepcopy__", [](const OpenMS::ProteinIdentification& self, nb::dict) { return OpenMS::ProteinIdentification(self); }, "memo"_a)
        .def_static("getAllNamesOfPeakMassType", []() { return OpenMS::ProteinIdentification::getAllNamesOfPeakMassType(); }, "Returns all peak mass type names known to OpenMS")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getHits", [](OpenMS::ProteinIdentification& self) -> std::vector<OpenMS::ProteinHit> & { return self.getHits(); }, nb::rv_policy::reference_internal, "Returns the protein hits")
        .def("insertHit", [](OpenMS::ProteinIdentification& self, const OpenMS::ProteinHit& input) { return self.insertHit(input); }, "input"_a, "Appends a protein hit")
        .def("insertHit", [](OpenMS::ProteinIdentification& self, OpenMS::ProteinHit& input) { return self.insertHit(input); }, "input"_a, "Appends a protein hit")
        .def("setHits", [](OpenMS::ProteinIdentification& self, const std::vector<OpenMS::ProteinHit>& hits) { return self.setHits(hits); }, "hits"_a, "Sets the protein hits")
        .def("getProteinGroups", [](OpenMS::ProteinIdentification& self) -> std::vector<OpenMS::ProteinIdentification::ProteinGroup> & { return self.getProteinGroups(); }, nb::rv_policy::reference_internal, "Returns the protein groups")
        .def("insertProteinGroup", [](OpenMS::ProteinIdentification& self, const OpenMS::ProteinIdentification::ProteinGroup& group) { return self.insertProteinGroup(group); }, "group"_a, "Appends a new protein group")
        .def("getIndistinguishableProteins", [](OpenMS::ProteinIdentification& self) -> std::vector<OpenMS::ProteinIdentification::ProteinGroup> & { return self.getIndistinguishableProteins(); }, nb::rv_policy::reference_internal, "Returns the indistinguishable proteins")
        .def("insertIndistinguishableProteins", [](OpenMS::ProteinIdentification& self, const OpenMS::ProteinIdentification::ProteinGroup& group) { return self.insertIndistinguishableProteins(group); }, "group"_a, "Appends new indistinguishable proteins")
        .def("getSignificanceThreshold", [](const OpenMS::ProteinIdentification& self) { return self.getSignificanceThreshold(); }, "Returns the protein significance threshold value")
        .def("setSignificanceThreshold", [](OpenMS::ProteinIdentification& self, double value) { return self.setSignificanceThreshold(value); }, "value"_a, "Sets the protein significance threshold value")
        .def("getScoreType", [](const OpenMS::ProteinIdentification& self) { return self.getScoreType(); }, "Returns the protein score type")
        .def("setScoreType", [](OpenMS::ProteinIdentification& self, const OpenMS::String& type) { return self.setScoreType(type); }, "type"_a, "Sets the protein score type")
        .def("isHigherScoreBetter", [](const OpenMS::ProteinIdentification& self) { return self.isHigherScoreBetter(); }, "Returns true if a higher score represents a better score")
        .def("setHigherScoreBetter", [](OpenMS::ProteinIdentification& self, bool higher_is_better) { return self.setHigherScoreBetter(higher_is_better); }, "higher_is_better"_a, "Sets the orientation of the score (is higher better?)")
        .def("sort", [](OpenMS::ProteinIdentification& self) { return self.sort(); }, "Sorts the protein hits according to their score")
        .def("computeCoverage", [](OpenMS::ProteinIdentification& self, const OpenMS::PeptideIdentificationList& pep_ids) { self.computeCoverage(pep_ids); }, "pep_ids"_a,
            R"doc(Compute the coverage (in percent) of all ProteinHits given PeptideHits.

Does not return anything but stores the coverage inside the ProteinHit objects.

:param pep_ids: List of PeptideIdentification objects
:raises Exception.MissingInformation: if ProteinHits do not have sequence information
)doc")
        .def("computeCoverage", [](OpenMS::ProteinIdentification& self, const OpenMS::ConsensusMap& cmap, bool use_unassigned_ids) { self.computeCoverage(cmap, use_unassigned_ids); }, "cmap"_a, "use_unassigned_ids"_a,
            R"doc(Compute the coverage (in percent) of all ProteinHits given a ConsensusMap.

Does not return anything but stores the coverage inside the ProteinHit objects.

:param cmap: ConsensusMap containing peptide identifications
:param use_unassigned_ids: Whether to also use unassigned peptide IDs
)doc")
        .def("setPrimaryMSRunPath", [](OpenMS::ProteinIdentification& self, const std::vector<OpenMS::String>& s, bool raw) { self.setPrimaryMSRunPath(s, raw); }, "s"_a, "raw"_a = false,
            R"doc(Set the file paths to the primary MS runs.

:param s: The file paths
:param raw: Store paths to the raw files (or equivalent) rather than mzMLs
)doc")
        .def("setPrimaryMSRunPath", [](OpenMS::ProteinIdentification& self, const std::vector<OpenMS::String>& s, OpenMS::MSExperiment& e) { self.setPrimaryMSRunPath(s, e); }, "s"_a, "e"_a,
            R"doc(Set the file path to the primary MS run but try to use the mzML annotated in the MSExperiment.

:param s: The file paths
:param e: MSExperiment to try to get the mzML path from
)doc")
        .def("getDateTime", [](const OpenMS::ProteinIdentification& self) -> const OpenMS::DateTime & { return self.getDateTime(); }, nb::rv_policy::reference_internal, "Returns the date of the protein identification run")
        .def("setDateTime", [](OpenMS::ProteinIdentification& self, const OpenMS::DateTime& date) { return self.setDateTime(date); }, "date"_a, "Sets the date of the protein identification run")
        .def("setSearchEngine", [](OpenMS::ProteinIdentification& self, const OpenMS::String& search_engine) { return self.setSearchEngine(search_engine); }, "search_engine"_a, "Sets the search engine type")
        .def("getSearchEngine", [](const OpenMS::ProteinIdentification& self) { return self.getSearchEngine(); }, "Returns the type of search engine used")
        .def("setSearchEngineVersion", [](OpenMS::ProteinIdentification& self, const OpenMS::String& search_engine_version) { return self.setSearchEngineVersion(search_engine_version); }, "search_engine_version"_a, "Sets the search engine version")
        .def("getSearchEngineVersion", [](const OpenMS::ProteinIdentification& self) { return self.getSearchEngineVersion(); }, "Returns the search engine version")
        .def("setSearchParameters", [](OpenMS::ProteinIdentification& self, const OpenMS::ProteinIdentification::SearchParameters& search_parameters) { return self.setSearchParameters(search_parameters); }, "search_parameters"_a, "Sets the search parameters")
        .def("setSearchParameters", [](OpenMS::ProteinIdentification& self, OpenMS::ProteinIdentification::SearchParameters& search_parameters) { return self.setSearchParameters(search_parameters); }, "search_parameters"_a, "Sets the search parameters")
        .def("getSearchParameters", [](OpenMS::ProteinIdentification& self) -> OpenMS::ProteinIdentification::SearchParameters & { return self.getSearchParameters(); }, nb::rv_policy::reference_internal, "Returns the search parameters")
        .def("getIdentifier", [](const OpenMS::ProteinIdentification& self) { return self.getIdentifier(); }, "Returns the identifier")
        .def("setIdentifier", [](OpenMS::ProteinIdentification& self, const OpenMS::String& id) { return self.setIdentifier(id); }, "id"_a, "Sets the identifier")
        .def("addPrimaryMSRunPath", [](OpenMS::ProteinIdentification& self, const OpenMS::String& s, bool raw) { return self.addPrimaryMSRunPath(s, raw); }, "s"_a, "raw"_a = false)
        .def("addPrimaryMSRunPath", [](OpenMS::ProteinIdentification& self, const std::vector<OpenMS::String>& s, bool raw) { return self.addPrimaryMSRunPath(s, raw); }, "s"_a, "raw"_a = false)
        .def("getPrimaryMSRunPath", [](const OpenMS::ProteinIdentification& self, bool raw) { std::vector<OpenMS::String> output; self.getPrimaryMSRunPath(output, raw); return output; }, "raw"_a = false)
        
        .def("__hash__", [](const OpenMS::ProteinIdentification& self) { return std::hash<OpenMS::ProteinIdentification>{}(self); })
        ;
    // PeakMassType enum nested under ProteinIdentification
    nb::enum_<OpenMS::ProteinIdentification::PeakMassType>(proteinidentification_class, "PeakMassType")
        .value("MONOISOTOPIC", OpenMS::ProteinIdentification::PeakMassType::MONOISOTOPIC)
        .value("AVERAGE", OpenMS::ProteinIdentification::PeakMassType::AVERAGE)
        .value("SIZE_OF_PEAKMASSTYPE", OpenMS::ProteinIdentification::PeakMassType::SIZE_OF_PEAKMASSTYPE)
        ;
    def_MetaInfoInterface<OpenMS::ProteinIdentification>(proteinidentification_class);

    // -----------------------------------------------------------------------
    // Publication
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedExperimentHelper::Publication, OpenMS::CVTermList>(m, "Publication", "CVTermList")
        .def(nb::init<>())
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        
        .def_rw("id", &OpenMS::TargetedExperimentHelper::Publication::id)
        ;

    // -----------------------------------------------------------------------
    // RangeBase
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::RangeBase>(m, "RangeBase", "OpenMS class RangeBase")
        .def(nb::init<>())
        .def(nb::init<double>())
        .def(nb::init<double, double>())
        .def(nb::init<const OpenMS::RangeBase &>())
        .def("__copy__", [](const OpenMS::RangeBase& self) { return OpenMS::RangeBase(self); })
        .def("__deepcopy__", [](const OpenMS::RangeBase& self, nb::dict) { return OpenMS::RangeBase(self); }, "memo"_a)
        .def("contains", [](const OpenMS::RangeBase& self, double value) { return self.contains(value); }, "value"_a, "Is value within [min, max]?")
        .def("contains", [](const OpenMS::RangeBase& self, const OpenMS::RangeBase& inner_range) { return self.contains(inner_range); }, "inner_range"_a, "Is value within [min, max]?")
        .def("setMin", [](OpenMS::RangeBase& self, double min) { return self.setMin(min); }, "min"_a)
        .def("setMax", [](OpenMS::RangeBase& self, double max) { return self.setMax(max); }, "max"_a)
        .def("getMin", [](const OpenMS::RangeBase& self) { return self.getMin(); })
        .def("getMax", [](const OpenMS::RangeBase& self) { return self.getMax(); })
        .def("extend", [](OpenMS::RangeBase& self, const OpenMS::RangeBase& other) { return self.extend(other); }, "other"_a, "Extend the range such that it includes the given @p value")
        .def("extend", [](OpenMS::RangeBase& self, double value) { return self.extend(value); }, "value"_a, "Extend the range such that it includes the given @p value")
        ;

    // -----------------------------------------------------------------------
    // RangeIntensity
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::RangeIntensity, OpenMS::RangeBase>(m, "RangeIntensity", "OpenMS class RangeIntensity")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::RangeIntensity &>())
        .def("__copy__", [](const OpenMS::RangeIntensity& self) { return OpenMS::RangeIntensity(self); })
        .def("__deepcopy__", [](const OpenMS::RangeIntensity& self, nb::dict) { return OpenMS::RangeIntensity(self); }, "memo"_a)
        .def("setMinIntensity", [](OpenMS::RangeIntensity& self, double min) { return self.setMinIntensity(min); }, "min"_a)
        .def("setMaxIntensity", [](OpenMS::RangeIntensity& self, double max) { return self.setMaxIntensity(max); }, "max"_a)
        .def("getMinIntensity", [](const OpenMS::RangeIntensity& self) { return self.getMinIntensity(); }, "Returns the minimum intensity")
        .def("getMaxIntensity", [](const OpenMS::RangeIntensity& self) { return self.getMaxIntensity(); }, "Returns the maximum intensity")
        .def("extendIntensity", [](OpenMS::RangeIntensity& self, double value) { return self.extendIntensity(value); }, "value"_a, "Extend the range such that it includes the given @p value")
        .def("containsIntensity", [](const OpenMS::RangeIntensity& self, double value) { return self.containsIntensity(value); }, "value"_a, "Is value within [min, max]?")
        .def("containsIntensity", [](const OpenMS::RangeIntensity& self, const OpenMS::RangeBase& inner_range) { return self.containsIntensity(inner_range); }, "inner_range"_a, "Is value within [min, max]?")
        .def("contains", [](const OpenMS::RangeIntensity& self, double value) { return self.contains(value); }, "value"_a, "Is value within [min, max]?")
        .def("setMin", [](OpenMS::RangeIntensity& self, double min) { return self.setMin(min); }, "min"_a)
        .def("setMax", [](OpenMS::RangeIntensity& self, double max) { return self.setMax(max); }, "max"_a)
        .def("getMin", [](const OpenMS::RangeIntensity& self) { return self.getMin(); })
        .def("getMax", [](const OpenMS::RangeIntensity& self) { return self.getMax(); })
        .def("extend", [](OpenMS::RangeIntensity& self, const OpenMS::RangeBase& other) { return self.extend(other); }, "other"_a, "Extend the range such that it includes the given @p value")
        ;

    // -----------------------------------------------------------------------
    // RangeMZ
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::RangeMZ, OpenMS::RangeBase>(m, "RangeMZ", "OpenMS class RangeMZ")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::RangeMZ &>())
        .def("__copy__", [](const OpenMS::RangeMZ& self) { return OpenMS::RangeMZ(self); })
        .def("__deepcopy__", [](const OpenMS::RangeMZ& self, nb::dict) { return OpenMS::RangeMZ(self); }, "memo"_a)
        .def("setMinMZ", [](OpenMS::RangeMZ& self, double min) { return self.setMinMZ(min); }, "min"_a)
        .def("setMaxMZ", [](OpenMS::RangeMZ& self, double max) { return self.setMaxMZ(max); }, "max"_a)
        .def("getMinMZ", [](const OpenMS::RangeMZ& self) { return self.getMinMZ(); })
        .def("getMaxMZ", [](const OpenMS::RangeMZ& self) { return self.getMaxMZ(); })
        .def("extendMZ", [](OpenMS::RangeMZ& self, double value) { return self.extendMZ(value); }, "value"_a, "Extend the range such that it includes the given @p value")
        .def("containsMZ", [](const OpenMS::RangeMZ& self, double value) { return self.containsMZ(value); }, "value"_a, "Is value within [min, max]?")
        .def("containsMZ", [](const OpenMS::RangeMZ& self, const OpenMS::RangeBase& inner_range) { return self.containsMZ(inner_range); }, "inner_range"_a, "Is value within [min, max]?")
        .def("contains", [](const OpenMS::RangeMZ& self, double value) { return self.contains(value); }, "value"_a, "Is value within [min, max]?")
        .def("setMin", [](OpenMS::RangeMZ& self, double min) { return self.setMin(min); }, "min"_a)
        .def("setMax", [](OpenMS::RangeMZ& self, double max) { return self.setMax(max); }, "max"_a)
        .def("getMin", [](const OpenMS::RangeMZ& self) { return self.getMin(); })
        .def("getMax", [](const OpenMS::RangeMZ& self) { return self.getMax(); })
        .def("extend", [](OpenMS::RangeMZ& self, const OpenMS::RangeBase& other) { return self.extend(other); }, "other"_a, "Extend the range such that it includes the given @p value")
        ;

    // -----------------------------------------------------------------------
    // RangeMobility
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::RangeMobility, OpenMS::RangeBase>(m, "RangeMobility", "OpenMS class RangeMobility")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::RangeMobility &>())
        .def("__copy__", [](const OpenMS::RangeMobility& self) { return OpenMS::RangeMobility(self); })
        .def("__deepcopy__", [](const OpenMS::RangeMobility& self, nb::dict) { return OpenMS::RangeMobility(self); }, "memo"_a)
        .def("setMinMobility", [](OpenMS::RangeMobility& self, double min) { return self.setMinMobility(min); }, "min"_a)
        .def("setMaxMobility", [](OpenMS::RangeMobility& self, double max) { return self.setMaxMobility(max); }, "max"_a)
        .def("getMinMobility", [](const OpenMS::RangeMobility& self) { return self.getMinMobility(); }, "Returns the minimum Mobility")
        .def("getMaxMobility", [](const OpenMS::RangeMobility& self) { return self.getMaxMobility(); }, "Returns the maximum Mobility")
        .def("extendMobility", [](OpenMS::RangeMobility& self, double value) { return self.extendMobility(value); }, "value"_a, "Extend the range such that it includes the given @p value")
        .def("containsMobility", [](const OpenMS::RangeMobility& self, double value) { return self.containsMobility(value); }, "value"_a, "Is value within [min, max]?")
        .def("containsMobility", [](const OpenMS::RangeMobility& self, const OpenMS::RangeBase& inner_range) { return self.containsMobility(inner_range); }, "inner_range"_a, "Is value within [min, max]?")
        .def("contains", [](const OpenMS::RangeMobility& self, double value) { return self.contains(value); }, "value"_a, "Is value within [min, max]?")
        .def("setMin", [](OpenMS::RangeMobility& self, double min) { return self.setMin(min); }, "min"_a)
        .def("setMax", [](OpenMS::RangeMobility& self, double max) { return self.setMax(max); }, "max"_a)
        .def("getMin", [](const OpenMS::RangeMobility& self) { return self.getMin(); })
        .def("getMax", [](const OpenMS::RangeMobility& self) { return self.getMax(); })
        .def("extend", [](OpenMS::RangeMobility& self, const OpenMS::RangeBase& other) { return self.extend(other); }, "other"_a, "Extend the range such that it includes the given @p value")
        ;

    // -----------------------------------------------------------------------
    // RangeRT
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::RangeRT, OpenMS::RangeBase>(m, "RangeRT", "OpenMS class RangeRT")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::RangeRT &>())
        .def("__copy__", [](const OpenMS::RangeRT& self) { return OpenMS::RangeRT(self); })
        .def("__deepcopy__", [](const OpenMS::RangeRT& self, nb::dict) { return OpenMS::RangeRT(self); }, "memo"_a)
        .def("setMinRT", [](OpenMS::RangeRT& self, double min) { return self.setMinRT(min); }, "min"_a)
        .def("setMaxRT", [](OpenMS::RangeRT& self, double max) { return self.setMaxRT(max); }, "max"_a)
        .def("getMinRT", [](const OpenMS::RangeRT& self) { return self.getMinRT(); })
        .def("getMaxRT", [](const OpenMS::RangeRT& self) { return self.getMaxRT(); })
        .def("extendRT", [](OpenMS::RangeRT& self, double value) { return self.extendRT(value); }, "value"_a, "Extend the range such that it includes the given @p value")
        .def("containsRT", [](const OpenMS::RangeRT& self, double value) { return self.containsRT(value); }, "value"_a, "Is value within [min, max]?")
        .def("containsRT", [](const OpenMS::RangeRT& self, const OpenMS::RangeBase& inner_range) { return self.containsRT(inner_range); }, "inner_range"_a, "Is value within [min, max]?")
        .def("contains", [](const OpenMS::RangeRT& self, double value) { return self.contains(value); }, "value"_a, "Is value within [min, max]?")
        .def("setMin", [](OpenMS::RangeRT& self, double min) { return self.setMin(min); }, "min"_a)
        .def("setMax", [](OpenMS::RangeRT& self, double max) { return self.setMax(max); }, "max"_a)
        .def("getMin", [](const OpenMS::RangeRT& self) { return self.getMin(); })
        .def("getMax", [](const OpenMS::RangeRT& self) { return self.getMax(); })
        .def("extend", [](OpenMS::RangeRT& self, const OpenMS::RangeBase& other) { return self.extend(other); }, "other"_a, "Extend the range such that it includes the given @p value")
        ;

    // -----------------------------------------------------------------------
    // ReactionMonitoringTransition
    // -----------------------------------------------------------------------
    auto reactionmonitoringtransition_class = nb::class_<OpenMS::ReactionMonitoringTransition, OpenMS::CVTermList>(m, "ReactionMonitoringTransition", 
        R"doc(
CVTermList

This class stores a SRM/MRM transition
This class is capable of representing a <Transition> tag in a TraML
document completely and contains all associated information
The default values for precursor m/z is 0.0 which indicates that it is
uninitialized
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ReactionMonitoringTransition &>())
        .def("__copy__", [](const OpenMS::ReactionMonitoringTransition& self) { return OpenMS::ReactionMonitoringTransition(self); })
        .def("__deepcopy__", [](const OpenMS::ReactionMonitoringTransition& self, nb::dict) { return OpenMS::ReactionMonitoringTransition(self); }, "memo"_a)
        .def("setName", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a)
        .def("getName", [](const OpenMS::ReactionMonitoringTransition& self) { return self.getName(); })
        .def("setNativeID", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::String& name) { return self.setNativeID(name); }, "name"_a)
        .def("getNativeID", [](const OpenMS::ReactionMonitoringTransition& self) { return self.getNativeID(); })
        .def("setPeptideRef", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::String& peptide_ref) { return self.setPeptideRef(peptide_ref); }, "peptide_ref"_a)
        .def("getPeptideRef", [](const OpenMS::ReactionMonitoringTransition& self) { return self.getPeptideRef(); })
        .def("setCompoundRef", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::String& compound_ref) { return self.setCompoundRef(compound_ref); }, "compound_ref"_a)
        .def("getCompoundRef", [](const OpenMS::ReactionMonitoringTransition& self) { return self.getCompoundRef(); })
        .def("setPrecursorMZ", [](OpenMS::ReactionMonitoringTransition& self, double mz) { return self.setPrecursorMZ(mz); }, "mz"_a, "Sets the precursor mz (Q1 value)")
        .def("getPrecursorMZ", [](const OpenMS::ReactionMonitoringTransition& self) { return self.getPrecursorMZ(); }, "Returns the precursor mz (Q1 value)")
        .def("hasPrecursorCVTerms", [](const OpenMS::ReactionMonitoringTransition& self) { return self.hasPrecursorCVTerms(); }, "Returns true if precursor CV Terms exist (means it is safe to call getPrecursorCVTermList)")
        .def("setPrecursorCVTermList", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::CVTermList& list) { return self.setPrecursorCVTermList(list); }, "list"_a, "Sets a list of precursor CV Terms")
        .def("addPrecursorCVTerm", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::CVTerm& cv_term) { return self.addPrecursorCVTerm(cv_term); }, "cv_term"_a, "Adds precursor CV Term")
        .def("getPrecursorCVTermList", [](const OpenMS::ReactionMonitoringTransition& self) -> const OpenMS::CVTermList & { return self.getPrecursorCVTermList(); }, nb::rv_policy::reference_internal, "Obtains the list of CV Terms for the precursor")
        .def("setProductMZ", [](OpenMS::ReactionMonitoringTransition& self, double mz) { return self.setProductMZ(mz); }, "mz"_a)
        .def("getProductMZ", [](const OpenMS::ReactionMonitoringTransition& self) { return self.getProductMZ(); })
        .def("setProduct", [](OpenMS::ReactionMonitoringTransition& self, OpenMS::ReactionMonitoringTransition::Product product) { self.setProduct(product); }, "product"_a, "Sets the product ion")
        .def("getProduct", [](const OpenMS::ReactionMonitoringTransition& self) -> const OpenMS::ReactionMonitoringTransition::Product& { return self.getProduct(); }, nb::rv_policy::reference_internal, "Returns the product ion")
        .def("getProductChargeState", [](const OpenMS::ReactionMonitoringTransition& self) { return self.getProductChargeState(); }, "Returns the charge state of the product")
        .def("isProductChargeStateSet", [](const OpenMS::ReactionMonitoringTransition& self) { return self.isProductChargeStateSet(); }, "Returns true if charge state of product is already set")
        .def("addProductCVTerm", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::CVTerm& cv_term) { return self.addProductCVTerm(cv_term); }, "cv_term"_a)
        .def("getIntermediateProducts", [](const OpenMS::ReactionMonitoringTransition& self) -> const std::vector<OpenMS::TargetedExperimentHelper::TraMLProduct> & { return self.getIntermediateProducts(); }, nb::rv_policy::reference_internal)
        .def("addIntermediateProduct", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::TargetedExperimentHelper::TraMLProduct& product) { return self.addIntermediateProduct(product); }, "product"_a)
        .def("setIntermediateProducts", [](OpenMS::ReactionMonitoringTransition& self, const std::vector<OpenMS::TargetedExperimentHelper::TraMLProduct>& products) { return self.setIntermediateProducts(products); }, "products"_a)
        .def("hasPrediction", [](const OpenMS::ReactionMonitoringTransition& self) { return self.hasPrediction(); }, "Returns true if a Prediction object exists (means it is safe to call getPrediction)")
        .def("setPrediction", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::ReactionMonitoringTransition::Prediction& prediction) { self.setPrediction(prediction); }, "prediction"_a, "Sets the Prediction object")
        .def("getPrediction", [](const OpenMS::ReactionMonitoringTransition& self) -> const OpenMS::ReactionMonitoringTransition::Prediction& { return self.getPrediction(); }, nb::rv_policy::reference_internal, "Returns the Prediction object")
        .def("addPredictionTerm", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::CVTerm& prediction) { return self.addPredictionTerm(prediction); }, "prediction"_a, "Adds prediction term")
        .def("getDecoyTransitionType", [](const OpenMS::ReactionMonitoringTransition& self) { return self.getDecoyTransitionType(); }, "Returns the type of transition (target or decoy)")
        .def("setDecoyTransitionType", [](OpenMS::ReactionMonitoringTransition& self, const OpenMS::ReactionMonitoringTransition::DecoyTransitionType& d) { return self.setDecoyTransitionType(d); }, "d"_a, "Sets the type of transition (target or decoy)")
        .def("getLibraryIntensity", [](const OpenMS::ReactionMonitoringTransition& self) { return self.getLibraryIntensity(); }, "Returns the library intensity (ion count or normalized ion count from a spectral library)")
        .def("setLibraryIntensity", [](OpenMS::ReactionMonitoringTransition& self, double intensity) { return self.setLibraryIntensity(intensity); }, "intensity"_a, "Sets the library intensity (ion count or normalized ion count from a spectral library)")
        .def("setRetentionTime", [](OpenMS::ReactionMonitoringTransition& self, OpenMS::TargetedExperimentHelper::RetentionTime rt) { return self.setRetentionTime(rt); }, "rt"_a)
        .def("getRetentionTime", [](const OpenMS::ReactionMonitoringTransition& self) -> const OpenMS::TargetedExperimentHelper::RetentionTime & { return self.getRetentionTime(); }, nb::rv_policy::reference_internal)
        .def("isDetectingTransition", [](const OpenMS::ReactionMonitoringTransition& self) { return self.isDetectingTransition(); })
        .def("setDetectingTransition", [](OpenMS::ReactionMonitoringTransition& self, bool val) { return self.setDetectingTransition(val); }, "val"_a)
        .def("isIdentifyingTransition", [](const OpenMS::ReactionMonitoringTransition& self) { return self.isIdentifyingTransition(); })
        .def("setIdentifyingTransition", [](OpenMS::ReactionMonitoringTransition& self, bool val) { return self.setIdentifyingTransition(val); }, "val"_a)
        .def("isQuantifyingTransition", [](const OpenMS::ReactionMonitoringTransition& self) { return self.isQuantifyingTransition(); })
        .def("setQuantifyingTransition", [](OpenMS::ReactionMonitoringTransition& self, bool val) { return self.setQuantifyingTransition(val); }, "val"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        
        .def("__hash__", [](const OpenMS::ReactionMonitoringTransition& self) { return std::hash<OpenMS::ReactionMonitoringTransition>{}(self); })
        .def("__repr__", [](const OpenMS::ReactionMonitoringTransition& self) {
            std::ostringstream oss;
            oss << std::fixed;
            oss << "ReactionMonitoringTransition(";
            bool first = true;
            auto nid = self.getNativeID();
            if (!nid.empty()) {
                oss << "id='" << nid << "'";
                first = false;
            }
            if (!first) oss << ", ";
            oss << std::setprecision(4) << "precursor_mz=" << self.getPrecursorMZ();
            oss << ", product_mz=" << self.getProductMZ();
            auto pref = self.getPeptideRef();
            if (!pref.empty()) {
                oss << ", peptide_ref='" << pref << "'";
            }
            if (self.getLibraryIntensity() > 0) {
                oss << std::setprecision(1) << ", library_intensity=" << self.getLibraryIntensity();
            }
            oss << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::ReactionMonitoringTransition& self) {
            // Delegate to __repr__ via Python
            return nb::cast(self).attr("__repr__")();
        })
        ;
    // DecoyTransitionType enum nested under ReactionMonitoringTransition
    nb::enum_<OpenMS::ReactionMonitoringTransition::DecoyTransitionType>(reactionmonitoringtransition_class, "DecoyTransitionType", "Whether a transition is target or decoy")
        .value("UNKNOWN", OpenMS::ReactionMonitoringTransition::DecoyTransitionType::UNKNOWN)
        .value("TARGET", OpenMS::ReactionMonitoringTransition::DecoyTransitionType::TARGET)
        .value("DECOY", OpenMS::ReactionMonitoringTransition::DecoyTransitionType::DECOY)

        .export_values();

    // -----------------------------------------------------------------------
    // RetentionTime
    // -----------------------------------------------------------------------
    auto retentiontime_class = nb::class_<OpenMS::TargetedExperimentHelper::RetentionTime, OpenMS::CVTermListInterface>(m, "RetentionTime",
        R"doc(
Represents a retention time entry for targeted experiment compounds and peptides
CVTermList
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TargetedExperimentHelper::RetentionTime &>())
        .def("__copy__", [](const OpenMS::TargetedExperimentHelper::RetentionTime& self) { return OpenMS::TargetedExperimentHelper::RetentionTime(self); })
        .def("__deepcopy__", [](const OpenMS::TargetedExperimentHelper::RetentionTime& self, nb::dict) { return OpenMS::TargetedExperimentHelper::RetentionTime(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("isRTset", [](const OpenMS::TargetedExperimentHelper::RetentionTime& self) { return self.isRTset(); })
        .def("setRT", [](OpenMS::TargetedExperimentHelper::RetentionTime& self, double rt) { return self.setRT(rt); }, "rt"_a)
        .def("getRT", [](const OpenMS::TargetedExperimentHelper::RetentionTime& self) { return self.getRT(); })
        .def(nb::self != nb::self)
        
        .def_rw("software_ref", &OpenMS::TargetedExperimentHelper::RetentionTime::software_ref)
        .def_rw("retention_time_unit", &OpenMS::TargetedExperimentHelper::RetentionTime::retention_time_unit)
        .def_rw("retention_time_type", &OpenMS::TargetedExperimentHelper::RetentionTime::retention_time_type)
        ;
    // RTUnit enum nested under RetentionTime
    nb::enum_<OpenMS::TargetedExperimentHelper::RetentionTime::RTUnit>(retentiontime_class, "RTUnit")
        .value("SECOND", OpenMS::TargetedExperimentHelper::RetentionTime::RTUnit::SECOND)
        .value("MINUTE", OpenMS::TargetedExperimentHelper::RetentionTime::RTUnit::MINUTE)
        .value("UNKNOWN", OpenMS::TargetedExperimentHelper::RetentionTime::RTUnit::UNKNOWN)
        .value("SIZE_OF_RTUNIT", OpenMS::TargetedExperimentHelper::RetentionTime::RTUnit::SIZE_OF_RTUNIT)
        .export_values();
    // RTType enum nested under RetentionTime
    nb::enum_<OpenMS::TargetedExperimentHelper::RetentionTime::RTType>(retentiontime_class, "RTType")
        .value("LOCAL", OpenMS::TargetedExperimentHelper::RetentionTime::RTType::LOCAL)
        .value("NORMALIZED", OpenMS::TargetedExperimentHelper::RetentionTime::RTType::NORMALIZED)
        .value("PREDICTED", OpenMS::TargetedExperimentHelper::RetentionTime::RTType::PREDICTED)
        .value("HPINS", OpenMS::TargetedExperimentHelper::RetentionTime::RTType::HPINS)
        .value("IRT", OpenMS::TargetedExperimentHelper::RetentionTime::RTType::IRT)
        .value("UNKNOWN", OpenMS::TargetedExperimentHelper::RetentionTime::RTType::UNKNOWN)
        .value("SIZE_OF_RTTYPE", OpenMS::TargetedExperimentHelper::RetentionTime::RTType::SIZE_OF_RTTYPE)
        .export_values();

    // -----------------------------------------------------------------------
    // Sample
    // -----------------------------------------------------------------------
    auto sample_class = nb::class_<OpenMS::Sample, OpenMS::MetaInfoInterface>(m, "Sample", 
        R"doc(
Meta information about the sample
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Sample &>())
        .def("__copy__", [](const OpenMS::Sample& self) { return OpenMS::Sample(self); })
        .def("__deepcopy__", [](const OpenMS::Sample& self, nb::dict) { return OpenMS::Sample(self); }, "memo"_a)
        .def_static("getAllNamesOfSampleState", []() { return OpenMS::Sample::getAllNamesOfSampleState(); }, "Returns all sample state names known to OpenMS")
        .def_static("sampleStateToString", [](OpenMS::Sample::SampleState state) { return OpenMS::Sample::sampleStateToString(state); }, "state"_a, "Convert a SampleState enum to string. Throws Exception::InvalidValue if state is SIZE_OF_SAMPLESTATE")
        .def_static("toSampleState", [](const OpenMS::String& name) { return OpenMS::Sample::toSampleState(name); }, "name"_a, "Convert a string to SampleState enum. Throws Exception::InvalidValue if name is not found")
        .def(nb::self == nb::self)
        .def("getName", [](const OpenMS::Sample& self) { return self.getName(); })
        .def("setName", [](OpenMS::Sample& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a)
        .def("getOrganism", [](const OpenMS::Sample& self) { return self.getOrganism(); })
        .def("setOrganism", [](OpenMS::Sample& self, const OpenMS::String& organism) { return self.setOrganism(organism); }, "organism"_a)
        .def("getNumber", [](const OpenMS::Sample& self) { return self.getNumber(); }, "Returns the sample number")
        .def("setNumber", [](OpenMS::Sample& self, const OpenMS::String& number) { return self.setNumber(number); }, "number"_a, "Sets the sample number (e.g. sample ID)")
        .def("getComment", [](const OpenMS::Sample& self) { return self.getComment(); }, 
            R"doc(
Returns the comment (default "")
)doc")
        .def("setComment", [](OpenMS::Sample& self, const OpenMS::String& comment) { return self.setComment(comment); }, "comment"_a, "Sets the comment (may contain newline characters)")
        .def("getState", [](const OpenMS::Sample& self) { return self.getState(); }, "Returns the state of aggregation (default SAMPLENULL)")
        .def("setState", [](OpenMS::Sample& self, OpenMS::Sample::SampleState state) { return self.setState(state); }, "state"_a, "Sets the state of aggregation")
        .def("getMass", [](const OpenMS::Sample& self) { return self.getMass(); }, "Returns the mass (in gram) (default 0.0)")
        .def("setMass", [](OpenMS::Sample& self, double mass) { return self.setMass(mass); }, "mass"_a, "Sets the mass (in gram)")
        .def("getVolume", [](const OpenMS::Sample& self) { return self.getVolume(); }, "Returns the volume (in ml) (default 0.0)")
        .def("setVolume", [](OpenMS::Sample& self, double volume) { return self.setVolume(volume); }, "volume"_a, "Sets the volume (in ml)")
        .def("getConcentration", [](const OpenMS::Sample& self) { return self.getConcentration(); }, "Returns the concentration (in g/l) (default 0.0)")
        .def("setConcentration", [](OpenMS::Sample& self, double concentration) { return self.setConcentration(concentration); }, "concentration"_a, "Sets the concentration (in g/l)")
        .def("getSubsamples", [](OpenMS::Sample& self) -> std::vector<OpenMS::Sample> & { return self.getSubsamples(); }, nb::rv_policy::reference_internal, "Returns a reference to the vector of subsamples that were combined to create this sample")
        .def("setSubsamples", [](OpenMS::Sample& self, const std::vector<OpenMS::Sample>& subsamples) { return self.setSubsamples(subsamples); }, "subsamples"_a, "Sets the vector of subsamples that were combined to create this sample")
        .def(nb::self != nb::self)
        
        ;
    // SampleState enum nested under Sample
    nb::enum_<OpenMS::Sample::SampleState>(sample_class, "SampleState")
        .value("SAMPLENULL", OpenMS::Sample::SampleState::SAMPLENULL)
        .value("SOLID", OpenMS::Sample::SampleState::SOLID)
        .value("LIQUID", OpenMS::Sample::SampleState::LIQUID)
        .value("GAS", OpenMS::Sample::SampleState::GAS)
        .value("SOLUTION", OpenMS::Sample::SampleState::SOLUTION)
        .value("EMULSION", OpenMS::Sample::SampleState::EMULSION)
        .value("SUSPENSION", OpenMS::Sample::SampleState::SUSPENSION)
        .value("SIZE_OF_SAMPLESTATE", OpenMS::Sample::SampleState::SIZE_OF_SAMPLESTATE)
        ;

    // -----------------------------------------------------------------------
    // ScanWindow
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ScanWindow, OpenMS::MetaInfoInterface>(m, "ScanWindow", "OpenMS class ScanWindow")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ScanWindow &>())
        .def("__copy__", [](const OpenMS::ScanWindow& self) { return OpenMS::ScanWindow(self); })
        .def("__deepcopy__", [](const OpenMS::ScanWindow& self, nb::dict) { return OpenMS::ScanWindow(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        
        .def_rw("begin", &OpenMS::ScanWindow::begin)
        .def_rw("end", &OpenMS::ScanWindow::end)
        ;

    // -----------------------------------------------------------------------
    // SearchParameters
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProteinIdentification::SearchParameters, OpenMS::MetaInfoInterface>(m, "SearchParameters", "OpenMS class SearchParameters")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ProteinIdentification::SearchParameters &>())
        .def("__copy__", [](const OpenMS::ProteinIdentification::SearchParameters& self) { return OpenMS::ProteinIdentification::SearchParameters(self); })
        .def("__deepcopy__", [](const OpenMS::ProteinIdentification::SearchParameters& self, nb::dict) { return OpenMS::ProteinIdentification::SearchParameters(self); }, "memo"_a)
        .def_rw("db", &OpenMS::ProteinIdentification::SearchParameters::db)
        .def_rw("db_version", &OpenMS::ProteinIdentification::SearchParameters::db_version)
        .def_rw("taxonomy", &OpenMS::ProteinIdentification::SearchParameters::taxonomy)
        .def_rw("charges", &OpenMS::ProteinIdentification::SearchParameters::charges)
        .def_rw("mass_type", &OpenMS::ProteinIdentification::SearchParameters::mass_type)
        .def_rw("fixed_modifications", &OpenMS::ProteinIdentification::SearchParameters::fixed_modifications)
        .def_rw("variable_modifications", &OpenMS::ProteinIdentification::SearchParameters::variable_modifications)
        .def_rw("missed_cleavages", &OpenMS::ProteinIdentification::SearchParameters::missed_cleavages)
        .def_rw("fragment_mass_tolerance", &OpenMS::ProteinIdentification::SearchParameters::fragment_mass_tolerance)
        .def_rw("fragment_mass_tolerance_ppm", &OpenMS::ProteinIdentification::SearchParameters::fragment_mass_tolerance_ppm)
        .def_rw("precursor_mass_tolerance", &OpenMS::ProteinIdentification::SearchParameters::precursor_mass_tolerance)
        .def_rw("precursor_mass_tolerance_ppm", &OpenMS::ProteinIdentification::SearchParameters::precursor_mass_tolerance_ppm)
        .def_rw("digestion_enzyme", &OpenMS::ProteinIdentification::SearchParameters::digestion_enzyme)
        ;

    // -----------------------------------------------------------------------
    // Software
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Software, OpenMS::CVTermList>(m, "Software", "Description of the software used for processing")
        .def(nb::init<const OpenMS::Software &>())
        .def(nb::init<OpenMS::String, OpenMS::String>())
        .def("__copy__", [](const OpenMS::Software& self) { return OpenMS::Software(self); })
        .def("__deepcopy__", [](const OpenMS::Software& self, nb::dict) { return OpenMS::Software(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("getName", [](const OpenMS::Software& self) { return self.getName(); }, "Returns the name of the software")
        .def("setName", [](OpenMS::Software& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the software")
        .def("getVersion", [](const OpenMS::Software& self) { return self.getVersion(); }, "Returns the software version")
        .def("setVersion", [](OpenMS::Software& self, const OpenMS::String& version) { return self.setVersion(version); }, "version"_a, "Sets the software version")
        .def("__hash__", [](const OpenMS::Software& self) { return std::hash<OpenMS::Software>{}(self); })

        .def(nb::init<>(), "Default constructor")
        ;

    // -----------------------------------------------------------------------
    // SourceFile
    // -----------------------------------------------------------------------
    auto sourcefile_class = nb::class_<OpenMS::SourceFile, OpenMS::CVTermList>(m, "SourceFile", "OpenMS class SourceFile")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SourceFile &>())
        .def("__copy__", [](const OpenMS::SourceFile& self) { return OpenMS::SourceFile(self); })
        .def("__deepcopy__", [](const OpenMS::SourceFile& self, nb::dict) { return OpenMS::SourceFile(self); }, "memo"_a)
        .def_static("getAllNamesOfChecksumType", []() { return OpenMS::SourceFile::getAllNamesOfChecksumType(); }, "Returns all checksum type names known to OpenMS")
        .def("getNameOfFile", [](const OpenMS::SourceFile& self) { return self.getNameOfFile(); }, "Returns the file name")
        .def("setNameOfFile", [](OpenMS::SourceFile& self, const OpenMS::String& name_of_file) { return self.setNameOfFile(name_of_file); }, "name_of_file"_a, "Sets the file name")
        .def("getPathToFile", [](const OpenMS::SourceFile& self) { return self.getPathToFile(); }, "Returns the file path")
        .def("setPathToFile", [](OpenMS::SourceFile& self, const OpenMS::String& path_path_to_file) { return self.setPathToFile(path_path_to_file); }, "path_path_to_file"_a, "Sets the file path")
        .def("getFileSize", [](const OpenMS::SourceFile& self) { return self.getFileSize(); }, "Returns the file size in MB")
        .def("setFileSize", [](OpenMS::SourceFile& self, float file_size) { return self.setFileSize(file_size); }, "file_size"_a, "Sets the file size in MB")
        .def("getFileType", [](const OpenMS::SourceFile& self) { return self.getFileType(); }, "Returns the file type")
        .def("setFileType", [](OpenMS::SourceFile& self, const OpenMS::String& file_type) { return self.setFileType(file_type); }, "file_type"_a, "Sets the file type")
        .def("getChecksum", [](const OpenMS::SourceFile& self) { return self.getChecksum(); }, "Returns the file's checksum")
        .def("setChecksum", [](OpenMS::SourceFile& self, const OpenMS::String& checksum, OpenMS::SourceFile::ChecksumType type) { return self.setChecksum(checksum, type); }, "checksum"_a, "type"_a, "Sets the file's checksum")
        .def("getChecksumType", [](const OpenMS::SourceFile& self) { return self.getChecksumType(); }, "Returns the checksum type")
        .def("getNativeIDType", [](const OpenMS::SourceFile& self) { return self.getNativeIDType(); }, "Returns the native ID type of the spectra")
        .def("setNativeIDType", [](OpenMS::SourceFile& self, const OpenMS::String& type) { return self.setNativeIDType(type); }, "type"_a, "Sets the native ID type of the spectra")
        .def("getNativeIDTypeAccession", [](const OpenMS::SourceFile& self) { return self.getNativeIDTypeAccession(); }, "Returns the nativeID of the spectra")
        .def("setNativeIDTypeAccession", [](OpenMS::SourceFile& self, const OpenMS::String& accesssion) { return self.setNativeIDTypeAccession(accesssion); }, "accesssion"_a, "Sets the native ID of the spectra")
        ;
    // ChecksumType enum nested under SourceFile
    nb::enum_<OpenMS::SourceFile::ChecksumType>(sourcefile_class, "ChecksumType", "Type of file checksum")
        .value("UNKNOWN_CHECKSUM", OpenMS::SourceFile::ChecksumType::UNKNOWN_CHECKSUM)
        .value("SHA1", OpenMS::SourceFile::ChecksumType::SHA1)
        .value("MD5", OpenMS::SourceFile::ChecksumType::MD5)
        .value("SIZE_OF_CHECKSUMTYPE", OpenMS::SourceFile::ChecksumType::SIZE_OF_CHECKSUMTYPE)
        ;

    // -----------------------------------------------------------------------
    // SpectrumRanges (BaseType of SpectrumRangeManager)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumRangeManager::BaseType>(m, "SpectrumRanges",
        R"doc(
Range container for a single MS level holding m/z, intensity, RT, and mobility ranges.

Returned by :meth:`SpectrumRangeManager.byMSLevel` to provide per-MS-level range queries.
)doc")
        .def("getMinRT", [](const OpenMS::SpectrumRangeManager::BaseType& self) { return self.getMinRT(); }, "Get the minimum RT value")
        .def("getMaxRT", [](const OpenMS::SpectrumRangeManager::BaseType& self) { return self.getMaxRT(); }, "Get the maximum RT value")
        .def("getMinMZ", [](const OpenMS::SpectrumRangeManager::BaseType& self) { return self.getMinMZ(); }, "Get the minimum m/z value")
        .def("getMaxMZ", [](const OpenMS::SpectrumRangeManager::BaseType& self) { return self.getMaxMZ(); }, "Get the maximum m/z value")
        .def("getMinIntensity", [](const OpenMS::SpectrumRangeManager::BaseType& self) { return self.getMinIntensity(); }, "Get the minimum intensity value")
        .def("getMaxIntensity", [](const OpenMS::SpectrumRangeManager::BaseType& self) { return self.getMaxIntensity(); }, "Get the maximum intensity value")
        .def("getMinMobility", [](const OpenMS::SpectrumRangeManager::BaseType& self) { return self.getMinMobility(); }, "Get the minimum mobility value")
        .def("getMaxMobility", [](const OpenMS::SpectrumRangeManager::BaseType& self) { return self.getMaxMobility(); }, "Get the maximum mobility value")
        ;

    // -----------------------------------------------------------------------
    // SpectrumRangeManager
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumRangeManager>(m, "SpectrumRangeManager",
        R"doc(
Advanced range manager for MS spectra with separate ranges for each MS level
This class extends the basic RangeManager to provide separate range tracking for different MS levels
(MS1, MS2, etc.). It manages four types of ranges:
- m/z (mass-to-charge ratio)
- intensity
- retention time (RT)
- ion mobility
A global range is tracked for all MS levels, and additional ranges are maintained for each specific MS level.
This allows for efficient querying of ranges for specific MS levels, which is useful for visualization,
filtering, and processing operations that need to work with specific MS levels.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SpectrumRangeManager &>())
        .def("__copy__", [](const OpenMS::SpectrumRangeManager& self) { return OpenMS::SpectrumRangeManager(self); })
        .def("__deepcopy__", [](const OpenMS::SpectrumRangeManager& self, nb::dict) { return OpenMS::SpectrumRangeManager(self); }, "memo"_a)
        .def("clearRanges", [](OpenMS::SpectrumRangeManager& self) { return self.clearRanges(); })
        .def("getMSLevels", [](const OpenMS::SpectrumRangeManager& self) { return self.getMSLevels(); })
        .def("extendRT", [](OpenMS::SpectrumRangeManager& self, double rt, unsigned int ms_level) { return self.extendRT(rt, ms_level); }, "rt"_a, "ms_level"_a = 0)
        .def("extendMZ", [](OpenMS::SpectrumRangeManager& self, double mz, unsigned int ms_level) { return self.extendMZ(mz, ms_level); }, "mz"_a, "ms_level"_a = 0)
        .def("extendUnsafe", [](OpenMS::SpectrumRangeManager& self, const OpenMS::MSSpectrum& spectrum, unsigned int ms_level) { return self.extendUnsafe(spectrum, ms_level); }, "spectrum"_a, "ms_level"_a = 0)
        .def("getMinRT", [](const OpenMS::SpectrumRangeManager& self) { return self.getMinRT(); }, "Get the minimum RT value")
        .def("getMaxRT", [](const OpenMS::SpectrumRangeManager& self) { return self.getMaxRT(); }, "Get the maximum RT value")
        .def("getMinMZ", [](const OpenMS::SpectrumRangeManager& self) { return self.getMinMZ(); }, "Get the minimum m/z value")
        .def("getMaxMZ", [](const OpenMS::SpectrumRangeManager& self) { return self.getMaxMZ(); }, "Get the maximum m/z value")
        .def("getMinIntensity", [](const OpenMS::SpectrumRangeManager& self) { return self.getMinIntensity(); }, "Get the minimum intensity value")
        .def("getMaxIntensity", [](const OpenMS::SpectrumRangeManager& self) { return self.getMaxIntensity(); }, "Get the maximum intensity value")
        .def("getMinMobility", [](const OpenMS::SpectrumRangeManager& self) { return self.getMinMobility(); }, "Get the minimum mobility value")
        .def("getMaxMobility", [](const OpenMS::SpectrumRangeManager& self) { return self.getMaxMobility(); }, "Get the maximum mobility value")
        .def("byMSLevel", [](const OpenMS::SpectrumRangeManager& self, unsigned int ms_level) -> OpenMS::SpectrumRangeManager::BaseType { return self.byMSLevel(ms_level); },
            "ms_level"_a,
            R"doc(Get ranges for a specific MS level.

Returns a SpectrumRanges object containing the m/z, RT, intensity, and mobility
ranges for the specified MS level.

:param ms_level: The MS level (e.g., 1 for MS1, 2 for MS2)
:raises RuntimeError: If no ranges exist for the specified MS level
)doc")
        ;

    // -----------------------------------------------------------------------
    // SpectrumSettings
    // -----------------------------------------------------------------------
    auto spectrumsettings_class = nb::class_<OpenMS::SpectrumSettings, OpenMS::MetaInfoInterface>(m, "SpectrumSettings", 
        R"doc(
Representation of 1D spectrum settings
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SpectrumSettings &>())
        .def("__copy__", [](const OpenMS::SpectrumSettings& self) { return OpenMS::SpectrumSettings(self); })
        .def("__deepcopy__", [](const OpenMS::SpectrumSettings& self, nb::dict) { return OpenMS::SpectrumSettings(self); }, "memo"_a)
        .def_static("getAllNamesOfSpectrumType", []() { return OpenMS::SpectrumSettings::getAllNamesOfSpectrumType(); }, "Returns all spectrum type names known to OpenMS")
        .def_static("spectrumTypeToString", [](OpenMS::SpectrumSettings::SpectrumType type) { return OpenMS::SpectrumSettings::spectrumTypeToString(type); }, "type"_a, "Convert a SpectrumType enum to String. Throws Exception::InvalidValue if type is SIZE_OF_SPECTRUMTYPE")
        .def_static("toSpectrumType", [](const OpenMS::String& name) { return OpenMS::SpectrumSettings::toSpectrumType(name); }, "name"_a, "Convert a string to SpectrumType enum. Throws Exception::InvalidValue if name is not a valid spectrum type")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("unify", [](OpenMS::SpectrumSettings& self, const OpenMS::SpectrumSettings& rhs) { return self.unify(rhs); }, "rhs"_a)
        .def("getType", [](const OpenMS::SpectrumSettings& self) { return self.getType(); }, "Returns the spectrum type (centroided (PEAKS) or profile data (RAW))")
        .def("setType", [](OpenMS::SpectrumSettings& self, OpenMS::SpectrumSettings::SpectrumType type) { return self.setType(type); }, "type"_a, "Sets the spectrum type")
        .def("getNativeID", [](const OpenMS::SpectrumSettings& self) { return self.getNativeID(); }, "Returns the native identifier for the spectrum, used by the acquisition software")
        .def("setNativeID", [](OpenMS::SpectrumSettings& self, const OpenMS::String& native_id) { return self.setNativeID(native_id); }, "native_id"_a, "Sets the native identifier for the spectrum, used by the acquisition software")
        .def("getComment", [](const OpenMS::SpectrumSettings& self) { return self.getComment(); }, "Returns the free-text comment")
        .def("setComment", [](OpenMS::SpectrumSettings& self, const OpenMS::String& comment) { return self.setComment(comment); }, "comment"_a, "Sets the free-text comment")
        .def("getInstrumentSettings", [](OpenMS::SpectrumSettings& self) -> OpenMS::InstrumentSettings & { return self.getInstrumentSettings(); }, nb::rv_policy::reference_internal, "Returns a const reference to the instrument settings of the current spectrum")
        .def("setInstrumentSettings", [](OpenMS::SpectrumSettings& self, const OpenMS::InstrumentSettings& instrument_settings) { return self.setInstrumentSettings(instrument_settings); }, "instrument_settings"_a, "Sets the instrument settings of the current spectrum")
        .def("getAcquisitionInfo", [](OpenMS::SpectrumSettings& self) -> OpenMS::AcquisitionInfo & { return self.getAcquisitionInfo(); }, nb::rv_policy::reference_internal, "Returns a const reference to the acquisition info")
        .def("setAcquisitionInfo", [](OpenMS::SpectrumSettings& self, const OpenMS::AcquisitionInfo& acquisition_info) { return self.setAcquisitionInfo(acquisition_info); }, "acquisition_info"_a, "Sets the acquisition info")
        .def("getSourceFile", [](OpenMS::SpectrumSettings& self) -> OpenMS::SourceFile & { return self.getSourceFile(); }, nb::rv_policy::reference_internal, "Returns a const reference to the source file")
        .def("setSourceFile", [](OpenMS::SpectrumSettings& self, const OpenMS::SourceFile& source_file) { return self.setSourceFile(source_file); }, "source_file"_a, "Sets the source file")
        .def("getPrecursors", [](OpenMS::SpectrumSettings& self) -> std::vector<OpenMS::Precursor> & { return self.getPrecursors(); }, nb::rv_policy::reference_internal, "Returns a const reference to the precursors")
        .def("setPrecursors", [](OpenMS::SpectrumSettings& self, const std::vector<OpenMS::Precursor>& precursors) { return self.setPrecursors(precursors); }, "precursors"_a, "Sets the precursors")
        .def("getProducts", [](OpenMS::SpectrumSettings& self) -> std::vector<OpenMS::Product> & { return self.getProducts(); }, nb::rv_policy::reference_internal, "Returns a const reference to the products")
        .def("setProducts", [](OpenMS::SpectrumSettings& self, const std::vector<OpenMS::Product>& products) { return self.setProducts(products); }, "products"_a, "Sets the products")
        .def("setDataProcessing", [](OpenMS::SpectrumSettings& self, const std::vector<std::shared_ptr<OpenMS::DataProcessing>>& data_processing) { return self.setDataProcessing(data_processing); }, "data_processing"_a)
        .def("getDataProcessing", [](OpenMS::SpectrumSettings& self) -> std::vector<std::shared_ptr<OpenMS::DataProcessing>> & { return self.getDataProcessing(); }, nb::rv_policy::reference_internal)
        
        .def("__hash__", [](const OpenMS::SpectrumSettings& self) { return std::hash<OpenMS::SpectrumSettings>{}(self); })
        ;
    // SpectrumType enum nested under SpectrumSettings
    nb::enum_<OpenMS::SpectrumSettings::SpectrumType>(spectrumsettings_class, "SpectrumType")
        .value("UNKNOWN", OpenMS::SpectrumSettings::SpectrumType::UNKNOWN)
        .value("CENTROID", OpenMS::SpectrumSettings::SpectrumType::CENTROID)
        .value("PROFILE", OpenMS::SpectrumSettings::SpectrumType::PROFILE)
        .value("SIZE_OF_SPECTRUMTYPE", OpenMS::SpectrumSettings::SpectrumType::SIZE_OF_SPECTRUMTYPE)
        ;

    // -----------------------------------------------------------------------
    // StringDataArray
    // -----------------------------------------------------------------------
    auto stringdataarray_class = nb::class_<OpenMS::DataArrays::StringDataArray>(m, "StringDataArray", 
        R"doc(
MetaInfoDescription

The representation of extra string data attached to a spectrum or chromatogram.
Commonly used for storing ion annotation names or other per-peak string annotations.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DataArrays::StringDataArray &>())
        .def("__copy__", [](const OpenMS::DataArrays::StringDataArray& self) { return OpenMS::DataArrays::StringDataArray(self); })
        .def("__deepcopy__", [](const OpenMS::DataArrays::StringDataArray& self, nb::dict) { return OpenMS::DataArrays::StringDataArray(self); }, "memo"_a)
        .def(nb::self != nb::self)
        .def(nb::self == nb::self)
        .def("getName", [](const OpenMS::DataArrays::StringDataArray& self) { return self.getName(); }, "Returns the name of the peak annotations")
        .def("setName", [](OpenMS::DataArrays::StringDataArray& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the peak annotations")
        .def("setDataProcessing", [](OpenMS::DataArrays::StringDataArray& self, const std::vector<std::shared_ptr<OpenMS::DataProcessing>>& data_processing) { return self.setDataProcessing(data_processing); }, "data_processing"_a, "Sets the description of the applied processing")
        
        .def("__len__", [](OpenMS::DataArrays::StringDataArray& self) { return self.size(); })
        .def("size", [](const OpenMS::DataArrays::StringDataArray& self) { return self.size(); }, "Returns the number of elements")
        .def("__getitem__", [](OpenMS::DataArrays::StringDataArray& self, size_t i) -> OpenMS::String& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)
        .def("__setitem__", [](OpenMS::DataArrays::StringDataArray& self, size_t i, const OpenMS::String& val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        })

        .def(nb::init<>(), "Default constructor")
        .def(nb::init<const OpenMS::DataArrays::StringDataArray&>(), "Copy constructor")
        .def("__copy__", [](const OpenMS::DataArrays::StringDataArray& self) { return OpenMS::DataArrays::StringDataArray(self); })
        .def("__deepcopy__", [](const OpenMS::DataArrays::StringDataArray& self, nb::dict) { return OpenMS::DataArrays::StringDataArray(self); }, "memo"_a)

        .def("get_data", [](const OpenMS::DataArrays::StringDataArray& self) {
            std::vector<OpenMS::String> arr(self.begin(), self.end());
            return arr;
        }, "Returns a copy of the data as a list")

        .def("set_data", [](OpenMS::DataArrays::StringDataArray& self, std::vector<OpenMS::String> arr) {
            self.assign(arr.begin(), arr.end());
        }, "data"_a, "Set data from a list")

        .def("push_back", [](OpenMS::DataArrays::StringDataArray& self, const OpenMS::String& val) {
            self.push_back(val);
        }, "val"_a, "Add a value to the array")

        .def("resize", [](OpenMS::DataArrays::StringDataArray& self, size_t n) {
            self.resize(n);
        }, "n"_a, "Resize the array")

        .def("clear", [](OpenMS::DataArrays::StringDataArray& self) {
            self.clear();
        }, "Clear the array")
        .def("getDataProcessing", [](const OpenMS::DataArrays::StringDataArray& self) { return self.getDataProcessing(); }, "Returns the data processing steps")
        .def("__repr__", [](const OpenMS::DataArrays::StringDataArray& self) {
            return "StringDataArray(name='" + self.getName() + "', size=" + std::to_string(self.size()) + ")";
        })
        .def("__str__", [](const OpenMS::DataArrays::StringDataArray& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;
    def_MetaInfoInterface<OpenMS::DataArrays::StringDataArray>(stringdataarray_class);

    // -----------------------------------------------------------------------
    // TargetedExperiment_Interpretation
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedExperimentHelper::Interpretation, OpenMS::CVTermListInterface>(m, "TargetedExperiment_Interpretation",
        R"doc(
Product ion interpretation

Stores information about the ion type, ordinal, and rank for MS product ions.
CVTermListInterface
)doc")
        .def(nb::init<>())
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def_rw("ordinal", &OpenMS::TargetedExperimentHelper::Interpretation::ordinal)
        .def_rw("rank", &OpenMS::TargetedExperimentHelper::Interpretation::rank)
        .def_rw("iontype", &OpenMS::TargetedExperimentHelper::Interpretation::iontype)
        ;

    // -----------------------------------------------------------------------
    // TraMLProduct
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedExperimentHelper::TraMLProduct, OpenMS::CVTermListInterface>(m, "TraMLProduct", "OpenMS class TraMLProduct")
        .def(nb::init<>())
        .def(nb::self == nb::self)
        .def("setChargeState", [](OpenMS::TargetedExperimentHelper::TraMLProduct& self, int charge) { return self.setChargeState(charge); }, "charge"_a)
        .def("hasCharge", [](const OpenMS::TargetedExperimentHelper::TraMLProduct& self) { return self.hasCharge(); })
        .def("getChargeState", [](const OpenMS::TargetedExperimentHelper::TraMLProduct& self) { return self.getChargeState(); })
        .def("getMZ", [](const OpenMS::TargetedExperimentHelper::TraMLProduct& self) { return self.getMZ(); })
        .def("setMZ", [](OpenMS::TargetedExperimentHelper::TraMLProduct& self, double mz) { return self.setMZ(mz); }, "mz"_a)
        .def("getConfigurationList", [](const OpenMS::TargetedExperimentHelper::TraMLProduct& self) -> const std::vector<OpenMS::TargetedExperimentHelper::Configuration> & { return self.getConfigurationList(); }, nb::rv_policy::reference_internal)
        .def("addConfiguration", [](OpenMS::TargetedExperimentHelper::TraMLProduct& self, const OpenMS::TargetedExperimentHelper::Configuration& configuration) { return self.addConfiguration(configuration); }, "configuration"_a)
        .def("getInterpretationList", [](const OpenMS::TargetedExperimentHelper::TraMLProduct& self) -> const std::vector<OpenMS::TargetedExperimentHelper::Interpretation> & { return self.getInterpretationList(); }, nb::rv_policy::reference_internal)
        .def("addInterpretation", [](OpenMS::TargetedExperimentHelper::TraMLProduct& self, const OpenMS::TargetedExperimentHelper::Interpretation& interpretation) { return self.addInterpretation(interpretation); }, "interpretation"_a)
        .def("resetInterpretations", [](OpenMS::TargetedExperimentHelper::TraMLProduct& self) { return self.resetInterpretations(); })
        .def(nb::self != nb::self)
        
        ;

    // -----------------------------------------------------------------------
    // UniqueIdInterface
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::UniqueIdInterface>(m, "UniqueIdInterface", 
        R"doc(
A base class defining a common interface for all classes having a
unique id
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::UniqueIdInterface &>())
        .def("__copy__", [](const OpenMS::UniqueIdInterface& self) { return OpenMS::UniqueIdInterface(self); })
        .def("__deepcopy__", [](const OpenMS::UniqueIdInterface& self, nb::dict) { return OpenMS::UniqueIdInterface(self); }, "memo"_a)
        .def_static("isValid", [](size_t unique_id) { return OpenMS::UniqueIdInterface::isValid(unique_id); }, "unique_id"_a, "Returns true if the unique_id is valid, false otherwise")
        .def("getUniqueId", [](const OpenMS::UniqueIdInterface& self) { return self.getUniqueId(); }, "Returns the unique id")
        .def("clearUniqueId", [](OpenMS::UniqueIdInterface& self) { return self.clearUniqueId(); }, "Clear the unique id. The new unique id will be invalid.  Returns 1 if the unique id was changed, 0 otherwise")
        .def("hasValidUniqueId", [](const OpenMS::UniqueIdInterface& self) { return self.hasValidUniqueId(); }, "Returns whether the unique id is valid.  Returns 1 if the unique id is valid, 0 otherwise")
        .def("hasInvalidUniqueId", [](const OpenMS::UniqueIdInterface& self) { return self.hasInvalidUniqueId(); }, "Returns whether the unique id is invalid.  Returns 1 if the unique id is invalid, 0 otherwise")
        .def("ensureUniqueId", [](OpenMS::UniqueIdInterface& self) { return self.ensureUniqueId(); }, "Assigns a valid unique id, but only if the present one is invalid.  Returns 1 if the unique id was changed, 0 otherwise")
        .def("setUniqueId", [](OpenMS::UniqueIdInterface& self, size_t rhs) { return self.setUniqueId(rhs); }, "rhs"_a, "Assigns a new, valid unique id.  Always returns 1")
        .def("setUniqueId", [](OpenMS::UniqueIdInterface& self, const OpenMS::String& rhs) { return self.setUniqueId(rhs); }, "rhs"_a, "Assigns a new, valid unique id.  Always returns 1")
        ;

    // -----------------------------------------------------------------------
    // ConsensusMap
    // -----------------------------------------------------------------------
    auto consensusmap_class = nb::class_<OpenMS::ConsensusMap>(m, "ConsensusMap", 
        R"doc(
UniqueIdInterface
DocumentIdentifier
RangeManagerRtMzInt
MetaInfoInterface

A container for consensus elements.
A ConsensusMap is a container holding 2-dimensional consensus elements
(ConsensusFeature) which in turn represent analytes that have been
quantified across multiple LC-MS/MS experiments. Each analyte in a
ConsensusFeature is linked to its original LC-MS/MS run, the links are
maintained by the ConsensusMap class.
The map is implemented as a vector of elements of type ConsensusFeature.
To be consistent, all maps who are referenced by ConsensusFeature objects
(through a unique id) need to be registered in this class.
This class supports direct iteration in Python.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ConsensusMap &>())
        .def("__copy__", [](const OpenMS::ConsensusMap& self) { return OpenMS::ConsensusMap(self); })
        .def("__deepcopy__", [](const OpenMS::ConsensusMap& self, nb::dict) { return OpenMS::ConsensusMap(self); }, "memo"_a)
        .def(nb::init<size_t>())
        .def("appendRows", [](OpenMS::ConsensusMap& self, const OpenMS::ConsensusMap& rhs) -> OpenMS::ConsensusMap & { return self.appendRows(rhs); }, "rhs"_a, nb::rv_policy::reference_internal, "Add consensus map entries as new rows")
        .def("appendColumns", [](OpenMS::ConsensusMap& self, const OpenMS::ConsensusMap& rhs) -> OpenMS::ConsensusMap & { return self.appendColumns(rhs); }, "rhs"_a, nb::rv_policy::reference_internal, "Add consensus map entries as new columns")
        .def("clear", [](OpenMS::ConsensusMap& self, bool clear_meta_data) { return self.clear(clear_meta_data); }, "clear_meta_data"_a = true, "Clears all data and meta data")
        .def("getExperimentType", [](const OpenMS::ConsensusMap& self) { return self.getExperimentType(); }, "Non-mutable access to the experiment type")
        .def("setExperimentType", [](OpenMS::ConsensusMap& self, const OpenMS::String& experiment_type) { return self.setExperimentType(experiment_type); }, "experiment_type"_a, "Mutable access to the experiment type")
        .def("sortByIntensity", [](OpenMS::ConsensusMap& self, bool reverse) { return self.sortByIntensity(reverse); }, "reverse"_a = false, "Sorts the peaks according to ascending intensity.")
        .def("sortByRT", [](OpenMS::ConsensusMap& self) { return self.sortByRT(); }, "Sorts the peaks according to RT position")
        .def("sortByMZ", [](OpenMS::ConsensusMap& self) { return self.sortByMZ(); }, "Sorts the peaks according to m/z position")
        .def("sortByPosition", [](OpenMS::ConsensusMap& self) { return self.sortByPosition(); }, "Lexicographically sorts the peaks by their position (First RT then m/z)")
        .def("sortByQuality", [](OpenMS::ConsensusMap& self, bool reverse) { return self.sortByQuality(reverse); }, "reverse"_a = false, "Sorts the peaks according to ascending quality.")
        .def("sortBySize", [](OpenMS::ConsensusMap& self) { return self.sortBySize(); }, "Sorts with respect to the size (number of elements)")
        .def("sortByMaps", [](OpenMS::ConsensusMap& self) { return self.sortByMaps(); }, "Sorts with respect to the sets of maps covered by the consensus features (lexicographically)")
        .def("sortPeptideIdentificationsByMapIndex", [](OpenMS::ConsensusMap& self) { return self.sortPeptideIdentificationsByMapIndex(); }, "Sorts PeptideIdentifications of consensus features with respect to their map index.")
        .def("updateRanges", [](OpenMS::ConsensusMap& self) { return self.updateRanges(); }, "Updates the RT, m/z, and intensity ranges based on contained consensus features")
        .def("clearRanges", [](OpenMS::ConsensusMap& self) { return self.clearRanges(); }, "Clear all ranges")
        .def("empty", [](const OpenMS::ConsensusMap& self) { return self.empty(); }, "Returns True if the map is empty")
        .def("reserve", [](OpenMS::ConsensusMap& self, size_t n) { return self.reserve(n); }, "n"_a, "Reserves space for n consensus features")
        .def("getMinRT", [](const OpenMS::ConsensusMap& self) { return self.getMinRT(); }, "Get the minimum RT value")
        .def("getMaxRT", [](const OpenMS::ConsensusMap& self) { return self.getMaxRT(); }, "Get the maximum RT value")
        .def("getMinMZ", [](const OpenMS::ConsensusMap& self) { return self.getMinMZ(); }, "Get the minimum m/z value")
        .def("getMaxMZ", [](const OpenMS::ConsensusMap& self) { return self.getMaxMZ(); }, "Get the maximum m/z value")
        .def("getMinIntensity", [](const OpenMS::ConsensusMap& self) { return self.getMinIntensity(); }, "Get the minimum intensity value")
        .def("getMaxIntensity", [](const OpenMS::ConsensusMap& self) { return self.getMaxIntensity(); }, "Get the maximum intensity value")
        .def("getProteinIdentifications", [](OpenMS::ConsensusMap& self) -> std::vector<OpenMS::ProteinIdentification> & { return self.getProteinIdentifications(); }, nb::rv_policy::reference_internal, "Returns the protein identification runs stored in this map")
        .def("setProteinIdentifications", [](OpenMS::ConsensusMap& self, const std::vector<OpenMS::ProteinIdentification>& protein_identifications) { return self.setProteinIdentifications(protein_identifications); }, "protein_identifications"_a, "Sets the protein identifications")
        .def("setProteinIdentifications", [](OpenMS::ConsensusMap& self, std::vector<OpenMS::ProteinIdentification>& protein_identifications) { return self.setProteinIdentifications(protein_identifications); }, "protein_identifications"_a, "Sets the protein identifications")
        .def("getUnassignedPeptideIdentifications", [](OpenMS::ConsensusMap& self) -> OpenMS::PeptideIdentificationList & { return self.getUnassignedPeptideIdentifications(); }, nb::rv_policy::reference_internal, "Returns peptide identifications that are not assigned to any consensus feature")
        .def("setUnassignedPeptideIdentifications", [](OpenMS::ConsensusMap& self, const OpenMS::PeptideIdentificationList& unassigned_peptide_identifications) { return self.setUnassignedPeptideIdentifications(unassigned_peptide_identifications); }, "unassigned_peptide_identifications"_a, "Sets the unassigned PeptideIdentificationList")
        .def("getDataProcessing", [](OpenMS::ConsensusMap& self) -> std::vector<OpenMS::DataProcessing> & { return self.getDataProcessing(); }, nb::rv_policy::reference_internal, "Returns a const reference to the description of the applied data processing")
        .def("setDataProcessing", [](OpenMS::ConsensusMap& self, const std::vector<OpenMS::DataProcessing>& processing_method) { return self.setDataProcessing(processing_method); }, "processing_method"_a, "Sets the description of the applied data processing")
        .def("setPrimaryMSRunPath", [](OpenMS::ConsensusMap& self, const std::vector<OpenMS::String>& s) { return self.setPrimaryMSRunPath(s); }, "s"_a, "Sets the file paths to the primary MS run (stored in ColumnHeaders)")
        .def("setPrimaryMSRunPath", [](OpenMS::ConsensusMap& self, const std::vector<OpenMS::String>& s, OpenMS::MSExperiment& e) { return self.setPrimaryMSRunPath(s, e); }, "s"_a, "e"_a)
        .def("getPrimaryMSRunPath", [](const OpenMS::ConsensusMap& self) { std::vector<OpenMS::String> toFill; self.getPrimaryMSRunPath(toFill); return toFill; }, "Returns the MS run path (stored in ColumnHeaders)")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        
        .def("__iter__", [](OpenMS::ConsensusMap& self) { return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::ConsensusMap>(), "ConsensusMap_iter", self.begin(), self.end()); })
        .def("__len__", [](OpenMS::ConsensusMap& self) { return self.size(); })
        .def("__getitem__", [](OpenMS::ConsensusMap& self, size_t i) -> OpenMS::ConsensusFeature & {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)
        .def("__setitem__", [](OpenMS::ConsensusMap& self, size_t i, const OpenMS::ConsensusFeature& val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        }, "i"_a, "val"_a, "Sets consensus feature at index i")

        .def("getColumnHeaders", [](const OpenMS::ConsensusMap& self) {
            return self.getColumnHeaders();
        }, "Returns the column headers")

        .def("setColumnHeaders", [](OpenMS::ConsensusMap& self, const std::map<OpenMS::UInt64, OpenMS::ConsensusMap::ColumnHeader>& headers) {
            self.setColumnHeaders(headers);
        }, "headers"_a, "Sets the column headers")

        .def("size", [](const OpenMS::ConsensusMap& self) {
            return self.size();
        }, "Returns the number of consensus features")

        .def("push_back", [](OpenMS::ConsensusMap& self, const OpenMS::ConsensusFeature& f) {
            self.push_back(f);
        }, "feature"_a, "Add a consensus feature to the map")
        .def("append", [](OpenMS::ConsensusMap& self, const OpenMS::ConsensusFeature& f) {
            self.push_back(f);
        }, "feature"_a, "Add a consensus feature to the map (alias for push_back)")
        .def("extend", [](OpenMS::ConsensusMap& self, const nb::list& items) {
            for (auto item : items) {
                self.push_back(nb::cast<OpenMS::ConsensusFeature>(item));
            }
        }, "items"_a, "Add multiple consensus features from a list")
        .def("extend", [](OpenMS::ConsensusMap& self, const OpenMS::ConsensusMap& other) {
            for (size_t i = 0; i < other.size(); ++i) {
                self.push_back(other[i]);
            }
        }, "other"_a, "Add features from another ConsensusMap")
        .def("setUniqueIds", [](OpenMS::ConsensusMap& self) {
            self.setUniqueId();
            self.applyMemberFunction(&OpenMS::UniqueIdInterface::setUniqueId);
        }, "Sets unique IDs on the map and all its child consensus features")
        .def("__repr__", [](const OpenMS::ConsensusMap& self) {
            return "ConsensusMap(num_consensus_features=" + std::to_string(self.size()) + ")";
        })
        .def("__str__", [](const OpenMS::ConsensusMap& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;
    def_MetaInfoInterface<OpenMS::ConsensusMap>(consensusmap_class);
    def_DocumentIdentifier<OpenMS::ConsensusMap>(consensusmap_class);
    def_UniqueIdInterface<OpenMS::ConsensusMap>(consensusmap_class);

    // -----------------------------------------------------------------------
    // FeatureHandle
    // -----------------------------------------------------------------------
    auto featurehandle_class = nb::class_<OpenMS::FeatureHandle>(m, "FeatureHandle", 
        R"doc(
Representation of a Peak2D, RichPeak2D or Feature
Peak2D
UniqueIdInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<size_t, OpenMS::Peak2D, size_t>())
        .def(nb::init<size_t, OpenMS::BaseFeature>())
        .def(nb::init<const OpenMS::FeatureHandle &>())
        .def("__copy__", [](const OpenMS::FeatureHandle& self) { return OpenMS::FeatureHandle(self); })
        .def("__deepcopy__", [](const OpenMS::FeatureHandle& self, nb::dict) { return OpenMS::FeatureHandle(self); }, "memo"_a)
        .def("getMapIndex", [](const OpenMS::FeatureHandle& self) { return self.getMapIndex(); }, "Returns the map index")
        .def("setMapIndex", [](OpenMS::FeatureHandle& self, size_t i) { return self.setMapIndex(i); }, "i"_a, "Sets the map index")
        .def("setCharge", [](OpenMS::FeatureHandle& self, int charge) { return self.setCharge(charge); }, "charge"_a, "Sets the charge")
        .def("getCharge", [](const OpenMS::FeatureHandle& self) { return self.getCharge(); }, "Returns the charge")
        .def("setWidth", [](OpenMS::FeatureHandle& self, float width) { return self.setWidth(width); }, "width"_a, "Sets the width (FWHM)")
        .def("getWidth", [](const OpenMS::FeatureHandle& self) { return self.getWidth(); }, "Returns the width (FWHM)")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getIntensity", [](const OpenMS::FeatureHandle& self) { return self.getIntensity(); }, "Returns the data point intensity (height)")
        .def("setIntensity", [](OpenMS::FeatureHandle& self, float intensity) { return self.setIntensity(intensity); }, "intensity"_a, "Sets the data point intensity (height)")
        .def("getMZ", [](const OpenMS::FeatureHandle& self) { return self.getMZ(); }, "Returns the m/z coordinate (index 1)")
        .def("setMZ", [](OpenMS::FeatureHandle& self, double coordinate) { return self.setMZ(coordinate); }, "coordinate"_a, "Sets the m/z coordinate (index 1)")
        .def("getRT", [](const OpenMS::FeatureHandle& self) { return self.getRT(); }, "Returns the RT coordinate (index 0)")
        .def("setRT", [](OpenMS::FeatureHandle& self, double coordinate) { return self.setRT(coordinate); }, "coordinate"_a, "Sets the RT coordinate (index 0)")
        .def("getRT", [](OpenMS::FeatureHandle& self) { return self.getRT(); }, "Returns the retention time")
        .def("setRT", [](OpenMS::FeatureHandle& self, double rt) { self.setRT(rt); }, "rt"_a, "Sets the retention time")
        .def("getMZ", [](OpenMS::FeatureHandle& self) { return self.getMZ(); }, "Returns the m/z")
        .def("setMZ", [](OpenMS::FeatureHandle& self, double mz) { self.setMZ(mz); }, "mz"_a, "Sets the m/z")
        .def("getIntensity", [](OpenMS::FeatureHandle& self) { return self.getIntensity(); }, "Returns the intensity")
        .def("setIntensity", [](OpenMS::FeatureHandle& self, float intensity) { self.setIntensity(intensity); }, "intensity"_a, "Sets the intensity")
        .def("__hash__", [](const OpenMS::FeatureHandle& self) { return std::hash<OpenMS::FeatureHandle>{}(self); })
        ;
    def_UniqueIdInterface<OpenMS::FeatureHandle>(featurehandle_class);

    // -----------------------------------------------------------------------
    // FeatureMap
    // -----------------------------------------------------------------------
    auto featuremap_class = nb::class_<OpenMS::FeatureMap>(m, "FeatureMap", 
        R"doc(
UniqueIdInterface
DocumentIdentifier
RangeManagerRtMzInt
MetaInfoInterface

A container for LC-MS features with metadata and identification information
FeatureMap is one of the core data structures in OpenMS for storing detected features
from LC-MS experiments. A feature represents a detected chemical entity (peptide, protein,
metabolite, etc.) with its elution profile and mass information.
Key capabilities:
- Store and manage Feature objects (detected analytes)
- Associate protein and peptide identifications with features
- Sort features by various criteria (RT, m/z, intensity, quality)
- Store experimental metadata and data processing information
- Support direct iteration and indexing in Python
Example usage:
.. code-block:: python
feature_map = oms.FeatureMap()
# Add a feature
feature = oms.Feature()
feature.setRT(1234.5)
feature.setMZ(445.678)
feature.setIntensity(100000.0)
feature_map.push_back(feature)
# Access features
print(f"Number of features: {feature_map.size()}")
first_feature = feature_map[0]
# Sort by RT
feature_map.sortByRT()
# Iterate over features
for feat in feature_map:
print(f"RT: {feat.getRT()}, m/z: {feat.getMZ()}")
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FeatureMap &>())
        .def("__copy__", [](const OpenMS::FeatureMap& self) { return OpenMS::FeatureMap(self); })
        .def("__deepcopy__", [](const OpenMS::FeatureMap& self, nb::dict) { return OpenMS::FeatureMap(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def(nb::self + nb::self)
        .def("sortByIntensity", [](OpenMS::FeatureMap& self, bool reverse) { return self.sortByIntensity(reverse); }, "reverse"_a = false, 
            R"doc(
Sorts features by ascending intensity
After sorting, features can be accessed in order from lowest to highest intensity
)doc")
        .def("sortByPosition", [](OpenMS::FeatureMap& self) { return self.sortByPosition(); }, 
            R"doc(
Sorts features by intensity with optional reverse order
:param reverse: If True, sorts in descending order (highest to lowest intensity)
)doc")
        .def("sortByRT", [](OpenMS::FeatureMap& self) { return self.sortByRT(); }, 
            R"doc(
Sorts features by position using lexicographical comparison
Compares RT first, then m/z for features with the same RT
)doc")
        .def("sortByMZ", [](OpenMS::FeatureMap& self) { return self.sortByMZ(); }, 
            R"doc(
Sorts features by retention time (RT) in ascending order
This is useful for time-based analysis or visualization
)doc")
        .def("sortByOverallQuality", [](OpenMS::FeatureMap& self, bool reverse) { return self.sortByOverallQuality(reverse); }, "reverse"_a = false, 
            R"doc(
Sorts features by mass-to-charge ratio (m/z) in ascending order
Useful for mass-based grouping or analysis
)doc")
        .def("updateRanges", [](OpenMS::FeatureMap& self) { return self.updateRanges(); }, "Updates the RT, m/z, and intensity ranges based on contained features")
        .def("clearRanges", [](OpenMS::FeatureMap& self) { return self.clearRanges(); }, "Clear all ranges")
        .def("getMinRT", [](const OpenMS::FeatureMap& self) { return self.getMinRT(); }, "Get the minimum RT value")
        .def("getMaxRT", [](const OpenMS::FeatureMap& self) { return self.getMaxRT(); }, "Get the maximum RT value")
        .def("getMinMZ", [](const OpenMS::FeatureMap& self) { return self.getMinMZ(); }, "Get the minimum m/z value")
        .def("getMaxMZ", [](const OpenMS::FeatureMap& self) { return self.getMaxMZ(); }, "Get the maximum m/z value")
        .def("getMinIntensity", [](const OpenMS::FeatureMap& self) { return self.getMinIntensity(); }, "Get the minimum intensity value")
        .def("getMaxIntensity", [](const OpenMS::FeatureMap& self) { return self.getMaxIntensity(); }, "Get the maximum intensity value")
        .def("swapFeaturesOnly", [](OpenMS::FeatureMap& self, OpenMS::FeatureMap& from) { return self.swapFeaturesOnly(from); }, "from"_a, "Swaps the feature content (plus its range information) of this map")
        .def("swap", [](OpenMS::FeatureMap& self, OpenMS::FeatureMap& from) { return self.swap(from); }, "from"_a, 
            R"doc(
Sorts features by overall quality score in ascending order
Higher quality scores indicate better feature detection confidence
)doc")
        .def("getProteinIdentifications", [](OpenMS::FeatureMap& self) -> std::vector<OpenMS::ProteinIdentification> & { return self.getProteinIdentifications(); }, nb::rv_policy::reference_internal)
        .def("setProteinIdentifications", [](OpenMS::FeatureMap& self, const std::vector<OpenMS::ProteinIdentification>& protein_identifications) { return self.setProteinIdentifications(protein_identifications); }, "protein_identifications"_a, 
            R"doc(
Returns the protein identification runs stored in this map
:return: Protein identification data from database searches
Protein identifications contain metadata about search parameters and protein hits
)doc")
        .def("getUnassignedPeptideIdentifications", [](OpenMS::FeatureMap& self) -> OpenMS::PeptideIdentificationList & { return self.getUnassignedPeptideIdentifications(); }, nb::rv_policy::reference_internal, 
            R"doc(
Sets the protein identifications for this map
:param protein_ids: Protein identification results to associate with this map
)doc")
        .def("setUnassignedPeptideIdentifications", [](OpenMS::FeatureMap& self, const OpenMS::PeptideIdentificationList& unassigned_peptide_identifications) { return self.setUnassignedPeptideIdentifications(unassigned_peptide_identifications); }, "unassigned_peptide_identifications"_a, 
            R"doc(
Returns peptide identifications that are not assigned to any feature
:return: Unassigned peptide identification results
These are peptide IDs that could not be matched to features, possibly due to feature detection issues or filtering
)doc")
        .def("getDataProcessing", [](OpenMS::FeatureMap& self) -> std::vector<OpenMS::DataProcessing> & { return self.getDataProcessing(); }, nb::rv_policy::reference_internal)
        .def("setDataProcessing", [](OpenMS::FeatureMap& self, const std::vector<OpenMS::DataProcessing>& processing_method) { return self.setDataProcessing(processing_method); }, "processing_method"_a, "Sets the description of the applied data processing")
        .def("setPrimaryMSRunPath", [](OpenMS::FeatureMap& self, const std::vector<OpenMS::String>& s) { return self.setPrimaryMSRunPath(s); }, "s"_a, "Sets the file path to the primary MS run (usually the mzML file obtained after data conversion from raw files)")
        .def("setPrimaryMSRunPath", [](OpenMS::FeatureMap& self, const std::vector<OpenMS::String>& s, OpenMS::MSExperiment& e) { return self.setPrimaryMSRunPath(s, e); }, "s"_a, "e"_a, "Sets the file path to the primary MS run using the mzML annotated in the MSExperiment argument `e`")
        .def("getPrimaryMSRunPath", [](const OpenMS::FeatureMap& self) { std::vector<OpenMS::String> toFill; self.getPrimaryMSRunPath(toFill); return toFill; }, "Returns the file path to the first MS run")
        .def("clear", [](OpenMS::FeatureMap& self, bool clear_meta_data) { return self.clear(clear_meta_data); }, "clear_meta_data"_a = true, 
            R"doc(
Clears all feature data and metadata
After calling this, the map will be empty (size() returns 0)
)doc")
        
        .def("__iter__", [](OpenMS::FeatureMap& self) { return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::FeatureMap>(), "FeatureMap_iter", self.begin(), self.end()); })
        .def("__len__", [](OpenMS::FeatureMap& self) { return self.size(); })
        .def("__getitem__", [](OpenMS::FeatureMap& self, size_t i) -> OpenMS::Feature & {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, nb::rv_policy::reference_internal)
        .def("__setitem__", [](OpenMS::FeatureMap& self, size_t i, const OpenMS::Feature& val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        }, "i"_a, "val"_a, "Sets feature at index i")
        .def("__iadd__", [](OpenMS::FeatureMap& self, const OpenMS::FeatureMap& other) -> OpenMS::FeatureMap& { self += other; return self; }, "other"_a, "Appends all features from another FeatureMap")

        .def("size", [](const OpenMS::FeatureMap& self) {
            return self.size();
        }, "Returns the number of features")

        .def("push_back", [](OpenMS::FeatureMap& self, const OpenMS::Feature& f) {
            self.push_back(f);
        }, "feature"_a, "Add a feature to the map")
        .def("append", [](OpenMS::FeatureMap& self, const OpenMS::Feature& f) {
            self.push_back(f);
        }, "feature"_a, "Add a feature to the map (alias for push_back)")
        .def("extend", [](OpenMS::FeatureMap& self, const nb::list& items) {
            for (auto item : items) {
                self.push_back(nb::cast<OpenMS::Feature>(item));
            }
        }, "items"_a, "Add multiple features from a list")
        .def("extend", [](OpenMS::FeatureMap& self, const OpenMS::FeatureMap& other) {
            for (const auto& f : other) {
                self.push_back(f);
            }
        }, "other"_a, "Add all features from another FeatureMap")
        .def("setUniqueIds", [](OpenMS::FeatureMap& self) {
            self.setUniqueId();
            self.applyMemberFunction(&OpenMS::UniqueIdInterface::setUniqueId);
        }, "Sets unique IDs on the map and all its child features")
        .def("__repr__", [](const OpenMS::FeatureMap& self) {
            return "FeatureMap(num_features=" + std::to_string(self.size()) + ")";
        })
        .def("__str__", [](const OpenMS::FeatureMap& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;
    def_MetaInfoInterface<OpenMS::FeatureMap>(featuremap_class);
    def_DocumentIdentifier<OpenMS::FeatureMap>(featuremap_class);
    def_UniqueIdInterface<OpenMS::FeatureMap>(featuremap_class);

    // -----------------------------------------------------------------------
    // RichPeak2D
    // -----------------------------------------------------------------------
    auto richpeak2d_class = nb::class_<OpenMS::RichPeak2D>(m, "RichPeak2D", 
        R"doc(
A 2-dimensional raw data point or peak with meta information
Peak2D
UniqueIdInterface
MetaInfoInterface
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::RichPeak2D &>())
        .def("__copy__", [](const OpenMS::RichPeak2D& self) { return OpenMS::RichPeak2D(self); })
        .def("__deepcopy__", [](const OpenMS::RichPeak2D& self, nb::dict) { return OpenMS::RichPeak2D(self); }, "memo"_a)
        .def(nb::init<OpenMS::Peak2D>())
        .def(nb::init<OpenMS::DPosition<2>, float>())
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getIntensity", [](const OpenMS::RichPeak2D& self) { return self.getIntensity(); }, "Returns the data point intensity (height)")
        .def("setIntensity", [](OpenMS::RichPeak2D& self, float intensity) { return self.setIntensity(intensity); }, "intensity"_a, "Sets the data point intensity (height)")
        .def("getMZ", [](const OpenMS::RichPeak2D& self) { return self.getMZ(); }, "Returns the m/z coordinate (index 1)")
        .def("setMZ", [](OpenMS::RichPeak2D& self, double coordinate) { return self.setMZ(coordinate); }, "coordinate"_a, "Sets the m/z coordinate (index 1)")
        .def("getRT", [](const OpenMS::RichPeak2D& self) { return self.getRT(); }, "Returns the RT coordinate (index 0)")
        .def("setRT", [](OpenMS::RichPeak2D& self, double coordinate) { return self.setRT(coordinate); }, "coordinate"_a, "Sets the RT coordinate (index 0)")
        
        .def("getRT", [](OpenMS::RichPeak2D& self) { return self.getRT(); }, "Returns the retention time")
        .def("setRT", [](OpenMS::RichPeak2D& self, double rt) { self.setRT(rt); }, "rt"_a, "Sets the retention time")
        .def("getMZ", [](OpenMS::RichPeak2D& self) { return self.getMZ(); }, "Returns the m/z")
        .def("setMZ", [](OpenMS::RichPeak2D& self, double mz) { self.setMZ(mz); }, "mz"_a, "Sets the m/z")
        .def("getIntensity", [](OpenMS::RichPeak2D& self) { return self.getIntensity(); }, "Returns the intensity")
        .def("setIntensity", [](OpenMS::RichPeak2D& self, float intensity) { self.setIntensity(intensity); }, "intensity"_a, "Sets the intensity")
        ;
    def_MetaInfoInterface<OpenMS::RichPeak2D>(richpeak2d_class);
    def_UniqueIdInterface<OpenMS::RichPeak2D>(richpeak2d_class);

    // -----------------------------------------------------------------------
    // BaseFeature
    // -----------------------------------------------------------------------
    auto basefeature_class = nb::class_<OpenMS::BaseFeature, OpenMS::RichPeak2D>(m, "BaseFeature", 
        R"doc(
A basic LC-MS feature
UniqueIdInterface
RichPeak2D
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::BaseFeature &>())
        .def("__copy__", [](const OpenMS::BaseFeature& self) { return OpenMS::BaseFeature(self); })
        .def("__deepcopy__", [](const OpenMS::BaseFeature& self, nb::dict) { return OpenMS::BaseFeature(self); }, "memo"_a)
        .def(nb::init<const OpenMS::BaseFeature &, size_t>())
        .def(nb::init<OpenMS::Peak2D>())
        .def(nb::init<OpenMS::RichPeak2D>())
        .def(nb::init<OpenMS::FeatureHandle>())
        .def("getQuality", [](const OpenMS::BaseFeature& self) { return self.getQuality(); }, "Returns the overall quality")
        .def("setQuality", [](OpenMS::BaseFeature& self, float q) { return self.setQuality(q); }, "q"_a, "Sets the overall quality")
        .def("getWidth", [](const OpenMS::BaseFeature& self) { return self.getWidth(); }, "Returns the features width (full width at half max, FWHM)")
        .def("setWidth", [](OpenMS::BaseFeature& self, float fwhm) { return self.setWidth(fwhm); }, "fwhm"_a, "Sets the width of the feature (FWHM)")
        .def("getCharge", [](const OpenMS::BaseFeature& self) { return self.getCharge(); }, "Returns the charge state")
        .def("setCharge", [](OpenMS::BaseFeature& self, const int& ch) { return self.setCharge(ch); }, "ch"_a, "Sets the charge state")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getPeptideIdentifications", [](OpenMS::BaseFeature& self) -> OpenMS::PeptideIdentificationList & { return self.getPeptideIdentifications(); }, nb::rv_policy::reference_internal, "Returns the PeptideIdentification vector")
        .def("setPeptideIdentifications", [](OpenMS::BaseFeature& self, const OpenMS::PeptideIdentificationList& peptides) { return self.setPeptideIdentifications(peptides); }, "peptides"_a, "Sets the PeptideIdentification vector")
        .def("getAnnotationState", [](const OpenMS::BaseFeature& self) { return self.getAnnotationState(); }, "State of peptide identifications attached to this feature. If one ID has multiple hits, the output depends on the top-hit only")
        
        ;
    // AnnotationState enum nested under BaseFeature
    nb::enum_<OpenMS::BaseFeature::AnnotationState>(basefeature_class, "AnnotationState", "State of peptide identification annotation for a feature")
        .value("FEATURE_ID_NONE", OpenMS::BaseFeature::AnnotationState::FEATURE_ID_NONE)
        .value("FEATURE_ID_SINGLE", OpenMS::BaseFeature::AnnotationState::FEATURE_ID_SINGLE)
        .value("FEATURE_ID_MULTIPLE_SAME", OpenMS::BaseFeature::AnnotationState::FEATURE_ID_MULTIPLE_SAME)
        .value("FEATURE_ID_MULTIPLE_DIVERGENT", OpenMS::BaseFeature::AnnotationState::FEATURE_ID_MULTIPLE_DIVERGENT)
        .value("SIZE_OF_ANNOTATIONSTATE", OpenMS::BaseFeature::AnnotationState::SIZE_OF_ANNOTATIONSTATE)

        .export_values();

    // -----------------------------------------------------------------------
    // ConsensusFeature
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusFeature, OpenMS::BaseFeature>(m, "ConsensusFeature", 
        R"doc(
UniqueIdInterface
BaseFeature

A consensus feature spanning multiple LC-MS/MS experiments.
A ConsensusFeature represents analytes that have been
quantified across multiple LC-MS/MS experiments. Each analyte in a
ConsensusFeature is linked to its original LC-MS/MS run through a
unique identifier.
Get access to the underlying features through getFeatureList()
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ConsensusFeature &>())
        .def("__copy__", [](const OpenMS::ConsensusFeature& self) { return OpenMS::ConsensusFeature(self); })
        .def("__deepcopy__", [](const OpenMS::ConsensusFeature& self, nb::dict) { return OpenMS::ConsensusFeature(self); }, "memo"_a)
        .def(nb::init<OpenMS::BaseFeature>())
        .def(nb::init<size_t, OpenMS::Peak2D, size_t>())
        .def(nb::init<size_t, OpenMS::BaseFeature>())
        .def("insert", [](OpenMS::ConsensusFeature& self, const OpenMS::ConsensusFeature& cf) { return self.insert(cf); }, "cf"_a, "Inserts a feature handle from a peak with map index and element index")
        .def("insert", [](OpenMS::ConsensusFeature& self, OpenMS::ConsensusFeature& cf) { return self.insert(cf); }, "cf"_a, "Inserts a feature handle from a peak with map index and element index")
        .def("insert", [](OpenMS::ConsensusFeature& self, const OpenMS::FeatureHandle& handle) { return self.insert(handle); }, "handle"_a, "Inserts a feature handle from a peak with map index and element index")
        .def("insert", [](OpenMS::ConsensusFeature& self, OpenMS::FeatureHandle& handle) { return self.insert(handle); }, "handle"_a, "Inserts a feature handle from a peak with map index and element index")
        .def("insert", [](OpenMS::ConsensusFeature& self, const std::set<OpenMS::FeatureHandle, OpenMS::FeatureHandle::IndexLess>& handle_set) { return self.insert(handle_set); }, "handle_set"_a, "Inserts a feature handle from a peak with map index and element index")
        .def("insert", [](OpenMS::ConsensusFeature& self, std::set<OpenMS::FeatureHandle, OpenMS::FeatureHandle::IndexLess>& handle_set) { return self.insert(handle_set); }, "handle_set"_a, "Inserts a feature handle from a peak with map index and element index")
        .def("insert", [](OpenMS::ConsensusFeature& self, size_t map_index, const OpenMS::Peak2D& element, size_t element_index) { return self.insert(map_index, element, element_index); }, "map_index"_a, "element"_a, "element_index"_a, "Inserts a feature handle from a peak with map index and element index")
        .def("insert", [](OpenMS::ConsensusFeature& self, size_t map_index, const OpenMS::BaseFeature& element) { return self.insert(map_index, element); }, "map_index"_a, "element"_a, "Inserts a feature handle from a base feature with map index")
        .def("getFeatureList", [](const OpenMS::ConsensusFeature& self) { return self.getFeatureList(); }, "Returns a list of all contained feature handles")
        .def("computeConsensus", [](OpenMS::ConsensusFeature& self) { return self.computeConsensus(); }, "Computes and updates the consensus position, intensity, and charge")
        .def("computeMonoisotopicConsensus", [](OpenMS::ConsensusFeature& self) { return self.computeMonoisotopicConsensus(); }, "Computes and updates the consensus position, intensity, and charge")
        .def("addRatio", [](OpenMS::ConsensusFeature& self, const OpenMS::ConsensusFeature::Ratio& r) { return self.addRatio(r); }, "r"_a, "Connects a ratio to the ConsensusFeature.")
        .def("setRatios", [](OpenMS::ConsensusFeature& self, std::vector<OpenMS::ConsensusFeature::Ratio>& rs) { return self.setRatios(rs); }, "rs"_a, "Connects the ratios to the ConsensusFeature.")
        .def("getRatios", [](OpenMS::ConsensusFeature& self) -> std::vector<OpenMS::ConsensusFeature::Ratio> & { return self.getRatios(); }, nb::rv_policy::reference_internal, "Get the ratio vector.")
        .def("size", [](const OpenMS::ConsensusFeature& self) { return self.size(); }, "Returns the number of feature handles in this consensus feature")
        .def("clear", [](OpenMS::ConsensusFeature& self) { return self.clear(); }, "Clears all feature handles from this consensus feature")
        .def("empty", [](const OpenMS::ConsensusFeature& self) { return self.empty(); }, "Returns True if this consensus feature contains no feature handles")
        .def("getQuality", [](const OpenMS::ConsensusFeature& self) { return self.getQuality(); }, "Returns the overall quality")
        .def("setQuality", [](OpenMS::ConsensusFeature& self, float q) { return self.setQuality(q); }, "q"_a, "Sets the overall quality")
        .def("getWidth", [](const OpenMS::ConsensusFeature& self) { return self.getWidth(); }, "Returns the features width (full width at half max, FWHM)")
        .def("setWidth", [](OpenMS::ConsensusFeature& self, float fwhm) { return self.setWidth(fwhm); }, "fwhm"_a, "Sets the width of the feature (FWHM)")
        .def("getCharge", [](const OpenMS::ConsensusFeature& self) { return self.getCharge(); }, "Returns the charge state")
        .def("setCharge", [](OpenMS::ConsensusFeature& self, const int& ch) { return self.setCharge(ch); }, "ch"_a, "Sets the charge state")
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("getPeptideIdentifications", [](const OpenMS::ConsensusFeature& self) -> const OpenMS::PeptideIdentificationList & { return self.getPeptideIdentifications(); }, nb::rv_policy::reference_internal, "Returns the PeptideIdentification vector")
        .def("setPeptideIdentifications", [](OpenMS::ConsensusFeature& self, const OpenMS::PeptideIdentificationList& peptides) { return self.setPeptideIdentifications(peptides); }, "peptides"_a, "Sets the PeptideIdentification vector")
        .def("getAnnotationState", [](const OpenMS::ConsensusFeature& self) { return self.getAnnotationState(); }, "State of peptide identifications attached to this feature. If one ID has multiple hits, the output depends on the top-hit only")
        .def("computeDechargeConsensus", [](OpenMS::ConsensusFeature& self, const OpenMS::FeatureMap& fm, bool intensity_weighted_averaging) { return self.computeDechargeConsensus(fm, intensity_weighted_averaging); }, "fm"_a, "intensity_weighted_averaging"_a = false, "Computes and updates the consensus position, intensity, and charge using decharge grouping")

        .def("__len__", [](OpenMS::ConsensusFeature& self) { return self.size(); })
        .def("__repr__", [](const OpenMS::ConsensusFeature& self) {
            std::ostringstream os;
            os << "ConsensusFeature(rt=" << self.getRT() << ", mz=" << self.getMZ()
               << ", intensity=" << self.getIntensity() << ", charge=" << self.getCharge()
               << ", num_features=" << self.size() << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::ConsensusFeature& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;

    // -----------------------------------------------------------------------
    // Feature
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Feature, OpenMS::BaseFeature>(m, "Feature", 
        R"doc(
UniqueIdInterface
RichPeak2D

An LC-MS feature representing a detected analyte signal
The Feature class represents a two-dimensional (RT and m/z) signal from an analyte
in LC-MS data. It is one of the core data structures in OpenMS for representing
detected peaks or compounds.
A Feature stores:
- Position: retention time (RT) and mass-to-charge ratio (m/z)
- Intensity: the signal strength (typically total ion count)
- Quality metrics: scores indicating detection confidence
- Charge state: the charge of the ion
- Convex hulls: the 2D area occupied by the feature in RT-m/z space
- Peptide identifications: for identified peptides (optional)
- Subordinate features: for isotopic peaks or related signals
By convention, the feature's position is at the maximum of the elution profile
(RT dimension) and at the monoisotopic peak (m/z dimension).
Example usage:
.. code-block:: python
feature = oms.Feature()
feature.setRT(1234.5)  # Set retention time in seconds
feature.setMZ(445.678)  # Set m/z value
feature.setIntensity(100000.0)  # Set intensity
feature.setCharge(2)  # Set charge state
feature.setOverallQuality(0.95)  # Set quality score (0-1)
# Access the values
print(f"RT: {feature.getRT()}, m/z: {feature.getMZ()}, charge: {feature.getCharge()}")
)doc")
        .def(nb::init<>())
        .def(nb::init<OpenMS::BaseFeature>())
        .def(nb::init<const OpenMS::Feature &>())
        .def("__copy__", [](const OpenMS::Feature& self) { return OpenMS::Feature(self); })
        .def("__deepcopy__", [](const OpenMS::Feature& self, nb::dict) { return OpenMS::Feature(self); }, "memo"_a)
        .def("getOverallQuality", [](const OpenMS::Feature& self) { return self.getOverallQuality(); },
            R"doc(
Returns the overall quality score of the feature
:return: Overall quality score (typically 0-1, where 1 is highest quality)
This score represents the overall confidence in the feature detection
)doc")
        .def("setOverallQuality", [](OpenMS::Feature& self, float q) { return self.setOverallQuality(q); }, "q"_a,
            R"doc(
Sets the overall quality score of the feature
:param q: Overall quality score (typically 0-1, where 1 is highest quality)
)doc")
        .def("getQuality", [](const OpenMS::Feature& self, size_t index) { return self.getQuality(index); }, "index"_a,
            R"doc(
Returns the quality score in a specific dimension
:param index: The dimension index (0 for RT, 1 for m/z)
:return: Quality score for the specified dimension (typically 0-1 range)
)doc")
        .def("setQuality", [](OpenMS::Feature& self, size_t index, float q) { return self.setQuality(index, q); }, "index"_a, "q"_a,
            R"doc(
Sets the quality score for a specific dimension
:param index: The dimension index (0 for RT, 1 for m/z)
:param q: Quality score to set (typically 0-1 range)
)doc")
        .def("getConvexHulls", [](OpenMS::Feature& self) -> std::vector<OpenMS::ConvexHull2D> & { return self.getConvexHulls(); }, nb::rv_policy::reference_internal,
            R"doc(
Returns the convex hulls of individual mass traces
:return: List of convex hulls, one for each isotopic mass trace
Each isotopic peak typically has its own convex hull in RT-m/z space
)doc")
        .def("setConvexHulls", [](OpenMS::Feature& self, const std::vector<OpenMS::ConvexHull2D>& hulls) { return self.setConvexHulls(hulls); }, "hulls"_a,
            R"doc(
Sets the convex hulls of individual mass traces
:param hulls: List of convex hulls to associate with this feature
)doc")
        .def("getConvexHull", [](const OpenMS::Feature& self) -> OpenMS::ConvexHull2D & { return self.getConvexHull(); }, nb::rv_policy::reference_internal,
            R"doc(
Returns the overall convex hull of the feature
:return: The overall 2D convex hull encompassing all mass traces
This is the union of all individual mass trace convex hulls
)doc")
        .def("encloses", [](const OpenMS::Feature& self, double rt, double mz) { return self.encloses(rt, mz); }, "rt"_a, "mz"_a,
            R"doc(
Checks if the feature's convex hulls enclose a given position
:param rt: Retention time in seconds
:param mz: Mass-to-charge ratio
:return: True if the position (rt, mz) is within the feature's convex hulls, False otherwise
This uses the feature's convex hull representation to determine spatial containment
)doc")
        .def(nb::self == nb::self)
        .def("getSubordinates", [](OpenMS::Feature& self) -> std::vector<OpenMS::Feature> & { return self.getSubordinates(); }, nb::rv_policy::reference_internal,
            R"doc(
Returns subordinate features (e.g., isotopic peaks)
:return: List of subordinate features associated with this feature
Subordinate features often represent individual isotopic peaks of the same compound
)doc")
        .def("setSubordinates", [](OpenMS::Feature& self, const std::vector<OpenMS::Feature>& rhs) { return self.setSubordinates(rhs); }, "rhs"_a,
            R"doc(
Sets the subordinate features
:param rhs: List of subordinate features to associate with this feature
)doc")
        .def("getWidth", [](const OpenMS::Feature& self) { return self.getWidth(); },
            R"doc(
Returns the width (FWHM) of the feature in RT dimension
:return: Full Width at Half Maximum (FWHM) in seconds
Represents the elution peak width
)doc")
        .def("setWidth", [](OpenMS::Feature& self, float fwhm) { return self.setWidth(fwhm); }, "fwhm"_a,
            R"doc(
Sets the width (FWHM) of the feature in RT dimension
:param fwhm: Full Width at Half Maximum in seconds
)doc")
        .def("getCharge", [](const OpenMS::Feature& self) { return self.getCharge(); },
            R"doc(
Returns the charge state of the feature
:return: Charge state (e.g., 2 for doubly charged ions, 0 if unknown)
)doc")
        .def("setCharge", [](OpenMS::Feature& self, const int& ch) { return self.setCharge(ch); }, "ch"_a,
            R"doc(
Sets the charge state of the feature
:param ch: Charge state (e.g., 2 for doubly charged ions)
)doc")
        .def(nb::self != nb::self)
        .def("getPeptideIdentifications", [](const OpenMS::Feature& self) -> const OpenMS::PeptideIdentificationList & { return self.getPeptideIdentifications(); }, nb::rv_policy::reference_internal,
            R"doc(
Returns the peptide identifications associated with this feature
:return: List of peptide identifications from database search
Only relevant for peptide features. Contains results from peptide identification tools
)doc")
        .def("setPeptideIdentifications", [](OpenMS::Feature& self, const OpenMS::PeptideIdentificationList& peptides) { return self.setPeptideIdentifications(peptides); }, "peptides"_a,
            R"doc(
Sets the peptide identifications associated with this feature
:param peptides: List of peptide identifications to associate with this feature
)doc")
        .def("getAnnotationState", [](const OpenMS::Feature& self) { return self.getAnnotationState(); },
            R"doc(
Returns the annotation state of the feature
:return: Enum indicating the annotation status of this feature
)doc")
        
        .def("__copy__", [](const OpenMS::Feature& self) { return OpenMS::Feature(self); })
        .def("__deepcopy__", [](const OpenMS::Feature& self, nb::dict) { return OpenMS::Feature(self); }, "memo"_a)
        .def("__repr__", [](const OpenMS::Feature& self) {
            std::ostringstream os;
            os << "Feature(rt=" << self.getRT() << ", mz=" << self.getMZ()
               << ", intensity=" << self.getIntensity() << ", charge=" << self.getCharge()
               << ", quality=" << self.getOverallQuality() << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::Feature& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;


    // SpectrumHelper is a namespace-level helper class with only static methods
    struct SpectrumHelper_Dummy {};
    // -----------------------------------------------------------------------
    // SpectrumHelper
    // -----------------------------------------------------------------------
    nb::class_<SpectrumHelper_Dummy>(m, "SpectrumHelper", "OpenMS class SpectrumHelper")
        .def_static("removePeaks", [](OpenMS::MSSpectrum& spec, double min_mz, double max_mz) {
            OpenMS::removePeaks(spec, min_mz, max_mz);
        }, "spectrum"_a, "min_mz"_a, "max_mz"_a, "Remove peaks outside the given m/z range")
        .def_static("subtractMinimumIntensity", [](OpenMS::MSSpectrum& spec) {
            OpenMS::subtractMinimumIntensity(spec);
        }, "spectrum"_a, "Subtract the minimum intensity from all peaks")
        ;


    // -----------------------------------------------------------------------
    // MRMFeature
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFeature, OpenMS::Feature>(m, "MRMFeature", "OpenMS class MRMFeature")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMFeature&>())
        .def("__copy__", [](const OpenMS::MRMFeature& self) { return OpenMS::MRMFeature(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFeature& self, nb::dict) { return OpenMS::MRMFeature(self); }, "memo"_a)
        .def("getScores", [](const OpenMS::MRMFeature& self) { return self.getScores(); })
        .def("setScores", &OpenMS::MRMFeature::setScores, "s"_a)
        .def("getFeature", [](const OpenMS::MRMFeature& self, const OpenMS::String& key) { return self.getFeature(key); }, "key"_a)
        .def("addFeature", [](OpenMS::MRMFeature& self, const OpenMS::Feature& f, const OpenMS::String& key) { self.addFeature(f, key); }, "f"_a, "key"_a)
        .def("getFeatures", [](const OpenMS::MRMFeature& self) { return self.getFeatures(); })
        .def("getFeatureIDs", [](OpenMS::MRMFeature& self, nb::list& output) {
            std::vector<OpenMS::String> result;
            self.getFeatureIDs(result);
            for (const auto& s : result) output.append(nb::str(s.c_str()));
        }, "output"_a, "Get feature IDs (fills the provided list)")
        .def("getFeatureIDs", [](OpenMS::MRMFeature& self) {
            std::vector<OpenMS::String> result;
            self.getFeatureIDs(result);
            return result;
        })
        .def("getPrecursorFeature", [](const OpenMS::MRMFeature& self, const OpenMS::String& key) { return self.getPrecursorFeature(key); }, "key"_a)
        .def("addPrecursorFeature", [](OpenMS::MRMFeature& self, const OpenMS::Feature& f, const OpenMS::String& key) { self.addPrecursorFeature(f, key); }, "f"_a, "key"_a)
        .def("getPrecursorFeatureIDs", [](OpenMS::MRMFeature& self) {
            std::vector<OpenMS::String> result;
            self.getPrecursorFeatureIDs(result);
            return result;
        })
        .def("__eq__", [](const OpenMS::MRMFeature& a, const OpenMS::MRMFeature& b) { return a == b; })
        .def("__ne__", [](const OpenMS::MRMFeature& a, const OpenMS::MRMFeature& b) { return a != b; })
        .def("__repr__", [](const OpenMS::MRMFeature& self) {
            std::ostringstream os;
            std::vector<OpenMS::String> ids;
            const_cast<OpenMS::MRMFeature&>(self).getFeatureIDs(ids);
            os << "MRMFeature(rt=" << self.getRT() << ", mz=" << self.getMZ()
               << ", intensity=" << self.getIntensity() << ", charge=" << self.getCharge()
               << ", num_transitions=" << ids.size() << ")";
            return os.str();
        })
        .def("__str__", [](const OpenMS::MRMFeature& self) {
            return nb::cast(self).attr("__repr__")();
        })
        ;

    // -----------------------------------------------------------------------
    // __static_* module-level wrappers for SpectrumSettings
    // -----------------------------------------------------------------------
    m.def("__static_SpectrumSettings_spectrumTypeToString", [](OpenMS::SpectrumSettings::SpectrumType type) -> OpenMS::String { return OpenMS::SpectrumSettings::spectrumTypeToString(type); }, "type"_a);
    m.def("__static_SpectrumSettings_toSpectrumType", [](const OpenMS::String& name) -> OpenMS::SpectrumSettings::SpectrumType { return OpenMS::SpectrumSettings::toSpectrumType(name); }, "name"_a);

    // -----------------------------------------------------------------------
    // __static_* module-level wrappers for SpectrumHelper
    // -----------------------------------------------------------------------
    m.def("__static_SpectrumHelper_removePeaks", [](OpenMS::MSSpectrum& spec, double min_mz, double max_mz) -> void { OpenMS::removePeaks(spec, min_mz, max_mz); }, "spectrum"_a, "min_mz"_a, "max_mz"_a);
    m.def("__static_SpectrumHelper_subtractMinimumIntensity", [](OpenMS::MSSpectrum& spec) -> void { OpenMS::subtractMinimumIntensity(spec); }, "spectrum"_a);

}
