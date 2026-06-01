// pyOpenMS nanobind bindings
// Domain: format

#include "all_casters.h"
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/FORMAT/AbsoluteQuantitationStandardsFile.h>
#include <OpenMS/FORMAT/Base64.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/CachedMzML.h>
#include <OpenMS/FORMAT/ChromeleonFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/EDTAFile.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/FLASHDeconvFeatureFile.h>
#include <OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/GNPSMetaValueFile.h>
#include <OpenMS/FORMAT/GNPSQuantificationFile.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>
#include <OpenMS/FORMAT/IBSpectraFile.h>
#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/FORMAT/InspectInfile.h>
#include <OpenMS/FORMAT/InspectOutfile.h>
#include <OpenMS/FORMAT/KroenikFile.h>
#include <OpenMS/FORMAT/MRMFile.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>
#include <OpenMS/FORMAT/MSstatsFile.h>
#include <OpenMS/FORMAT/MsInspectFile.h>
#include <OpenMS/FORMAT/MzQCFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTabM.h>
#include <OpenMS/FORMAT/MzTabMFile.h>
#include <OpenMS/FORMAT/OMSSACSVFile.h>
#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/PEFFFile.h>
#include <OpenMS/FORMAT/ParamCTDFile.h>
#include <OpenMS/FORMAT/ParquetFilter.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/FORMAT/PercolatorOutfile.h>
#include <OpenMS/FORMAT/SequestInfile.h>
#include <OpenMS/FORMAT/SequestOutfile.h>
#include <OpenMS/FORMAT/SqMassFile.h>
#include <OpenMS/FORMAT/TargetedDataFileLoader.h>
#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteHandler.h>
#include <OpenMS/FORMAT/MSPGenericFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessSqMass.h>
#include <OpenMS/FORMAT/XICParquetFile.h>
#include <OpenMS/FORMAT/XIMParquetFile.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>
#include <OpenMS/FORMAT/FeatureMapArrowIO.h>
#include <OpenMS/FORMAT/ConsensusMapArrowIO.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
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

NB_MODULE(_pyopenms_format, m) {
    m.doc() = "pyOpenMS format bindings";

    // -----------------------------------------------------------------------
    // AbsoluteQuantitationStandardsFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::AbsoluteQuantitationStandardsFile>(m, "AbsoluteQuantitationStandardsFile", "Load files containing runConcentration data")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::AbsoluteQuantitationStandardsFile& self) { return OpenMS::AbsoluteQuantitationStandardsFile(self); })
        .def("__deepcopy__", [](const OpenMS::AbsoluteQuantitationStandardsFile& self, nb::dict) { return OpenMS::AbsoluteQuantitationStandardsFile(self); }, "memo"_a)
        .def("load", [](const OpenMS::AbsoluteQuantitationStandardsFile& self, const OpenMS::String& filename, std::vector<OpenMS::AbsoluteQuantitationStandards::runConcentration>& run_concentrations) { self.load(filename, run_concentrations); }, "filename"_a, "run_concentrations"_a)
        ;

    // -----------------------------------------------------------------------
    // Base64
    // -----------------------------------------------------------------------
    auto base64_class = nb::class_<OpenMS::Base64>(m, "Base64", "Class to encode and decode Base64")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::Base64& self) { return OpenMS::Base64(self); })
        .def("__deepcopy__", [](const OpenMS::Base64& self, nb::dict) { return OpenMS::Base64(self); }, "memo"_a)
        .def_static("encodeStrings", [](const std::vector<OpenMS::String>& in, bool zlib_compression, bool append_null_byte) { OpenMS::String out; OpenMS::Base64::encodeStrings(in, out, zlib_compression, append_null_byte); return out; }, "in"_a, "zlib_compression"_a, "append_null_byte"_a, "Encodes a vector of strings to a Base64 string")
        .def_static("decodeStrings", [](const OpenMS::String& in, bool zlib_compression) { std::vector<OpenMS::String> out; OpenMS::Base64::decodeStrings(in, out, zlib_compression); return out; }, "in"_a, "zlib_compression"_a, "Decodes a Base64 string to a vector of (null-terminated) strings")
        .def_static("encode64", [](std::vector<double> in, OpenMS::Base64::ByteOrder byte_order, bool zlib_compression) {
            OpenMS::String out;
            OpenMS::Base64::encode(in, byte_order, out, zlib_compression);
            return out;
        }, "in"_a, "byte_order"_a, "zlib_compression"_a = false, "Encodes a vector of 64-bit floats (doubles) to a Base64 string")
        .def_static("decode64", [](const OpenMS::String& in, OpenMS::Base64::ByteOrder byte_order, bool zlib_compression) {
            std::vector<double> out;
            OpenMS::Base64::decode(in, byte_order, out, zlib_compression);
            return out;
        }, "in"_a, "byte_order"_a, "zlib_compression"_a = false, "Decodes a Base64 string to a vector of 64-bit floats (doubles)")
        .def_static("encode32", [](std::vector<float> in, OpenMS::Base64::ByteOrder byte_order, bool zlib_compression) {
            OpenMS::String out;
            OpenMS::Base64::encode(in, byte_order, out, zlib_compression);
            return out;
        }, "in"_a, "byte_order"_a, "zlib_compression"_a = false, "Encodes a vector of 32-bit floats to a Base64 string")
        .def_static("decode32", [](const OpenMS::String& in, OpenMS::Base64::ByteOrder byte_order, bool zlib_compression) {
            std::vector<float> out;
            OpenMS::Base64::decode(in, byte_order, out, zlib_compression);
            return out;
        }, "in"_a, "byte_order"_a, "zlib_compression"_a = false, "Decodes a Base64 string to a vector of 32-bit floats")
        .def_static("encodeIntegers", [](std::vector<int32_t> in, OpenMS::Base64::ByteOrder byte_order, bool zlib_compression) {
            OpenMS::String out;
            OpenMS::Base64::encodeIntegers(in, byte_order, out, zlib_compression);
            return out;
        }, "in"_a, "byte_order"_a, "zlib_compression"_a = false, "Encodes a vector of 32-bit integers to a Base64 string")
        .def_static("decodeIntegers", [](const OpenMS::String& in, OpenMS::Base64::ByteOrder byte_order, bool zlib_compression) {
            std::vector<int32_t> out;
            OpenMS::Base64::decodeIntegers(in, byte_order, out, zlib_compression);
            return out;
        }, "in"_a, "byte_order"_a, "zlib_compression"_a = false, "Decodes a Base64 string to a vector of 32-bit integers")
        ;
    // ByteOrder enum nested under Base64
    nb::enum_<OpenMS::Base64::ByteOrder>(base64_class, "ByteOrder", nb::is_arithmetic())
        .value("BYTEORDER_BIGENDIAN", OpenMS::Base64::ByteOrder::BYTEORDER_BIGENDIAN)
        .value("BYTEORDER_LITTLEENDIAN", OpenMS::Base64::ByteOrder::BYTEORDER_LITTLEENDIAN)
        .export_values();

#ifdef WITH_OPENTIMS
    // -----------------------------------------------------------------------
    // BrukerTimsFile
    // -----------------------------------------------------------------------
    auto brukertimsfile_class = nb::class_<OpenMS::BrukerTimsFile>(m, "BrukerTimsFile",
        R"doc(
Reader for Bruker TimsTOF .d directories via opentims.

Supports DDA-PASEF, DIA-PASEF, and raw frame-level 4D access.
Ion mobility data is stored in VSSC (1/K0) units using IM_PEAK format
for MS1 and DIA MS2, and scalar drift times for DDA MS2.

This class is only available when OpenMS is built with
``WITH_OPENTIMS=ON`` (the default). The Bruker SDK
(``timsdata.dll`` / ``libtimsdata.so``) is not bundled with pyOpenMS;
supply it via the ``OPENMS_BRUKER_SDK_PATH`` environment variable or
``Config.bruker_sdk_path``.

Use ``hasattr(pyopenms, "BrukerTimsFile")`` to feature-detect at
runtime.
)doc");
    brukertimsfile_class
        .def(nb::init<>())
        .def("load", [](OpenMS::BrukerTimsFile& self, const OpenMS::String& path) {
            OpenMS::MSExperiment exp; self.load(path, exp); return exp;
        }, "path"_a, "Load a .d directory into a new MSExperiment using default configuration")
        .def("load", [](OpenMS::BrukerTimsFile& self, const OpenMS::String& path, const OpenMS::BrukerTimsFile::Config& config) {
            OpenMS::MSExperiment exp; self.load(path, exp, config); return exp;
        }, "path"_a, "config"_a, "Load a .d directory into a new MSExperiment with explicit configuration")
        .def("readDIAMetadata", [](OpenMS::BrukerTimsFile& self, const OpenMS::String& path) {
            OpenMS::ExperimentalSettings settings;
            auto meta = self.readDIAMetadata(path, settings);
            return nb::make_tuple(meta, settings);
        }, "path"_a, "Read DIA SWATH boundaries and spectrum counts (no peak data). Returns (DIAStreamingMetadata, ExperimentalSettings).")
        .def("readDIAMetadata", [](OpenMS::BrukerTimsFile& self, const OpenMS::String& path, const OpenMS::BrukerTimsFile::Config& config) {
            OpenMS::ExperimentalSettings settings;
            auto meta = self.readDIAMetadata(path, settings, config);
            return nb::make_tuple(meta, settings);
        }, "path"_a, "config"_a, "Read DIA SWATH boundaries and spectrum counts with explicit configuration. Returns (DIAStreamingMetadata, ExperimentalSettings).")
        ;
    def_ProgressLogger<OpenMS::BrukerTimsFile>(brukertimsfile_class);

    // Nested Config struct
    auto brukertimsfile_config = nb::class_<OpenMS::BrukerTimsFile::Config>(brukertimsfile_class, "Config",
        "Processing and export configuration for BrukerTimsFile.load() and readDIAMetadata()")
        .def(nb::init<>())
        .def_rw("calibration_tolerance", &OpenMS::BrukerTimsFile::Config::calibration_tolerance, "m/z recalibration tolerance in Da (0 = default 0.1 Da)")
        .def_rw("calibrate", &OpenMS::BrukerTimsFile::Config::calibrate, "Enable m/z recalibration (off by default; may fail on some datasets)")
        .def_rw("load_ms1", &OpenMS::BrukerTimsFile::Config::load_ms1, "Load MS1 spectra. Disable for MS2-only workflows.")
        .def_rw("ms1_centroid_mz_ppm", &OpenMS::BrukerTimsFile::Config::ms1_centroid_mz_ppm, "MS1 IM-centroiding m/z tolerance in ppm (0 = disabled, suggested: 5.0)")
        .def_rw("ms1_centroid_im_pct", &OpenMS::BrukerTimsFile::Config::ms1_centroid_im_pct, "MS1 IM-centroiding ion mobility tolerance in percent (0 = disabled, suggested: 3.0)")
        .def_rw("ms1_centroid_max_peaks", &OpenMS::BrukerTimsFile::Config::ms1_centroid_max_peaks, "Upper bound on centroided peaks per MS1 spectrum")
        .def_rw("dia_ms2_n_neighbors", &OpenMS::BrukerTimsFile::Config::dia_ms2_n_neighbors, "DIA MS2 frame aggregation: number of adjacent frames on each side (0 = disabled)")
        .def_rw("dia_ms2_min_support", &OpenMS::BrukerTimsFile::Config::dia_ms2_min_support, "DIA MS2 denoising: minimum occupied neighbors in 3x3 (m/z x IM) grid")
        .def_rw("dia_ms2_centroid", &OpenMS::BrukerTimsFile::Config::dia_ms2_centroid, "DIA MS2 2D peak picking: Gaussian smoothing + local maxima detection")
        .def_rw("ms1_n_neighbors", &OpenMS::BrukerTimsFile::Config::ms1_n_neighbors, "MS1 frame aggregation: adjacent MS1 frames on each side (0 = disabled)")
        .def_rw("ms1_min_support", &OpenMS::BrukerTimsFile::Config::ms1_min_support, "MS1 denoising after aggregation: min occupied 3x3 neighbors (0 = disabled)")
        .def_rw("ms1_max_rt_distance_sec", &OpenMS::BrukerTimsFile::Config::ms1_max_rt_distance_sec, "Cap RT distance (s) between neighbor MS1 frame and center frame during aggregation (0 = no cap)")
        .def_rw("export_mode", &OpenMS::BrukerTimsFile::Config::export_mode, "AUTO detects DDA vs DIA; SPECTRUM forces per-precursor; FRAME returns raw 4D frames")
        .def_rw("tims_calibration_strategy", &OpenMS::BrukerTimsFile::Config::tims_calibration_strategy, "Strategy for converting TIMS scan indices to 1/K0 values")
        .def_rw("pressure_compensation", &OpenMS::BrukerTimsFile::Config::pressure_compensation, "Pressure compensation strategy (only effective with BRUKER_SDK calibration)")
        .def_rw("bruker_sdk_path", &OpenMS::BrukerTimsFile::Config::bruker_sdk_path, "Path to Bruker SDK library (empty = discover from OPENMS_BRUKER_SDK_PATH env var)")
        ;

    nb::enum_<OpenMS::BrukerTimsFile::Config::ExportMode>(brukertimsfile_config, "ExportMode")
        .value("AUTO", OpenMS::BrukerTimsFile::Config::ExportMode::AUTO)
        .value("SPECTRUM", OpenMS::BrukerTimsFile::Config::ExportMode::SPECTRUM)
        .value("FRAME", OpenMS::BrukerTimsFile::Config::ExportMode::FRAME)
        ;

    nb::enum_<OpenMS::BrukerTimsFile::Config::TimsCalibrationStrategy>(brukertimsfile_config, "TimsCalibrationStrategy")
        .value("AUTO", OpenMS::BrukerTimsFile::Config::TimsCalibrationStrategy::AUTO)
        .value("BRUKER_SDK", OpenMS::BrukerTimsFile::Config::TimsCalibrationStrategy::BRUKER_SDK)
        .value("RATIONAL", OpenMS::BrukerTimsFile::Config::TimsCalibrationStrategy::RATIONAL)
        .value("LINEAR", OpenMS::BrukerTimsFile::Config::TimsCalibrationStrategy::LINEAR)
        ;

    nb::enum_<OpenMS::BrukerTimsFile::Config::PressureCompensation>(brukertimsfile_config, "PressureCompensation")
        .value("NONE", OpenMS::BrukerTimsFile::Config::PressureCompensation::NONE)
        .value("GLOBAL", OpenMS::BrukerTimsFile::Config::PressureCompensation::GLOBAL)
        .value("PER_FRAME", OpenMS::BrukerTimsFile::Config::PressureCompensation::PER_FRAME)
        ;

    // Nested DIAStreamingMetadata struct
    nb::class_<OpenMS::BrukerTimsFile::DIAStreamingMetadata>(brukertimsfile_class, "DIAStreamingMetadata",
        "Metadata describing DIA SWATH layout for streaming consumers")
        .def(nb::init<>())
        .def_rw("boundaries", &OpenMS::BrukerTimsFile::DIAStreamingMetadata::boundaries, "One SwathMap per DIA window (MS2 only)")
        .def_rw("nr_ms1_spectra", &OpenMS::BrukerTimsFile::DIAStreamingMetadata::nr_ms1_spectra, "Number of MS1 frames")
        .def_rw("nr_ms2_spectra", &OpenMS::BrukerTimsFile::DIAStreamingMetadata::nr_ms2_spectra, "Per-window spectrum counts (parallel to boundaries)")
        ;
#endif // WITH_OPENTIMS

    // -----------------------------------------------------------------------
    // CachedSwathFileConsumer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CachedSwathFileConsumer>(m, "CachedSwathFileConsumer", "FullSwathFileConsumer")
        .def("__copy__", [](const OpenMS::CachedSwathFileConsumer& self) { return OpenMS::CachedSwathFileConsumer(self); })
        .def("__deepcopy__", [](const OpenMS::CachedSwathFileConsumer& self, nb::dict) { return OpenMS::CachedSwathFileConsumer(self); }, "memo"_a)
        .def(nb::init<OpenMS::String, OpenMS::String, size_t, std::vector<int>>())
        .def(nb::init<std::vector<OpenSwath::SwathMap>, OpenMS::String, OpenMS::String, size_t, std::vector<int>>())
        .def("setExpectedSize", [](OpenMS::CachedSwathFileConsumer& self, size_t p0, size_t p1) { return self.setExpectedSize(p0, p1); })
        .def("setExperimentalSettings", [](OpenMS::CachedSwathFileConsumer& self, const OpenMS::ExperimentalSettings& exp) { return self.setExperimentalSettings(exp); }, "exp"_a)
        .def("retrieveSwathMaps", [](OpenMS::CachedSwathFileConsumer& self) { std::vector<OpenSwath::SwathMap> maps; self.retrieveSwathMaps(maps); return maps; })
        .def("consumeChromatogram", [](OpenMS::CachedSwathFileConsumer& self, OpenMS::MSChromatogram& p0) { return self.consumeChromatogram(p0); })
        .def("consumeSpectrum", [](OpenMS::CachedSwathFileConsumer& self, OpenMS::MSSpectrum& s) { return self.consumeSpectrum(s); }, "s"_a)
        ;

    // -----------------------------------------------------------------------
    // CachedmzML
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CachedmzML>(m, "CachedmzML", 
        R"doc(
An class that uses on-disk caching to read and write spectra and
chromatograms
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::CachedmzML &>())
        .def(nb::init<OpenMS::String>())
        .def("__copy__", [](const OpenMS::CachedmzML& self) { return OpenMS::CachedmzML(self); })
        .def("__deepcopy__", [](const OpenMS::CachedmzML& self, nb::dict) { return OpenMS::CachedmzML(self); }, "memo"_a)
        .def("getSpectrum", [](OpenMS::CachedmzML& self, size_t id) { return self.getSpectrum(id); }, "id"_a)
        .def("getChromatogram", [](OpenMS::CachedmzML& self, size_t id) { return self.getChromatogram(id); }, "id"_a)
        .def("getNrSpectra", [](const OpenMS::CachedmzML& self) { return self.getNrSpectra(); })
        .def("getNrChromatograms", [](const OpenMS::CachedmzML& self) { return self.getNrChromatograms(); })
        .def("getMetaData", [](const OpenMS::CachedmzML& self) -> const OpenMS::MSExperiment & { return self.getMetaData(); }, nb::rv_policy::reference_internal)
        .def_static("store", [](const OpenMS::String& filename, const OpenMS::MSExperiment& map) { return OpenMS::CachedmzML::store(filename, map); }, "filename"_a, "map"_a)

        .def_static("load", [](const OpenMS::String& filename, OpenMS::CachedmzML& cached) {
            OpenMS::CachedmzML::load(filename, cached);
        }, "filename"_a, "cached"_a, "Load a cached mzML file into a CachedmzML object (backward-compatible)")
        ;

    // -----------------------------------------------------------------------
    // ChromeleonFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ChromeleonFile>(m, "ChromeleonFile", "Load Chromeleon HPLC text file and save it into a `MSExperiment`")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::ChromeleonFile& self) { return OpenMS::ChromeleonFile(self); })
        .def("__deepcopy__", [](const OpenMS::ChromeleonFile& self, nb::dict) { return OpenMS::ChromeleonFile(self); }, "memo"_a)
        .def("load", [](const OpenMS::ChromeleonFile& self, const OpenMS::String& filename) { OpenMS::MSExperiment experiment; self.load(filename, experiment); return experiment; }, "filename"_a, "Load the file's data and metadata, and save it into a `MSExperiment`")
        ;

    // -----------------------------------------------------------------------
    // DTAFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DTAFile>(m, "DTAFile", "File adapter for DTA files")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::DTAFile& self) { return OpenMS::DTAFile(self); })
        .def("__deepcopy__", [](const OpenMS::DTAFile& self, nb::dict) { return OpenMS::DTAFile(self); }, "memo"_a)
        .def("load", [](OpenMS::DTAFile& self, const OpenMS::String& filename) {
            OpenMS::MSSpectrum spectrum;
            self.load(filename, spectrum);
            return spectrum;
        }, "filename"_a, "Loads a DTA file into an MSSpectrum")
        .def("store", [](const OpenMS::DTAFile& self, const OpenMS::String& filename, const OpenMS::MSSpectrum& spectrum) {
            self.store(filename, spectrum);
        }, "filename"_a, "spectrum"_a, "Stores an MSSpectrum to a DTA file")
        ;

    // -----------------------------------------------------------------------
    // EDTAFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::EDTAFile>(m, "EDTAFile", "File adapter for Enhanced DTA files")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::EDTAFile& self) { return OpenMS::EDTAFile(self); })
        .def("__deepcopy__", [](const OpenMS::EDTAFile& self, nb::dict) { return OpenMS::EDTAFile(self); }, "memo"_a)
        .def("load", [](OpenMS::EDTAFile& self, const OpenMS::String& filename) { OpenMS::ConsensusMap consensus_map; self.load(filename, consensus_map); return consensus_map; }, "filename"_a)
        .def("store", [](const OpenMS::EDTAFile& self, const OpenMS::String& filename, const OpenMS::ConsensusMap& map) { return self.store(filename, map); }, "filename"_a, "map"_a)
        .def("store", [](const OpenMS::EDTAFile& self, const OpenMS::String& filename, const OpenMS::FeatureMap& map) { return self.store(filename, map); }, "filename"_a, "map"_a)
        ;

    // -----------------------------------------------------------------------
    // ExperimentalDesignFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ExperimentalDesignFile>(m, "ExperimentalDesignFile", 
        R"doc(
Load an experimental design from a TSV file. (see ExperimentalDesign
for details on the supported format)
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ExperimentalDesignFile &>())
        .def("__copy__", [](const OpenMS::ExperimentalDesignFile& self) { return OpenMS::ExperimentalDesignFile(self); })
        .def("__deepcopy__", [](const OpenMS::ExperimentalDesignFile& self, nb::dict) { return OpenMS::ExperimentalDesignFile(self); }, "memo"_a)
        .def_static("load", [](const OpenMS::String& tsv_file, bool require_spectra_files) {
            return OpenMS::ExperimentalDesignFile::load(tsv_file, require_spectra_files);
        }, "tsv_file"_a, "require_spectra_files"_a, "Load an experimental design from a TSV file")
        ;

    // -----------------------------------------------------------------------
    // FLASHDeconvFeatureFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FLASHDeconvFeatureFile>(m, "FLASHDeconvFeatureFile", 
        R"doc(
FLASHDeconv feature level output (.tsv, .ms1ft for Promex, .feature for TopPIC) file formats.
This class provides static methods for writing mass feature data.
Note: Methods taking std::ostream are not directly exposed. Use file-based workflows.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FLASHDeconvFeatureFile &>())
        .def("__copy__", [](const OpenMS::FLASHDeconvFeatureFile& self) { return OpenMS::FLASHDeconvFeatureFile(self); })
        .def("__deepcopy__", [](const OpenMS::FLASHDeconvFeatureFile& self, nb::dict) { return OpenMS::FLASHDeconvFeatureFile(self); }, "memo"_a)
        ;

    // -----------------------------------------------------------------------
    // FLASHDeconvSpectrumFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FLASHDeconvSpectrumFile>(m, "FLASHDeconvSpectrumFile", 
        R"doc(
FLASHDeconv Spectrum level output (.tsv, .msalign for TopPIC) file formats.
This class provides static methods for writing deconvolved spectrum data.
Note: Methods taking std::ostream are not directly exposed. Use file-based workflows.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FLASHDeconvSpectrumFile &>())
        .def("__copy__", [](const OpenMS::FLASHDeconvSpectrumFile& self) { return OpenMS::FLASHDeconvSpectrumFile(self); })
        .def("__deepcopy__", [](const OpenMS::FLASHDeconvSpectrumFile& self, nb::dict) { return OpenMS::FLASHDeconvSpectrumFile(self); }, "memo"_a)
        .def_static("writeMzML", [](const OpenMS::MSExperiment& map, std::vector<OpenMS::DeconvolvedSpectrum>& deconvolved_spectra, const OpenMS::String& deconvolved_mzML_file, const OpenMS::String& annotated_mzML_file, int mzml_charge, std::vector<double> tols) { return OpenMS::FLASHDeconvSpectrumFile::writeMzML(map, deconvolved_spectra, deconvolved_mzML_file, annotated_mzML_file, mzml_charge, tols); }, "map"_a, "deconvolved_spectra"_a, "deconvolved_mzML_file"_a, "annotated_mzML_file"_a, "mzml_charge"_a, "tols"_a)
        ;

    // -----------------------------------------------------------------------
    // FeatureFileOptions
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureFileOptions>(m, "FeatureFileOptions", "Options for loading files containing features")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::FeatureFileOptions& self) { return OpenMS::FeatureFileOptions(self); })
        .def("__deepcopy__", [](const OpenMS::FeatureFileOptions& self, nb::dict) { return OpenMS::FeatureFileOptions(self); }, "memo"_a)
        .def("setLoadConvexHull", [](OpenMS::FeatureFileOptions& self, bool convex) { return self.setLoadConvexHull(convex); }, "convex"_a, "Sets whether or not to load convex hull")
        .def("getLoadConvexHull", [](const OpenMS::FeatureFileOptions& self) { return self.getLoadConvexHull(); }, "Returns whether or not to load convex hull")
        .def("setLoadSubordinates", [](OpenMS::FeatureFileOptions& self, bool sub) { return self.setLoadSubordinates(sub); }, "sub"_a, "Sets whether or not load subordinates")
        .def("getLoadSubordinates", [](const OpenMS::FeatureFileOptions& self) { return self.getLoadSubordinates(); }, "Returns whether or not to load subordinates")
        .def("setMetadataOnly", [](OpenMS::FeatureFileOptions& self, bool only) { return self.setMetadataOnly(only); }, "only"_a, "Sets whether or not to load only meta data")
        .def("getMetadataOnly", [](const OpenMS::FeatureFileOptions& self) { return self.getMetadataOnly(); }, "Returns whether or not to load only meta data")
        .def("setSizeOnly", [](OpenMS::FeatureFileOptions& self, bool only) { return self.setSizeOnly(only); }, "only"_a, "Sets whether or not to load only feature count")
        .def("getSizeOnly", [](const OpenMS::FeatureFileOptions& self) { return self.getSizeOnly(); }, "Returns whether or not to load only meta data")
        .def("setRTRange", [](OpenMS::FeatureFileOptions& self, const OpenMS::DRange<1>& range) { return self.setRTRange(range); }, "range"_a, "Restricts the range of RT values for peaks to load")
        .def("hasRTRange", [](const OpenMS::FeatureFileOptions& self) { return self.hasRTRange(); }, "Returns true if an RT range has been set")
        .def("getRTRange", [](const OpenMS::FeatureFileOptions& self) -> const OpenMS::DRange<1> & { return self.getRTRange(); }, nb::rv_policy::reference_internal, "Returns the RT range")
        .def("setMZRange", [](OpenMS::FeatureFileOptions& self, const OpenMS::DRange<1>& range) { return self.setMZRange(range); }, "range"_a, "Restricts the range of MZ values for peaks to load")
        .def("hasMZRange", [](const OpenMS::FeatureFileOptions& self) { return self.hasMZRange(); }, "Returns true if an MZ range has been set")
        .def("getMZRange", [](const OpenMS::FeatureFileOptions& self) -> const OpenMS::DRange<1> & { return self.getMZRange(); }, nb::rv_policy::reference_internal, "Returns the MZ range")
        .def("setIntensityRange", [](OpenMS::FeatureFileOptions& self, const OpenMS::DRange<1>& range) { return self.setIntensityRange(range); }, "range"_a, "Restricts the range of intensity values for peaks to load")
        .def("hasIntensityRange", [](const OpenMS::FeatureFileOptions& self) { return self.hasIntensityRange(); }, "Returns true if an intensity range has been set")
        .def("getIntensityRange", [](const OpenMS::FeatureFileOptions& self) -> const OpenMS::DRange<1> & { return self.getIntensityRange(); }, nb::rv_policy::reference_internal, "Returns the intensity range")
        ;

    // -----------------------------------------------------------------------
    // FileHandler
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FileHandler>(m, "FileHandler", 
        R"doc(
FileHandler().loadExperiment("test.mzML", exp)
)doc")
        .def(nb::init<>())
        .def_static("getTypeByFileName", [](const OpenMS::String& filename) { return OpenMS::FileHandler::getTypeByFileName(filename); }, "filename"_a, 
            R"doc(
Determines the file type based on the file name and/or content
:param filename: Path to the file
:returns: Integer representation of the file type
)doc")
        .def_static("hasValidExtension", [](const OpenMS::String& filename, OpenMS::FileTypes::Type type) { return OpenMS::FileHandler::hasValidExtension(filename, type); }, "filename"_a, "type"_a,
            R"doc(
Checks whether the file has a valid extension for the given type
:param filename: Path to the file
:param type: The expected file type
:returns: True if the extension matches the file type
)doc")
        .def_static("stripExtension", [](const OpenMS::String& filename) { return OpenMS::FileHandler::stripExtension(filename); }, "filename"_a,
            R"doc(
Returns the file name without the extension
:param filename: Path to the file
:returns: File path without extension
)doc")
        .def_static("swapExtension", [](const OpenMS::String& filename, OpenMS::FileTypes::Type new_type) { return OpenMS::FileHandler::swapExtension(filename, new_type); }, "filename"_a, "new_type"_a,
            R"doc(
Replaces the file extension with the extension for the given type
:param filename: Path to the file
:param new_type: The new file type whose extension to use
:returns: File path with the new extension
)doc")
        .def_static("getTypeByContent", [](const OpenMS::String& filename) { return OpenMS::FileHandler::getTypeByContent(filename); }, "filename"_a,
            R"doc(
Determines the file type based on the file content
:param filename: Path to the file
:returns: The file type based on file content analysis
)doc")
        .def_static("isSupported", [](OpenMS::FileTypes::Type type) { return OpenMS::FileHandler::isSupported(type); }, "type"_a,
            R"doc(
Checks whether the given file type is supported
:param type: The file type to check
:returns: True if the file type is supported
)doc")
        .def("getOptions", [](OpenMS::FileHandler& self) -> OpenMS::PeakFileOptions & { return self.getOptions(); }, nb::rv_policy::reference_internal, "Access to the options for loading/storing")
        .def("setOptions", [](OpenMS::FileHandler& self, const OpenMS::PeakFileOptions& p0) { return self.setOptions(p0); }, "Sets options for loading/storing")
        .def("loadFeatures", [](OpenMS::FileHandler& self, const OpenMS::String& filename) { OpenMS::FeatureMap map; self.loadFeatures(filename, map); return map; }, "filename"_a,
            R"doc(
Loads features from a file into a FeatureMap (auto-detect file type)
:param filename: The name of the file to load features from
:returns: FeatureMap containing the loaded features
:raises:
Exception: FileNotFound is thrown if the file could not be opened
)doc")
        .def("loadFeatures", [](OpenMS::FileHandler& self, const OpenMS::String& filename, std::vector<OpenMS::FileTypes::Type> allowed_types, OpenMS::ProgressLogger::LogType log) { OpenMS::FeatureMap map; self.loadFeatures(filename, map, allowed_types, log); return map; }, "filename"_a, "allowed_types"_a, "log"_a,
            R"doc(
Loads features from a file into a FeatureMap
:param filename: The name of the file to load features from
:param allowed_types: List of allowed file types
:param log: Progress logging mode
:returns: FeatureMap containing the loaded features
:raises:
Exception: FileNotFound is thrown if the file could not be opened
)doc")
        .def_static("computeFileHash", [](const OpenMS::String& filename) { return OpenMS::FileHandler::computeFileHash(filename); }, "filename"_a,
            R"doc(
Computes a SHA-1 hash of the file content
:param filename: Path to the file
:returns: SHA-1 hash string of the file content
)doc")

        .def("loadExperiment", [](OpenMS::FileHandler& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp) {
            nb::gil_scoped_release release;
            self.loadExperiment(filename, exp);
        }, "filename"_a, "exp"_a, "Load experiment from file")
        .def("loadExperiment", [](OpenMS::FileHandler& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log,
             bool rewrite_source_file, bool compute_hash) {
            nb::gil_scoped_release release;
            self.loadExperiment(filename, exp, allowed_types, log, rewrite_source_file, compute_hash);
        }, "filename"_a, "exp"_a, "allowed_types"_a, "log"_a, "rewrite_source_file"_a = false, "compute_hash"_a = false, "Load experiment with options")

        .def("storeExperiment", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::MSExperiment& exp) {
            nb::gil_scoped_release release;
            self.storeExperiment(filename, exp);
        }, "filename"_a, "exp"_a, "Store experiment to file")
        .def("storeExperiment", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::MSExperiment& exp,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log) {
            nb::gil_scoped_release release;
            self.storeExperiment(filename, exp, allowed_types, log);
        }, "filename"_a, "exp"_a, "allowed_types"_a, "log"_a, "Store experiment with options")

        .def_static("getType", [](const OpenMS::String& filename) {
            return OpenMS::FileHandler::getType(filename);
        }, "filename"_a, "Determine the file type from the file name")

        .def_static("getConsistentOutputfileType", [](const OpenMS::String& output_filename, const OpenMS::String& requested_type) {
            return OpenMS::FileHandler::getConsistentOutputfileType(output_filename, requested_type);
        }, "output_filename"_a, "requested_type"_a,
           "Checks consistency of output file type from filename and requested type. Returns consistent type or UNKNOWN on conflict")

        .def("getFeatOptions", [](OpenMS::FileHandler& self) -> OpenMS::FeatureFileOptions& { return self.getFeatOptions(); }, nb::rv_policy::reference_internal,
            "Mutable access to the feature file options for loading/storing")
        .def("setFeatOptions", [](OpenMS::FileHandler& self, const OpenMS::FeatureFileOptions& opts) { self.setFeatOptions(opts); }, "options"_a,
            "Set feature file options for loading/storing")

        .def("loadSpectrum", [](OpenMS::FileHandler& self, const OpenMS::String& filename) {
            OpenMS::MSSpectrum spec;
            self.loadSpectrum(filename, spec);
            return spec;
        }, "filename"_a, "Load a single MSSpectrum from a file")
        .def("loadSpectrum", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const std::vector<OpenMS::FileTypes::Type>& allowed_types) {
            OpenMS::MSSpectrum spec;
            self.loadSpectrum(filename, spec, allowed_types);
            return spec;
        }, "filename"_a, "allowed_types"_a, "Load a single MSSpectrum from a file with allowed types filter")

        .def("storeSpectrum", [](OpenMS::FileHandler& self, const OpenMS::String& filename, OpenMS::MSSpectrum& spec) {
            self.storeSpectrum(filename, spec);
        }, "filename"_a, "spec"_a, "Store a single MSSpectrum to a file")
        .def("storeSpectrum", [](OpenMS::FileHandler& self, const OpenMS::String& filename, OpenMS::MSSpectrum& spec, const std::vector<OpenMS::FileTypes::Type>& allowed_types) {
            self.storeSpectrum(filename, spec, allowed_types);
        }, "filename"_a, "spec"_a, "allowed_types"_a, "Store a single MSSpectrum to a file with allowed types filter")

        .def("storeFeatures", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::FeatureMap& map) {
            self.storeFeatures(filename, map);
        }, "filename"_a, "map"_a, "Store a FeatureMap to a file")
        .def("storeFeatures", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::FeatureMap& map,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log) {
            self.storeFeatures(filename, map, allowed_types, log);
        }, "filename"_a, "map"_a, "allowed_types"_a, "log"_a, "Store a FeatureMap to a file with options")

        .def("loadConsensusFeatures", [](OpenMS::FileHandler& self, const OpenMS::String& filename) {
            OpenMS::ConsensusMap map;
            self.loadConsensusFeatures(filename, map);
            return map;
        }, "filename"_a, "Load a ConsensusMap from a file")
        .def("loadConsensusFeatures", [](OpenMS::FileHandler& self, const OpenMS::String& filename,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log) {
            OpenMS::ConsensusMap map;
            self.loadConsensusFeatures(filename, map, allowed_types, log);
            return map;
        }, "filename"_a, "allowed_types"_a, "log"_a, "Load a ConsensusMap from a file with options")

        .def("storeConsensusFeatures", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::ConsensusMap& map) {
            self.storeConsensusFeatures(filename, map);
        }, "filename"_a, "map"_a, "Store a ConsensusMap to a file")
        .def("storeConsensusFeatures", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::ConsensusMap& map,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log) {
            self.storeConsensusFeatures(filename, map, allowed_types, log);
        }, "filename"_a, "map"_a, "allowed_types"_a, "log"_a, "Store a ConsensusMap to a file with options")

        .def("loadIdentifications", [](OpenMS::FileHandler& self, const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            self.loadIdentifications(filename, proteins, peptides);
            return nb::make_tuple(proteins, peptides);
        }, "filename"_a, "Load identifications from file. Returns (list[ProteinIdentification], PeptideIdentificationList)")
        .def("loadIdentifications", [](OpenMS::FileHandler& self, const OpenMS::String& filename,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log) {
            std::vector<OpenMS::ProteinIdentification> proteins;
            OpenMS::PeptideIdentificationList peptides;
            self.loadIdentifications(filename, proteins, peptides, allowed_types, log);
            return nb::make_tuple(proteins, peptides);
        }, "filename"_a, "allowed_types"_a, "log"_a, "Load identifications from file with options. Returns (list[ProteinIdentification], PeptideIdentificationList)")

        .def("storeIdentifications", [](OpenMS::FileHandler& self, const OpenMS::String& filename,
             const std::vector<OpenMS::ProteinIdentification>& proteins, const OpenMS::PeptideIdentificationList& peptides) {
            self.storeIdentifications(filename, proteins, peptides);
        }, "filename"_a, "proteins"_a, "peptides"_a, "Store identifications to file")
        .def("storeIdentifications", [](OpenMS::FileHandler& self, const OpenMS::String& filename,
             const std::vector<OpenMS::ProteinIdentification>& proteins, const OpenMS::PeptideIdentificationList& peptides,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log) {
            self.storeIdentifications(filename, proteins, peptides, allowed_types, log);
        }, "filename"_a, "proteins"_a, "peptides"_a, "allowed_types"_a, "log"_a, "Store identifications to file with options")

        .def("loadTransitions", [](OpenMS::FileHandler& self, const OpenMS::String& filename) {
            OpenMS::TargetedExperiment library;
            self.loadTransitions(filename, library);
            return library;
        }, "filename"_a, "Load targeted experiment transitions from file")
        .def("loadTransitions", [](OpenMS::FileHandler& self, const OpenMS::String& filename,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log) {
            OpenMS::TargetedExperiment library;
            self.loadTransitions(filename, library, allowed_types, log);
            return library;
        }, "filename"_a, "allowed_types"_a, "log"_a, "Load targeted experiment transitions from file with options")

        .def("storeTransitions", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::TargetedExperiment& library) {
            self.storeTransitions(filename, library);
        }, "filename"_a, "library"_a, "Store targeted experiment transitions to file")
        .def("storeTransitions", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::TargetedExperiment& library,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types, OpenMS::ProgressLogger::LogType log) {
            self.storeTransitions(filename, library, allowed_types, log);
        }, "filename"_a, "library"_a, "allowed_types"_a, "log"_a, "Store targeted experiment transitions to file with options")

        .def("loadTransformations", [](OpenMS::FileHandler& self, const OpenMS::String& filename, bool fit_model) {
            OpenMS::TransformationDescription map;
            self.loadTransformations(filename, map, fit_model);
            return map;
        }, "filename"_a, "fit_model"_a = true, "Load transformation description from file")
        .def("loadTransformations", [](OpenMS::FileHandler& self, const OpenMS::String& filename, bool fit_model,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types) {
            OpenMS::TransformationDescription map;
            self.loadTransformations(filename, map, fit_model, allowed_types);
            return map;
        }, "filename"_a, "fit_model"_a, "allowed_types"_a, "Load transformation description from file with allowed types filter")

        .def("storeTransformations", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::TransformationDescription& map) {
            self.storeTransformations(filename, map);
        }, "filename"_a, "map"_a, "Store transformation description to file")
        .def("storeTransformations", [](OpenMS::FileHandler& self, const OpenMS::String& filename, const OpenMS::TransformationDescription& map,
             const std::vector<OpenMS::FileTypes::Type>& allowed_types) {
            self.storeTransformations(filename, map, allowed_types);
        }, "filename"_a, "map"_a, "allowed_types"_a, "Store transformation description to file with allowed types filter")
        ;

    // -----------------------------------------------------------------------
    // FileTypes
    // -----------------------------------------------------------------------
    auto filetypes_class = nb::class_<OpenMS::FileTypes>(m, "FileTypes", "OpenMS class FileTypes")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::FileTypes &>())
        .def("__copy__", [](const OpenMS::FileTypes& self) { return OpenMS::FileTypes(self); })
        .def("__deepcopy__", [](const OpenMS::FileTypes& self, nb::dict) { return OpenMS::FileTypes(self); }, "memo"_a)
        .def_static("typeToName", [](OpenMS::FileTypes::Type type) { return OpenMS::FileTypes::typeToName(type); }, "type"_a, "Returns the name/extension of the type")
        .def_static("typeToDescription", [](OpenMS::FileTypes::Type type) { return OpenMS::FileTypes::typeToDescription(type); }, "type"_a, "Returns the human-readable explanation of the type")
        .def_static("nameToType", [](const OpenMS::String& name) { return OpenMS::FileTypes::nameToType(name); }, "name"_a)
        .def_static("typeToMZML", [](OpenMS::FileTypes::Type type) { return OpenMS::FileTypes::typeToMZML(type); }, "type"_a, "Returns the mzML name")
        ;
    // FileType enum nested under FileTypes
    nb::enum_<OpenMS::FileTypes::Type>(filetypes_class, "FileType", "Enum for known file types", nb::is_arithmetic())
        .value("UNKNOWN", OpenMS::FileTypes::Type::UNKNOWN)
        .value("DTA", OpenMS::FileTypes::Type::DTA)
        .value("DTA2D", OpenMS::FileTypes::Type::DTA2D)
        .value("MZDATA", OpenMS::FileTypes::Type::MZDATA)
        .value("MZXML", OpenMS::FileTypes::Type::MZXML)
        .value("FEATUREXML", OpenMS::FileTypes::Type::FEATUREXML)
        .value("IDXML", OpenMS::FileTypes::Type::IDXML)
        .value("CONSENSUSXML", OpenMS::FileTypes::Type::CONSENSUSXML)
        .value("MGF", OpenMS::FileTypes::Type::MGF)
        .value("INI", OpenMS::FileTypes::Type::INI)
        .value("TOPPAS", OpenMS::FileTypes::Type::TOPPAS)
        .value("TRANSFORMATIONXML", OpenMS::FileTypes::Type::TRANSFORMATIONXML)
        .value("MZML", OpenMS::FileTypes::Type::MZML)
        .value("CACHEDMZML", OpenMS::FileTypes::Type::CACHEDMZML)
        .value("MS2", OpenMS::FileTypes::Type::MS2)
        .value("PEPXML", OpenMS::FileTypes::Type::PEPXML)
        .value("PROTXML", OpenMS::FileTypes::Type::PROTXML)
        .value("MZIDENTML", OpenMS::FileTypes::Type::MZIDENTML)
        .value("QCML", OpenMS::FileTypes::Type::QCML)
        .value("MZQC", OpenMS::FileTypes::Type::MZQC)
        .value("GELML", OpenMS::FileTypes::Type::GELML)
        .value("TRAML", OpenMS::FileTypes::Type::TRAML)
        .value("MSP", OpenMS::FileTypes::Type::MSP)
        .value("OMSSAXML", OpenMS::FileTypes::Type::OMSSAXML)
        .value("MASCOTXML", OpenMS::FileTypes::Type::MASCOTXML)
        .value("PNG", OpenMS::FileTypes::Type::PNG)
        .value("XMASS", OpenMS::FileTypes::Type::XMASS)
        .value("TSV", OpenMS::FileTypes::Type::TSV)
        .value("MZTAB", OpenMS::FileTypes::Type::MZTAB)
        .value("PEPLIST", OpenMS::FileTypes::Type::PEPLIST)
        .value("HARDKLOER", OpenMS::FileTypes::Type::HARDKLOER)
        .value("KROENIK", OpenMS::FileTypes::Type::KROENIK)
        .value("FASTA", OpenMS::FileTypes::Type::FASTA)
        .value("EDTA", OpenMS::FileTypes::Type::EDTA)
        .value("CSV", OpenMS::FileTypes::Type::CSV)
        .value("TXT", OpenMS::FileTypes::Type::TXT)
        .value("OBO", OpenMS::FileTypes::Type::OBO)
        .value("HTML", OpenMS::FileTypes::Type::HTML)
        .value("ANALYSISXML", OpenMS::FileTypes::Type::ANALYSISXML)
        .value("XSD", OpenMS::FileTypes::Type::XSD)
        .value("PSQ", OpenMS::FileTypes::Type::PSQ)
        .value("MRM", OpenMS::FileTypes::Type::MRM)
        .value("SQMASS", OpenMS::FileTypes::Type::SQMASS)
        .value("PQP", OpenMS::FileTypes::Type::PQP)
        .value("MS", OpenMS::FileTypes::Type::MS)
        .value("OSW", OpenMS::FileTypes::Type::OSW)
        .value("PSMS", OpenMS::FileTypes::Type::PSMS)
        .value("PIN", OpenMS::FileTypes::Type::PIN)
        .value("PARAMXML", OpenMS::FileTypes::Type::PARAMXML)
        .value("SPLIB", OpenMS::FileTypes::Type::SPLIB)
        .value("NOVOR", OpenMS::FileTypes::Type::NOVOR)
        .value("XQUESTXML", OpenMS::FileTypes::Type::XQUESTXML)
        .value("SPECXML", OpenMS::FileTypes::Type::SPECXML)
        .value("JSON", OpenMS::FileTypes::Type::JSON)
        .value("RAW", OpenMS::FileTypes::Type::RAW)
        .value("OMS", OpenMS::FileTypes::Type::OMS)
        .value("EXE", OpenMS::FileTypes::Type::EXE)
        .value("XML", OpenMS::FileTypes::Type::XML)
        .value("BZ2", OpenMS::FileTypes::Type::BZ2)
        .value("GZ", OpenMS::FileTypes::Type::GZ)
        .value("PARQUET", OpenMS::FileTypes::Type::PARQUET)
        .value("SIZE_OF_TYPE", OpenMS::FileTypes::Type::SIZE_OF_TYPE)

        .export_values();

    // -----------------------------------------------------------------------
    // GNPSMetaValueFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::GNPSMetaValueFile>(m, "GNPSMetaValueFile", "OpenMS class GNPSMetaValueFile")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::GNPSMetaValueFile &>())
        .def("__copy__", [](const OpenMS::GNPSMetaValueFile& self) { return OpenMS::GNPSMetaValueFile(self); })
        .def("__deepcopy__", [](const OpenMS::GNPSMetaValueFile& self, nb::dict) { return OpenMS::GNPSMetaValueFile(self); }, "memo"_a)
        .def_static("store", [](const OpenMS::ConsensusMap& consensus_map, const OpenMS::String& output_file) { return OpenMS::GNPSMetaValueFile::store(consensus_map, output_file); }, "consensus_map"_a, "output_file"_a)
        ;

    // -----------------------------------------------------------------------
    // GNPSQuantificationFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::GNPSQuantificationFile>(m, "GNPSQuantificationFile", "OpenMS class GNPSQuantificationFile")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::GNPSQuantificationFile &>())
        .def("__copy__", [](const OpenMS::GNPSQuantificationFile& self) { return OpenMS::GNPSQuantificationFile(self); })
        .def("__deepcopy__", [](const OpenMS::GNPSQuantificationFile& self, nb::dict) { return OpenMS::GNPSQuantificationFile(self); }, "memo"_a)
        .def_static("store", [](const OpenMS::ConsensusMap& consensus_map, const OpenMS::String& output_file) { return OpenMS::GNPSQuantificationFile::store(consensus_map, output_file); }, "consensus_map"_a, "output_file"_a)
        ;

    // -----------------------------------------------------------------------
    // IBSpectraFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IBSpectraFile>(m, "IBSpectraFile", 
        R"doc(
Implements the export of consensusmaps into the IBSpectra format used
by isobar to load quantification results
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IBSpectraFile &>())
        .def("__copy__", [](const OpenMS::IBSpectraFile& self) { return OpenMS::IBSpectraFile(self); })
        .def("__deepcopy__", [](const OpenMS::IBSpectraFile& self, nb::dict) { return OpenMS::IBSpectraFile(self); }, "memo"_a)
        .def("store", [](OpenMS::IBSpectraFile& self, const OpenMS::String& filename, const OpenMS::ConsensusMap& cm) { return self.store(filename, cm); }, "filename"_a, "cm"_a)
        ;

    // -----------------------------------------------------------------------
    // IndexedMzMLDecoder
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IndexedMzMLDecoder>(m, "IndexedMzMLDecoder", 
        R"doc(
A class to analyze indexedmzML files and extract the offsets of individual tags
Specifically, this class allows one to extract the offsets of the <indexList>
tag and of all <spectrum> and <chromatogram> tag using the indices found at
the end of the indexedmzML XML structure
While findIndexListOffset tries extracts the offset of the indexList tag from
the last 1024 bytes of the file, this offset allows the function parseOffsets
to extract all elements contained in the <indexList> tag and thus get access
to all spectra and chromatogram offsets
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::IndexedMzMLDecoder &>())
        .def("__copy__", [](const OpenMS::IndexedMzMLDecoder& self) { return OpenMS::IndexedMzMLDecoder(self); })
        .def("__deepcopy__", [](const OpenMS::IndexedMzMLDecoder& self, nb::dict) { return OpenMS::IndexedMzMLDecoder(self); }, "memo"_a)
        .def("findIndexListOffset", [](OpenMS::IndexedMzMLDecoder& self, const OpenMS::String& filename, int buffersize) { return self.findIndexListOffset(filename, buffersize); }, "filename"_a, "buffersize"_a = 1023)
        ;

    // -----------------------------------------------------------------------
    // IndexedMzMLFileLoader
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::IndexedMzMLFileLoader>(m, "IndexedMzMLFileLoader", "A class to load an indexedmzML file")
        .def(nb::init<>())
        .def("getOptions", [](OpenMS::IndexedMzMLFileLoader& self) -> OpenMS::PeakFileOptions & { return self.getOptions(); }, nb::rv_policy::reference_internal, "Returns the options for loading/storing")
        .def("setOptions", [](OpenMS::IndexedMzMLFileLoader& self, const OpenMS::PeakFileOptions& p0) { return self.setOptions(p0); }, "Returns the options for loading/storing")
        .def("load", [](OpenMS::IndexedMzMLFileLoader& self, const OpenMS::String& filename, OpenMS::OnDiscMSExperiment& exp) { nb::gil_scoped_release release; return self.load(filename, exp); }, "filename"_a, "exp"_a)
        .def("store", [](OpenMS::IndexedMzMLFileLoader& self, const OpenMS::String& filename, OpenMS::OnDiscMSExperiment& exp) { nb::gil_scoped_release release; return self.store(filename, exp); }, "filename"_a, "exp"_a,
            R"doc(
Load a file\n
Tries to parse the file, success needs to be checked with the return value
:param filename: Filename determines where the file is located
:param exp: Object which will contain the data after the call
:return: Indicates whether parsing was successful (if it is false, the file most likely was not an mzML or not indexed)
)doc")
        .def("store", [](OpenMS::IndexedMzMLFileLoader& self, const OpenMS::String& filename, OpenMS::MSExperiment& exp) { nb::gil_scoped_release release; return self.store(filename, exp); }, "filename"_a, "exp"_a,
            R"doc(
Load a file\n
Tries to parse the file, success needs to be checked with the return value
:param filename: Filename determines where the file is located
:param exp: Object which will contain the data after the call
:return: Indicates whether parsing was successful (if it is false, the file most likely was not an mzML or not indexed)
)doc")
        ;

    // -----------------------------------------------------------------------
    // IndexedMzMLHandler
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Internal::IndexedMzMLHandler>(m, "IndexedMzMLHandler", "A low-level class to read an indexedmzML file")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Internal::IndexedMzMLHandler &>())
        .def(nb::init<OpenMS::String>())
        .def("__copy__", [](const OpenMS::Internal::IndexedMzMLHandler& self) { return OpenMS::Internal::IndexedMzMLHandler(self); })
        .def("__deepcopy__", [](const OpenMS::Internal::IndexedMzMLHandler& self, nb::dict) { return OpenMS::Internal::IndexedMzMLHandler(self); }, "memo"_a)
        .def("openFile", [](OpenMS::Internal::IndexedMzMLHandler& self, const OpenMS::String& filename) { return self.openFile(filename); }, "filename"_a)
        .def("getParsingSuccess", [](const OpenMS::Internal::IndexedMzMLHandler& self) { return self.getParsingSuccess(); })
        .def("getNrSpectra", [](const OpenMS::Internal::IndexedMzMLHandler& self) { return self.getNrSpectra(); })
        .def("getNrChromatograms", [](const OpenMS::Internal::IndexedMzMLHandler& self) { return self.getNrChromatograms(); })
        .def("getSpectrumById", [](OpenMS::Internal::IndexedMzMLHandler& self, int id) { return self.getSpectrumById(id); }, "id"_a)
        .def("getMSSpectrumById", [](OpenMS::Internal::IndexedMzMLHandler& self, int id) { return self.getMSSpectrumById(id); }, "id"_a)
        .def("getMSSpectrumByNativeId", [](OpenMS::Internal::IndexedMzMLHandler& self, const OpenMS::String& id) { OpenMS::MSSpectrum s; self.getMSSpectrumByNativeId(id, s); return s; }, "id"_a)
        .def("getMSSpectrumById", [](OpenMS::Internal::IndexedMzMLHandler& self, int id, OpenMS::MSSpectrum& s) { return self.getMSSpectrumById(id, s); }, "id"_a, "s"_a)
        .def("getChromatogramById", [](OpenMS::Internal::IndexedMzMLHandler& self, int id) { return self.getChromatogramById(id); }, "id"_a)
        .def("getMSChromatogramById", [](OpenMS::Internal::IndexedMzMLHandler& self, int id) { return self.getMSChromatogramById(id); }, "id"_a)
        .def("getMSChromatogramByNativeId", [](OpenMS::Internal::IndexedMzMLHandler& self, const OpenMS::String& id) { OpenMS::MSChromatogram c; self.getMSChromatogramByNativeId(id, c); return c; }, "id"_a)
        .def("getMSChromatogramById", [](OpenMS::Internal::IndexedMzMLHandler& self, int id, OpenMS::MSChromatogram& c) { return self.getMSChromatogramById(id, c); }, "id"_a, "c"_a)
        .def("setSkipXMLChecks", [](OpenMS::Internal::IndexedMzMLHandler& self, bool skip) { return self.setSkipXMLChecks(skip); }, "skip"_a)
        ;

    // -----------------------------------------------------------------------
    // InspectInfile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::InspectInfile>(m, "InspectInfile", "Inspect input file adapter")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::InspectInfile &>())
        .def("__copy__", [](const OpenMS::InspectInfile& self) { return OpenMS::InspectInfile(self); })
        .def("__deepcopy__", [](const OpenMS::InspectInfile& self, nb::dict) { return OpenMS::InspectInfile(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("store", [](OpenMS::InspectInfile& self, const OpenMS::String& filename) { return self.store(filename); }, "filename"_a, "Stores the experiment data in an Inspect input file that can be used as input for Inspect shell execution")
        .def("handlePTMs", [](OpenMS::InspectInfile& self, const OpenMS::String& modification_line, const OpenMS::String& modifications_filename, bool monoisotopic) { return self.handlePTMs(modification_line, modifications_filename, monoisotopic); }, "modification_line"_a, "modifications_filename"_a, "monoisotopic"_a)
        .def("getSpectra", [](const OpenMS::InspectInfile& self) { return self.getSpectra(); }, "Specifies a spectrum file to search")
        .def("setSpectra", [](OpenMS::InspectInfile& self, const OpenMS::String& spectra) { return self.setSpectra(spectra); }, "spectra"_a, "Specifies a spectrum file to search")
        .def("getDb", [](const OpenMS::InspectInfile& self) { return self.getDb(); }, "Specifies the name of a database (.trie file) to search")
        .def("setDb", [](OpenMS::InspectInfile& self, const OpenMS::String& db) { return self.setDb(db); }, "db"_a, "Specifies the name of a database (.trie file) to search")
        .def("getEnzyme", [](const OpenMS::InspectInfile& self) { return self.getEnzyme(); }, 
            R"doc(
Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values
)doc")
        .def("setEnzyme", [](OpenMS::InspectInfile& self, const OpenMS::String& enzyme) { return self.setEnzyme(enzyme); }, "enzyme"_a, 
            R"doc(
Specifies the name of a enzyme. "Trypsin", "None", and "Chymotrypsin" are the available values
)doc")
        .def("getModificationsPerPeptide", [](const OpenMS::InspectInfile& self) { return self.getModificationsPerPeptide(); }, "Number of PTMs permitted in a single peptide")
        .def("setModificationsPerPeptide", [](OpenMS::InspectInfile& self, int modifications_per_peptide) { return self.setModificationsPerPeptide(modifications_per_peptide); }, "modifications_per_peptide"_a, "Number of PTMs permitted in a single peptide")
        .def("getBlind", [](const OpenMS::InspectInfile& self) { return self.getBlind(); }, "Run inspect in a blind mode")
        .def("setBlind", [](OpenMS::InspectInfile& self, unsigned int blind) { return self.setBlind(blind); }, "blind"_a, "Run inspect in a blind mode")
        .def("getMaxPTMsize", [](const OpenMS::InspectInfile& self) { return self.getMaxPTMsize(); }, "The maximum modification size (in Da) to consider in a blind search")
        .def("setMaxPTMsize", [](OpenMS::InspectInfile& self, float maxptmsize) { return self.setMaxPTMsize(maxptmsize); }, "maxptmsize"_a, "The maximum modification size (in Da) to consider in a blind search")
        .def("getPrecursorMassTolerance", [](const OpenMS::InspectInfile& self) { return self.getPrecursorMassTolerance(); }, "Specifies the parent mass tolerance, in Daltons")
        .def("setPrecursorMassTolerance", [](OpenMS::InspectInfile& self, float precursor_mass_tolerance) { return self.setPrecursorMassTolerance(precursor_mass_tolerance); }, "precursor_mass_tolerance"_a, "Specifies the parent mass tolerance, in Daltons")
        .def("getPeakMassTolerance", [](const OpenMS::InspectInfile& self) { return self.getPeakMassTolerance(); }, "How far b and y peaks can be shifted from their expected masses.")
        .def("setPeakMassTolerance", [](OpenMS::InspectInfile& self, float peak_mass_tolerance) { return self.setPeakMassTolerance(peak_mass_tolerance); }, "peak_mass_tolerance"_a, "How far b and y peaks can be shifted from their expected masses")
        .def("getMulticharge", [](const OpenMS::InspectInfile& self) { return self.getMulticharge(); }, "If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible")
        .def("setMulticharge", [](OpenMS::InspectInfile& self, unsigned int multicharge) { return self.setMulticharge(multicharge); }, "multicharge"_a, "If set to true, attempt to guess the precursor charge and mass, and consider multiple charge states if feasible")
        .def("getInstrument", [](const OpenMS::InspectInfile& self) { return self.getInstrument(); }, "If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass")
        .def("setInstrument", [](OpenMS::InspectInfile& self, const OpenMS::String& instrument) { return self.setInstrument(instrument); }, "instrument"_a, "If set to QTOF, uses a QTOF-derived fragmentation model, and does not attempt to correct the parent mass")
        .def("getTagCount", [](const OpenMS::InspectInfile& self) { return self.getTagCount(); }, "Number of tags to generate")
        .def("setTagCount", [](OpenMS::InspectInfile& self, int TagCount) { return self.setTagCount(TagCount); }, "TagCount"_a, "Number of tags to generate")
        .def("getModifications", [](const OpenMS::InspectInfile& self) { return self.getModifications(); }, "Returns the modifications")
        ;

    // -----------------------------------------------------------------------
    // InspectOutfile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::InspectOutfile>(m, "InspectOutfile", "Representation of an Inspect outfile")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::InspectOutfile &>())
        .def("__copy__", [](const OpenMS::InspectOutfile& self) { return OpenMS::InspectOutfile(self); })
        .def("__deepcopy__", [](const OpenMS::InspectOutfile& self, nb::dict) { return OpenMS::InspectOutfile(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("load", [](OpenMS::InspectOutfile& self, const OpenMS::String& result_filename, OpenMS::PeptideIdentificationList& peptide_identifications, OpenMS::ProteinIdentification& protein_identification, double p_value_threshold, const OpenMS::String& database_filename) { return self.load(result_filename, peptide_identifications, protein_identification, p_value_threshold, database_filename); }, "result_filename"_a, "peptide_identifications"_a, "protein_identification"_a, "p_value_threshold"_a, "database_filename"_a = "")
        .def("getWantedRecords", [](OpenMS::InspectOutfile& self, const OpenMS::String& result_filename, double p_value_threshold) { return self.getWantedRecords(result_filename, p_value_threshold); }, "result_filename"_a, "p_value_threshold"_a, 
            R"doc(
Load the results of an Inspect search
:param result_filename: Input parameter which is the file name of the input file
:param peptide_identifications: Output parameter which holds the peptide identifications from the given file
:param protein_identification: Output parameter which holds the protein identifications from the given file
:param p_value_threshold:
:param database_filename:
:raises:
Exception: FileNotFound is thrown if the given file could not be found
:raises:
Exception: ParseError is thrown if the given file could not be parsed
:raises:
Exception: FileEmpty is thrown if the given file is empty
)doc")
        .def("compressTrieDB", [](OpenMS::InspectOutfile& self, const OpenMS::String& database_filename, const OpenMS::String& index_filename, std::vector<size_t>& wanted_records, const OpenMS::String& snd_database_filename, const OpenMS::String& snd_index_filename, bool append) { return self.compressTrieDB(database_filename, index_filename, wanted_records, snd_database_filename, snd_index_filename, append); }, "database_filename"_a, "index_filename"_a, "wanted_records"_a, "snd_database_filename"_a, "snd_index_filename"_a, "append"_a = false, "Generates a trie database from another one, using the wanted records only")
        .def("generateTrieDB", [](OpenMS::InspectOutfile& self, const OpenMS::String& source_database_filename, const OpenMS::String& database_filename, const OpenMS::String& index_filename, bool append, const OpenMS::String& species) { return self.generateTrieDB(source_database_filename, database_filename, index_filename, append, species); }, "source_database_filename"_a, "database_filename"_a, "index_filename"_a, "append"_a = false, "species"_a = "", "Generates a trie database from a given one (the type of database is determined by getLabels)")
        .def("getACAndACType", [](OpenMS::InspectOutfile& self, OpenMS::String line) { OpenMS::String accession; OpenMS::String accession_type; self.getACAndACType(line, accession, accession_type); return std::make_tuple(accession, accession_type); }, "line"_a, "Retrieve the accession type and accession number from a protein description line")
        .def("getLabels", [](OpenMS::InspectOutfile& self, const OpenMS::String& source_database_filename) { OpenMS::String ac_label, sequence_start_label, sequence_end_label, comment_label, species_label; self.getLabels(source_database_filename, ac_label, sequence_start_label, sequence_end_label, comment_label, species_label); return std::make_tuple(ac_label, sequence_start_label, sequence_end_label, comment_label, species_label); }, "source_database_filename"_a, "Retrieve the labels of a given database (at the moment FASTA and Swissprot)")
        .def("getSequences", [](OpenMS::InspectOutfile& self, const OpenMS::String& database_filename, const std::map<size_t, size_t>& wanted_records) { std::vector<OpenMS::String> sequences; self.getSequences(database_filename, wanted_records, sequences); return sequences; }, "database_filename"_a, "wanted_records"_a, "Retrieve sequences from a trie database")
        .def("getExperiment", [](OpenMS::InspectOutfile& self, const OpenMS::String& in_filename) { OpenMS::MSExperiment exp; OpenMS::String type; self.getExperiment(exp, type, in_filename); return std::make_tuple(exp, type); }, "in_filename"_a, "Get the experiment from a file")
        .def("getSearchEngineAndVersion", [](OpenMS::InspectOutfile& self, const OpenMS::String& cmd_output, OpenMS::ProteinIdentification& protein_identification) { return self.getSearchEngineAndVersion(cmd_output, protein_identification); }, "cmd_output"_a, "protein_identification"_a, "Get the search engine and its version from the output of the InsPecT executable without parameters. Returns true on success, false otherwise")
        .def("readOutHeader", [](OpenMS::InspectOutfile& self, const OpenMS::String& filename, const OpenMS::String& header_line) {
            OpenMS::Int spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column;
            OpenMS::Size number_of_columns;
            self.readOutHeader(filename, header_line, spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column, number_of_columns);
            return nb::make_tuple(spectrum_file_column, scan_column, peptide_column, protein_column, charge_column, MQ_score_column, p_value_column, record_number_column, DB_file_pos_column, spec_file_pos_column, number_of_columns);
        }, "filename"_a, "header_line"_a, "Reads the header of an inspect output file and returns column indices as a tuple")
        ;

    // -----------------------------------------------------------------------
    // KroenikFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::KroenikFile>(m, "KroenikFile", 
        R"doc(
File adapter for Kroenik (HardKloer sibling) files
The first line is the header and contains the column names:
File,  First Scan,  Last Scan,  Num of Scans,  Charge,  Monoisotopic Mass,  Base Isotope Peak,  Best Intensity,  Summed Intensity,  First RTime,  Last RTime,  Best RTime,  Best Correlation,  Modifications
Every subsequent line is a feature
All properties in the file are converted to Feature properties, whereas "First Scan", "Last Scan", "Num of Scans" and "Modifications" are stored as
metavalues with the following names "FirstScan", "LastScan", "NumOfScans" and "AveragineModifications"
The width in m/z of the overall convex hull of each feature is set to 3 Th in lack of a value provided by the Kroenik file
)doc")
        .def(nb::init<>())
        .def("load", [](OpenMS::KroenikFile& self, const OpenMS::String& filename) { OpenMS::FeatureMap feature_map; self.load(filename, feature_map); return feature_map; }, "filename"_a)
        .def("store", [](const OpenMS::KroenikFile& self, const OpenMS::String& filename, const OpenMS::FeatureMap& fmap) { self.store(filename, fmap); }, "filename"_a, "fmap"_a, "Stores a FeatureMap to a Kroenik file")
        ;

    // -----------------------------------------------------------------------
    // MRMFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MRMFile>(m, "MRMFile", "Minimal SRM/MRM file loader returning a single SwathMap wrapping the chromatogram container.")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MRMFile &>())
        .def("__copy__", [](const OpenMS::MRMFile& self) { return OpenMS::MRMFile(self); })
        .def("__deepcopy__", [](const OpenMS::MRMFile& self, nb::dict) { return OpenMS::MRMFile(self); }, "memo"_a)
        .def_static("loadMzML", [](const OpenMS::String& file, const OpenMS::String& tmp, std::shared_ptr<OpenMS::ExperimentalSettings>& exp_meta) { return OpenMS::MRMFile::loadMzML(file, tmp, exp_meta); }, "file"_a, "tmp"_a, "exp_meta"_a)
        ;

    // -----------------------------------------------------------------------
    // MSNumpressCoder
    // -----------------------------------------------------------------------
    auto msnumpresscoder_class = nb::class_<OpenMS::MSNumpressCoder>(m, "MSNumpressCoder", "Class to encode and decode data encoded with MSNumpress")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MSNumpressCoder& self) { return OpenMS::MSNumpressCoder(self); })
        .def("__deepcopy__", [](const OpenMS::MSNumpressCoder& self, nb::dict) { return OpenMS::MSNumpressCoder(self); }, "memo"_a)

        .def("encodeNP", [](OpenMS::MSNumpressCoder& self, std::vector<double> in, nb::object out_obj, bool zlib_compression, OpenMS::MSNumpressCoder::NumpressConfig config) {
            OpenMS::String result;
            self.encodeNP(in, result, zlib_compression, config);
            // Write result back to the output String-like object
            if (nb::hasattr(out_obj, "_value")) {
                out_obj.attr("_value") = std::string(result);
            }
        }, "in"_a, "result"_a, "zlib_compression"_a, "config"_a, "Encode vector of doubles to Base64 numpress string")

        .def("decodeNP", [](OpenMS::MSNumpressCoder& self, nb::object in_obj, nb::list out, bool zlib_compression, OpenMS::MSNumpressCoder::NumpressConfig config) {
            std::string in_str;
            if (nb::isinstance<nb::bytes>(in_obj)) {
                auto b = nb::cast<nb::bytes>(in_obj);
                in_str = std::string(b.c_str(), b.size());
            } else {
                in_str = nb::cast<std::string>(in_obj);
            }
            std::vector<double> result;
            self.decodeNP(in_str, result, zlib_compression, config);
            for (double v : result) out.append(v);
        }, "in"_a, "out"_a, "zlib_compression"_a, "config"_a, "Decode Base64 numpress string to vector of doubles")

        .def("encodeNPRaw", [](OpenMS::MSNumpressCoder& self, std::vector<double> in, nb::object out_obj, OpenMS::MSNumpressCoder::NumpressConfig config) {
            OpenMS::String result;
            self.encodeNPRaw(in, result, config);
            if (nb::hasattr(out_obj, "_value")) {
                // Raw encoding may contain null bytes - use the full size
                out_obj.attr("_value") = nb::bytes(result.c_str(), result.size());
            }
        }, "in"_a, "result"_a, "config"_a, "Encode vector of doubles to raw numpress byte array")

        .def("decodeNPRaw", [](OpenMS::MSNumpressCoder& self, nb::object in_obj, nb::list out, OpenMS::MSNumpressCoder::NumpressConfig config) {
            std::string in_str;
            if (nb::isinstance<nb::bytes>(in_obj)) {
                auto b = nb::cast<nb::bytes>(in_obj);
                in_str = std::string(b.c_str(), b.size());
            } else {
                in_str = nb::cast<std::string>(in_obj);
            }
            std::vector<double> result;
            self.decodeNPRaw(in_str, result, config);
            for (double v : result) out.append(v);
        }, "in"_a, "out"_a, "config"_a, "Decode raw numpress byte array to vector of doubles")
        ;
    // NumpressCompression enum nested under MSNumpressCoder
    nb::enum_<OpenMS::MSNumpressCoder::NumpressCompression>(msnumpresscoder_class, "NumpressCompression", nb::is_arithmetic())
        .value("NONE", OpenMS::MSNumpressCoder::NumpressCompression::NONE)
        .value("LINEAR", OpenMS::MSNumpressCoder::NumpressCompression::LINEAR)
        .value("PIC", OpenMS::MSNumpressCoder::NumpressCompression::PIC)
        .value("SLOF", OpenMS::MSNumpressCoder::NumpressCompression::SLOF)
        .value("SIZE_OF_NUMPRESSCOMPRESSION", OpenMS::MSNumpressCoder::NumpressCompression::SIZE_OF_NUMPRESSCOMPRESSION)
        .export_values();


    // -----------------------------------------------------------------------
    // NumpressConfig
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MSNumpressCoder::NumpressConfig>(m, "NumpressConfig", "OpenMS class NumpressConfig")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MSNumpressCoder::NumpressConfig& self) { return OpenMS::MSNumpressCoder::NumpressConfig(self); })
        .def("__deepcopy__", [](const OpenMS::MSNumpressCoder::NumpressConfig& self, nb::dict) { return OpenMS::MSNumpressCoder::NumpressConfig(self); }, "memo"_a)
        .def_rw("numpressFixedPoint", &OpenMS::MSNumpressCoder::NumpressConfig::numpressFixedPoint)
        .def_rw("numpressErrorTolerance", &OpenMS::MSNumpressCoder::NumpressConfig::numpressErrorTolerance)
        .def_rw("np_compression", &OpenMS::MSNumpressCoder::NumpressConfig::np_compression)
        .def_rw("estimate_fixed_point", &OpenMS::MSNumpressCoder::NumpressConfig::estimate_fixed_point)
        .def_rw("linear_fp_mass_acc", &OpenMS::MSNumpressCoder::NumpressConfig::linear_fp_mass_acc)
        .def("setCompression", [](OpenMS::MSNumpressCoder::NumpressConfig& self, const std::string& compression) { self.setCompression(compression); }, "compression"_a, "Sets compression from string ('none', 'linear', 'pic', 'slof')")
        ;

    // -----------------------------------------------------------------------
    // MSstatsFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MSstatsFile>(m, "MSstatsFile", "File adapter for MSstats files")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MSstatsFile& self) { return OpenMS::MSstatsFile(self); })
        .def("__deepcopy__", [](const OpenMS::MSstatsFile& self, nb::dict) { return OpenMS::MSstatsFile(self); }, "memo"_a)
        .def("storeLFQ", [](OpenMS::MSstatsFile& self, const OpenMS::String& filename, const OpenMS::ConsensusMap& consensus_map, const OpenMS::ExperimentalDesign& design, const OpenMS::StringList& reannotate_filenames, const bool is_isotope_label_type, const OpenMS::String& bioreplicate, const OpenMS::String& condition, const OpenMS::String& retention_time_summarization_method, const bool remove_shared_peptides) { self.storeLFQ(filename, consensus_map, design, reannotate_filenames, is_isotope_label_type, bioreplicate, condition, retention_time_summarization_method, remove_shared_peptides); }, "filename"_a, "consensus_map"_a, "design"_a, "reannotate_filenames"_a, "is_isotope_label_type"_a, "bioreplicate"_a, "condition"_a, "retention_time_summarization_method"_a, "remove_shared_peptides"_a = true)
        .def("storeISO", [](OpenMS::MSstatsFile& self, const OpenMS::String& filename, const OpenMS::ConsensusMap& consensus_map, const OpenMS::ExperimentalDesign& design, const OpenMS::StringList& reannotate_filenames, const OpenMS::String& bioreplicate, const OpenMS::String& condition, const OpenMS::String& mixture, const OpenMS::String& retention_time_summarization_method, const bool remove_shared_peptides) { self.storeISO(filename, consensus_map, design, reannotate_filenames, bioreplicate, condition, mixture, retention_time_summarization_method, remove_shared_peptides); }, "filename"_a, "consensus_map"_a, "design"_a, "reannotate_filenames"_a, "bioreplicate"_a, "condition"_a, "mixture"_a, "retention_time_summarization_method"_a, "remove_shared_peptides"_a = true)
        ;

    // -----------------------------------------------------------------------
    // MsInspectFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MsInspectFile>(m, "MsInspectFile", "File adapter for MsInspect files")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MsInspectFile& self) { return OpenMS::MsInspectFile(self); })
        .def("__deepcopy__", [](const OpenMS::MsInspectFile& self, nb::dict) { return OpenMS::MsInspectFile(self); }, "memo"_a)
        .def("load", [](OpenMS::MsInspectFile& self, const OpenMS::String& filename, OpenMS::FeatureMap& feature_map) { self.load(filename, feature_map); }, "filename"_a, "feature_map"_a)
        .def("store", [](const OpenMS::MsInspectFile& self, const OpenMS::String& filename, const OpenMS::MSSpectrum& spectrum) { self.store(filename, spectrum); }, "filename"_a, "spectrum"_a)
        ;

    // -----------------------------------------------------------------------
    // MzMLSpectrumDecoder
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzMLSpectrumDecoder>(m, "MzMLSpectrumDecoder", 
        R"doc(
A class to decode input strings that contain an mzML chromatogram or spectrum tag
It uses xercesc to parse a string containing either a exactly one mzML
spectrum or chromatogram (from <chromatogram> to </chromatogram> or
<spectrum> to </spectrum> tag). It returns the data contained in the
binaryDataArray for Intensity / mass-to-charge or Intensity / time
)doc")
        .def(nb::init<bool>())
        .def("__copy__", [](const OpenMS::MzMLSpectrumDecoder& self) { return OpenMS::MzMLSpectrumDecoder(self); })
        .def("__deepcopy__", [](const OpenMS::MzMLSpectrumDecoder& self, nb::dict) { return OpenMS::MzMLSpectrumDecoder(self); }, "memo"_a)
        .def("domParseSpectrum", [](OpenMS::MzMLSpectrumDecoder& self, const OpenMS::String& in, std::shared_ptr<OpenMS::Interfaces::Spectrum>& sptr) { return self.domParseSpectrum(in, sptr); }, "in"_a, "sptr"_a, 
            R"doc(
Extract data from a string which contains a full mzML chromatogram
Extracts data from the input string which is expected to contain exactly
one <chromatogram> tag (from <chromatogram> to </chromatogram>). This
function will extract the contained binaryDataArray and provide the
result as Chromatogram
:param in: Input string containing the raw XML
:param cptr: Resulting chromatogram
)doc")
        .def("domParseSpectrum", [](OpenMS::MzMLSpectrumDecoder& self, const OpenMS::String& in, OpenMS::MSSpectrum& s) { return self.domParseSpectrum(in, s); }, "in"_a, "s"_a, 
            R"doc(
Extract data from a string which contains a full mzML chromatogram
Extracts data from the input string which is expected to contain exactly
one <chromatogram> tag (from <chromatogram> to </chromatogram>). This
function will extract the contained binaryDataArray and provide the
result as Chromatogram
:param in: Input string containing the raw XML
:param cptr: Resulting chromatogram
)doc")
        .def("domParseChromatogram", [](OpenMS::MzMLSpectrumDecoder& self, const OpenMS::String& in, OpenMS::MSChromatogram& c) { return self.domParseChromatogram(in, c); }, "in"_a, "c"_a)
        .def("domParseChromatogram", [](OpenMS::MzMLSpectrumDecoder& self, const OpenMS::String& in, std::shared_ptr<OpenMS::Interfaces::Chromatogram>& cptr) { return self.domParseChromatogram(in, cptr); }, "in"_a, "cptr"_a)
        .def("setSkipXMLChecks", [](OpenMS::MzMLSpectrumDecoder& self, bool only) { return self.setSkipXMLChecks(only); }, "only"_a, "Whether to skip some XML checks (e.g. removing whitespace inside base64 arrays) and be fast instead")
        ;

    // -----------------------------------------------------------------------
    // MzMLSwathFileConsumer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzMLSwathFileConsumer>(m, "MzMLSwathFileConsumer", "FullSwathFileConsumer")
        .def(nb::init<OpenMS::String, OpenMS::String, size_t, std::vector<int>>())
        .def(nb::init<std::vector<OpenSwath::SwathMap>, OpenMS::String, OpenMS::String, size_t, std::vector<int>>())
        .def("setExpectedSize", [](OpenMS::MzMLSwathFileConsumer& self, size_t p0, size_t p1) { return self.setExpectedSize(p0, p1); })
        .def("setExperimentalSettings", [](OpenMS::MzMLSwathFileConsumer& self, const OpenMS::ExperimentalSettings& exp) { return self.setExperimentalSettings(exp); }, "exp"_a)
        .def("retrieveSwathMaps", [](OpenMS::MzMLSwathFileConsumer& self) { std::vector<OpenSwath::SwathMap> maps; self.retrieveSwathMaps(maps); return maps; })
        .def("consumeChromatogram", [](OpenMS::MzMLSwathFileConsumer& self, OpenMS::MSChromatogram& p0) { return self.consumeChromatogram(p0); })
        .def("consumeSpectrum", [](OpenMS::MzMLSwathFileConsumer& self, OpenMS::MSSpectrum& s) { return self.consumeSpectrum(s); }, "s"_a)
        ;

    // -----------------------------------------------------------------------
    // MzQCFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzQCFile>(m, "MzQCFile", 
        R"doc(
File adapter for mzQC files used to load and store mzQC files
This class collects the data for the mzQC File
)doc")
        .def(nb::init<>())
        .def("store", [](const OpenMS::MzQCFile& self, const OpenMS::String& input_file, const OpenMS::String& output_file, const OpenMS::MSExperiment& exp, const OpenMS::String& contact_name, const OpenMS::String& contact_address, const OpenMS::String& description, const OpenMS::String& label, const OpenMS::FeatureMap& feature_map, std::vector<OpenMS::ProteinIdentification>& prot_ids, OpenMS::PeptideIdentificationList& pep_ids) { self.store(input_file, output_file, exp, contact_name, contact_address, description, label, feature_map, prot_ids, pep_ids); }, "input_file"_a, "output_file"_a, "exp"_a, "contact_name"_a, "contact_address"_a, "description"_a, "label"_a, "feature_map"_a, "prot_ids"_a, "pep_ids"_a)
        ;

    // -----------------------------------------------------------------------
    // MzTab
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzTab>(m, "MzTab", 
        R"doc(
Data model of MzTab files
Please see the official MzTab specification at https://code.google.com/p/mztab/
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MzTab& self) { return OpenMS::MzTab(self); })
        .def("__deepcopy__", [](const OpenMS::MzTab& self, nb::dict) { return OpenMS::MzTab(self); }, "memo"_a)
        ;

    // -----------------------------------------------------------------------
    // MzTabFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzTabFile>(m, "MzTabFile", "File adapter for MzTab files")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MzTabFile& self) { return OpenMS::MzTabFile(self); })
        .def("__deepcopy__", [](const OpenMS::MzTabFile& self, nb::dict) { return OpenMS::MzTabFile(self); }, "memo"_a)
        .def("store", [](const OpenMS::MzTabFile& self, const OpenMS::String& filename, const OpenMS::MzTab& mz_tab) { return self.store(filename, mz_tab); }, "filename"_a, "mz_tab"_a, "Stores MzTab file")
        .def("store", [](OpenMS::MzTabFile& self, const OpenMS::String& filename, const std::vector<OpenMS::ProteinIdentification>& protein_identifications, const OpenMS::PeptideIdentificationList& peptide_identifications, bool first_run_inference_only, bool export_empty_pep_ids, bool export_all_psms, const OpenMS::String& title) { return self.store(filename, protein_identifications, peptide_identifications, first_run_inference_only, export_empty_pep_ids, export_all_psms, title); }, "filename"_a, "protein_identifications"_a, "peptide_identifications"_a, "first_run_inference_only"_a, "export_empty_pep_ids"_a = false, "export_all_psms"_a = false, "title"_a = "ID export from OpenMS", "Stores MzTab file")
        .def("load", [](OpenMS::MzTabFile& self, const OpenMS::String& filename) { OpenMS::MzTab mz_tab; self.load(filename, mz_tab); return mz_tab; }, "filename"_a, "Loads MzTab file")
        ;

    // -----------------------------------------------------------------------
    // MzTabM
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzTabM>(m, "MzTabM", 
        R"doc(
Data model of MzTabM files
Please see the official MzTabM specification at https://github.com/HUPO-PSI/mzTab/tree/master/specification_document-releases/2_0-Metabolomics-Release
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MzTabM& self) { return OpenMS::MzTabM(self); })
        .def("__deepcopy__", [](const OpenMS::MzTabM& self, nb::dict) { return OpenMS::MzTabM(self); }, "memo"_a)
        .def_static("exportFeatureMapToMzTabM", [](const OpenMS::FeatureMap& feature_map) { return OpenMS::MzTabM::exportFeatureMapToMzTabM(feature_map); }, "feature_map"_a, "Export FeatureMap with Identifications to MzTabM")
        ;

    // -----------------------------------------------------------------------
    // MzTabMFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::MzTabMFile>(m, "MzTabMFile", "File adapter for MzTab-M files")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::MzTabMFile& self) { return OpenMS::MzTabMFile(self); })
        .def("__deepcopy__", [](const OpenMS::MzTabMFile& self, nb::dict) { return OpenMS::MzTabMFile(self); }, "memo"_a)
        .def("store", [](const OpenMS::MzTabMFile& self, const OpenMS::String& filename, const OpenMS::MzTabM& mztab_m) { return self.store(filename, mztab_m); }, "filename"_a, "mztab_m"_a, "Store MzTabM file")
        ;

    // -----------------------------------------------------------------------
    // NoopMSDataWritingConsumer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::NoopMSDataWritingConsumer>(m, "NoopMSDataWritingConsumer", 
        R"doc(
Consumer class that perform no operation
This is sometimes necessary to fulfill the requirement of passing an
valid MSDataWritingConsumer object or pointer but no operation is
required
)doc")
        .def(nb::init<OpenMS::String>())
        .def("setExperimentalSettings", [](OpenMS::NoopMSDataWritingConsumer& self, const OpenMS::ExperimentalSettings& p0) { return self.setExperimentalSettings(p0); })
        .def("consumeSpectrum", [](OpenMS::NoopMSDataWritingConsumer& self, OpenMS::MSSpectrum& p0) { return self.consumeSpectrum(p0); })
        .def("consumeChromatogram", [](OpenMS::NoopMSDataWritingConsumer& self, OpenMS::MSChromatogram& p0) { return self.consumeChromatogram(p0); })
        .def("setExpectedSize", [](OpenMS::NoopMSDataWritingConsumer& self, size_t expectedSpectra, size_t expectedChromatograms) { return self.setExpectedSize(expectedSpectra, expectedChromatograms); }, "expectedSpectra"_a, "expectedChromatograms"_a)
        .def("addDataProcessing", [](OpenMS::NoopMSDataWritingConsumer& self, OpenMS::DataProcessing d) { return self.addDataProcessing(d); }, "d"_a)
        .def("getNrSpectraWritten", [](OpenMS::NoopMSDataWritingConsumer& self) { return self.getNrSpectraWritten(); })
        .def("getNrChromatogramsWritten", [](OpenMS::NoopMSDataWritingConsumer& self) { return self.getNrChromatogramsWritten(); })
        ;

    // -----------------------------------------------------------------------
    // OMSSACSVFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::OMSSACSVFile>(m, "OMSSACSVFile",
        R"doc(
File adapter for OMSSACSV files
The files contain the results of the OMSSA algorithm in a comma separated manner. This file adapter is able to
load the data from such a file into the structures of OpenMS
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::OMSSACSVFile& self) { return OpenMS::OMSSACSVFile(self); })
        .def("__deepcopy__", [](const OpenMS::OMSSACSVFile& self, nb::dict) { return OpenMS::OMSSACSVFile(self); }, "memo"_a)
        .def("load", [](const OpenMS::OMSSACSVFile& self, const OpenMS::String& filename) { OpenMS::ProteinIdentification protein_identification; OpenMS::PeptideIdentificationList id_data; self.load(filename, protein_identification, id_data); return std::make_tuple(protein_identification, id_data); }, "filename"_a)
        ;

    // -----------------------------------------------------------------------
    // PEFFCustomKeyDef
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PEFFCustomKeyDef>(m, "PEFFCustomKeyDef", 
        R"doc(
Represents a custom key definition from the PEFF header
Attributes:
key_name: Custom key name
description: Description of the key
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PEFFCustomKeyDef &>())
        .def("__copy__", [](const OpenMS::PEFFCustomKeyDef& self) { return OpenMS::PEFFCustomKeyDef(self); })
        .def("__deepcopy__", [](const OpenMS::PEFFCustomKeyDef& self, nb::dict) { return OpenMS::PEFFCustomKeyDef(self); }, "memo"_a)
        .def_rw("key_name", &OpenMS::PEFFCustomKeyDef::key_name)
        .def_rw("description", &OpenMS::PEFFCustomKeyDef::description)
        .def_rw("concept_curie", &OpenMS::PEFFCustomKeyDef::concept_curie)
        .def_rw("regexp", &OpenMS::PEFFCustomKeyDef::regexp)
        .def_rw("field_names", &OpenMS::PEFFCustomKeyDef::field_names)
        .def_rw("field_types", &OpenMS::PEFFCustomKeyDef::field_types)
        ;

    // -----------------------------------------------------------------------
    // PEFFDatabaseMetadata
    // -----------------------------------------------------------------------
    auto peffdatabasemetadata_class = nb::class_<OpenMS::PEFFDatabaseMetadata>(m, "PEFFDatabaseMetadata", 
        R"doc(
Metadata from a PEFF database header section
Attributes:
version: PEFF version
db_name: Database name
prefix: Database prefix
db_description: Database description
is_decoy: Whether this is a decoy database
db_sources: List of database sources
db_version: Database version
number_of_entries: Number of entries
sequence_type: Sequence type (AA or NA)
general_comments: List of general comments
conversion: Conversion notes
specific_keys: Custom key definitions
specific_values: Custom value type definitions
has_annotation_identifiers: Whether entries use annotation identifiers
is_proteoform_db: Whether this is a proteoform database
unrecognized_keys: Unrecognized header keys preserved for round-trip
)doc")
        .def(nb::init<>())
        .def_rw("version", &OpenMS::PEFFDatabaseMetadata::version)
        .def_rw("db_name", &OpenMS::PEFFDatabaseMetadata::db_name)
        .def_rw("prefix", &OpenMS::PEFFDatabaseMetadata::prefix)
        .def_rw("db_description", &OpenMS::PEFFDatabaseMetadata::db_description)
        .def_rw("is_decoy", &OpenMS::PEFFDatabaseMetadata::is_decoy)
        .def_rw("db_sources", &OpenMS::PEFFDatabaseMetadata::db_sources)
        .def_rw("db_version", &OpenMS::PEFFDatabaseMetadata::db_version)
        .def_rw("number_of_entries", &OpenMS::PEFFDatabaseMetadata::number_of_entries)
        .def_rw("sequence_type", &OpenMS::PEFFDatabaseMetadata::sequence_type)
        .def_rw("general_comments", &OpenMS::PEFFDatabaseMetadata::general_comments)
        .def_rw("conversion", &OpenMS::PEFFDatabaseMetadata::conversion)
        .def_rw("specific_keys", &OpenMS::PEFFDatabaseMetadata::specific_keys)
        .def_rw("specific_values", &OpenMS::PEFFDatabaseMetadata::specific_values)
        .def_rw("optional_tag_defs", &OpenMS::PEFFDatabaseMetadata::optional_tag_defs)
        .def_rw("custom_key_defs", &OpenMS::PEFFDatabaseMetadata::custom_key_defs)
        .def_rw("has_annotation_identifiers", &OpenMS::PEFFDatabaseMetadata::has_annotation_identifiers)
        .def_rw("is_proteoform_db", &OpenMS::PEFFDatabaseMetadata::is_proteoform_db)
        .def_rw("unrecognized_keys", &OpenMS::PEFFDatabaseMetadata::unrecognized_keys)
        ;
    // SequenceType enum nested under PEFFDatabaseMetadata
    nb::enum_<OpenMS::PEFFDatabaseMetadata::SequenceType>(peffdatabasemetadata_class, "SequenceType", nb::is_arithmetic())
        .value("AA", OpenMS::PEFFDatabaseMetadata::SequenceType::AA)
        .value("NA", OpenMS::PEFFDatabaseMetadata::SequenceType::NA)
        ;

    // -----------------------------------------------------------------------
    // PEFFDisulfideBond
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PEFFDisulfideBond>(m, "PEFFDisulfideBond", 
        R"doc(
Represents a disulfide bond annotation in PEFF
Attributes:
id1: First cysteine reference (AnnotationIdentifier)
id2: Second cysteine reference (AnnotationIdentifier)
optional_tag: Optional tag (e.g., "between chains")
annotation_id: Optional annotation identifier (UInt, max value = not set)
)doc")
        .def(nb::init<>())
        .def(nb::init<OpenMS::String, OpenMS::String, OpenMS::String, unsigned int>())
        .def_rw("id1", &OpenMS::PEFFDisulfideBond::id1)
        .def_rw("id2", &OpenMS::PEFFDisulfideBond::id2)
        .def_rw("optional_tag", &OpenMS::PEFFDisulfideBond::optional_tag)
        .def_rw("annotation_id", &OpenMS::PEFFDisulfideBond::annotation_id)
        ;

    // -----------------------------------------------------------------------
    // PEFFEntry
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PEFFEntry>(m, "PEFFEntry", 
        R"doc(
Represents a single entry in a PEFF file with all annotations.
Each entry corresponds to one description line and sequence.
The description line format is: >Prefix:DbUniqueId \\key=value ...
Attributes:
prefix: Database prefix from description line
identifier: Protein identifier (Prefix:DbUniqueId)
sequence: Amino acid sequence
protein_names: List of protein names
gene_name: Gene name
ncbi_tax_id: NCBI taxonomy ID
taxonomy_name: Taxonomy name
sequence_length: Sequence length
modifications: List of modifications
simple_variants: List of simple variants
complex_variants: List of complex variants
processed_regions: List of processed regions
proteoforms: List of proteoforms in ProForma notation
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PEFFEntry &>())
        .def("__copy__", [](const OpenMS::PEFFEntry& self) { return OpenMS::PEFFEntry(self); })
        .def("__deepcopy__", [](const OpenMS::PEFFEntry& self, nb::dict) { return OpenMS::PEFFEntry(self); }, "memo"_a)
        .def("toFASTAEntry", [](const OpenMS::PEFFEntry& self) { return self.toFASTAEntry(); })
        .def_static("fromFASTAEntry", [](const OpenMS::FASTAFile::FASTAEntry& fasta) { return OpenMS::PEFFEntry::fromFASTAEntry(fasta); }, "fasta"_a)
        .def("getSequence", [](const OpenMS::PEFFEntry& self) { return self.getSequence(); }, "Convert to a FASTAEntry (loses PEFF-specific annotations)")
        .def("getModifiedSequence", [](const OpenMS::PEFFEntry& self) { return self.getModifiedSequence(); }, "Get the base AASequence for this entry (unmodified sequence)")
        .def("getVariantSequences", [](const OpenMS::PEFFEntry& self, bool include_complex) { std::vector<std::basic_string<char>> descriptions; std::vector<OpenMS::AASequence> sequences; self.getVariantSequences(descriptions, sequences, include_complex); return std::make_tuple(descriptions, sequences); }, "include_complex"_a, 
            R"doc(
Get an AASequence with all annotated modifications applied.
Modifications with unknown positions (position == 0) are skipped.
)doc")
        .def("getProcessedSequence", [](const OpenMS::PEFFEntry& self, const OpenMS::String& region_accession) { return self.getProcessedSequence(region_accession); }, "region_accession"_a, 
            R"doc(
Get processed sequence (e.g., mature protein without signal peptide).
Uses default region accession "PEFF:0001021" (signal peptide).
)doc")
        .def_rw("prefix", &OpenMS::PEFFEntry::prefix)
        .def_rw("identifier", &OpenMS::PEFFEntry::identifier)
        .def_rw("sequence", &OpenMS::PEFFEntry::sequence)
        .def_rw("protein_names", &OpenMS::PEFFEntry::protein_names)
        .def_rw("gene_name", &OpenMS::PEFFEntry::gene_name)
        .def_rw("ncbi_tax_id", &OpenMS::PEFFEntry::ncbi_tax_id)
        .def_rw("taxonomy_name", &OpenMS::PEFFEntry::taxonomy_name)
        .def_rw("sequence_length", &OpenMS::PEFFEntry::sequence_length)
        .def_rw("sequence_version", &OpenMS::PEFFEntry::sequence_version)
        .def_rw("entry_version", &OpenMS::PEFFEntry::entry_version)
        .def_rw("protein_existence", &OpenMS::PEFFEntry::protein_existence)
        .def_rw("db_unique_id", &OpenMS::PEFFEntry::db_unique_id)
        .def_rw("entry_id", &OpenMS::PEFFEntry::entry_id)
        .def_rw("alt_accessions", &OpenMS::PEFFEntry::alt_accessions)
        .def_rw("modifications", &OpenMS::PEFFEntry::modifications)
        .def_rw("simple_variants", &OpenMS::PEFFEntry::simple_variants)
        .def_rw("complex_variants", &OpenMS::PEFFEntry::complex_variants)
        .def_rw("processed_regions", &OpenMS::PEFFEntry::processed_regions)
        .def_rw("disulfide_bonds", &OpenMS::PEFFEntry::disulfide_bonds)
        .def_rw("proteoforms", &OpenMS::PEFFEntry::proteoforms)
        .def_rw("custom_annotations", &OpenMS::PEFFEntry::custom_annotations)
        .def("digestWithVariants", [](const OpenMS::PEFFEntry& self, const OpenMS::ProteaseDigestion& digestor, OpenMS::Size min_length, OpenMS::Size max_length, bool include_reference, bool include_variants, bool include_modifications) {
            std::vector<std::string> descriptions;
            std::vector<OpenMS::AASequence> sequences;
            self.digestWithVariants(digestor, descriptions, sequences, min_length, max_length, include_reference, include_variants, include_modifications);
            return nb::make_tuple(descriptions, sequences);
        }, "digestor"_a, "min_length"_a = 6, "max_length"_a = 40, "include_reference"_a = true, "include_variants"_a = true, "include_modifications"_a = false, "Digest with variants. Returns (descriptions, sequences)")
        .def("generatePeptides", [](const OpenMS::PEFFEntry& self, const OpenMS::ProteaseDigestion& digestor, const std::vector<std::string>& fixed_mods, const std::vector<std::string>& variable_mods, OpenMS::Size max_variable_mods_per_peptide, OpenMS::Size min_length, OpenMS::Size max_length, bool include_reference, bool include_peff_variants, bool include_peff_modifications) {
            std::vector<std::string> descriptions;
            std::vector<OpenMS::AASequence> sequences;
            self.generatePeptides(digestor, descriptions, sequences, fixed_mods, variable_mods, max_variable_mods_per_peptide, min_length, max_length, include_reference, include_peff_variants, include_peff_modifications);
            return nb::make_tuple(descriptions, sequences);
        }, "digestor"_a, "fixed_mods"_a = std::vector<std::string>{}, "variable_mods"_a = std::vector<std::string>{}, "max_variable_mods_per_peptide"_a = 2, "min_length"_a = 6, "max_length"_a = 40, "include_reference"_a = true, "include_peff_variants"_a = true, "include_peff_modifications"_a = true, "Generate peptides with PEFF annotations. Returns (descriptions, sequences)")
        ;

    // -----------------------------------------------------------------------
    // PEFFModification
    // -----------------------------------------------------------------------
    auto peffmodification_class = nb::class_<OpenMS::PEFFModification>(m, "PEFFModification", 
        R"doc(
Represents a PEFF modification annotation
Attributes:
position: 1-based position (0 = unknown)
accession: Modification accession (MOD:xxxxx, UNIMOD:xx, or custom)
name: Human-readable name
optional_tag: Optional tag (last component of annotation tuple)
annotation_id: Optional annotation identifier (UInt, max value = not set)
type: Modification type (PSI_MOD, UNIMOD, or GENERIC)
)doc")
        .def(nb::init<>())
        .def(nb::init<size_t, OpenMS::String, OpenMS::String, OpenMS::String, unsigned int>())
        .def_rw("position", &OpenMS::PEFFModification::position)
        .def_rw("accession", &OpenMS::PEFFModification::accession)
        .def_rw("name", &OpenMS::PEFFModification::name)
        .def_rw("optional_tag", &OpenMS::PEFFModification::optional_tag)
        .def_rw("annotation_id", &OpenMS::PEFFModification::annotation_id)
        .def_rw("type", &OpenMS::PEFFModification::type)
        ;
    // ModificationType enum nested under PEFFModification
    nb::enum_<OpenMS::PEFFModification::Type>(peffmodification_class, "ModificationType", nb::is_arithmetic())
        .value("PSI_MOD", OpenMS::PEFFModification::Type::PSI_MOD)
        .value("UNIMOD", OpenMS::PEFFModification::Type::UNIMOD)
        .value("GENERIC", OpenMS::PEFFModification::Type::GENERIC)
        ;

    // -----------------------------------------------------------------------
    // PEFFProcessedRegion
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PEFFProcessedRegion>(m, "PEFFProcessedRegion", 
        R"doc(
Represents a PEFF processed region (signal peptide, etc.)
Attributes:
start_position: 1-based start position
end_position: 1-based end position
accession: PEFF CV accession (e.g., PEFF:0001021)
name: Optional name
optional_tag: Optional tag (last component of annotation tuple)
annotation_id: Optional annotation identifier (UInt, max value = not set)
)doc")
        .def(nb::init<>())
        .def(nb::init<size_t, size_t, OpenMS::String, OpenMS::String, OpenMS::String, unsigned int>())
        .def_rw("start_position", &OpenMS::PEFFProcessedRegion::start_position)
        .def_rw("end_position", &OpenMS::PEFFProcessedRegion::end_position)
        .def_rw("accession", &OpenMS::PEFFProcessedRegion::accession)
        .def_rw("name", &OpenMS::PEFFProcessedRegion::name)
        .def_rw("optional_tag", &OpenMS::PEFFProcessedRegion::optional_tag)
        .def_rw("annotation_id", &OpenMS::PEFFProcessedRegion::annotation_id)
        ;

    // -----------------------------------------------------------------------
    // PEFFVariantComplex
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PEFFVariantComplex>(m, "PEFFVariantComplex", 
        R"doc(
Represents a complex PEFF variant (multi-residue change)
Attributes:
start_position: 1-based start position
end_position: 1-based end position
replacement: Replacement sequence (empty = deletion)
optional_tag: Optional tag (last component of annotation tuple)
annotation_id: Optional annotation identifier (UInt, max value = not set)
)doc")
        .def(nb::init<>())
        .def(nb::init<size_t, size_t, OpenMS::String, OpenMS::String, unsigned int>())
        .def_rw("start_position", &OpenMS::PEFFVariantComplex::start_position)
        .def_rw("end_position", &OpenMS::PEFFVariantComplex::end_position)
        .def_rw("replacement", &OpenMS::PEFFVariantComplex::replacement)
        .def_rw("optional_tag", &OpenMS::PEFFVariantComplex::optional_tag)
        .def_rw("annotation_id", &OpenMS::PEFFVariantComplex::annotation_id)
        ;

    // -----------------------------------------------------------------------
    // PEFFVariantSimple
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PEFFVariantSimple>(m, "PEFFVariantSimple", 
        R"doc(
Represents a simple PEFF variant (single amino acid substitution)
Attributes:
position: 1-based position
variant_aa: Variant amino acid character
optional_tag: Optional tag (last component of annotation tuple)
annotation_id: Optional annotation identifier (UInt, max value = not set)
)doc")
        .def(nb::init<>())
        .def(nb::init<size_t, char, OpenMS::String, unsigned int>())
        .def_rw("position", &OpenMS::PEFFVariantSimple::position)
        .def_rw("variant_aa", &OpenMS::PEFFVariantSimple::variant_aa)
        .def_rw("optional_tag", &OpenMS::PEFFVariantSimple::optional_tag)
        .def_rw("annotation_id", &OpenMS::PEFFVariantSimple::annotation_id)
        ;

    // -----------------------------------------------------------------------
    // ParamCTDFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ParamCTDFile>(m, "ParamCTDFile", "Serializes a Param object to/from a CTD (Common Tool Description) file")
        .def(nb::init<>())
        .def("store", [](const OpenMS::ParamCTDFile& self, const OpenMS::String& filename, const OpenMS::Param& param, const OpenMS::ToolInfo& tool_info) { return self.store(filename, param, tool_info); }, "filename"_a, "param"_a, "tool_info"_a)
        ;

    // -----------------------------------------------------------------------
    // ParquetFilter
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ParquetFilter>(m, "ParquetFilter", "OpenMS class ParquetFilter")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ParquetFilter &>())
        .def("__copy__", [](const OpenMS::ParquetFilter& self) { return OpenMS::ParquetFilter(self); })
        .def("__deepcopy__", [](const OpenMS::ParquetFilter& self, nb::dict) { return OpenMS::ParquetFilter(self); }, "memo"_a)
        .def("andNext", [](OpenMS::ParquetFilter& self) -> OpenMS::ParquetFilter & { return self.andNext(); }, nb::rv_policy::reference_internal)
        .def("orNext", [](OpenMS::ParquetFilter& self) -> OpenMS::ParquetFilter & { return self.orNext(); }, nb::rv_policy::reference_internal)
        .def("eq", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilter & { return self.eq(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("ne", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilter & { return self.ne(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("lt", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilter & { return self.lt(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("le", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilter & { return self.le(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("gt", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilter & { return self.gt(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("ge", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilter & { return self.ge(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("eq", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilter & { return self.eq(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("ne", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilter & { return self.ne(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("lt", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilter & { return self.lt(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("le", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilter & { return self.le(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("gt", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilter & { return self.gt(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("ge", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilter & { return self.ge(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("empty", [](const OpenMS::ParquetFilter& self) { return self.empty(); })
        .def("in_", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, const std::vector<OpenMS::Int64>& values) -> OpenMS::ParquetFilter & { return self.in(column, values); }, "column"_a, "values"_a, nb::rv_policy::reference_internal)
        .def("in_", [](OpenMS::ParquetFilter& self, const OpenMS::String& column, const std::vector<OpenMS::String>& values) -> OpenMS::ParquetFilter & { return self.in(column, values); }, "column"_a, "values"_a, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // ParquetFilterBuilder
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ParquetFilterBuilder>(m, "ParquetFilterBuilder", "OpenMS class ParquetFilterBuilder")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ParquetFilterBuilder &>())
        .def("__copy__", [](const OpenMS::ParquetFilterBuilder& self) { return OpenMS::ParquetFilterBuilder(self); })
        .def("__deepcopy__", [](const OpenMS::ParquetFilterBuilder& self, nb::dict) { return OpenMS::ParquetFilterBuilder(self); }, "memo"_a)
        .def("andNext", [](OpenMS::ParquetFilterBuilder& self) -> OpenMS::ParquetFilterBuilder & { return self.andNext(); }, nb::rv_policy::reference_internal)
        .def("orNext", [](OpenMS::ParquetFilterBuilder& self) -> OpenMS::ParquetFilterBuilder & { return self.orNext(); }, nb::rv_policy::reference_internal)
        .def("eq", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilterBuilder & { return self.eq(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("ne", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilterBuilder & { return self.ne(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("lt", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilterBuilder & { return self.lt(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("le", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilterBuilder & { return self.le(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("gt", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilterBuilder & { return self.gt(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("ge", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, long value) -> OpenMS::ParquetFilterBuilder & { return self.ge(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("eq", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilterBuilder & { return self.eq(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("ne", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilterBuilder & { return self.ne(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("lt", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilterBuilder & { return self.lt(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("le", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilterBuilder & { return self.le(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("gt", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilterBuilder & { return self.gt(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("ge", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, const OpenMS::String& value) -> OpenMS::ParquetFilterBuilder & { return self.ge(column, value); }, "column"_a, "value"_a, nb::rv_policy::reference_internal)
        .def("filter", [](const OpenMS::ParquetFilterBuilder& self) -> const OpenMS::ParquetFilter & { return self.filter(); }, nb::rv_policy::reference_internal)
        .def("empty", [](const OpenMS::ParquetFilterBuilder& self) { return self.empty(); })
        .def("in_", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, const std::vector<OpenMS::Int64>& values) -> OpenMS::ParquetFilterBuilder & { return self.in(column, values); }, "column"_a, "values"_a, nb::rv_policy::reference_internal)
        .def("in_", [](OpenMS::ParquetFilterBuilder& self, const OpenMS::String& column, const std::vector<OpenMS::String>& values) -> OpenMS::ParquetFilterBuilder & { return self.in(column, values); }, "column"_a, "values"_a, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // PeakFileOptions
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakFileOptions>(m, "PeakFileOptions", "Options for loading files containing peak data")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PeakFileOptions &>())
        .def("__copy__", [](const OpenMS::PeakFileOptions& self) { return OpenMS::PeakFileOptions(self); })
        .def("__deepcopy__", [](const OpenMS::PeakFileOptions& self, nb::dict) { return OpenMS::PeakFileOptions(self); }, "memo"_a)
        .def("setMetadataOnly", [](OpenMS::PeakFileOptions& self, bool only) { return self.setMetadataOnly(only); }, "only"_a, "Sets whether or not to load only meta data")
        .def("getMetadataOnly", [](const OpenMS::PeakFileOptions& self) { return self.getMetadataOnly(); }, "Returns whether or not to load only meta data")
        .def("setForceMQCompatability", [](OpenMS::PeakFileOptions& self, bool forceMQ) { return self.setForceMQCompatability(forceMQ); }, "forceMQ"_a, "[mzXML only!]Returns Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)")
        .def("getForceMQCompatability", [](const OpenMS::PeakFileOptions& self) { return self.getForceMQCompatability(); }, "[mzXML only!]Returns Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)")
        .def("setForceTPPCompatability", [](OpenMS::PeakFileOptions& self, bool forceTPP) { return self.setForceTPPCompatability(forceTPP); }, "forceTPP"_a, 
            R"doc(
[ mzML only!]Returns Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
)doc")
        .def("getForceTPPCompatability", [](const OpenMS::PeakFileOptions& self) { return self.getForceTPPCompatability(); }, 
            R"doc(
[mzML only!]Returns Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
)doc")
        .def("setWriteSupplementalData", [](OpenMS::PeakFileOptions& self, bool write) { return self.setWriteSupplementalData(write); }, "write"_a, "Sets whether or not to write supplemental peak data in MzData files")
        .def("getWriteSupplementalData", [](const OpenMS::PeakFileOptions& self) { return self.getWriteSupplementalData(); }, "Returns whether or not to write supplemental peak data in MzData files")
        .def("setRTRange", [](OpenMS::PeakFileOptions& self, const OpenMS::DRange<1>& range) { return self.setRTRange(range); }, "range"_a, "Restricts the range of RT values for peaks to load")
        .def("hasRTRange", [](const OpenMS::PeakFileOptions& self) { return self.hasRTRange(); }, "Returns true if an RT range has been set")
        .def("getRTRange", [](const OpenMS::PeakFileOptions& self) -> const OpenMS::DRange<1> & { return self.getRTRange(); }, nb::rv_policy::reference_internal, "Returns the RT range")
        .def("setMZRange", [](OpenMS::PeakFileOptions& self, const OpenMS::DRange<1>& range) { return self.setMZRange(range); }, "range"_a, "Restricts the range of MZ values for peaks to load")
        .def("hasMZRange", [](const OpenMS::PeakFileOptions& self) { return self.hasMZRange(); }, "Returns true if an MZ range has been set")
        .def("getMZRange", [](const OpenMS::PeakFileOptions& self) -> const OpenMS::DRange<1> & { return self.getMZRange(); }, nb::rv_policy::reference_internal, "Returns the MZ range")
        .def("setIntensityRange", [](OpenMS::PeakFileOptions& self, const OpenMS::DRange<1>& range) { return self.setIntensityRange(range); }, "range"_a, "Restricts the range of intensity values for peaks to load")
        .def("hasIntensityRange", [](const OpenMS::PeakFileOptions& self) { return self.hasIntensityRange(); }, "Returns true if an intensity range has been set")
        .def("getIntensityRange", [](const OpenMS::PeakFileOptions& self) -> const OpenMS::DRange<1> & { return self.getIntensityRange(); }, nb::rv_policy::reference_internal, "Returns the intensity range")
        .def("setPrecursorMZRange", [](OpenMS::PeakFileOptions& self, const OpenMS::DRange<1>& range) { return self.setPrecursorMZRange(range); }, "range"_a, "Restricts the range of precursor m/z values for MS2+ spectra to load")
        .def("hasPrecursorMZRange", [](const OpenMS::PeakFileOptions& self) { return self.hasPrecursorMZRange(); }, "Returns true if a precursor m/z range has been set")
        .def("getPrecursorMZRange", [](const OpenMS::PeakFileOptions& self) -> const OpenMS::DRange<1> & { return self.getPrecursorMZRange(); }, nb::rv_policy::reference_internal, "Returns the precursor m/z range")
        .def("setMSLevels", [](OpenMS::PeakFileOptions& self, const std::vector<int>& levels) { return self.setMSLevels(levels); }, "levels"_a, "Sets the desired MS levels for peaks to load")
        .def("addMSLevel", [](OpenMS::PeakFileOptions& self, int level) { return self.addMSLevel(level); }, "level"_a, "Adds a desired MS level for peaks to load")
        .def("clearMSLevels", [](OpenMS::PeakFileOptions& self) { return self.clearMSLevels(); }, "Clears the MS levels")
        .def("hasMSLevels", [](const OpenMS::PeakFileOptions& self) { return self.hasMSLevels(); }, "Returns true, if MS levels have been set")
        .def("containsMSLevel", [](const OpenMS::PeakFileOptions& self, int level) { return self.containsMSLevel(level); }, "level"_a, "Returns true, if MS level `level` has been set")
        .def("getMSLevels", [](const OpenMS::PeakFileOptions& self) -> const std::vector<int> & { return self.getMSLevels(); }, nb::rv_policy::reference_internal, "Returns the set MS levels")
        .def("setCompression", [](OpenMS::PeakFileOptions& self, bool compress) { return self.setCompression(compress); }, "compress"_a, "Sets if data should be compressed when writing")
        .def("getCompression", [](const OpenMS::PeakFileOptions& self) { return self.getCompression(); }, "Returns true, if data should be compressed when writing")
        .def("setFillData", [](OpenMS::PeakFileOptions& self, bool only) { return self.setFillData(only); }, "only"_a, "Sets whether to fill the actual data into the container (spectrum/chromatogram)")
        .def("getFillData", [](const OpenMS::PeakFileOptions& self) { return self.getFillData(); }, "Returns whether to fill the actual data into the container (spectrum/chromatogram)")
        .def("setSkipXMLChecks", [](OpenMS::PeakFileOptions& self, bool only) { return self.setSkipXMLChecks(only); }, "only"_a, "Sets whether to skip some XML checks and be fast instead")
        .def("getSkipXMLChecks", [](const OpenMS::PeakFileOptions& self) { return self.getSkipXMLChecks(); }, "Returns whether to skip some XML checks and be fast instead")
        .def("setSortSpectraByMZ", [](OpenMS::PeakFileOptions& self, bool sort) { return self.setSortSpectraByMZ(sort); }, "sort"_a, "Sets whether or not to sort peaks in spectra")
        .def("getSortSpectraByMZ", [](const OpenMS::PeakFileOptions& self) { return self.getSortSpectraByMZ(); }, "Returns whether or not peaks in spectra should be sorted")
        .def("setSortChromatogramsByRT", [](OpenMS::PeakFileOptions& self, bool sort) { return self.setSortChromatogramsByRT(sort); }, "sort"_a, "Sets whether or not to sort peaks in chromatograms")
        .def("getSortChromatogramsByRT", [](const OpenMS::PeakFileOptions& self) { return self.getSortChromatogramsByRT(); }, "Returns whether or not peaks in chromatograms should be sorted")
        .def("setMz32Bit", [](OpenMS::PeakFileOptions& self, bool mz_32_bit) { return self.setMz32Bit(mz_32_bit); }, "mz_32_bit"_a, "Sets if mz-data and rt-data should be stored with 32bit or 64bit precision")
        .def("getMz32Bit", [](const OpenMS::PeakFileOptions& self) { return self.getMz32Bit(); }, "Returns true, if mz-data and rt-data should be stored with 32bit precision")
        .def("setIntensity32Bit", [](OpenMS::PeakFileOptions& self, bool int_32_bit) { return self.setIntensity32Bit(int_32_bit); }, "int_32_bit"_a, "Sets if intensity data should be stored with 32bit or 64bit precision")
        .def("getIntensity32Bit", [](const OpenMS::PeakFileOptions& self) { return self.getIntensity32Bit(); }, "Returns true, if intensity data should be stored with 32bit precision")
        .def("getWriteIndex", [](const OpenMS::PeakFileOptions& self) { return self.getWriteIndex(); }, "Returns whether to write an index at the end of the file (e.g. indexedmzML file format)")
        .def("setWriteIndex", [](OpenMS::PeakFileOptions& self, bool write_index) { return self.setWriteIndex(write_index); }, "write_index"_a, "Returns whether to write an index at the end of the file (e.g. indexedmzML file format)")
        .def("getNumpressConfigurationMassTime", [](const OpenMS::PeakFileOptions& self) { return self.getNumpressConfigurationMassTime(); }, "Sets numpress configuration options for m/z or rt dimension")
        .def("setNumpressConfigurationMassTime", [](OpenMS::PeakFileOptions& self, OpenMS::MSNumpressCoder::NumpressConfig config) { return self.setNumpressConfigurationMassTime(config); }, "config"_a, "Returns numpress configuration options for m/z or rt dimension")
        .def("getNumpressConfigurationIntensity", [](const OpenMS::PeakFileOptions& self) { return self.getNumpressConfigurationIntensity(); }, "Sets numpress configuration options for intensity dimension")
        .def("setNumpressConfigurationIntensity", [](OpenMS::PeakFileOptions& self, OpenMS::MSNumpressCoder::NumpressConfig config) { return self.setNumpressConfigurationIntensity(config); }, "config"_a, "Returns numpress configuration options for intensity dimension")
        .def("getNumpressConfigurationFloatDataArray", [](const OpenMS::PeakFileOptions& self) { return self.getNumpressConfigurationFloatDataArray(); }, "Sets numpress configuration options for float data arrays")
        .def("setNumpressConfigurationFloatDataArray", [](OpenMS::PeakFileOptions& self, OpenMS::MSNumpressCoder::NumpressConfig config) { return self.setNumpressConfigurationFloatDataArray(config); }, "config"_a, "Returns numpress configuration options for float data arrays")
        .def("getMaxDataPoolSize", [](const OpenMS::PeakFileOptions& self) { return self.getMaxDataPoolSize(); }, "Returns maximal size of the data pool")
        .def("setMaxDataPoolSize", [](OpenMS::PeakFileOptions& self, size_t size) { return self.setMaxDataPoolSize(size); }, "size"_a, "Sets maximal size of the data pool")
        .def("hasFilters", [](const OpenMS::PeakFileOptions& self) { return self.hasFilters(); })
        ;


    // -----------------------------------------------------------------------
    // DRange1
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DRange<1>>(m, "DRange1", "OpenMS class DRange1")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DRange<1>&>())
        .def("__copy__", [](const OpenMS::DRange<1>& self) { return OpenMS::DRange<1>(self); })
        .def("__deepcopy__", [](const OpenMS::DRange<1>& self, nb::dict) { return OpenMS::DRange<1>(self); }, "memo"_a)
        .def("__init__", [](OpenMS::DRange<1>* self, double min_val, double max_val) {
            new (self) OpenMS::DRange<1>(OpenMS::DPosition<1>(min_val), OpenMS::DPosition<1>(max_val));
        }, "min"_a, "max"_a)
        .def("__eq__", [](const OpenMS::DRange<1>& self, const OpenMS::DRange<1>& other) {
            return self == other;
        })
        .def("encloses", [](const OpenMS::DRange<1>& self, double pos) {
            return self.encloses(OpenMS::DPosition<1>(pos));
        }, "position"_a, "Check if position is within this range")
        .def("united", [](const OpenMS::DRange<1>& self, const OpenMS::DRange<1>& other) {
            return self.united(other);
        }, "other"_a, "Returns the union of this range with another")
        .def("isIntersected", [](const OpenMS::DRange<1>& self, const OpenMS::DRange<1>& other) {
            return self.isIntersected(other);
        }, "other"_a, "Check if two ranges intersect")
        .def("isEmpty", [](const OpenMS::DRange<1>& self) { return self.isEmpty(); })
        .def("minX", [](const OpenMS::DRange<1>& self) { return self.minPosition()[0]; })
        .def("maxX", [](const OpenMS::DRange<1>& self) { return self.maxPosition()[0]; })
        .def("setMinX", [](OpenMS::DRange<1>& self, double c) {
            self.min_[0] = c;
        }, "c"_a)
        .def("setMaxX", [](OpenMS::DRange<1>& self, double c) {
            self.max_[0] = c;
        }, "c"_a)
        .def("__repr__", [](const OpenMS::DRange<1>& self) {
            std::ostringstream oss;
            oss << std::showpoint << "DRange1(" << self.minPosition()[0] << ", " << self.maxPosition()[0] << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::DRange<1>& self) {
            std::ostringstream oss;
            oss << std::showpoint << "DRange1(" << self.minPosition()[0] << ", " << self.maxPosition()[0] << ")";
            return oss.str();
        })
        ;

    // -----------------------------------------------------------------------
    // DRangeIntersection (DRange::DRangeIntersection)
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::DRange<1>::DRangeIntersection>(m, "DRangeIntersection",
        "Result of range intersection check", nb::is_arithmetic())
        .value("Disjoint", OpenMS::DRange<1>::Disjoint)
        .value("Intersects", OpenMS::DRange<1>::Intersects)
        .value("Inside", OpenMS::DRange<1>::Inside)

        ;

    // -----------------------------------------------------------------------
    // DRange2
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::DRange<2>>(m, "DRange2", "Two-dimensional range")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::DRange<2>&>())
        .def("__copy__", [](const OpenMS::DRange<2>& self) { return OpenMS::DRange<2>(self); })
        .def("__deepcopy__", [](const OpenMS::DRange<2>& self, nb::dict) { return OpenMS::DRange<2>(self); }, "memo"_a)
        .def("__init__", [](OpenMS::DRange<2>* self, double minx, double miny, double maxx, double maxy) {
            new (self) OpenMS::DRange<2>(minx, miny, maxx, maxy);
        }, "min_x"_a, "min_y"_a, "max_x"_a, "max_y"_a)
        .def("__eq__", [](const OpenMS::DRange<2>& self, const OpenMS::DRange<2>& other) { return self == other; })
        .def("encloses", [](const OpenMS::DRange<2>& self, double x, double y) {
            return self.encloses(OpenMS::DPosition<2>(x, y));
        }, "x"_a, "y"_a, "Check if position is within this range")
        .def("united", [](const OpenMS::DRange<2>& self, const OpenMS::DRange<2>& other) { return self.united(other); }, "other"_a, "Returns the union of this range with another")
        .def("isIntersected", [](const OpenMS::DRange<2>& self, const OpenMS::DRange<2>& other) { return self.isIntersected(other); }, "other"_a, "Check if two ranges intersect")
        .def("isEmpty", [](const OpenMS::DRange<2>& self) { return self.isEmpty(); }, "Returns true if the range is empty")
        .def("minX", [](const OpenMS::DRange<2>& self) { return self.minPosition()[0]; })
        .def("minY", [](const OpenMS::DRange<2>& self) { return self.minPosition()[1]; })
        .def("maxX", [](const OpenMS::DRange<2>& self) { return self.maxPosition()[0]; })
        .def("maxY", [](const OpenMS::DRange<2>& self) { return self.maxPosition()[1]; })
        .def("setMinX", [](OpenMS::DRange<2>& self, double c) { self.min_[0] = c; }, "c"_a)
        .def("setMinY", [](OpenMS::DRange<2>& self, double c) { self.min_[1] = c; }, "c"_a)
        .def("setMaxX", [](OpenMS::DRange<2>& self, double c) { self.max_[0] = c; }, "c"_a)
        .def("setMaxY", [](OpenMS::DRange<2>& self, double c) { self.max_[1] = c; }, "c"_a)
        .def("__repr__", [](const OpenMS::DRange<2>& self) {
            std::ostringstream oss;
            oss << std::showpoint << "DRange2((" << self.minPosition()[0] << ", " << self.minPosition()[1]
                << "), (" << self.maxPosition()[0] << ", " << self.maxPosition()[1] << "))";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::DRange<2>& self) {
            std::ostringstream oss;
            oss << std::showpoint << "DRange2((" << self.minPosition()[0] << ", " << self.minPosition()[1]
                << "), (" << self.maxPosition()[0] << ", " << self.maxPosition()[1] << "))";
            return oss.str();
        })
        ;

    // -----------------------------------------------------------------------
    // PeakTypeEstimator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PeakTypeEstimator>(m, "PeakTypeEstimator", "Estimates if the data of a spectrum is raw data or peak data")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PeakTypeEstimator &>())
        .def("__copy__", [](const OpenMS::PeakTypeEstimator& self) { return OpenMS::PeakTypeEstimator(self); })
        .def("__deepcopy__", [](const OpenMS::PeakTypeEstimator& self, nb::dict) { return OpenMS::PeakTypeEstimator(self); }, "memo"_a)
        .def("estimateType", [](const OpenMS::PeakTypeEstimator& /*self*/, const OpenMS::MSSpectrum& spec) {
            return OpenMS::PeakTypeEstimator::estimateType(spec.begin(), spec.end());
        }, "spectrum"_a, "Estimates the peak type of the given spectrum")
        ;

    // -----------------------------------------------------------------------
    // PercolatorInfile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PercolatorInfile>(m, "PercolatorInfile", "Class for storing Percolator tab-delimited input files")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::PercolatorInfile &>())
        .def("__copy__", [](const OpenMS::PercolatorInfile& self) { return OpenMS::PercolatorInfile(self); })
        .def("__deepcopy__", [](const OpenMS::PercolatorInfile& self, nb::dict) { return OpenMS::PercolatorInfile(self); }, "memo"_a)
        .def_static("store", [](const OpenMS::String& pin_file, const OpenMS::PeptideIdentificationList& peptide_ids, const std::vector<OpenMS::String>& feature_set, const OpenMS::String& enz, int min_charge, int max_charge) { return OpenMS::PercolatorInfile::store(pin_file, peptide_ids, feature_set, enz, min_charge, max_charge); }, "pin_file"_a, "peptide_ids"_a, "feature_set"_a, "enz"_a, "min_charge"_a, "max_charge"_a)
        ;

    // -----------------------------------------------------------------------
    // PercolatorOutfile
    // -----------------------------------------------------------------------
    auto percolatoroutfile_class = nb::class_<OpenMS::PercolatorOutfile>(m, "PercolatorOutfile", 
        R"doc(
Class for reading Percolator tab-delimited output files
For PSM-level output, the file extension should be ".psms"
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::PercolatorOutfile& self) { return OpenMS::PercolatorOutfile(self); })
        .def("__deepcopy__", [](const OpenMS::PercolatorOutfile& self, nb::dict) { return OpenMS::PercolatorOutfile(self); }, "memo"_a)
        .def_static("getScoreType", [](OpenMS::String score_type_name) { return OpenMS::PercolatorOutfile::getScoreType(score_type_name); }, "score_type_name"_a, "Returns a score type given its name")
        .def("load", [](OpenMS::PercolatorOutfile& self, const OpenMS::String& filename, OpenMS::ProteinIdentification& proteins, OpenMS::PeptideIdentificationList& peptides, OpenMS::SpectrumMetaDataLookup& lookup, OpenMS::PercolatorOutfile::ScoreType output_score) { return self.load(filename, proteins, peptides, lookup, output_score); }, "filename"_a, "proteins"_a, "peptides"_a, "lookup"_a, "output_score"_a, "Loads a Percolator output file")
        ;
    // ScoreType enum nested under PercolatorOutfile
    nb::enum_<OpenMS::PercolatorOutfile::ScoreType>(percolatoroutfile_class, "ScoreType", nb::is_arithmetic())
        .value("QVALUE", OpenMS::PercolatorOutfile::ScoreType::QVALUE)
        .value("POSTERRPROB", OpenMS::PercolatorOutfile::ScoreType::POSTERRPROB)
        .value("SCORE", OpenMS::PercolatorOutfile::ScoreType::SCORE)
        .value("SIZE_OF_SCORETYPE", OpenMS::PercolatorOutfile::ScoreType::SIZE_OF_SCORETYPE)
        ;

    // -----------------------------------------------------------------------
    // PlainMSDataWritingConsumer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::PlainMSDataWritingConsumer>(m, "PlainMSDataWritingConsumer", "Consumer class that writes MS data to disk using the mzML format")
        .def(nb::init<OpenMS::String>())
        .def("setExperimentalSettings", [](OpenMS::PlainMSDataWritingConsumer& self, const OpenMS::ExperimentalSettings& exp) { return self.setExperimentalSettings(exp); }, "exp"_a)
        .def("setExpectedSize", [](OpenMS::PlainMSDataWritingConsumer& self, size_t expectedSpectra, size_t expectedChromatograms) { return self.setExpectedSize(expectedSpectra, expectedChromatograms); }, "expectedSpectra"_a, "expectedChromatograms"_a, 
            R"doc(
Set experimental settings for the whole file
:param exp: Experimental settings to be used for this file (from this and the first spectrum/chromatogram, the class will deduce most of the header of the mzML file)
)doc")
        .def("consumeSpectrum", [](OpenMS::PlainMSDataWritingConsumer& self, OpenMS::MSSpectrum& s) { return self.consumeSpectrum(s); }, "s"_a)
        .def("consumeChromatogram", [](OpenMS::PlainMSDataWritingConsumer& self, OpenMS::MSChromatogram& c) { return self.consumeChromatogram(c); }, "c"_a)
        .def("addDataProcessing", [](OpenMS::PlainMSDataWritingConsumer& self, OpenMS::DataProcessing d) { return self.addDataProcessing(d); }, "d"_a, 
            R"doc(
Set expected size of spectra and chromatograms to be written
These numbers will be written in the spectrumList and chromatogramList
tag in the mzML file. Therefore, these will contain wrong numbers if
the expected size is not set correctly
:param expectedSpectra: Number of spectra expected
:param expectedChromatograms: Number of chromatograms expected
)doc")
        .def("getNrSpectraWritten", [](OpenMS::PlainMSDataWritingConsumer& self) { return self.getNrSpectraWritten(); }, "Returns the number of spectra written")
        .def("getNrChromatogramsWritten", [](OpenMS::PlainMSDataWritingConsumer& self) { return self.getNrChromatogramsWritten(); }, "Returns the number of chromatograms written")
        .def("setOptions", [](OpenMS::PlainMSDataWritingConsumer& self, const OpenMS::PeakFileOptions& opt) { return self.setOptions(opt); }, "opt"_a)
        .def("getOptions", [](OpenMS::PlainMSDataWritingConsumer& self) -> OpenMS::PeakFileOptions & { return self.getOptions(); }, nb::rv_policy::reference_internal)
        ;

    // -----------------------------------------------------------------------
    // RegularSwathFileConsumer
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::RegularSwathFileConsumer>(m, "RegularSwathFileConsumer", 
        R"doc(
Consumes spectra from a SWATH experiment stored in a single file.
Groups MS2 spectra by their precursor m/z into SWATH windows
FullSwathFileConsumer
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::RegularSwathFileConsumer& self) { return OpenMS::RegularSwathFileConsumer(self); })
        .def("__deepcopy__", [](const OpenMS::RegularSwathFileConsumer& self, nb::dict) { return OpenMS::RegularSwathFileConsumer(self); }, "memo"_a)
        .def(nb::init<std::vector<OpenSwath::SwathMap>>())
        .def("setExpectedSize", [](OpenMS::RegularSwathFileConsumer& self, size_t p0, size_t p1) { return self.setExpectedSize(p0, p1); })
        .def("setExperimentalSettings", [](OpenMS::RegularSwathFileConsumer& self, const OpenMS::ExperimentalSettings& exp) { return self.setExperimentalSettings(exp); }, "exp"_a)
        .def("retrieveSwathMaps", [](OpenMS::RegularSwathFileConsumer& self) { std::vector<OpenSwath::SwathMap> maps; self.retrieveSwathMaps(maps); return maps; })
        .def("consumeChromatogram", [](OpenMS::RegularSwathFileConsumer& self, OpenMS::MSChromatogram& p0) { return self.consumeChromatogram(p0); })
        .def("consumeSpectrum", [](OpenMS::RegularSwathFileConsumer& self, OpenMS::MSSpectrum& s) { return self.consumeSpectrum(s); }, "s"_a)
        ;

    // -----------------------------------------------------------------------
    // SemanticValidator_CVTerm
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Internal::SemanticValidator::CVTerm>(m, "SemanticValidator_CVTerm", "OpenMS class SemanticValidator_CVTerm")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::Internal::SemanticValidator::CVTerm &>())
        .def("__copy__", [](const OpenMS::Internal::SemanticValidator::CVTerm& self) { return OpenMS::Internal::SemanticValidator::CVTerm(self); })
        .def("__deepcopy__", [](const OpenMS::Internal::SemanticValidator::CVTerm& self, nb::dict) { return OpenMS::Internal::SemanticValidator::CVTerm(self); }, "memo"_a)
        .def_rw("accession", &OpenMS::Internal::SemanticValidator::CVTerm::accession)
        .def_rw("name", &OpenMS::Internal::SemanticValidator::CVTerm::name)
        .def_rw("value", &OpenMS::Internal::SemanticValidator::CVTerm::value)
        .def_rw("has_value", &OpenMS::Internal::SemanticValidator::CVTerm::has_value)
        .def_rw("unit_accession", &OpenMS::Internal::SemanticValidator::CVTerm::unit_accession)
        .def_rw("has_unit_accession", &OpenMS::Internal::SemanticValidator::CVTerm::has_unit_accession)
        .def_rw("unit_name", &OpenMS::Internal::SemanticValidator::CVTerm::unit_name)
        .def_rw("has_unit_name", &OpenMS::Internal::SemanticValidator::CVTerm::has_unit_name)
        ;

    // -----------------------------------------------------------------------
    // SequestInfile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SequestInfile>(m, "SequestInfile", "OpenMS class SequestInfile")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SequestInfile &>())
        .def("__copy__", [](const OpenMS::SequestInfile& self) { return OpenMS::SequestInfile(self); })
        .def("__deepcopy__", [](const OpenMS::SequestInfile& self, nb::dict) { return OpenMS::SequestInfile(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("store", [](OpenMS::SequestInfile& self, const OpenMS::String& filename) { return self.store(filename); }, "filename"_a)
        .def("getEnzymeInfoAsString", [](const OpenMS::SequestInfile& self) { return self.getEnzymeInfoAsString(); }, "Returns the enzyme list as a string")
        .def("getDatabase", [](const OpenMS::SequestInfile& self) { return self.getDatabase(); }, "Returns the used database")
        .def("setDatabase", [](OpenMS::SequestInfile& self, const OpenMS::String& database) { return self.setDatabase(database); }, "database"_a, "Sets the used database")
        .def("getNeutralLossesForIons", [](const OpenMS::SequestInfile& self) { return self.getNeutralLossesForIons(); }, "Returns whether neutral losses are considered for the a-, b- and y-ions")
        .def("setNeutralLossesForIons", [](OpenMS::SequestInfile& self, const OpenMS::String& neutral_losses_for_ions) { return self.setNeutralLossesForIons(neutral_losses_for_ions); }, "neutral_losses_for_ions"_a, "Sets whether neutral losses are considered for the a-, b- and y-ions")
        .def("getIonSeriesWeights", [](const OpenMS::SequestInfile& self) { return self.getIonSeriesWeights(); }, "Returns the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series")
        .def("setIonSeriesWeights", [](OpenMS::SequestInfile& self, const OpenMS::String& ion_series_weights) { return self.setIonSeriesWeights(ion_series_weights); }, "ion_series_weights"_a, "Sets the weights for the a-, b-, c-, d-, v-, w-, x-, y- and z-ion series")
        .def("getPartialSequence", [](const OpenMS::SequestInfile& self) { return self.getPartialSequence(); }, "Returns the partial sequences (space delimited) that have to occur in the theoretical spectra")
        .def("setPartialSequence", [](OpenMS::SequestInfile& self, const OpenMS::String& partial_sequence) { return self.setPartialSequence(partial_sequence); }, "partial_sequence"_a, "Sets the partial sequences (space delimited) that have to occur in the theoretical spectra")
        .def("getSequenceHeaderFilter", [](const OpenMS::SequestInfile& self) { return self.getSequenceHeaderFilter(); }, "Returns the sequences (space delimited) that have to occur, or be absent (preceded by a tilde) in the header of a protein to be considered")
        .def("setSequenceHeaderFilter", [](OpenMS::SequestInfile& self, const OpenMS::String& sequence_header_filter) { return self.setSequenceHeaderFilter(sequence_header_filter); }, "sequence_header_filter"_a, "Sets the sequences (space delimited) that have to occur, or be absent (preceded by a tilde) in the header of a protein to be considered")
        .def("getProteinMassFilter", [](const OpenMS::SequestInfile& self) { return self.getProteinMassFilter(); }, "Returns the protein mass filter (either min and max mass, or mass and tolerance value in percent)")
        .def("setProteinMassFilter", [](OpenMS::SequestInfile& self, const OpenMS::String& protein_mass_filter) { return self.setProteinMassFilter(protein_mass_filter); }, "protein_mass_filter"_a, "Sets the protein mass filter (either min and max mass, or mass and tolerance value in percent)")
        .def("getPeakMassTolerance", [](const OpenMS::SequestInfile& self) { return self.getPeakMassTolerance(); }, "Returns the peak mass tolerance")
        .def("setPeakMassTolerance", [](OpenMS::SequestInfile& self, float peak_mass_tolerance) { return self.setPeakMassTolerance(peak_mass_tolerance); }, "peak_mass_tolerance"_a, "Sets the peak mass tolerance")
        .def("getPrecursorMassTolerance", [](const OpenMS::SequestInfile& self) { return self.getPrecursorMassTolerance(); }, "Returns the precursor mass tolerance")
        .def("setPrecursorMassTolerance", [](OpenMS::SequestInfile& self, float precursor_mass_tolerance) { return self.setPrecursorMassTolerance(precursor_mass_tolerance); }, "precursor_mass_tolerance"_a, "Sets the precursor mass tolerance")
        .def("getMatchPeakTolerance", [](const OpenMS::SequestInfile& self) { return self.getMatchPeakTolerance(); }, "Returns the match peak tolerance")
        .def("setMatchPeakTolerance", [](OpenMS::SequestInfile& self, float match_peak_tolerance) { return self.setMatchPeakTolerance(match_peak_tolerance); }, "match_peak_tolerance"_a, "Sets the match peak tolerance")
        .def("getIonCutoffPercentage", [](const OpenMS::SequestInfile& self) { return self.getIonCutoffPercentage(); }, "Returns the the cutoff of the ratio matching theoretical peaks/theoretical peaks")
        .def("setIonCutoffPercentage", [](OpenMS::SequestInfile& self, float ion_cutoff_percentage) { return self.setIonCutoffPercentage(ion_cutoff_percentage); }, "ion_cutoff_percentage"_a, "Sets the ion cutoff of the ratio matching theoretical peaks/theoretical peaks")
        .def("getPeptideMassUnit", [](const OpenMS::SequestInfile& self) { return self.getPeptideMassUnit(); }, "Returns the peptide mass unit")
        .def("setPeptideMassUnit", [](OpenMS::SequestInfile& self, size_t peptide_mass_unit) { return self.setPeptideMassUnit(peptide_mass_unit); }, "peptide_mass_unit"_a, "Sets the peptide mass unit")
        .def("getOutputLines", [](const OpenMS::SequestInfile& self) { return self.getOutputLines(); }, "Returns the number of peptides to be displayed")
        .def("setOutputLines", [](OpenMS::SequestInfile& self, size_t output_lines) { return self.setOutputLines(output_lines); }, "output_lines"_a, "Sets the number of peptides to be displayed")
        .def("getEnzymeNumber", [](const OpenMS::SequestInfile& self) { return self.getEnzymeNumber(); }, "Returns the enzyme used for cleavage (by means of the number from a list of enzymes)")
        .def("getEnzymeName", [](const OpenMS::SequestInfile& self) { return self.getEnzymeName(); }, "Returns the enzyme used for cleavage")
        .def("setEnzyme", [](OpenMS::SequestInfile& self, const OpenMS::String& enzyme_name) { return self.setEnzyme(enzyme_name); }, "enzyme_name"_a, "Sets the enzyme used for cleavage (by means of the number from a list of enzymes)")
        .def("getMaxAAPerModPerPeptide", [](const OpenMS::SequestInfile& self) { return self.getMaxAAPerModPerPeptide(); }, "Returns the maximum number of amino acids containing the same modification in a peptide")
        .def("setMaxAAPerModPerPeptide", [](OpenMS::SequestInfile& self, size_t max_aa_per_mod_per_peptide) { return self.setMaxAAPerModPerPeptide(max_aa_per_mod_per_peptide); }, "max_aa_per_mod_per_peptide"_a, "Sets the maximum number of amino acids containing the same modification in a peptide")
        .def("getMaxModsPerPeptide", [](const OpenMS::SequestInfile& self) { return self.getMaxModsPerPeptide(); }, "Returns the maximum number of modifications that are allowed in a peptide")
        .def("setMaxModsPerPeptide", [](OpenMS::SequestInfile& self, size_t max_mods_per_peptide) { return self.setMaxModsPerPeptide(max_mods_per_peptide); }, "max_mods_per_peptide"_a, "Sets the maximum number of modifications that are allowed in a peptide")
        .def("getNucleotideReadingFrame", [](const OpenMS::SequestInfile& self) { return self.getNucleotideReadingFrame(); }, "Returns the nucleotide reading frame")
        .def("setNucleotideReadingFrame", [](OpenMS::SequestInfile& self, size_t nucleotide_reading_frame) { return self.setNucleotideReadingFrame(nucleotide_reading_frame); }, "nucleotide_reading_frame"_a, "Sets the nucleotide reading frame")
        .def("getMaxInternalCleavageSites", [](const OpenMS::SequestInfile& self) { return self.getMaxInternalCleavageSites(); }, "Returns the maximum number of internal cleavage sites")
        .def("setMaxInternalCleavageSites", [](OpenMS::SequestInfile& self, size_t max_internal_cleavage_sites) { return self.setMaxInternalCleavageSites(max_internal_cleavage_sites); }, "max_internal_cleavage_sites"_a, "Sets the maximum number of internal cleavage sites")
        .def("getMatchPeakCount", [](const OpenMS::SequestInfile& self) { return self.getMatchPeakCount(); }, "Returns the number of top abundant peaks to match with theoretical ones")
        .def("setMatchPeakCount", [](OpenMS::SequestInfile& self, size_t match_peak_count) { return self.setMatchPeakCount(match_peak_count); }, "match_peak_count"_a, "Sets the number of top abundant peaks to with theoretical ones")
        .def("getMatchPeakAllowedError", [](const OpenMS::SequestInfile& self) { return self.getMatchPeakAllowedError(); }, "Returns the number of top abundant peaks that are allowed not to match with a theoretical peak")
        .def("setMatchPeakAllowedError", [](OpenMS::SequestInfile& self, size_t match_peak_allowed_error) { return self.setMatchPeakAllowedError(match_peak_allowed_error); }, "match_peak_allowed_error"_a, "Sets the number of top abundant peaks that are allowed not to match with a theoretical peak")
        .def("getShowFragmentIons", [](const OpenMS::SequestInfile& self) { return self.getShowFragmentIons(); }, "Returns whether fragment ions shall be displayed")
        .def("setShowFragmentIons", [](OpenMS::SequestInfile& self, bool show_fragments) { return self.setShowFragmentIons(show_fragments); }, "show_fragments"_a, "Sets whether fragment ions shall be displayed")
        .def("getPrintDuplicateReferences", [](const OpenMS::SequestInfile& self) { return self.getPrintDuplicateReferences(); }, "Returns whether all proteins containing a found peptide should be displayed")
        .def("setPrintDuplicateReferences", [](OpenMS::SequestInfile& self, bool print_duplicate_references) { return self.setPrintDuplicateReferences(print_duplicate_references); }, "print_duplicate_references"_a, "Sets whether all proteins containing a found peptide should be displayed")
        .def("getRemovePrecursorNearPeaks", [](const OpenMS::SequestInfile& self) { return self.getRemovePrecursorNearPeaks(); }, "Returns whether peaks near (15 amu) the precursor peak are removed")
        .def("setRemovePrecursorNearPeaks", [](OpenMS::SequestInfile& self, bool remove_precursor_near_peaks) { return self.setRemovePrecursorNearPeaks(remove_precursor_near_peaks); }, "remove_precursor_near_peaks"_a, "Sets whether peaks near (15 amu) the precursor peak are removed")
        .def("getMassTypeParent", [](const OpenMS::SequestInfile& self) { return self.getMassTypeParent(); }, "Returns the mass type of the parent (0 - monoisotopic, 1 - average mass)")
        .def("setMassTypeParent", [](OpenMS::SequestInfile& self, bool mass_type_parent) { return self.setMassTypeParent(mass_type_parent); }, "mass_type_parent"_a, "Sets the mass type of the parent (0 - monoisotopic, 1 - average mass)")
        .def("getMassTypeFragment", [](const OpenMS::SequestInfile& self) { return self.getMassTypeFragment(); }, "Returns the mass type of the fragments (0 - monoisotopic, 1 - average mass)")
        .def("setMassTypeFragment", [](OpenMS::SequestInfile& self, bool mass_type_fragment) { return self.setMassTypeFragment(mass_type_fragment); }, "mass_type_fragment"_a, "Sets the mass type of the fragments (0 - monoisotopic, 1 - average mass)")
        .def("getNormalizeXcorr", [](const OpenMS::SequestInfile& self) { return self.getNormalizeXcorr(); }, "Returns whether normalized xcorr values are displayed")
        .def("setNormalizeXcorr", [](OpenMS::SequestInfile& self, bool normalize_xcorr) { return self.setNormalizeXcorr(normalize_xcorr); }, "normalize_xcorr"_a, "Sets whether normalized xcorr values are displayed")
        .def("getResiduesInUpperCase", [](const OpenMS::SequestInfile& self) { return self.getResiduesInUpperCase(); }, "Returns whether residues are in upper case")
        .def("setResiduesInUpperCase", [](OpenMS::SequestInfile& self, bool residues_in_upper_case) { return self.setResiduesInUpperCase(residues_in_upper_case); }, "residues_in_upper_case"_a, "Sets whether residues are in upper case")
        .def("addEnzymeInfo", [](OpenMS::SequestInfile& self, std::vector<OpenMS::String>& enzyme_info) { return self.addEnzymeInfo(enzyme_info); }, "enzyme_info"_a, "Adds an enzyme to the list and sets is as used")
        .def("handlePTMs", [](OpenMS::SequestInfile& self, const OpenMS::String& modification_line, const OpenMS::String& modifications_filename, bool monoisotopic) { return self.handlePTMs(modification_line, modifications_filename, monoisotopic); }, "modification_line"_a, "modifications_filename"_a, "monoisotopic"_a)
        .def("getModifications", [](const OpenMS::SequestInfile& self) -> std::map<OpenMS::String, std::vector<OpenMS::String>> { return self.getModifications(); }, "Returns the modifications (the modification names map to the affected residues, the mass change and the type)")
        ;

    // -----------------------------------------------------------------------
    // SequestOutfile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SequestOutfile>(m, "SequestOutfile", "Representation of a Sequest output file")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SequestOutfile &>())
        .def("__copy__", [](const OpenMS::SequestOutfile& self) { return OpenMS::SequestOutfile(self); })
        .def("__deepcopy__", [](const OpenMS::SequestOutfile& self, nb::dict) { return OpenMS::SequestOutfile(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def("load", [](OpenMS::SequestOutfile& self, const OpenMS::String& result_filename, OpenMS::PeptideIdentificationList& peptide_identifications, OpenMS::ProteinIdentification& protein_identification, double p_value_threshold, const OpenMS::String& database, bool ignore_proteins_per_peptide) { std::vector<double> pvalues; self.load(result_filename, peptide_identifications, protein_identification, p_value_threshold, pvalues, database, ignore_proteins_per_peptide); return pvalues; }, "result_filename"_a, "peptide_identifications"_a, "protein_identification"_a, "p_value_threshold"_a, "database"_a, "ignore_proteins_per_peptide"_a)
        .def("getColumns", [](OpenMS::SequestOutfile& self, const OpenMS::String& line, size_t number_of_columns, size_t reference_column) { std::vector<OpenMS::String> substrings; self.getColumns(line, substrings, number_of_columns, reference_column); return substrings; }, "line"_a, "number_of_columns"_a, "reference_column"_a, "Retrieves columns from a Sequest outfile line")
        .def("getACAndACType", [](OpenMS::SequestOutfile& self, OpenMS::String line) { OpenMS::String accession; OpenMS::String accession_type; self.getACAndACType(line, accession, accession_type); return std::make_tuple(accession, accession_type); }, "line"_a, "Retrieves the accession type and accession number from a protein description line")
        .def("getSequences", [](OpenMS::SequestOutfile& self, const OpenMS::String& database_filename, const std::map<OpenMS::String, size_t>& ac_position_map) {
            std::vector<OpenMS::String> sequences;
            std::vector<std::pair<OpenMS::String, size_t>> found;
            std::map<OpenMS::String, size_t> not_found;
            self.getSequences(database_filename, ac_position_map, sequences, found, not_found);
            return nb::make_tuple(sequences, found, not_found);
        }, "database_filename"_a, "ac_position_map"_a, "Retrieves sequences from a FASTA database. Returns (sequences, found, not_found)")
        ;

    // -----------------------------------------------------------------------
    // SiriusFragmentAnnotation
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SiriusFragmentAnnotation>(m, "SiriusFragmentAnnotation", "OpenMS class SiriusFragmentAnnotation")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::SiriusFragmentAnnotation &>())
        .def("__copy__", [](const OpenMS::SiriusFragmentAnnotation& self) { return OpenMS::SiriusFragmentAnnotation(self); })
        .def("__deepcopy__", [](const OpenMS::SiriusFragmentAnnotation& self, nb::dict) { return OpenMS::SiriusFragmentAnnotation(self); }, "memo"_a)
        .def_static("extractAndResolveSiriusAnnotations", [](const std::vector<OpenMS::String>& sirius_workspace_subdirs, double score_threshold, bool use_exact_mass, bool decoy_generation) { return OpenMS::SiriusFragmentAnnotation::extractAndResolveSiriusAnnotations(sirius_workspace_subdirs, score_threshold, use_exact_mass, decoy_generation); }, "sirius_workspace_subdirs"_a, "score_threshold"_a, "use_exact_mass"_a, "decoy_generation"_a)
        .def_static("extractAnnotationsFromSiriusFile", [](const OpenMS::String& path_to_sirius_workspace, size_t max_rank, bool decoy, bool use_exact_mass) { return OpenMS::SiriusFragmentAnnotation::extractAnnotationsFromSiriusFile(path_to_sirius_workspace, max_rank, decoy, use_exact_mass); }, "path_to_sirius_workspace"_a, "max_rank"_a, "decoy"_a, "use_exact_mass"_a)
        .def_static("extract_columnname_to_columnindex", [](const OpenMS::CsvFile& csvfile) { return OpenMS::SiriusFragmentAnnotation::extract_columnname_to_columnindex(csvfile); }, "csvfile"_a)
        ;

    // -----------------------------------------------------------------------
    // SiriusFragmentAnnotation_SiriusTargetDecoySpectra
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SiriusFragmentAnnotation::SiriusTargetDecoySpectra>(m, "SiriusFragmentAnnotation_SiriusTargetDecoySpectra", "OpenMS class SiriusFragmentAnnotation_SiriusTargetDecoySpectra")
        .def(nb::init<>())
        .def(nb::init<OpenMS::MSSpectrum, OpenMS::MSSpectrum>())
        .def("__copy__", [](const OpenMS::SiriusFragmentAnnotation::SiriusTargetDecoySpectra& self) { return OpenMS::SiriusFragmentAnnotation::SiriusTargetDecoySpectra(self); })
        .def("__deepcopy__", [](const OpenMS::SiriusFragmentAnnotation::SiriusTargetDecoySpectra& self, nb::dict) { return OpenMS::SiriusFragmentAnnotation::SiriusTargetDecoySpectra(self); }, "memo"_a)
        .def_rw("target", &OpenMS::SiriusFragmentAnnotation::SiriusTargetDecoySpectra::target)
        .def_rw("decoy", &OpenMS::SiriusFragmentAnnotation::SiriusTargetDecoySpectra::decoy)
        ;

    // -----------------------------------------------------------------------
    // SqMassConfig (SqMassFile::SqMassConfig)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SqMassFile::SqMassConfig>(m, "SqMassConfig",
        "Configuration for SqMassFile write operations")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::SqMassFile::SqMassConfig& self) { return OpenMS::SqMassFile::SqMassConfig(self); })
        .def("__deepcopy__", [](const OpenMS::SqMassFile::SqMassConfig& self, nb::dict) { return OpenMS::SqMassFile::SqMassConfig(self); }, "memo"_a)
        .def_rw("write_full_meta", &OpenMS::SqMassFile::SqMassConfig::write_full_meta)
        .def_rw("use_lossy_numpress", &OpenMS::SqMassFile::SqMassConfig::use_lossy_numpress)
        .def_rw("linear_fp_mass_acc", &OpenMS::SqMassFile::SqMassConfig::linear_fp_mass_acc)
        ;

    // -----------------------------------------------------------------------
    // SqMassFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SqMassFile>(m, "SqMassFile", 
        R"doc(
An class that uses on-disk SQLite database to read and write spectra and chromatograms
This class provides functions to read and write spectra and chromatograms
to disk using a SQLite database and store them in sqMass format. This
allows users to access, select and filter spectra and chromatograms
on-demand even in a large collection of data
)doc")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::SqMassFile& self) { return OpenMS::SqMassFile(self); })
        .def("__deepcopy__", [](const OpenMS::SqMassFile& self, nb::dict) { return OpenMS::SqMassFile(self); }, "memo"_a)
        .def("load", [](const OpenMS::SqMassFile& self, const OpenMS::String& filename) { OpenMS::MSExperiment map; self.load(filename, map); return map; }, "filename"_a, "Read / Write a complete mass spectrometric experiment")
        .def("store", [](const OpenMS::SqMassFile& self, const OpenMS::String& filename, const OpenMS::MSExperiment& map) { return self.store(filename, map); }, "filename"_a, "map"_a, "Store an MSExperiment in sqMass format")
        .def("setConfig", [](OpenMS::SqMassFile& self, const OpenMS::SqMassFile::SqMassConfig& config) { return self.setConfig(config); }, "config"_a)
        ;

    // -----------------------------------------------------------------------
    // TargetedDataFileLoader
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::TargetedDataFileLoader>(m, "TargetedDataFileLoader", 
        R"doc(
Dispatcher that detects whether an mzML contains spectra (SWATH/DIA)
or chromatograms only (SRM/MRM) and forwards to the appropriate loader.
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::TargetedDataFileLoader &>())
        .def("__copy__", [](const OpenMS::TargetedDataFileLoader& self) { return OpenMS::TargetedDataFileLoader(self); })
        .def("__deepcopy__", [](const OpenMS::TargetedDataFileLoader& self, nb::dict) { return OpenMS::TargetedDataFileLoader(self); }, "memo"_a)
        .def_static("loadFile", [](const OpenMS::String& file, const OpenMS::String& tmp, std::shared_ptr<OpenMS::ExperimentalSettings>& exp_meta, const OpenMS::String& readoptions, OpenMS::Interfaces::IMSDataConsumer * plugin_consumer) { return OpenMS::TargetedDataFileLoader::loadFile(file, tmp, exp_meta, readoptions, plugin_consumer); }, "file"_a, "tmp"_a, "exp_meta"_a, "readoptions"_a, "plugin_consumer"_a)
        ;

    // -----------------------------------------------------------------------
    // ToolInfo
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ToolInfo>(m, "ToolInfo", "OpenMS class ToolInfo")
        .def(nb::init<const OpenMS::ToolInfo &>())
        .def("__copy__", [](const OpenMS::ToolInfo& self) { return OpenMS::ToolInfo(self); })
        .def("__deepcopy__", [](const OpenMS::ToolInfo& self, nb::dict) { return OpenMS::ToolInfo(self); }, "memo"_a)
        .def_rw("version_", &OpenMS::ToolInfo::version_)
        .def_rw("name_", &OpenMS::ToolInfo::name_)
        .def_rw("docurl_", &OpenMS::ToolInfo::docurl_)
        .def_rw("category_", &OpenMS::ToolInfo::category_)
        .def_rw("description_", &OpenMS::ToolInfo::description_)
        .def_rw("citations_", &OpenMS::ToolInfo::citations_)
        ;

    // -----------------------------------------------------------------------
    // XICParquetFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::XICParquetFile>(m, "XICParquetFile", "OpenMS class XICParquetFile")
        .def(nb::init<const OpenMS::XICParquetFile &>())
        .def(nb::init<OpenMS::String>())
        .def(nb::init<std::vector<OpenMS::String>>())
        .def("__copy__", [](const OpenMS::XICParquetFile& self) { return OpenMS::XICParquetFile(self); })
        .def("__deepcopy__", [](const OpenMS::XICParquetFile& self, nb::dict) { return OpenMS::XICParquetFile(self); }, "memo"_a)
        .def("getFilename", [](const OpenMS::XICParquetFile& self) { return self.getFilename(); }, "Reader for multiple OpenSWATH chromatogram Parquet files (.xic).")
        .def("getFilenames", [](const OpenMS::XICParquetFile& self) -> const std::vector<OpenMS::String> & { return self.getFilenames(); }, nb::rv_policy::reference_internal)

        .def("getColumns", [](const OpenMS::XICParquetFile& self) {
            std::vector<OpenMS::String> columns;
            self.getColumns(columns);
            nb::list result;
            for (const auto& col : columns) {
                result.append(nb::str(col.c_str()));
            }
            return result;
        }, "Return parquet schema column names as a list")

        .def("getRuns", [](const OpenMS::XICParquetFile& self) {
            std::vector<OpenMS::XICParquetFile::XICRunInfo> runs;
            self.getRuns(runs);
            nb::list run_ids, source_files;
            for (const auto& r : runs) {
                run_ids.append(r.run_id);
                source_files.append(nb::str(r.source_file.c_str()));
            }
            nb::dict result;
            result["run_id"] = run_ids;
            result["source_file"] = source_files;
            return result;
        }, "Return unique run metadata as a dict")

        .def("getAnalytes", [](const OpenMS::XICParquetFile& self, bool nest_transitions, nb::object columns_obj) {
            std::vector<OpenMS::String> columns;
            std::set<std::string> requested_cols;
            bool filter_cols = !columns_obj.is_none();
            if (filter_cols) {
                for (auto item : columns_obj) {
                    std::string col = nb::cast<std::string>(item);
                    columns.push_back(col);
                    // Store lowercase version for result dict filtering
                    std::string lower_col = col;
                    std::transform(lower_col.begin(), lower_col.end(), lower_col.begin(), ::tolower);
                    requested_cols.insert(lower_col);
                }
            }
            std::vector<OpenMS::XICParquetFile::XICAnalyte> analytes;
            self.getAnalytes(analytes, columns, nest_transitions);

            nb::list precursor_id_list, modified_sequence_list, precursor_charge_list, precursor_decoy_list;
            nb::list transition_id_list, product_charge_list, transition_ordinal_list;
            nb::list detecting_transition_list, product_decoy_list, transition_type_list, annotation_list;

            auto want = [&](const std::string& col) { return !filter_cols || requested_cols.count(col) > 0; };

            for (const auto& a : analytes) {
                if (want("precursor_id")) precursor_id_list.append(a.has_precursor_id ? nb::cast(a.precursor_id) : nb::none());
                if (want("modified_sequence")) modified_sequence_list.append(nb::str(a.modified_sequence.c_str()));
                if (want("precursor_charge")) precursor_charge_list.append(a.has_precursor_charge ? nb::cast(a.precursor_charge) : nb::none());
                if (want("precursor_decoy")) precursor_decoy_list.append(a.has_precursor_decoy ? nb::cast(a.precursor_decoy) : nb::none());

                if (!nest_transitions) {
                    if (want("transition_id")) transition_id_list.append(a.has_transition_id ? nb::cast(a.transition_id) : nb::none());
                    if (want("product_charge")) product_charge_list.append(a.has_product_charge ? nb::cast(a.product_charge) : nb::none());
                    if (want("transition_ordinal")) transition_ordinal_list.append(a.has_transition_ordinal ? nb::cast(a.transition_ordinal) : nb::none());
                    if (want("detecting_transition")) detecting_transition_list.append(a.has_detecting_transition ? nb::cast(a.detecting_transition) : nb::none());
                    if (want("product_decoy")) product_decoy_list.append(a.has_product_decoy ? nb::cast(a.product_decoy) : nb::none());
                    if (want("transition_type")) transition_type_list.append(nb::str(a.transition_type.c_str()));
                    if (want("annotation")) annotation_list.append(nb::str(a.annotation.c_str()));
                } else {
                    nb::list t_ids, p_charges, t_ordinals, d_transitions, p_decoys, t_types, annots;
                    for (auto v : a.transition_ids) t_ids.append(v >= 0 ? nb::cast(v) : nb::none());
                    for (auto v : a.product_charges) p_charges.append(v >= 0 ? nb::cast(v) : nb::none());
                    for (auto v : a.transition_ordinals) t_ordinals.append(v >= 0 ? nb::cast(v) : nb::none());
                    for (auto v : a.detecting_transitions) d_transitions.append(v >= 0 ? nb::cast(v) : nb::none());
                    for (auto v : a.product_decoys) p_decoys.append(v >= 0 ? nb::cast(v) : nb::none());
                    for (const auto& s : a.transition_types) t_types.append(nb::str(s.c_str()));
                    for (const auto& s : a.annotations) annots.append(nb::str(s.c_str()));
                    if (want("transition_id")) transition_id_list.append(t_ids);
                    if (want("product_charge")) product_charge_list.append(p_charges);
                    if (want("transition_ordinal")) transition_ordinal_list.append(t_ordinals);
                    if (want("detecting_transition")) detecting_transition_list.append(d_transitions);
                    if (want("product_decoy")) product_decoy_list.append(p_decoys);
                    if (want("transition_type")) transition_type_list.append(t_types);
                    if (want("annotation")) annotation_list.append(annots);
                }
            }

            nb::dict result;
            if (want("precursor_id")) result["precursor_id"] = precursor_id_list;
            if (want("modified_sequence")) result["modified_sequence"] = modified_sequence_list;
            if (want("precursor_charge")) result["precursor_charge"] = precursor_charge_list;
            if (want("precursor_decoy")) result["precursor_decoy"] = precursor_decoy_list;
            if (want("transition_id")) result["transition_id"] = transition_id_list;
            if (want("product_charge")) result["product_charge"] = product_charge_list;
            if (want("transition_ordinal")) result["transition_ordinal"] = transition_ordinal_list;
            if (want("detecting_transition")) result["detecting_transition"] = detecting_transition_list;
            if (want("product_decoy")) result["product_decoy"] = product_decoy_list;
            if (want("transition_type")) result["transition_type"] = transition_type_list;
            if (want("annotation")) result["annotation"] = annotation_list;
            return result;
        }, "nest_transitions"_a = true, "columns"_a = nb::none(),
           "Return unique analyte metadata as a dict")

        .def("getChromatograms", [](const OpenMS::XICParquetFile& self,
                                    int64_t precursor_id, int64_t transition_id,
                                    const std::string& modified_sequence,
                                    int64_t precursor_charge, int64_t product_charge,
                                    int64_t ms_level, int64_t run_id,
                                    const std::string& filter, bool explode) {
            std::vector<OpenMS::XICParquetFile::XICChromatogram> chroms;
            self.getChromatograms(chroms, precursor_id, transition_id,
                                  OpenMS::String(modified_sequence),
                                  precursor_charge, product_charge,
                                  ms_level, run_id, OpenMS::String(filter));

            nb::list run_id_list, source_file_list, ms_level_list;
            nb::list precursor_id_list, transition_id_list, modified_sequence_list;
            nb::list precursor_charge_list, product_charge_list, detecting_transition_list;
            nb::list precursor_decoy_list, product_decoy_list, transition_ordinal_list;
            nb::list transition_type_list, annotation_list, rt_list, intensity_list;

            for (const auto& c : chroms) {
                if (explode) {
                    if (c.rt.empty()) continue;
                    for (size_t j = 0; j < c.rt.size(); ++j) {
                        run_id_list.append(c.run_id);
                        source_file_list.append(nb::str(c.source_file.c_str()));
                        ms_level_list.append(c.ms_level);
                        precursor_id_list.append(c.has_precursor_id ? nb::cast(c.precursor_id) : nb::none());
                        transition_id_list.append(c.has_transition_id ? nb::cast(c.transition_id) : nb::none());
                        modified_sequence_list.append(nb::str(c.modified_sequence.c_str()));
                        precursor_charge_list.append(c.has_precursor_charge ? nb::cast(c.precursor_charge) : nb::none());
                        product_charge_list.append(c.has_product_charge ? nb::cast(c.product_charge) : nb::none());
                        detecting_transition_list.append(c.has_detecting_transition ? nb::cast(c.detecting_transition) : nb::none());
                        precursor_decoy_list.append(c.has_precursor_decoy ? nb::cast(c.precursor_decoy) : nb::none());
                        product_decoy_list.append(c.has_product_decoy ? nb::cast(c.product_decoy) : nb::none());
                        transition_ordinal_list.append(c.has_transition_ordinal ? nb::cast(c.transition_ordinal) : nb::none());
                        transition_type_list.append(nb::str(c.transition_type.c_str()));
                        annotation_list.append(nb::str(c.annotation.c_str()));
                        rt_list.append(c.rt[j]);
                        intensity_list.append(c.intensity[j]);
                    }
                } else {
                    run_id_list.append(c.run_id);
                    source_file_list.append(nb::str(c.source_file.c_str()));
                    ms_level_list.append(c.ms_level);
                    precursor_id_list.append(c.has_precursor_id ? nb::cast(c.precursor_id) : nb::none());
                    transition_id_list.append(c.has_transition_id ? nb::cast(c.transition_id) : nb::none());
                    modified_sequence_list.append(nb::str(c.modified_sequence.c_str()));
                    precursor_charge_list.append(c.has_precursor_charge ? nb::cast(c.precursor_charge) : nb::none());
                    product_charge_list.append(c.has_product_charge ? nb::cast(c.product_charge) : nb::none());
                    detecting_transition_list.append(c.has_detecting_transition ? nb::cast(c.detecting_transition) : nb::none());
                    precursor_decoy_list.append(c.has_precursor_decoy ? nb::cast(c.precursor_decoy) : nb::none());
                    product_decoy_list.append(c.has_product_decoy ? nb::cast(c.product_decoy) : nb::none());
                    transition_ordinal_list.append(c.has_transition_ordinal ? nb::cast(c.transition_ordinal) : nb::none());
                    transition_type_list.append(nb::str(c.transition_type.c_str()));
                    annotation_list.append(nb::str(c.annotation.c_str()));
                    nb::list rt_vals, int_vals;
                    for (auto v : c.rt) rt_vals.append(v);
                    for (auto v : c.intensity) int_vals.append(v);
                    rt_list.append(rt_vals);
                    intensity_list.append(int_vals);
                }
            }

            nb::dict result;
            result["run_id"] = run_id_list;
            result["source_file"] = source_file_list;
            result["ms_level"] = ms_level_list;
            result["precursor_id"] = precursor_id_list;
            result["transition_id"] = transition_id_list;
            result["modified_sequence"] = modified_sequence_list;
            result["precursor_charge"] = precursor_charge_list;
            result["product_charge"] = product_charge_list;
            result["detecting_transition"] = detecting_transition_list;
            result["precursor_decoy"] = precursor_decoy_list;
            result["product_decoy"] = product_decoy_list;
            result["transition_ordinal"] = transition_ordinal_list;
            result["transition_type"] = transition_type_list;
            result["annotation"] = annotation_list;
            result["rt"] = rt_list;
            result["intensity"] = intensity_list;
            return result;
        }, "precursor_id"_a = -1, "transition_id"_a = -1, "modified_sequence"_a = "",
           "precursor_charge"_a = -1, "product_charge"_a = -1, "ms_level"_a = -1,
           "run_id"_a = -1, "filter"_a = "", "explode"_a = false,
           "Return chromatogram data as a dict")
        .def("__repr__", [](const OpenMS::XICParquetFile& self) {
            const auto& files = self.getFilenames();
            return "XICParquetFile(n_files=" + std::to_string(files.size()) + ")";
        })
        .def("__str__", [](const OpenMS::XICParquetFile& self) { return nb::cast(self).attr("__repr__")(); })
        ;

    // -----------------------------------------------------------------------
    // XIMParquetFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::XIMParquetFile>(m, "XIMParquetFile", "OpenMS class XIMParquetFile")
        .def(nb::init<const OpenMS::XIMParquetFile &>())
        .def(nb::init<OpenMS::String>())
        .def(nb::init<std::vector<OpenMS::String>>())
        .def("__copy__", [](const OpenMS::XIMParquetFile& self) { return OpenMS::XIMParquetFile(self); })
        .def("__deepcopy__", [](const OpenMS::XIMParquetFile& self, nb::dict) { return OpenMS::XIMParquetFile(self); }, "memo"_a)
        .def("getFilename", [](const OpenMS::XIMParquetFile& self) { return self.getFilename(); }, "Reader for multiple OpenSWATH mobilogram Parquet files (.xim).")
        .def("getFilenames", [](const OpenMS::XIMParquetFile& self) -> const std::vector<OpenMS::String> & { return self.getFilenames(); }, nb::rv_policy::reference_internal)

        .def("getColumns", [](const OpenMS::XIMParquetFile& self) {
            std::vector<OpenMS::String> columns;
            self.getColumns(columns);
            nb::list result;
            for (const auto& col : columns) {
                result.append(nb::str(col.c_str()));
            }
            return result;
        }, "Return parquet schema column names as a list")

        .def("getRuns", [](const OpenMS::XIMParquetFile& self) {
            std::vector<OpenMS::XIMParquetFile::XIMRunInfo> runs;
            self.getRuns(runs);
            nb::list run_ids, source_files;
            for (const auto& r : runs) {
                run_ids.append(r.run_id);
                source_files.append(nb::str(r.source_file.c_str()));
            }
            nb::dict result;
            result["run_id"] = run_ids;
            result["source_file"] = source_files;
            return result;
        }, "Return unique run metadata as a dict")

        .def("getAnalytes", [](const OpenMS::XIMParquetFile& self, bool nest_transitions, nb::object columns_obj) {
            std::vector<OpenMS::String> columns;
            std::set<std::string> requested_cols;
            bool filter_cols = !columns_obj.is_none();
            if (filter_cols) {
                for (auto item : columns_obj) {
                    std::string col = nb::cast<std::string>(item);
                    columns.push_back(col);
                    std::string lower_col = col;
                    std::transform(lower_col.begin(), lower_col.end(), lower_col.begin(), ::tolower);
                    requested_cols.insert(lower_col);
                }
            }
            std::vector<OpenMS::XIMParquetFile::XIMAnalyte> analytes;
            self.getAnalytes(analytes, columns, nest_transitions);

            nb::list precursor_id_list, modified_sequence_list, precursor_charge_list, precursor_decoy_list;
            nb::list transition_id_list, product_charge_list, transition_ordinal_list;
            nb::list detecting_transition_list, product_decoy_list, transition_type_list, annotation_list;

            auto want = [&](const std::string& col) { return !filter_cols || requested_cols.count(col) > 0; };

            for (const auto& a : analytes) {
                if (want("precursor_id")) precursor_id_list.append(a.has_precursor_id ? nb::cast(a.precursor_id) : nb::none());
                if (want("modified_sequence")) modified_sequence_list.append(nb::str(a.modified_sequence.c_str()));
                if (want("precursor_charge")) precursor_charge_list.append(a.has_precursor_charge ? nb::cast(a.precursor_charge) : nb::none());
                if (want("precursor_decoy")) precursor_decoy_list.append(a.has_precursor_decoy ? nb::cast(a.precursor_decoy) : nb::none());

                if (!nest_transitions) {
                    if (want("transition_id")) transition_id_list.append(a.has_transition_id ? nb::cast(a.transition_id) : nb::none());
                    if (want("product_charge")) product_charge_list.append(a.has_product_charge ? nb::cast(a.product_charge) : nb::none());
                    if (want("transition_ordinal")) transition_ordinal_list.append(a.has_transition_ordinal ? nb::cast(a.transition_ordinal) : nb::none());
                    if (want("detecting_transition")) detecting_transition_list.append(a.has_detecting_transition ? nb::cast(a.detecting_transition) : nb::none());
                    if (want("product_decoy")) product_decoy_list.append(a.has_product_decoy ? nb::cast(a.product_decoy) : nb::none());
                    if (want("transition_type")) transition_type_list.append(nb::str(a.transition_type.c_str()));
                    if (want("annotation")) annotation_list.append(nb::str(a.annotation.c_str()));
                } else {
                    nb::list t_ids, p_charges, t_ordinals, d_transitions, p_decoys, t_types, annots;
                    for (auto v : a.transition_ids) t_ids.append(v >= 0 ? nb::cast(v) : nb::none());
                    for (auto v : a.product_charges) p_charges.append(v >= 0 ? nb::cast(v) : nb::none());
                    for (auto v : a.transition_ordinals) t_ordinals.append(v >= 0 ? nb::cast(v) : nb::none());
                    for (auto v : a.detecting_transitions) d_transitions.append(v >= 0 ? nb::cast(v) : nb::none());
                    for (auto v : a.product_decoys) p_decoys.append(v >= 0 ? nb::cast(v) : nb::none());
                    for (const auto& s : a.transition_types) t_types.append(nb::str(s.c_str()));
                    for (const auto& s : a.annotations) annots.append(nb::str(s.c_str()));
                    if (want("transition_id")) transition_id_list.append(t_ids);
                    if (want("product_charge")) product_charge_list.append(p_charges);
                    if (want("transition_ordinal")) transition_ordinal_list.append(t_ordinals);
                    if (want("detecting_transition")) detecting_transition_list.append(d_transitions);
                    if (want("product_decoy")) product_decoy_list.append(p_decoys);
                    if (want("transition_type")) transition_type_list.append(t_types);
                    if (want("annotation")) annotation_list.append(annots);
                }
            }

            nb::dict result;
            if (want("precursor_id")) result["precursor_id"] = precursor_id_list;
            if (want("modified_sequence")) result["modified_sequence"] = modified_sequence_list;
            if (want("precursor_charge")) result["precursor_charge"] = precursor_charge_list;
            if (want("precursor_decoy")) result["precursor_decoy"] = precursor_decoy_list;
            if (want("transition_id")) result["transition_id"] = transition_id_list;
            if (want("product_charge")) result["product_charge"] = product_charge_list;
            if (want("transition_ordinal")) result["transition_ordinal"] = transition_ordinal_list;
            if (want("detecting_transition")) result["detecting_transition"] = detecting_transition_list;
            if (want("product_decoy")) result["product_decoy"] = product_decoy_list;
            if (want("transition_type")) result["transition_type"] = transition_type_list;
            if (want("annotation")) result["annotation"] = annotation_list;
            return result;
        }, "nest_transitions"_a = true, "columns"_a = nb::none(),
           "Return unique analyte metadata as a dict")

        .def("getMobilograms", [](const OpenMS::XIMParquetFile& self,
                                    int64_t precursor_id, int64_t transition_id,
                                    const std::string& modified_sequence,
                                    int64_t precursor_charge, int64_t product_charge,
                                    int64_t ms_level, int64_t run_id,
                                    const std::string& mobilogram_type, int64_t feature_id, double feature_rt,
                                    const std::string& filter, bool explode) {
            std::vector<OpenMS::XIMParquetFile::XIMMobilogram> mobilos;
            self.getMobilograms(mobilos, precursor_id, transition_id,
                                  OpenMS::String(modified_sequence),
                                  precursor_charge, product_charge,
                                  ms_level, run_id, OpenMS::String(mobilogram_type), feature_id, feature_rt, OpenMS::String(filter));

            nb::list run_id_list, source_file_list, ms_level_list;
            nb::list precursor_id_list, transition_id_list, modified_sequence_list;
            nb::list precursor_charge_list, product_charge_list, detecting_transition_list;
            nb::list precursor_decoy_list, product_decoy_list, transition_ordinal_list;
            nb::list transition_type_list, annotation_list, mobility_list, intensity_list;
            nb::list mobilogram_type_list, feature_id_list, feature_rt_list;

            for (const auto& c : mobilos) {
                if (explode) {
                    if (c.mobility.empty()) continue;
                    if (c.mobility.size() != c.intensity.size()) {
                        throw std::runtime_error("XIMParquetFile: mobility/intensity length mismatch");
                    }
                    for (size_t j = 0; j < c.mobility.size(); ++j) {
                        run_id_list.append(c.run_id);
                        source_file_list.append(nb::str(c.source_file.c_str()));
                        ms_level_list.append(c.ms_level);
                        precursor_id_list.append(c.has_precursor_id ? nb::cast(c.precursor_id) : nb::none());
                        transition_id_list.append(c.has_transition_id ? nb::cast(c.transition_id) : nb::none());
                        modified_sequence_list.append(nb::str(c.modified_sequence.c_str()));
                        precursor_charge_list.append(c.has_precursor_charge ? nb::cast(c.precursor_charge) : nb::none());
                        product_charge_list.append(c.has_product_charge ? nb::cast(c.product_charge) : nb::none());
                        detecting_transition_list.append(c.has_detecting_transition ? nb::cast(c.detecting_transition) : nb::none());
                        precursor_decoy_list.append(c.has_precursor_decoy ? nb::cast(c.precursor_decoy) : nb::none());
                        product_decoy_list.append(c.has_product_decoy ? nb::cast(c.product_decoy) : nb::none());
                        transition_ordinal_list.append(c.has_transition_ordinal ? nb::cast(c.transition_ordinal) : nb::none());
                        transition_type_list.append(nb::str(c.transition_type.c_str()));
                        annotation_list.append(nb::str(c.annotation.c_str()));
                        mobilogram_type_list.append(nb::str(c.mobilogram_type.c_str()));
                        feature_id_list.append(c.has_feature_id ? nb::cast(c.feature_id) : nb::none());
                        feature_rt_list.append(c.has_feature_rt ? nb::cast(c.feature_rt) : nb::none());
                        mobility_list.append(c.mobility[j]);
                        intensity_list.append(c.intensity[j]);
                    }
                } else {
                    run_id_list.append(c.run_id);
                    source_file_list.append(nb::str(c.source_file.c_str()));
                    ms_level_list.append(c.ms_level);
                    precursor_id_list.append(c.has_precursor_id ? nb::cast(c.precursor_id) : nb::none());
                    transition_id_list.append(c.has_transition_id ? nb::cast(c.transition_id) : nb::none());
                    modified_sequence_list.append(nb::str(c.modified_sequence.c_str()));
                    precursor_charge_list.append(c.has_precursor_charge ? nb::cast(c.precursor_charge) : nb::none());
                    product_charge_list.append(c.has_product_charge ? nb::cast(c.product_charge) : nb::none());
                    detecting_transition_list.append(c.has_detecting_transition ? nb::cast(c.detecting_transition) : nb::none());
                    precursor_decoy_list.append(c.has_precursor_decoy ? nb::cast(c.precursor_decoy) : nb::none());
                    product_decoy_list.append(c.has_product_decoy ? nb::cast(c.product_decoy) : nb::none());
                    transition_ordinal_list.append(c.has_transition_ordinal ? nb::cast(c.transition_ordinal) : nb::none());
                    transition_type_list.append(nb::str(c.transition_type.c_str()));
                    annotation_list.append(nb::str(c.annotation.c_str()));
                    mobilogram_type_list.append(nb::str(c.mobilogram_type.c_str()));
                    feature_id_list.append(c.has_feature_id ? nb::cast(c.feature_id) : nb::none());
                    feature_rt_list.append(c.has_feature_rt ? nb::cast(c.feature_rt) : nb::none());
                    nb::list mobility_vals, int_vals;
                    for (auto v : c.mobility) mobility_vals.append(v);
                    for (auto v : c.intensity) int_vals.append(v);
                    mobility_list.append(mobility_vals);
                    intensity_list.append(int_vals);
                }
            }

            nb::dict result;
            result["run_id"] = run_id_list;
            result["source_file"] = source_file_list;
            result["ms_level"] = ms_level_list;
            result["precursor_id"] = precursor_id_list;
            result["transition_id"] = transition_id_list;
            result["modified_sequence"] = modified_sequence_list;
            result["precursor_charge"] = precursor_charge_list;
            result["product_charge"] = product_charge_list;
            result["detecting_transition"] = detecting_transition_list;
            result["precursor_decoy"] = precursor_decoy_list;
            result["product_decoy"] = product_decoy_list;
            result["transition_ordinal"] = transition_ordinal_list;
            result["transition_type"] = transition_type_list;
            result["annotation"] = annotation_list;
            result["mobilogram_type"] = mobilogram_type_list;
            result["feature_id"] = feature_id_list;
            result["feature_rt"] = feature_rt_list;
            result["mobility"] = mobility_list;
            result["intensity"] = intensity_list;
            return result;
        }, "precursor_id"_a = -1, "transition_id"_a = -1, "modified_sequence"_a = "",
           "precursor_charge"_a = -1, "product_charge"_a = -1, "ms_level"_a = -1,
           "run_id"_a = -1, "mobilogram_type"_a = "", "feature_id"_a = -1, "feature_rt"_a = -1.0, "filter"_a = "", "explode"_a = false,
           "Return mobilogram data as a dict")
        .def("__repr__", [](const OpenMS::XIMParquetFile& self) {
            const auto& files = self.getFilenames();
            return "XIMParquetFile(n_files=" + std::to_string(files.size()) + ")";
        })
        .def("__str__", [](const OpenMS::XIMParquetFile& self) { return nb::cast(self).attr("__repr__")(); })
        ;

    // -----------------------------------------------------------------------
    // ParquetWriteConfig
    // -----------------------------------------------------------------------
    auto parquetwriteconfig_class = nb::class_<OpenMS::ParquetWriteConfig>(m, "ParquetWriteConfig",
        "Configuration for Parquet file writing (compression, row group size, etc.)")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ParquetWriteConfig&>())
        .def("__copy__", [](const OpenMS::ParquetWriteConfig& self) { return OpenMS::ParquetWriteConfig(self); })
        .def("__deepcopy__", [](const OpenMS::ParquetWriteConfig& self, nb::dict) { return OpenMS::ParquetWriteConfig(self); }, "memo"_a)
        .def_rw("compression", &OpenMS::ParquetWriteConfig::compression)
        .def_rw("compression_level", &OpenMS::ParquetWriteConfig::compression_level)
        .def_rw("row_group_size", &OpenMS::ParquetWriteConfig::row_group_size)
        .def_rw("write_statistics", &OpenMS::ParquetWriteConfig::write_statistics)
        .def_rw("data_page_size", &OpenMS::ParquetWriteConfig::data_page_size)
        ;

    nb::enum_<OpenMS::ParquetWriteConfig::Compression>(parquetwriteconfig_class, "Compression", nb::is_arithmetic())
        .value("NONE", OpenMS::ParquetWriteConfig::Compression::NONE)
        .value("SNAPPY", OpenMS::ParquetWriteConfig::Compression::SNAPPY)
        .value("GZIP", OpenMS::ParquetWriteConfig::Compression::GZIP)
        .value("LZ4", OpenMS::ParquetWriteConfig::Compression::LZ4)
        .value("ZSTD", OpenMS::ParquetWriteConfig::Compression::ZSTD)
        .export_values();

    // -----------------------------------------------------------------------
    // QPXFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::QPXFile>(m, "QPXFile",
        "Export PSM data to Apache Arrow/Parquet format following QPX PSM schema")
        .def(nb::init<>())
        .def_static("exportToParquet",
            static_cast<bool (*)(const std::vector<OpenMS::ProteinIdentification>&,
                                 const OpenMS::PeptideIdentificationList&,
                                 const OpenMS::String&,
                                 bool,
                                 const OpenMS::ParquetWriteConfig&)>(&OpenMS::QPXFile::exportToParquet),
            "protein_identifications"_a, "peptide_identifications"_a,
            "filename"_a, "export_all_psms"_a = false,
            "config"_a = OpenMS::ParquetWriteConfig{},
            "Export PSM data to Parquet file. Returns True on success")
        ;

    // -----------------------------------------------------------------------
    // FeatureMapArrowIO
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::FeatureMapArrowIO>(m, "FeatureMapArrowIO",
        "Import and export FeatureMap data to/from Apache Arrow/Parquet format")
        .def(nb::init<>())
        .def_static("exportToParquet", &OpenMS::FeatureMapArrowIO::exportToParquet,
            "feature_map"_a, "directory"_a,
            "config"_a = OpenMS::ParquetWriteConfig{},
            "Export FeatureMap to a directory of Parquet files. Returns True on success")
        .def_static("importFromParquet", &OpenMS::FeatureMapArrowIO::importFromParquet,
            "directory"_a, "feature_map"_a,
            "Import FeatureMap from a directory of Parquet files. Returns True on success")
        ;

    // -----------------------------------------------------------------------
    // ConsensusMapArrowIO
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ConsensusMapArrowIO>(m, "ConsensusMapArrowIO",
        "Import and export ConsensusMap data to/from Apache Arrow/Parquet format")
        .def(nb::init<>())
        .def_static("exportToParquet", &OpenMS::ConsensusMapArrowIO::exportToParquet,
            "cmap"_a, "directory"_a,
            "config"_a = OpenMS::ParquetWriteConfig{},
            "Export ConsensusMap to a directory of Parquet files. Returns True on success")
        .def_static("importFromParquet", &OpenMS::ConsensusMapArrowIO::importFromParquet,
            "directory"_a, "cmap"_a,
            "Import ConsensusMap from a directory of Parquet files. Returns True on success")
        ;

    // -----------------------------------------------------------------------
    // ProteinIdentificationArrowIO
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ProteinIdentificationArrowIO>(m, "ProteinIdentificationArrowIO",
        "Import and export ProteinIdentification data to/from Apache Arrow/Parquet format")
        .def(nb::init<>())
        .def_static("exportProteinsToParquet", &OpenMS::ProteinIdentificationArrowIO::exportProteinsToParquet,
            "protein_identifications"_a, "filename"_a,
            "config"_a = OpenMS::ParquetWriteConfig{},
            "Export protein hits to Parquet file. Returns True on success")
        .def_static("exportProteinGroupsToParquet", &OpenMS::ProteinIdentificationArrowIO::exportProteinGroupsToParquet,
            "protein_identifications"_a, "filename"_a,
            "config"_a = OpenMS::ParquetWriteConfig{},
            "Export protein groups to Parquet file. Returns True on success")
        .def_static("exportSearchParamsToParquet", &OpenMS::ProteinIdentificationArrowIO::exportSearchParamsToParquet,
            "protein_identifications"_a, "filename"_a,
            "config"_a = OpenMS::ParquetWriteConfig{},
            "Export search parameters to Parquet file. Returns True on success")
        // Import methods: std::vector<ProteinIdentification>& is an output param.
        // Since vectors go through nanobind's STL type caster (creates copies),
        // we must use lambdas that return the modified vector.
        .def_static("importFromParquet", [](const OpenMS::String& proteins_filename,
            const OpenMS::String& protein_groups_filename,
            const OpenMS::String& search_params_filename) {
            std::vector<OpenMS::ProteinIdentification> prot_ids;
            bool ok = OpenMS::ProteinIdentificationArrowIO::importFromParquet(
                proteins_filename, protein_groups_filename, search_params_filename, prot_ids);
            if (!ok) throw std::runtime_error("Failed to import ProteinIdentifications from Parquet");
            return prot_ids;
        }, "proteins_filename"_a, "protein_groups_filename"_a, "search_params_filename"_a,
            "Import all three Parquet files and return reconstructed ProteinIdentifications")
        .def_static("importSearchParamsFromParquet", [](const OpenMS::String& filename) {
            std::vector<OpenMS::ProteinIdentification> prot_ids;
            bool ok = OpenMS::ProteinIdentificationArrowIO::importSearchParamsFromParquet(filename, prot_ids);
            if (!ok) throw std::runtime_error("Failed to import search parameters from Parquet");
            return prot_ids;
        }, "filename"_a,
            "Import search parameters from Parquet file. Returns list of ProteinIdentifications")
        .def_static("importProteinsFromParquet", [](const OpenMS::String& filename,
            std::vector<OpenMS::ProteinIdentification> prot_ids) {
            bool ok = OpenMS::ProteinIdentificationArrowIO::importProteinsFromParquet(filename, prot_ids);
            if (!ok) throw std::runtime_error("Failed to import proteins from Parquet");
            return prot_ids;
        }, "filename"_a, "protein_identifications"_a,
            "Import protein hits from Parquet file into existing ProteinIdentifications. Returns updated list")
        .def_static("importProteinGroupsFromParquet", [](const OpenMS::String& filename,
            std::vector<OpenMS::ProteinIdentification> prot_ids) {
            bool ok = OpenMS::ProteinIdentificationArrowIO::importProteinGroupsFromParquet(filename, prot_ids);
            if (!ok) throw std::runtime_error("Failed to import protein groups from Parquet");
            return prot_ids;
        }, "filename"_a, "protein_identifications"_a,
            "Import protein groups from Parquet file into existing ProteinIdentifications. Returns updated list")
        ;

    // -----------------------------------------------------------------------
    // ControlledVocabulary
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ControlledVocabulary>(m, "ControlledVocabulary",
        "Representation of a controlled vocabulary")
        .def(nb::init<>())
        .def("__copy__", [](const OpenMS::ControlledVocabulary& self) { return OpenMS::ControlledVocabulary(self); })
        .def("__deepcopy__", [](const OpenMS::ControlledVocabulary& self, nb::dict) { return OpenMS::ControlledVocabulary(self); }, "memo"_a)
        .def("loadFromOBO", [](OpenMS::ControlledVocabulary& self, const OpenMS::String& name, const OpenMS::String& filename) { self.loadFromOBO(name, filename); }, "name"_a, "filename"_a, "Loads the CV from an OBO file")
        .def("name", [](const OpenMS::ControlledVocabulary& self) { return self.name(); })
        .def("exists", [](const OpenMS::ControlledVocabulary& self, const OpenMS::String& id) { return self.exists(id); }, "id"_a)
        .def("hasTermWithName", [](const OpenMS::ControlledVocabulary& self, const OpenMS::String& name) { return self.hasTermWithName(name); }, "name"_a)
        .def("isChildOf", [](const OpenMS::ControlledVocabulary& self, const OpenMS::String& child_id, const OpenMS::String& parent_id) { return self.isChildOf(child_id, parent_id); }, "child_id"_a, "parent_id"_a)
        .def("getTerm", [](const OpenMS::ControlledVocabulary& self, const OpenMS::String& id) {
            const auto& term = self.getTerm(id);
            nb::dict d;
            d["id"] = nb::str(term.id.c_str());
            d["name"] = nb::str(term.name.c_str());
            d["description"] = nb::str(term.description.c_str());
            return d;
        }, "id"_a, "Returns the term with the given id as a dict")
        .def("getTermByName", [](const OpenMS::ControlledVocabulary& self, const OpenMS::String& name, const OpenMS::String& desc) {
            const auto& term = self.getTermByName(name, desc);
            nb::dict d;
            d["id"] = nb::str(term.id.c_str());
            d["name"] = nb::str(term.name.c_str());
            d["description"] = nb::str(term.description.c_str());
            return d;
        }, "name"_a, "desc"_a = "", "Returns the term with the given name as a dict")
        .def("getAllChildTerms", [](const OpenMS::ControlledVocabulary& self, const OpenMS::String& parent_id) { std::set<OpenMS::String> terms; self.getAllChildTerms(terms, parent_id); return terms; }, "parent_id"_a, "Returns all child terms of the given parent term")
        ;

    // -----------------------------------------------------------------------
    // CVTerm_ControlledVocabulary (ControlledVocabulary::CVTerm)
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::ControlledVocabulary::CVTerm>(m, "CVTerm_ControlledVocabulary",
        "Representation of a CV term in a controlled vocabulary")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::ControlledVocabulary::CVTerm&>())
        .def("__copy__", [](const OpenMS::ControlledVocabulary::CVTerm& self) { return OpenMS::ControlledVocabulary::CVTerm(self); })
        .def("__deepcopy__", [](const OpenMS::ControlledVocabulary::CVTerm& self, nb::dict) { return OpenMS::ControlledVocabulary::CVTerm(self); }, "memo"_a)
        .def_rw("name", &OpenMS::ControlledVocabulary::CVTerm::name)
        .def_rw("id", &OpenMS::ControlledVocabulary::CVTerm::id)
        .def_rw("parents", &OpenMS::ControlledVocabulary::CVTerm::parents)
        .def_rw("children", &OpenMS::ControlledVocabulary::CVTerm::children)
        .def_rw("obsolete", &OpenMS::ControlledVocabulary::CVTerm::obsolete)
        .def_rw("description", &OpenMS::ControlledVocabulary::CVTerm::description)
        .def_rw("synonyms", &OpenMS::ControlledVocabulary::CVTerm::synonyms)
        .def_rw("unparsed", &OpenMS::ControlledVocabulary::CVTerm::unparsed)
        .def_rw("xref_type", &OpenMS::ControlledVocabulary::CVTerm::xref_type)
        .def_rw("xref_binary", &OpenMS::ControlledVocabulary::CVTerm::xref_binary)
        .def_rw("units", &OpenMS::ControlledVocabulary::CVTerm::units)
        .def_static("getXRefTypeName", &OpenMS::ControlledVocabulary::CVTerm::getXRefTypeName, "type"_a,
            "Returns the name of the cross-reference type")
        .def_static("isHigherBetterScore", &OpenMS::ControlledVocabulary::CVTerm::isHigherBetterScore, "term"_a,
            "Returns whether a higher score is better for this term")
        .def("toXMLString", [](const OpenMS::ControlledVocabulary::CVTerm& self, const OpenMS::String& ref, const OpenMS::String& value) {
            return self.toXMLString(ref, value);
        }, "ref"_a, "value"_a = "", "Returns the XML representation of this term")
        ;
    m.def("__static_CVTerm_ControlledVocabulary_getXRefTypeName", [](OpenMS::ControlledVocabulary::CVTerm::XRefType type) -> OpenMS::String { return OpenMS::ControlledVocabulary::CVTerm::getXRefTypeName(type); }, "type"_a);
    m.def("__static_CVTerm_ControlledVocabulary_isHigherBetterScore", [](OpenMS::ControlledVocabulary::CVTerm term) -> bool { return OpenMS::ControlledVocabulary::CVTerm::isHigherBetterScore(term); }, "term"_a);

    // -----------------------------------------------------------------------
    // XRefType_CVTerm_ControlledVocabulary (ControlledVocabulary::CVTerm::XRefType)
    // -----------------------------------------------------------------------
    nb::enum_<OpenMS::ControlledVocabulary::CVTerm::XRefType>(m, "XRefType_CVTerm_ControlledVocabulary",
        "Cross-reference type for CV terms", nb::is_arithmetic())
        .value("XSD_STRING", OpenMS::ControlledVocabulary::CVTerm::XRefType::XSD_STRING)
        .value("XSD_INTEGER", OpenMS::ControlledVocabulary::CVTerm::XRefType::XSD_INTEGER)
        .value("XSD_DECIMAL", OpenMS::ControlledVocabulary::CVTerm::XRefType::XSD_DECIMAL)
        .value("XSD_NEGATIVE_INTEGER", OpenMS::ControlledVocabulary::CVTerm::XRefType::XSD_NEGATIVE_INTEGER)
        .value("XSD_POSITIVE_INTEGER", OpenMS::ControlledVocabulary::CVTerm::XRefType::XSD_POSITIVE_INTEGER)
        .value("XSD_NON_NEGATIVE_INTEGER", OpenMS::ControlledVocabulary::CVTerm::XRefType::XSD_NON_NEGATIVE_INTEGER)
        .value("XSD_NON_POSITIVE_INTEGER", OpenMS::ControlledVocabulary::CVTerm::XRefType::XSD_NON_POSITIVE_INTEGER)
        .value("XSD_BOOLEAN", OpenMS::ControlledVocabulary::CVTerm::XRefType::XSD_BOOLEAN)
        .value("XSD_DATE", OpenMS::ControlledVocabulary::CVTerm::XRefType::XSD_DATE)
        .value("XSD_ANYURI", OpenMS::ControlledVocabulary::CVTerm::XRefType::XSD_ANYURI)
        .value("NONE", OpenMS::ControlledVocabulary::CVTerm::XRefType::NONE)
        ;

    // -----------------------------------------------------------------------
    // CVMappingFile
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::CVMappingFile>(m, "CVMappingFile",
        "Used to load CvMapping files")
        .def(nb::init<>())
        .def("load", [](OpenMS::CVMappingFile& self, const OpenMS::String& filename, OpenMS::CVMappings& cv_mappings, bool strip_namespaces) { self.load(filename, cv_mappings, strip_namespaces); }, "filename"_a, "cv_mappings"_a, "strip_namespaces"_a = false, "Loads CvMappings from the given file")
        ;

    // -----------------------------------------------------------------------
    // SemanticValidator
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Internal::SemanticValidator>(m, "SemanticValidator",
        "Semantically validates XML files using CVMappings and a ControlledVocabulary")
        .def(nb::init<const OpenMS::CVMappings &, const OpenMS::ControlledVocabulary &>(), "mapping"_a, "cv"_a)
        .def("validate", [](OpenMS::Internal::SemanticValidator& self, const OpenMS::String& filename) {
            OpenMS::StringList errors, warnings;
            bool result = self.validate(filename, errors, warnings);
            return nb::make_tuple(result, errors, warnings);
        }, "filename"_a, "Validates an XML file. Returns (success, errors, warnings)")
        .def("setCheckTermValueTypes", [](OpenMS::Internal::SemanticValidator& self, bool check) { self.setCheckTermValueTypes(check); }, "check"_a)
        .def("setCheckUnits", [](OpenMS::Internal::SemanticValidator& self, bool check) { self.setCheckUnits(check); }, "check"_a)
        .def("setTag", [](OpenMS::Internal::SemanticValidator& self, const OpenMS::String& tag) { self.setTag(tag); }, "tag"_a)
        .def("setAccessionAttribute", [](OpenMS::Internal::SemanticValidator& self, const OpenMS::String& accession) { self.setAccessionAttribute(accession); }, "accession"_a)
        .def("setNameAttribute", [](OpenMS::Internal::SemanticValidator& self, const OpenMS::String& name) { self.setNameAttribute(name); }, "name"_a)
        .def("setValueAttribute", [](OpenMS::Internal::SemanticValidator& self, const OpenMS::String& value) { self.setValueAttribute(value); }, "value"_a)
        .def("setUnitAccessionAttribute", [](OpenMS::Internal::SemanticValidator& self, const OpenMS::String& accession) { self.setUnitAccessionAttribute(accession); }, "accession"_a, "Sets the name of the unit accession attribute")
        .def("setUnitNameAttribute", [](OpenMS::Internal::SemanticValidator& self, const OpenMS::String& name) { self.setUnitNameAttribute(name); }, "name"_a, "Sets the name of the unit name attribute")
        .def("locateTerm", [](const OpenMS::Internal::SemanticValidator& self,
                const OpenMS::String& path,
                const OpenMS::Internal::SemanticValidator::CVTerm& parsed_term) {
            return self.locateTerm(path, parsed_term);
        }, "path"_a, "parsed_term"_a, "Locate a CV term in the mapping file")
        ;


    // -----------------------------------------------------------------------
    // MzMLSqliteHandler
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::Internal::MzMLSqliteHandler>(m, "MzMLSqliteHandler",
        "SQLite-based handler for mzML data (sqMass format)")
        .def("__copy__", [](const OpenMS::Internal::MzMLSqliteHandler& self) { return OpenMS::Internal::MzMLSqliteHandler(self); })
        .def("__deepcopy__", [](const OpenMS::Internal::MzMLSqliteHandler& self, nb::dict) { return OpenMS::Internal::MzMLSqliteHandler(self); }, "memo"_a)
        .def(nb::init<const OpenMS::String&, const OpenMS::UInt64>(), "filename"_a, "run_id"_a)
        .def("readExperiment", [](const OpenMS::Internal::MzMLSqliteHandler& self, OpenMS::MSExperiment& exp, bool meta_only) { self.readExperiment(exp, meta_only); }, "exp"_a, "meta_only"_a = false)
        .def("writeExperiment", [](OpenMS::Internal::MzMLSqliteHandler& self, const OpenMS::MSExperiment& exp) { self.writeExperiment(exp); }, "exp"_a)
        .def("createTables", [](OpenMS::Internal::MzMLSqliteHandler& self) { self.createTables(); })
        .def("writeSpectra", [](OpenMS::Internal::MzMLSqliteHandler& self, const std::vector<OpenMS::MSSpectrum>& spectra) { self.writeSpectra(spectra); }, "spectra"_a)
        .def("writeChromatograms", [](OpenMS::Internal::MzMLSqliteHandler& self, const std::vector<OpenMS::MSChromatogram>& chroms) { self.writeChromatograms(chroms); }, "chromatograms"_a)
        .def("writeRunLevelInformation", [](OpenMS::Internal::MzMLSqliteHandler& self, const OpenMS::MSExperiment& exp, bool write_full_meta) { self.writeRunLevelInformation(exp, write_full_meta); }, "exp"_a, "write_full_meta"_a)
        .def("readSpectra", [](const OpenMS::Internal::MzMLSqliteHandler& self, std::vector<OpenMS::MSSpectrum>& spectra, const std::vector<int>& indices, bool meta_only) { self.readSpectra(spectra, indices, meta_only); }, "spectra"_a, "indices"_a, "meta_only"_a = false)
        .def("readChromatograms", [](const OpenMS::Internal::MzMLSqliteHandler& self, std::vector<OpenMS::MSChromatogram>& chroms, const std::vector<int>& indices, bool meta_only) { self.readChromatograms(chroms, indices, meta_only); }, "chromatograms"_a, "indices"_a, "meta_only"_a = false)
        .def("getNrSpectra", [](const OpenMS::Internal::MzMLSqliteHandler& self) { return self.getNrSpectra(); })
        .def("getNrChromatograms", [](const OpenMS::Internal::MzMLSqliteHandler& self) { return self.getNrChromatograms(); })
        .def("getRunID", [](const OpenMS::Internal::MzMLSqliteHandler& self) { return self.getRunID(); })
        .def("getSpectraIndicesbyRT", [](const OpenMS::Internal::MzMLSqliteHandler& self, double RT, double deltaRT, const std::vector<int>& indices) { return self.getSpectraIndicesbyRT(RT, deltaRT, indices); }, "RT"_a, "deltaRT"_a, "indices"_a)
        .def("setConfig", [](OpenMS::Internal::MzMLSqliteHandler& self, bool write_full_meta, bool use_lossy, double linear_abs, int batch_size) { self.setConfig(write_full_meta, use_lossy, linear_abs, batch_size); }, "write_full_meta"_a, "use_lossy_compression"_a, "linear_abs_mass_acc"_a, "sql_batch_size"_a = 500)
        ;

    // -----------------------------------------------------------------------
    // MSPGenericFile
    // -----------------------------------------------------------------------
    // Note: DefaultParamHandler base not specified because it's in _pyopenms_misc
    // which loads after _pyopenms_format. Methods added directly instead.
    auto mspgenericfile_class = nb::class_<OpenMS::MSPGenericFile>(m, "MSPGenericFile",
        "MSP spectral library file reader/writer")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::String&, OpenMS::MSExperiment&>(), "filename"_a, "library"_a)
        .def("__copy__", [](const OpenMS::MSPGenericFile& self) { return OpenMS::MSPGenericFile(self); })
        .def("__deepcopy__", [](const OpenMS::MSPGenericFile& self, nb::dict) { return OpenMS::MSPGenericFile(self); }, "memo"_a)
        .def("load", [](OpenMS::MSPGenericFile& self, const OpenMS::String& filename, OpenMS::MSExperiment& library) { self.load(filename, library); }, "filename"_a, "library"_a)
        .def("store", [](const OpenMS::MSPGenericFile& self, const OpenMS::String& filename, const OpenMS::MSExperiment& library) { self.store(filename, library); }, "filename"_a, "library"_a)
        .def("getDefaultParameters", [](OpenMS::MSPGenericFile& self) { OpenMS::Param p; self.getDefaultParameters(p); return p; }, "Gets the class' default parameters")
        ;
    def_DefaultParamHandler<OpenMS::MSPGenericFile>(mspgenericfile_class);

    // -----------------------------------------------------------------------
    // SpectrumAccessSqMass
    // -----------------------------------------------------------------------
    nb::class_<OpenMS::SpectrumAccessSqMass>(m, "SpectrumAccessSqMass",
        "SQL-based spectrum access via sqMass files")
        .def(nb::init<const OpenMS::Internal::MzMLSqliteHandler&>(), "handler"_a)
        .def("__copy__", [](const OpenMS::SpectrumAccessSqMass& self) { return OpenMS::SpectrumAccessSqMass(self); })
        .def("__deepcopy__", [](const OpenMS::SpectrumAccessSqMass& self, nb::dict) { return OpenMS::SpectrumAccessSqMass(self); }, "memo"_a)
        .def(nb::init<const OpenMS::Internal::MzMLSqliteHandler&, const std::vector<int>&>(), "handler"_a, "indices"_a)
        .def("getNrSpectra", [](const OpenMS::SpectrumAccessSqMass& self) { return self.getNrSpectra(); })
        .def("getNrChromatograms", [](const OpenMS::SpectrumAccessSqMass& self) { return self.getNrChromatograms(); })
        .def("getSpectrumById", [](OpenMS::SpectrumAccessSqMass& self, int id) { return self.getSpectrumById(id); }, "id"_a)
        .def("getChromatogramById", [](OpenMS::SpectrumAccessSqMass& self, int id) { return self.getChromatogramById(id); }, "id"_a)
        .def("getChromatogramNativeID", [](const OpenMS::SpectrumAccessSqMass& self, int id) { return self.getChromatogramNativeID(id); }, "id"_a, "Returns the native ID of the chromatogram with the given id")
        .def("getSpectraByRT", [](const OpenMS::SpectrumAccessSqMass& self, double RT, double deltaRT) { return self.getSpectraByRT(RT, deltaRT); }, "RT"_a, "deltaRT"_a, "Returns indices of spectra within RT +/- deltaRT")
        ;

    // Free function alias for backward compatibility
    m.def("fromFASTAEntry", [](const OpenMS::FASTAFile::FASTAEntry& fasta) {
        return OpenMS::PEFFEntry::fromFASTAEntry(fasta);
    }, "fasta"_a, "Create a PEFFEntry from a FASTAEntry");

}
