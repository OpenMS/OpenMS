// pyOpenMS nanobind bindings
// Domain: spectrum

#include "all_casters.h"
#include <type_traits>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Product.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>
#include "binding_utils.h"
#include <iomanip>
#include <sstream>

namespace nb = nanobind;
using namespace nb::literals;

// Helper: ensure an object is a C-contiguous numpy array of type T.
template <typename T>
nb::ndarray<nb::numpy, T, nb::ndim<1>, nb::c_contig> as_numpy_array(nb::object obj) {
    nb::ndarray<nb::numpy, T, nb::ndim<1>, nb::c_contig> arr;
    if (nb::try_cast(obj, arr)) return arr;
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

NB_MODULE(_pyopenms_spectrum, m) {
    m.doc() = "pyOpenMS spectrum bindings";

    // ABI guards for zero-copy structured array access (get_peaks_struct dtype depends on these)
    static_assert(std::is_standard_layout_v<OpenMS::Peak1D>,
        "Peak1D must be standard-layout for zero-copy struct views (guarantees member order matches dtype)");
    static_assert(sizeof(OpenMS::Peak1D) == 16,
        "Peak1D must be 16 bytes for zero-copy structured array access");
    static_assert(std::is_same_v<OpenMS::Peak1D::CoordinateType, double>,
        "Peak1D::CoordinateType must be double (dtype assumes float64 for position)");
    static_assert(std::is_same_v<OpenMS::Peak1D::IntensityType, float>,
        "Peak1D::IntensityType must be float (dtype assumes float32 for intensity)");

    // -----------------------------------------------------------------------
    // MSSpectrum
    // -----------------------------------------------------------------------
    auto msspectrum_class = nb::class_<OpenMS::MSSpectrum>(m, "MSSpectrum",
        R"doc(
SpectrumSettings
RangeManagerMzInt

The representation of a 1D spectrum.
Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
Iterations yields access to underlying peak objects but is slower
Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
See help(SpectrumSettings) for information about meta-information
Usage:

.. code-block:: python

  # Access data from an existing spectrum
  ms_level = spectrum.getMSLevel()
  rt = spectrum.getRT()
  mz, intensities = spectrum.get_peaks()

  # Create a new spectrum from scratch
  from pyopenms import *
  spectrum = MSSpectrum()
  spectrum.setDriftTime(25) # 25 ms
  spectrum.setRT(205.2) # 205.2 s
  spectrum.setMSLevel(3) # MS3
  p = Precursor()
  p.setIsolationWindowLowerOffset(1.5)
  p.setIsolationWindowUpperOffset(1.5)
  p.setMZ(600) # isolation at 600 +/- 1.5 Th
  p.setActivationEnergy(40) # 40 eV
  p.setCharge(4) # 4+ ion
  spectrum.setPrecursors( [p] )
  # Add raw data to spectrum
  spectrum.set_peaks( ([401.5], [900]) )
  # Additional data arrays / peak annotations
  fda = FloatDataArray()
  fda.setName("Signal to Noise Array")
  fda.push_back(15)
  sda = StringDataArray()
  sda.setName("Peak annotation")
  sda.push_back("y15++")
  spectrum.setFloatDataArrays( [fda] )
  spectrum.setStringDataArrays( [sda] )
  # Add spectrum to MSExperiment
  exp = MSExperiment()
  exp.addSpectrum(spectrum)
  # Add second spectrum and store as mzML file
  spectrum2 = MSSpectrum()
  spectrum2.set_peaks( ([1, 2], [1, 2]) )
  exp.addSpectrum(spectrum2)
  MzMLFile().store("testfile.mzML", exp)
)doc")
        .def(nb::init<>())
        .def(nb::init<std::initializer_list<OpenMS::Peak1D>>())
        .def(nb::init<const OpenMS::MSSpectrum &>())
        .def("__copy__", [](const OpenMS::MSSpectrum& self) { return OpenMS::MSSpectrum(self); })
        .def("__deepcopy__", [](const OpenMS::MSSpectrum& self, nb::dict) { return OpenMS::MSSpectrum(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("updateRanges", [](OpenMS::MSSpectrum& self) { return self.updateRanges(); }, "Recalculates the m/z and intensity ranges of the spectrum")
        .def("getRT", [](const OpenMS::MSSpectrum& self) { return self.getRT(); }, "Returns the absolute retention time (in seconds)")
        .def("setRT", [](OpenMS::MSSpectrum& self, double rt) { return self.setRT(rt); }, "rt"_a, "Sets the absolute retention time (in seconds)")
        .def("getDriftTime", [](const OpenMS::MSSpectrum& self) { return self.getDriftTime(); }, "Returns the drift time (-1 if not set)")
        .def("setDriftTime", [](OpenMS::MSSpectrum& self, double dt) { return self.setDriftTime(dt); }, "dt"_a, "Sets the drift time (-1 if not set)")
        .def("getDriftTimeUnit", [](const OpenMS::MSSpectrum& self) { return self.getDriftTimeUnit(); }, "Returns the ion mobility drift time unit")
        .def("getDriftTimeUnitAsString", [](const OpenMS::MSSpectrum& self) { return self.getDriftTimeUnitAsString(); }, "Returns the ion mobility drift time unit as string")
        .def("setDriftTimeUnit", [](OpenMS::MSSpectrum& self, OpenMS::DriftTimeUnit dt) { return self.setDriftTimeUnit(dt); }, "dt"_a, "Sets the ion mobility drift time unit")
        .def("getMSLevel", [](const OpenMS::MSSpectrum& self) { return self.getMSLevel(); }, "Returns the MS level")
        .def("setMSLevel", [](OpenMS::MSSpectrum& self, unsigned int ms_level) { return self.setMSLevel(ms_level); }, "ms_level"_a, "Sets the MS level")
        .def("getName", [](const OpenMS::MSSpectrum& self) { return self.getName(); }, "Returns the name of the spectrum")
        .def("setName", [](OpenMS::MSSpectrum& self, const OpenMS::String& name) { return self.setName(name); }, "name"_a, "Sets the name of the spectrum")
        .def("sortByIntensity", [](OpenMS::MSSpectrum& self, bool reverse) { return self.sortByIntensity(reverse); }, "reverse"_a = false, "Sorts the peaks by intensity (ascending if reverse is False, descending if True)")
        .def("sortByPosition", [](OpenMS::MSSpectrum& self) { return self.sortByPosition(); }, "Sorts the peaks by m/z position")
        .def("isSorted", [](const OpenMS::MSSpectrum& self) { return self.isSorted(); }, "Returns True if the spectrum is sorted by m/z")
        .def("findNearest", [](const OpenMS::MSSpectrum& self, double mz) { return self.findNearest(mz); }, "mz"_a, "Returns the index of the closest peak in m/z")
        .def("findNearest", [](const OpenMS::MSSpectrum& self, double mz, double tolerance) { return self.findNearest(mz, tolerance); }, "mz"_a, "tolerance"_a, "Returns the index of the closest peak in the provided +/- m/z tolerance window (-1 if none match)")
        .def("findNearest", [](const OpenMS::MSSpectrum& self, double mz, double tolerance_left, double tolerance_right) { return self.findNearest(mz, tolerance_left, tolerance_right); }, "mz"_a, "tolerance_left"_a, "tolerance_right"_a, "Returns the index of the closest peak in the provided abs. m/z tolerance window to the left and right (-1 if none match)")
        .def("containsIMData", [](const OpenMS::MSSpectrum& self) { return self.containsIMData(); }, "Returns whether the spectrum contains ion mobility data")
        .def("clear", [](OpenMS::MSSpectrum& self, bool clear_meta_data) { return self.clear(clear_meta_data); }, "clear_meta_data"_a, "Clears all data (and meta data if clear_meta_data is True)")
        .def("getType", [](const OpenMS::MSSpectrum& self, bool query_data) { return self.getType(query_data); }, "query_data"_a = false, "Returns the spectrum type (centroided, profile or unknown). If SpectrumSettings and DataProcessing information are not sufficient and query_data is True, the data will be queried (potentially expensive)")
        .def_static("getAllNamesOfSpectrumType", []() { return OpenMS::MSSpectrum::getAllNamesOfSpectrumType(); }, "Returns all spectrum type names known to OpenMS")
        .def_static("spectrumTypeToString", [](OpenMS::SpectrumSettings::SpectrumType type) { return OpenMS::SpectrumSettings::spectrumTypeToString(type); }, "type"_a, "Convert a SpectrumType enum to String. Throws Exception::InvalidValue if type is SIZE_OF_SPECTRUMTYPE")
        .def_static("toSpectrumType", [](const OpenMS::String& name) { return OpenMS::SpectrumSettings::toSpectrumType(name); }, "name"_a, "Convert a string to SpectrumType enum. Throws Exception::InvalidValue if name is not a valid spectrum type")
        .def("unify", [](OpenMS::MSSpectrum& self, const OpenMS::SpectrumSettings& rhs) { return self.unify(rhs); }, "rhs"_a)
        .def("setType", [](OpenMS::MSSpectrum& self, OpenMS::SpectrumSettings::SpectrumType type) { return self.setType(type); }, "type"_a, "Sets the spectrum type")
        .def("setIMFormat", [](OpenMS::MSSpectrum& self, const OpenMS::IMFormat& im_type) { return self.setIMFormat(im_type); }, "im_type"_a, "Sets the ion mobility format")
        .def("getIMFormat", [](const OpenMS::MSSpectrum& self) { return self.getIMFormat(); }, "Returns the ion mobility format")
        .def("setIMPeakType", [](OpenMS::MSSpectrum& self, OpenMS::IMPeakType pt) { self.setIMPeakType(pt); }, "pt"_a, "Sets the IM peak type (profile/centroided)")
        .def("getIMPeakType", [](const OpenMS::MSSpectrum& self) { return self.getIMPeakType(); }, "Returns the IM peak type")
        .def("getNativeID", [](const OpenMS::MSSpectrum& self) { return self.getNativeID(); }, "Returns the native identifier for the spectrum, used by the acquisition software")
        .def("setNativeID", [](OpenMS::MSSpectrum& self, const OpenMS::String& native_id) { return self.setNativeID(native_id); }, "native_id"_a, "Sets the native identifier for the spectrum, used by the acquisition software")
        .def("getComment", [](const OpenMS::MSSpectrum& self) { return self.getComment(); }, "Returns the free-text comment")
        .def("setComment", [](OpenMS::MSSpectrum& self, const OpenMS::String& comment) { return self.setComment(comment); }, "comment"_a, "Sets the free-text comment")
        .def("getInstrumentSettings", [](const OpenMS::MSSpectrum& self) -> const OpenMS::InstrumentSettings & { return self.getInstrumentSettings(); }, nb::rv_policy::reference_internal, "Returns a const reference to the instrument settings of the current spectrum")
        .def("setInstrumentSettings", [](OpenMS::MSSpectrum& self, const OpenMS::InstrumentSettings& instrument_settings) { return self.setInstrumentSettings(instrument_settings); }, "instrument_settings"_a, "Sets the instrument settings of the current spectrum")
        .def("getAcquisitionInfo", [](const OpenMS::MSSpectrum& self) -> const OpenMS::AcquisitionInfo & { return self.getAcquisitionInfo(); }, nb::rv_policy::reference_internal, "Returns a const reference to the acquisition info")
        .def("setAcquisitionInfo", [](OpenMS::MSSpectrum& self, const OpenMS::AcquisitionInfo& acquisition_info) { return self.setAcquisitionInfo(acquisition_info); }, "acquisition_info"_a, "Sets the acquisition info")
        .def("getSourceFile", [](const OpenMS::MSSpectrum& self) -> const OpenMS::SourceFile & { return self.getSourceFile(); }, nb::rv_policy::reference_internal, "Returns a const reference to the source file")
        .def("setSourceFile", [](OpenMS::MSSpectrum& self, const OpenMS::SourceFile& source_file) { return self.setSourceFile(source_file); }, "source_file"_a, "Sets the source file")
        .def("getPrecursors", [](const OpenMS::MSSpectrum& self) -> const std::vector<OpenMS::Precursor> & { return self.getPrecursors(); }, nb::rv_policy::reference_internal, "Returns a const reference to the precursors")
        .def("setPrecursors", [](OpenMS::MSSpectrum& self, const std::vector<OpenMS::Precursor>& precursors) { return self.setPrecursors(precursors); }, "precursors"_a, "Sets the precursors")
        .def("getProducts", [](const OpenMS::MSSpectrum& self) -> const std::vector<OpenMS::Product> & { return self.getProducts(); }, nb::rv_policy::reference_internal, "Returns a const reference to the products")
        .def("setProducts", [](OpenMS::MSSpectrum& self, const std::vector<OpenMS::Product>& products) { return self.setProducts(products); }, "products"_a, "Sets the products")
        .def("setDataProcessing", [](OpenMS::MSSpectrum& self, const std::vector<std::shared_ptr<OpenMS::DataProcessing>>& data_processing) { return self.setDataProcessing(data_processing); }, "data_processing"_a)
        .def("getDataProcessing", [](OpenMS::MSSpectrum& self) -> std::vector<std::shared_ptr<OpenMS::DataProcessing>> & { return self.getDataProcessing(); }, nb::rv_policy::reference_internal)

        .def("__iter__", [](OpenMS::MSSpectrum& self) { return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::MSSpectrum>(), "MSSpectrum_iter", self.begin(), self.end()); })
        .def("__len__", [](OpenMS::MSSpectrum& self) { return self.size(); })

        .def("getIMData", [](const OpenMS::MSSpectrum& self) {
            auto result = self.getIMData();
            return nb::make_tuple((int)result.first, (int)result.second);
        }, "Returns (index, drift_time_unit) for ion mobility data")

        .def("_get_peaks_view", [](nb::object self_obj) {
            auto& self = nb::cast<OpenMS::MSSpectrum&>(self_obj);
            // Cast to a raw byte pointer
            uint8_t* data_ptr = self.empty() ? nullptr : reinterpret_cast<uint8_t*>(self.data());

            // Shape is total number of peaks * size of one peak (16 bytes)
            size_t shape[1] = { self.size() * sizeof(OpenMS::Peak1D) };

            // Return as a 1D NumPy array of unsigned 8-bit integers (bytes)
            return nb::ndarray<nb::numpy, uint8_t, nb::c_contig>(
                data_ptr,
                1,
                shape,
                self_obj
            );
        },
        nb::rv_policy::reference_internal,
        "Returns a raw byte view of the underlying Peak1D array (AoS layout).")

        .def("get_peaks_struct",
            [](nb::object self_obj) -> nb::object {
                auto& self = nb::cast<OpenMS::MSSpectrum&>(self_obj);
                size_t n = self.size();
                auto np = nb::module_::import_("numpy");
                nb::dict dtype_dict;
                // Derive dtype from C++ layout (validated by static_asserts at module init)
                constexpr size_t pos_offset = 0; // standard-layout: first member at offset 0
                constexpr size_t int_offset = sizeof(OpenMS::Peak1D::PositionType);
                dtype_dict["names"] = nb::make_tuple("mz", "intensity");
                dtype_dict["formats"] = nb::make_tuple(np.attr("float64"), np.attr("float32"));
                dtype_dict["offsets"] = nb::make_tuple(pos_offset, int_offset);
                dtype_dict["itemsize"] = sizeof(OpenMS::Peak1D);
                auto py_dtype = np.attr("dtype")(dtype_dict);
                if (n == 0) {
                    return np.attr("empty")(0, py_dtype);
                }
                uint8_t* data_ptr = reinterpret_cast<uint8_t*>(self.data());
                size_t byte_shape[1] = { n * sizeof(OpenMS::Peak1D) };
                auto raw = nb::ndarray<nb::numpy, uint8_t, nb::c_contig>(
                    data_ptr, 1, byte_shape, self_obj
                );
                return np.attr("frombuffer")(raw, py_dtype);
            },
            nb::rv_policy::reference_internal,
            "Returns zero-copy structured array with fields 'mz' (float64) and 'intensity' (float32)."
        )

        .def("get_peaks", [](const OpenMS::MSSpectrum& self) {
            // Return (mz_array, intensity_array) as numpy arrays
            // mz as float64 (double), intensity as float32 (float) matching C++ storage
            // Single allocation + single capsule to reduce overhead for small arrays
            const size_t n = self.size();
            const size_t mz_bytes = n * sizeof(double);
            const size_t int_bytes = n * sizeof(float);
            char* buf = new char[mz_bytes + int_bytes];
            double* mz_data = reinterpret_cast<double*>(buf);
            float* int_data = reinterpret_cast<float*>(buf + mz_bytes);
            for (size_t i = 0; i < n; ++i) {
                mz_data[i] = self[i].getMZ();
                int_data[i] = self[i].getIntensity();
            }
            nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<char*>(p); });
            auto mz_arr = nb::ndarray<nb::numpy, double, nb::ndim<1>>(mz_data, {n}, owner);
            auto int_arr = nb::ndarray<nb::numpy, float, nb::ndim<1>>(int_data, {n}, owner);
            return nb::make_tuple(mz_arr, int_arr);
        }, "Returns a tuple of (mz_array, intensity_array) as numpy arrays")

        .def("set_peaks", [](OpenMS::MSSpectrum& self, nb::object mz_obj, nb::object int_obj) {
            // Fast path: direct pointer access from numpy arrays (no intermediate vector copy)
            // mz is double, intensity is float matching Peak1D storage
            auto mz_arr = as_numpy_array<double>(mz_obj);
            auto int_arr = as_numpy_array<float>(int_obj);
            const size_t n = mz_arr.shape(0);
            if (int_arr.shape(0) != n) {
                throw std::runtime_error("mz and intensity arrays must have same length");
            }
            self.resize(n);
            const double* mz_ptr = static_cast<const double*>(mz_arr.data());
            const float* int_ptr = static_cast<const float*>(int_arr.data());
            for (size_t i = 0; i < n; ++i) {
                self[i].setMZ(mz_ptr[i]);
                self[i].setIntensity(int_ptr[i]);
            }
        }, "mz"_a, "intensity"_a, "Set peaks from mz and intensity arrays")
        .def("set_peaks", [](OpenMS::MSSpectrum& self, nb::object peaks_seq) {
            if (nb::len(peaks_seq) != 2) {
                throw std::runtime_error("set_peaks sequence must contain exactly 2 arrays (mz, intensity)");
            }
            auto mz_arr = as_numpy_array<double>(peaks_seq[0]);
            auto int_arr = as_numpy_array<float>(peaks_seq[1]);
            const size_t n = mz_arr.shape(0);
            if (int_arr.shape(0) != n) {
                throw std::runtime_error("mz and intensity arrays must have same length");
            }
            self.resize(n);
            const double* mz_ptr = static_cast<const double*>(mz_arr.data());
            const float* int_ptr = static_cast<const float*>(int_arr.data());
            for (size_t i = 0; i < n; ++i) {
                self[i].setMZ(mz_ptr[i]);
                self[i].setIntensity(int_ptr[i]);
            }
        }, "peaks"_a, "Set peaks from a tuple of (mz_array, intensity_array)")

        .def("push_back", [](OpenMS::MSSpectrum& self, const OpenMS::Peak1D& p) {
            self.push_back(p);
        }, "peak"_a, "Add a peak to the spectrum")

        .def("size", [](const OpenMS::MSSpectrum& self) {
            return self.size();
        }, "Returns the number of peaks")

        .def("__getitem__", [](const OpenMS::MSSpectrum& self, size_t i) {
            if (i >= self.size()) throw nb::index_error();
            return self[i];  // Return by value (copy)
        }, "i"_a, "Returns a copy of the peak at index i")
        .def("__setitem__", [](OpenMS::MSSpectrum& self, size_t i, const OpenMS::Peak1D& val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        }, "i"_a, "val"_a, "Sets peak at index i")
        .def("findHighestInWindow", [](const OpenMS::MSSpectrum& self, double mz, double tolerance_left, double tolerance_right) { return self.findHighestInWindow(mz, tolerance_left, tolerance_right); }, "mz"_a, "tolerance_left"_a, "tolerance_right"_a, "Returns the index of the highest peak in the provided abs. m/z tolerance window (-1 if none match)")
        .def("select", [](OpenMS::MSSpectrum& self, const std::vector<size_t>& indices) -> OpenMS::MSSpectrum& { return self.select(indices); }, nb::rv_policy::reference_internal, "indices"_a, "Selects peaks by indices, removing all others")

        .def("getMinMZ", &OpenMS::MSSpectrum::getMinMZ, "Returns minimum m/z value")

        .def("getMaxMZ", &OpenMS::MSSpectrum::getMaxMZ, "Returns maximum m/z value")

        .def("getMinIntensity", &OpenMS::MSSpectrum::getMinIntensity, "Returns minimum intensity value")

        .def("getMaxIntensity", &OpenMS::MSSpectrum::getMaxIntensity, "Returns maximum intensity value")

        .def("clearRanges", &OpenMS::MSSpectrum::clearRanges, "Resets all range dimensions")

        .def("get_mz_array", [](const OpenMS::MSSpectrum& self) {
            size_t n = self.size();
            double* data = new double[n];
            for (size_t i = 0; i < n; ++i) data[i] = self[i].getMZ();
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            return nb::ndarray<nb::numpy, double, nb::ndim<1>>(data, {n}, owner);
        }, "Returns m/z values as numpy array")

        .def("get_intensity_array", [](const OpenMS::MSSpectrum& self) {
            size_t n = self.size();
            float* data = new float[n];
            for (size_t i = 0; i < n; ++i) data[i] = self[i].getIntensity();
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(data, {n}, owner);
        }, "Returns intensity values as numpy array")

        .def("get_drift_time_array", [](const OpenMS::MSSpectrum& self) -> std::optional<nb::ndarray<nb::numpy, float, nb::ndim<1>>> {
            if (!self.containsIMData()) return std::nullopt;
            auto [im_index, im_unit] = self.getIMData();
            const auto& arr = self.getFloatDataArrays()[im_index];
            size_t n = arr.size();
            float* data = new float[n];
            std::copy(arr.begin(), arr.end(), data);
            nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(data, {n}, owner);
        }, "Returns drift time array if ion mobility data exists, else None")

        .def("get_drift_time_array_view", [](nb::object self_obj) -> std::optional<nb::ndarray<nb::numpy, float, nb::ndim<1>>> {
            // Writable zero-copy view into the IM float data array
            auto& self = nb::cast<OpenMS::MSSpectrum&>(self_obj);
            if (!self.containsIMData()) return std::nullopt;
            auto [im_index, im_unit] = self.getIMData();
            auto& arr = self.getFloatDataArrays()[im_index];
            float* data_ptr = arr.empty() ? nullptr : arr.data();
            return nb::ndarray<nb::numpy, float, nb::ndim<1>>(
                data_ptr, {arr.size()}, self_obj
            );
        }, "Returns writable view of drift time array if ion mobility data exists, else None")

        .def("getFloatDataArrays", [](OpenMS::MSSpectrum& self) -> std::vector<OpenMS::DataArrays::FloatDataArray>& {
            return self.getFloatDataArrays();
        }, nb::rv_policy::reference_internal, "Returns the float data arrays")

        .def("setFloatDataArrays", [](OpenMS::MSSpectrum& self, const std::vector<OpenMS::DataArrays::FloatDataArray>& arrays) {
            self.setFloatDataArrays(arrays);
        }, "arrays"_a, "Set the float data arrays")

        .def("getIntegerDataArrays", [](OpenMS::MSSpectrum& self) -> std::vector<OpenMS::DataArrays::IntegerDataArray>& {
            return self.getIntegerDataArrays();
        }, nb::rv_policy::reference_internal, "Returns the integer data arrays")

        .def("setIntegerDataArrays", [](OpenMS::MSSpectrum& self, const std::vector<OpenMS::DataArrays::IntegerDataArray>& arrays) {
            self.setIntegerDataArrays(arrays);
        }, "arrays"_a, "Set the integer data arrays")

        .def("getStringDataArrays", [](OpenMS::MSSpectrum& self) -> std::vector<OpenMS::DataArrays::StringDataArray>& {
            return self.getStringDataArrays();
        }, nb::rv_policy::reference_internal, "Returns the string data arrays")

        .def("setStringDataArrays", [](OpenMS::MSSpectrum& self, const std::vector<OpenMS::DataArrays::StringDataArray>& arrays) {
            self.setStringDataArrays(arrays);
        }, "arrays"_a, "Set the string data arrays")

        .def("get_drift_time_unit", [](const OpenMS::MSSpectrum& self) -> std::optional<OpenMS::DriftTimeUnit> {
            if (!self.containsIMData()) return std::nullopt;
            return self.getDriftTimeUnit();
        }, "Returns drift time unit if ion mobility data exists, else None")

        .def("calculateTIC", [](const OpenMS::MSSpectrum& self) -> double {
            return static_cast<double>(self.calculateTIC());
        }, "Returns the total ion current (sum of all peak intensities)")

        .def("reserve", [](OpenMS::MSSpectrum& self, size_t n) {
            self.reserve(n);
        }, "n"_a, "Reserves space for n peaks in the underlying container")

        .def("resize", [](OpenMS::MSSpectrum& self, size_t n) {
            self.resize(n);
        }, "n"_a, "Resizes the spectrum to contain n peaks")

        .def("rasterizeIMFrame", [](const OpenMS::MSSpectrum& self,
                nb::ndarray<float, nb::ndim<2>, nb::device::cpu> output,
                double min_im, double max_im,
                double min_mz, double max_mz,
                const std::string& aggregation) {
            // Check for C-contiguous array (required for correct pointer arithmetic)
            if (output.stride(1) != 1 ||
                output.stride(0) != static_cast<int64_t>(output.shape(1))) {
                throw std::invalid_argument("Output array must be C-contiguous. "
                    "Use numpy.ascontiguousarray() to convert.");
            }
            const size_t mz_bins = output.shape(0);
            const size_t im_bins = output.shape(1);
            float* output_ptr = output.data();
            OpenMS::MSSpectrum::RasterAggregation agg_mode;
            if (aggregation == "sum" || aggregation == "SUM")
                agg_mode = OpenMS::MSSpectrum::RasterAggregation::SUM;
            else if (aggregation == "max" || aggregation == "MAX")
                agg_mode = OpenMS::MSSpectrum::RasterAggregation::MAX;
            else
                throw std::invalid_argument("Invalid aggregation mode '" + aggregation + "'. Must be 'sum' or 'max'.");
            self.rasterizeIMFrame(output_ptr, im_bins, mz_bins, min_im, max_im, min_mz, max_mz, agg_mode);
        }, "output"_a, "min_im"_a, "max_im"_a, "min_mz"_a, "max_mz"_a, "aggregation"_a = "sum",
        R"doc(Rasterize an ion mobility frame into a 2D intensity matrix (IM vs m/z).

Creates a 2D heatmap representation by binning peak intensities into a regular grid.
Designed for spectra in IM_PEAK format where each peak has an associated ion mobility value.

Parameters
----------
output : numpy.ndarray
    Pre-allocated 2D float32 C-contiguous array with shape (mz_bins, im_bins).
    Will be filled in-place (zero-initialize before calling).
min_im : float
    Minimum ion mobility value for the output range.
max_im : float
    Maximum ion mobility value for the output range.
min_mz : float
    Minimum m/z value for the output range.
max_mz : float
    Maximum m/z value for the output range.
aggregation : str, optional
    Aggregation mode: "sum" (default) or "max".
)doc")
        .def("__repr__", [](const OpenMS::MSSpectrum& self) {
            std::ostringstream oss;
            oss << "MSSpectrum(n_peaks=" << self.size()
                << ", rt=" << std::fixed << std::setprecision(2) << self.getRT() << "s"
                << ", ms_level=" << self.getMSLevel() << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::MSSpectrum& self) { return nb::cast(self).attr("__repr__")(); })
        ;
    def_MetaInfoInterface<OpenMS::MSSpectrum>(msspectrum_class);
    // MSSpectrumRasterAggregation enum nested under MSSpectrum
    nb::enum_<OpenMS::MSSpectrum::RasterAggregation>(msspectrum_class, "MSSpectrumRasterAggregation")
        .value("SUM", OpenMS::MSSpectrum::RasterAggregation::SUM)
        .value("MAX", OpenMS::MSSpectrum::RasterAggregation::MAX)
        ;

} // NB_MODULE
