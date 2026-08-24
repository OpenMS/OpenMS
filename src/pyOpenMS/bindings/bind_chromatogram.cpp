// pyOpenMS nanobind bindings
// Domain: chromatogram

#include "all_casters.h"
#include <type_traits>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/METADATA/ChromatogramSettings.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/Product.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>
#include "binding_utils.h"
#include "peak_layout.h"
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

// ABI guards for zero-copy structured array access (peaks_struct dtype depends on these).
// The static_assert is what instantiates PeakLayout, so the guards run even if the dtype below
// is refactored away; it also restates the offsets peaks_struct publishes to Python.
using ChromatogramPeakLayout = pyopenms::PeakLayout<OpenMS::ChromatogramPeak>;
static_assert(ChromatogramPeakLayout::position_offset == 0 && ChromatogramPeakLayout::intensity_offset == 8,
    "ChromatogramPeak's structured dtype is documented as rt (float64) at 0, intensity (float32) at 8");

NB_MODULE(_pyopenms_chromatogram, m) {
    m.doc() = "pyOpenMS chromatogram bindings";

    // -----------------------------------------------------------------------
    // MSChromatogram
    // -----------------------------------------------------------------------
    auto mschromatogram_class = nb::class_<OpenMS::MSChromatogram>(m, "MSChromatogram",
        R"doc(
ChromatogramSettings
RangeManagerRtInt

The representation of a chromatogram.
Raw data access is proved by `get_peaks` and `set_peaks`, which yields numpy arrays
Indexing and iteration yield copies of the peaks; write changes back with chrom[i] = peak
Extra data arrays can be accessed through getFloatDataArrays / getIntegerDataArrays / getStringDataArrays
See help(ChromatogramSettings) for information about meta-information
Usage:
.. code-block:: python
precursor = chromatogram.getPrecursor()
product = chromatogram.getProduct()
rt, intensities = chromatogram.get_peaks()
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MSChromatogram &>())
        .def("__copy__", [](const OpenMS::MSChromatogram& self) { return OpenMS::MSChromatogram(self); })
        .def("__deepcopy__", [](const OpenMS::MSChromatogram& self, nb::dict) { return OpenMS::MSChromatogram(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("updateRanges", [](OpenMS::MSChromatogram& self) { return self.updateRanges(); }, "Recalculates the RT and intensity ranges of the chromatogram")
        .def("clearRanges", [](OpenMS::MSChromatogram& self) { return self.clearRanges(); }, "Clear all ranges")
        .def("resize", [](OpenMS::MSChromatogram& self, size_t new_size) { return self.resize(new_size); }, "new_size"_a, "Resize the peak array")
        .def("reserve", [](OpenMS::MSChromatogram& self, size_t new_size) { return self.reserve(new_size); }, "new_size"_a, "Reserve space for peaks")
        .def("getMinRT", [](const OpenMS::MSChromatogram& self) { return self.getMinRT(); }, "Get the minimum RT value")
        .def("getMaxRT", [](const OpenMS::MSChromatogram& self) { return self.getMaxRT(); }, "Get the maximum RT value")
        .def("getMinIntensity", [](const OpenMS::MSChromatogram& self) { return self.getMinIntensity(); }, "Get the minimum intensity value")
        .def("getMaxIntensity", [](const OpenMS::MSChromatogram& self) { return self.getMaxIntensity(); }, "Get the maximum intensity value")
        .def("getName", [](const OpenMS::MSChromatogram& self) { return self.getName(); }, "Returns the name")
        .def("setName", [](OpenMS::MSChromatogram& self, const std::string& name) { return self.setName(name); }, "name"_a, "Sets the name")
        .def("getMZ", [](const OpenMS::MSChromatogram& self) { return self.getMZ(); }, "Returns the mz of the product entry, makes sense especially for MRM scans")
        .def("sortByIntensity", [](OpenMS::MSChromatogram& self, bool reverse) { return self.sortByIntensity(reverse); }, "reverse"_a = false,
            R"doc(
Lexicographically sorts the peaks by their intensity.
Sorts the peaks according to ascending intensity (lowest to highest) when reverse is False.
When reverse is True, sorts by descending intensity (highest to lowest).
Meta data arrays will be sorted accordingly
)doc")
        .def("sortByPosition", [](OpenMS::MSChromatogram& self) { return self.sortByPosition(); },
            R"doc(
Lexicographically sorts the peaks by their position (RT).
Sorts the peaks according to ascending RT. Meta data arrays will be sorted accordingly
)doc")
        .def("isSorted", [](const OpenMS::MSChromatogram& self) { return self.isSorted(); }, "Checks if all peaks are sorted with respect to ascending RT")
        .def("findNearest", [](const OpenMS::MSChromatogram& self, double rt) { return self.findNearest(rt); }, "rt"_a,
            R"doc(
Binary search for the peak nearest to a specific RT.
Returns the index of the peak. The chromatogram must be sorted with respect to RT, otherwise
the result is undefined. Raises an exception if the chromatogram is empty.
)doc")
        .def("clear", [](OpenMS::MSChromatogram& self, bool clear_meta_data) { return self.clear(clear_meta_data); }, "clear_meta_data"_a,
            R"doc(
Clears all data and meta data.
Always deletes all peaks, associated data arrays (float, integer, string), and ranges.
If clear_meta_data is True, also deletes the descriptive meta data (ChromatogramSettings, name).
)doc")
        .def("getNativeID", [](const OpenMS::MSChromatogram& self) { return self.getNativeID(); }, "Returns the native identifier for the spectrum, used by the acquisition software.")
        .def("setNativeID", [](OpenMS::MSChromatogram& self, const std::string& native_id) { return self.setNativeID(native_id); }, "native_id"_a, "Sets the native identifier for the spectrum, used by the acquisition software.")
        .def("getComment", [](const OpenMS::MSChromatogram& self) { return self.getComment(); }, "Returns the free-text comment")
        .def("setComment", [](OpenMS::MSChromatogram& self, const std::string& comment) { return self.setComment(comment); }, "comment"_a, "Sets the free-text comment")
        .def("getInstrumentSettings", [](const OpenMS::MSChromatogram& self) -> OpenMS::InstrumentSettings { return self.getInstrumentSettings(); }, "Returns the instrument settings of the current spectrum")
        .def("setInstrumentSettings", [](OpenMS::MSChromatogram& self, const OpenMS::InstrumentSettings& instrument_settings) { return self.setInstrumentSettings(instrument_settings); }, "instrument_settings"_a, "Sets the instrument settings of the current spectrum")
        .def("getAcquisitionInfo", [](const OpenMS::MSChromatogram& self) -> OpenMS::AcquisitionInfo { return self.getAcquisitionInfo(); }, "Returns the acquisition info")
        .def("setAcquisitionInfo", [](OpenMS::MSChromatogram& self, const OpenMS::AcquisitionInfo& acquisition_info) { return self.setAcquisitionInfo(acquisition_info); }, "acquisition_info"_a, "Sets the acquisition info")
        .def("getSourceFile", [](const OpenMS::MSChromatogram& self) -> OpenMS::SourceFile { return self.getSourceFile(); }, "Returns the source file")
        .def("setSourceFile", [](OpenMS::MSChromatogram& self, const OpenMS::SourceFile& source_file) { return self.setSourceFile(source_file); }, "source_file"_a, "Sets the source file")
        .def("getPrecursor", [](const OpenMS::MSChromatogram& self) -> OpenMS::Precursor { return self.getPrecursor(); }, "Returns the precursors")
        .def("setPrecursor", [](OpenMS::MSChromatogram& self, const OpenMS::Precursor& precursor) { return self.setPrecursor(precursor); }, "precursor"_a, "Sets the precursors")
        .def("getProduct", [](const OpenMS::MSChromatogram& self) -> OpenMS::Product { return self.getProduct(); }, "Returns the product ion")
        .def("setProduct", [](OpenMS::MSChromatogram& self, const OpenMS::Product& product) { return self.setProduct(product); }, "product"_a, "Sets the product ion")
        .def("getChromatogramType", [](const OpenMS::MSChromatogram& self) { return self.getChromatogramType(); }, "Get the chromatogram type")
        .def("setChromatogramType", [](OpenMS::MSChromatogram& self, OpenMS::ChromatogramSettings::ChromatogramType type) { return self.setChromatogramType(type); }, "type"_a, "Sets the chromatogram type")
        .def("setDataProcessing", [](OpenMS::MSChromatogram& self, const std::vector<std::shared_ptr<OpenMS::DataProcessing>>& data_processing) { return self.setDataProcessing(data_processing); }, "data_processing"_a, "Sets the description of the applied processing")
        .def("getDataProcessing", [](OpenMS::MSChromatogram& self) -> std::vector<std::shared_ptr<OpenMS::DataProcessing>> { return self.getDataProcessing(); }, "Returns the description of the applied processing")

        .def("__iter__", [](OpenMS::MSChromatogram& self) { return nb::make_iterator<nb::rv_policy::copy>(nb::type<OpenMS::MSChromatogram>(), "MSChromatogram_iter", self.begin(), self.end()); }, nb::keep_alive<0, 1>())
        .def("__len__", [](OpenMS::MSChromatogram& self) { return self.size(); })
        .def("__getitem__", [](const OpenMS::MSChromatogram& self, size_t i) -> OpenMS::ChromatogramPeak {
            if (i >= self.size()) throw nb::index_error();
            return self[i];  // by value: element access yields an owned copy
        }, "i"_a, "Returns a copy of the peak at index i")
        .def("__setitem__", [](OpenMS::MSChromatogram& self, size_t i, const OpenMS::ChromatogramPeak& val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        }, "i"_a, "val"_a, "Sets peak at index i")

        .def("get_peaks", [](const OpenMS::MSChromatogram& self) {
            const size_t n = self.size();
            const size_t rt_bytes = n * sizeof(double);
            const size_t int_bytes = n * sizeof(float);
            char* buf = new char[rt_bytes + int_bytes];
            double* rt_data = reinterpret_cast<double*>(buf);
            float* int_data = reinterpret_cast<float*>(buf + rt_bytes);
            for (size_t i = 0; i < n; ++i) {
                rt_data[i] = self[i].getRT();
                int_data[i] = self[i].getIntensity();
            }
            nb::capsule owner(buf, [](void* p) noexcept { delete[] static_cast<char*>(p); });
            auto rt_arr = nb::ndarray<nb::numpy, double, nb::ndim<1>>(rt_data, {n}, owner);
            auto int_arr = nb::ndarray<nb::numpy, float, nb::ndim<1>>(int_data, {n}, owner);
            return nb::make_tuple(rt_arr, int_arr);
        }, "Returns a tuple of (rt_array, intensity_array) as numpy arrays")
    
    
    
        .def("_get_peaks_view", [](nb::object self_obj) {
            auto& self = nb::cast<OpenMS::MSChromatogram&>(self_obj);
            uint8_t* data_ptr = self.empty() ? nullptr : reinterpret_cast<uint8_t*>(&self[0]);
            size_t shape[1] = { self.size() * sizeof(OpenMS::ChromatogramPeak) };
            return nb::ndarray<nb::numpy, uint8_t, nb::c_contig>(
                data_ptr,
                1,
                shape,
                self_obj
            );
        },
        nb::rv_policy::reference_internal,
        "Returns a raw byte view of the underlying ChromatogramPeak array (AoS layout).")

        .def("peaks_struct",
            [](nb::object self_obj) -> nb::object {
                auto& self = nb::cast<OpenMS::MSChromatogram&>(self_obj);
                size_t n = self.size();
                auto np = nb::module_::import_("numpy");

                // Offsets are read off the C++ layout, not reconstructed from it
                nb::dict dtype_dict;
                dtype_dict["names"] = nb::make_tuple("rt", "intensity");
                dtype_dict["formats"] = nb::make_tuple(np.attr("float64"), np.attr("float32"));
                dtype_dict["offsets"] = nb::make_tuple(ChromatogramPeakLayout::position_offset, ChromatogramPeakLayout::intensity_offset);
                dtype_dict["itemsize"] = ChromatogramPeakLayout::itemsize;
                auto py_dtype = np.attr("dtype")(dtype_dict);

                if (n == 0) {
                    return np.attr("empty")(0, py_dtype);
                }

                uint8_t* data_ptr = reinterpret_cast<uint8_t*>(&self[0]);
                size_t byte_shape[1] = { n * sizeof(OpenMS::ChromatogramPeak) };
                auto raw = nb::ndarray<nb::numpy, uint8_t, nb::c_contig>(
                    data_ptr,
                    1,
                    byte_shape,
                    self_obj
                );

                return np.attr("frombuffer")(raw, py_dtype);
            },
            nb::rv_policy::reference_internal,
            "Returns zero-copy structured array with fields 'rt' (float64) and 'intensity' (float32)."
        )
    
        .def("set_peaks", [](OpenMS::MSChromatogram& self, nb::object rt_obj, nb::object int_obj, const std::string& metadata) {
            // set_peaks(peaks, "clear") binds here rather than to the sequence overload below,
            // because this one is registered first and a 2-tuple is a perfectly good nb::object.
            // Without this it would fail deep inside as_numpy_array with an unrelated message.
            if (nb::isinstance<nb::str>(int_obj) || nb::isinstance<nb::bytes>(int_obj)) {
                throw std::invalid_argument("set_peaks() received a string as the intensity argument. "
                                            "Did you mean set_peaks(peaks, metadata=...)?");
            }
            const PeakMetadataPolicy policy = parsePeakMetadataPolicy(metadata);
            auto rt_arr = as_numpy_array<double>(rt_obj);
            auto int_arr = as_numpy_array<float>(int_obj);
            const size_t n = rt_arr.shape(0);
            if (int_arr.shape(0) != n) {
                throw std::runtime_error("rt and intensity arrays must have same length");
            }
            const double* rt_ptr = static_cast<const double*>(rt_arr.data());
            const float* int_ptr = static_cast<const float*>(int_arr.data());
            writePeaksWithPolicy(self, n, policy, "chromatogram",
                                 [rt_ptr, int_ptr](OpenMS::MSChromatogram& chrom, size_t count) {
                                     for (size_t i = 0; i < count; ++i) {
                                         chrom[i].setRT(rt_ptr[i]);
                                         chrom[i].setIntensity(int_ptr[i]);
                                     }
                                 });
        }, "rt"_a, "intensity"_a, "metadata"_a = "error",
           "Set peaks from rt and intensity arrays" PYOPENMS_SET_PEAKS_METADATA_DOC)
        .def("set_peaks", [](OpenMS::MSChromatogram& self, nb::object peaks_seq, const std::string& metadata) {
            const PeakMetadataPolicy policy = parsePeakMetadataPolicy(metadata);
            if (nb::len(peaks_seq) != 2) {
                throw std::runtime_error("set_peaks sequence must contain exactly 2 arrays (rt, intensity)");
            }
            auto rt_arr = as_numpy_array<double>(peaks_seq[0]);
            auto int_arr = as_numpy_array<float>(peaks_seq[1]);
            const size_t n = rt_arr.shape(0);
            if (int_arr.shape(0) != n) {
                throw std::runtime_error("rt and intensity arrays must have same length");
            }
            const double* rt_ptr = static_cast<const double*>(rt_arr.data());
            const float* int_ptr = static_cast<const float*>(int_arr.data());
            writePeaksWithPolicy(self, n, policy, "chromatogram",
                                 [rt_ptr, int_ptr](OpenMS::MSChromatogram& chrom, size_t count) {
                                     for (size_t i = 0; i < count; ++i) {
                                         chrom[i].setRT(rt_ptr[i]);
                                         chrom[i].setIntensity(int_ptr[i]);
                                     }
                                 });
        }, "peaks"_a, nb::kw_only(), "metadata"_a = "error",
           "Set peaks from a tuple/list of (rt_array, intensity_array)" PYOPENMS_SET_PEAKS_METADATA_DOC)

        .def("size", [](const OpenMS::MSChromatogram& self) {
            return self.size();
        }, "Returns the number of peaks")

        .def("push_back", [](OpenMS::MSChromatogram& self, const OpenMS::ChromatogramPeak& peak) {
            self.push_back(peak);
        }, "peak"_a, "Append a peak")

        .def("_float_data_array_count", [](const OpenMS::MSChromatogram& self) { return self.getFloatDataArrays().size(); })
        .def("float_data_array_view", [](OpenMS::MSChromatogram& self, size_t i) -> OpenMS::DataArrays::FloatDataArray& {
            if (i >= self.getFloatDataArrays().size()) throw nb::index_error();
            return self.getFloatDataArrays()[i];
        }, nb::rv_policy::reference_internal, "i"_a,
            "Returns a live view of the float data array at index i. Chain .data_view() on it for zero-copy numpy access into this object's storage. The view aliases this object's storage: edits through it are visible immediately, and it stays valid only until the data array list is resized or reordered. The parent object is kept alive automatically. For an owned copy use getFloatDataArrays()[i].")
        .def("getFloatDataArrays", [](OpenMS::MSChromatogram& self) -> std::vector<OpenMS::DataArrays::FloatDataArray> {
            return self.getFloatDataArrays();
        }, "Returns the float data arrays")

        .def("setFloatDataArrays", [](OpenMS::MSChromatogram& self, const std::vector<OpenMS::DataArrays::FloatDataArray>& arrays) {
            self.setFloatDataArrays(arrays);
        }, "arrays"_a, "Set the float data arrays")

        .def("_integer_data_array_count", [](const OpenMS::MSChromatogram& self) { return self.getIntegerDataArrays().size(); })
        .def("integer_data_array_view", [](OpenMS::MSChromatogram& self, size_t i) -> OpenMS::DataArrays::IntegerDataArray& {
            if (i >= self.getIntegerDataArrays().size()) throw nb::index_error();
            return self.getIntegerDataArrays()[i];
        }, nb::rv_policy::reference_internal, "i"_a,
            "Returns a live view of the integer data array at index i. Chain .data_view() on it for zero-copy numpy access into this object's storage. The view aliases this object's storage: edits through it are visible immediately, and it stays valid only until the data array list is resized or reordered. The parent object is kept alive automatically. For an owned copy use getIntegerDataArrays()[i].")
        .def("getIntegerDataArrays", [](OpenMS::MSChromatogram& self) -> std::vector<OpenMS::DataArrays::IntegerDataArray> {
            return self.getIntegerDataArrays();
        }, "Returns the integer data arrays")

        .def("setIntegerDataArrays", [](OpenMS::MSChromatogram& self, const std::vector<OpenMS::DataArrays::IntegerDataArray>& arrays) {
            self.setIntegerDataArrays(arrays);
        }, "arrays"_a, "Set the integer data arrays")

        .def("_string_data_array_count", [](const OpenMS::MSChromatogram& self) { return self.getStringDataArrays().size(); })
        .def("string_data_array_view", [](OpenMS::MSChromatogram& self, size_t i) -> OpenMS::DataArrays::StringDataArray& {
            if (i >= self.getStringDataArrays().size()) throw nb::index_error();
            return self.getStringDataArrays()[i];
        }, nb::rv_policy::reference_internal, "i"_a,
            "Returns a live view of the string data array at index i. Chain .data_view() on it for zero-copy numpy access into this object's storage. The view aliases this object's storage: edits through it are visible immediately, and it stays valid only until the data array list is resized or reordered. The parent object is kept alive automatically. For an owned copy use getStringDataArrays()[i].")
        .def("getStringDataArrays", [](OpenMS::MSChromatogram& self) -> std::vector<OpenMS::DataArrays::StringDataArray> {
            return self.getStringDataArrays();
        }, "Returns the string data arrays")

        .def("setStringDataArrays", [](OpenMS::MSChromatogram& self, const std::vector<OpenMS::DataArrays::StringDataArray>& arrays) {
            self.setStringDataArrays(arrays);
        }, "arrays"_a, "Set the string data arrays")
        .def("__repr__", [](const OpenMS::MSChromatogram& self) {
            std::ostringstream oss;
            oss << "MSChromatogram(native_id='" << std::string(self.getNativeID())
                << "', num_peaks=" << self.size() << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::MSChromatogram& self) { return nb::cast(self).attr("__repr__")(); })
        ;
    def_MetaInfoInterface<OpenMS::MSChromatogram>(mschromatogram_class);

} // NB_MODULE
