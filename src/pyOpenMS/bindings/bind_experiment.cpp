// pyOpenMS nanobind bindings
// Domain: experiment

#include "all_casters.h"
#include <type_traits>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/SpectrumRangeManager.h>
#include <OpenMS/KERNEL/ChromatogramRangeManager.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/HPLC.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <nanobind/make_iterator.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/operators.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/vector.h>
#include "binding_utils.h"
#include <sstream>
#include <unordered_set>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(_pyopenms_experiment, m) {
    m.doc() = "pyOpenMS experiment bindings";

    // -----------------------------------------------------------------------
    // MSExperiment
    // -----------------------------------------------------------------------
    auto msexperiment_class = nb::class_<OpenMS::MSExperiment, OpenMS::ExperimentalSettings>(m, "MSExperiment",
        R"doc(
ExperimentalSettings

In-Memory representation of a mass spectrometry experiment.
Contains the data and metadata of an experiment performed with an MS (or
HPLC and MS). This representation of an MS experiment is organized as list
of spectra and chromatograms and provides an in-memory representation of
popular mass-spectrometric file formats such as mzXML or mzML. The
meta-data associated with an experiment is contained in
ExperimentalSettings (by inheritance) while the raw data (as well as
spectra and chromatogram level meta data) is stored in objects of type
MSSpectrum and MSChromatogram, which are accessible through the getSpectrum
and getChromatogram functions.
Spectra can be accessed by direct iteration or by getSpectrum(),
while chromatograms are accessed through getChromatogram().
See help(ExperimentalSettings) for information about meta-data.
Usage:
.. code-block:: python
exp = MSExperiment()
MzMLFile().load(path_to_file, exp)
for spectrum in exp:
print(spectrum.size()) # prints number of peaks
mz, intensities = spectrum.get_peaks()
)doc")
        .def(nb::init<>())
        .def(nb::init<const OpenMS::MSExperiment &>())
        .def("__copy__", [](const OpenMS::MSExperiment& self) { return OpenMS::MSExperiment(self); })
        .def("__deepcopy__", [](const OpenMS::MSExperiment& self, nb::dict) { return OpenMS::MSExperiment(self); }, "memo"_a)
        .def(nb::self == nb::self)
        .def(nb::self != nb::self)
        .def("size", [](const OpenMS::MSExperiment& self) { return self.size(); }, "Returns the number of spectra")
        .def("resize", [](OpenMS::MSExperiment& self, size_t n) { return self.resize(n); }, "n"_a, "Resizes the spectrum container to the specified number of spectra")
        .def("empty", [](const OpenMS::MSExperiment& self) { return self.empty(); }, "Returns True if the experiment contains no spectra")
        .def("reserve", [](OpenMS::MSExperiment& self, size_t n) { return self.reserve(n); }, "n"_a, "Reserves space for the specified number of spectra")
        .def("begin", [](const OpenMS::MSExperiment& self) { return self.begin(); })
        .def("end", [](const OpenMS::MSExperiment& self) { return self.end(); })
        .def("reserveSpaceSpectra", [](OpenMS::MSExperiment& self, size_t s) { return self.reserveSpaceSpectra(s); }, "s"_a, "Reserves space for the specified number of spectra")
        .def("reserveSpaceChromatograms", [](OpenMS::MSExperiment& self, size_t s) { return self.reserveSpaceChromatograms(s); }, "s"_a, "Reserves space for the specified number of chromatograms")
        .def("aggregateFromMatrix", [](const OpenMS::MSExperiment& self, const OpenMS::Matrix<double>& ranges, unsigned int ms_level, const OpenMS::String& mz_agg) { return self.aggregateFromMatrix(ranges, ms_level, mz_agg); }, "ranges"_a, "ms_level"_a, "mz_agg"_a, "Aggregates intensity values for multiple m/z and RT ranges specified in a matrix")
        .def("extractXICsFromMatrix", [](const OpenMS::MSExperiment& self, const OpenMS::Matrix<double>& ranges, unsigned int ms_level, const OpenMS::String& mz_agg) { return self.extractXICsFromMatrix(ranges, ms_level, mz_agg); }, "ranges"_a, "ms_level"_a, "mz_agg"_a, "Extracts XIC chromatograms for multiple m/z and RT ranges specified in a matrix")
        .def("clearRanges", [](OpenMS::MSExperiment& self) { return self.clearRanges(); }, "Clear all ranges in all range managers")
        .def("getMinRT", [](const OpenMS::MSExperiment& self) { return self.getMinRT(); }, "Get the minimum RT value from the combined ranges")
        .def("getMaxRT", [](const OpenMS::MSExperiment& self) { return self.getMaxRT(); }, "Get the maximum RT value from the combined ranges")
        .def("getMinMZ", [](const OpenMS::MSExperiment& self) { return self.getMinMZ(); }, "Get the minimum m/z value from the combined ranges")
        .def("getMaxMZ", [](const OpenMS::MSExperiment& self) { return self.getMaxMZ(); }, "Get the maximum m/z value from the combined ranges")
        .def("getMinIntensity", [](const OpenMS::MSExperiment& self) { return self.getMinIntensity(); }, "Get the minimum intensity value from the combined ranges")
        .def("getMaxIntensity", [](const OpenMS::MSExperiment& self) { return self.getMaxIntensity(); }, "Get the maximum intensity value from the combined ranges")
        .def("getMinMobility", [](const OpenMS::MSExperiment& self) { return self.getMinMobility(); }, "Get the minimum mobility value from the combined ranges")
        .def("getMaxMobility", [](const OpenMS::MSExperiment& self) { return self.getMaxMobility(); }, "Get the maximum mobility value from the combined ranges")
        .def("updateRanges", [](OpenMS::MSExperiment& self) { return self.updateRanges(); }, "Recalculate global ranges for both spectra and chromatrograms after changes to the data has been made.")
        .def("getSize", [](const OpenMS::MSExperiment& self) { return self.getSize(); }, "Returns the total number of peaks")
        .def("sortSpectra", [](OpenMS::MSExperiment& self, bool sort_mz) { nb::gil_scoped_release release; return self.sortSpectra(sort_mz); }, "sort_mz"_a = true, "Sorts spectra by RT. If sort_mz=True also sort each peak in a spectrum by m/z")
        .def("sortChromatograms", [](OpenMS::MSExperiment& self, bool sort_rt) { nb::gil_scoped_release release; return self.sortChromatograms(sort_rt); }, "sort_rt"_a = true, "Sorts chromatograms by m/z. If sort_rt=True also sort each chromatogram RT")
        .def("isSorted", [](const OpenMS::MSExperiment& self, bool check_mz) { return self.isSorted(check_mz); }, "check_mz"_a = true, "Checks if all spectra are sorted with respect to ascending RT")
        .def("reset", [](OpenMS::MSExperiment& self) { return self.reset(); }, "Clears all data and meta data")
        .def("clearMetaDataArrays", [](OpenMS::MSExperiment& self) { return self.clearMetaDataArrays(); }, "Clears the meta data arrays of all contained spectra")
        .def("getExperimentalSettings", [](OpenMS::MSExperiment& self) -> OpenMS::ExperimentalSettings & { return self.getExperimentalSettings(); }, nb::rv_policy::reference_internal, "Returns the meta information of this experiment")
        .def("getPrimaryMSRunPath", [](const OpenMS::MSExperiment& self) { std::vector<OpenMS::String> toFill; self.getPrimaryMSRunPath(toFill); return toFill; }, "References to the first MS file(s) after conversions. Used to trace results back to original data.")
        .def("getPrecursorSpectrum", [](const OpenMS::MSExperiment& self, int zero_based_index) { return self.getPrecursorSpectrum(zero_based_index); }, "zero_based_index"_a, "Returns the index of the precursor spectrum for spectrum at index @p zero_based_index")
        .def("swap", [](OpenMS::MSExperiment& self, OpenMS::MSExperiment& from) { return self.swap(from); }, "from"_a, "Swaps the content of this experiment with another")
        .def("setSpectra", [](OpenMS::MSExperiment& self, const std::vector<OpenMS::MSSpectrum>& spectra) { return self.setSpectra(spectra); }, "spectra"_a, "Sets the spectrum list")
        .def("setSpectra", [](OpenMS::MSExperiment& self, std::vector<OpenMS::MSSpectrum>& spectra) { return self.setSpectra(spectra); }, "spectra"_a, "Sets the spectrum list")
        .def("addSpectrum", [](OpenMS::MSExperiment& self, const OpenMS::MSSpectrum& spectrum) { return self.addSpectrum(spectrum); }, "spectrum"_a, "Adds a spectrum to the experiment")
        .def("addSpectrum", [](OpenMS::MSExperiment& self, OpenMS::MSSpectrum& spectrum) { return self.addSpectrum(spectrum); }, "spectrum"_a, "Adds a spectrum to the experiment")
        .def("getSpectra", [](OpenMS::MSExperiment& self) -> std::vector<OpenMS::MSSpectrum> & { return self.getSpectra(); }, nb::rv_policy::reference_internal, "Returns the list of spectra")
        .def("setChromatograms", [](OpenMS::MSExperiment& self, const std::vector<OpenMS::MSChromatogram>& chromatograms) { return self.setChromatograms(chromatograms); }, "chromatograms"_a, "Sets the chromatogram list")
        .def("setChromatograms", [](OpenMS::MSExperiment& self, std::vector<OpenMS::MSChromatogram>& chromatograms) { return self.setChromatograms(chromatograms); }, "chromatograms"_a, "Sets the chromatogram list")
        .def("addChromatogram", [](OpenMS::MSExperiment& self, const OpenMS::MSChromatogram& chromatogram) { return self.addChromatogram(chromatogram); }, "chromatogram"_a, "Adds a chromatogram to the experiment")
        .def("addChromatogram", [](OpenMS::MSExperiment& self, OpenMS::MSChromatogram& chrom) { return self.addChromatogram(chrom); }, "chrom"_a, "Adds a chromatogram to the experiment")
        .def("getChromatograms", [](OpenMS::MSExperiment& self) -> std::vector<OpenMS::MSChromatogram> & { return self.getChromatograms(); }, nb::rv_policy::reference_internal, "Returns the list of chromatograms")
        .def("getChromatogram", [](OpenMS::MSExperiment& self, size_t id) -> OpenMS::MSChromatogram& { return self.getChromatogram(id); }, nb::rv_policy::reference_internal, "id"_a, "Returns a single chromatogram by index")
        .def("getNrSpectra", [](const OpenMS::MSExperiment& self) { return self.getNrSpectra(); }, "Returns the number of MS spectra")
        .def("getNrChromatograms", [](const OpenMS::MSExperiment& self) { return self.getNrChromatograms(); }, "Returns the number of chromatograms")
        .def("getMSLevels", [](const OpenMS::MSExperiment& self) { return self.getMSLevels(); }, "Returns a sorted list of unique MS levels in the experiment")
        .def("calculateTIC", [](const OpenMS::MSExperiment& self, float rt_bin_size, unsigned int ms_level) { return self.calculateTIC(rt_bin_size, ms_level); }, "rt_bin_size"_a = 0, "ms_level"_a = 1, "Returns the total ion chromatogram")
        .def("clear", [](OpenMS::MSExperiment& self, bool clear_meta_data) { return self.clear(clear_meta_data); }, "clear_meta_data"_a, "Clear all spectra data and meta data (if called with True)")
        .def("spectrumRanges", [](const OpenMS::MSExperiment& self) -> const OpenMS::SpectrumRangeManager & { return self.spectrumRanges(); }, nb::rv_policy::reference_internal, "Returns a reference to the spectrum range manager")
        .def("chromatogramRanges", [](const OpenMS::MSExperiment& self) -> const OpenMS::ChromatogramRangeManager & { return self.chromatogramRanges(); }, nb::rv_policy::reference_internal, "Returns a reference to the chromatogram range manager")
        .def("getSample", [](const OpenMS::MSExperiment& self) -> const OpenMS::Sample & { return self.getSample(); }, nb::rv_policy::reference_internal, "Returns a reference to the sample description")
        .def("setSample", [](OpenMS::MSExperiment& self, const OpenMS::Sample& sample) { return self.setSample(sample); }, "sample"_a, "Sets the sample description")
        .def("getSourceFiles", [](const OpenMS::MSExperiment& self) -> const std::vector<OpenMS::SourceFile> & { return self.getSourceFiles(); }, nb::rv_policy::reference_internal, "Returns a reference to the source data file")
        .def("setSourceFiles", [](OpenMS::MSExperiment& self, const std::vector<OpenMS::SourceFile>& source_files) { return self.setSourceFiles(source_files); }, "source_files"_a, "Sets the source data file")
        .def("getContacts", [](const OpenMS::MSExperiment& self) -> const std::vector<OpenMS::ContactPerson> & { return self.getContacts(); }, nb::rv_policy::reference_internal, "Returns a reference to the list of contact persons")
        .def("setContacts", [](OpenMS::MSExperiment& self, const std::vector<OpenMS::ContactPerson>& contacts) { return self.setContacts(contacts); }, "contacts"_a, "Sets the list of contact persons")
        .def("getInstrument", [](const OpenMS::MSExperiment& self) -> const OpenMS::Instrument & { return self.getInstrument(); }, nb::rv_policy::reference_internal, "Returns a reference to the MS instrument description")
        .def("setInstrument", [](OpenMS::MSExperiment& self, const OpenMS::Instrument& instrument) { return self.setInstrument(instrument); }, "instrument"_a, "Sets the MS instrument description")
        .def("getHPLC", [](const OpenMS::MSExperiment& self) -> const OpenMS::HPLC & { return self.getHPLC(); }, nb::rv_policy::reference_internal, "Returns a reference to the description of the HPLC run")
        .def("setHPLC", [](OpenMS::MSExperiment& self, const OpenMS::HPLC& hplc) { return self.setHPLC(hplc); }, "hplc"_a, "Sets the description of the HPLC run")
        .def("getDateTime", [](const OpenMS::MSExperiment& self) -> const OpenMS::DateTime & { return self.getDateTime(); }, nb::rv_policy::reference_internal, "Returns the date the experiment was performed")
        .def("setDateTime", [](OpenMS::MSExperiment& self, const OpenMS::DateTime& date) { return self.setDateTime(date); }, "date"_a, "Sets the date the experiment was performed")
        .def("getComment", [](const OpenMS::MSExperiment& self) { return self.getComment(); }, "Returns the free-text comment")
        .def("setComment", [](OpenMS::MSExperiment& self, const OpenMS::String& comment) { return self.setComment(comment); }, "comment"_a, "Sets the free-text comment")
        .def("getFractionIdentifier", [](const OpenMS::MSExperiment& self) { return self.getFractionIdentifier(); }, "Returns fraction identifier")
        .def("setFractionIdentifier", [](OpenMS::MSExperiment& self, const OpenMS::String& fraction_identifier) { return self.setFractionIdentifier(fraction_identifier); }, "fraction_identifier"_a, "Sets the fraction identifier")

        .def("__iter__", [](OpenMS::MSExperiment& self) { return nb::make_iterator<nb::rv_policy::reference_internal>(nb::type<OpenMS::MSExperiment>(), "MSExperiment_iter", self.begin(), self.end()); })
        .def("__len__", [](OpenMS::MSExperiment& self) { return self.size(); })

        .def("__getitem__", [](OpenMS::MSExperiment& self, size_t i) -> OpenMS::MSSpectrum& {
            if (i >= self.size()) throw nb::index_error();
            return self[i];
        }, "i"_a, nb::rv_policy::reference_internal)
        .def("__setitem__", [](OpenMS::MSExperiment& self, size_t i, const OpenMS::MSSpectrum& val) {
            if (i >= self.size()) throw nb::index_error();
            self[i] = val;
        }, "i"_a, "val"_a, "Sets spectrum at index i")
        .def("getSpectrum", [](OpenMS::MSExperiment& self, size_t id) -> OpenMS::MSSpectrum& { return self.getSpectrum(id); }, nb::rv_policy::reference_internal, "id"_a, "Returns a single spectrum by index")

        .def("get2DPeakData", [](const OpenMS::MSExperiment& self,
                                double min_rt, double max_rt,
                                double min_mz, double max_mz,
                                size_t ms_level) {
            std::vector<float> rt, mz, intensity;
            self.get2DPeakData(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, intensity);
            const size_t n = rt.size();

            std::unique_ptr<float[]> rt_uptr(new float[n]);
            std::unique_ptr<float[]> mz_uptr(new float[n]);
            std::unique_ptr<float[]> int_uptr(new float[n]);
            std::copy(rt.begin(), rt.end(), rt_uptr.get());
            std::copy(mz.begin(), mz.end(), mz_uptr.get());
            std::copy(intensity.begin(), intensity.end(), int_uptr.get());

            float* rt_data = rt_uptr.release();
            float* mz_data = mz_uptr.release();
            float* int_data = int_uptr.release();

            nb::capsule rt_owner(rt_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            nb::capsule mz_owner(mz_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            nb::capsule int_owner(int_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });

            return nb::make_tuple(
                nb::ndarray<nb::numpy, float, nb::ndim<1>>(rt_data, {n}, rt_owner),
                nb::ndarray<nb::numpy, float, nb::ndim<1>>(mz_data, {n}, mz_owner),
                nb::ndarray<nb::numpy, float, nb::ndim<1>>(int_data, {n}, int_owner)
            );
        }, "min_rt"_a, "max_rt"_a, "min_mz"_a, "max_mz"_a, "ms_level"_a,
           "Retrieves peak data in the given mz-rt range as (rt, mz, intensity) numpy arrays")

        .def("get2DPeakDataIM", [](const OpenMS::MSExperiment& self,
                                   double min_rt, double max_rt,
                                   double min_mz, double max_mz,
                                   size_t ms_level) {
            std::vector<float> rt, mz, intensity, ion_mobility;
            self.get2DPeakDataIM(min_rt, max_rt, min_mz, max_mz, ms_level, rt, mz, intensity, ion_mobility);
            const size_t n = rt.size();

            std::unique_ptr<float[]> rt_uptr(new float[n]);
            std::unique_ptr<float[]> mz_uptr(new float[n]);
            std::unique_ptr<float[]> int_uptr(new float[n]);
            std::unique_ptr<float[]> im_uptr(new float[n]);
            std::copy(rt.begin(), rt.end(), rt_uptr.get());
            std::copy(mz.begin(), mz.end(), mz_uptr.get());
            std::copy(intensity.begin(), intensity.end(), int_uptr.get());
            std::copy(ion_mobility.begin(), ion_mobility.end(), im_uptr.get());

            float* rt_data = rt_uptr.release();
            float* mz_data = mz_uptr.release();
            float* int_data = int_uptr.release();
            float* im_data = im_uptr.release();

            nb::capsule rt_owner(rt_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            nb::capsule mz_owner(mz_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            nb::capsule int_owner(int_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            nb::capsule im_owner(im_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });

            return nb::make_tuple(
                nb::ndarray<nb::numpy, float, nb::ndim<1>>(rt_data, {n}, rt_owner),
                nb::ndarray<nb::numpy, float, nb::ndim<1>>(mz_data, {n}, mz_owner),
                nb::ndarray<nb::numpy, float, nb::ndim<1>>(int_data, {n}, int_owner),
                nb::ndarray<nb::numpy, float, nb::ndim<1>>(im_data, {n}, im_owner)
            );
        }, "min_rt"_a, "max_rt"_a, "min_mz"_a, "max_mz"_a, "ms_level"_a,
           "Retrieves peak data with ion mobility in the given mz-rt range as (rt, mz, intensity, ion_mobility) numpy arrays")

        .def("get2DPeakDataLong", [](const OpenMS::MSExperiment& self,
                                     double min_rt, double max_rt,
                                     double min_mz, double max_mz,
                                     size_t ms_level) {
            // Returns (rt, mz, intensity) as flat numpy arrays with full precision
            // (float64 for rt/mz, float32 for intensity). Replaces the Python addon.

            // Pass 1: count total matching peaks
            size_t total_peaks = 0;
            for (size_t i = 0; i < self.size(); ++i) {
                const auto& spec = self[i];
                if (spec.getMSLevel() != ms_level) continue;
                double rt = spec.getRT();
                if (rt < min_rt || rt > max_rt) continue;
                for (size_t j = 0; j < spec.size(); ++j) {
                    double mz = spec[j].getMZ();
                    if (mz >= min_mz && mz <= max_mz) ++total_peaks;
                }
            }

            // Allocate output arrays
            auto rt_uptr = std::make_unique<double[]>(total_peaks);
            auto mz_uptr = std::make_unique<double[]>(total_peaks);
            auto int_uptr = std::make_unique<float[]>(total_peaks);

            // Pass 2: fill arrays
            size_t pos = 0;
            for (size_t i = 0; i < self.size(); ++i) {
                const auto& spec = self[i];
                if (spec.getMSLevel() != ms_level) continue;
                double rt = spec.getRT();
                if (rt < min_rt || rt > max_rt) continue;
                for (size_t j = 0; j < spec.size(); ++j) {
                    double mz = spec[j].getMZ();
                    if (mz >= min_mz && mz <= max_mz) {
                        rt_uptr[pos] = rt;
                        mz_uptr[pos] = mz;
                        int_uptr[pos] = spec[j].getIntensity();
                        ++pos;
                    }
                }
            }

            double* rt_data = rt_uptr.release();
            double* mz_data = mz_uptr.release();
            float* int_data = int_uptr.release();

            nb::capsule rt_owner(rt_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            nb::capsule mz_owner(mz_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            nb::capsule int_owner(int_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });

            return nb::make_tuple(
                nb::ndarray<nb::numpy, double, nb::ndim<1>>(rt_data, {total_peaks}, rt_owner),
                nb::ndarray<nb::numpy, double, nb::ndim<1>>(mz_data, {total_peaks}, mz_owner),
                nb::ndarray<nb::numpy, float, nb::ndim<1>>(int_data, {total_peaks}, int_owner)
            );
        }, "min_rt"_a, "max_rt"_a, "min_mz"_a, "max_mz"_a, "ms_level"_a,
           "Retrieves peak data in the given mz-rt range as (rt, mz, intensity) numpy arrays with full precision")

        .def("get2DPeakDataSemiWide", [](const OpenMS::MSExperiment& self,
                                         std::vector<int> ms_levels) {
            // Returns (rt, ms_level, mz_flat, intensity_flat, offsets) where
            // mz_flat/intensity_flat are contiguous buffers and offsets[i]:offsets[i+1]
            // gives the slice for spectrum i. Use np.split(mz_flat, offsets[1:-1])
            // on the Python side to get zero-copy views per spectrum.
            std::unordered_set<int> level_set(ms_levels.begin(), ms_levels.end());

            // Pass 1: count matching spectra and total peaks
            size_t n_spectra = 0;
            size_t total_peaks = 0;
            for (size_t i = 0; i < self.size(); ++i) {
                const auto& spec = self[i];
                if (level_set.count(spec.getMSLevel())) {
                    ++n_spectra;
                    total_peaks += spec.size();
                }
            }

            // Allocate output arrays
            auto rt_uptr = std::make_unique<float[]>(n_spectra);
            auto ms_uptr = std::make_unique<uint8_t[]>(n_spectra);
            auto off_uptr = std::make_unique<int64_t[]>(n_spectra + 1);
            auto mz_uptr = std::make_unique<double[]>(total_peaks);
            auto int_uptr = std::make_unique<float[]>(total_peaks);

            // Pass 2: fill arrays
            size_t spec_idx = 0;
            size_t peak_pos = 0;
            for (size_t i = 0; i < self.size(); ++i) {
                const auto& spec = self[i];
                if (!level_set.count(spec.getMSLevel())) continue;

                rt_uptr[spec_idx] = static_cast<float>(spec.getRT());
                ms_uptr[spec_idx] = static_cast<uint8_t>(spec.getMSLevel());
                off_uptr[spec_idx] = static_cast<int64_t>(peak_pos);

                for (size_t j = 0; j < spec.size(); ++j) {
                    mz_uptr[peak_pos] = spec[j].getMZ();
                    int_uptr[peak_pos] = spec[j].getIntensity();
                    ++peak_pos;
                }
                ++spec_idx;
            }
            off_uptr[n_spectra] = static_cast<int64_t>(peak_pos);

            // Transfer ownership to capsules
            float* rt_data = rt_uptr.release();
            uint8_t* ms_data = ms_uptr.release();
            int64_t* off_data = off_uptr.release();
            double* mz_data = mz_uptr.release();
            float* int_data = int_uptr.release();

            nb::capsule rt_owner(rt_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });
            nb::capsule ms_owner(ms_data, [](void* p) noexcept { delete[] static_cast<uint8_t*>(p); });
            nb::capsule off_owner(off_data, [](void* p) noexcept { delete[] static_cast<int64_t*>(p); });
            nb::capsule mz_owner(mz_data, [](void* p) noexcept { delete[] static_cast<double*>(p); });
            nb::capsule int_owner(int_data, [](void* p) noexcept { delete[] static_cast<float*>(p); });

            return nb::make_tuple(
                nb::ndarray<nb::numpy, float, nb::ndim<1>>(rt_data, {n_spectra}, rt_owner),
                nb::ndarray<nb::numpy, uint8_t, nb::ndim<1>>(ms_data, {n_spectra}, ms_owner),
                nb::ndarray<nb::numpy, double, nb::ndim<1>>(mz_data, {total_peaks}, mz_owner),
                nb::ndarray<nb::numpy, float, nb::ndim<1>>(int_data, {total_peaks}, int_owner),
                nb::ndarray<nb::numpy, int64_t, nb::ndim<1>>(off_data, {n_spectra + 1}, off_owner)
            );
        }, "ms_levels"_a,
           "Retrieves peak data as flat arrays + offsets for efficient semi-wide DataFrame construction. "
           "Returns (rt, ms_level, mz_flat, intensity_flat, offsets) numpy arrays.")

        .def("rasterizeRTMZ", [](OpenMS::MSExperiment& self,
                                nb::ndarray<float, nb::ndim<2>, nb::device::cpu> output,
                                double min_rt, double max_rt,
                                double min_mz, double max_mz,
                                unsigned int ms_level,
                                const std::string& aggregation) {
            // Check for C-contiguous array (legacy pyOpenMS compatibility)
            // DLTensor strides are in elements, not bytes
            // C-contiguous: stride(1)=1, stride(0)=shape(1)
            if (output.stride(1) != 1 ||
                output.stride(0) != static_cast<int64_t>(output.shape(1))) {
                throw std::invalid_argument("Output array must be C-contiguous. "
                    "Use numpy.ascontiguousarray() to convert.");
            }

            // Output array has shape [mz_bins, rt_bins]
            size_t mz_bins = output.shape(0);
            size_t rt_bins = output.shape(1);
            float* output_ptr = output.data();

            OpenMS::MSExperiment::RasterAggregation agg_mode;
            if (aggregation == "sum" || aggregation == "SUM") {
                agg_mode = OpenMS::MSExperiment::RasterAggregation::SUM;
            } else if (aggregation == "max" || aggregation == "MAX") {
                agg_mode = OpenMS::MSExperiment::RasterAggregation::MAX;
            } else {
                throw std::invalid_argument("Invalid aggregation mode '" + aggregation + "'. Must be 'sum' or 'max'.");
            }

            self.rasterizeRTMZ(output_ptr, rt_bins, mz_bins, min_rt, max_rt, min_mz, max_mz, ms_level, agg_mode);
        }, "output"_a, "min_rt"_a, "max_rt"_a, "min_mz"_a, "max_mz"_a, "ms_level"_a, "aggregation"_a = "sum",
           "Rasterize peak data into a 2D intensity matrix. Output shape is [mz_bins, rt_bins].")
        .def("__repr__", [](const OpenMS::MSExperiment& self) {
            std::ostringstream oss;
            oss << "MSExperiment(num_spectra=" << self.getNrSpectra()
                << ", num_chromatograms=" << self.getNrChromatograms() << ")";
            return oss.str();
        })
        .def("__str__", [](const OpenMS::MSExperiment& self) { return nb::cast(self).attr("__repr__")(); })
        ;
    // MSExperimentRasterAggregation enum nested under MSExperiment
    nb::enum_<OpenMS::MSExperiment::RasterAggregation>(msexperiment_class, "MSExperimentRasterAggregation")
        .value("SUM", OpenMS::MSExperiment::RasterAggregation::SUM)
        .value("MAX", OpenMS::MSExperiment::RasterAggregation::MAX)
        ;

} // NB_MODULE
