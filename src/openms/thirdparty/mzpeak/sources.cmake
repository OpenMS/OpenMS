# src/openms/thirdparty/mzpeak/sources.cmake
# The .cpp files of the vendored mzPeak library (the lib_sources list of its meson.build) plus the Boost.JSON
# header-only instantiation. Compiled into the OpenMS_mzpeak static library (CMakeLists.txt here).

set(OpenMS_mzpeak_sources
  src/data/array_index.cpp
  src/data/signals.cpp
  src/index.cpp
  src/io/directory.cpp
  src/io/file.cpp
  src/io/zip.cpp
  src/metadata/scans.cpp
  src/metadata/spectrum.cpp
  src/metadata/table.cpp
  src/open.cpp
  src/run_metadata.cpp
  src/schema/buffer_format.cpp
  src/schema/cv.cpp
  src/schema/data_kind.cpp
  src/schema/entity_type.cpp
  src/schema/file.cpp
  src/schema/group.cpp
  src/schema/psi/array_type.cpp
  src/schema/psi/chunk_encoding.cpp
  src/schema/psi/data_type.cpp
  src/chromatogram.cpp
  src/chromatograms.cpp
  src/schema/psi/transform.cpp
  src/spectra.cpp
  src/spectrum.cpp
  src/wavelength_spectrum.cpp
  src/wavelength_spectra.cpp
  src/util/arrow.cpp
  src/util/decode_pool.cpp
  src/util/decoders.cpp
  src/util/executor.cpp
  src/util/manager.cpp
  src/util/parquet.cpp
  src/util/sha512.cpp
  src/util/json_writer.cpp
  src/util/key_runs.cpp
  src/util/metadata_model.cpp
  src/util/null_fill.cpp
  src/util/numpress.cpp
  src/util/parquet_writer.cpp
  src/util/planner.cpp
  src/util/row_group_cache.cpp
  src/util/projection.cpp
  src/util/query.cpp
  src/util/slice.cpp
  src/util/types.cpp
  src/writer.cpp
  src/boost_json_src.cpp
)

list(TRANSFORM OpenMS_mzpeak_sources PREPEND "${CMAKE_CURRENT_LIST_DIR}/")
