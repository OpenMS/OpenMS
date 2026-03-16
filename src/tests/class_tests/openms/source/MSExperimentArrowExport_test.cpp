// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>
///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/type.h>
#include <arrow/c/bridge.h>

using namespace OpenMS;
using namespace std;

/////////////////////////////////////////////////////////////
// Helper functions to convert C Data Interface to Arrow Table
/////////////////////////////////////////////////////////////

std::shared_ptr<arrow::Table> exportSpectraViaInterface(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config)
{
  ArrowSchema schema;
  ArrowArray array;

  if (!MSExperimentArrowExport::exportSpectraToArrowCDataInterface(exp, config, &schema, &array))
  {
    return nullptr;
  }

  // Import as RecordBatch
  auto batch_result = arrow::ImportRecordBatch(&array, &schema);
  if (!batch_result.ok())
  {
    return nullptr;
  }

  // Convert to Table
  return arrow::Table::FromRecordBatches({batch_result.ValueOrDie()}).ValueOrDie();
}

std::shared_ptr<arrow::Table> exportChromatogramsViaInterface(
  const MSExperiment& exp,
  const ArrowChromatogramExportConfig& config)
{
  ArrowSchema schema;
  ArrowArray array;

  if (!MSExperimentArrowExport::exportChromatogramsToArrowCDataInterface(exp, config, &schema, &array))
  {
    return nullptr;
  }

  // Import as RecordBatch
  auto batch_result = arrow::ImportRecordBatch(&array, &schema);
  if (!batch_result.ok())
  {
    return nullptr;
  }

  // Convert to Table
  return arrow::Table::FromRecordBatches({batch_result.ValueOrDie()}).ValueOrDie();
}

/////////////////////////////////////////////////////////////
// Helper function to create a test MSExperiment
/////////////////////////////////////////////////////////////

MSExperiment createTestExperiment()
{
  MSExperiment exp;

  // Create MS1 spectrum with 3 peaks
  MSSpectrum ms1;
  ms1.setRT(100.0);
  ms1.setMSLevel(1);
  ms1.setNativeID("spectrum=0");
  ms1.push_back(Peak1D(100.0, 1000.0));
  ms1.push_back(Peak1D(200.0, 2000.0));
  ms1.push_back(Peak1D(300.0, 3000.0));
  exp.addSpectrum(ms1);

  // Create MS2 spectrum with precursor info
  MSSpectrum ms2;
  ms2.setRT(110.0);
  ms2.setMSLevel(2);
  ms2.setNativeID("spectrum=1");

  Precursor prec;
  prec.setMZ(500.0);
  prec.setCharge(2);
  prec.setIntensity(50000.0);
  prec.setIsolationWindowLowerOffset(1.5);
  prec.setIsolationWindowUpperOffset(1.5);
  ms2.setPrecursors({prec});

  ms2.push_back(Peak1D(150.0, 500.0));
  ms2.push_back(Peak1D(250.0, 750.0));
  exp.addSpectrum(ms2);

  // Create another MS1 spectrum
  MSSpectrum ms1_2;
  ms1_2.setRT(200.0);
  ms1_2.setMSLevel(1);
  ms1_2.setNativeID("spectrum=2");
  ms1_2.push_back(Peak1D(150.0, 1500.0));
  ms1_2.push_back(Peak1D(250.0, 2500.0));
  exp.addSpectrum(ms1_2);

  return exp;
}

MSExperiment createEmptyExperiment()
{
  return MSExperiment();
}

/////////////////////////////////////////////////////////////

START_TEST(MSExperimentArrowExport, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(exportSpectraToArrowCDataInterface - empty experiment)
{
  MSExperiment exp = createEmptyExperiment();
  ArrowSpectraExportConfig config;

  auto table = exportSpectraViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
  // Schema should still be present
  TEST_EQUAL(table->num_columns() > 0, true)
}
END_SECTION

START_SECTION(exportSpectraToArrowCDataInterface - long format basic)
{
  MSExperiment exp = createTestExperiment();
  ArrowSpectraExportConfig config;
  config.format = ArrowExportFormat::Long;

  auto table = exportSpectraViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  // Total peaks: 3 (MS1) + 2 (MS2) + 2 (MS1) = 7
  TEST_EQUAL(table->num_rows(), 7)

  // Check schema
  auto schema = table->schema();
  TEST_EQUAL(schema->GetFieldIndex("mz") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("intensity") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("rt") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("spectrum_index") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("ms_level") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("native_id") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("precursor_mz") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("precursor_charge") >= 0, true)

  // Check mz column type
  auto mz_field = schema->GetFieldByName("mz");
  TEST_EQUAL(mz_field->type()->id(), arrow::Type::DOUBLE)

  // Check intensity column type
  auto intensity_field = schema->GetFieldByName("intensity");
  TEST_EQUAL(intensity_field->type()->id(), arrow::Type::FLOAT)

  // Check ms_level column type
  auto ms_level_field = schema->GetFieldByName("ms_level");
  TEST_EQUAL(ms_level_field->type()->id(), arrow::Type::UINT8)
}
END_SECTION

START_SECTION(exportSpectraToArrowCDataInterface - long format with MS level filter)
{
  MSExperiment exp = createTestExperiment();
  ArrowSpectraExportConfig config;
  config.format = ArrowExportFormat::Long;
  config.ms_levels = {1}; // Only MS1

  auto table = exportSpectraViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  // Only MS1 peaks: 3 + 2 = 5
  TEST_EQUAL(table->num_rows(), 5)
}
END_SECTION

START_SECTION(exportSpectraToArrowCDataInterface - long format with RT filter)
{
  MSExperiment exp = createTestExperiment();
  ArrowSpectraExportConfig config;
  config.format = ArrowExportFormat::Long;
  config.min_rt = 100.0;
  config.max_rt = 110.0;

  auto table = exportSpectraViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  // Spectra in RT range: ms1 (3 peaks) + ms2 (2 peaks) = 5
  TEST_EQUAL(table->num_rows(), 5)
}
END_SECTION

START_SECTION(exportSpectraToArrowCDataInterface - long format with m/z filter)
{
  MSExperiment exp = createTestExperiment();
  ArrowSpectraExportConfig config;
  config.format = ArrowExportFormat::Long;
  config.min_mz = 150.0;
  config.max_mz = 250.0;

  auto table = exportSpectraViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  // Peaks in m/z range: 200.0 (ms1), 150.0+250.0 (ms2), 150.0+250.0 (ms1_2) = 5
  TEST_EQUAL(table->num_rows(), 5)
}
END_SECTION

START_SECTION(exportSpectraToArrowCDataInterface - long format with column selection)
{
  MSExperiment exp = createTestExperiment();
  ArrowSpectraExportConfig config;
  config.format = ArrowExportFormat::Long;
  config.columns = {"mz", "intensity", "rt"};

  auto table = exportSpectraViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_columns(), 3)
  TEST_EQUAL(table->schema()->GetFieldIndex("mz") >= 0, true)
  TEST_EQUAL(table->schema()->GetFieldIndex("intensity") >= 0, true)
  TEST_EQUAL(table->schema()->GetFieldIndex("rt") >= 0, true)
  TEST_EQUAL(table->schema()->GetFieldIndex("spectrum_index"), -1)
}
END_SECTION

START_SECTION(exportSpectraToArrowCDataInterface - long format no precursor info)
{
  MSExperiment exp = createTestExperiment();
  ArrowSpectraExportConfig config;
  config.format = ArrowExportFormat::Long;
  config.include_precursor_info = false;

  auto table = exportSpectraViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  auto schema = table->schema();
  TEST_EQUAL(schema->GetFieldIndex("precursor_mz"), -1)
  TEST_EQUAL(schema->GetFieldIndex("precursor_charge"), -1)
  TEST_EQUAL(schema->GetFieldIndex("isolation_lower"), -1)
}
END_SECTION

START_SECTION(exportSpectraToArrowCDataInterface - semi-wide format basic)
{
  MSExperiment exp = createTestExperiment();
  ArrowSpectraExportConfig config;
  config.format = ArrowExportFormat::SemiWide;

  auto table = exportSpectraViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  // One row per spectrum
  TEST_EQUAL(table->num_rows(), 3)

  // Check schema - mz and intensity should be list types
  auto schema = table->schema();
  auto mz_field = schema->GetFieldByName("mz");
  TEST_EQUAL(mz_field->type()->id(), arrow::Type::LIST)

  auto intensity_field = schema->GetFieldByName("intensity");
  TEST_EQUAL(intensity_field->type()->id(), arrow::Type::LIST)

  // Scalar columns should be regular types
  auto rt_field = schema->GetFieldByName("rt");
  TEST_EQUAL(rt_field->type()->id(), arrow::Type::FLOAT)
}
END_SECTION

START_SECTION(exportSpectraToArrowCDataInterface - semi-wide format with MS level filter)
{
  MSExperiment exp = createTestExperiment();
  ArrowSpectraExportConfig config;
  config.format = ArrowExportFormat::SemiWide;
  config.ms_levels = {2}; // Only MS2

  auto table = exportSpectraViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1)
}
END_SECTION

START_SECTION(getSpectraArrowColumnNames - long format)
{
  MSExperiment exp = createTestExperiment();
  ArrowSpectraExportConfig config;
  config.format = ArrowExportFormat::Long;

  auto columns = MSExperimentArrowExport::getSpectraArrowColumnNames(exp, config);

  // Check required columns are present
  TEST_EQUAL(std::find(columns.begin(), columns.end(), "mz") != columns.end(), true)
  TEST_EQUAL(std::find(columns.begin(), columns.end(), "intensity") != columns.end(), true)
  TEST_EQUAL(std::find(columns.begin(), columns.end(), "rt") != columns.end(), true)
  TEST_EQUAL(std::find(columns.begin(), columns.end(), "spectrum_index") != columns.end(), true)
  TEST_EQUAL(std::find(columns.begin(), columns.end(), "ms_level") != columns.end(), true)
  TEST_EQUAL(std::find(columns.begin(), columns.end(), "precursor_mz") != columns.end(), true)
}
END_SECTION

START_SECTION(getSpectraArrowColumnNames - with column filter)
{
  MSExperiment exp = createTestExperiment();
  ArrowSpectraExportConfig config;
  config.columns = {"mz", "intensity"};

  auto columns = MSExperimentArrowExport::getSpectraArrowColumnNames(exp, config);

  TEST_EQUAL(columns.size(), 2)
  TEST_EQUAL(std::find(columns.begin(), columns.end(), "mz") != columns.end(), true)
  TEST_EQUAL(std::find(columns.begin(), columns.end(), "intensity") != columns.end(), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Chromatogram export tests
/////////////////////////////////////////////////////////////

START_SECTION(exportChromatogramsToArrowCDataInterface - empty chromatograms)
{
  MSExperiment exp;
  ArrowChromatogramExportConfig config;

  auto table = exportChromatogramsViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 0)
}
END_SECTION

START_SECTION(exportChromatogramsToArrowCDataInterface - long format)
{
  MSExperiment exp;

  // Create a simple chromatogram
  MSChromatogram chrom;
  chrom.setNativeID("TIC");
  chrom.getPrecursor().setMZ(500.0);
  chrom.getProduct().setMZ(200.0);
  chrom.push_back(ChromatogramPeak(10.0, 100.0));
  chrom.push_back(ChromatogramPeak(20.0, 200.0));
  chrom.push_back(ChromatogramPeak(30.0, 150.0));
  exp.addChromatogram(chrom);

  // Add another chromatogram
  MSChromatogram chrom2;
  chrom2.setNativeID("XIC_1");
  chrom2.getPrecursor().setMZ(600.0);
  chrom2.getProduct().setMZ(300.0);
  chrom2.push_back(ChromatogramPeak(15.0, 50.0));
  chrom2.push_back(ChromatogramPeak(25.0, 75.0));
  exp.addChromatogram(chrom2);

  ArrowChromatogramExportConfig config;
  config.format = ArrowExportFormat::Long;

  auto table = exportChromatogramsViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  // Total points: 3 + 2 = 5
  TEST_EQUAL(table->num_rows(), 5)

  auto schema = table->schema();
  TEST_EQUAL(schema->GetFieldIndex("rt") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("intensity") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("chromatogram_index") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("precursor_mz") >= 0, true)
  TEST_EQUAL(schema->GetFieldIndex("product_mz") >= 0, true)
}
END_SECTION

START_SECTION(exportChromatogramsToArrowCDataInterface - semi-wide format)
{
  MSExperiment exp;

  MSChromatogram chrom;
  chrom.setNativeID("TIC");
  chrom.push_back(ChromatogramPeak(10.0, 100.0));
  chrom.push_back(ChromatogramPeak(20.0, 200.0));
  exp.addChromatogram(chrom);

  ArrowChromatogramExportConfig config;
  config.format = ArrowExportFormat::SemiWide;

  auto table = exportChromatogramsViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  TEST_EQUAL(table->num_rows(), 1) // One row per chromatogram

  auto schema = table->schema();
  auto rt_field = schema->GetFieldByName("rt");
  TEST_EQUAL(rt_field->type()->id(), arrow::Type::LIST)
}
END_SECTION

START_SECTION(exportChromatogramsToArrowCDataInterface - with RT filter)
{
  MSExperiment exp;

  MSChromatogram chrom;
  chrom.setNativeID("TIC");
  chrom.push_back(ChromatogramPeak(10.0, 100.0));
  chrom.push_back(ChromatogramPeak(20.0, 200.0));
  chrom.push_back(ChromatogramPeak(30.0, 150.0));
  exp.addChromatogram(chrom);

  ArrowChromatogramExportConfig config;
  config.format = ArrowExportFormat::Long;
  config.min_rt = 15.0;
  config.max_rt = 25.0;

  auto table = exportChromatogramsViaInterface(exp, config);

  TEST_NOT_EQUAL(table, nullptr)
  // Only point at RT=20.0 should be included
  TEST_EQUAL(table->num_rows(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Config struct tests
/////////////////////////////////////////////////////////////

START_SECTION(ArrowSpectraExportConfig - default values)
{
  ArrowSpectraExportConfig config;

  TEST_EQUAL(config.format == ArrowExportFormat::Long, true)
  TEST_EQUAL(config.ms_levels.empty(), true)
  TEST_REAL_SIMILAR(config.min_rt, 0.0)
  TEST_REAL_SIMILAR(config.max_rt, 0.0)
  TEST_REAL_SIMILAR(config.min_mz, 0.0)
  TEST_REAL_SIMILAR(config.max_mz, 0.0)
  TEST_EQUAL(config.columns.empty(), true)
  TEST_EQUAL(config.include_precursor_info, true)
  TEST_EQUAL(config.include_ion_mobility, true)
}
END_SECTION

START_SECTION(ArrowChromatogramExportConfig - default values)
{
  ArrowChromatogramExportConfig config;

  TEST_EQUAL(config.format == ArrowExportFormat::Long, true)
  TEST_REAL_SIMILAR(config.min_rt, 0.0)
  TEST_REAL_SIMILAR(config.max_rt, 0.0)
  TEST_EQUAL(config.columns.empty(), true)
}
END_SECTION

START_SECTION(ParquetWriteConfig - default values)
{
  ParquetWriteConfig config;

  TEST_EQUAL(config.compression == ParquetWriteConfig::Compression::ZSTD, true)
  TEST_EQUAL(config.compression_level, 3)
  TEST_EQUAL(config.row_group_size, 128 * 1024 * 1024)
  TEST_EQUAL(config.write_statistics, true)
  TEST_EQUAL(config.data_page_size, 1024 * 1024)
}
END_SECTION

START_SECTION(exportSpectraToParquet - basic export)
{
  MSExperiment exp = createTestExperiment();

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  bool success = MSExperimentArrowExport::exportSpectraToParquet(exp, filename);
  TEST_EQUAL(success, true)

  // Verify file was created
  TEST_EQUAL(File::exists(filename), true)
  TEST_EQUAL(File::empty(filename), false)
}
END_SECTION

START_SECTION(exportSpectraToParquet - with compression options)
{
  MSExperiment exp = createTestExperiment();

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  // Test with different compression settings
  ParquetWriteConfig pq_config;
  pq_config.compression = ParquetWriteConfig::Compression::ZSTD;
  pq_config.compression_level = 9;
  pq_config.row_group_size = 1024 * 1024; // 1MB for testing

  bool success = MSExperimentArrowExport::exportSpectraToParquet(exp, filename, ArrowSpectraExportConfig{}, pq_config);
  TEST_EQUAL(success, true)

  // Verify file was created
  TEST_EQUAL(File::exists(filename), true)
}
END_SECTION

START_SECTION(exportSpectraToParquet - with filtering)
{
  MSExperiment exp = createTestExperiment();

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  // Export only MS2
  ArrowSpectraExportConfig config;
  config.ms_levels = {2};

  bool success = MSExperimentArrowExport::exportSpectraToParquet(exp, filename, config);
  TEST_EQUAL(success, true)

  // Verify file was created
  TEST_EQUAL(File::exists(filename), true)
}
END_SECTION

START_SECTION(exportSpectraToParquet - empty experiment)
{
  MSExperiment exp;

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  bool success = MSExperimentArrowExport::exportSpectraToParquet(exp, filename);
  TEST_EQUAL(success, true)

  // Verify file was created (even for empty data)
  TEST_EQUAL(File::exists(filename), true)
}
END_SECTION

START_SECTION(exportChromatogramsToParquet - basic export)
{
  MSExperiment exp;

  // Create SRM/MRM chromatograms with precursor and product m/z
  MSChromatogram chrom1;
  chrom1.setNativeID("TIC");
  chrom1.getPrecursor().setMZ(500.0);
  chrom1.getProduct().setMZ(200.0);
  chrom1.push_back(ChromatogramPeak(10.0, 100.0));
  chrom1.push_back(ChromatogramPeak(20.0, 200.0));
  chrom1.push_back(ChromatogramPeak(30.0, 150.0));
  exp.addChromatogram(chrom1);

  MSChromatogram chrom2;
  chrom2.setNativeID("XIC_1");
  chrom2.getPrecursor().setMZ(600.0);
  chrom2.getProduct().setMZ(300.0);
  chrom2.push_back(ChromatogramPeak(15.0, 50.0));
  chrom2.push_back(ChromatogramPeak(25.0, 75.0));
  exp.addChromatogram(chrom2);

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  bool success = MSExperimentArrowExport::exportChromatogramsToParquet(exp, filename);
  TEST_EQUAL(success, true)

  // Verify file was created
  TEST_EQUAL(File::exists(filename), true)
  TEST_EQUAL(File::empty(filename), false)
}
END_SECTION

START_SECTION(exportChromatogramsToParquet - with compression options)
{
  MSExperiment exp;

  // Create SRM/MRM chromatogram
  MSChromatogram chrom;
  chrom.setNativeID("SRM_transition_1");
  chrom.getPrecursor().setMZ(450.0);
  chrom.getProduct().setMZ(175.0);
  chrom.push_back(ChromatogramPeak(5.0, 500.0));
  chrom.push_back(ChromatogramPeak(10.0, 1000.0));
  chrom.push_back(ChromatogramPeak(15.0, 750.0));
  exp.addChromatogram(chrom);

  String filename;
  NEW_TMP_FILE(filename);
  filename += ".parquet";

  // Test with SNAPPY compression
  ParquetWriteConfig pq_config;
  pq_config.compression = ParquetWriteConfig::Compression::SNAPPY;

  bool success = MSExperimentArrowExport::exportChromatogramsToParquet(exp, filename, ArrowChromatogramExportConfig{}, pq_config);
  TEST_EQUAL(success, true)

  // Verify file was created
  TEST_EQUAL(File::exists(filename), true)
}
END_SECTION

END_TEST
