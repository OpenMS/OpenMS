// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>

#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

using namespace OpenMS;
using namespace std;

START_TEST(ChromatogramExtractor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChromatogramExtractor* ptr = nullptr;
ChromatogramExtractor* nullPointer = nullptr;

START_SECTION(ChromatogramExtractor())
{
	ptr = new ChromatogramExtractor();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ChromatogramExtractor())
{
  delete ptr;
}
END_SECTION

START_SECTION(void extractChromatograms(const OpenSwath::SpectrumAccessPtr input, std::vector< OpenSwath::ChromatogramPtr > &output, std::vector< ExtractionCoordinates > extraction_coordinates, double mz_extraction_window, bool ppm, String filter))
{
  NOT_TESTABLE // is tested in ChromatogramExtractorAlgorithm
}
END_SECTION

START_SECTION(void prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms, std::vector< ExtractionCoordinates > & coordinates, OpenMS::TargetedExperiment & transition_exp, const double rt_extraction_window, const bool ms1) const)
{
  TargetedExperiment transitions;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.TraML"), transitions);
  TargetedExperiment transitions_;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.TraML"), transitions_);
  double rt_extraction_window = 1.0;

  // Test transitions
  {
    std::vector< OpenSwath::ChromatogramPtr > output_chromatograms;
    std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
    ChromatogramExtractor extractor;
    extractor.prepare_coordinates(output_chromatograms, coordinates, transitions, rt_extraction_window, false);

    TEST_TRUE(transitions == transitions_)
    TEST_EQUAL(output_chromatograms.size(), coordinates.size())
    TEST_EQUAL(coordinates.size(), 3)
    TEST_EQUAL(coordinates[0].mz, 618.31)
    TEST_EQUAL(coordinates[1].mz, 628.45)
    TEST_EQUAL(coordinates[2].mz, 654.38)

    TEST_REAL_SIMILAR(coordinates[0].rt_start, 1.5)
    TEST_REAL_SIMILAR(coordinates[1].rt_start, 43.5)
    TEST_REAL_SIMILAR(coordinates[2].rt_start, 43.5)

    TEST_REAL_SIMILAR(coordinates[0].rt_end, 2.5)
    TEST_REAL_SIMILAR(coordinates[1].rt_end, 44.5)
    TEST_REAL_SIMILAR(coordinates[2].rt_end, 44.5)

    // Note: they are ordered according to m/z
    TEST_EQUAL(coordinates[0].id, "tr3")
    TEST_EQUAL(coordinates[1].id, "tr1")
    TEST_EQUAL(coordinates[2].id, "tr2")
  }

  // Test peptides
  {
    std::vector< OpenSwath::ChromatogramPtr > output_chromatograms;
    std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
    ChromatogramExtractor extractor;
    extractor.prepare_coordinates(output_chromatograms, coordinates, transitions, rt_extraction_window, true);

    TEST_TRUE(transitions == transitions_)
    TEST_EQUAL(output_chromatograms.size(), coordinates.size())
    TEST_EQUAL(coordinates.size(), 2)
    TEST_EQUAL(coordinates[0].mz, 500)
    TEST_EQUAL(coordinates[1].mz, 501)

    TEST_REAL_SIMILAR(coordinates[0].rt_start, 43.5)
    TEST_REAL_SIMILAR(coordinates[1].rt_start, 1.5)

    TEST_REAL_SIMILAR(coordinates[0].rt_end, 44.5)
    TEST_REAL_SIMILAR(coordinates[1].rt_end, 2.5)

    TEST_EQUAL(coordinates[0].id, "tr_gr1_Precursor_i0")
    TEST_EQUAL(coordinates[1].id, "tr_gr2_Precursor_i0")

    TEST_EQUAL(OpenSwathHelper::computeTransitionGroupId(coordinates[0].id), "tr_gr1")
    TEST_EQUAL(OpenSwathHelper::computeTransitionGroupId(coordinates[1].id), "tr_gr2")
  }

}
END_SECTION

START_SECTION((template < typename TransitionExpT > static void return_chromatogram(std::vector< OpenSwath::ChromatogramPtr > &chromatograms, std::vector< ExtractionCoordinates > &coordinates, TransitionExpT &transition_exp_used, SpectrumSettings settings, std::vector< OpenMS::MSChromatogram > &output_chromatograms, bool ms1)))
{
  double extract_window = 0.05;
  double ppm = false;
  double rt_extraction_window = -1;
  String extraction_function = "tophat";

  TargetedExperiment transitions;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.TraML"), transitions);

  boost::shared_ptr<PeakMap > exp(new PeakMap);
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.mzML"), *exp);
  OpenSwath::SpectrumAccessPtr expptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);

  std::vector< OpenMS::MSChromatogram > chromatograms;
  {
    std::vector< OpenSwath::ChromatogramPtr > output_chromatograms;
    std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
    ChromatogramExtractor extractor;
    extractor.prepare_coordinates(output_chromatograms, coordinates, transitions, rt_extraction_window, false);

    extractor.extractChromatograms(expptr, output_chromatograms, coordinates, 
        extract_window, ppm, extraction_function);

    extractor.return_chromatogram(output_chromatograms, coordinates, transitions, (*exp)[0], chromatograms, false);
  }

  TEST_EQUAL(chromatograms.size(), 3)
  TEST_EQUAL(chromatograms[0].getChromatogramType(), ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
  TEST_REAL_SIMILAR(chromatograms[0].getProduct().getMZ(), 618.31)
  TEST_EQUAL(chromatograms[0].getPrecursor().metaValueExists("peptide_sequence"), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

