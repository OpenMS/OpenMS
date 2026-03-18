// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Andreas Bertsch, Lars Nilse $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SpectraMerger, "$Id$")

/////////////////////////////////////////////////////////////

SpectraMerger* e_ptr = nullptr;
SpectraMerger* e_nullPointer = nullptr;
START_SECTION((SpectraMerger()))
	e_ptr = new SpectraMerger;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~SpectraMerger()))
	delete e_ptr;
END_SECTION

e_ptr = new SpectraMerger();

START_SECTION((SpectraMerger(const SpectraMerger& source)))
	SpectraMerger copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
END_SECTION

START_SECTION((SpectraMerger& operator=(const SpectraMerger& source)))
	SpectraMerger copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
END_SECTION

START_SECTION((template < typename MapType > void mergeSpectraBlockWise(MapType &exp)))
	PeakMap exp, exp2;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SpectraMerger_input_2.mzML"), exp);
	TEST_EQUAL(exp.size(), 144)

  exp2 = exp;

	SpectraMerger merger;
  Param p;
  p.setValue("mz_binning_width", 0.0001, "Max m/z distance of two peaks to be merged.", {"advanced"});
  p.setValue("mz_binning_width_unit", "Da", "Unit in which the distance between two peaks is given.", {"advanced"});

  p.setValue("block_method:rt_block_size", 5);
  p.setValue("block_method:ms_levels", ListUtils::create<Int>("1"));
  merger.setParameters(p);
  merger.mergeSpectraBlockWise(exp);
  TEST_EQUAL(exp.size(), 130);
  exp=exp2;

  p.setValue("block_method:rt_block_size", 4);
  p.setValue("block_method:ms_levels", ListUtils::create<Int>("2"));
  merger.setParameters(p);
  merger.mergeSpectraBlockWise(exp);
  TEST_EQUAL(exp.size(), 50);
  TEST_REAL_SIMILAR(exp[0].getRT(),201.0275)
  TEST_REAL_SIMILAR(exp[1].getRT(),204.34075)
  TEST_EQUAL(exp[1].getMSLevel(), 2);
  TEST_EQUAL(exp[2].getMSLevel(), 1);
  exp=exp2;

  p.setValue("block_method:rt_block_size", 4);
  p.setValue("block_method:ms_levels", ListUtils::create<Int>("1,2"));
  merger.setParameters(p);
  merger.mergeSpectraBlockWise(exp);
  TEST_EQUAL(exp.size(), 37);

END_SECTION

START_SECTION((template < typename MapType > void mergeSpectraPrecursors(MapType &exp)))
	PeakMap exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SpectraMerger_input_precursor.mzML"), exp);

	SpectraMerger merger;
	TEST_EQUAL(exp.size(), 17)

  Param p;
  p.setValue("mz_binning_width", 0.3, "Max m/z distance of two peaks to be merged.", {"advanced"});

  p.setValue("mz_binning_width_unit", "Da", "Unit in which the distance between two peaks is given.", {"advanced"});

  // same precursor MS/MS merging
 	p.setValue("precursor_method:mz_tolerance", 10e-5, "Max m/z distance of the precursor entries of two spectra to be merged in [Da].");
  p.setValue("precursor_method:rt_tolerance", 5.0, "Max RT distance of the precursor entries of two spectra to be merged in [s].");
  merger.setParameters(p);
  merger.mergeSpectraPrecursors(exp);

  PeakMap exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SpectraMerger_output_precursor.mzML"), exp2);

	TEST_EQUAL(exp.size(), exp2.size());
  ABORT_IF(exp.size() != exp2.size());

  for (Size i=0;i<exp.size();++i)
  {
    TEST_EQUAL(exp[i].size(), exp2[i].size())
    TEST_EQUAL(exp[i].getMSLevel (), exp2[i].getMSLevel ())
  }

END_SECTION

START_SECTION((bool areMassesMatched(double mz1, double mz2, double tol_ppm, int max_c)))
  SpectraMerger merger;
  bool non_matched = merger.areMassesMatched(100, 1000, 10, 5);
  bool matched = merger.areMassesMatched(1000, 1000.001, 10, 5);

  TEST_EQUAL(non_matched, false);
  TEST_EQUAL(matched, true);
END_SECTION

START_SECTION((template < typename MapType > void averageGaussian(MapType &exp)))
	PeakMap exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SpectraMerger_input_3.mzML"), exp);    // profile mode

	SpectraMerger merger;
	TEST_EQUAL(exp.size(), 28)

  Param p;
  p.setValue("mz_binning_width", 0.0001, "Max m/z distance of two peaks to be merged.", {"advanced"});
  p.setValue("mz_binning_width_unit", "Da", "Unit in which the distance between two peaks is given.", {"advanced"});

  // same precursor MS/MS merging
 	p.setValue("average_gaussian:spectrum_type", "automatic", "Spectrum type of the MS level to be averaged");
 	p.setValue("average_gaussian:ms_level", 1, "Average spectra of this level. All other spectra remain unchanged.");
  p.setValue("average_gaussian:rt_FWHM", 5.0, "FWHM of Gauss curve in seconds to be averaged over.");
  p.setValue("average_gaussian:cutoff", 0.01, "Intensity cutoff for Gaussian. The Gaussian RT profile decreases from 1 at its apex to 0 at infinity. Spectra for which the intensity of the Gaussian drops below the cutoff do not contribute to the average.", {"advanced"});
  merger.setParameters(p);
  merger.average(exp,"gaussian");

  PeakMap exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SpectraMerger_output_3.mzML"), exp2);

	TEST_EQUAL(exp.size(), exp2.size());
  ABORT_IF(exp.size() != exp2.size());

  for (Size i=0;i<exp.size();++i)
  {
    TEST_EQUAL(exp[i].size(), exp2[i].size())
    TEST_EQUAL(exp[i].getMSLevel (), exp2[i].getMSLevel ())
  }

END_SECTION

START_SECTION((template < typename MapType > void averageGaussian(MapType &exp)))
	PeakMap exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SpectraMerger_input_4.mzML"), exp);    // centroid mode

	SpectraMerger merger;
	TEST_EQUAL(exp.size(), 28)

  Param p;
  p.setValue("mz_binning_width", 0.0001, "Max m/z distance of two peaks to be merged.", {"advanced"});
  p.setValue("mz_binning_width_unit", "Da", "Unit in which the distance between two peaks is given.", {"advanced"});

  // same precursor MS/MS merging
 	p.setValue("average_gaussian:spectrum_type", "automatic", "Spectrum type of the MS level to be averaged");
 	p.setValue("average_gaussian:ms_level", 1, "Average spectra of this level. All other spectra remain unchanged.");
  p.setValue("average_gaussian:rt_FWHM", 5.0, "FWHM of Gauss curve in seconds to be averaged over.");
  p.setValue("average_gaussian:cutoff", 0.01, "Intensity cutoff for Gaussian. The Gaussian RT profile decreases from 1 at its apex to 0 at infinity. Spectra for which the intensity of the Gaussian drops below the cutoff do not contribute to the average.", {"advanced"});
  merger.setParameters(p);
  merger.average(exp,"gaussian");

  PeakMap exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SpectraMerger_output_4.mzML"), exp2);

	TEST_EQUAL(exp.size(), exp2.size());
  ABORT_IF(exp.size() != exp2.size());

  for (Size i=0;i<exp.size();++i)
  {
    TEST_EQUAL(exp[i].size(), exp2[i].size())
    TEST_EQUAL(exp[i].getMSLevel (), exp2[i].getMSLevel ())
  }

END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
