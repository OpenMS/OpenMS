// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <fstream>
#include <QtCore/QDir>
#include <QtCore/QString>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SiriusFragmentAnnotation, "$Id$")

/////////////////////////////////////////////////////////////

SiriusFragmentAnnotation* fa_ptr = nullptr;
SiriusFragmentAnnotation* fa_null = nullptr;

START_SECTION(SiriusFragmentAnnotation())
{
    fa_ptr = new SiriusFragmentAnnotation;
    TEST_NOT_EQUAL(fa_ptr, fa_null);
}
END_SECTION

START_SECTION(~SiriusFragmentAnnotation())
{
    delete fa_ptr;
}
END_SECTION

// mz  intensity   rel.intensity   exactmass   explanation
// 123.000363  65857.38    3.36    122.999604  C7H3Cl

// test function 
START_SECTION(static void extractSiriusFragmentAnnotationMapping(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill, bool use_exact_mass))
{
    String test_path = OPENMS_GET_TEST_DATA_PATH("SiriusFragmentAnnotation_test");
    MSSpectrum annotated_msspectrum = SiriusFragmentAnnotation::extractAnnotationsFromSiriusFile(test_path, 1, false, false)[0];

    TEST_STRING_SIMILAR(annotated_msspectrum.getNativeID(), "sample=1 period=1 cycle=676 experiment=4|sample=1 period=1 cycle=677 experiment=5|sample=1 period=1 cycle=678 experiment=3");
    TEST_EQUAL(annotated_msspectrum.getMSLevel(), 2);

    TEST_EQUAL(annotated_msspectrum.empty(), false);
    TEST_REAL_SIMILAR(annotated_msspectrum[0].getMZ(), 70.040098);
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("peak_mz"), "mz");
    TEST_STRING_SIMILAR(annotated_msspectrum.getFloatDataArrays()[0].getName(), "exact_mass");
    TEST_REAL_SIMILAR(annotated_msspectrum.getFloatDataArrays()[0][0], 70.040098);
    TEST_STRING_SIMILAR(annotated_msspectrum.getStringDataArrays()[0][0], "C2H3N3");
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("annotated_sumformula"), "C15H17ClN4");
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("annotated_adduct"), "[M+H]+");
    TEST_REAL_SIMILAR(annotated_msspectrum.getMetaValue("decoy"), 0);
}
END_SECTION

// test exact mass output 
START_SECTION(static void extractSiriusFragmentAnnotationMapping(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill, bool use_exact_mass))
{
    String test_path = OPENMS_GET_TEST_DATA_PATH("SiriusFragmentAnnotation_test");
    MSSpectrum annotated_msspectrum = SiriusFragmentAnnotation::extractAnnotationsFromSiriusFile(test_path, 1, false, true)[0];

    TEST_STRING_SIMILAR(annotated_msspectrum.getNativeID(), "sample=1 period=1 cycle=676 experiment=4|sample=1 period=1 cycle=677 experiment=5|sample=1 period=1 cycle=678 experiment=3");
    TEST_EQUAL(annotated_msspectrum.getMSLevel(), 2);

    TEST_EQUAL(annotated_msspectrum.empty(), false);
    TEST_REAL_SIMILAR(annotated_msspectrum[0].getMZ(), 70.040098)
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("peak_mz"), "exact_mass");
    TEST_STRING_SIMILAR(annotated_msspectrum.getFloatDataArrays()[0].getName(), "mz");
    TEST_REAL_SIMILAR(annotated_msspectrum.getFloatDataArrays()[0][0], 70.040098);
    TEST_STRING_SIMILAR(annotated_msspectrum.getStringDataArrays()[0][0], "C2H3N3");
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("annotated_sumformula"), "C15H17ClN4");
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("annotated_adduct"), "[M+H]+");
    TEST_REAL_SIMILAR(annotated_msspectrum.getMetaValue("decoy"), 0);
}
END_SECTION

// test decoy extraction
START_SECTION(static void extractSiriusDecoyAnnotationMapping(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill))
{
    String test_path = OPENMS_GET_TEST_DATA_PATH("SiriusFragmentAnnotation_test");
    MSSpectrum decoy_msspectrum = SiriusFragmentAnnotation::extractAnnotationsFromSiriusFile(test_path, 1, true, false)[0];

    TEST_STRING_SIMILAR(decoy_msspectrum.getNativeID(), "sample=1 period=1 cycle=676 experiment=4|sample=1 period=1 cycle=677 experiment=5|sample=1 period=1 cycle=678 experiment=3");
    TEST_EQUAL(decoy_msspectrum.getMSLevel(), 2);

    TEST_EQUAL(decoy_msspectrum.empty(), false);
    TEST_REAL_SIMILAR(decoy_msspectrum[0].getMZ(), 53.013424);
    TEST_STRING_SIMILAR(decoy_msspectrum.getMetaValue("peak_mz"), "mz");
    TEST_STRING_SIMILAR(decoy_msspectrum.getStringDataArrays()[0][0], "C2N2");
    TEST_STRING_SIMILAR(decoy_msspectrum.getMetaValue("annotated_sumformula"), "C15H17ClN4");
    TEST_STRING_SIMILAR(decoy_msspectrum.getMetaValue("annotated_adduct"), "[M+H]+");
    TEST_REAL_SIMILAR(decoy_msspectrum.getMetaValue("decoy"), 1);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
