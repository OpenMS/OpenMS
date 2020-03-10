// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    MSSpectrum annotated_msspectrum;

    SiriusFragmentAnnotation::extractSiriusFragmentAnnotationMapping(test_path, annotated_msspectrum);
    TEST_STRING_SIMILAR(annotated_msspectrum.getNativeID(), "sample=1 period=1 cycle=2056 experiment=3");
    TEST_EQUAL(annotated_msspectrum.getMSLevel(), 2);

    TEST_EQUAL(annotated_msspectrum.empty(), false);
    TEST_REAL_SIMILAR(annotated_msspectrum[0].getMZ(), 123.000363);
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("peak_mz"), "mz");
    TEST_STRING_SIMILAR(annotated_msspectrum.getFloatDataArrays()[0].getName(), "exact_mass");
    TEST_REAL_SIMILAR(annotated_msspectrum.getFloatDataArrays()[0][0], 122.999604);
    TEST_STRING_SIMILAR(annotated_msspectrum.getStringDataArrays()[0][0], "C7H3Cl");
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("annotated_sumformula"), "C12H11Cl2N3O2");
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("annotated_adduct"), "M+H+");
}
END_SECTION

// test exact mass output 
START_SECTION(static void extractSiriusFragmentAnnotationMapping(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill, bool use_exact_mass))
{
    String test_path = OPENMS_GET_TEST_DATA_PATH("SiriusFragmentAnnotation_test");
    MSSpectrum annotated_msspectrum;

    SiriusFragmentAnnotation::extractSiriusFragmentAnnotationMapping(test_path, annotated_msspectrum, true);    
    TEST_STRING_SIMILAR(annotated_msspectrum.getNativeID(), "sample=1 period=1 cycle=2056 experiment=3");
    TEST_EQUAL(annotated_msspectrum.getMSLevel(), 2);

    TEST_EQUAL(annotated_msspectrum.empty(), false);
    TEST_REAL_SIMILAR(annotated_msspectrum[0].getMZ(), 122.999604)
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("peak_mz"), "exact_mass");
    TEST_STRING_SIMILAR(annotated_msspectrum.getFloatDataArrays()[0].getName(), "mz");
    TEST_REAL_SIMILAR(annotated_msspectrum.getFloatDataArrays()[0][0], 123.000363);
    TEST_STRING_SIMILAR(annotated_msspectrum.getStringDataArrays()[0][0], "C7H3Cl");
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("annotated_sumformula"), "C12H11Cl2N3O2");
    TEST_STRING_SIMILAR(annotated_msspectrum.getMetaValue("annotated_adduct"), "M+H+");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
