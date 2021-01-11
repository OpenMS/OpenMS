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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
///////////////////////////
#include <OpenMS/QC/FWHM.h>

#include <OpenMS/KERNEL/FeatureMap.h>

///////////////////////////

START_TEST(FWHM, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
using namespace OpenMS;
using namespace std;

FWHM* ptr = nullptr;
FWHM* nullPointer = nullptr;
START_SECTION(MzCalibration())
ptr = new FWHM();
TEST_NOT_EQUAL(ptr, nullPointer);
END_SECTION

START_SECTION(~FWHM())
delete ptr;
END_SECTION


START_SECTION(void compute(FeatureMap& features))
{
  Feature f;
  PeptideIdentification pi;
  pi.getHits().push_back(PeptideHit(1.0, 1, 3, AASequence::fromString("KKK")));
  f.getPeptideIdentifications().push_back(pi);
  f.setMetaValue("FWHM", 123.4);
  FeatureMap fm;
  fm.push_back(f);
  f.clearMetaInfo();
  f.setMetaValue("model_FWHM", 98.1);
  fm.push_back(f); 
  FWHM fw;
  fw.compute(fm);
  TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getMetaValue("FWHM"), 123.4)
  TEST_EQUAL(fm[1].getPeptideIdentifications()[0].getMetaValue("FWHM"), 98.1)
}
END_SECTION

START_SECTION(QCBase::Status requires() const override)
{
  FWHM fw;
  TEST_EQUAL(fw.requires() == (QCBase::Status() | QCBase::Requires::POSTFDRFEAT), true);
}
END_SECTION

START_SECTION(const String & getName() const)
{
  TEST_EQUAL(FWHM().getName(), "FWHM");
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
