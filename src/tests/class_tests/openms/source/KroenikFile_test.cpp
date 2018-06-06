// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#include <OpenMS/FORMAT/KroenikFile.h>
///////////////////////////

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;
using namespace std;

START_TEST(KroenikFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

KroenikFile* ptr = nullptr;
KroenikFile* null_ptr = nullptr;
START_SECTION(KroenikFile())
{
	ptr = new KroenikFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~KroenikFile())
{
	delete ptr;
}
END_SECTION


START_SECTION((template < typename FeatureMapType > void load(const String &filename, FeatureMapType &feature_map)))
{
  KroenikFile f;
  FeatureMap fm;
  f.load(OPENMS_GET_TEST_DATA_PATH("KroenikFile_test_1.krf"), fm);
  TEST_EQUAL(fm.size(),3)
  ABORT_IF(fm.size()!=3)
  TEST_EQUAL(fm[0].getRT(), 63.2)
  TEST_REAL_SIMILAR(fm[0].getMZ(), 1002.11)
  TEST_EQUAL(fm[0].getIntensity(), 999999)
  TEST_EQUAL(fm[0].getCharge(), 1)
  TEST_EQUAL(String(fm[0].getMetaValue("AveragineModifications")), String("Carbamido"))
  TEST_EQUAL(fm[1].getRT(),  62.2)
  TEST_REAL_SIMILAR(fm[1].getMZ(), 252.057	)
  TEST_EQUAL(fm[1].getIntensity(), 9999)
  TEST_EQUAL(fm[1].getCharge(), 2)
  TEST_EQUAL(String(fm[1].getMetaValue("AveragineModifications")), String("Carbamido2"))

  

  TEST_EXCEPTION(Exception::ParseError, f.load(OPENMS_GET_TEST_DATA_PATH("KroenikFile_test_2.krf"), fm));
  
  TEST_EXCEPTION(Exception::FileNotFound, f.load(OPENMS_GET_TEST_DATA_PATH("KroenikFile_test_2_doesnotexist.edta"), fm));
}
END_SECTION

START_SECTION((template < typename SpectrumType > void store(const String &filename, const SpectrumType &spectrum) const ))
{
  KroenikFile f;
  MSSpectrum spec;
  TEST_EXCEPTION(Exception::NotImplemented, f.store("bla", spec))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



