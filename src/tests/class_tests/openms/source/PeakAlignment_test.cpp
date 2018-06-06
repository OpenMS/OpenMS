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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/COMPARISON/SPECTRA/PeakAlignment.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakAlignment* ptr = nullptr;
PeakAlignment* nullPointer = nullptr;
START_SECTION(PeakAlignment())
{
	ptr = new PeakAlignment();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~PeakAlignment())
{
	delete ptr;
}
END_SECTION

START_SECTION((PeakAlignment(const PeakAlignment &source)))
{
	ptr = new PeakAlignment();
	PeakAlignment copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((PeakAlignment& operator=(const PeakAlignment &source)))
{
	PeakAlignment copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((double operator()(const PeakSpectrum &spec1, const PeakSpectrum &spec2) const ))
{
	PeakAlignment pa;
	PeakSpectrum s1, s2;
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s2);
	s2.pop_back();
	double score = pa(s1, s2);
	TEST_REAL_SIMILAR(score, 0.997477)

  // Test empty spectra - they should return zero
  PeakSpectrum empty_spectrum;
	score = pa(empty_spectrum, s2);
	TEST_REAL_SIMILAR(score, 0.0)

	score = pa(s1, empty_spectrum);
	TEST_REAL_SIMILAR(score, 0.0)
}
END_SECTION

START_SECTION((double operator()(const PeakSpectrum &spec) const ))
{
  PeakSpectrum s1;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
  double score = (*ptr)(s1);
  TEST_REAL_SIMILAR(score, 1);
}
END_SECTION

START_SECTION((vector< pair<Size,Size> > getAlignmentTraceback(const PeakSpectrum &spec1, const PeakSpectrum &spec2) const ))
{
	PeakAlignment pa;
	PeakSpectrum s1, s2;
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s2);
	vector< pair<Size,Size> > result, tester;
	result = pa.getAlignmentTraceback(s1,s2);
	for (Size i = 0; i < 127; ++i)
	{
		tester.push_back(pair<Size,Size>(i,i));
	}
	TEST_EQUAL(tester.size(),result.size())
	for (Size i = 0; i < tester.size(); ++i)
	{
		TEST_EQUAL(tester.at(i).first,result.at(i).first)
	}
}
END_SECTION

START_SECTION((static PeakSpectrumCompareFunctor* create()))
{
	PeakSpectrumCompareFunctor* psf = PeakAlignment::create();
	PeakAlignment pa;
	TEST_EQUAL(psf->getParameters(), pa.getParameters())
	TEST_EQUAL(psf->getName(), pa.getName())
}
END_SECTION

START_SECTION((static const String getProductName()))
{
	TEST_EQUAL(ptr->getProductName(), "PeakAlignment")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



