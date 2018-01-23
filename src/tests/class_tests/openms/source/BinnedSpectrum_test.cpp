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
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/FORMAT/DTAFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BinnedSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BinnedSpectrum* ptr = nullptr;
BinnedSpectrum* nullPointer = nullptr;
START_SECTION(BinnedSpectrum())
{
	ptr = new BinnedSpectrum();
	TEST_NOT_EQUAL(ptr, nullPointer)
  TEST_EXCEPTION(BinnedSpectrum::NoSpectrumIntegrated,ptr->setBinning();)
}
END_SECTION

START_SECTION(~BinnedSpectrum())
{
	delete ptr;
}
END_SECTION

  BinnedSpectrum* bs1;
  DTAFile dtafile;
  PeakSpectrum s1;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);

START_SECTION((BinnedSpectrum(float size, UInt spread, PeakSpectrum ps)))
{
  bs1 = new BinnedSpectrum(1.5,2,s1);
  TEST_NOT_EQUAL(bs1,nullPointer)
}
END_SECTION

START_SECTION((BinnedSpectrum(const BinnedSpectrum &source)))
{
  BinnedSpectrum copy(*bs1);
  TEST_EQUAL(copy.getRawSpectrum().getName(), bs1->getRawSpectrum().getName());
  TEST_EQUAL(copy.getBinSize(), bs1->getBinSize());
  TEST_EQUAL((UInt)copy.getRawSpectrum().getPrecursors()[0].getMZ(),(UInt)bs1->getRawSpectrum().getPrecursors()[0].getMZ());
}
END_SECTION

START_SECTION((BinnedSpectrum& operator=(const BinnedSpectrum &source)))
{
  BinnedSpectrum copy(*bs1);
  bs1 = new BinnedSpectrum(1.5,2,s1);
  TEST_EQUAL(copy.getRawSpectrum().getName(), bs1->getRawSpectrum().getName());
  TEST_EQUAL(copy.getBinSize(), bs1->getBinSize());
  TEST_EQUAL((UInt)copy.getRawSpectrum().getPrecursors()[0].getMZ(),(UInt)bs1->getRawSpectrum().getPrecursors()[0].getMZ());
}
END_SECTION

START_SECTION((BinnedSpectrum& operator=(const PeakSpectrum &source)))
{
  bs1 = new BinnedSpectrum();
  *bs1 = s1;
  TEST_EQUAL(bs1->getRawSpectrum().getPrecursors()[0].getMZ(),s1.getPrecursors()[0].getMZ());
  bs1->setBinSize(1.5);
  bs1->setBinSpread(2);
}
END_SECTION

START_SECTION((const PeakSpectrum& getRawSpectrum() const))
{
  bs1 = new BinnedSpectrum(1.5,2,s1);
  TEST_EQUAL(bs1->getRawSpectrum(),s1)
}
END_SECTION

START_SECTION((bool operator==(const BinnedSpectrum &rhs) const ))
{
	BinnedSpectrum copy = *bs1;
	TEST_EQUAL((*bs1==copy),true)
}
END_SECTION

START_SECTION((bool operator!=(const BinnedSpectrum &rhs) const ))
{
	BinnedSpectrum copy = *bs1;
	TEST_EQUAL((*bs1!=copy),false)
}
END_SECTION

START_SECTION((bool operator==(const PeakSpectrum &rhs) const ))
{
	TEST_EQUAL((*bs1==s1),true)
}
END_SECTION

START_SECTION((bool operator!=(const PeakSpectrum &rhs) const ))
{
	BinnedSpectrum copy = *bs1;
	TEST_EQUAL((*bs1!=s1),false)
}
END_SECTION

START_SECTION((double getBinSize() const ))
{
	TEST_EQUAL(bs1->getBinSize(),1.5)
}
END_SECTION

START_SECTION((UInt getBinSpread() const ))
{
	TEST_EQUAL(bs1->getBinSpread(),2)
}
END_SECTION

START_SECTION((UInt getBinNumber() const ))
{
	TEST_EQUAL(bs1->getBinNumber(),659)
}
END_SECTION

START_SECTION((UInt getFilledBinNumber() const ))
{
	TEST_EQUAL(bs1->getFilledBinNumber(),347)
}
END_SECTION

START_SECTION((const SparseVector<float>& getBins() const))
{
	TEST_EQUAL(bs1->getBins().at(658),501645)
}
END_SECTION

START_SECTION((SparseVector<float>& getBins()))
{
	TEST_EQUAL(bs1->getBins().at(658),501645)
}
END_SECTION

START_SECTION((const_bin_iterator begin() const ))
{
	UInt c(0);
	for (BinnedSpectrum::const_bin_iterator it1 = bs1->begin(); it1 != bs1->end(); ++it1)
	{
		++c;
	}
	TEST_EQUAL(bs1->getBinNumber(),c)
}
END_SECTION

START_SECTION((const_bin_iterator end() const ))
{
	NOT_TESTABLE
	//tested above
}
END_SECTION

START_SECTION((bin_iterator begin()))
{
	UInt c(0);
	for (BinnedSpectrum::bin_iterator it1 = bs1->begin(); it1 != bs1->end(); ++it1)
	{
		++c;
	}
	TEST_EQUAL(bs1->getBinNumber(),c)
}
END_SECTION

START_SECTION((bin_iterator end()))
{
  NOT_TESTABLE
	//tested above
}
END_SECTION

START_SECTION((void setBinSize(double s)))
{
	TEST_EQUAL(bs1->getBinSize(),1.5)
}
END_SECTION

START_SECTION((void setBinSpread(UInt s)))
{
	TEST_EQUAL(bs1->getBinSpread(),2)
}
END_SECTION

START_SECTION((void setBinning()))
{
	NOT_TESTABLE
	//tested within another test
}
END_SECTION

START_SECTION((bool checkCompliance(const BinnedSpectrum &bs) const ))
{
	TEST_EQUAL(bs1->checkCompliance(BinnedSpectrum()),false)
}
END_SECTION

START_SECTION(([BinnedSpectrum::NoSpectrumIntegrated] NoSpectrumIntegrated(const char *file, int line, const char *function, const char *message="BinnedSpectrum hasn't got a PeakSpectrum to base on yet")))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([BinnedSpectrum::NoSpectrumIntegrated] virtual ~NoSpectrumIntegrated()))
{
  NOT_TESTABLE
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



