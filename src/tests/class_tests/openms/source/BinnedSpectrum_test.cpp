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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
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

START_SECTION(~BinnedSpectrum())
{
  delete ptr;
}
END_SECTION

BinnedSpectrum* bs1;
DTAFile dtafile;
PeakSpectrum s1;
DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);

START_SECTION((BinnedSpectrum(const PeakSpectrum & ps, float size, UInt spread, float offset)))
{
  bs1 = new BinnedSpectrum(s1, 1.5, false, 2, 0.0);
  TEST_NOT_EQUAL(bs1, nullPointer)
}
END_SECTION

START_SECTION((BinnedSpectrum(const BinnedSpectrum &source)))
{
  BinnedSpectrum copy(*bs1);
  TEST_EQUAL(copy.getBinSize(), bs1->getBinSize());
  TEST_EQUAL(copy.getPrecursors().size(), 1);
  TEST_EQUAL(bs1->getPrecursors().size(), 1);
  TEST_EQUAL((UInt)copy.getPrecursors()[0].getMZ(),(UInt)bs1->getPrecursors()[0].getMZ());
}
END_SECTION

START_SECTION((BinnedSpectrum& operator=(const BinnedSpectrum &source)))
{
  BinnedSpectrum copy(*bs1);
  bs1 = new BinnedSpectrum(s1, 1.5, false, 2, 0.0);
  TEST_EQUAL(copy.getBinSize(), bs1->getBinSize());
  TEST_EQUAL((UInt)copy.getPrecursors()[0].getMZ(),(UInt)bs1->getPrecursors()[0].getMZ());
}
END_SECTION

START_SECTION((bool operator==(const BinnedSpectrum &rhs) const ))
{
  BinnedSpectrum copy = *bs1;
  TEST_EQUAL((*bs1 == copy), true)
}
END_SECTION

START_SECTION((bool operator!=(const BinnedSpectrum &rhs) const ))
{
  BinnedSpectrum copy = *bs1;
  TEST_EQUAL((*bs1!=copy),false)
}
END_SECTION

START_SECTION((float getBinSize() const ))
{
  TEST_EQUAL(bs1->getBinSize(),1.5)
}
END_SECTION

START_SECTION((UInt getBinSpread() const ))
{
  TEST_EQUAL(bs1->getBinSpread(),2)
}
END_SECTION

START_SECTION((SparseVectorIndexType getBinIndex(double mz) const))
{
  bs1 = new BinnedSpectrum(s1, 10, true, 0, 0.0); // 10 ppm bins
  TEST_EQUAL(bs1->getBinIndex(1.0), 0);
  TEST_EQUAL(bs1->getBinIndex(10.0), 230259);
  TEST_EQUAL(bs1->getBinIndex(100.0), 460519);
  TEST_EQUAL(bs1->getBinIndex(1000.0), 690778);
}
END_SECTION

START_SECTION((float getBinLowerMZ(size_t i) const))
{
  bs1 = new BinnedSpectrum(s1, 10, true, 0, 0); // 10 ppm bins
  TEST_REAL_SIMILAR(bs1->getBinLowerMZ(0), 1); // m/z = 1 corresponds to lowest index
  TEST_REAL_SIMILAR(bs1->getBinLowerMZ(1), 1 + 10 * 1e-6); // (1 + 10 ppm)
  TEST_REAL_SIMILAR(bs1->getBinLowerMZ(1000), pow(1 + 10 * 1e-6, 1000));
  TEST_REAL_SIMILAR(bs1->getBinLowerMZ(bs1->getBinIndex(1.0)), 1);
  TEST_REAL_SIMILAR(bs1->getBinLowerMZ(bs1->getBinIndex(10.0)), 10.0);
  TEST_REAL_SIMILAR(bs1->getBinLowerMZ(bs1->getBinIndex(100.0)), 100.0);
  TEST_REAL_SIMILAR(bs1->getBinLowerMZ(bs1->getBinIndex(1000.0)), 1000.0);

  BinnedSpectrum* bs2 = new BinnedSpectrum(s1, 1.0, false, 0, 0.5); // 1.0 m/z bins with 0.5 offset
  // offset ensures that floats close to nominal masses fall into same bin 
  TEST_EQUAL(bs2->getBinIndex(999.99), bs2->getBinIndex(1000.01));
  TEST_EQUAL(bs2->getBinIndex(99.99), bs2->getBinIndex(100.01));
  TEST_EQUAL(bs2->getBinIndex(9.99), bs2->getBinIndex(10.01));
  TEST_EQUAL(bs2->getBinIndex(0.99), bs2->getBinIndex(1.01));
  TEST_EQUAL(bs2->getBinIndex(0), bs2->getBinIndex(0.01));
  TEST_REAL_SIMILAR(bs2->getBinLowerMZ(0), -0.5); // because of offset, bin starts at -0.5
  TEST_REAL_SIMILAR(bs2->getBinLowerMZ(1), 0.5);
  TEST_REAL_SIMILAR(bs2->getBinLowerMZ(1000), 999.5);
  TEST_REAL_SIMILAR(bs2->getBinLowerMZ(bs2->getBinIndex(0.5)), 0.5);
  TEST_REAL_SIMILAR(bs2->getBinLowerMZ(bs2->getBinIndex(9.5)), 9.5);
  TEST_REAL_SIMILAR(bs2->getBinLowerMZ(bs2->getBinIndex(99.5)), 99.5);
  TEST_REAL_SIMILAR(bs2->getBinLowerMZ(bs2->getBinIndex(999.5)), 999.5);
}
END_SECTION

START_SECTION((const SparseVectorType& getBins() const))
{
  bs1 = new BinnedSpectrum(s1, 1.5, false, 2, 0);
  // count non-zero elements before access
  TEST_EQUAL(bs1->getBins().nonZeros(), 347)

  // access by bin index
  TEST_EQUAL(bs1->getBins().coeffRef(658), 501645)

  // check if number of non-zero elements is still the same
  TEST_EQUAL(bs1->getBins().nonZeros(), 347)

  // some additional tests for the underlying Eigen SparseVector
  UInt c = 0;
  for (BinnedSpectrum::SparseVectorIteratorType it(bs1->getBins()); it; ++it) { ++c; }
  TEST_EQUAL(bs1->getBins().nonZeros(), c)
}
END_SECTION

START_SECTION((SparseVectorType& getBins()))
{
  TEST_EQUAL(bs1->getBins().coeffRef(658),501645)
}
END_SECTION

START_SECTION((void setBinning()))
{
  NOT_TESTABLE
  //tested within another test
}
END_SECTION

// static
START_SECTION((bool BinnedSpectrum::isCompatible(const BinnedSpectrum& a, const BinnedSpectrum& b)))
{
  BinnedSpectrum bs2(s1, 1.234, false, 2, 0.0);
  TEST_EQUAL(BinnedSpectrum::isCompatible(*bs1, bs2), false)
  TEST_EQUAL(BinnedSpectrum::isCompatible(*bs1, *bs1), true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


