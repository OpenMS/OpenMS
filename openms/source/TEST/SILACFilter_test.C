// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SILACFilter, "$Id$")

std::vector<DoubleReal> mass_separations;
mass_separations.push_back(4);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((SILACFilter(std::vector< DoubleReal > mass_separations, Int charge, DoubleReal model_deviation, Int isotopes_per_peptide, DoubleReal intensity_cutoff, DoubleReal intensity_correlation, bool allow_missing_peaks)))
{
  SILACFilter f(mass_separations, 1, 1, 3, 0, 0, false);
  TEST_EQUAL(f.getCharge(), 1);
}
END_SECTION

START_SECTION((std::vector<DoubleReal> getPeakPositions()))
{
  SILACFilter f(mass_separations, 1, 1, 3, 0, 0, false);
  // XXX: Segfaults
  // TEST_EQUAL(f.getPeakPositions().size(), 0);
}
END_SECTION

START_SECTION((const std::vector<DoubleReal>& getExpectedMzShifts()))
{
  const UInt peaks_per_peptide = 3;
  SILACFilter f(mass_separations, 1, 1, peaks_per_peptide, 0, 0, false);
  TEST_EQUAL(f.getExpectedMzShifts().size(), (mass_separations.size() + 1) * peaks_per_peptide);
}
END_SECTION

START_SECTION((std::vector<SILACPattern>& getElements()))
{
  SILACFilter f(mass_separations, 1, 1, 3, 0, 0, false);
  TEST_EQUAL(f.getElements().size(), 0);
}
END_SECTION

START_SECTION((Int getCharge()))
{
  SILACFilter f(mass_separations, 1, 1, 3, 0, 0, false);
  TEST_EQUAL(f.getCharge(), 1);
}
END_SECTION

START_SECTION((std::vector<DoubleReal>& getMassSeparations()))
{
  SILACFilter f(mass_separations, 1, 1, 3, 0, 0, false);
  TEST_EQUAL(f.getMassSeparations() == mass_separations, true);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

