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
#include <OpenMS/CHEMISTRY/WeightWrapper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(WeightWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

WeightWrapper* ptr = nullptr;
WeightWrapper* null_ptr = nullptr;
START_SECTION(WeightWrapper())
{
  ptr = new WeightWrapper();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~WeightWrapper())
{
  delete ptr;
}
END_SECTION

START_SECTION((WeightWrapper(const WEIGHTMODE weight_mode)))
{
  WeightWrapper ww(WeightWrapper::MONO);
  TEST_EQUAL(ww.getWeightMode(), WeightWrapper::MONO)
  WeightWrapper ww2(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww2.getWeightMode(), WeightWrapper::AVERAGE)
}
END_SECTION

START_SECTION((WeightWrapper(const WeightWrapper &source)))
{
  WeightWrapper ww(WeightWrapper::AVERAGE);
  WeightWrapper ww2(ww);

  TEST_EQUAL(ww.getWeightMode(), ww2.getWeightMode())
}
END_SECTION

START_SECTION((void setWeightMode(const WEIGHTMODE mode)))
{
  WeightWrapper ww;
  TEST_EXCEPTION(Exception::IllegalArgument, ww.setWeightMode(WeightWrapper::SIZE_OF_WEIGHTMODE))
  ww.setWeightMode(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww.getWeightMode(), WeightWrapper::AVERAGE)
}
END_SECTION

START_SECTION((WEIGHTMODE getWeightMode() const ))
{
  WeightWrapper ww;
  TEST_EQUAL(ww.getWeightMode(), WeightWrapper::MONO)
}
END_SECTION

START_SECTION((double getWeight(const AASequence &aa) const ))
{
  WeightWrapper ww;
  AASequence aa= AASequence::fromString("DFINAGER");
  TEST_EQUAL(ww.getWeight(aa), aa.getMonoWeight())
  WeightWrapper ww2(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww2.getWeight(aa), aa.getAverageWeight())
}
END_SECTION

START_SECTION((double getWeight(const EmpiricalFormula &ef) const ))
{
  WeightWrapper ww;
  EmpiricalFormula aa("C12H544");
  TEST_EQUAL(ww.getWeight(aa), aa.getMonoWeight())
  WeightWrapper ww2(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww2.getWeight(aa), aa.getAverageWeight())
}
END_SECTION

START_SECTION((double getWeight(const Residue &r, Residue::ResidueType res_type=Residue::Full) const ))
{
  WeightWrapper ww;
  Residue aa("L", "LEU", "L", EmpiricalFormula("C454H33"));
  TEST_EQUAL(ww.getWeight(aa), aa.getMonoWeight())
  WeightWrapper ww2(WeightWrapper::AVERAGE);
  TEST_EQUAL(ww2.getWeight(aa), aa.getAverageWeight())
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



