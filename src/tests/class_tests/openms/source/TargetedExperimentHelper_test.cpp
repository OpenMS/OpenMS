// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hannes Roest$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;
using namespace OpenMS::TargetedExperimentHelper;

START_TEST(TargetedExperimentHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(TargetedExperimentHelper::Configuration())
{
  // Ensure that Configuration has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Configuration(std::declval<TargetedExperimentHelper::Configuration&&>())), true)

  auto ptr = new TargetedExperimentHelper::Configuration();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::CV())
{
  // Ensure that CV has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::CV(std::declval<TargetedExperimentHelper::CV&&>())), true)

  auto ptr = new TargetedExperimentHelper::CV("", "", "", "");
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Protein())
{
  // Ensure that Protein has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Protein(std::declval<TargetedExperimentHelper::Protein&&>())), true)

  auto ptr = new TargetedExperimentHelper::Protein();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::RetentionTime())
{
  // Ensure that RetentionTime has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::RetentionTime(std::declval<TargetedExperimentHelper::RetentionTime&&>())), true)

  auto ptr = new TargetedExperimentHelper::RetentionTime();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::PeptideCompound())
{
  // Ensure that PeptideCompound has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::PeptideCompound(std::declval<TargetedExperimentHelper::PeptideCompound&&>())), true)

  auto ptr = new TargetedExperimentHelper::PeptideCompound();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Peptide())
{
  // Ensure that Peptide has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Peptide(std::declval<TargetedExperimentHelper::Peptide&&>())), true)

  auto ptr = new TargetedExperimentHelper::Peptide();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Compound())
{
  // Ensure that Compound has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Compound(std::declval<TargetedExperimentHelper::Compound&&>())), true)

  auto ptr = new TargetedExperimentHelper::Compound();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Contact())
{
  // Ensure that Contact has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Contact(std::declval<TargetedExperimentHelper::Contact&&>())), true)

  auto ptr = new TargetedExperimentHelper::Contact();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Publication())
{
  // Ensure that Publication has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Publication(std::declval<TargetedExperimentHelper::Publication&&>())), true)

  auto ptr = new TargetedExperimentHelper::Publication();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Instrument())
{
  // Ensure that Instrument has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Instrument(std::declval<TargetedExperimentHelper::Instrument&&>())), true)

  auto ptr = new TargetedExperimentHelper::Instrument();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Prediction())
{
  // Ensure that Prediction has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Prediction(std::declval<TargetedExperimentHelper::Prediction&&>())), true)

  auto ptr = new TargetedExperimentHelper::Prediction();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Interpretation())
{
  // Ensure that Interpretation has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Interpretation(std::declval<TargetedExperimentHelper::Interpretation&&>())), true)

  auto ptr = new TargetedExperimentHelper::Interpretation();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::TraMLProduct())
{
  // Ensure that TraMLProduct has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::TraMLProduct(std::declval<TargetedExperimentHelper::TraMLProduct&&>())), true)

  auto ptr = new TargetedExperimentHelper::TraMLProduct();
  TEST_EQUAL(ptr != nullptr, true)
  delete ptr;
}
END_SECTION

TargetedExperimentHelper::Peptide* ptr = nullptr;
TargetedExperimentHelper::Peptide* null_ptr = nullptr;

START_SECTION(TargetedExperimentHelper::Peptide())
{
  ptr = new TargetedExperimentHelper::Peptide();
  TEST_NOT_EQUAL(ptr, null_ptr)

  // Ensure that Peptide has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(Peptide(std::declval<Peptide&&>())), true)
}
END_SECTION

START_SECTION(~TargetedExperimentHelper::Peptide())
{
  delete ptr;
}
END_SECTION

START_SECTION((TargetedExperiment::RetentionTime))
{
  TargetedExperimentHelper::RetentionTime rt;
    
  TEST_EQUAL(rt.isRTset(), false)

  rt.setRT(5.0);
  TEST_EQUAL(rt.isRTset(), true)
  TEST_REAL_SIMILAR (rt.getRT(), 5.0)

}
END_SECTION

START_SECTION((TargetedExperiment::Peptide))
{
  TargetedExperimentHelper::Peptide p;
    
  TEST_EQUAL(p.hasRetentionTime(), false)
  TEST_EQUAL(p.rts.size(), 0)

  // add a RT
  TargetedExperimentHelper::RetentionTime rt;
  rt.setRT(5.0);
  rt.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
  rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;
  p.rts.push_back(rt);

  // test the RT methods
  TEST_EQUAL(p.rts.size(), 1)
  TEST_EQUAL(p.rts[0] == rt, true)
  TEST_EQUAL(p.rts[0].retention_time_unit == TargetedExperimentHelper::RetentionTime::RTUnit::SECOND, true)
  TEST_EQUAL(p.rts[0].retention_time_type == TargetedExperimentHelper::RetentionTime::RTType::PREDICTED, true)
  TEST_REAL_SIMILAR(p.rts[0].getRT(), 5.0)

  // test the Peptide methods
  TEST_EQUAL(p.hasRetentionTime(), true)
  TEST_REAL_SIMILAR(p.getRetentionTime(), 5.0)
  TEST_EQUAL(p.getRetentionTimeUnit() == TargetedExperimentHelper::RetentionTime::RTUnit::SECOND, true)
  TEST_EQUAL(p.getRetentionTimeType() == TargetedExperimentHelper::RetentionTime::RTType::PREDICTED, true)

  TEST_EQUAL(p.getPeptideGroupLabel(), "")
  p.setPeptideGroupLabel("test1");
  TEST_EQUAL(p.getPeptideGroupLabel(), "test1")

  TEST_EQUAL(p.hasCharge(), false)
  p.setChargeState(-1);
  TEST_EQUAL(p.getChargeState(), -1)
  p.setChargeState(2);
  TEST_EQUAL(p.getChargeState(), 2)
}
END_SECTION

START_SECTION((TargetedExperiment::Compound))
{
  TargetedExperimentHelper::Compound p;
    
  TEST_EQUAL(p.hasRetentionTime(), false)
  TEST_EQUAL(p.rts.size(), 0)

  // add a RT
  TargetedExperimentHelper::RetentionTime rt;
  rt.setRT(5.0);
  rt.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
  rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;
  p.rts.push_back(rt);

  // test the RT methods
  TEST_EQUAL(p.rts.size(), 1)
  TEST_EQUAL(p.rts[0] == rt, true)
  TEST_EQUAL(p.rts[0].retention_time_unit == TargetedExperimentHelper::RetentionTime::RTUnit::SECOND, true)
  TEST_EQUAL(p.rts[0].retention_time_type == TargetedExperimentHelper::RetentionTime::RTType::PREDICTED, true)
  TEST_REAL_SIMILAR(p.rts[0].getRT(), 5.0)

  // test the Compound methods
  TEST_EQUAL(p.hasRetentionTime(), true)
  TEST_REAL_SIMILAR(p.getRetentionTime(), 5.0)
  TEST_EQUAL(p.getRetentionTimeUnit() == TargetedExperimentHelper::RetentionTime::RTUnit::SECOND, true)
  TEST_EQUAL(p.getRetentionTimeType() == TargetedExperimentHelper::RetentionTime::RTType::PREDICTED, true)

  // TEST_EQUAL(p.getPeptideGroupLabel(), "")
  // p.setPeptideGroupLabel("test1");
  // TEST_EQUAL(p.getPeptideGroupLabel(), "test1")

  TEST_EQUAL(p.hasCharge(), false)
  p.setChargeState(-1);
  TEST_EQUAL(p.getChargeState(), -1)
  p.setChargeState(2);
  TEST_EQUAL(p.getChargeState(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



