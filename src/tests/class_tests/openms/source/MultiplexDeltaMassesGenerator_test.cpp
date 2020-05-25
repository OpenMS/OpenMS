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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MultiplexDeltaMassesGenerator.h>

using namespace OpenMS;

START_TEST(MultiplexDeltaMassesGenerator, "$Id$")

std::map<String, double> label_mass_shift;
label_mass_shift.insert(std::make_pair("Arg6", 6.0201290268));
label_mass_shift.insert(std::make_pair("Arg10", 10.008268600));
label_mass_shift.insert(std::make_pair("Lys4", 4.0251069836));
label_mass_shift.insert(std::make_pair("Lys8", 8.0141988132));

// triple SILAC
String labels = "[][Lys4,Arg6][Lys8,Arg10]";
int missed_cleavages = 1;

MultiplexDeltaMassesGenerator* nullPointer = nullptr;
MultiplexDeltaMassesGenerator* ptr;

START_SECTION(MultiplexDeltaMassesGenerator(String labels, int missed_cleavages, std::map<String,double> label_mass_shift))
    MultiplexDeltaMassesGenerator list(labels, missed_cleavages, label_mass_shift);
    TEST_EQUAL(list.getDeltaMassesList().size(), 5);
    ptr = new MultiplexDeltaMassesGenerator(labels, missed_cleavages, label_mass_shift);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexDeltaMassesGenerator list(labels, missed_cleavages, label_mass_shift);

START_SECTION(std::vector<MultiplexDeltaMasses> getDeltaMassesList())
  std::vector<MultiplexDeltaMasses> masses = list.getDeltaMassesList();
  TEST_EQUAL(masses.size(), 5);
  TEST_REAL_SIMILAR(masses[2].getDeltaMasses()[1].delta_mass, 8.0502139672);
  TEST_REAL_SIMILAR(masses[4].getDeltaMasses()[2].delta_mass, 20.0165372);
END_SECTION

START_SECTION(std::vector<std::vector<String> > MultiplexDeltaMassesGenerator::getSamplesLabelsList())
  std::vector<std::vector<String> > samples_labels = list.getSamplesLabelsList();
  TEST_EQUAL(samples_labels.size(), 3);
  TEST_EQUAL(samples_labels[1][0], "Lys4");
  TEST_EQUAL(samples_labels[2][1], "Arg10");
END_SECTION

START_SECTION(void generateKnockoutDeltaMasses())
  list.generateKnockoutDeltaMasses();
  std::vector<MultiplexDeltaMasses> masses_knockout = list.getDeltaMassesList();
  TEST_EQUAL(masses_knockout.size(), 21);
  TEST_REAL_SIMILAR(masses_knockout[6].getDeltaMasses()[1].delta_mass, 8.0141988132);
  TEST_REAL_SIMILAR(masses_knockout[19].getDeltaMasses()[1].delta_mass, 20.0165372);
  TEST_EQUAL(masses_knockout[20].getDeltaMasses().size(), 1);
END_SECTION

START_SECTION(String MultiplexDeltaMassesGenerator::getLabelShort(String label))
  TEST_EQUAL(list.getLabelShort("Label:13C(6)15N(2)"), "Lys8");
  TEST_EQUAL(list.getLabelShort("Dimethyl:2H(6)13C(2)"), "Dimethyl8");
END_SECTION

START_SECTION(String MultiplexDeltaMassesGenerator::getLabelLong(String label))
  TEST_EQUAL(list.getLabelLong("Lys8"), "Label:13C(6)15N(2)");
  TEST_EQUAL(list.getLabelLong("Dimethyl8"), "Dimethyl:2H(6)13C(2)");
END_SECTION

START_SECTION(MultiplexDeltaMasses::LabelSet MultiplexDeltaMassesGenerator::extractLabelSet(AASequence sequence))
  AASequence sequence = AASequence::fromString("LAPITSDPTEAAAVGAVEASFK(Label:13C(6)15N(2))");
  MultiplexDeltaMasses::LabelSet label_set = list.extractLabelSet(sequence);
  TEST_EQUAL(label_set.size(), 1);
  TEST_EQUAL(*(label_set.begin()), "Lys8");
END_SECTION

END_TEST
