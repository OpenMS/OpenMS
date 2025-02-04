// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabet.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetTextParser.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

class IMSAlphabetParserImpl
    : public IMSAlphabetParser<>
{
private:
  ContainerType elements_;
public:
  ContainerType& getElements() override
  {
    return elements_;
  }

  void parse(std::istream& ) override
  {
    // ignore istream, just enter something into the map
    elements_.insert(std::make_pair("hydrogen", 1.0));
    elements_.insert(std::make_pair("oxygen", 16.0));
  }
};


START_TEST(IMSAlphabet, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IMSElement hydrogen("hydrogen", 1.0);
IMSElement oxygen("oxygen", 16.0);
IMSElement nitrogen("nitrogen", 14.0);

vector<IMSElement> elements;
elements.push_back(hydrogen);
elements.push_back(oxygen);
elements.push_back(nitrogen);

double carbon_mass = 12.0;
IMSElement carbon("carbon", carbon_mass);

IMSAlphabet* ptr = nullptr;
IMSAlphabet* null_ptr = nullptr;
START_SECTION(IMSAlphabet())
{
  ptr = new IMSAlphabet();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IMSAlphabet())
{
  delete ptr;
}
END_SECTION

START_SECTION((IMSAlphabet(const container &elements)))
{
  ptr = new IMSAlphabet(elements);
  TEST_NOT_EQUAL(ptr, null_ptr)

  TEST_EQUAL(ptr->size(), 3)
  TEST_EQUAL(ptr->getName(0), "hydrogen")

  delete ptr;
}
END_SECTION

// test instance for the following tests
IMSAlphabet alphabet(elements);

START_SECTION((IMSAlphabet(const IMSAlphabet &alphabet)))
{
  IMSAlphabet alphabet_copy(alphabet);

  TEST_EQUAL(alphabet_copy.size(), 3)
  TEST_EQUAL(alphabet_copy.getName(0), "hydrogen")
}
END_SECTION

START_SECTION((const element_type& getElement(const name_type &name) const ))
{
  TEST_EQUAL(alphabet.getElement("hydrogen"), hydrogen)
  TEST_EQUAL(alphabet.getElement("oxygen"), oxygen)
  TEST_EQUAL(alphabet.getElement("nitrogen"), nitrogen)
}
END_SECTION

START_SECTION((const element_type& getElement(size_type index) const ))
{
  TEST_EQUAL(alphabet.getElement(0), hydrogen)
  TEST_EQUAL(alphabet.getElement(1), oxygen)
  TEST_EQUAL(alphabet.getElement(2), nitrogen)

  TEST_EXCEPTION(Exception::InvalidValue, alphabet.getElement("nitrogen2"))
}
END_SECTION

START_SECTION((void setElement(const name_type &name, mass_type mass, bool forced=false)))
{
  IMSAlphabet alphabet_copy(alphabet);

  alphabet_copy.setElement("hydrogen", 2.0, false);
  TEST_EQUAL(alphabet_copy.size(), 3)
  TEST_EQUAL(alphabet_copy.getMass("hydrogen"), 2.0)

  // this should not change the alphabet, since there is no
  // element named carbon
  alphabet_copy.setElement("carbon", carbon_mass, false);
  TEST_EQUAL(alphabet_copy.size(), 3)

  alphabet_copy.setElement("carbon", carbon_mass, true);
  TEST_EQUAL(alphabet_copy.size(), 4)
  TEST_EQUAL(alphabet_copy.getMass("carbon"), 12)
}
END_SECTION

START_SECTION((const name_type& getName(size_type index) const ))
{
  TEST_STRING_EQUAL(alphabet.getName(0),"hydrogen")
  TEST_STRING_EQUAL(alphabet.getName(1),"oxygen")
  TEST_STRING_EQUAL(alphabet.getName(2),"nitrogen")
}
END_SECTION

START_SECTION((mass_type getMass(const name_type &name) const ))
{
  TEST_EQUAL(alphabet.getMass("hydrogen"), hydrogen.getMass())
  TEST_EQUAL(alphabet.getMass("oxygen"), oxygen.getMass())
  TEST_EQUAL(alphabet.getMass("nitrogen"), nitrogen.getMass())

  TEST_EXCEPTION(Exception::InvalidValue, alphabet.getMass("nitrogen2"))
}
END_SECTION

START_SECTION((mass_type getMass(size_type index) const ))
{
  TEST_EQUAL(alphabet.getMass(0), hydrogen.getMass())
  TEST_EQUAL(alphabet.getMass(1), oxygen.getMass())
  TEST_EQUAL(alphabet.getMass(2), nitrogen.getMass())
}
END_SECTION

START_SECTION((masses_type getMasses(size_type isotope_index=0) const ))
{
  IMSAlphabet::masses_type masses = alphabet.getMasses();
  TEST_EQUAL(masses.size(), 3)
  ABORT_IF(masses.size() != 3)

  TEST_EQUAL(masses[0], hydrogen.getMass())
  TEST_EQUAL(masses[1], oxygen.getMass())
  TEST_EQUAL(masses[2], nitrogen.getMass())
}
END_SECTION

START_SECTION((masses_type getAverageMasses() const ))
{
  IMSAlphabet::masses_type average_masses = alphabet.getAverageMasses();

  TEST_EQUAL(average_masses.size(), 3)
  ABORT_IF(average_masses.size() != 3)

  TEST_EQUAL(average_masses[0], hydrogen.getAverageMass())
  TEST_EQUAL(average_masses[1], oxygen.getAverageMass())
  TEST_EQUAL(average_masses[2], nitrogen.getAverageMass())
}
END_SECTION

START_SECTION((bool hasName(const name_type &name) const ))
{
  TEST_EQUAL(alphabet.hasName("nitrogen"), true)
  TEST_EQUAL(alphabet.hasName("oxygen"), true)
  TEST_EQUAL(alphabet.hasName("oxygen2"), false)
}
END_SECTION

START_SECTION((void push_back(const name_type &name, mass_type value)))
{
  IMSAlphabet alphabet_copy(alphabet);
  alphabet_copy.push_back("carbon", carbon_mass);
  TEST_EQUAL(alphabet_copy.size(), 4)
  ABORT_IF(alphabet_copy.size() != 4)

  TEST_EQUAL(alphabet_copy.getElement(3).getName(), "carbon")
  TEST_EQUAL(alphabet_copy.hasName("carbon"),true)
  TEST_EQUAL(alphabet_copy.getMass(3), carbon_mass)
}
END_SECTION

START_SECTION((void push_back(const element_type &element)))
{
  IMSAlphabet alphabet_copy(alphabet);
  alphabet_copy.push_back(carbon);
  TEST_EQUAL(alphabet_copy.size(), 4)
  ABORT_IF(alphabet_copy.size() != 4)

  TEST_EQUAL(alphabet_copy.getElement(3).getName(), "carbon")
  TEST_EQUAL(alphabet_copy.hasName("carbon"),true)
  TEST_EQUAL(alphabet_copy.getMass(3), carbon.getMass())
}
END_SECTION

START_SECTION((void clear()))
{
  IMSAlphabet alphabet_copy(alphabet);
  alphabet_copy.clear();
  TEST_EQUAL(alphabet_copy.size(), 0)
  TEST_EQUAL(alphabet_copy.hasName("oxygen"), false)
}
END_SECTION

START_SECTION((virtual void sortByNames()))
{
  IMSAlphabet alphabet_copy(alphabet);
  alphabet_copy.sortByNames();
  TEST_EQUAL(alphabet_copy.size(), 3)
  TEST_EQUAL(alphabet_copy.getElement(0), hydrogen)
  TEST_EQUAL(alphabet_copy.getElement(1), nitrogen)
  TEST_EQUAL(alphabet_copy.getElement(2), oxygen)
}
END_SECTION

START_SECTION((virtual void sortByValues()))
{
  IMSAlphabet alphabet_copy(alphabet);
  alphabet_copy.push_back(carbon);
  alphabet_copy.sortByValues();
  TEST_EQUAL(alphabet_copy.size(), 4)
  TEST_EQUAL(alphabet_copy.getElement(0), hydrogen)
  TEST_EQUAL(alphabet_copy.getElement(1), carbon)
  TEST_EQUAL(alphabet_copy.getElement(2), nitrogen)
  TEST_EQUAL(alphabet_copy.getElement(3), oxygen)
}
END_SECTION

START_SECTION((virtual void load(const std::string &fname)))
{
  NOT_TESTABLE // will be tested in virtual void load(const std::string &fname, IMSAlphabetParser<> *parser)
}
END_SECTION

START_SECTION((virtual void load(const std::string &fname, IMSAlphabetParser<> *parser)))
{
  IMSAlphabet load_copy(alphabet);
  TEST_EQUAL(load_copy.size(), 3)

  TEST_EXCEPTION(Exception::IOException, load_copy.load(""))

  // this should not clear the alphabet
  TEST_EQUAL(load_copy.size(), 3)

  // create file with minimal content
  String filename;
  NEW_TMP_FILE(filename)

  // just create the file
  ofstream of;
  of.open(filename.c_str());
  of << "# a comment which should be ignored" << std::endl;
  of << "hydrogen\t1.0" << std::endl;
  of << "oxygen\t16.0" << std::endl;
  of << "nitrogen\t14.0" << std::endl;
  of.close();

  IMSAlphabetTextParser parser;
  load_copy.load(filename, parser);

  TEST_EQUAL(load_copy.size(), 3)
  TEST_EQUAL(load_copy.hasName("hydrogen"), true)
  TEST_REAL_SIMILAR(load_copy.getMass("hydrogen"), 1.0)
  TEST_EQUAL(load_copy.hasName("nitrogen"), true)
  TEST_REAL_SIMILAR(load_copy.getMass("nitrogen"), 14.0)
  TEST_EQUAL(load_copy.hasName("oxygen"), true)
  TEST_REAL_SIMILAR(load_copy.getMass("oxygen"), 16.0)
}
END_SECTION

START_SECTION((size_type size() const ))
{
  TEST_EQUAL(alphabet.size(),3)

  IMSAlphabet alphabet_copy(alphabet);
  TEST_EQUAL(alphabet_copy.size(),3)

  alphabet_copy.push_back(carbon);
  TEST_EQUAL(alphabet_copy.size(), 4)
}
END_SECTION

START_SECTION((bool erase(const name_type &name)))
{
  IMSAlphabet alphabet_copy(alphabet);
  TEST_EQUAL(alphabet_copy.size(),3)

  TEST_EQUAL(alphabet_copy.erase("hydrogen"), true)
  TEST_EQUAL(alphabet_copy.size(),2)
  TEST_EQUAL(alphabet_copy.erase("hydrogen"), false)

  TEST_EQUAL(alphabet_copy.erase("oxygen"), true)
  TEST_EQUAL(alphabet_copy.size(),1)

  TEST_EQUAL(alphabet_copy.erase("nitrogen"), true)
  TEST_EQUAL(alphabet_copy.size(),0)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



