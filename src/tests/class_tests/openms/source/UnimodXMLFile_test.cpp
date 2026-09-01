// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <cmath>
#include <vector>

///////////////////////////

START_TEST(UnimodXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

UnimodXMLFile xml_file;
UnimodXMLFile* ptr;
UnimodXMLFile* nullPointer = nullptr;

START_SECTION((UnimodXMLFile()))
	ptr = new UnimodXMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~UnimodXMLFile())
	delete ptr;
END_SECTION

ptr = new UnimodXMLFile();

START_SECTION(void load(const std::string& filename, vector<ResidueModification*>& modifications))
	vector<ResidueModification*> modifications;
	ptr->load("CHEMISTRY/unimod.xml", modifications);

	//cerr << "#modifications read: " << modifications.size() << endl;
	//for (vector<ResidueModification>::const_iterator it = modifications.begin(); it != modifications.end(); ++it)
	//{
	//	cerr << it->getTitle() << "\t" << it->getFullName() << "\t" << it->getAllowedPositionName() << "\t" << it->getSite() << "\t" << it->getClassification() << "\t" << it->getComposition() << "\t" << it->getMonoMass() << endl;
	//}

	TEST_EQUAL(modifications.size() > 1, true)

  // cleanup
  for (Size k = 0; k < modifications.size(); k++)
  {
    delete modifications[k];
  }
END_SECTION

delete ptr;

START_SECTION([EXTRA] diff formula and diff mono mass agree for every entry)
{
  // The composition of a modification is given by the <element> children of its <delta> block.
  // <Ignore> blocks (reagent fragments a search must not consider) carry <element> children as
  // well; folding those into the formula made ICAT-D (UNIMOD:12) and ICAT-D:2H(8) (UNIMOD:13)
  // roughly 2250 Da heavier than their delta mass. getDiffMonoMass() is read from the mono_mass
  // attribute and stayed correct, so only formula-derived results - AASequence::getFormula() and
  // every isotope pattern computed from it - were wrong, silently.
  // See https://github.com/OpenMS/OpenMS/issues/10030
  UnimodXMLFile file;
  vector<ResidueModification*> modifications;
  file.load("CHEMISTRY/unimod.xml", modifications);

  // OpenMS' element masses differ marginally from those UniMod used for its mono_mass attribute;
  // the largest deviation over all entries is ~3e-5 Da, so 1e-3 Da leaves ample headroom while
  // still catching a stray element.
  const double tolerance = 1e-3;
  Size mismatches = 0;
  for (const ResidueModification* mod : modifications)
  {
    const double from_formula = mod->getDiffFormula().getMonoWeight();
    if (std::fabs(from_formula - mod->getDiffMonoMass()) > tolerance)
    {
      ++mismatches;
      STATUS(mod->getUniModAccession() << " (" << mod->getId() << "): diff formula "
             << mod->getDiffFormula().toString() << " weighs " << from_formula
             << " but diff mono mass is " << mod->getDiffMonoMass());
    }
  }
  TEST_EQUAL(mismatches, 0)

  // the two entries the scoping bug hit, checked by name
  for (const ResidueModification* mod : modifications)
  {
    if (mod->getId() == "ICAT-D")
    {
      TEST_EQUAL(mod->getDiffFormula() == EmpiricalFormula("C20H34N4O5S"), true)
      TEST_REAL_SIMILAR(mod->getDiffMonoMass(), 442.224991)
    }
    else if (mod->getId() == "ICAT-D:2H(8)")
    {
      TEST_EQUAL(mod->getDiffFormula() == EmpiricalFormula("C20H26(2)H8N4O5S"), true)
      TEST_REAL_SIMILAR(mod->getDiffMonoMass(), 450.275205)
    }
  }

  // cleanup
  for (Size k = 0; k < modifications.size(); k++)
  {
    delete modifications[k];
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
