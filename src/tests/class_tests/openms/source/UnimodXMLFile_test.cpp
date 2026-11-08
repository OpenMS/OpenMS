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
#include <OpenMS/SYSTEM/File.h>

#include <cmath>
#include <fstream>
#include <iterator>
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

START_SECTION([EXTRA] neutral loss formulas and neutral loss masses agree for every entry)
{
  // A <NeutralLoss> carries its own mono_mass/avge_mass. Those attributes were never read - the
  // <delta> masses were pushed instead, and since <specificity> (which encloses <NeutralLoss>)
  // always precedes <delta>, every neutral loss mass in the DB was 0. The mass vectors were also
  // written to the modification once per <specificity> rather than per site, so the last site's
  // vector was copied over all of them and could not even match the per-site formula vector in
  // length. See https://github.com/OpenMS/OpenMS/issues/10030
  UnimodXMLFile file;
  vector<ResidueModification*> modifications;
  file.load("CHEMISTRY/unimod.xml", modifications);

  const double mono_tolerance = 1e-3;  // as above: OpenMS vs. UniMod element masses
  const double avg_tolerance = 0.05;   // average atomic weights differ more; worst observed is ~0.006 Da
  Size with_neutral_loss = 0;
  Size length_mismatches = 0;
  Size mass_mismatches = 0;
  for (const ResidueModification* mod : modifications)
  {
    if (!mod->hasNeutralLoss()) continue;
    ++with_neutral_loss;

    const vector<EmpiricalFormula>& formulas = mod->getNeutralLossDiffFormulas();
    const vector<double> mono_masses = mod->getNeutralLossMonoMasses();
    const vector<double> avg_masses = mod->getNeutralLossAverageMasses();
    if (formulas.size() != mono_masses.size() || formulas.size() != avg_masses.size())
    {
      ++length_mismatches;
      STATUS(mod->getId() << " (" << mod->getOrigin() << "): " << formulas.size() << " neutral loss formula(s) but "
             << mono_masses.size() << " mono and " << avg_masses.size() << " average mass(es)");
      continue;
    }
    for (Size i = 0; i < formulas.size(); ++i)
    {
      if (std::fabs(formulas[i].getMonoWeight() - mono_masses[i]) > mono_tolerance ||
          std::fabs(formulas[i].getAverageWeight() - avg_masses[i]) > avg_tolerance)
      {
        ++mass_mismatches;
        STATUS(mod->getId() << " (" << mod->getOrigin() << "): neutral loss " << formulas[i].toString()
               << " weighs " << formulas[i].getMonoWeight() << " / " << formulas[i].getAverageWeight()
               << " but its masses are " << mono_masses[i] << " / " << avg_masses[i]);
      }
    }
  }
  TEST_EQUAL(with_neutral_loss > 0, true)
  TEST_EQUAL(length_mismatches, 0)
  TEST_EQUAL(mass_mismatches, 0)

  // cleanup
  for (Size k = 0; k < modifications.size(); k++)
  {
    delete modifications[k];
  }
}
END_SECTION

START_SECTION([EXTRA] a file parsed in slices gives the modifications of the file parsed as a whole)
{
  // unimod.xml is large enough to be parsed in slices; a comment in its <umod:modifications> section makes the
  // same content parse as a whole (see UnimodXMLFile::loadInSlices_())
  std::string text;
  {
    std::ifstream in(File::find("CHEMISTRY/unimod.xml"), std::ios::binary);
    text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  }
  const std::string open_tag = "<umod:modifications>";
  const size_t open = text.find(open_tag);
  TEST_EQUAL(open != std::string::npos, true)
  auto write = [](const std::string& content)
  {
    std::string file;
    NEW_TMP_FILE(file)
    std::ofstream(file, std::ios::binary) << content;
    return file;
  };
  // loads 'content' (in slices if the loader can) and the same with a comment in the section (as a whole) and
  // compares; returns the number of modifications
  auto compare = [&](const std::string& content)
  {
    std::string whole = content;
    whole.insert(whole.rfind(open_tag) + open_tag.size(), "<!-- parsed as a whole -->"); // the real section
    vector<ResidueModification*> sliced, as_whole;
    UnimodXMLFile().load(write(content), sliced);
    UnimodXMLFile().load(write(whole), as_whole);
    TEST_EQUAL(sliced.size(), as_whole.size())
    Size differences = 0;
    for (Size i = 0; i < std::min(sliced.size(), as_whole.size()); ++i)
    {
      if (!(*sliced[i] == *as_whole[i]) || sliced[i]->getFullId() != as_whole[i]->getFullId() ||
          sliced[i]->getNeutralLossMonoMasses() != as_whole[i]->getNeutralLossMonoMasses()) ++differences;
    }
    TEST_EQUAL(differences, 0)
    const Size n = sliced.size();
    for (auto* m : sliced) delete m;
    for (auto* m : as_whole) delete m;
    return n;
  };
  const Size n_unimod = compare(text);
  TEST_EQUAL(n_unimod > 1000, true)

  // a comment before the root element that looks like a large section: it is not read as modifications
  {
    std::string fake = "<!-- " + open_tag + std::string(600 * 1024, ' ');
    for (int i = 0; i < 8; ++i) fake += "<umod:mod title=\"x\"></umod:mod>" + std::string(70 * 1024, ' ');
    fake += "</umod:modifications> -->\n";
    std::string content = text;
    const size_t root = content.find("<umod:unimod");
    content.insert(root, fake);
    TEST_EQUAL(compare(content), n_unimod)
  }
  // a neutral loss after a specificity (its masses carry over into the next modification when read as a whole)
  {
    std::string content = text;
    const size_t p = content.find("</umod:specificity>", content.size() / 2) + 19;
    content.insert(p, "<umod:NeutralLoss mono_mass=\"1.5\" avge_mass=\"1.5\" flag=\"false\"><umod:element symbol=\"H\" number=\"1\"/></umod:NeutralLoss>");
    compare(content);
  }
  // an unknown position in single quotes, with blanks around '='
  {
    std::string content = text;
    const size_t p = content.find("position=\"Anywhere\"", content.size() / 2);
    content.replace(p, 19, "position = 'Somewhere'");
    compare(content);
  }

  // an error in one slice: the whole file is parsed again and reports it
  std::string broken = text;
  const size_t middle = broken.find("</umod:mod>", broken.size() / 2);
  broken.replace(middle, 11, "</umod:mox>");
  vector<ResidueModification*> none;
  TEST_EXCEPTION(Exception::ParseError, UnimodXMLFile().load(write(broken), none))
  for (auto* m : none) delete m;

  // a nested, unclosed <umod:modifications> in the section (malformed): an error as with the whole file
  {
    std::string nested = text;
    const size_t p = nested.find("</umod:mod>", nested.size() / 3) + 11;
    nested.insert(p, open_tag);
    vector<ResidueModification*> mods;
    TEST_EXCEPTION(Exception::ParseError, UnimodXMLFile().load(write(nested), mods))
    for (auto* m : mods) delete m;
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
