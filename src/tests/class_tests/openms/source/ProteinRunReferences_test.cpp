// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/ProteinRunReferences.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(ProteinRunReferences, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

auto run = [](const std::string& identifier)
{
  ProteinIdentification r;
  r.setIdentifier(identifier);
  return r;
};

auto peptide = [](const std::string& identifier)
{
  PeptideIdentification p;
  p.setIdentifier(identifier);
  return p;
};

const std::string missing_message =
  "Peptide identification has no matching protein run: 'other'. Every peptide identification needs the protein "
  "identification run (search run) with its identifier, which may have no protein hits.";

START_SECTION(static std::string missingRunMessage(const std::string& identifier))
{
  TEST_STRING_EQUAL(ProteinRunReferences::missingRunMessage("other"), missing_message)
}
END_SECTION

START_SECTION(static void check(const std::vector<ProteinIdentification>& runs, const std::string& identifier))
{
  const std::vector<ProteinIdentification> runs = {run("search"), run("")};
  ProteinRunReferences::check(runs, "search");
  ProteinRunReferences::check(runs, ""); // an empty identifier is an identifier like any other
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ProteinRunReferences::check(runs, "other"), missing_message)
  TEST_EXCEPTION(Exception::InvalidParameter, ProteinRunReferences::check(std::vector<ProteinIdentification>(), ""))
}
END_SECTION

START_SECTION(static void check(const std::vector<ProteinIdentification>& runs, const PeptideIdentificationList& peptides))
{
  // a run without protein hits (e.g. peptidomics or de novo) is a run
  const std::vector<ProteinIdentification> runs = {run("search")};
  PeptideIdentificationList peptides;
  ProteinRunReferences::check(std::vector<ProteinIdentification>(), peptides); // nothing to check
  peptides.push_back(peptide("search"));
  ProteinRunReferences::check(runs, peptides);
  peptides.push_back(peptide("other"));
  peptides.push_back(peptide("third"));
  // the first peptide identification without its run is reported
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ProteinRunReferences::check(runs, peptides), missing_message)
  TEST_EXCEPTION(Exception::InvalidParameter, ProteinRunReferences::check(std::vector<ProteinIdentification>(), peptides))
}
END_SECTION

START_SECTION(static void check(const std::vector<ProteinIdentification>& runs, const std::vector<const PeptideIdentification*>& peptides))
{
  const std::vector<ProteinIdentification> runs = {run("search")};
  const PeptideIdentification good = peptide("search"), bad = peptide("other");
  ProteinRunReferences::check(runs, std::vector<const PeptideIdentification*>{&good});
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ProteinRunReferences::check(runs, std::vector<const PeptideIdentification*>{&good, &bad}), missing_message)
}
END_SECTION

START_SECTION(static void checkProteinAccessions(const std::vector<ProteinIdentification>& runs, const std::vector<const PeptideIdentification*>& peptides))
{
  auto with_protein = [](const std::string& identifier, const std::string& accession)
  {
    PeptideIdentification p;
    p.setIdentifier(identifier);
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDE"));
    hit.setCharge(2);
    hit.addPeptideEvidence(PeptideEvidence(accession, 0, 6, '-', '-'));
    hit.addPeptideEvidence(PeptideEvidence("", 0, 6, '-', '-')); // empty accessions are not written
    p.insertHit(hit);
    return p;
  };
  ProteinIdentification search = run("search"), other = run("other");
  search.insertHit(ProteinHit(0.0, 1, "PROT_A", ""));
  other.insertHit(ProteinHit(0.0, 1, "PROT_B", ""));
  const std::vector<ProteinIdentification> runs = {search, other};
  const auto first = with_protein("search", "PROT_A"), second = with_protein("other", "PROT_B");
  ProteinRunReferences::checkProteinAccessions(runs, {&first, &second});

  // a protein of another run
  const auto cross_run = with_protein("search", "PROT_B");
  TEST_EXCEPTION_WITH_MESSAGE(Exception::ElementNotFound, ProteinRunReferences::checkProteinAccessions(runs, {&first, &second, &cross_run}),
    "the element 'No accession PROT_B found in run 'search' for PSM PEPTIDE_2. Every protein of a peptide evidence needs to be a protein hit of the peptide identification's run.' could not be found")

  // a protein that is no protein hit of any run
  const auto unknown = with_protein("search", "PROT_C");
  TEST_EXCEPTION(Exception::ElementNotFound, ProteinRunReferences::checkProteinAccessions(runs, {&first, &unknown}))

  // a peptide identification without its run (one pass checks both)
  const auto no_run = with_protein("third", "PROT_A");
  TEST_EXCEPTION(Exception::InvalidParameter, ProteinRunReferences::checkProteinAccessions(runs, {&first, &second, &no_run}))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
