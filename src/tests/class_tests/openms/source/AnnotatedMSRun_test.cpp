// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Voigt $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/METADATA/AnnotatedMSRun.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <unordered_set>
#include <unordered_map>

START_TEST(AnnotatedMSRun, "$Id$")

using namespace OpenMS;

// Default constructor
AnnotatedMSRun* ptr = nullptr;
AnnotatedMSRun* nullPointer = nullptr;

START_SECTION((AnnotatedMSRun()))
  ptr = new AnnotatedMSRun();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~AnnotatedMSRun()))
  delete ptr;
END_SECTION

START_SECTION((explicit AnnotatedMSRun(MSExperiment&& experiment)))
  MSExperiment exp;
  MSSpectrum spec;
  spec.setRT(42.0);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);
  
  AnnotatedMSRun annotated_data(std::move(exp));
  TEST_EQUAL(annotated_data.getMSExperiment().size(), 1)
  TEST_REAL_SIMILAR(annotated_data.getMSExperiment()[0].getRT(), 42.0)
END_SECTION

START_SECTION((ProteinIdentification& getProteinIdentifications()))
  AnnotatedMSRun annotated_data;
  
  auto& prot_id = annotated_data.getProteinIdentifications();
  prot_id.resize(1);
  prot_id[0].setIdentifier("Test");
  TEST_EQUAL(annotated_data.getProteinIdentifications()[0].getIdentifier(), "Test")
END_SECTION

START_SECTION((const ProteinIdentification& getProteinIdentifications() const))
  AnnotatedMSRun annotated_data;
  auto& prot_id = annotated_data.getProteinIdentifications();
  prot_id.resize(1);
  prot_id[0].setIdentifier("Test");
  
  const AnnotatedMSRun& const_data = annotated_data;
  TEST_EQUAL(const_data.getProteinIdentifications()[0].getIdentifier(), "Test")
END_SECTION

START_SECTION((PeptideIdentification& getPeptideIdentification(size_t index)))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec1, spec2;
  exp.addSpectrum(spec1);
  exp.addSpectrum(spec2);
  annotated_data.getMSExperiment() = std::move(exp);
  
  // Resize peptide identifications to match spectra
  annotated_data.getPeptideIdentifications().resize(2);
  
  // Add data to the first peptide identification
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  annotated_data.getPeptideIdentifications()[0].insertHit(hit);
  
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPTIDE")
END_SECTION

START_SECTION((const PeptideIdentification& getPeptideIdentification(size_t index) const))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec1, spec2;
  exp.addSpectrum(spec1);
  exp.addSpectrum(spec2);
  annotated_data.getMSExperiment() = std::move(exp);
  
  // Resize peptide identifications to match spectra
  annotated_data.getPeptideIdentifications().resize(2);
  
  // Add data to the first peptide identification
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  annotated_data.getPeptideIdentifications()[0].insertHit(hit);
  
  const AnnotatedMSRun& const_data = annotated_data;
  TEST_EQUAL(const_data.getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(const_data.getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPTIDE")
END_SECTION


START_SECTION((PeptideIdentificationList& getPeptideIdentifications()))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec1, spec2;
  exp.addSpectrum(spec1);
  exp.addSpectrum(spec2);
  annotated_data.getMSExperiment() = std::move(exp);
  
  // Resize peptide identifications to match spectra
  annotated_data.getPeptideIdentifications().resize(2);
  
  // Add data to the peptide identifications
  PeptideHit hit1, hit2;
  hit1.setSequence(AASequence("PEPTIDER"));
  hit2.setSequence(AASequence("PEPTIDAR"));
  annotated_data.getPeptideIdentifications()[0].insertHit(hit1);
  annotated_data.getPeptideIdentifications()[1].insertHit(hit2);
  
  TEST_EQUAL(annotated_data.getPeptideIdentifications().size(), 2)
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[1].getHits().size(), 1)
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPTIDER")
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[1].getHits()[0].getSequence().toString(), "PEPTIDAR")
END_SECTION

START_SECTION((const PeptideIdentificationList& getPeptideIdentifications() const))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec1, spec2;
  exp.addSpectrum(spec1);
  exp.addSpectrum(spec2);
  annotated_data.getMSExperiment() = std::move(exp);
  
  // Resize peptide identifications to match spectra
  annotated_data.getPeptideIdentifications().resize(2);
  
  // Add data to the peptide identifications
  PeptideHit hit1, hit2;
  hit1.setSequence(AASequence::fromString("PEPTIDER"));
  hit2.setSequence(AASequence::fromString("PEPTIDAR"));
  annotated_data.getPeptideIdentifications()[0].insertHit(hit1);
  annotated_data.getPeptideIdentifications()[1].insertHit(hit2);
  
  const AnnotatedMSRun& const_data = annotated_data;
  TEST_EQUAL(const_data.getPeptideIdentifications().size(), 2)
  TEST_EQUAL(const_data.getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(const_data.getPeptideIdentifications()[1].getHits().size(), 1)
  TEST_EQUAL(const_data.getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPTIDER")
  TEST_EQUAL(const_data.getPeptideIdentifications()[1].getHits()[0].getSequence().toString(), "PEPTIDAR")
END_SECTION


START_SECTION((void setPeptideIdentification(PeptideIdentification&& id, size_t index)))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec1, spec2;
  exp.addSpectrum(spec1);
  exp.addSpectrum(spec2);
  annotated_data.getMSExperiment() = std::move(exp);
  
  // Resize peptide identifications to match spectra
  annotated_data.getPeptideIdentifications().resize(2);
  
  // Create a peptide identification
  PeptideIdentification pep_id;
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  pep_id.insertHit(hit);
  
  // Set the peptide identification
  annotated_data.getPeptideIdentifications()[0] = pep_id;
  
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPTIDE")
END_SECTION


START_SECTION((void setPeptideIdentifications(PeptideIdentificationList&& ids)))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec1, spec2;
  exp.addSpectrum(spec1);
  exp.addSpectrum(spec2);
  annotated_data.getMSExperiment() = std::move(exp);
  
  // Create a vector of peptide identifications
  PeptideIdentificationList pep_ids;
  PeptideIdentification pep_id1, pep_id2;
  PeptideHit hit1, hit2;
  hit1.setSequence(AASequence::fromString("PEPTIDER"));
  hit2.setSequence(AASequence::fromString("PEPTIDAR"));
  pep_id1.insertHit(hit1);
  pep_id2.insertHit(hit2);
  pep_ids.push_back(pep_id1);
  pep_ids.push_back(pep_id2);
  
  // Set all peptide identifications
  annotated_data.setPeptideIdentifications(std::move(pep_ids));
  
  TEST_EQUAL(annotated_data.getPeptideIdentifications().size(), 2)
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[1].getHits().size(), 1)
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[0].getHits()[0].getSequence().toString(), "PEPTIDER")
  TEST_EQUAL(annotated_data.getPeptideIdentifications()[1].getHits()[0].getSequence().toString(), "PEPTIDAR")
END_SECTION


START_SECTION((void clearAllPeptideIdentifications()))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec1, spec2;
  exp.addSpectrum(spec1);
  exp.addSpectrum(spec2);
  annotated_data.getMSExperiment() = std::move(exp);
  
  // Resize peptide identifications to match spectra
  annotated_data.getPeptideIdentifications().resize(2);
  
  // Add data to the peptide identifications
  PeptideHit hit1, hit2;
  hit1.setSequence(AASequence::fromString("PEPTIDER"));
  hit2.setSequence(AASequence::fromString("PEPTIDAR"));
  annotated_data.getPeptideIdentifications()[0].insertHit(hit1);
  annotated_data.getPeptideIdentifications()[1].insertHit(hit2);
  
  // Clear all peptide identifications
  annotated_data.getPeptideIdentifications().clear();
  
  TEST_EQUAL(annotated_data.getPeptideIdentifications().size(), 0)
END_SECTION

START_SECTION((MSExperiment& getMSExperiment()))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec;
  spec.setRT(42.0);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);
  
  annotated_data.getMSExperiment() = std::move(exp);
  TEST_EQUAL(annotated_data.getMSExperiment().size(), 1)
  TEST_REAL_SIMILAR(annotated_data.getMSExperiment()[0].getRT(), 42.0)
END_SECTION

START_SECTION((const MSExperiment& getMSExperiment() const))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec;
  spec.setRT(42.0);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);
  
  annotated_data.getMSExperiment() = std::move(exp);
  
  const AnnotatedMSRun& const_data = annotated_data;
  TEST_EQUAL(const_data.getMSExperiment().size(), 1)
  TEST_REAL_SIMILAR(const_data.getMSExperiment()[0].getRT(), 42.0)
END_SECTION

START_SECTION((Iterator functionality))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec1, spec2;
  spec1.setRT(10.0);
  spec2.setRT(20.0);
  exp.addSpectrum(spec1);
  exp.addSpectrum(spec2);
  annotated_data.getMSExperiment() = std::move(exp);
  
  // Resize peptide identifications to match spectra
  annotated_data.getPeptideIdentifications().resize(2);
  
  // Add data to the peptide identifications
  PeptideHit hit1, hit2;
  hit1.setSequence(AASequence::fromString("PEPTIDER"));
  hit2.setSequence(AASequence::fromString("PEPTIDAR"));
  annotated_data.getPeptideIdentifications()[0].insertHit(hit1);
  annotated_data.getPeptideIdentifications()[1].insertHit(hit2);
  
  // Test iterator functionality
  size_t count = 0;
  for (auto [spectrum, peptide_id] : annotated_data)
  {
    if (count == 0)
    {
      TEST_REAL_SIMILAR(spectrum.getRT(), 10.0)
      TEST_EQUAL(peptide_id.getHits().size(), 1)
      TEST_EQUAL(peptide_id.getHits()[0].getSequence().toString(), "PEPTIDER")
    }
    else if (count == 1)
    {
      TEST_REAL_SIMILAR(spectrum.getRT(), 20.0)
      TEST_EQUAL(peptide_id.getHits().size(), 1)
      TEST_EQUAL(peptide_id.getHits()[0].getSequence().toString(), "PEPTIDAR")
    }
    count++;
  }
  TEST_EQUAL(count, 2)
END_SECTION

START_SECTION((Operator[] functionality))
  AnnotatedMSRun annotated_data;
  MSExperiment exp;
  MSSpectrum spec1, spec2;
  spec1.setRT(10.0);
  spec2.setRT(20.0);
  exp.addSpectrum(spec1);
  exp.addSpectrum(spec2);
  annotated_data.getMSExperiment() = std::move(exp);
  
  // Resize peptide identifications to match spectra
  annotated_data.getPeptideIdentifications().resize(2);
  
  // Add data to the peptide identifications
  PeptideHit hit1, hit2;
  hit1.setSequence(AASequence::fromString("PEPTIDER"));
  hit2.setSequence(AASequence::fromString("PEPTIDAR"));
  annotated_data.getPeptideIdentifications()[0].insertHit(hit1);
  annotated_data.getPeptideIdentifications()[1].insertHit(hit2);
  
  // Test operator[] functionality
  auto [spectrum, peptide_id] = annotated_data[0];
  TEST_REAL_SIMILAR(spectrum.getRT(), 10.0)
  TEST_EQUAL(peptide_id.getHits().size(), 1)
  TEST_EQUAL(peptide_id.getHits()[0].getSequence().toString(), "PEPTIDER")
  
  auto [spectrum2, peptide_id2] = annotated_data[1];
  TEST_REAL_SIMILAR(spectrum2.getRT(), 20.0)
  TEST_EQUAL(peptide_id2.getHits().size(), 1)
  TEST_EQUAL(peptide_id2.getHits()[0].getSequence().toString(), "PEPTIDAR")
  
  // Test const operator[] functionality
  const AnnotatedMSRun& const_data = annotated_data;
  auto [const_spectrum, const_peptide_id] = const_data[0];
  TEST_REAL_SIMILAR(const_spectrum.getRT(), 10.0)
  TEST_EQUAL(const_peptide_id.getHits().size(), 1)
  TEST_EQUAL(const_peptide_id.getHits()[0].getSequence().toString(), "PEPTIDER")
END_SECTION

START_SECTION((bool operator==(const AnnotatedMSRun& rhs) const))
  // Create two identical AnnotatedMSRun objects
  AnnotatedMSRun run1, run2;

  // Set up MSExperiment
  MSExperiment exp;
  MSSpectrum spec;
  spec.setRT(42.0);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  run1.setMSExperiment(exp);
  run2.setMSExperiment(exp);

  // Set up peptide identifications
  PeptideIdentificationList pep_ids;
  PeptideIdentification pep_id;
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  pep_id.insertHit(hit);
  pep_ids.push_back(pep_id);

  run1.setPeptideIdentifications(pep_ids);
  run2.setPeptideIdentifications(pep_ids);

  // Set up protein identifications
  std::vector<ProteinIdentification> prot_ids;
  ProteinIdentification prot_id;
  prot_id.setIdentifier("TestProtein");
  prot_ids.push_back(prot_id);

  run1.setProteinIdentifications(prot_ids);
  run2.setProteinIdentifications(prot_ids);

  // Test equality
  TEST_EQUAL(run1 == run2, true)
  TEST_EQUAL(run1 != run2, false)

  // Test inequality - modify run2
  run2.getMSExperiment()[0].setRT(100.0);
  TEST_EQUAL(run1 == run2, false)
  TEST_EQUAL(run1 != run2, true)
END_SECTION

START_SECTION(([EXTRA] std::hash<AnnotatedMSRun>))
  // Test 1: Equal objects have equal hashes
  AnnotatedMSRun run1, run2;

  // Set up MSExperiment
  MSExperiment exp;
  MSSpectrum spec;
  spec.setRT(42.0);
  spec.setMSLevel(2);
  exp.addSpectrum(spec);

  run1.setMSExperiment(exp);
  run2.setMSExperiment(exp);

  // Set up peptide identifications
  PeptideIdentificationList pep_ids;
  PeptideIdentification pep_id;
  pep_id.setIdentifier("test_id");
  pep_id.setScoreType("test_score");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDE"));
  pep_id.insertHit(hit);
  pep_ids.push_back(pep_id);

  run1.setPeptideIdentifications(pep_ids);
  run2.setPeptideIdentifications(pep_ids);

  // Set up protein identifications
  std::vector<ProteinIdentification> prot_ids;
  ProteinIdentification prot_id;
  prot_id.setIdentifier("TestProtein");
  prot_id.setSearchEngine("TestEngine");
  prot_ids.push_back(prot_id);

  run1.setProteinIdentifications(prot_ids);
  run2.setProteinIdentifications(prot_ids);

  TEST_EQUAL(run1 == run2, true)
  TEST_EQUAL(std::hash<AnnotatedMSRun>{}(run1), std::hash<AnnotatedMSRun>{}(run2))

  // Test 2: Different objects have different hashes (likely, not guaranteed)
  AnnotatedMSRun run3;
  MSExperiment exp3;
  MSSpectrum spec3;
  spec3.setRT(100.0);  // Different RT
  spec3.setMSLevel(1); // Different MS level
  exp3.addSpectrum(spec3);
  run3.setMSExperiment(exp3);

  TEST_EQUAL(run1 == run3, false)
  TEST_NOT_EQUAL(std::hash<AnnotatedMSRun>{}(run1), std::hash<AnnotatedMSRun>{}(run3))

  // Test 3: Test use in unordered_set
  std::unordered_set<AnnotatedMSRun> run_set;
  run_set.insert(run1);
  run_set.insert(run2);  // Should not increase size (equal to run1)
  run_set.insert(run3);

  TEST_EQUAL(run_set.size(), 2)
  TEST_EQUAL(run_set.count(run1), 1)
  TEST_EQUAL(run_set.count(run3), 1)

  // Test 4: Test use in unordered_map
  std::unordered_map<AnnotatedMSRun, std::string> run_map;
  run_map[run1] = "first";
  run_map[run3] = "third";

  TEST_EQUAL(run_map.size(), 2)
  TEST_EQUAL(run_map[run1], "first")
  TEST_EQUAL(run_map[run2], "first")  // run2 == run1, so should get same value
  TEST_EQUAL(run_map[run3], "third")

  // Test 5: Empty AnnotatedMSRun has consistent hash
  AnnotatedMSRun empty1, empty2;
  TEST_EQUAL(empty1 == empty2, true)
  TEST_EQUAL(std::hash<AnnotatedMSRun>{}(empty1), std::hash<AnnotatedMSRun>{}(empty2))
END_SECTION

END_TEST