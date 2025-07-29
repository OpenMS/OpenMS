// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileTypes.h>

#include <boost/assign/std/vector.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TransitionTSVFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TransitionTSVFile* ptr = nullptr;
TransitionTSVFile* nullPointer = nullptr;

START_SECTION(TransitionTSVFile())
{
  ptr = new TransitionTSVFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~TransitionTSVFile())
{
  delete ptr;
}
END_SECTION

START_SECTION( void convertTargetedExperimentToTSV(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
{
  // see TOPP tool test
  NOT_TESTABLE
}
END_SECTION

START_SECTION( void convertTSVToTargetedExperiment(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
{
  // see TOPP tool test
  NOT_TESTABLE
}
END_SECTION

START_SECTION( void validateTargetedExperiment(OpenMS::TargetedExperiment & targeted_exp))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION( void createOligo_(const String & oligo_name, OpenMS::TargetedExperiment::Oligo & oligo))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION( void convertTSVToTargetedExperiment(const char * filename, FileTypes::Type filetype, OpenSwath::LightTargetedExperiment & targeted_exp))
{
  // Create a minimal nucleotide TSV content
  String tsv_content =
    "PrecursorMz\tProductMz\tLibraryIntensity\tNormalizedRetentionTime\tNuctideSequence\tFullNuctideName\tOligoName\tTransitionGroupId\tTransitionId\n"
    "500.25\t400.15\t1000.0\t25.5\tAUGCGAUU\tAUG[Cm]CGAUU\ttest_oligo\tnuc_group_1\ttransition_1\n";
  
  // Write to temporary file
  String temp_filename = File::getTemporaryFile();
  std::ofstream temp_file(temp_filename.c_str());
  temp_file << tsv_content;
  temp_file.close();
  
  // Test reading
  TransitionTSVFile tsv_file;
  TargetedExperiment exp;
  
  tsv_file.convertTSVToTargetedExperiment(temp_filename.c_str(), FileTypes::TSV, exp);
  
  // Verify the experiment contains nucleotide data
  TEST_EQUAL(exp.getTransitions().size(), 1)
  TEST_EQUAL(exp.getNuctides().size(), 1)
  TEST_EQUAL(exp.getOligos().size(), 1)
  
  // Test nucleotide properties
  const TargetedExperiment::Nuctide& nuctide = exp.getNuctides()[0];
  TEST_EQUAL(nuctide.id, "nuc_group_1")
  TEST_EQUAL(nuctide.getSequence(), "AUGCGAUU")
  
  // Test oligo properties
  const TargetedExperiment::Oligo& oligo = exp.getOligos()[0];
  TEST_EQUAL(oligo.id, "test_oligo")
  
  // Test transition properties
  const ReactionMonitoringTransition& transition = exp.getTransitions()[0];
  TEST_REAL_SIMILAR(transition.getPrecursorMZ(), 500.25)
  TEST_REAL_SIMILAR(transition.getProductMZ(), 400.15)
  TEST_EQUAL(transition.getTransGroupRef(), "nuc_group_1")
  
  // Clean up
  File::remove(temp_filename);
}
END_SECTION

START_SECTION( void convertTargetedExperimentToTSV(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
{
  // Create a TargetedExperiment with nucleotide data
  TargetedExperiment exp;
  
  // Create nucleotide
  TargetedExperiment::Nuctide nuctide;
  nuctide.id = "nuc_group_test";
  nuctide.setSequence("AUGCGAUU");
  nuctide.setMetaValue("full_nuctide_name", "AUG[Cm]CGAUU");
  
  // Add retention time
  TargetedExperimentHelper::RetentionTime rt;
  rt.setRT(30.5);
  rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::IRT;
  nuctide.rts.push_back(rt);
  
  exp.addNuctide(nuctide);
  
  // Create oligo
  TargetedExperiment::Oligo oligo;
  oligo.id = "test_oligo_write";
  exp.addOligo(oligo);
  
  // Create transition
  ReactionMonitoringTransition transition;
  transition.setNativeID("transition_test");
  transition.setPrecursorMZ(500.25);
  transition.setProductMZ(400.15);
  transition.setLibraryIntensity(1500.0);
  transition.setTransGroupRef("nuc_group_test", OpenSwath::TransType::NUCTIDE);
  
  exp.addTransition(transition);
  
  // Write to TSV
  String temp_filename = File::getTemporaryFile();
  TransitionTSVFile tsv_file;
  tsv_file.convertTargetedExperimentToTSV(temp_filename.c_str(), exp);
  
  // Read back and verify
  std::ifstream temp_file(temp_filename.c_str());
  String line;
  
  // Read header
  std::getline(temp_file, line);
  TEST_EQUAL(line.hasSubstring("NuctideSequence"), true)
  TEST_EQUAL(line.hasSubstring("OligoName"), true)
  
  // Read data line
  std::getline(temp_file, line);
  TEST_EQUAL(line.hasSubstring("AUGCGAUU"), true)
  TEST_EQUAL(line.hasSubstring("500.25"), true)
  // Check for ProductMZ value with floating point tolerance
  TEST_EQUAL(line.hasSubstring("400.14") || line.hasSubstring("400.15"), true)
  
  temp_file.close();
  
  // Clean up
  File::remove(temp_filename);

  // Create TSV with both peptide and nucleotide data
  String tsv_content =
    "PrecursorMz\tProductMz\tLibraryIntensity\tNormalizedRetentionTime\tPeptideSequence\tFullPeptideName\tProteinId\tNuctideSequence\tFullNuctideName\tOligoName\tTransitionGroupId\tTransitionId\tCompoundName\n"
    "500.25\t400.15\t1000.0\t25.5\tPEPTIDE\tPEPTIDE\ttest_protein\t\t\t\tpep_group_1\ttransition_1\t\n"
    "600.35\t500.25\t1200.0\t30.5\t\t\t\tAUGCGAUU\tAUG[Cm]CGAUU\ttest_oligo\tnuc_group_1\ttransition_2\t\n";
  
  String temp_filename = File::getTemporaryFile();
  std::ofstream temp_out_file(temp_filename.c_str());
  temp_out_file << tsv_content;
  temp_out_file.close();

  TransitionTSVFile tsv_file;
  TargetedExperiment exp;
  
  tsv_file.convertTSVToTargetedExperiment(temp_filename.c_str(), FileTypes::TSV, exp);
  
  // Verify mixed content
  TEST_EQUAL(exp.getTransitions().size(), 2)
  TEST_EQUAL(exp.getPeptides().size(), 1)
  TEST_EQUAL(exp.getNuctides().size(), 1)
  TEST_EQUAL(exp.getProteins().size(), 1)
  TEST_EQUAL(exp.getOligos().size(), 1)
  
  // Test peptide transition
  const ReactionMonitoringTransition& pep_transition = exp.getTransitions()[0];
  TEST_EQUAL(pep_transition.getTransGroupRef(), "pep_group_1")

  // Test nucleotide transition
  const ReactionMonitoringTransition& nuc_transition = exp.getTransitions()[1];
  TEST_EQUAL(nuc_transition.getTransGroupRef(), "nuc_group_1")
  
  File::remove(temp_filename);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



