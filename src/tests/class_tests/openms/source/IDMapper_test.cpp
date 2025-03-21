// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <iostream>

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

///////////////////////////

using namespace OpenMS;
class IDMapper2 : public IDMapper
{
  public:
    double getAbsoluteMZTolerance2_(const double mz)
    {
      return getAbsoluteMZTolerance_(mz);
    }

    bool isMatch2_(const double rt_distance, const double mz_theoretical, const double mz_observed)
    {
      return isMatch_(rt_distance, mz_theoretical, mz_observed);
    }

};

START_TEST(IDMapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


using namespace std;

IDMapper* ptr = nullptr;
IDMapper* nullPointer = nullptr;
START_SECTION((IDMapper()))
  ptr = new IDMapper();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~IDMapper()))
  delete ptr;
END_SECTION

START_SECTION((IDMapper(const IDMapper& cp)))
{
  IDMapper mapper;
  Param p = mapper.getParameters();
  p.setValue("rt_tolerance", 0.5);
  p.setValue("mz_tolerance", 0.05);
  p.setValue("mz_measure","ppm");
  mapper.setParameters(p);
  IDMapper m2(mapper);
  TEST_EQUAL(m2.getParameters(), p);
}
END_SECTION


START_SECTION((IDMapper& operator = (const IDMapper& rhs)))
  IDMapper mapper;
  Param p = mapper.getParameters();
  p.setValue("rt_tolerance", 0.5);
  p.setValue("mz_tolerance", 0.05);
  p.setValue("mz_measure", "ppm");
  mapper.setParameters(p);
  IDMapper m2=mapper;
  TEST_EQUAL(m2.getParameters(), p);
END_SECTION

START_SECTION((template <typename PeakType> void annotate(MSExperiment<PeakType>& map, FeatureMap fmap, const bool clear_ids = false, const bool mapMS1 = false)))
  // create id
  FeatureMap fm;
  Feature f;
  f.setMZ(900.0);
  f.setRT(9.0);
  std::vector< PeptideIdentification > pids;
  PeptideIdentification pid;
  pid.setIdentifier("myID");
  pid.setHits(std::vector<PeptideHit>(4));
  pids.push_back(pid); // without MZ&RT for PID (take feature instead)
  pid.setMZ(800.0);
  pid.setRT(9.05);
  pids.push_back(pid); // with MZ&RT from PID
  f.setPeptideIdentifications(pids);
  fm.push_back(f);
  std::vector< ProteinIdentification > prids(2);
  fm.setProteinIdentifications(prids);

  // create experiment
  PeakMap experiment;
  MSSpectrum spectrum;
  Precursor precursor;
  precursor.setMZ(0);
  spectrum.setRT(8.9);
  experiment.addSpectrum(spectrum);
  experiment[0].getPrecursors().push_back(precursor);
  precursor.setMZ(20);
  spectrum.setRT(9.1);
  experiment.addSpectrum(spectrum);
  experiment[1].getPrecursors().push_back(precursor);
  precursor.setMZ(11);
  spectrum.setRT(12.0);
  experiment.addSpectrum(spectrum);
  experiment[2].getPrecursors().push_back(precursor);

  // map
  IDMapper mapper;
  Param p = mapper.getParameters();
  p.setValue("rt_tolerance", 0.3);
  p.setValue("mz_tolerance", 0.05);
  p.setValue("mz_measure", "Da");
  p.setValue("ignore_charge", "true");
  mapper.setParameters(p);


  mapper.annotate(experiment, fm, true, true);

  //test
  TEST_EQUAL(experiment.getProteinIdentifications().size(), 2)
  //scan 1
  TEST_EQUAL(experiment[0].getPeptideIdentifications().size(), 2)
  //scan 2
  TEST_EQUAL(experiment[1].getPeptideIdentifications().size(), 2)
  ABORT_IF(experiment[1].getPeptideIdentifications().size() != 2)
  TEST_EQUAL(experiment[1].getPeptideIdentifications()[0].getHits().size(), 4)
  TEST_EQUAL(experiment[1].getPeptideIdentifications()[0].getMZ(), 900.0)
  TEST_EQUAL(experiment[1].getPeptideIdentifications()[1].getHits().size(), 4)
  TEST_EQUAL(experiment[1].getPeptideIdentifications()[1].getMZ(), 800.0)
  //scan 3
  TEST_EQUAL(experiment[2].getPeptideIdentifications().size(), 0)

  mapper.annotate(experiment, fm, true, false); // no MS1 mapping. MZ threshold never fulfilled
  //test
  TEST_EQUAL(experiment.getProteinIdentifications().size(), 2)
  //scan 1
  TEST_EQUAL(experiment[0].getPeptideIdentifications().size(), 0)
  //scan 2
  TEST_EQUAL(experiment[1].getPeptideIdentifications().size(), 0)
  //scan 3
  TEST_EQUAL(experiment[2].getPeptideIdentifications().size(), 0)

END_SECTION

START_SECTION((template <typename PeakType> void annotate(MSExperiment<PeakType>& map, const std::vector<PeptideIdentification>& peptide_ids, const std::vector<ProteinIdentification>& protein_ids, const bool clear_ids = false, const bool mapMS1 = false)))
  // load id
  vector<PeptideIdentification> identifications;
  vector<ProteinIdentification> protein_identifications;
  String document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_1.idXML"), protein_identifications, identifications, document_id);

  vector<PeptideIdentification> identifications2;
  vector<ProteinIdentification> protein_identifications2;
  String document_id2;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_7.idXML"), protein_identifications2, identifications2, document_id2);

  TEST_EQUAL(identifications.size(),3)
  TEST_EQUAL(identifications[0].getHits().size(), 2)
  TEST_EQUAL(identifications[1].getHits().size(), 1)
  TEST_EQUAL(identifications[2].getHits().size(), 2)
  TEST_EQUAL(protein_identifications.size(),1)
  TEST_EQUAL(protein_identifications[0].getHits().size(), 2)

  TEST_EQUAL(identifications2.size(),3)
  TEST_EQUAL(identifications2[0].getHits().size(), 2)
  TEST_EQUAL(identifications2[1].getHits().size(), 1)
  TEST_EQUAL(identifications2[2].getHits().size(), 2)
  TEST_EQUAL(protein_identifications2.size(),1)
  TEST_EQUAL(protein_identifications2[0].getHits().size(), 2)

  // TEST RT MAPPING
  // create experiment
  PeakMap experiment;
  MSSpectrum spectrum;
  Precursor precursor;
  precursor.setMZ(0);
  spectrum.setRT(60);
  experiment.addSpectrum(spectrum);
  experiment[0].getPrecursors().push_back(precursor);
  precursor.setMZ(20);
  spectrum.setRT(181);
  experiment.addSpectrum(spectrum);
  experiment[1].getPrecursors().push_back(precursor);
  precursor.setMZ(11);
  spectrum.setRT(120.0001);
  experiment.addSpectrum(spectrum);
  experiment[2].getPrecursors().push_back(precursor);

  //map
  IDMapper mapper;
  Param p = mapper.getParameters();
  p.setValue("rt_tolerance", 0.5);
  p.setValue("mz_tolerance", 0.05);
  p.setValue("mz_measure","Da");
  p.setValue("ignore_charge", "true");
  mapper.setParameters(p);

  mapper.annotate(experiment, identifications, protein_identifications);

  //test
  TEST_EQUAL(experiment.getProteinIdentifications().size(), 1)
  TEST_EQUAL(experiment.getProteinIdentifications()[0].getHits().size(),2)
  TEST_EQUAL(experiment.getProteinIdentifications()[0].getHits()[0].getAccession(),"ABCDE")
  TEST_EQUAL(experiment.getProteinIdentifications()[0].getHits()[1].getAccession(),"FGHIJ")
  //scan 1
  TEST_EQUAL(experiment[0].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(experiment[0].getPeptideIdentifications()[0].getHits().size(), 2)
  TEST_EQUAL(experiment[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("LHASGITVTEIPVTATNFK"))
  TEST_EQUAL(experiment[0].getPeptideIdentifications()[0].getHits()[1].getSequence(), AASequence::fromString("MRSLGYVAVISAVATDTDK"))
  //scan 2
  TEST_EQUAL(experiment[1].getPeptideIdentifications().size(), 0)
  //scan 3
  TEST_EQUAL(experiment[2].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(experiment[2].getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(experiment[2].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("HSKLSAK"))

  //-----------------------------------------------------------------------------------
  // TEST NATIVE_ID MAPPING

  // create experiment
  PeakMap experiment2;
  MSSpectrum spectrum2;
  Precursor precursor2;
  precursor2.setMZ(0);
  spectrum2.setRT(60);
  spectrum2.setNativeID("spectrum=1234");
  experiment2.addSpectrum(spectrum2);
  experiment2[0].getPrecursors().push_back(precursor2);
  precursor2.setMZ(20);
  spectrum2.setRT(181);
  spectrum2.setNativeID("spectrum=6666");
  experiment2.addSpectrum(spectrum2);
  experiment2[1].getPrecursors().push_back(precursor2);
  precursor2.setMZ(11);
  spectrum2.setRT(120.0001);
  spectrum2.setNativeID("spectrum=4321");
  experiment2.addSpectrum(spectrum2);
  experiment2[2].getPrecursors().push_back(precursor2);

  IDMapper mapper2;
  Param p2 = mapper2.getParameters();
  p2.setValue("rt_tolerance", 0.5);
  p2.setValue("mz_tolerance", 0.05);
  p2.setValue("mz_measure","Da");
  p2.setValue("ignore_charge", "true");
  mapper2.setParameters(p2);

  mapper2.annotate(experiment2, identifications2, protein_identifications2);

  //test
  TEST_EQUAL(experiment2.getProteinIdentifications().size(), 1)
  TEST_EQUAL(experiment2.getProteinIdentifications()[0].getHits().size(),2)
  TEST_EQUAL(experiment2.getProteinIdentifications()[0].getHits()[0].getAccession(),"ABCDE")
  TEST_EQUAL(experiment2.getProteinIdentifications()[0].getHits()[1].getAccession(),"FGHIJ")
  //scan 1
  TEST_EQUAL(experiment2[0].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(experiment2[0].getPeptideIdentifications()[0].getHits().size(), 2)
  TEST_EQUAL(experiment2[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("LHASGITVTEIPVTATNFK"))
  TEST_EQUAL(experiment2[0].getPeptideIdentifications()[0].getHits()[1].getSequence(), AASequence::fromString("MRSLGYVAVISAVATDTDK"))
  //scan 2
  TEST_EQUAL(experiment2[1].getPeptideIdentifications().size(), 0)
  //scan 3
  TEST_EQUAL(experiment2[2].getPeptideIdentifications().size(), 1)
  TEST_EQUAL(experiment2[2].getPeptideIdentifications()[0].getHits().size(), 1)
  TEST_EQUAL(experiment2[2].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("HSKLSAK"))

END_SECTION


START_SECTION((template < typename FeatureType > void annotate(FeatureMap< FeatureType > &map, const std::vector< PeptideIdentification > &ids, const std::vector< ProteinIdentification > &protein_ids, bool use_centroid_rt=false, bool use_centroid_mz=false)))
{
  //load id data
  vector<PeptideIdentification> identifications;
  vector<ProteinIdentification> protein_identifications;
  String document_id;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_2.idXML"), protein_identifications, identifications, document_id);

  //--------------------------------------------------------------------------------------
  //TEST MAPPING TO CONVEX HULLS
  FeatureMap fm;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_2.featureXML"), fm);

  IDMapper mapper;
  Param p = mapper.getParameters();
  p.setValue("rt_tolerance", 0.0);
  p.setValue("mz_tolerance", 0.0);
  p.setValue("mz_measure","Da");
  p.setValue("ignore_charge", "true");
  mapper.setParameters(p);

  mapper.annotate(fm,identifications,protein_identifications);

  //test protein ids
  TEST_EQUAL(fm.getProteinIdentifications().size(),1)
  TEST_EQUAL(fm.getProteinIdentifications()[0].getHits().size(),2)
  TEST_EQUAL(fm.getProteinIdentifications()[0].getHits()[0].getAccession(),"ABCDE")
  TEST_EQUAL(fm.getProteinIdentifications()[0].getHits()[1].getAccession(),"FGHIJ")

  //test peptide ids
  TEST_EQUAL(fm[0].getPeptideIdentifications().size(),7)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits().size(),1)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[1].getHits().size(),1)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[2].getHits().size(),1)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[3].getHits().size(),1)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[4].getHits().size(),1)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[5].getHits().size(),1)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[6].getHits().size(),1)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits()[0].getSequence(),AASequence::fromString("A"))
  TEST_EQUAL(fm[0].getPeptideIdentifications()[1].getHits()[0].getSequence(),AASequence::fromString("K"))
  TEST_EQUAL(fm[0].getPeptideIdentifications()[2].getHits()[0].getSequence(),AASequence::fromString("C"))
  TEST_EQUAL(fm[0].getPeptideIdentifications()[3].getHits()[0].getSequence(),AASequence::fromString("D"))
  TEST_EQUAL(fm[0].getPeptideIdentifications()[4].getHits()[0].getSequence(),AASequence::fromString("E"))
  TEST_EQUAL(fm[0].getPeptideIdentifications()[5].getHits()[0].getSequence(),AASequence::fromString("F"))
  TEST_EQUAL(fm[0].getPeptideIdentifications()[6].getHits()[0].getSequence(),AASequence::fromString("I"))

  //test unassigned peptide ids
  TEST_EQUAL(fm.getUnassignedPeptideIdentifications().size(),3)
  TEST_EQUAL(fm.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence(),AASequence::fromString("G"))
  TEST_EQUAL(fm.getUnassignedPeptideIdentifications()[1].getHits()[0].getSequence(),AASequence::fromString("H"))
  TEST_EQUAL(fm.getUnassignedPeptideIdentifications()[2].getHits()[0].getSequence(),AASequence::fromString("L"))

  //--------------------------------------------------------------------------------------
  //TEST MAPPING TO CENTROIDS
  FeatureMap fm2;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_2.featureXML"), fm2);
  p.setValue("rt_tolerance", 4.0);
  p.setValue("mz_tolerance", 1.5);
  p.setValue("mz_measure","Da");
  p.setValue("ignore_charge", "true");
  mapper.setParameters(p);

  mapper.annotate(fm2,identifications,protein_identifications, true, true);

  //test protein ids
  TEST_EQUAL(fm2.getProteinIdentifications().size(),1)
  TEST_EQUAL(fm2.getProteinIdentifications()[0].getHits().size(),2)
  TEST_EQUAL(fm2.getProteinIdentifications()[0].getHits()[0].getAccession(),"ABCDE")
  TEST_EQUAL(fm2.getProteinIdentifications()[0].getHits()[1].getAccession(),"FGHIJ")

  //test peptide ids
  TEST_EQUAL(fm2[0].getPeptideIdentifications().size(),2)
  TEST_EQUAL(fm2[0].getPeptideIdentifications()[0].getHits().size(),1)
  TEST_EQUAL(fm2[0].getPeptideIdentifications()[1].getHits().size(),1)
  TEST_EQUAL(fm2[0].getPeptideIdentifications()[0].getHits()[0].getSequence(),AASequence::fromString("A"))
  TEST_EQUAL(fm2[0].getPeptideIdentifications()[1].getHits()[0].getSequence(),AASequence::fromString("K"))

  //test unassigned peptide ids
  TEST_EQUAL(fm2.getUnassignedPeptideIdentifications().size(),8)
  TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence(),AASequence::fromString("C"))
  TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[1].getHits()[0].getSequence(),AASequence::fromString("D"))
  TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[2].getHits()[0].getSequence(),AASequence::fromString("E"))
  TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[3].getHits()[0].getSequence(),AASequence::fromString("F"))
  TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[4].getHits()[0].getSequence(),AASequence::fromString("G"))
  TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[5].getHits()[0].getSequence(),AASequence::fromString("H"))
  TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[6].getHits()[0].getSequence(),AASequence::fromString("I"))
  TEST_EQUAL(fm2.getUnassignedPeptideIdentifications()[7].getHits()[0].getSequence(),AASequence::fromString("L"))

  // ******* test charge-specific matching *******

  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_2.featureXML"), fm);

  p.setValue("rt_tolerance", 0.0);
  p.setValue("mz_tolerance", 0.0);
  p.setValue("mz_measure", "Da");
  p.setValue("ignore_charge", "false");
  mapper.setParameters(p);

  mapper.annotate(fm, identifications, protein_identifications);

  //test protein ids
  TEST_EQUAL(fm.getProteinIdentifications().size(), 1)
  TEST_EQUAL(fm.getProteinIdentifications()[0].getHits().size(), 2)
  TEST_EQUAL(fm.getProteinIdentifications()[0].getHits()[0].getAccession(),
             "ABCDE")
  TEST_EQUAL(fm.getProteinIdentifications()[0].getHits()[1].getAccession(),
             "FGHIJ")

  //test peptide ids
  TEST_EQUAL(fm[0].getPeptideIdentifications().size(), 3)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits().size(),1)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[1].getHits().size(),1)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[2].getHits().size(),1)
  TEST_EQUAL(fm[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("A"))
  TEST_EQUAL(fm[0].getPeptideIdentifications()[1].getHits()[0].getSequence(), AASequence::fromString("K"))
  TEST_EQUAL(fm[0].getPeptideIdentifications()[2].getHits()[0].getSequence(), AASequence::fromString("C"))

  //test unassigned peptide ids
  TEST_EQUAL(fm.getUnassignedPeptideIdentifications().size(), 7)


  // ******* PPM test *******
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_4.idXML"), protein_identifications, identifications);

  FeatureMap fm_ppm;
  FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_4.featureXML"), fm_ppm);
  p.setValue("rt_tolerance", 4.0);
  p.setValue("mz_tolerance", 3.0);
  p.setValue("mz_measure","ppm");
  p.setValue("ignore_charge", "true");
  mapper.setParameters(p);

  mapper.annotate(fm_ppm,identifications,protein_identifications);

  //test peptide ids
  TEST_EQUAL(fm_ppm[0].getPeptideIdentifications().size(),1)
  TEST_EQUAL(fm_ppm[0].getPeptideIdentifications()[0].getHits().size(),2)
  TEST_EQUAL(fm_ppm[0].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("LHASGITVTEIPVTATNFK"))

  TEST_EQUAL(fm_ppm[1].getPeptideIdentifications().size(),0)

  TEST_EQUAL(fm_ppm[2].getPeptideIdentifications().size(),1)
  TEST_EQUAL(fm_ppm[2].getPeptideIdentifications()[0].getHits().size(),1)
  TEST_EQUAL(fm_ppm[2].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("HSKLSAK"))

  TEST_EQUAL(fm_ppm[3].getPeptideIdentifications().size(),0)

  TEST_EQUAL(fm_ppm[4].getPeptideIdentifications().size(),1)
  TEST_EQUAL(fm_ppm[4].getPeptideIdentifications()[0].getHits().size(),2)
  TEST_EQUAL(fm_ppm[4].getPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("RASNSPQDPQSATAHSFR"))

  TEST_EQUAL(fm_ppm[5].getPeptideIdentifications().size(),0)


  TEST_EQUAL(fm_ppm.getUnassignedPeptideIdentifications().size(),2)
  TEST_EQUAL(fm_ppm.getUnassignedPeptideIdentifications()[0].getHits()[0].getSequence(), AASequence::fromString("DEAD"))
  TEST_EQUAL(fm_ppm.getUnassignedPeptideIdentifications()[0].getHits()[1].getSequence(), AASequence::fromString("DEADA"))
  TEST_EQUAL(fm_ppm.getUnassignedPeptideIdentifications()[1].getHits()[0].getSequence(), AASequence::fromString("DEADAA"))
  TEST_EQUAL(fm_ppm.getUnassignedPeptideIdentifications()[1].getHits()[1].getSequence(), AASequence::fromString("DEADAAA"))
}
END_SECTION


START_SECTION((void annotate(ConsensusMap& map, const std::vector<PeptideIdentification>& ids, const std::vector<ProteinIdentification>& protein_ids, bool measure_from_subelements=false)))
{
  IDMapper mapper;
  Param p = mapper.getParameters();
  p.setValue("mz_tolerance", 0.01);
  p.setValue("mz_measure","Da");
  p.setValue("ignore_charge", "true");
  mapper.setParameters(p);

  TOLERANCE_ABSOLUTE(0.01);

  std::vector<ProteinIdentification> protein_ids;
  std::vector<ProteinIdentification> protein_ids2;
  std::vector<PeptideIdentification> peptide_ids;
  std::vector<PeptideIdentification> peptide_ids2;
  String document_id;
  String document_id2;
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_3.idXML"), protein_ids, peptide_ids, document_id);
  IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IDMapper_5.idXML"), protein_ids2, peptide_ids2, document_id2);

  ConsensusXMLFile cons_file;

  {
    std::string tmp_filename;
    NEW_TMP_FILE(tmp_filename);
    ConsensusMap cons_map;
    cons_file.load(OPENMS_GET_TEST_DATA_PATH("IDMapper_3.consensusXML"), cons_map);
    mapper.annotate(cons_map, peptide_ids, protein_ids);
    cons_file.store(tmp_filename,cons_map);
    WHITELIST("<?xml-stylesheet, date=");
    TEST_FILE_SIMILAR(tmp_filename,OPENMS_GET_TEST_DATA_PATH("IDMapper_3_out1.consensusXML"));
  }

  {
    std::string tmp_filename;
    NEW_TMP_FILE(tmp_filename);
    ConsensusMap cons_map;
    cons_file.load(OPENMS_GET_TEST_DATA_PATH("IDMapper_3.consensusXML"), cons_map);
    mapper.annotate(cons_map, peptide_ids, protein_ids, true);
    cons_file.store(tmp_filename,cons_map);
    WHITELIST("<?xml-stylesheet, date=");
    TEST_FILE_SIMILAR(tmp_filename,OPENMS_GET_TEST_DATA_PATH("IDMapper_3_out2.consensusXML"));
  }

  {
    IDMapper mapper5;
    Param p5 = mapper5.getParameters();
    p5.setValue("rt_tolerance", 20.0);
    p5.setValue("mz_tolerance", 20.0);
    p5.setValue("mz_measure","ppm");
    p5.setValue("ignore_charge", "true");
    p5.setValue("consensus:use_subelements", "true");
    p5.setValue("consensus:annotate_ids_with_subelements", "true");
    mapper5.setParameters(p5);

    std::string tmp_filename;
    NEW_TMP_FILE(tmp_filename);
    ConsensusMap cons_map;
    cons_file.load(OPENMS_GET_TEST_DATA_PATH("IDMapper_5.consensusXML"), cons_map);
    mapper5.annotate(cons_map, peptide_ids2, protein_ids2, true, true);
    cons_file.store(tmp_filename,cons_map);
    WHITELIST("<?xml-stylesheet, date=");
    TEST_FILE_SIMILAR(tmp_filename,OPENMS_GET_TEST_DATA_PATH("IDMapper_5_out1.consensusXML"));
  }

  // check charge-specific matching:
  {
    ConsensusMap cm;
    cm.resize(1);
    cm[0].setRT(4101.48);
    cm[0].setMZ(117.1);
    cm[0].setCharge(2);

    mapper.annotate(cm, peptide_ids, protein_ids);

    TEST_EQUAL(cm[0].getPeptideIdentifications().size(), 1);
    TEST_EQUAL(cm[0].getPeptideIdentifications()[0].getHits()[0].getSequence(),
    AASequence::fromString("ACSF"));
    TEST_EQUAL(cm.getUnassignedPeptideIdentifications().size(),
    peptide_ids.size() - 1);

    cm[0].getPeptideIdentifications().clear();
    cm.getUnassignedPeptideIdentifications().clear();
    p.setValue("ignore_charge", "false");
    mapper.setParameters(p);
    mapper.annotate(cm, peptide_ids, protein_ids);
    TEST_EQUAL(cm[0].getPeptideIdentifications().size(), 0);
    TEST_EQUAL(cm.getUnassignedPeptideIdentifications().size(),
               peptide_ids.size());
  }

  // annotation of precursors without id
  IDMapper mapper6;
  p = mapper6.getParameters();
  p.setValue("mz_tolerance", 0.01);
  p.setValue("mz_measure","Da");
  p.setValue("ignore_charge", "true");
  mapper6.setParameters(p);

  TOLERANCE_ABSOLUTE(0.01);

  PeakMap experiment;
  MSSpectrum spectrum;

  // match exactly to the first 10 consensusXML centroids 
  double mzs[10] = { 426.849, 405.85, 506.815, 484.83, 496.244, 430.212, 446.081, 453.233, 400.172, 437.227 }; 
  double rts[10] = { 306.58, 306.58, 312.738, 312.738, 3112.53, 3840.95, 3849.22, 3870.67, 3880.9, 3892.26}; 

  for (Size i = 0; i != 10; ++i)
  {
    vector<Precursor> precursors;
    Precursor prec;
    prec.setMZ(mzs[i]);
    precursors.push_back(prec);
    spectrum.setRT(rts[i]);
    spectrum.setPrecursors(precursors);
    experiment.addSpectrum(spectrum);
  }

  {
    std::string tmp_filename;
    NEW_TMP_FILE(tmp_filename);
    ConsensusMap cons_map;
    cons_file.load(OPENMS_GET_TEST_DATA_PATH("IDMapper_3.consensusXML"), cons_map);
    mapper6.annotate(cons_map, vector<PeptideIdentification>(), vector<ProteinIdentification>(), false, false, experiment);
    cons_file.store(tmp_filename, cons_map);
    WHITELIST("<?xml-stylesheet, date=");
    TEST_FILE_SIMILAR(tmp_filename, OPENMS_GET_TEST_DATA_PATH("IDMapper_6_out1.consensusXML"));
  }

  experiment.clear(true);


  // only 5 should be in the 0.01 Da tolerance (every second entry is too much off)
  double mzs_5_mismatch[10] = { 426.85899, 405, 506.815, 484.85, 496.244, 430, 446.081, 453, 400.172, 437.239 }; 

  for (Size i = 0; i != 10; ++i)
  {
    vector<Precursor> precursors;
    Precursor prec;
    prec.setMZ(mzs_5_mismatch[i]);
    precursors.push_back(prec);
    spectrum.setRT(rts[i]);
    spectrum.setPrecursors(precursors);
    experiment.addSpectrum(spectrum);
  }

  {
    std::string tmp_filename;
    NEW_TMP_FILE(tmp_filename);
    ConsensusMap cons_map;
    cons_file.load(OPENMS_GET_TEST_DATA_PATH("IDMapper_3.consensusXML"), cons_map);
    mapper6.annotate(cons_map, vector<PeptideIdentification>(), vector<ProteinIdentification>(), false, false, experiment);
    cons_file.store(tmp_filename, cons_map);
    WHITELIST("<?xml-stylesheet, date=");
    TEST_FILE_SIMILAR(tmp_filename, OPENMS_GET_TEST_DATA_PATH("IDMapper_6_out2.consensusXML"));
  }

  // check mappings of multiple precursors to one consensus feature
  experiment.clear(true);
  double rts_multiple[5] = { 306.58, 305.58, 307.58, 304.58, 308.58 };
  for (Size i = 0; i != 5; ++i)
  {
    vector<Precursor> precursors;
    Precursor prec;
    prec.setMZ(426.849);
    precursors.push_back(prec);
    spectrum.setRT(rts_multiple[i]);
    spectrum.setPrecursors(precursors);
    experiment.addSpectrum(spectrum);
  }

  {
    std::string tmp_filename;
    NEW_TMP_FILE(tmp_filename);
    ConsensusMap cons_map;
    cons_file.load(OPENMS_GET_TEST_DATA_PATH("IDMapper_3.consensusXML"), cons_map);
    mapper6.annotate(cons_map, vector<PeptideIdentification>(), vector<ProteinIdentification>(), false, false, experiment);
    cons_file.store(tmp_filename, cons_map);
    WHITELIST("<?xml-stylesheet, date=");
    TEST_FILE_SIMILAR(tmp_filename, OPENMS_GET_TEST_DATA_PATH("IDMapper_6_out3.consensusXML"));
  }
}
END_SECTION

START_SECTION([EXTRA] double getAbsoluteMZTolerance_(const double mz) const)
  IDMapper2 mapper;
  Param p = mapper.getParameters();
  p.setValue("mz_tolerance", 1.0);
  mapper.setParameters(p);
  TEST_REAL_SIMILAR(mapper.getAbsoluteMZTolerance2_(1000), 0.001)
  p.setValue("mz_tolerance", 3.0);
  mapper.setParameters(p);
  TEST_REAL_SIMILAR(mapper.getAbsoluteMZTolerance2_(1000), 0.003)
  p.setValue("mz_measure","Da");
  mapper.setParameters(p);
  TEST_REAL_SIMILAR(mapper.getAbsoluteMZTolerance2_(1000), 3)
END_SECTION

START_SECTION([EXTRA] bool isMatch_(const double rt_distance, const double mz_theoretical, const double mz_observed) const)
  IDMapper2 mapper;
  TEST_EQUAL(mapper.isMatch2_(1, 1000, 1000.001), true)
  Param p = mapper.getParameters();
  p.setValue("mz_tolerance", 3.0);
  mapper.setParameters(p);
  TEST_EQUAL(mapper.isMatch2_(4, 1000, 1000.0028), true)
  TEST_EQUAL(mapper.isMatch2_(4, 1000, 1000.004), false)
  TEST_EQUAL(mapper.isMatch2_(4, 1000, 999.9972), true)
  TEST_EQUAL(mapper.isMatch2_(4, 1000, 999.996), false)
  p.setValue("mz_measure","Da");
  mapper.setParameters(p);
  TEST_EQUAL(mapper.isMatch2_(5, 999, 1002), true)
  TEST_EQUAL(mapper.isMatch2_(5, 999, 1002.1), false)
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
