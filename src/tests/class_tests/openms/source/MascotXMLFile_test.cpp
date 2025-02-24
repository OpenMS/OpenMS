// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>


///////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

///////////////////////////

START_TEST(MascotXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MascotXMLFile xml_file;
MascotXMLFile* ptr;
ProteinIdentification protein_identification;
vector<PeptideIdentification> peptide_identifications;
vector<PeptideIdentification> peptide_identifications2;
DateTime date;
PeptideHit peptide_hit;
vector<String> references;

date.set("2006-03-09 11:31:52");

MascotXMLFile* nullPointer = nullptr;
START_SECTION((MascotXMLFile()))
  ptr = new MascotXMLFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
  delete ptr;
END_SECTION

START_SECTION((static void initializeLookup(SpectrumMetaDataLookup& lookup, PeakMap& experiment, const String& scan_regex = "")))
{
  PeakMap exp;
  exp.getSpectra().resize(1);
  SpectrumMetaDataLookup lookup;
  xml_file.initializeLookup(lookup, exp);
  TEST_EQUAL(lookup.empty(), false);
}
END_SECTION

START_SECTION((void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, SpectrumMetaDataLookup& lookup)))
{
  SpectrumMetaDataLookup lookup;
  xml_file.load(OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_1.mascotXML"),
                protein_identification, peptide_identifications, lookup);

  {
    ProteinIdentification::SearchParameters search_parameters = protein_identification.getSearchParameters();
    TEST_EQUAL(search_parameters.missed_cleavages, 1);
    TEST_EQUAL(search_parameters.taxonomy, ". . Eukaryota (eucaryotes)");
    TEST_EQUAL(search_parameters.mass_type, ProteinIdentification::AVERAGE);
    TEST_EQUAL(search_parameters.db, "MSDB_chordata");
    TEST_EQUAL(search_parameters.db_version, "MSDB_chordata_20070910.fasta");
    TEST_EQUAL(search_parameters.fragment_mass_tolerance, 0.2);
    TEST_EQUAL(search_parameters.precursor_mass_tolerance, 1.4);
    TEST_EQUAL(search_parameters.fragment_mass_tolerance_ppm, false);
    TEST_EQUAL(search_parameters.precursor_mass_tolerance_ppm, false);
    TEST_EQUAL(search_parameters.charges, "1+, 2+ and 3+");
    TEST_EQUAL(search_parameters.fixed_modifications.size(), 4);
    TEST_EQUAL(search_parameters.fixed_modifications[0], "Carboxymethyl (C)");
    TEST_EQUAL(search_parameters.fixed_modifications[1], "Deamidated (N)");
    TEST_EQUAL(search_parameters.fixed_modifications[2], "Deamidated (Q)");
    TEST_EQUAL(search_parameters.fixed_modifications[3], "Guanidinyl (K)");
    TEST_EQUAL(search_parameters.variable_modifications.size(), 3);
    TEST_EQUAL(search_parameters.variable_modifications[0], "Acetyl (Protein N-term)");
    TEST_EQUAL(search_parameters.variable_modifications[1], "Biotin (K)");
    TEST_EQUAL(search_parameters.variable_modifications[2], "Carbamyl (K)");
    TEST_EQUAL(peptide_identifications.size(), 3);
    TOLERANCE_ABSOLUTE(0.0001);
    TEST_REAL_SIMILAR(peptide_identifications[0].getMZ(), 789.83);
    TEST_REAL_SIMILAR(peptide_identifications[1].getMZ(), 135.29);
    TEST_REAL_SIMILAR(peptide_identifications[2].getMZ(), 982.58);
    TOLERANCE_ABSOLUTE(0.00001);
    TEST_EQUAL(protein_identification.getHits().size(), 2);
    TEST_EQUAL(protein_identification.getHits()[0].getAccession(), "AAN17824");
    TEST_EQUAL(protein_identification.getHits()[1].getAccession(), "GN1736");
    TEST_REAL_SIMILAR(protein_identification.getHits()[0].getScore(), 619);
    TEST_REAL_SIMILAR(protein_identification.getHits()[1].getScore(), 293);
    TEST_EQUAL(protein_identification.getScoreType(), "Mascot");
    TEST_EQUAL(protein_identification.getDateTime().get(), "2006-03-09 11:31:52");

    TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 31.8621);
    TEST_EQUAL(peptide_identifications[0].getHits().size(), 2);

    peptide_hit = peptide_identifications[0].getHits()[0];
    set<String> ref_set = peptide_hit.extractProteinAccessionsSet();
    vector<String> references(ref_set.begin(), ref_set.end());
    TEST_EQUAL(references.size(), 2);
    TEST_EQUAL(references[0], "AAN17824");
    TEST_EQUAL(references[1], "GN1736");
    peptide_hit = peptide_identifications[0].getHits()[1];
    ref_set = peptide_hit.extractProteinAccessionsSet();
    references = vector<String>(ref_set.begin(), ref_set.end());
    TEST_EQUAL(references.size(), 1);
    TEST_EQUAL(references[0], "AAN17824");
    peptide_hit = peptide_identifications[1].getHits()[0];
    ref_set = peptide_hit.extractProteinAccessionsSet();
    references = vector<String>(ref_set.begin(), ref_set.end());
    TEST_EQUAL(references.size(), 1);
    TEST_EQUAL(references[0], "GN1736");

    TEST_EQUAL(peptide_identifications[1].getHits().size(), 1);
    TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), 33.85);
    TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[1].getScore(), 33.12);
    TEST_REAL_SIMILAR(peptide_identifications[1].getHits()[0].getScore(), 43.9);
    TEST_EQUAL(peptide_identifications[0].getScoreType(), "Mascot");
    TEST_EQUAL(peptide_identifications[1].getScoreType(), "Mascot");
    TEST_EQUAL(protein_identification.getDateTime() == date, true);
    TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), AASequence::fromString("LHASGITVTEIPVTATN(MOD:00565)FK(MOD:00445)"));
    TEST_EQUAL(peptide_identifications[0].getHits()[1].getSequence(), AASequence::fromString("MRSLGYVAVISAVATDTDK(MOD:00445)"));
    TEST_EQUAL(peptide_identifications[1].getHits()[0].getSequence(), AASequence::fromString("HSK(MOD:00445)LSAK(MOD:00445)"));

    String identifier = protein_identification.getIdentifier();
    TEST_EQUAL(!identifier.empty(), true);
    for (Size i = 0; i < peptide_identifications.size(); ++i)
    {
      TEST_EQUAL(identifier, peptide_identifications[i].getIdentifier())
    }
  }


  /// for new MascotXML 2.1 as used by Mascot Server 2.3
  xml_file.load(OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_2.mascotXML"),
                protein_identification, peptide_identifications, lookup);
  {
    ProteinIdentification::SearchParameters search_parameters = protein_identification.getSearchParameters();
    TEST_EQUAL(search_parameters.missed_cleavages, 7);
    TEST_EQUAL(search_parameters.taxonomy, "All entries");
    TEST_EQUAL(search_parameters.mass_type, ProteinIdentification::MONOISOTOPIC);
    TEST_EQUAL(search_parameters.db, "IPI_human");
    TEST_EQUAL(search_parameters.db_version, "ipi.HUMAN.v3.61.fasta");
    TEST_EQUAL(search_parameters.fragment_mass_tolerance, 0.3);
    TEST_EQUAL(search_parameters.precursor_mass_tolerance, 3);
    TEST_EQUAL(search_parameters.fragment_mass_tolerance_ppm, false);
    TEST_EQUAL(search_parameters.precursor_mass_tolerance_ppm, false);
    TEST_EQUAL(search_parameters.charges, "");
    TEST_EQUAL(search_parameters.fixed_modifications.size(), 1);
    TEST_EQUAL(search_parameters.fixed_modifications[0], "Carbamidomethyl (C)");
    TEST_EQUAL(search_parameters.variable_modifications.size(), 3);
    TEST_EQUAL(search_parameters.variable_modifications[0], "Oxidation (M)");
    TEST_EQUAL(search_parameters.variable_modifications[1], "Acetyl (N-term)");
    TEST_EQUAL(search_parameters.variable_modifications[2], "Phospho (Y)");
  // not necessarily equal to numQueries as some hits might not be contained, e.g. peptide's might start with <peptide rank="10"...> so 9 peptides are missing
  // thus empty peptides are removed (see MascotXMLFile.cpp::load() ) after the handler() call
    TEST_EQUAL(peptide_identifications.size(), 1112);
    TOLERANCE_ABSOLUTE(0.0001);
    TEST_REAL_SIMILAR(peptide_identifications[0].getMZ(), 304.6967);
    TEST_REAL_SIMILAR(peptide_identifications[1].getMZ(), 314.1815);
    TEST_REAL_SIMILAR(peptide_identifications[1111].getMZ(), 583.7948);
    TOLERANCE_ABSOLUTE(0.00001);
    TEST_EQUAL(protein_identification.getHits().size(), 66);
    TEST_EQUAL(protein_identification.getHits()[0].getAccession(), "IPI00745872");
    TEST_EQUAL(protein_identification.getHits()[1].getAccession(), "IPI00908876");
    TEST_REAL_SIMILAR(protein_identification.getHits()[0].getScore(), 122);
    TEST_REAL_SIMILAR(protein_identification.getHits()[1].getScore(), 122);
    TEST_EQUAL(protein_identification.getScoreType(), "Mascot");
    TEST_EQUAL(protein_identification.getDateTime().get(), "2011-06-24 19:34:54");

    TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 5);
    TEST_EQUAL(peptide_identifications[0].getHits().size(), 1);

    peptide_hit = peptide_identifications[0].getHits()[0];
    vector<PeptideEvidence> pes = peptide_hit.getPeptideEvidences();
    TEST_EQUAL(pes.size(), 0);
    pes = peptide_identifications[34].getHits()[0].getPeptideEvidences();
    set<String> accessions = peptide_identifications[34].getHits()[0].extractProteinAccessionsSet();
    references = vector<String>(accessions.begin(), accessions.end()); // corresponds to <peptide query="35" ...>
    ABORT_IF(references.size() != 5);
    TEST_EQUAL(references[0], "IPI00022434");
    TEST_EQUAL(references[1], "IPI00384697");
    TEST_EQUAL(references[2], "IPI00745872");
    TEST_EQUAL(references[3], "IPI00878517");
    TEST_EQUAL(references[4], "IPI00908876");

    TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), 5.34);
    TEST_REAL_SIMILAR(peptide_identifications[49].getHits()[0].getScore(), 14.83);
    TEST_REAL_SIMILAR(peptide_identifications[49].getHits()[1].getScore(), 17.5);
    TEST_EQUAL(peptide_identifications[0].getScoreType(), "Mascot");
    TEST_EQUAL(peptide_identifications[1].getScoreType(), "Mascot");
    TEST_EQUAL(protein_identification.getDateTime().get() == "2011-06-24 19:34:54", true);
    TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), AASequence::fromString("VVFIK"));
    TEST_EQUAL(peptide_identifications[49].getHits()[0].getSequence(), AASequence::fromString("LASYLDK"));
    TEST_EQUAL(peptide_identifications[49].getHits()[1].getSequence(), AASequence::fromString("(Acetyl)AAFESDK"));
  //for (int i=520;i<540;++i) std::cerr << "i: " << i << " " << peptide_identifications[i].getHits()[0].getSequence() << "\n";
    TEST_EQUAL(peptide_identifications[522].getHits()[0].getSequence(), AASequence::fromString("(Acetyl)GALM(Oxidation)NEIQAAK"));
    TEST_EQUAL(peptide_identifications[67].getHits()[0].getSequence(), AASequence::fromString("SHY(Phospho)GGSR"));

    String identifier = protein_identification.getIdentifier();
    TEST_EQUAL(!identifier.empty(), true);
    for (Size i = 0; i < peptide_identifications.size(); ++i)
    {
      TEST_EQUAL(identifier, peptide_identifications[i].getIdentifier())
    }
  }

  xml_file.load(OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_3.mascotXML"),
                protein_identification, peptide_identifications, lookup);
  {
    std::vector<ProteinIdentification> pids;
    pids.push_back(protein_identification);
    String filename;
    NEW_TMP_FILE(filename)
    IdXMLFile().store(filename, pids, peptide_identifications);
    FuzzyStringComparator fuzzy;
    fuzzy.setWhitelist(ListUtils::create<String>("<?xml-stylesheet"));
    fuzzy.setAcceptableAbsolute(0.0001);
    bool result = fuzzy.compareFiles(filename, OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_out_3.idXML"));
    TEST_EQUAL(result, true);
  }
}
END_SECTION

START_SECTION((void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, std::map<String, std::vector<AASequence> >& peptides, SpectrumMetaDataLookup& lookup)))
  std::map<String, vector<AASequence> > modified_peptides;
  AASequence aa_sequence_1;
  AASequence aa_sequence_2;
  AASequence aa_sequence_3;
  vector<AASequence> temp;

  aa_sequence_1 = AASequence::fromString("LHASGITVTEIPVTATNFK");
  aa_sequence_1.setModification(16, "Deamidated");
  aa_sequence_2 = AASequence::fromString("MRSLGYVAVISAVATDTDK");
  aa_sequence_2.setModification(2, "Phospho");
  aa_sequence_3 = AASequence::fromString("HSKLSAK");
  aa_sequence_3.setModification(4, "Phospho");
  temp.push_back(aa_sequence_1);
  temp.push_back(aa_sequence_2);
  modified_peptides.insert(make_pair("789.83", temp));
  temp.clear();
  temp.push_back(aa_sequence_3);
  modified_peptides.insert(make_pair("135.29", temp));

  SpectrumMetaDataLookup lookup;
  xml_file.load(OPENMS_GET_TEST_DATA_PATH("MascotXMLFile_test_1.mascotXML"),
                protein_identification, peptide_identifications, 
                modified_peptides, lookup);

  TEST_EQUAL(peptide_identifications.size(), 3)
  TOLERANCE_ABSOLUTE(0.0001)
  TEST_REAL_SIMILAR(peptide_identifications[0].getMZ(), 789.83)
  TEST_REAL_SIMILAR(peptide_identifications[1].getMZ(), 135.29)
  TEST_REAL_SIMILAR(peptide_identifications[2].getMZ(), 982.58)
  TOLERANCE_ABSOLUTE(0.00001)
  TEST_EQUAL(protein_identification.getHits().size(), 2)
  TEST_EQUAL(protein_identification.getHits()[0].getAccession(), "AAN17824")
  TEST_EQUAL(protein_identification.getHits()[1].getAccession(), "GN1736")
  TEST_REAL_SIMILAR(protein_identification.getHits()[0].getScore(), 619)
  TEST_REAL_SIMILAR(protein_identification.getHits()[1].getScore(), 293)
  TEST_EQUAL(protein_identification.getScoreType(), "Mascot")
  TEST_EQUAL(protein_identification.getDateTime().get(), "2006-03-09 11:31:52")

  TEST_REAL_SIMILAR(peptide_identifications[0].getSignificanceThreshold(), 31.8621)
  TEST_EQUAL(peptide_identifications[0].getHits().size(), 2)

  peptide_hit = peptide_identifications[0].getHits()[0];  
  set<String> accessions = peptide_hit.extractProteinAccessionsSet();
  references = vector<String>(accessions.begin(), accessions.end());
  TEST_EQUAL(references.size(), 2)
  TEST_EQUAL(references[0], "AAN17824")
  TEST_EQUAL(references[1], "GN1736")  
  peptide_hit = peptide_identifications[0].getHits()[1];
  accessions = peptide_hit.extractProteinAccessionsSet();
  references = vector<String>(accessions.begin(), accessions.end());
  TEST_EQUAL(references.size(), 1)
  TEST_EQUAL(references[0], "AAN17824")
  peptide_hit = peptide_identifications[1].getHits()[0];
  accessions = peptide_hit.extractProteinAccessionsSet();
  references = vector<String>(accessions.begin(), accessions.end());
  TEST_EQUAL(references.size(), 1)
  TEST_EQUAL(references[0], "GN1736")

  TEST_EQUAL(peptide_identifications[1].getHits().size(), 1)
  TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[0].getScore(), 33.85)
  TEST_REAL_SIMILAR(peptide_identifications[0].getHits()[1].getScore(), 33.12)
  TEST_REAL_SIMILAR(peptide_identifications[1].getHits()[0].getScore(), 43.9)
  TEST_EQUAL(peptide_identifications[0].getScoreType(), "Mascot")
  TEST_EQUAL(peptide_identifications[1].getScoreType(), "Mascot")
  TEST_EQUAL(protein_identification.getDateTime() == date, true)
  TEST_EQUAL(peptide_identifications[0].getHits()[0].getSequence(), aa_sequence_1)
  TEST_EQUAL(peptide_identifications[0].getHits()[1].getSequence(), aa_sequence_2)
  TEST_EQUAL(peptide_identifications[1].getHits()[0].getSequence(), aa_sequence_3)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
