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

#include <OpenMS/FORMAT/UniProtXMLFile.h>

#include <algorithm>
#include <vector>

///////////////////////////

START_TEST(UniProtXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

UniProtXMLFile* ptr = nullptr;
UniProtXMLFile* null_pointer = nullptr;

START_SECTION((UniProtXMLFile()))
  ptr = new UniProtXMLFile();
  TEST_NOT_EQUAL(ptr, null_pointer)
END_SECTION

START_SECTION((~UniProtXMLFile()))
  delete ptr;
END_SECTION

START_SECTION((void load(const std::string& filename, std::vector<UniProtEntry>& entries)))
  UniProtXMLFile xml;
  vector<UniProtEntry> entries;
  xml.load(OPENMS_GET_TEST_DATA_PATH("UniProtXMLFile_test.xml"), entries);

  // Two entries: one Swiss-Prot kitchen-sink, one minimal TrEMBL.
  TEST_EQUAL(entries.size(), 2)

  // ── Entry 1: Swiss-Prot ──────────────────────────────────────────
  const UniProtEntry& e1 = entries[0];
  TEST_EQUAL(e1.dataset, "Swiss-Prot")
  TEST_EQUAL(e1.entry_version, "42")
  TEST_EQUAL(e1.accession, "P00001")
  TEST_EQUAL(e1.alt_accessions.size(), 1)
  TEST_EQUAL(e1.alt_accessions[0], "P99999")
  TEST_EQUAL(e1.name, "TEST1_HUMAN")
  TEST_EQUAL(e1.full_name, "Test protein one")
  // The primary gene name must be captured even though a synonym appears first.
  TEST_EQUAL(e1.primary_gene, "PRIMGENE")
  TEST_EQUAL(e1.ncbi_tax_id, "9606")
  // The scientific organism name must be captured even though a common name appears first.
  TEST_EQUAL(e1.tax_name, "Homo sapiens")
  TEST_EQUAL(e1.protein_existence, "evidence at protein level")
  TEST_EQUAL(e1.sequence_version, "2")
  TEST_EQUAL(e1.sequence, "MSEQNCEMTFQIQRIYNCKDISFEAPNAPC")
  // Isoform sequences (inside <alternative products>) must be skipped, so
  // the canonical 30-residue sequence is the one we get.
  TEST_EQUAL(e1.sequence.size(), 30)

  // Four features: signal peptide, modified residue, disulfide bond, sequence variant.
  TEST_EQUAL(e1.features.size(), 4)

  const UniProtFeature& f_sig = e1.features[0];
  TEST_EQUAL(f_sig.type, "signal peptide")
  TEST_EQUAL(f_sig.has_range, true)
  TEST_EQUAL(f_sig.begin, 1)
  TEST_EQUAL(f_sig.end,   10)

  const UniProtFeature& f_mod = e1.features[1];
  TEST_EQUAL(f_mod.type, "modified residue")
  TEST_EQUAL(f_mod.description, "Phosphoserine")
  TEST_EQUAL(f_mod.has_position, true)
  TEST_EQUAL(f_mod.position, 20)

  const UniProtFeature& f_ss = e1.features[2];
  TEST_EQUAL(f_ss.type, "disulfide bond")
  TEST_EQUAL(f_ss.has_range, true)
  TEST_EQUAL(f_ss.begin, 15)
  TEST_EQUAL(f_ss.end,   25)

  const UniProtFeature& f_var = e1.features[3];
  TEST_EQUAL(f_var.type, "sequence variant")
  TEST_EQUAL(f_var.original, "A")
  TEST_EQUAL(f_var.variation, "V")
  TEST_EQUAL(f_var.has_position, true)
  TEST_EQUAL(f_var.position, 30)

  // ── Entry 2: TrEMBL ──────────────────────────────────────────────
  const UniProtEntry& e2 = entries[1];
  TEST_EQUAL(e2.dataset, "TrEMBL")
  TEST_EQUAL(e2.accession, "Q00002")
  TEST_EQUAL(e2.name, "TEST2_HUMAN")
  TEST_EQUAL(e2.sequence, "MSEQE")
  TEST_EQUAL(e2.features.size(), 0)
END_SECTION

START_SECTION((void loadStreaming(const std::string& filename, const std::function<void(UniProtEntry&&)>& callback)))
  UniProtXMLFile xml;
  size_t entry_count = 0;
  vector<string> accessions;
  xml.loadStreaming(OPENMS_GET_TEST_DATA_PATH("UniProtXMLFile_test.xml"),
                    [&](UniProtEntry&& e) {
                      ++entry_count;
                      accessions.push_back(e.accession);
                    });
  TEST_EQUAL(entry_count, 2)
  TEST_EQUAL(accessions.size(), 2)
  TEST_EQUAL(accessions[0], "P00001")
  TEST_EQUAL(accessions[1], "Q00002")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
