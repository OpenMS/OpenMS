// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <vector>
using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ProteaseDigestion, "$Id$")

/////////////////////////////////////////////////////////////
ProteaseDigestion* pd_ptr = nullptr;
ProteaseDigestion* pd_null = nullptr;

START_SECTION(([EXTRA] ProteaseDigestion()))
    pd_ptr = new ProteaseDigestion;
    TEST_NOT_EQUAL(pd_ptr, pd_null)
END_SECTION

START_SECTION([EXTRA] ~ProteaseDigestion())
    delete pd_ptr;
END_SECTION

START_SECTION(([EXTRA] ProteaseDigestion(const ProteaseDigestion& rhs)))
    ProteaseDigestion pd;
    pd.setMissedCleavages(1234);
    pd.setEnzyme("no cleavage");
    pd.setSpecificity(EnzymaticDigestion::SPEC_SEMI);

    ProteaseDigestion pd2(pd);

    TEST_EQUAL(pd.getMissedCleavages(), pd2.getMissedCleavages());
    TEST_EQUAL(pd.getEnzymeName(), pd2.getEnzymeName());
    TEST_EQUAL(pd.getSpecificity(), pd2.getSpecificity());
END_SECTION

START_SECTION(([EXTRA] ProteaseDigestion& operator=(const ProteaseDigestion& rhs)))
    ProteaseDigestion pd;
    pd.setMissedCleavages(1234);
    pd.setEnzyme("no cleavage");
    pd.setSpecificity(EnzymaticDigestion::SPEC_SEMI);

    ProteaseDigestion pd2;
    pd2 = pd;

    TEST_EQUAL(pd.getMissedCleavages(), pd2.getMissedCleavages());
    TEST_EQUAL(pd.getEnzymeName(), pd2.getEnzymeName());
    TEST_EQUAL(pd.getSpecificity(), pd2.getSpecificity());
END_SECTION

START_SECTION((void setEnzyme(const String& enzyme_name)))
    ProteaseDigestion pd;
    pd.setEnzyme("Trypsin");
    TEST_EQUAL(pd.getEnzymeName(), "Trypsin");
    pd.setEnzyme("Trypsin/P");
    TEST_EQUAL(pd.getEnzymeName(), "Trypsin/P");
END_SECTION

START_SECTION((Size peptideCount(const AASequence& protein)))
    ProteaseDigestion pd;
    for (int i = 0; i < 2; ++i) // common cases for Trypsin and Trypsin_P
    {
      if (i == 0)
      {
        pd.setEnzyme("Trypsin");
      }
      else if (i == 1)
      {
        pd.setEnzyme("Trypsin/P");
      }
      pd.setMissedCleavages(0);
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACDE")), 1)
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACKDE")), 2)
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACRDE")), 2)
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("ARCRDRE")), 4)
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("RKR")), 3)
      pd.setMissedCleavages(1);
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACDE")), 1)
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACRDE")), 3)
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("ARCDRE")), 5)
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("RKR")), 5)
      pd.setMissedCleavages(3);
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACDE")), 1)
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACRDE")), 3)
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("ARCDRE")), 6)
      TEST_EQUAL(pd.peptideCount(AASequence::fromString("RKR")), 6)
    }
    // special cases:
    pd.setMissedCleavages(0);
    pd.setEnzyme("Trypsin");
    TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACKPDE")), 1)
    TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACRPDE")), 1)
    TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACKPDERA")), 2)
    TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACRPDEKA")), 2)
    pd.setEnzyme("Trypsin/P");
    TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACKPDE")), 2)
    TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACRPDE")), 2)
    TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACKPDERA")), 3)
    TEST_EQUAL(pd.peptideCount(AASequence::fromString("ACRPDEKA")), 3)
END_SECTION

START_SECTION((Size digest(const AASequence& protein, std::vector<AASequence>& output, Size min_length = 1, Size max_length = 0) const))
    ProteaseDigestion pd;
    vector<AASequence> out;

    pd.digest(AASequence::fromString("ACDE"), out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].toString(), "ACDE")

    pd.digest(AASequence::fromString("ACKDE"), out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].toString(), "ACK")
    TEST_EQUAL(out[1].toString(), "DE")

    pd.digest(AASequence::fromString("ACRDE"), out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].toString(), "ACR")
    TEST_EQUAL(out[1].toString(), "DE")

    pd.digest(AASequence::fromString("ACKPDE"), out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].toString(), "ACKPDE")

    pd.digest(AASequence::fromString("ACRPDE"), out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].toString(), "ACRPDE")

    pd.digest(AASequence::fromString("ARCRDRE"), out);
    TEST_EQUAL(out.size(), 4)
    TEST_EQUAL(out[0].toString(), "AR")
    TEST_EQUAL(out[1].toString(), "CR")
    TEST_EQUAL(out[2].toString(), "DR")
    TEST_EQUAL(out[3].toString(), "E")

    pd.digest(AASequence::fromString("RKR"), out);
    TEST_EQUAL(out.size(), 3)
    TEST_EQUAL(out[0].toString(), "R")
    TEST_EQUAL(out[1].toString(), "K")
    TEST_EQUAL(out[2].toString(), "R")

    pd.setMissedCleavages(1);

    pd.digest(AASequence::fromString("ACDE"), out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].toString(), "ACDE")

    pd.digest(AASequence::fromString("ACRDE"), out);
    TEST_EQUAL(out.size(), 3)
    TEST_EQUAL(out[0].toString(), "ACR")
    TEST_EQUAL(out[1].toString(), "DE")
    TEST_EQUAL(out[2].toString(), "ACRDE")

    pd.digest(AASequence::fromString("ARCDRE"), out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].toString(), "AR")
    TEST_EQUAL(out[1].toString(), "CDR")
    TEST_EQUAL(out[2].toString(), "E")
    TEST_EQUAL(out[3].toString(), "ARCDR")
    TEST_EQUAL(out[4].toString(), "CDRE")
    Size discarded = pd.digest(AASequence::fromString("ARCDRE"), out, 3, 4);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].toString(), "CDR")
    TEST_EQUAL(out[1].toString(), "CDRE")
    TEST_EQUAL(discarded, 3)

    pd.digest(AASequence::fromString("RKR"), out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].toString(), "R")
    TEST_EQUAL(out[1].toString(), "K")
    TEST_EQUAL(out[2].toString(), "R")
    TEST_EQUAL(out[3].toString(), "RK")
    TEST_EQUAL(out[4].toString(), "KR")

    pd.digest(AASequence::fromString("(ICPL:2H(4))ARCDRE"), out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].toString(), ".(ICPL:2H(4))AR")
    TEST_EQUAL(out[1].toString(), "CDR")
    TEST_EQUAL(out[2].toString(), "E")
    TEST_EQUAL(out[3].toString(), ".(ICPL:2H(4))ARCDR")
    TEST_EQUAL(out[4].toString(), "CDRE")

    pd.digest(AASequence::fromString("ARCDRE.(Amidated)"), out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].toString(), "AR")
    TEST_EQUAL(out[1].toString(), "CDR")
    TEST_EQUAL(out[2].toString(), "E.(Amidated)")
    TEST_EQUAL(out[3].toString(), "ARCDR")
    TEST_EQUAL(out[4].toString(), "CDRE.(Amidated)")
    discarded = pd.digest(AASequence::fromString("ARCDRE.(Amidated)"), out, 3, 4);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].toString(), "CDR")
    TEST_EQUAL(out[1].toString(), "CDRE.(Amidated)")
    TEST_EQUAL(discarded, 3)

    // ------------------------
    // Trypsin/P
    // ------------------------
    pd.setMissedCleavages(0);
    pd.setEnzyme("Trypsin/P");
    pd.digest(AASequence::fromString("ACKPDE"), out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].toString(), "ACK")
    TEST_EQUAL(out[1].toString(), "PDE")

    pd.digest(AASequence::fromString("ACRPDE"), out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].toString(), "ACR")
    TEST_EQUAL(out[1].toString(), "PDE")

    // ------------------------
    // unspecific cleavage
    // ------------------------
    pd.setEnzyme("unspecific cleavage");
    pd.digest(AASequence::fromString("ABCDEFGHIJ"), out);
    TEST_EQUAL(out.size(), 11*10/2)
    pd.digest(AASequence::fromString("ABC"), out);
    TEST_EQUAL(out.size(), 4*3/2)
END_SECTION

START_SECTION((bool isValidProduct(const String& protein, int pep_pos, int pep_length, bool ignore_missed_cleavages, bool allow_nterm_protein_cleavage, bool allow_random_asp_pro_cleavage)))
    NOT_TESTABLE // tested by overload below
END_SECTION

START_SECTION((bool isValidProduct(const AASequence& protein, int pep_pos, int pep_length, bool ignore_missed_cleavages, bool allow_nterm_protein_cleavage, bool allow_random_asp_pro_cleavage)))
    ProteaseDigestion pd;
    pd.setEnzyme("Trypsin");
    pd.setSpecificity(EnzymaticDigestion::SPEC_FULL); // require both sides

    AASequence prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(pd.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(pd.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(pd.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(pd.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size

    TEST_EQUAL(pd.isValidProduct(prot, 0, 3), false); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 0, 19), false); // invalid C-term - followed by proline
    TEST_EQUAL(pd.isValidProduct(prot, 8, 3), false); // invalid C-term
    TEST_EQUAL(pd.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7), false); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    prot = AASequence::fromString("MBCDEFGKABCRAAAKAA"); // starts with Met - we assume the cleaved form without Met occurs in vivo
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7, true, true), true); // valid N-term (since protein starts with Met)
    TEST_EQUAL(pd.isValidProduct(prot, 2, 6, true, true), true); // valid N-term (since protein starts with Met and second AA maybe cleaved in XTandem)
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7, true, false), false); // invalid N-term (since Met cleavage is not allowed.)
    TEST_EQUAL(pd.isValidProduct(prot, 2, 6, true, false), false); // invalid N-term (since Met cleavage is not allowed.)
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    //################################################
    // same as above, just with other specificity

    pd.setSpecificity(EnzymaticDigestion::SPEC_SEMI); // require one special cleavage site
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(pd.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(pd.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(pd.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(pd.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size

    TEST_EQUAL(pd.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 0, 19), true); // invalid C-term - followed by proline
    TEST_EQUAL(pd.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(pd.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7), true); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    prot = AASequence::fromString("MBCDEFGKABCRAAAKAA"); // starts with Met - we assume the cleaved form without Met occurs in vivo
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7, true, true), true); // valid N-term (since protein starts with Met)
    TEST_EQUAL(pd.isValidProduct(prot, 2, 6, true, true), true); // valid N-term (since protein starts with Met and second AA maybe cleaved in XTandem)
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7, true, false), true); // invalid N-term (since Met cleavage is not allowed.)
    TEST_EQUAL(pd.isValidProduct(prot, 2, 6, true, false), true); // invalid N-term (since Met cleavage is not allowed.)
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    //################################################
    // same as above, just with other specificity

    pd.setSpecificity(EnzymaticDigestion::SPEC_NONE); // require no special cleavage site
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(pd.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(pd.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(pd.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(pd.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size

    TEST_EQUAL(pd.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 0, 19), true); // invalid C-term - followed by proline
    TEST_EQUAL(pd.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(pd.isValidProduct(prot, 3, 6), true); // invalid C+N-term
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7), true); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    // ------------------------
    // Trypsin/P
    // ------------------------
    pd.setEnzyme("Trypsin/P");
    pd.setSpecificity(EnzymaticDigestion::SPEC_FULL); // require both sides

    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(pd.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(pd.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(pd.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(pd.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size

    TEST_EQUAL(pd.isValidProduct(prot, 0, 3), false); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 0, 19), true); // valid C-term - followed by proline
    TEST_EQUAL(pd.isValidProduct(prot, 8, 3), false); // invalid C-term
    TEST_EQUAL(pd.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7), false); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    prot = AASequence::fromString("MBCDEFGKABCRAAAKAA"); // starts with Met - we assume the cleaved form without Met occurs in vivo
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7, true, true), true); // valid N-term (since protein starts with Met)
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7), false); // invalid N-term (since Met cleavage is not allowpd.)
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    // test with different missed cleavages when this is not ignored (ignore_missed_cleavages = false)
    //                                    |8  |12 |16|19
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB"); // 4 cleavages at {(0),8,12,16,19}
    pd.setMissedCleavages(0); // redundant, by default zero, should be zero
    TEST_EQUAL(pd.isValidProduct(prot, 8, 4, false, false), true);  //  valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 8, 8, false, false), false); //  invalid, fully-tryptic but with a missing cleavage
    pd.setMissedCleavages(1);
    TEST_EQUAL(pd.isValidProduct(prot, 8, 8, false, false), true);  //  valid, fully-tryptic with 1 missing cleavage (allow)
    TEST_EQUAL(pd.isValidProduct(prot, 8, 11, false, false), false);//  invalid, fully-tryptic but with 2 missing cleavages
    pd.setMissedCleavages(2);
    TEST_EQUAL(pd.isValidProduct(prot, 8, 11, false, false), true); //  valid, fully-tryptic with 2 missing cleavages
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, true, false), true);  //  boundary case, length of protein (no checking of MCs)
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), false);//  boundary case, this exceeds missing cleavages
    TEST_EQUAL(pd.isValidProduct(prot, 0, 19, false, false), false);//  start-boundary case, 2 allowed, 3 required
    pd.setMissedCleavages(3);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), false);//  boundary case, invalid: 3 allowed, 4 required
    TEST_EQUAL(pd.isValidProduct(prot, 0, 19, false, false), true); //  start-boundary case, 3 allowed, 3 required
    pd.setMissedCleavages(4); // maximum cleavages for this peptide
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), true); //  boundary case, accepted: 4 allowed, 4 required
    TEST_EQUAL(pd.isValidProduct(prot, 0, 19, false, false), true); //  start-boundary case, 4 allowed, 3 required
    pd.setMissedCleavages(5); // allow even more ...
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), true); //  boundary case, accepted: 5 allowed, 4 required
    pd.setMissedCleavages(0); // set back to default

    //################################################
    // same as above, just with other specificity

    pd.setSpecificity(EnzymaticDigestion::SPEC_SEMI); // require one special cleavage site
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(pd.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(pd.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(pd.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(pd.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size

    TEST_EQUAL(pd.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 0, 19), true); // valid C-term - followed by proline
    TEST_EQUAL(pd.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(pd.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7), true); // invalid N-term valid C-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    prot = AASequence::fromString("MBCDEFGKABCRAAAKAA"); // starts with Met - we assume the cleaved form without Met occurs in vivo
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7, true, true), true); // valid N-term (since protein starts with Met)
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7, true, false), true); // invalid N-term (since Met cleavage is not allowpd.)
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    // test with different missed cleavages when this is not ignored (ignore_missed_cleavages = false)
    //                                    |8  |12 |16|19
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB"); // 4 cleavages at {(0),8,12,16,19}
    pd.setMissedCleavages(0); // redundant, by default zero, should be zero
    TEST_EQUAL(pd.isValidProduct(prot, 8, 3, false, false), true);  //  valid semi-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 8, 5, false, false), false); //  invalid, semi-tryptic but with a missing cleavage
    pd.setMissedCleavages(1);
    TEST_EQUAL(pd.isValidProduct(prot, 8, 5, false, false), true);  //  valid, semi-tryptic with 1 missing cleavage (allow)
    TEST_EQUAL(pd.isValidProduct(prot, 8, 10, false, false), false);//  invalid, semi-tryptic but with 2 missing cleavages
    pd.setMissedCleavages(2);
    TEST_EQUAL(pd.isValidProduct(prot, 8, 10, false, false), true); //  valid, semi-tryptic with 2 missing cleavages
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, true, false), true);  //  boundary case, length of protein (no checking of MCs)
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), false);//  boundary case, this exceeds missing cleavages
    TEST_EQUAL(pd.isValidProduct(prot, 0, 18, false, false), false);//  start-boundary case, 2 allowed, 3 required
    pd.setMissedCleavages(3);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), false);//  boundary case, invalid: 3 allowed, 4 required
    TEST_EQUAL(pd.isValidProduct(prot, 0, 18, false, false), true); //  start-boundary case, 3 allowed, 3 required
    pd.setMissedCleavages(4); // maximum cleavages for this peptide
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), true); //  boundary case, accepted: 4 allowed, 4 required
    TEST_EQUAL(pd.isValidProduct(prot, 0, 18, false, false), true); //  start-boundary case, 4 allowed, 3 required
    pd.setMissedCleavages(5); // allow even more ...
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), true); //  boundary case, accepted: 5 allowed, 4 required
    pd.setMissedCleavages(0); // set back to default

    //################################################
    // same as above, just with other specificity

    pd.setSpecificity(EnzymaticDigestion::SPEC_NONE); // require no special cleavage site
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(pd.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(pd.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(pd.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(pd.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size

    TEST_EQUAL(pd.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 0, 19), true); // valid C-term - followed by proline
    TEST_EQUAL(pd.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(pd.isValidProduct(prot, 3, 6), true); // invalid C+N-term
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7), true); // invalid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    // test with different missed cleavages when this is not ignored (ignore_missed_cleavages = false)
    //                                    |8  |12 |16|19
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB"); // 4 cleavages at {(0),8,12,16,19}
    pd.setMissedCleavages(0); // redundant, by default zero, should be zero
    TEST_EQUAL(pd.isValidProduct(prot, 9, 2, false, false), true);  //  valid not-tryptic
    TEST_EQUAL(pd.isValidProduct(prot, 9, 5, false, false), false); //  invalid, not-tryptic but with a missing cleavage
    pd.setMissedCleavages(1);
    TEST_EQUAL(pd.isValidProduct(prot, 9, 5, false, false), true);  //  valid, not-tryptic with 1 missing cleavage (allow)
    TEST_EQUAL(pd.isValidProduct(prot, 9, 9, false, false), false); //  invalid, semi-tryptic but with 2 missing cleavages
    pd.setMissedCleavages(2);
    TEST_EQUAL(pd.isValidProduct(prot, 9, 9, false, false), true);  //  valid, semi-tryptic with 2 missing cleavages
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, true, false), true);  //  boundary case, length of protein (no checking of MCs)
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), false);//  boundary case, this exceeds missing cleavages
    pd.setMissedCleavages(3);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), false);//  boundary case, invalid: 3 allowed, 4 required
    pd.setMissedCleavages(4); // maximum cleavages for this peptide
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), true); //  boundary case, accepted: 4 allowed, 4 required
    pd.setMissedCleavages(5); // allow even more ...
    TEST_EQUAL(pd.isValidProduct(prot, 0, 24, false, false), true); //  boundary case, accepted: 5 allowed, 4 required
    pd.setMissedCleavages(0); // set back to default

    // tests with random Asp-Pro Cleavages
    pd.setSpecificity(EnzymaticDigestion::SPEC_SEMI);
    //                                  |6*  |11|14*|18 |22|25 |29* |34*
    prot = AASequence::fromString("MCABCDPEFGKACDPBCRAAAKAARPBBDPBBCDP"); // 4 real cleavages at {(0),11,18,22,25}

    pd.setMissedCleavages(0); // redundant, by default zero, should be zero
    TEST_EQUAL(pd.isValidProduct(prot, 0, 2), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 1, 2, true, true), true); // valid N-term, M cleavage
    TEST_EQUAL(pd.isValidProduct(prot, 2, 2, true, true), true); // valid N-term, MC cleavage
    TEST_EQUAL(pd.isValidProduct(prot, 3, 2, true, true), false); // invalid N- and C-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 3, true, true, true), true); // valid N-term (random D|P cleavage)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 3, true, true, false), false); // invalid C- and N-term (random D|P cleavage not allowed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 5, true, true, true), true); // valid fully-tryptic (N-term random D|P)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 6, true, true, true), true); // valid N-term (random D|P)
    TEST_EQUAL(pd.isValidProduct(prot, 29, 4, true, true, true), true); // valid fully-tryptic (N-term random D|P)
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // the whole thing
    TEST_EQUAL(pd.isValidProduct(prot, 1, prot.size()-1, true), true); // valid C-term (end), M cleavage not allowed
    TEST_EQUAL(pd.isValidProduct(prot, 2, prot.size()-2, true), true); // valid C-term (end), MC cleavage not allowed
    TEST_EQUAL(pd.isValidProduct(prot, 1, prot.size()-1, true, true), true); // valid C- and N-term (end, start M cleavage)
    TEST_EQUAL(pd.isValidProduct(prot, 2, prot.size()-2, true, true), true); // valid C- and N-term (end, start MC cleavage)


    TEST_EQUAL(pd.isValidProduct(prot, 6, 6, false, false, true), false);  //  valid N-term (D|P), one real cleavage site skipped (11)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, false, false, true), false);  //  valid termini (both D|P), one real cleavage site skipped (11)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 10, false, false, true), false);  //  valid N-term (D|P), skipped real cleavage (11) and D|P (14)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 12, false, false, true), false);  //  valid termini (D|P, real), skipped real cleavage (11) and D|P (14)
    TEST_EQUAL(pd.isValidProduct(prot, 11, 4, false, false, true), true);  //  valid N-term (real), one D|P cleavage site skipped (14)
    TEST_EQUAL(pd.isValidProduct(prot, 11, 7, false, false, true), true);  //  valid termini (both real), one D|P cleavage site skipped (14)
    TEST_EQUAL(pd.isValidProduct(prot, 18, 16, false, false, true), false);  //  valid termini (real, D|P), two real cleavage sites skipped (22,25), one D|P (29)
    TEST_EQUAL(pd.isValidProduct(prot, 18, 17, false, false, true), false);  //  valid termini (real, end), two real cleavage sites skipped (22,25), two D|P (35)

    pd.setMissedCleavages(1);

    TEST_EQUAL(pd.isValidProduct(prot, 6, 6, false, false, true), true);  //  valid N-term (D|P), one real cleavage site skipped (11)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, false, false, true), true);  //  valid termini (both D|P), one real cleavage site skipped (11)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 10, false, false, true), true);  //  valid N-term (D|P), skipped real cleavage (11) and D|P (14)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 12, false, false, true), true);  //  valid termini (D|P, real), skipped real cleavage (11) and D|P (14)
    TEST_EQUAL(pd.isValidProduct(prot, 11, 4, false, false, true), true);  //  valid N-term (real), one D|P cleavage site skipped (14)
    TEST_EQUAL(pd.isValidProduct(prot, 11, 7, false, false, true), true);  //  valid termini (both real), one D|P cleavage site skipped (14)
    TEST_EQUAL(pd.isValidProduct(prot, 18, 16, false, false, true), false);  //  valid termini (real, D|P), two real cleavage sites skipped (22,25), one D|P (29)
    TEST_EQUAL(pd.isValidProduct(prot, 18, 17, false, false, true), false);  //  valid termini (real, end), two real cleavage sites skipped (22,25), two D|P (35)
    
    pd.setMissedCleavages(2);

    TEST_EQUAL(pd.isValidProduct(prot, 6, 6, false, false, true), true);  //  valid N-term (D|P), one real cleavage site skipped (11)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, false, false, true), true);  //  valid termini (both D|P), one real cleavage site skipped (11)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 10, false, false, true), true);  //  valid N-term (D|P), skipped real cleavage (11) and D|P (14)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 12, false, false, true), true);  //  valid termini (D|P, real), skipped real cleavage (11) and D|P (14)
    TEST_EQUAL(pd.isValidProduct(prot, 11, 4, false, false, true), true);  //  valid N-term (real), one D|P cleavage site skipped (14)
    TEST_EQUAL(pd.isValidProduct(prot, 11, 7, false, false, true), true);  //  valid termini (both real), one D|P cleavage site skipped (14)
    TEST_EQUAL(pd.isValidProduct(prot, 18, 16, false, false, true), true);  //  valid termini (real, D|P), two real cleavage sites skipped (22,25), one D|P (29)
    TEST_EQUAL(pd.isValidProduct(prot, 18, 17, false, false, true), true);  //  valid termini (real, end), two real cleavage sites skipped (22,25), two D|P (35)
   
    // more than needed
    pd.setMissedCleavages(3);

    TEST_EQUAL(pd.isValidProduct(prot, 6, 6, false, false, true), true);  //  valid N-term (D|P), one real cleavage site skipped (11)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, false, false, true), true);  //  valid termini (both D|P), one real cleavage site skipped (11)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 10, false, false, true), true);  //  valid N-term (D|P), skipped real cleavage (11) and D|P (14)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 12, false, false, true), true);  //  valid termini (D|P, real), skipped real cleavage (11) and D|P (14)
    TEST_EQUAL(pd.isValidProduct(prot, 11, 4, false, false, true), true);  //  valid N-term (real), one D|P cleavage site skipped (14)
    TEST_EQUAL(pd.isValidProduct(prot, 11, 7, false, false, true), true);  //  valid termini (both real), one D|P cleavage site skipped (14)
    TEST_EQUAL(pd.isValidProduct(prot, 18, 16, false, false, true), true);  //  valid termini (real, D|P), two real cleavage sites skipped (22,25), one D|P (29)
    TEST_EQUAL(pd.isValidProduct(prot, 18, 17, false, false, true), true);  //  valid termini (real, end), two real cleavage sites skipped (22,25), two D|P (35)

    // ------------------------
    // glutamyl endopeptidase (since it cleaves after D as the random XTandem cleavage)
    // ------------------------
    pd.setEnzyme("glutamyl endopeptidase");
    pd.setSpecificity(EnzymaticDigestion::SPEC_SEMI);
    //                                  |6*  |11|14*|18 |22|25 |29* |34*
    prot = AASequence::fromString("MCABCDPLFGKACDPHCRAAAKAARPHHDPHHCDP"); // 4 real cleavages at {(0),11,18,22,25}

    pd.setMissedCleavages(0); // redundant, by default zero, should be zero
    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, true, true, true), true); // valid C and N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, true, true, false), true); // valid C and N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 5, true, true, true), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 5, true, true, false), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 23, true, true, true), true); // valid C and N-term (but 1 missed, ignored)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 23, true, true, false), true); // valid C and N-term (but 1 missed, ignored)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 28, true, true, true), true); // valid C and N-term (but 2 missed, ignored)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 28, true, true, false), true); // valid C and N-term (but 2 missed, ignored)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 23, false, true, true), false); // valid C and N-term (but 1 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 23, false, true, false), false); // valid C and N-term (but 1 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 28, false, true, true), false); // valid C and N-term (but 2 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 28, false, true, false), false); // valid C and N-term (but 2 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size()), true); // valid N-term start, C-term end (3 missed, ignored)
    TEST_EQUAL(pd.isValidProduct(prot, 1, prot.size()-1, true), true); // valid C-term end (3 missed, ignored)
    TEST_EQUAL(pd.isValidProduct(prot, 2, prot.size()-2, true), true); // valid C-term end (3 missed, ignored)
    TEST_EQUAL(pd.isValidProduct(prot, 0, prot.size(), false), false); // valid N-term start, C-term end (but 3 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 1, prot.size()-1, false), false); // valid C-term end (but 3 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 2, prot.size()-2, false), false); // valid C-term end (but 3 missed)

    pd.setMissedCleavages(1);

    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, false, true, true), true); // valid C and N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, false, true, false), true); // valid C and N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 5, false, true, true), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 5, false, true, false), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 23, false, true, true), true); // valid C and N-term (but 1 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 23, false, true, false), true); // valid C and N-term (but 1 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 28, false, true, true), false); // valid C and N-term (but 2 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 28, false, true, false), false); // valid C and N-term (but 2 missed)

    pd.setMissedCleavages(2);

    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, false, true, true), true); // valid C and N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, false, true, false), true); // valid C and N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 5, false, true, true), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 5, false, true, false), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 23, false, true, true), true); // valid C and N-term (but 1 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 23, false, true, false), true); // valid C and N-term (but 1 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 28, false, true, true), true); // valid C and N-term (but 2 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 28, false, true, false), true); // valid C and N-term (but 2 missed)

    // more than needed
    pd.setMissedCleavages(3);

    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, false, true, true), true); // valid C and N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 8, false, true, false), true); // valid C and N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 5, false, true, true), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 5, false, true, false), true); // valid N-term
    TEST_EQUAL(pd.isValidProduct(prot, 6, 23, false, true, true), true); // valid C and N-term (but 1 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 23, false, true, false), true); // valid C and N-term (but 1 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 28, false, true, true), true); // valid C and N-term (but 2 missed)
    TEST_EQUAL(pd.isValidProduct(prot, 6, 28, false, true, false), true); // valid C and N-term (but 2 missed)

    /// no cleavage:
    pd.setEnzyme("no cleavage");
    pd.setSpecificity(EnzymaticDigestion::SPEC_FULL); // require both sides

    prot = AASequence::fromString("MADPDE"); // something which start with 'M' to test nterm prot cleavage
    TEST_EQUAL(pd.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(pd.isValidProduct(prot, 1, 300), false); // invalid length
    TEST_EQUAL(pd.isValidProduct(prot, 3, 0), false);   // invalid size

    const bool with_nterm_prot_cleavage = true;
    const bool no_nterm_prot_cleavage = false;
    const bool with_rnd_asppro_cleavage = true;
    const bool no_rnd_asppro_cleavage = false;
    TEST_EQUAL(pd.isValidProduct(prot, 0, 6, true, with_nterm_prot_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 6, true, no_nterm_prot_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 1, 5, true, with_nterm_prot_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 1, 5, true, no_nterm_prot_cleavage), false);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 5, true, with_nterm_prot_cleavage), false);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 5, true, no_nterm_prot_cleavage), false);

    // asppro_cleavage should make no difference
    TEST_EQUAL(pd.isValidProduct(prot, 0, 6, true, with_nterm_prot_cleavage, with_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 6, true, no_nterm_prot_cleavage, no_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 1, 5, true, with_nterm_prot_cleavage, with_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 1, 5, true, no_nterm_prot_cleavage, no_rnd_asppro_cleavage), false);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 5, true, with_nterm_prot_cleavage, with_rnd_asppro_cleavage), false);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 5, true, no_nterm_prot_cleavage, no_rnd_asppro_cleavage), false);
    // ... unless we hit one
    TEST_EQUAL(pd.isValidProduct(prot, 1, 2, true, no_nterm_prot_cleavage, no_rnd_asppro_cleavage), false);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 3, true, with_nterm_prot_cleavage, with_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 3, true, no_nterm_prot_cleavage, no_rnd_asppro_cleavage), false);

    pd.setSpecificity(EnzymaticDigestion::SPEC_SEMI); // require one side

    prot = AASequence::fromString("MADPDE");            // something which starts with 'M' to test nterm prot cleavage and has a D/P cleavage
    TEST_EQUAL(pd.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(pd.isValidProduct(prot, 1, 300), false); // invalid length
    TEST_EQUAL(pd.isValidProduct(prot, 3, 0), false);   // invalid size

    TEST_EQUAL(pd.isValidProduct(prot, 0, 6, true, with_nterm_prot_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 6, true, no_nterm_prot_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 1, 5, true, with_nterm_prot_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 1, 5, true, no_nterm_prot_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 5, true, with_nterm_prot_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 5, true, no_nterm_prot_cleavage), true);

    // asppro_cleavage should make no difference
    TEST_EQUAL(pd.isValidProduct(prot, 0, 6, true, with_nterm_prot_cleavage, with_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 6, true, no_nterm_prot_cleavage, no_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 1, 5, true, with_nterm_prot_cleavage, with_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 1, 5, true, no_nterm_prot_cleavage, no_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 5, true, with_nterm_prot_cleavage, with_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 5, true, no_nterm_prot_cleavage, no_rnd_asppro_cleavage), true);
    // ... unless we hit one
    TEST_EQUAL(pd.isValidProduct(prot, 1, 2, true, no_nterm_prot_cleavage, no_rnd_asppro_cleavage), false);
    TEST_EQUAL(pd.isValidProduct(prot, 1, 2, true, no_nterm_prot_cleavage, with_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 3, true, with_nterm_prot_cleavage, with_rnd_asppro_cleavage), true);
    TEST_EQUAL(pd.isValidProduct(prot, 0, 3, true, no_nterm_prot_cleavage, no_rnd_asppro_cleavage), true);

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
