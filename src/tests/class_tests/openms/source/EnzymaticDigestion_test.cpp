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

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <vector>
using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(EnzymaticDigestion, "$Id$")

/////////////////////////////////////////////////////////////
EnzymaticDigestion* ed_ptr = nullptr;
EnzymaticDigestion* ed_null = nullptr;

START_SECTION((EnzymaticDigestion()))
    ed_ptr = new EnzymaticDigestion;
    TEST_NOT_EQUAL(ed_ptr, ed_null)
END_SECTION

START_SECTION([EXTRA] virtual ~EnzymaticDigestion())
  delete ed_ptr;
  NOT_TESTABLE
END_SECTION

START_SECTION(([EXTRA] EnzymaticDigestion(const EnzymaticDigestion& rhs)))
    EnzymaticDigestion ed;
    ed.setMissedCleavages(1234);
    ed.setEnzyme(ProteaseDB::getInstance()->getEnzyme("no cleavage"));
    ed.setSpecificity(EnzymaticDigestion::SPEC_SEMI);

    EnzymaticDigestion ed2(ed);

    TEST_EQUAL(ed.getMissedCleavages(), ed2.getMissedCleavages());
    TEST_EQUAL(ed.getEnzymeName(), ed2.getEnzymeName());
    TEST_EQUAL(ed.getSpecificity(), ed2.getSpecificity());
END_SECTION

START_SECTION(([EXTRA] EnzymaticDigestion(const EnzymaticDigestion& rhs)))
    EnzymaticDigestion ed;
    ed.setMissedCleavages(1234);
    ed.setEnzyme(ProteaseDB::getInstance()->getEnzyme("no cleavage"));
    ed.setSpecificity(EnzymaticDigestion::SPEC_SEMI);

    EnzymaticDigestion ed2(ed);

    TEST_EQUAL(ed.getMissedCleavages(), ed2.getMissedCleavages());
    TEST_EQUAL(ed.getEnzymeName(), ed2.getEnzymeName());
    TEST_EQUAL(ed.getSpecificity(), ed2.getSpecificity());
END_SECTION

START_SECTION((Size getMissedCleavages() const))
    TEST_EQUAL(EnzymaticDigestion().getMissedCleavages(), 0)
END_SECTION

START_SECTION((String getEnzymeName() const))
    TEST_EQUAL(EnzymaticDigestion().getEnzymeName(), "Trypsin")
END_SECTION

START_SECTION((void setMissedCleavages(Size missed_cleavages)))
    EnzymaticDigestion ed;
    ed.setMissedCleavages(5);
    TEST_EQUAL(ed.getMissedCleavages(), 5)
END_SECTION

START_SECTION((void setEnzyme(const DigestionEnzyme* enzyme)))
    EnzymaticDigestion ed;
    ed.setEnzyme(ProteaseDB::getInstance()->getEnzyme("Trypsin/P"));
    TEST_EQUAL(ed.getEnzymeName(), "Trypsin/P");
END_SECTION

START_SECTION((Specificity getSpecificity() const))
    EnzymaticDigestion ed;
    TEST_EQUAL(ed.getSpecificity(), EnzymaticDigestion::SPEC_FULL);
    ed.setSpecificity(EnzymaticDigestion::SPEC_NONE);
    TEST_EQUAL(ed.getSpecificity(), EnzymaticDigestion::SPEC_NONE);
    ed.setSpecificity(EnzymaticDigestion::SPEC_SEMI);
    TEST_EQUAL(ed.getSpecificity(), EnzymaticDigestion::SPEC_SEMI);
END_SECTION

START_SECTION((void setSpecificity(Specificity spec)))
    NOT_TESTABLE // tested above
END_SECTION

START_SECTION((static Specificity getSpecificityByName(const String& name)))
    TEST_EQUAL(EnzymaticDigestion::getSpecificityByName(EnzymaticDigestion::NamesOfSpecificity[2]), EnzymaticDigestion::SPEC_FULL);
    TEST_EQUAL(EnzymaticDigestion::getSpecificityByName(EnzymaticDigestion::NamesOfSpecificity[1]), EnzymaticDigestion::SPEC_SEMI);
    TEST_EQUAL(EnzymaticDigestion::getSpecificityByName(EnzymaticDigestion::NamesOfSpecificity[0]), EnzymaticDigestion::SPEC_NONE);
    TEST_EQUAL(EnzymaticDigestion::getSpecificityByName("DoesNotExist"), EnzymaticDigestion::SPEC_UNKNOWN);
END_SECTION

START_SECTION((Size digestUnmodified(const StringView sequence, std::vector<StringView>& output, Size min_length, Size max_length)))
{
    EnzymaticDigestion ed;
    vector<StringView> out;

    // end without cutting site
    std::string s = "ACDE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].getString(), s)

    // end with cutting site
    s = "ACDEK";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].getString(), "ACDEK")

    s = "ACKDE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].getString(), "ACK")
    TEST_EQUAL(out[1].getString(), "DE")

    s = "ACRDE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].getString(), "ACR")
    TEST_EQUAL(out[1].getString(), "DE")

    s = "ACKPDE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].getString(), "ACKPDE")

    s = "ACRPDE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].getString(), "ACRPDE")

    s = "ARCRDRE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 4)
    TEST_EQUAL(out[0].getString(), "AR")
    TEST_EQUAL(out[1].getString(), "CR")
    TEST_EQUAL(out[2].getString(), "DR")
    TEST_EQUAL(out[3].getString(), "E")

    s = "RKR";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 3)
    TEST_EQUAL(out[0].getString(), "R")
    TEST_EQUAL(out[1].getString(), "K")
    TEST_EQUAL(out[2].getString(), "R")

    ed.setMissedCleavages(1);

    s = "ACDE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].getString(), "ACDE")

    s = "ACRDE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 3)
    TEST_EQUAL(out[0].getString(), "ACR")
    TEST_EQUAL(out[1].getString(), "DE")
    TEST_EQUAL(out[2].getString(), "ACRDE")

    s = "ARCDRE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].getString(), "AR")
    TEST_EQUAL(out[1].getString(), "CDR")
    TEST_EQUAL(out[2].getString(), "E")
    TEST_EQUAL(out[3].getString(), "ARCDR")
    TEST_EQUAL(out[4].getString(), "CDRE")

    s = "ARCDRER";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].getString(), "AR")
    TEST_EQUAL(out[1].getString(), "CDR")
    TEST_EQUAL(out[2].getString(), "ER")
    TEST_EQUAL(out[3].getString(), "ARCDR")
    TEST_EQUAL(out[4].getString(), "CDRER")

    s = "RKR";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].getString(), "R")
    TEST_EQUAL(out[1].getString(), "K")
    TEST_EQUAL(out[2].getString(), "R")
    TEST_EQUAL(out[3].getString(), "RK")
    TEST_EQUAL(out[4].getString(), "KR")

    s = "(ICPL:2H(4))ARCDRE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].getString(), "(ICPL:2H(4))AR")
    TEST_EQUAL(out[1].getString(), "CDR")
    TEST_EQUAL(out[2].getString(), "E")
    TEST_EQUAL(out[3].getString(), "(ICPL:2H(4))ARCDR")
    TEST_EQUAL(out[4].getString(), "CDRE")

    s = "ARCDRE(Amidated)";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].getString(), "AR")
    TEST_EQUAL(out[1].getString(), "CDR")
    TEST_EQUAL(out[2].getString(), "E(Amidated)")
    TEST_EQUAL(out[3].getString(), "ARCDR")
    TEST_EQUAL(out[4].getString(), "CDRE(Amidated)")

    ed.setMissedCleavages(2);
    s = "RKR";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 6)
    TEST_EQUAL(out[0].getString(), "R")
    TEST_EQUAL(out[1].getString(), "K")
    TEST_EQUAL(out[2].getString(), "R")
    TEST_EQUAL(out[3].getString(), "RK")
    TEST_EQUAL(out[4].getString(), "KR")
    TEST_EQUAL(out[5].getString(), "RKR")

    // min size
    ed.digestUnmodified(s, out, 2);
    TEST_EQUAL(out.size(), 3)
    TEST_EQUAL(out[0].getString(), "RK")
    TEST_EQUAL(out[1].getString(), "KR")
    TEST_EQUAL(out[2].getString(), "RKR")

    ed.digestUnmodified(s, out, 3);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].getString(), "RKR")

    // max size
    ed.digestUnmodified(s, out, 2,2);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].getString(), "RK")
    TEST_EQUAL(out[1].getString(), "KR")

    // ------------------------
    // Trypsin/P
    // ------------------------
    ed.setMissedCleavages(0);
    ed.setEnzyme(ProteaseDB::getInstance()->getEnzyme("Trypsin/P"));
    s = "ACKPDE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].getString(), "ACK")
    TEST_EQUAL(out[1].getString(), "PDE")

    s = "ACRPDE";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].getString(), "ACR")
    TEST_EQUAL(out[1].getString(), "PDE")

    // ------------------------
    // unspecific cleavage
    // ------------------------
    s = "ABCDEFGHIJ";
    ed.setEnzyme(ProteaseDB::getInstance()->getEnzyme("unspecific cleavage"));
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 11 * 10 / 2)

    // digest with min/max length
    ed.digestUnmodified(s, out, 5, 6);
    for (auto & a : out)
    {
      TEST_EQUAL(a.getString().size() == 5 
      || a.getString().size() == 6, true)    
    }

    s = "ABC";
    ed.digestUnmodified(s, out);
    TEST_EQUAL(out.size(), 4 * 3 / 2);
}
END_SECTION

START_SECTION((bool isValidProduct(const String& sequence, int pos, int length, bool ignore_missed_cleavages)))
{
    EnzymaticDigestion ed;
    ed.setEnzyme(ProteaseDB::getInstance()->getEnzyme("Trypsin"));
    ed.setSpecificity(EnzymaticDigestion::SPEC_FULL); // require both sides

    String prot = "ABCDEFGKABCRAAAKAARPBBBB";
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct("", 10, 0), false); // invalid size

    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), false); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), false); // invalid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), false); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), false); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    //################################################
    // same as above, just with other specificity

    ed.setSpecificity(EnzymaticDigestion::SPEC_SEMI); // require one special cleavage site
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct("", 10, 0), false); // invalid size

    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), true); // invalid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    //################################################
    // same as above, just with other specificity

    ed.setSpecificity(EnzymaticDigestion::SPEC_NONE); // require no special cleavage site
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct("", 10, 0), false); // invalid size

    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), true); // invalid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), true); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    // ------------------------
    // Trypsin/P
    // ------------------------
    ed.setEnzyme(ProteaseDB::getInstance()->getEnzyme("Trypsin/P"));
    ed.setSpecificity(EnzymaticDigestion::SPEC_FULL); // require both sides

    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct("", 10, 0), false); // invalid size

    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), false); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), true); // valid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), false); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), false); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    // test with different missed cleavages when this is not ignored (ignore_missed_cleavages = false)
    //             |8  |12 |16|19
    prot = "ABCDEFGKABCRAAAKAARPBBBB"; // 4 cleavages at {(0),8,12,16,19}
    ed.setMissedCleavages(0); // redundant, by default zero, should be zero
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4, false), true);  //  valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8, false), false); //  invalid, fully-tryptic but with a missing cleavage
    ed.setMissedCleavages(1);
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8, false), true);  //  valid, fully-tryptic with 1 missing cleavage (allow)
    TEST_EQUAL(ed.isValidProduct(prot, 8, 11, false), false);//  invalid, fully-tryptic but with 2 missing cleavages
    ed.setMissedCleavages(2);
    TEST_EQUAL(ed.isValidProduct(prot, 8, 11, false), true); //  valid, fully-tryptic with 2 missing cleavages
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, true), true);  //  boundary case, length of protein (no checking of MCs)
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), false);//  boundary case, this exceeds missing cleavages
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19, false), false);//  start-boundary case, 2 allowed, 3 required
    ed.setMissedCleavages(3);
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), false);//  boundary case, invalid: 3 allowed, 4 required
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19, false), true); //  start-boundary case, 3 allowed, 3 required
    ed.setMissedCleavages(4); // maximum cleavages for this peptide
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), true); //  boundary case, accepted: 4 allowed, 4 required
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19, false), true); //  start-boundary case, 4 allowed, 3 required
    ed.setMissedCleavages(5); // allow even more ...
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), true); //  boundary case, accepted: 5 allowed, 4 required
    ed.setMissedCleavages(0); // set back to default

    //################################################
    // same as above, just with other specificity

    ed.setSpecificity(EnzymaticDigestion::SPEC_SEMI); // require one special cleavage site
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct("", 10, 0), false); // invalid size

    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), true); // valid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), true); // invalid N-term valid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    // test with different missed cleavages when this is not ignored (ignore_missed_cleavages = false)
    //             |8  |12 |16|19
    prot = "ABCDEFGKABCRAAAKAARPBBBB"; // 4 cleavages at {(0),8,12,16,19}
    ed.setMissedCleavages(0); // redundant, by default zero, should be zero
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3, false), true);  //  valid semi-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 5, false), false); //  invalid, semi-tryptic but with a missing cleavage
    ed.setMissedCleavages(1);
    TEST_EQUAL(ed.isValidProduct(prot, 8, 5, false), true);  //  valid, semi-tryptic with 1 missing cleavage (allow)
    TEST_EQUAL(ed.isValidProduct(prot, 8, 10, false), false);//  invalid, semi-tryptic but with 2 missing cleavages
    ed.setMissedCleavages(2);
    TEST_EQUAL(ed.isValidProduct(prot, 8, 10, false), true); //  valid, semi-tryptic with 2 missing cleavages
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, true), true);  //  boundary case, length of protein (no checking of MCs)
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), false);//  boundary case, this exceeds missing cleavages
    TEST_EQUAL(ed.isValidProduct(prot, 0, 18, false), false);//  start-boundary case, 2 allowed, 3 required
    ed.setMissedCleavages(3);
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), false);//  boundary case, invalid: 3 allowed, 4 required
    TEST_EQUAL(ed.isValidProduct(prot, 0, 18, false), true); //  start-boundary case, 3 allowed, 3 required
    ed.setMissedCleavages(4); // maximum cleavages for this peptide
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), true); //  boundary case, accepted: 4 allowed, 4 required
    TEST_EQUAL(ed.isValidProduct(prot, 0, 18, false), true); //  start-boundary case, 4 allowed, 3 required
    ed.setMissedCleavages(5); // allow even more ...
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), true); //  boundary case, accepted: 5 allowed, 4 required
    ed.setMissedCleavages(0); // set back to default

    //################################################
    // same as above, just with other specificity

    ed.setSpecificity(EnzymaticDigestion::SPEC_NONE); // require no special cleavage site
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct("", 10, 0), false); // invalid size

    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); // valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); // valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), true); // valid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), true); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing

    // test with different missed cleavages when this is not ignored (ignore_missed_cleavages = false)
    //             |8  |12 |16|19
    prot = "ABCDEFGKABCRAAAKAARPBBBB"; // 4 cleavages at {(0),8,12,16,19}
    ed.setMissedCleavages(0); // redundant, by default zero, should be zero
    TEST_EQUAL(ed.isValidProduct(prot, 9, 2, false), true);  //  valid not-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 9, 5, false), false); //  invalid, not-tryptic but with a missing cleavage
    ed.setMissedCleavages(1);
    TEST_EQUAL(ed.isValidProduct(prot, 9, 5, false), true);  //  valid, not-tryptic with 1 missing cleavage (allow)
    TEST_EQUAL(ed.isValidProduct(prot, 9, 9, false), false); //  invalid, semi-tryptic but with 2 missing cleavages
    ed.setMissedCleavages(2);
    TEST_EQUAL(ed.isValidProduct(prot, 9, 9, false), true);  //  valid, semi-tryptic with 2 missing cleavages
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, true), true);  //  boundary case, length of protein (no checking of MCs)
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), false);//  boundary case, this exceeds missing cleavages
    ed.setMissedCleavages(3);
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), false);//  boundary case, invalid: 3 allowed, 4 required
    ed.setMissedCleavages(4); // maximum cleavages for this peptide
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), true); //  boundary case, accepted: 4 allowed, 4 required
    ed.setMissedCleavages(5); // allow even more ...
    TEST_EQUAL(ed.isValidProduct(prot, 0, 24, false), true); //  boundary case, accepted: 5 allowed, 4 required
    ed.setMissedCleavages(0); // set back to default
}
END_SECTION

START_SECTION([EXTRA] Size countMissedCleavages_(const std::vector<int>& cleavage_positions, Size pep_start, Size pep_end) const)
  EnzymaticDigestion ed;
  ed.setMissedCleavages(2);
  TEST_EQUAL(ed.isValidProduct("KKKK", 0, 4, false), false); // has 3 MC's, should not be valid
  ed.setMissedCleavages(3);
  TEST_EQUAL(ed.isValidProduct("KKKK", 0, 4, false), true);  // has 3 MC's, should be valid
END_SECTION

START_SECTION(Size countInternalCleavageSites(const String& sequence) )
  EnzymaticDigestion ed;
  ed.setMissedCleavages(0); // setting max missed cleavages should not have any impact
  TEST_EQUAL(ed.countInternalCleavageSites("PEEKEEKEEPKEEPK"), 3); // has 3 internal cleavage sites
  ed.setMissedCleavages(2);
  TEST_EQUAL(ed.countInternalCleavageSites("PEEKEEKEEPKEEPK"), 3); // has 3 internal cleavage sites
  TEST_EQUAL(ed.countInternalCleavageSites("EEEEEEEEEEEEEEE"), 0); // has 0 internal cleavage sites
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
