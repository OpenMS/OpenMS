// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
ProteaseDigestion* pd_ptr = 0;
ProteaseDigestion* pd_null = 0;

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

START_SECTION((void digest(const AASequence &protein, std::vector<AASequence>&output) const))
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

START_SECTION((bool isValidProduct(const String& protein, Size pep_pos, Size pep_length, bool ignore_missed_cleavages, bool methionine_cleavage)))
    NOT_TESTABLE // tested by overload below
END_SECTION

START_SECTION((bool isValidProduct(const AASequence& protein, Size pep_pos, Size pep_length, bool ignore_missed_cleavages, bool methionine_cleavage)))
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
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7, true, false), false); // invalid N-term (since Met cleavage is not allowpd.)
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
    TEST_EQUAL(pd.isValidProduct(prot, 1, 7, true, false), true); // invalid N-term (since Met cleavage is not allowpd.)
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

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
