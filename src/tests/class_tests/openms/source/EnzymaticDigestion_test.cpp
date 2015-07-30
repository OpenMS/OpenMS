// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(EnzymaticDigestion, "$Id$")

/////////////////////////////////////////////////////////////
    EnzymaticDigestion * e_ptr = 0;
    EnzymaticDigestion* e_nullPointer = 0;
    START_SECTION((EnzymaticDigestion()))
    e_ptr = new EnzymaticDigestion;
    TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION([EXTRA] ~EnzymaticDigestion())
    delete e_ptr;
END_SECTION

START_SECTION((EnzymaticDigestion(const EnzymaticDigestion &rhs)))
    EnzymaticDigestion ed;
    ed.setMissedCleavages(1234);
    ed.setEnzyme("no cleavage");
    ed.setSpecificity(EnzymaticDigestion::SPEC_SEMI);
    
    EnzymaticDigestion ed2(ed);
    
    TEST_EQUAL(ed.getMissedCleavages(), ed2.getMissedCleavages());
    TEST_EQUAL(ed.getEnzymeName(), ed2.getEnzymeName());
    TEST_EQUAL(ed.getSpecificity(), ed2.getSpecificity());

END_SECTION

START_SECTION((EnzymaticDigestion & operator=(const EnzymaticDigestion &rhs)))
    EnzymaticDigestion ed;
    ed.setMissedCleavages(1234);
    ed.setEnzyme("no cleavage");
    ed.setSpecificity(EnzymaticDigestion::SPEC_SEMI);
    
    EnzymaticDigestion ed2 = ed;
    
    TEST_EQUAL(ed.getMissedCleavages(), ed2.getMissedCleavages());
    TEST_EQUAL(ed.getEnzymeName(), ed2.getEnzymeName());
    TEST_EQUAL(ed.getSpecificity(), ed2.getSpecificity());

END_SECTION

START_SECTION((SignedSize getMissedCleavages() const))
    TEST_EQUAL(EnzymaticDigestion().getMissedCleavages(), 0)
END_SECTION

START_SECTION((Enzyme getEnzymeName() const))
    TEST_EQUAL(EnzymaticDigestion().getEnzymeName(), "Trypsin")
END_SECTION

START_SECTION((void setMissedCleavages(SignedSize missed_cleavages)))
    EnzymaticDigestion ed;
    ed.setMissedCleavages(5);
    TEST_EQUAL(ed.getMissedCleavages(), 5)
END_SECTION

START_SECTION((void setEnzyme(const String enzyme_name)))
    EnzymaticDigestion ed;
    ed.setEnzyme("Trypsin");
    TEST_EQUAL(ed.getEnzymeName(), "Trypsin");
    ed.setEnzyme("Trypsin/P");
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

START_SECTION((static Specificity getSpecificityByName(const String &name)))
    TEST_EQUAL(EnzymaticDigestion::getSpecificityByName(EnzymaticDigestion::NamesOfSpecificity[0]), EnzymaticDigestion::SPEC_FULL);
    TEST_EQUAL(EnzymaticDigestion::getSpecificityByName(EnzymaticDigestion::NamesOfSpecificity[1]), EnzymaticDigestion::SPEC_SEMI);
    TEST_EQUAL(EnzymaticDigestion::getSpecificityByName(EnzymaticDigestion::NamesOfSpecificity[2]), EnzymaticDigestion::SPEC_NONE);
    TEST_EQUAL(EnzymaticDigestion::getSpecificityByName("DoesNotExist"), EnzymaticDigestion::SIZE_OF_SPECIFICITY);
END_SECTION

START_SECTION((Size peptideCount(const AASequence &protein)))
    EnzymaticDigestion ed;
    for (int i = 0; i < 2; ++i) // common cases for Trypsin and Trypsin_P
    {
      if (i == 0) ed.setEnzyme("Trypsin");
      else if (i == 1)
        ed.setEnzyme("Trypsin/P");
    
      ed.setMissedCleavages(0);
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACDE")), 1)
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACKDE")), 2)
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACRDE")), 2)
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("ARCRDRE")), 4)
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("RKR")), 3)
      ed.setMissedCleavages(1);
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACDE")), 1)
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACRDE")), 3)
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("ARCDRE")), 5)
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("RKR")), 5)
      ed.setMissedCleavages(3);
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACDE")), 1)
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACRDE")), 3)
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("ARCDRE")), 6)
      TEST_EQUAL(ed.peptideCount(AASequence::fromString("RKR")), 6)
    }
    // special cases:
    ed.setMissedCleavages(0);
    ed.setEnzyme("Trypsin");
    TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACKPDE")), 1)
    TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACRPDE")), 1)
    TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACKPDERA")), 2)
    TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACRPDEKA")), 2)
    ed.setEnzyme("Trypsin/P");
    TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACKPDE")), 2)
    TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACRPDE")), 2)
    TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACKPDERA")), 3)
    TEST_EQUAL(ed.peptideCount(AASequence::fromString("ACRPDEKA")), 3)
END_SECTION

START_SECTION((void digest(const AASequence &protein, std::vector<AASequence>&output) const))
    EnzymaticDigestion ed;
    vector<AASequence> out;
    
    ed.digest(AASequence::fromString("ACDE"), out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].toString(), "ACDE")
    
    ed.digest(AASequence::fromString("ACKDE"), out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].toString(), "ACK")
    TEST_EQUAL(out[1].toString(), "DE")
    
    ed.digest(AASequence::fromString("ACRDE"), out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].toString(), "ACR")
    TEST_EQUAL(out[1].toString(), "DE")
    
    ed.digest(AASequence::fromString("ACKPDE"), out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].toString(), "ACKPDE")
    
    ed.digest(AASequence::fromString("ACRPDE"), out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].toString(), "ACRPDE")
    
    ed.digest(AASequence::fromString("ARCRDRE"), out);
    TEST_EQUAL(out.size(), 4)
    TEST_EQUAL(out[0].toString(), "AR")
    TEST_EQUAL(out[1].toString(), "CR")
    TEST_EQUAL(out[2].toString(), "DR")
    TEST_EQUAL(out[3].toString(), "E")
    
    ed.digest(AASequence::fromString("RKR"), out);
    TEST_EQUAL(out.size(), 3)
    TEST_EQUAL(out[0].toString(), "R")
    TEST_EQUAL(out[1].toString(), "K")
    TEST_EQUAL(out[2].toString(), "R")
    
    ed.setMissedCleavages(1);
    
    ed.digest(AASequence::fromString("ACDE"), out);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].toString(), "ACDE")
    
    ed.digest(AASequence::fromString("ACRDE"), out);
    TEST_EQUAL(out.size(), 3)
    TEST_EQUAL(out[0].toString(), "ACR")
    TEST_EQUAL(out[1].toString(), "DE")
    TEST_EQUAL(out[2].toString(), "ACRDE")
    
    ed.digest(AASequence::fromString("ARCDRE"), out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].toString(), "AR")
    TEST_EQUAL(out[1].toString(), "CDR")
    TEST_EQUAL(out[2].toString(), "E")
    TEST_EQUAL(out[3].toString(), "ARCDR")
    TEST_EQUAL(out[4].toString(), "CDRE")
    
    ed.digest(AASequence::fromString("RKR"), out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].toString(), "R")
    TEST_EQUAL(out[1].toString(), "K")
    TEST_EQUAL(out[2].toString(), "R")
    TEST_EQUAL(out[3].toString(), "RK")
    TEST_EQUAL(out[4].toString(), "KR")
    
    
    ed.digest(AASequence::fromString("(ICPL:2H(4))ARCDRE"), out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].toString(), "(ICPL:2H(4))AR")
    TEST_EQUAL(out[1].toString(), "CDR")
    TEST_EQUAL(out[2].toString(), "E")
    TEST_EQUAL(out[3].toString(), "(ICPL:2H(4))ARCDR")
    TEST_EQUAL(out[4].toString(), "CDRE")
    
    ed.digest(AASequence::fromString("ARCDRE(Amidated)"), out);
    TEST_EQUAL(out.size(), 5)
    TEST_EQUAL(out[0].toString(), "AR")
    TEST_EQUAL(out[1].toString(), "CDR")
    TEST_EQUAL(out[2].toString(), "E(Amidated)")
    TEST_EQUAL(out[3].toString(), "ARCDR")
    TEST_EQUAL(out[4].toString(), "CDRE(Amidated)")
    
    // ------------------------
    // Trypsin/P
    // ------------------------
    ed.setMissedCleavages(0);
    ed.setEnzyme("Trypsin/P");
    ed.digest(AASequence::fromString("ACKPDE"), out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].toString(), "ACK")
    TEST_EQUAL(out[1].toString(), "PDE")
    
    ed.digest(AASequence::fromString("ACRPDE"), out);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].toString(), "ACR")
    TEST_EQUAL(out[1].toString(), "PDE")


END_SECTION


START_SECTION((bool isValidProduct(const AASequence &protein, Size pep_pos, Size pep_length, bool methionine_cleavage)))
    EnzymaticDigestion ed;
    ed.setEnzyme("Trypsin");
    ed.setSpecificity(EnzymaticDigestion::SPEC_FULL); // require both sides
    
    AASequence prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size
    
    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), false); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); //   valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), false); // invalid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), false); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), false); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing
    
    
    prot = AASequence::fromString("MBCDEFGKABCRAAAKAA"); // starts with Met - we assume the cleaved form without Met occurs in vivo
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7, true), true); // valid N-term (since protein starts with Met)
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7, false), false); // invalid N-term (since Met cleavage is not allowed.)
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing
    
    //################################################
    // same as above, just with other specificity
    
    ed.setSpecificity(EnzymaticDigestion::SPEC_SEMI); // require one special cleavage site
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size
    
    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); //   valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), true); // invalid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing
    
    prot = AASequence::fromString("MBCDEFGKABCRAAAKAA"); // starts with Met - we assume the cleaved form without Met occurs in vivo
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7, true), true); // valid N-term (since protein starts with Met)
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7, false), true); // invalid N-term (since Met cleavage is not allowed.)
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing
    
    //################################################
    // same as above, just with other specificity
    
    ed.setSpecificity(EnzymaticDigestion::SPEC_NONE); // require no special cleavage site
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size
    
    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); //   valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), true); // invalid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), true); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing
    
    // ------------------------
    // Trypsin/P
    // ------------------------
    ed.setEnzyme("Trypsin/P");
    ed.setSpecificity(EnzymaticDigestion::SPEC_FULL); // require both sides
    
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size
    
    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), false); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); //   valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), true); //   valid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), false); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), false); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing
    
    prot = AASequence::fromString("MBCDEFGKABCRAAAKAA"); // starts with Met - we assume the cleaved form without Met occurs in vivo
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7, true), true); // valid N-term (since protein starts with Met)
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), false); // invalid N-term (since Met cleavage is not allowed.)
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing
    
    //################################################
    // same as above, just with other specificity
    
    ed.setSpecificity(EnzymaticDigestion::SPEC_SEMI); // require one special cleavage site
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size
    
    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); //   valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), true); //   valid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), false); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), true); // invalid N-term valid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing
    
    prot = AASequence::fromString("MBCDEFGKABCRAAAKAA"); // starts with Met - we assume the cleaved form without Met occurs in vivo
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7, true), true); // valid N-term (since protein starts with Met)
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7, false), true); // invalid N-term (since Met cleavage is not allowed.)
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing
    
    //################################################
    // same as above, just with other specificity
    
    ed.setSpecificity(EnzymaticDigestion::SPEC_NONE); // require no special cleavage site
    prot = AASequence::fromString("ABCDEFGKABCRAAAKAARPBBBB");
    TEST_EQUAL(ed.isValidProduct(prot, 100, 3), false); // invalid position
    TEST_EQUAL(ed.isValidProduct(prot, 10, 300), false); // invalid length
    TEST_EQUAL(ed.isValidProduct(prot, 10, 0), false); // invalid size
    TEST_EQUAL(ed.isValidProduct(AASequence::fromString(""), 10, 0), false); // invalid size
    
    TEST_EQUAL(ed.isValidProduct(prot, 0, 3), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, 8), true); //   valid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 8, 4), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 8, 8), true); //   valid fully-tryptic
    TEST_EQUAL(ed.isValidProduct(prot, 0, 19), true); //   valid C-term - followed by proline
    TEST_EQUAL(ed.isValidProduct(prot, 8, 3), true); // invalid C-term
    TEST_EQUAL(ed.isValidProduct(prot, 3, 6), true); // invalid C+N-term
    TEST_EQUAL(ed.isValidProduct(prot, 1, 7), true); // invalid N-term
    TEST_EQUAL(ed.isValidProduct(prot, 0, prot.size()), true); // the whole thing

END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
