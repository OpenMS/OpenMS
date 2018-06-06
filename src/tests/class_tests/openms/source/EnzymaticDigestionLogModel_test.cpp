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
// $Maintainer: Chris Bielow, Xiao Liang $
// $Authors: Marc Sturm, Chris Bielow, Xiao Liang $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/EnzymaticDigestionLogModel.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(EnzymaticDigestionLogModel, "$Id$")

/////////////////////////////////////////////////////////////

EnzymaticDigestionLogModel * e_ptr = nullptr;
EnzymaticDigestionLogModel* e_nullPointer = nullptr;
START_SECTION((EnzymaticDigestionLogModel()))
e_ptr = new EnzymaticDigestionLogModel;
TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((EnzymaticDigestionLogModel(const EnzymaticDigestionLogModel &rhs)))
EnzymaticDigestionLogModel ed;
ed.setEnzyme("no cleavage");
ed.setLogThreshold(81231);

EnzymaticDigestionLogModel ed2(ed);

TEST_EQUAL(ed.getEnzymeName(), ed2.getEnzymeName());
TEST_EQUAL(ed.getLogThreshold(), ed2.getLogThreshold());

END_SECTION

START_SECTION((EnzymaticDigestionLogModel & operator=(const EnzymaticDigestionLogModel &rhs)))
EnzymaticDigestionLogModel ed;
ed.setEnzyme("no cleavage");
ed.setLogThreshold(81231);

EnzymaticDigestionLogModel ed2 = ed;

TEST_EQUAL(ed.getEnzymeName(), ed2.getEnzymeName());
TEST_EQUAL(ed.getLogThreshold(), ed2.getLogThreshold());

END_SECTION

START_SECTION((Enzyme getEnzymeName() const))
TEST_EQUAL(EnzymaticDigestionLogModel().getEnzymeName(), "Trypsin")
END_SECTION

START_SECTION((void setEnzyme(const String enzyme_name)))
EnzymaticDigestionLogModel ed;
ed.setEnzyme("Trypsin");
TEST_EQUAL(ed.getEnzymeName(), "Trypsin");
END_SECTION

START_SECTION((double getLogThreshold() const))
EnzymaticDigestionLogModel ed;
ed.setLogThreshold(1.234);
TEST_EQUAL(ed.getLogThreshold(), 1.234);
END_SECTION

START_SECTION((void setLogThreshold(double threshold)))
// TESTED ABOVE
NOT_TESTABLE
END_SECTION

START_SECTION((Size peptideCount(const AASequence &protein)))
EnzymaticDigestionLogModel ed;
// with log L model:
ed.setEnzyme("Trypsin");
TEST_EQUAL(ed.peptideCount(AASequence::fromString("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFDKCCAADDKEACFAVEGPKLVVSTQTALA")), 9 + 1 + 1)   // K R + 1
// with non-standard amino-acids "O" and "U":
TEST_EQUAL(ed.peptideCount(AASequence::fromString("IITAQVUDRPONAIYMTY")), 2);
END_SECTION

START_SECTION((void digest(const AASequence &protein, std::vector<AASequence>&output) const))
EnzymaticDigestionLogModel ed;
vector<AASequence> out;
// with log L model:
ed.digest(AASequence::fromString("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFDKCCAADDKEACFAVEGPKLVVSTQTALA"), out);
TEST_EQUAL(out.size(), 11)
TEST_EQUAL(out[0].toString(), "MKWVTFISLLLLFSSAYSRGVFRRDTHK")
TEST_EQUAL(out[1].toString(), "SEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLK")
TEST_EQUAL(out[2].toString(), "PDPNTLCDEFKADEKK")
TEST_EQUAL(out[3].toString(), "FWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQK")
TEST_EQUAL(out[4].toString(), "FPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDK")
TEST_EQUAL(out[5].toString(), "PLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRK")
TEST_EQUAL(out[6].toString(), "VPQVSTPTLVEVSRSLGK")
TEST_EQUAL(out[7].toString(), "VGTRCCTK")
TEST_EQUAL(out[8].toString(), "PESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRR")
TEST_EQUAL(out[9].toString(), "PCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHK")
TEST_EQUAL(out[10].toString(), "PKATEEQLKTVMENFVAFDKCCAADDKEACFAVEGPKLVVSTQTALA")

// ------------------------
// Trypsin/P
// ------------------------
ed.setEnzyme("Trypsin/P");
// .. log-model only for restrictive Trypsin (with P constraint)
TEST_EXCEPTION(Exception::InvalidParameter, ed.digest(AASequence::fromString("ANGER"), out));

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
