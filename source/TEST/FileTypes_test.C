// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Stephan Aiche $
// $Authors: Marc Sturm, Andreas Bertsch, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

///////////////////////////

START_TEST(FileHandler, "Id")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((static String typeToName(FileTypes::Type type)))
{
  TEST_EQUAL(FileTypes::typeToName(FileTypes::UNKNOWN), "unknown");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::DTA), "dta");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::DTA2D), "dta2d");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::MZDATA), "mzData");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::MZXML), "mzXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::MZML), "mzML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::FEATUREXML), "featureXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::IDXML), "idXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::CONSENSUSXML), "consensusXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::TRANSFORMATIONXML), "trafoXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::INI), "ini");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::TOPPAS), "toppas");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::PNG), "png");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::TXT), "txt");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::CSV), "csv");
}
END_SECTION

START_SECTION((static FileTypes::Type nameToType(const String &name)))
{
  TEST_EQUAL(FileTypes::UNKNOWN, FileTypes::nameToType("unknown"));
  TEST_EQUAL(FileTypes::DTA, FileTypes::nameToType("dta"));
  TEST_EQUAL(FileTypes::DTA2D, FileTypes::nameToType("dta2d"));
  TEST_EQUAL(FileTypes::MZDATA, FileTypes::nameToType("mzData"));
  TEST_EQUAL(FileTypes::MZXML, FileTypes::nameToType("mzXML"));
  TEST_EQUAL(FileTypes::FEATUREXML, FileTypes::nameToType("featureXML"));
  TEST_EQUAL(FileTypes::IDXML, FileTypes::nameToType("idXmL")); // case-insensitivity
  TEST_EQUAL(FileTypes::CONSENSUSXML, FileTypes::nameToType("consensusXML"));
  TEST_EQUAL(FileTypes::MGF, FileTypes::nameToType("mgf"));
  TEST_EQUAL(FileTypes::INI, FileTypes::nameToType("ini"));
  TEST_EQUAL(FileTypes::TOPPAS, FileTypes::nameToType("toppas"));
  TEST_EQUAL(FileTypes::TRANSFORMATIONXML, FileTypes::nameToType("trafoXML"));
  TEST_EQUAL(FileTypes::MZML, FileTypes::nameToType("mzML"));
  TEST_EQUAL(FileTypes::MS2, FileTypes::nameToType("ms2"));
  TEST_EQUAL(FileTypes::PEPXML, FileTypes::nameToType("pepXML"));
  TEST_EQUAL(FileTypes::PROTXML, FileTypes::nameToType("protXML"));
  TEST_EQUAL(FileTypes::MZIDENTML, FileTypes::nameToType("mzid"));
  TEST_EQUAL(FileTypes::GELML, FileTypes::nameToType("gelML"));
  TEST_EQUAL(FileTypes::TRAML, FileTypes::nameToType("traML"));
  TEST_EQUAL(FileTypes::MSP, FileTypes::nameToType("msp"));
  TEST_EQUAL(FileTypes::OMSSAXML, FileTypes::nameToType("omssaXML"));
  TEST_EQUAL(FileTypes::PNG, FileTypes::nameToType("png"));
  TEST_EQUAL(FileTypes::XMASS, FileTypes::nameToType("fid"));
  TEST_EQUAL(FileTypes::TSV, FileTypes::nameToType("tsv"));
  TEST_EQUAL(FileTypes::PEPLIST, FileTypes::nameToType("peplist"));
  TEST_EQUAL(FileTypes::HARDKLOER, FileTypes::nameToType("hardkloer"));
  TEST_EQUAL(FileTypes::KROENIK, FileTypes::nameToType("kroenik"));
  TEST_EQUAL(FileTypes::FASTA, FileTypes::nameToType("fasta"));
  TEST_EQUAL(FileTypes::EDTA, FileTypes::nameToType("edta"));
  TEST_EQUAL(FileTypes::CSV, FileTypes::nameToType("csv"));
  TEST_EQUAL(FileTypes::TXT, FileTypes::nameToType("txt"));
}
END_SECTION

END_TEST
