// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MzTab.h>
///////////////////////////

START_TEST(MzTab, "$Id$")

using namespace OpenMS;
using namespace std;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzTab* ptr = nullptr;
MzTab* null_ptr = nullptr;
START_SECTION(MzTab())
{
	ptr = new MzTab();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MzTab())
{
	delete ptr;
}
END_SECTION

START_SECTION(std::vector<String> getPSMOptionalColumnNames() const)
{
  MzTab mztab;
  MzTabPSMSectionRow row;
  MzTabPSMSectionRows rows;
  MzTabString s;
  MzTabOptionalColumnEntry e;

  // row 1 //////////////////////
  row.sequence.fromCellString("NDYKAPPQPAPGK");
  row.PSM_ID.fromCellString("38");
  row.accession.fromCellString("IPI:B1");
  row.unique.fromCellString("1");
  row.database.fromCellString("null");
  row.database_version.fromCellString("null");
  row.search_engine.fromCellString("[, , Percolator, ]");
  row.search_engine_score[0].fromCellString("51.9678841193106");

  e.first = "Percolator_score";
  s.fromCellString("0.359083");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "Percolator_qvalue";
  s.fromCellString("0.00649874");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "Percolator_PEP";
  s.fromCellString("0.0420992");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "search_engine_sequence";
  s.fromCellString("NDYKAPPQPAPGK");
  e.second = s;
  row.opt_.push_back(e);

  rows.push_back(row);

  // row 2 //////////////////////
  row.sequence.fromCellString("IRRS(Phospho)SFSSK");
  row.PSM_ID.fromCellString("39");
  row.accession.fromCellString("IPI:IPI00009899.4");
  row.unique.fromCellString("0");
  row.database.fromCellString("null");
  row.database_version.fromCellString("null");
  row.search_engine.fromCellString("[, , Percolator, ]");
  row.search_engine_score[0].fromCellString("9.55915773892318");

  e.first = "Percolator_score";
  s.fromCellString("0.157068");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "Percolator_qvalue";
  s.fromCellString("0.00774619");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "Percolator_PEP";
  s.fromCellString("0.0779777");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "search_engine_sequence";
  s.fromCellString("IRRSSFS(Phospho)SK");
  e.second = s;
  row.opt_.push_back(e);

  e.first = "AScore_1";
  s.fromCellString("3.64384830671351");
  e.second = s;
  row.opt_.push_back(e);
  rows.push_back(row);

  mztab.setPSMSectionRows(rows);

  // Tests ///////////////////////////////
  vector<String> optional_columns = mztab.getPSMOptionalColumnNames();

  TEST_EQUAL(mztab.getPSMSectionRows().size(),2)
  TEST_EQUAL(optional_columns.size(),5)
}
END_SECTION

START_SECTION(static void addMetaInfoToOptionalColumns(const std::set<String>& keys, std::vector<MzTabOptionalColumnEntry>& opt, const String& id, const MetaInfoInterface& meta))
{
  // keys will have spaces replaced with underscore
  std::set<String> keys{ "FWHM", "with space", "ppm_errors" };
  // values should remain as they are
  MetaInfoInterface meta;
  meta.setMetaValue("FWHM", 34.5);
  DataValue dv(DoubleList{ 0.5, 1.4, -2.0, 0.1 });
  meta.setMetaValue("ppm_errors", dv);
  std::vector<MzTabOptionalColumnEntry> opt;
  MzTab::addMetaInfoToOptionalColumns(keys, opt, "global", meta);
  TEST_EQUAL(opt.size(), 3)
  TEST_EQUAL(opt[0].first, "opt_global_FWHM");
  TEST_EQUAL(opt[1].first, "opt_global_ppm_errors");
  TEST_EQUAL(opt[2].first, "opt_global_with_space");
  TEST_EQUAL(opt[0].second.toCellString(), "34.5");
  TEST_EQUAL(opt[1].second.toCellString(), "[0.5, 1.4, -2.0, 0.1]");
  TEST_EQUAL(opt[2].second.toCellString(), "null");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
