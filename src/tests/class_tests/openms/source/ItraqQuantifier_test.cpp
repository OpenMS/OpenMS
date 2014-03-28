// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
///////////////////////////

//#include <OpenMS/DATASTRUCTURES/Matrix.h>

using namespace OpenMS;
using namespace std;

ConsensusFeature getCFWithIntensites(double v[])
{
  ConsensusFeature cf;
  BaseFeature bf0, bf1, bf2, bf3;
  bf0.setIntensity(v[0]);
  bf1.setIntensity(v[1]);
  bf2.setIntensity(v[2]);
  bf3.setIntensity(v[3]);
  cf.insert(0, bf0);cf.insert(1, bf1);cf.insert(2, bf2);cf.insert(3, bf3);
  cf.setIntensity(v[0]+v[1]+v[2]+v[3]);
  return cf;
}

START_TEST(ItraqQuantifier, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ItraqQuantifier* ptr = 0;
ItraqQuantifier* nullPointer = 0;
START_SECTION(ItraqQuantifier())
{
	ptr = new ItraqQuantifier();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ItraqQuantifier())
{
	delete ptr;
}
END_SECTION

START_SECTION((ItraqQuantifier(Int itraq_type)))
{
  ItraqQuantifier iq(ItraqQuantifier::EIGHTPLEX);
	TEST_EQUAL((String) iq.getParameters().getValue("isotope_correction")=="true", true);
	TEST_EQUAL((Int) iq.getParameters().getValue("channel_reference"), 114);
  ItraqQuantifier iq2(ItraqQuantifier::FOURPLEX);
	TEST_EQUAL((String) iq2.getParameters().getValue("isotope_correction")=="true", true);
	TEST_EQUAL((Int) iq2.getParameters().getValue("channel_reference"), 114);
}
END_SECTION

START_SECTION((ItraqQuantifier(Int itraq_type, const Param &param)))
{
	Param p;
	p.setValue("isotope_correction:4plex", ListUtils::create<String>("114:0/0.3/4/0 , 116:0.1/0.3/3/0.2"));
	ItraqQuantifier iq(ItraqQuantifier::FOURPLEX, p);
	TEST_EQUAL(iq.getParameters().getValue("isotope_correction:4plex") == ListUtils::create<String>( "114:0/0.3/4/0 , 116:0.1/0.3/3/0.2"), true)

	// this should go wrong
	p.setValue("isotope_correction:4plex", ListUtils::create<String>("114:0/0.3/0 , 116:0.1/0.3/3/0.2"));
	TEST_EXCEPTION(Exception::InvalidParameter, ItraqQuantifier iq2(ItraqQuantifier::FOURPLEX, p));

	// this should go wrong too
	p.setValue("isotope_correction:4plex", ListUtils::create<String>("113:0/0.3/0/0.3 , 116:0.1/0.3/3/0.2"));
	TEST_EXCEPTION(Exception::InvalidParameter, ItraqQuantifier iq2(ItraqQuantifier::FOURPLEX, p));
}
END_SECTION


START_SECTION((ItraqQuantifier(const ItraqQuantifier &cp)))
{
	Param p;
	p.setValue("isotope_correction:4plex", ListUtils::create<String>("114:0/0.3/4/0 , 116:0.1/0.3/3/0.2"));
	ItraqQuantifier iq(ItraqQuantifier::EIGHTPLEX, p);

	ItraqQuantifier iq_cp(iq);

	TEST_EQUAL(iq_cp.getParameters(), iq.getParameters());

}
END_SECTION

START_SECTION((ItraqQuantifier& operator=(const ItraqQuantifier &rhs)))
{
	Param p;
	p.setValue("isotope_correction:4plex", ListUtils::create<String>("114:0/0.3/4/0 , 116:0.1/0.3/3/0.2"));
	ItraqQuantifier iq(ItraqQuantifier::EIGHTPLEX, p);

	ItraqQuantifier iq_cp;
	iq_cp = iq;

	TEST_EQUAL(iq_cp.getParameters(), iq.getParameters());

}
END_SECTION


START_SECTION((void run(const ConsensusMap &consensus_map_in, ConsensusMap &consensus_map_out)))
{
	ConsensusXMLFile cm_file;
	ConsensusMap cm_in, cm_out;
	cm_file.load(OPENMS_GET_TEST_DATA_PATH("ItraqChannelExtractor.consensusXML"),cm_in);

	ItraqQuantifier iq;
	Param p;
	p.setValue("do_normalization", "true");
	iq.setParameters(p);
	iq.run(cm_in,cm_out);

	String cm_file_out;// = OPENMS_GET_TEST_DATA_PATH("ItraqQuantifier.consensusXML");
	NEW_TMP_FILE(cm_file_out);
	cm_file.store(cm_file_out,cm_out);

	WHITELIST("<?xml-stylesheet,id=");
	// WHITELIST("<?xml-stylesheet,consensusElement id=");
	TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("ItraqQuantifier.consensusXML"));
}
END_SECTION


START_SECTION((ItraqQuantifierStats getStats() const))
{
  /*
  // prep: for data generation
  ItraqConstants::IsotopeMatrices isotope_corrections_;
  isotope_corrections_.resize(2);
  isotope_corrections_[0].setMatrix<4,4>(ItraqConstants::ISOTOPECORRECTIONS_FOURPLEX);
  isotope_corrections_[1].setMatrix<8,4>(ItraqConstants::ISOTOPECORRECTIONS_EIGHTPLEX);
  Matrix<double> channel_frequency = ItraqConstants::translateIsotopeMatrix(0, isotope_corrections_);
  std::cerr << "matrix: \n\n" << channel_frequency << "\n\n";

  // some R code to  get the nnls and naive values faster


  require("nnls")
  ## correction matrix (as obtained from 'channel_frequency' matrix above)
  m = matrix(c(0.929, 0.02, 0, 0,
  0.059, 0.923,  0.03, 0.001,
  0.002, 0.056, 0.924,  0.04,
  0    , 0.001, 0.045, 0.923), ncol=4, nrow=4, byrow=T)
  ## 'true' intensities
  x1 = c(-1,100,100,100)
  ## observed intensities
  i = m %*% x1 ##   1.071  95.341  101.998  96.900

  ## naive and nnls solution
  n = solve(m) %*% i ##  -1       100      100        100
  nn = nnls(m, i)$x  ##  0.00000  99.91414 100.00375  99.99990

  d = n-nn
  sum(abs(d[2:4]))

  */
  ConsensusXMLFile cm_file;
  ConsensusMap cm_in, cm_out;
  cm_file.load(OPENMS_GET_TEST_DATA_PATH("ItraqChannelExtractor.consensusXML"), cm_in);

  cm_in.clear(false);

  ItraqQuantifier iq;
  Param p;
  p.setValue("isotope_correction", "true");
  p.setValue("do_normalization", "false");
  iq.setParameters(p);

  // first run (empty):
  iq.run(cm_in, cm_out);

  ItraqQuantifier::ItraqQuantifierStats stats = iq.getStats();
  TEST_EQUAL(stats.channel_count, 4)
  TEST_EQUAL(stats.iso_number_ms2_negative, 0)
  TEST_EQUAL(stats.iso_number_reporter_negative, 0)
  TEST_EQUAL(stats.iso_number_reporter_different, 0)
  TEST_REAL_SIMILAR(stats.iso_solution_different_intensity,  0)
  TEST_REAL_SIMILAR(stats.iso_total_intensity_negative, 0)
  TEST_EQUAL(stats.number_ms2_total, cm_in.size())
  TEST_EQUAL(stats.number_ms2_empty, 0)
  TEST_EQUAL(stats.empty_channels[114], 0)
  TEST_EQUAL(stats.empty_channels[115], 0)
  TEST_EQUAL(stats.empty_channels[116], 0)
  TEST_EQUAL(stats.empty_channels[117], 0)


  // add some target results
  double v1[4] = {1.071,  95.341,  101.998,  96.900}; // naive yields: {-1,100,100,100};  NNLS: {0.00000  99.91414 100.00375  99.99990}
  cm_in.push_back(getCFWithIntensites(v1));

  iq.run(cm_in, cm_out);

  stats = iq.getStats();
  TEST_EQUAL(stats.channel_count, 4)
  TEST_EQUAL(stats.iso_number_ms2_negative, 1)
  TEST_EQUAL(stats.iso_number_reporter_negative, 1)
  TEST_EQUAL(stats.iso_number_reporter_different, 3)
  TEST_REAL_SIMILAR(stats.iso_solution_different_intensity, 0.089703566418)
  TEST_REAL_SIMILAR(stats.iso_total_intensity_negative, 299.9178)
  TEST_EQUAL(stats.number_ms2_total, cm_in.size())
  TEST_EQUAL(stats.number_ms2_empty, 0)
  TEST_EQUAL(stats.empty_channels[114], 1)
  TEST_EQUAL(stats.empty_channels[115], 0)
  TEST_EQUAL(stats.empty_channels[116], 0)
  TEST_EQUAL(stats.empty_channels[117], 0)

  // change some more... (second run)
  double v2[4] = {0,0,0,0};
  cm_in.push_back(getCFWithIntensites(v2));

  iq.run(cm_in, cm_out);

  stats = iq.getStats();
  TEST_EQUAL(stats.channel_count, 4)
  TEST_EQUAL(stats.iso_number_ms2_negative, 1)
  TEST_EQUAL(stats.iso_number_reporter_negative, 1)
  TEST_EQUAL(stats.iso_number_reporter_different, 3)
  TEST_REAL_SIMILAR(stats.iso_solution_different_intensity, 0.089703566418)
  TEST_REAL_SIMILAR(stats.iso_total_intensity_negative, 299.9178)
  TEST_EQUAL(stats.number_ms2_total, cm_in.size())
  TEST_EQUAL(stats.number_ms2_empty, 1)
  TEST_EQUAL(stats.empty_channels[114], 2)
  TEST_EQUAL(stats.empty_channels[115], 1)
  TEST_EQUAL(stats.empty_channels[116], 1)
  TEST_EQUAL(stats.empty_channels[117], 1)

  // grab some warning, such that it does not show in output
  // Macro#1 (use local scope to hide 'ss_grab' to avoid double declaration)
  // params: log steam, default (to be temporarily disabled), message (to grab))
  // WATCH_LOG(Log_warn, cout, "msg...")
  {
  ostringstream ss_grab;
  Log_warn.flush();
  Log_warn.insert(ss_grab);
  Log_warn.remove(cout);

  p.setValue("isotope_correction", "false");
  iq.setParameters(p);
  iq.run(cm_in, cm_out);

  // test result of warning
  // Macro#2 (check msg and close local scope)
  Log_warn.flush();
  TEST_EQUAL(ss_grab.str(), "Warning: Due to deactivated isotope-correction labeling statistics will be based on raw intensities, which might give too optimistic results.\n")
  Log_warn.remove(ss_grab);
  Log_warn.insert(cout);
  }

  stats = iq.getStats();
  TEST_EQUAL(stats.channel_count, 4)
  TEST_EQUAL(stats.iso_number_ms2_negative, 0)
  TEST_EQUAL(stats.iso_number_reporter_negative, 0)
  TEST_EQUAL(stats.iso_number_reporter_different, 0)
  TEST_REAL_SIMILAR(stats.iso_solution_different_intensity, 0)
  TEST_REAL_SIMILAR(stats.iso_total_intensity_negative, 0)
  TEST_EQUAL(stats.number_ms2_total, cm_in.size())
  TEST_EQUAL(stats.number_ms2_empty, 1)
  TEST_EQUAL(stats.empty_channels[114], 1)
  TEST_EQUAL(stats.empty_channels[115], 1)
  TEST_EQUAL(stats.empty_channels[116], 1)
  TEST_EQUAL(stats.empty_channels[117], 1)

}
END_SECTION

START_SECTION([ItraqQuantifier::ItraqQuantifierStats] ItraqQuantifierStats())

  ItraqQuantifier::ItraqQuantifierStats stats;

  // ... this is an unimportant test, as values are filled during run() method. Test it there...
  TEST_EQUAL(stats.channel_count, 0)
  TEST_EQUAL(stats.number_ms2_empty, 0)
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



