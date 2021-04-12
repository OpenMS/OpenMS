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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>

///////////////////////////


namespace OpenMS
{
  class FeatureDeconvolutionTest
    : public FeatureDeconvolution 
  {
  public:
      /// List of adducts used to explain mass differences
      MassExplainer::AdductsType getPotentialAdducts()
      { return potential_adducts_;}
      /// labeling table
      Map<Size, String> getMapLabels()
      { return map_label_;}

      /// labeling table inverse
      Map<String, Size> getMapLabelInverse()
      { return map_label_inverse_;}

			/// status of intensity filter for edges
			bool isIntensityFilterEnabled()
      { return enable_intensity_filter_;}

			/// status of charge discovery
			CHARGEMODE getChargeMode()
      { return q_try_;}


  };

}

using namespace OpenMS;
using namespace std;

START_TEST(FeatureDeconvolution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureDeconvolution* ptr = nullptr;
FeatureDeconvolution* nullPointer = nullptr;
START_SECTION(FeatureDeconvolution())
	ptr = new FeatureDeconvolution();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~FeatureDeconvolution())
	delete ptr;
END_SECTION

START_SECTION([EXTRA](void updateMembers_()))
	FeatureDeconvolutionTest fdt;
	
  Param p;
	p.setValue("charge_min", 11, "minimal possible charge");
  p.setValue("charge_max", 13, "maximal possible charge");
  p.setValue("retention_max_diff", 1.0, "maximum allowed RT difference between any two features if their relation shall be determined");
	p.setValue("retention_max_diff_local", 2.0, "maxi");
  p.setValue("potential_adducts", ListUtils::create<String>("H:+:0.7,Na:+:0.1,(2)H4H-4:0:0.1:-2:heavy"), "Ad");
	fdt.setParameters(p);
  
  {
	MassExplainer::AdductsType adducts = fdt.getPotentialAdducts();
  Map<Size, String> map = fdt.getMapLabels();
  Map<String, Size> map_i = fdt.getMapLabelInverse();
  bool b_filter = fdt.isIntensityFilterEnabled();
  FeatureDeconvolution::CHARGEMODE cm = fdt.getChargeMode();

	TEST_EQUAL(adducts.size(), 3)
  TEST_EQUAL(adducts[0].getFormula(), "H1");
  TEST_EQUAL(adducts[0].getRTShift(), 0);
  TEST_EQUAL(adducts[0].getCharge(), 1);
  TEST_REAL_SIMILAR(adducts[0].getLogProb(), log(0.7));
  TEST_EQUAL(adducts[1].getFormula(), "Na1");
  TEST_EQUAL(adducts[1].getRTShift(), 0);
  TEST_EQUAL(adducts[1].getCharge(), 1);
  TEST_REAL_SIMILAR(adducts[1].getLogProb(), log(0.1));
  TEST_EQUAL(adducts[2].getFormula(), "(2)H4H-4");
  TEST_EQUAL(adducts[2].getRTShift(), -2);
  TEST_EQUAL(adducts[2].getCharge(), 0);
  TEST_REAL_SIMILAR(adducts[2].getLogProb(), log(0.1));
  TEST_EQUAL(cm, FeatureDeconvolution::QFROMFEATURE)
  TEST_EQUAL(map.size(), 2)
  TEST_EQUAL(map_i.size(), 2)
  TEST_EQUAL(map[0], "decharged features");
  TEST_EQUAL(map_i["decharged features"], 0);
  TEST_EQUAL(map[1], "heavy");
  TEST_EQUAL(map_i["heavy"], 1);
  TEST_EQUAL(b_filter, false)
  Param p_internal = fdt.getParameters();
  TEST_REAL_SIMILAR((double) p_internal.getValue("retention_max_diff"), 1.0);
  TEST_REAL_SIMILAR((double) p_internal.getValue("retention_max_diff_local"), 1.0);
  }

  // second param set
	p.setValue("charge_min", 11, "minimal possible charge");
  p.setValue("charge_max", 13, "maximal possible charge");
  p.setValue("q_try", "heuristic", "Try dif");
  p.setValue("potential_adducts", ListUtils::create<String>("H:+:0.9,Na:++:0.1"));
  p.setValue("retention_max_diff", 1.0, "maximum ");
	p.setValue("retention_max_diff_local", 1.0, "maxim");
  p.setValue("intensity_filter", "true", "Enable");
  p.setValue("default_map_label", "mylabel", "Label");
  p.setValue("retention_max_diff", 2.0, "maximum allowed RT difference between any two features if their relation shall be determined");
	p.setValue("retention_max_diff_local", 5.0, "maxi");

	fdt.setParameters(p);
  {
  MassExplainer::AdductsType adducts = fdt.getPotentialAdducts();
  Map<Size, String> map = fdt.getMapLabels();
  Map<String, Size> map_i = fdt.getMapLabelInverse();
  bool b_filter = fdt.isIntensityFilterEnabled();
  FeatureDeconvolution::CHARGEMODE cm = fdt.getChargeMode();

	TEST_EQUAL(adducts.size(), 2)
  TEST_EQUAL(adducts[0].getFormula(), "H1");
  TEST_EQUAL(adducts[0].getRTShift(), 0);
  TEST_EQUAL(adducts[0].getCharge(), 1);
  TEST_REAL_SIMILAR(adducts[0].getLogProb(), log(0.9));
  TEST_EQUAL(adducts[1].getFormula(), "Na1");
  TEST_EQUAL(adducts[1].getRTShift(), 0);
  TEST_EQUAL(adducts[1].getCharge(), 2);
  TEST_REAL_SIMILAR(adducts[1].getLogProb(), log(0.1));

  TEST_EQUAL(cm, FeatureDeconvolution::QHEURISTIC)
  TEST_EQUAL(map.size(), 1)
  TEST_EQUAL(map_i.size(), 1)
  TEST_EQUAL(map[0], "mylabel");
  TEST_EQUAL(map_i["mylabel"], 0);
  TEST_EQUAL(b_filter, true)
  Param p_internal = fdt.getParameters();
  TEST_REAL_SIMILAR((double) p_internal.getValue("retention_max_diff"), 2.0);
  TEST_REAL_SIMILAR((double) p_internal.getValue("retention_max_diff_local"), 2.0);

  }

END_SECTION

START_SECTION(FeatureDeconvolution(const FeatureDeconvolution &source))
	FeatureDeconvolution fd;
	Param p;
	p.setValue("charge_min", 11, "minimal possible charge");
  p.setValue("charge_max", 13, "maximal possible charge");
	fd.setParameters(p);
	FeatureDeconvolution fd2(fd);
	FeatureDeconvolution fd_untouched;
	
	TEST_EQUAL(fd2.getParameters(), fd.getParameters())
	TEST_NOT_EQUAL(fd2.getParameters(), fd_untouched.getParameters())

END_SECTION

START_SECTION(FeatureDeconvolution& operator=(const FeatureDeconvolution &source))
	FeatureDeconvolution fd;
	Param p;
	p.setValue("charge_min", 11, "minimal possible charge");
  p.setValue("charge_max", 13, "maximal possible charge");
	fd.setParameters(p);
	FeatureDeconvolution fd2 = fd;
	FeatureDeconvolution fd_untouched;
	
	TEST_EQUAL(fd2.getParameters(), fd.getParameters())
	TEST_NOT_EQUAL(fd2.getParameters(), fd_untouched.getParameters())
END_SECTION


START_SECTION(void compute(const FeatureMapType &fm_in, FeatureMapType &fm_out, ConsensusMap &cons_map, ConsensusMap &cons_map_p))
//_CrtSetDbgFlag(_CrtSetDbgFlag(0)|_CRTDBG_CHECK_ALWAYS_DF);

	FeatureDeconvolution fd;
        Param p;
        p.setValue("potential_adducts", ListUtils::create<String>("H:+:0.7,Na:+:0.1,(2)H4H-4:0:0.1:-2:heavy"), "Ad");
	p.setValue("mass_max_diff", 0.1);
	fd.setParameters(p);

	FeatureMap fm_in, fm_out;
	ConsensusMap cm, cm2;
	FeatureXMLFile fl;
	fl.load(OPENMS_GET_TEST_DATA_PATH("FeatureDeconvolution_easy_input.featureXML"), fm_in);
	fd.compute(fm_in, fm_out, cm, cm2);

	String out_file;
	NEW_TMP_FILE(out_file)
	ConsensusXMLFile c1;

	c1.store(out_file,cm);

  WHITELIST("xml-stylesheet");
  // WHITELIST("xml-stylesheet,consensusElement id=");
	// WHITELIST("xml-stylesheet,map id,consensusElement id=");
	TEST_FILE_SIMILAR(out_file, OPENMS_GET_TEST_DATA_PATH("FeatureDeconvolution_easy_output.consensusXML"));

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


