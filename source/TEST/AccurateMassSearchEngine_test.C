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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzTab.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AccurateMassSearchEngine, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AccurateMassSearchEngine* ptr = 0;
AccurateMassSearchEngine* null_ptr = 0;
START_SECTION(AccurateMassSearchEngine())
{
    ptr = new AccurateMassSearchEngine();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~AccurateMassSearchEngine())
{
    delete ptr;
}
END_SECTION


AccurateMassSearchEngine ams_pos;
AccurateMassSearchEngine ams_neg;

DoubleReal query_mass_pos(308.09);
DoubleReal query_mass_neg(306.08);
// Param ams_param;
// ams_param.setValue("");

String id_list_pos[] = {"C10H17N3O6S", "C15H16O7", "C14H14N2OS2", "C16H15NO4", "C17H11N5", "C10H14NO6P", "C14H12O4", "C7H6O2"};
String id_list_neg[] = {/*"C32H28O11",*/ "C17H17Cl2N", "C10H13N5O5", "C6H14O6S2"};
                        // 588.16316173 this cannot be since its +2, but we restrict to +1 (or -1)
                        //                   305.073804963  283.091668551   246.02317956
START_SECTION((void queryByMass(const DoubleReal& observed_mass, const Int& observed_charge, std::vector<AccurateMassSearchResult>& results)))
{
    std::vector<AccurateMassSearchResult> hmdb_results_pos, hmdb_results_neg;
    ams_pos.queryByMass(query_mass_pos, 1.0, hmdb_results_pos);

    Size id_list_pos_length(sizeof(id_list_pos)/sizeof(id_list_pos[0]));
    Size id_list_neg_length(sizeof(id_list_neg)/sizeof(id_list_neg[0]));


    TEST_EQUAL(hmdb_results_pos.size(), id_list_pos_length)

    if (hmdb_results_pos.size() == id_list_pos_length)
    {
        for (Size i = 0; i < id_list_pos_length; ++i)
        {
            TEST_STRING_EQUAL(hmdb_results_pos[i].getFormulaString(), id_list_pos[i])
            // std::cout << hmdb_results_pos[i].getFormulaString() << std::endl;
        }
    }
    ams_pos.queryByMass(query_mass_neg, -1.0, hmdb_results_neg);  // this is not 100% correct, since we are still searching with positive adducts. 
    // However, it does not matter if +z or -z, since any FF will just give +1, even in negative mode

    TEST_EQUAL(hmdb_results_neg.size(), id_list_neg_length)

    if (hmdb_results_neg.size() == id_list_neg_length)
    {
        for (Size i = 0; i < hmdb_results_neg.size(); ++i)
        {
            TEST_STRING_EQUAL(hmdb_results_neg[i].getFormulaString(), id_list_neg[i])
            // std::cout << hmdb_results_neg[i].getFormulaString() << std::endl;
        }
    }
}
END_SECTION

Feature test_feat;
test_feat.setRT(300.0);
test_feat.setMZ(399.33486);
test_feat.setIntensity(100.0);
test_feat.setMetaValue("num_of_masstraces", 3);
test_feat.setCharge(1.0);

test_feat.setMetaValue("masstrace_intensity_0", 100.0);
test_feat.setMetaValue("masstrace_intensity_1", 26.1);
test_feat.setMetaValue("masstrace_intensity_2", 4.0);

AccurateMassSearchEngine ams_feat_test;

String feat_query_pos[] = {"C23H45NO4", "C20H37NO3", "C22H41NO"};

START_SECTION((void queryByFeature(const Feature& feature, const Size& feature_index, std::vector<AccurateMassSearchResult>& results)))
{
    std::vector<AccurateMassSearchResult> results;

    ams_feat_test.queryByFeature(test_feat, 0, results);

    TEST_EQUAL(results.size(), 3)

    for (Size i = 0; i < results.size(); ++i)
    {
        TEST_REAL_SIMILAR(results[i].getObservedRT(), 300.0)
        TEST_REAL_SIMILAR(results[i].getObservedIntensity(), 100.0)
    }

    Size feat_query_size(sizeof(feat_query_pos)/sizeof(feat_query_pos[0]));

    if (results.size() == feat_query_size)
    {
        for (Size i = 0; i < feat_query_size; ++i)
        {
            TEST_STRING_EQUAL(results[i].getFormulaString(), feat_query_pos[i])
        }
    }
}
END_SECTION

ConsensusFeature cons_feat;
cons_feat.setRT(300.0);
cons_feat.setMZ(399.33486);
cons_feat.setIntensity(100.0);
cons_feat.setCharge(1.0);

FeatureHandle fh1, fh2, fh3;
fh1.setRT(300.0);
fh1.setMZ(399.33485);
fh1.setIntensity(100.0);
fh1.setCharge(1.0);
fh1.setMapIndex(0);

fh2.setRT(310.0);
fh2.setMZ(399.33486);
fh2.setIntensity(300.0);
fh2.setCharge(1.0);
fh2.setMapIndex(1);

fh3.setRT(290.0);
fh3.setMZ(399.33487);
fh3.setIntensity(500.0);
fh3.setCharge(1.0);
fh3.setMapIndex(2);

cons_feat.insert(fh1);
cons_feat.insert(fh2);
cons_feat.insert(fh3);
cons_feat.computeConsensus();


START_SECTION((void queryByConsensusFeature(const ConsensusFeature& cfeat, const Size& cf_index, const Size& number_of_maps, std::vector<AccurateMassSearchResult>& results)))
{
    std::vector<AccurateMassSearchResult> results;

    ams_feat_test.queryByConsensusFeature(cons_feat, 0, 3, results);

    TEST_EQUAL(results.size(), 3)

    for (Size i = 0; i < results.size(); ++i)
    {
        TEST_REAL_SIMILAR(results[i].getObservedRT(), 300.0)
        TEST_REAL_SIMILAR(results[i].getObservedIntensity(), 0.0)
    }

    // std::cout << cons_feat.getMZ() << " " << results.size() << std::endl;

    for (Size i = 0; i < results.size(); ++i)
    {
        std::vector<DoubleReal> indiv_ints = results[i].getIndividualIntensities();
        TEST_EQUAL(indiv_ints.size(), 3)

        if (indiv_ints.size() == 3)
        {
            TEST_REAL_SIMILAR(indiv_ints[0], fh1.getIntensity());
            TEST_REAL_SIMILAR(indiv_ints[1], fh2.getIntensity());
            TEST_REAL_SIMILAR(indiv_ints[2], fh3.getIntensity());
        }
    }

    Size feat_query_size(sizeof(feat_query_pos)/sizeof(feat_query_pos[0]));

    if (results.size() == feat_query_size)
    {
        for (Size i = 0; i < feat_query_size; ++i)
        {
            TEST_STRING_EQUAL(results[i].getFormulaString(), feat_query_pos[i])
        }
    }
}
END_SECTION

FeatureMap<> exp_fm;
FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.featureXML"), exp_fm);
MzTab test_mztab;

Param ams_param;
ams_param.setValue("isotopic_similarity", "true");

String fm_id_list[] = {"C17H25ClO2","C10H19N3O4S","C10H19N3O4S","C10H19N3O4S","C10H19N3O4S","C18H15NO2","C18H15NO2","C8H16N2O4S","C8H16N2O4S","C10H9NO","C8H8N4","C10H9NO","C8H8N4","C8H18NO2","C8H18NO2","C10H8O2","C10H8O2","C10H8O2","C17H20N2S","C17H20N2S","C17H20N2S","C17H20N2S"};
String fm_id_filt_list[] = {"C17H25ClO2","C10H19N3O4S","C10H19N3O4S","C10H19N3O4S","C10H19N3O4S","C18H15NO2","C18H15NO2","C8H16N2O4S","C8H16N2O4S","C10H9NO","C8H8N4","C10H9NO","C8H18NO2","C8H18NO2","C10H8O2","C10H8O2","C10H8O2","C17H20N2S","C17H20N2S","C17H20N2S","C17H20N2S"};


START_SECTION((void run(const FeatureMap<> &, MzTab &)))
{
    ams_feat_test.run(exp_fm, test_mztab);

    MzTabSmallMoleculeSectionData sms = test_mztab.getSmallMoleculeSectionData();
    MzTabSmallMoleculeSectionRows sm_rows = sms["AccMassSearch"];

    Size fm_id_list_size = sizeof(fm_id_list)/sizeof(fm_id_list[0]);
    TEST_EQUAL(sm_rows.size(), fm_id_list_size);

    if (sm_rows.size() == fm_id_list_size)
    {
        for (Size i = 0; i < sm_rows.size(); ++i)
        {
            String sm_formula = sm_rows[i].chemical_formula.get();
            //std::cout << "\"" << sm_formula << "\",";
            TEST_STRING_EQUAL(sm_formula, fm_id_list[i]);
        }
    }

    ams_feat_test.setParameters(ams_param);

    // now with isotope filtering
    ams_feat_test.run(exp_fm, test_mztab);

    sms = test_mztab.getSmallMoleculeSectionData();
    sm_rows = sms["AccMassSearch"];

    Size fm_id_filt_list_size = sizeof(fm_id_filt_list)/sizeof(fm_id_filt_list[0]);
    TEST_EQUAL(sm_rows.size(), fm_id_filt_list_size);

    if (sm_rows.size() == fm_id_filt_list_size)
    {
        for (Size i = 0; i < sm_rows.size(); ++i)
        {
            String sm_formula = sm_rows[i].chemical_formula.get();
            // std::cout << "\"" << sm_formula << "\",";
            TEST_STRING_EQUAL(sm_formula, fm_id_filt_list[i]);
        }
    }
}
END_SECTION

ConsensusMap exp_cm;
ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("AccurateMassSearchEngine_input1.consensusXML"), exp_cm);
MzTab test_mztab2;

String cons_id_list[] = {"C27H41NO6","C27H36O6","C27H36O6","C17H22O2","C17H22O2","C17H22O2","C17H22O2","C17H22O2","C17H22O2","C17H22O2","C17H22O2","C26H52NO7P","C26H52NO7P","C27H40O6","C27H40O6","C15H21NO3","C26H46O6","C12H24N2O4","C31H48O2S2","C11H20N4O2","C60H86O18","C48H86O18P2","C48H86O18P2","C48H86O18P2","C48H86O18P2","C26H54NO7P","C26H54NO7P","C24H40N8O4","C27H42O6","C13H23NO3","C10H14N2O6","C10H14N2O6","C10H14N2O6","C9H18N2O4S","C9H18N2O4S"};

START_SECTION((void run(const ConsensusMap &, MzTab &)))
{
    ams_feat_test.run(exp_cm, test_mztab2);

    MzTabSmallMoleculeSectionData sms = test_mztab2.getSmallMoleculeSectionData();
    MzTabSmallMoleculeSectionRows sm_rows = sms["AccMassSearch"];

    Size cons_id_list_size = sizeof(cons_id_list)/sizeof(cons_id_list[0]);

    TEST_EQUAL(sm_rows.size(), cons_id_list_size);

    if (sm_rows.size() == cons_id_list_size)
    {
        for (Size i = 0; i < sm_rows.size(); ++i)
        {
            String sm_formula = sm_rows[i].chemical_formula.get();
            // std::cout << "\"" << sm_formula << "\",";
            TEST_STRING_EQUAL(sm_formula, cons_id_list[i]);
        }
    }
}
END_SECTION


START_SECTION((const String& getInternalIonMode()))
{
  AccurateMassSearchEngine ams;
  Param p;
  p.setValue("ionization_mode","auto");
  ams.setParameters(p);
  TEST_EQUAL(ams.getInternalIonMode(), "auto")
  std::vector<AccurateMassSearchResult> hmdb_results_pos;
  // query a mass in auto mode is invalid (auto mode requires a feature/cf map)
  TEST_EXCEPTION(Exception::InvalidParameter, ams.queryByMass(1234, 1, hmdb_results_pos))

  p.setValue("ionization_mode","negative");
  ams.setParameters(p);
  TEST_EQUAL(ams.getInternalIonMode(), "negative")
  p.setValue("ionization_mode","positive");
  ams.setParameters(p);
  TEST_EQUAL(ams.getInternalIonMode(), "positive")


}
END_SECTION


START_SECTION([EXTRA] template <typename MAPTYPE> void resolveAutoMode_(const MAPTYPE& map))
{

  FeatureMap<> fm_p = exp_fm;
  AccurateMassSearchEngine ams;
  MzTab mzt;
  Param p;
  p.setValue("ionization_mode","auto");
  ams.setParameters(p);


  TEST_EXCEPTION(Exception::InvalidParameter, ams.run(fm_p, mzt)); // has no scan_polarity meta value

  fm_p[0].setMetaValue("scan_polarity", "positive");          // should run ok
  ams.run(fm_p, mzt);

  fm_p[0].setMetaValue("scan_polarity", "negative");          // should run ok
  ams.run(fm_p, mzt);

  fm_p[0].setMetaValue("scan_polarity", "something;somethingelse");
  TEST_EXCEPTION(Exception::InvalidParameter, ams.run(fm_p, mzt)); // not unique

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
