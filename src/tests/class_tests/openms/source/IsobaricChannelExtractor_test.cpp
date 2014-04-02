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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(IsobaricChannelExtractor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsobaricChannelExtractor * ptr = 0;
IsobaricChannelExtractor* null_ptr = 0;
IsobaricQuantitationMethod* q_method = new ItraqFourPlexQuantitationMethod();

START_SECTION((IsobaricChannelExtractor(const IsobaricQuantitationMethod * const quant_method)))
{
  ptr = new IsobaricChannelExtractor(q_method);
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IsobaricChannelExtractor())
{
  delete ptr;
}

END_SECTION

START_SECTION((IsobaricChannelExtractor(const IsobaricChannelExtractor &other)))
{
  IsobaricChannelExtractor ice(q_method);
  Param p = ice.getParameters();
  p.setValue("select_activation", "");

  ice.setParameters(p);

  IsobaricChannelExtractor ice2(ice);
  TEST_EQUAL(ice2.getParameters(), p)
}

END_SECTION

START_SECTION((IsobaricChannelExtractor& operator=(const IsobaricChannelExtractor& rhs)))
{
  IsobaricChannelExtractor ice(q_method);
  Param p = ice.getParameters();
  p.setValue("reporter_mass_shift", 0.3);
  ice.setParameters(p);

  IsobaricChannelExtractor ice2(q_method);
  ice2 = ice;
  TEST_EQUAL(ice2.getParameters(), p)
}

END_SECTION

START_SECTION((void extractChannels(const MSExperiment<Peak1D>&ms_exp_data, ConsensusMap & consensus_map)))
{
  {
    // load test data
    MSExperiment<Peak1D> exp;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp);

    // add some more information to the quant method
    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "");
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ice.extractChannels(exp, cm_out);

    // check channel meta information
    TEST_EQUAL(cm_out.getFileDescriptions().size(), 4)
    ABORT_IF(cm_out.getFileDescriptions().size() != 4)

    TEST_EQUAL(cm_out.getFileDescriptions()[0].label, "itraq4plex_114")
    TEST_EQUAL(cm_out.getFileDescriptions()[0].getMetaValue("channel_name"), 114)
    TEST_EQUAL(cm_out.getFileDescriptions()[0].getMetaValue("channel_id"), 0)
    TEST_EQUAL(cm_out.getFileDescriptions()[0].getMetaValue("channel_description"), "ref")
    TEST_EQUAL(cm_out.getFileDescriptions()[0].getMetaValue("channel_center"), 114.1112)

    TEST_EQUAL(cm_out.getFileDescriptions()[1].label, "itraq4plex_115")
    TEST_EQUAL(cm_out.getFileDescriptions()[1].getMetaValue("channel_name"), 115)
    TEST_EQUAL(cm_out.getFileDescriptions()[1].getMetaValue("channel_id"), 1)
    TEST_EQUAL(cm_out.getFileDescriptions()[1].getMetaValue("channel_description"), "something")
    TEST_EQUAL(cm_out.getFileDescriptions()[1].getMetaValue("channel_center"), 115.1082)

    TEST_EQUAL(cm_out.getFileDescriptions()[2].label, "itraq4plex_116")
    TEST_EQUAL(cm_out.getFileDescriptions()[2].getMetaValue("channel_name"), 116)
    TEST_EQUAL(cm_out.getFileDescriptions()[2].getMetaValue("channel_id"), 2)
    TEST_EQUAL(cm_out.getFileDescriptions()[2].getMetaValue("channel_description"), "else")
    TEST_EQUAL(cm_out.getFileDescriptions()[2].getMetaValue("channel_center"), 116.1116)

    TEST_EQUAL(cm_out.getFileDescriptions()[3].label, "itraq4plex_117")
    TEST_EQUAL(cm_out.getFileDescriptions()[3].getMetaValue("channel_name"), 117)
    TEST_EQUAL(cm_out.getFileDescriptions()[3].getMetaValue("channel_id"), 3)
    TEST_EQUAL(cm_out.getFileDescriptions()[3].getMetaValue("channel_description"), "")
    TEST_EQUAL(cm_out.getFileDescriptions()[3].getMetaValue("channel_center"), 117.1149)

    // compare results
    TEST_EQUAL(cm_out.size(), 5)
    ABORT_IF(cm_out.size() != 5)
    ConsensusFeature::iterator cf_it;

    TEST_EQUAL(cm_out[0].size(), 4)
    TEST_EQUAL(cm_out[0].getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=2")
    TEST_REAL_SIMILAR(cm_out[0].getMetaValue("precursor_intensity"), 5251952.5)
    TEST_REAL_SIMILAR(cm_out[0].getMetaValue("precursor_charge"), 2)
    TEST_REAL_SIMILAR(cm_out[0].getIntensity(), 1490501.21)
    cf_it = cm_out[0].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 643005.56)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 458708.97)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 182238.38)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 206543.3)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[0].end())


    TEST_EQUAL(cm_out[1].size(), 4)
    TEST_EQUAL(cm_out[1].getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=4")
    TEST_REAL_SIMILAR(cm_out[1].getMetaValue("precursor_intensity"), 7365030)
    TEST_REAL_SIMILAR(cm_out[1].getMetaValue("precursor_charge"), 3)
    TEST_REAL_SIMILAR(cm_out[1].getIntensity(), 2358063.25)
    cf_it = cm_out[1].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 851248.38)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 875994.77)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 322173.1)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 308647)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[1].end())

    TEST_EQUAL(cm_out[2].size(), 4)
    TEST_EQUAL(cm_out[2].getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=6")
    TEST_REAL_SIMILAR(cm_out[2].getMetaValue("precursor_intensity"), 6835636)
    TEST_REAL_SIMILAR(cm_out[2].getMetaValue("precursor_charge"), 3)
    TEST_REAL_SIMILAR(cm_out[2].getIntensity(), 2623415.33)
    cf_it = cm_out[2].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 898583.7)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 977466.23)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 406220.4)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 341145)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[2].end())

    TEST_EQUAL(cm_out[3].size(), 4)
    TEST_EQUAL(cm_out[3].getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=8")
    TEST_REAL_SIMILAR(cm_out[3].getMetaValue("precursor_intensity"), 6762358)
    TEST_REAL_SIMILAR(cm_out[3].getMetaValue("precursor_charge"), 3)
    TEST_REAL_SIMILAR(cm_out[3].getIntensity(), 1692679.37)
    cf_it = cm_out[3].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 593009)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 661448.27)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 249740.1)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 188482)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[3].end())

    TEST_EQUAL(cm_out[4].size(), 4)
    TEST_EQUAL(cm_out[4].getMetaValue("scan_id"), "controllerType=0 controllerNumber=1 scan=10")
    TEST_REAL_SIMILAR(cm_out[4].getMetaValue("precursor_intensity"), 5464634.5)
    TEST_REAL_SIMILAR(cm_out[4].getMetaValue("precursor_charge"), 2)
    TEST_REAL_SIMILAR(cm_out[4].getIntensity(), 1746368)
    cf_it = cm_out[4].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 648863)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 632090)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 229391)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 236024)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[4].end())
  }
  { // test -> keep_unannotated_precursor
    // load test data
    MSExperiment<Peak1D> exp;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_7.mzML"), exp);

    // add some more information to the quant method
    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "");
    p.setValue("keep_unannotated_precursor", "false");
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ice.extractChannels(exp, cm_out);

    TEST_EQUAL(cm_out.size(), 4)
    ABORT_IF(cm_out.size() != 4)
    TEST_EQUAL(((double)cm_out[0].getMetaValue("precursor_intensity")) != 0.0, true)

    p.setValue("keep_unannotated_precursor", "true");
    ice.setParameters(p);

    ConsensusMap cm_out_w_unannotated;
    ice.extractChannels(exp, cm_out_w_unannotated);

    TEST_EQUAL(cm_out_w_unannotated.size(), 5)
    ABORT_IF(cm_out_w_unannotated.size() != 5)
    TEST_REAL_SIMILAR(cm_out_w_unannotated[0].getMetaValue("precursor_intensity"), 0.0)
  }
  {
    // load test data
    MSExperiment<Peak1D> exp;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp);

    // add some more information to the quant method
    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "");
    p.setValue("min_precursor_intensity", 5300000.0);
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ice.extractChannels(exp, cm_out);

    // compare results
    TEST_EQUAL(cm_out.size(), 4)
    ABORT_IF(cm_out.size() != 4)
    for(ConsensusMap::Iterator cf = cm_out.begin(); cf != cm_out.end(); ++cf)
    {
      double prec_intensity = cf->getMetaValue("precursor_intensity");
      TEST_EQUAL(prec_intensity > 5300000.0, true)
    }
  }
  {
    // load test data
    MSExperiment<Peak1D> exp;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp);

    // add some more information to the quant method
    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "");
    p.setValue("min_reporter_intensity", 200000.0);
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ConsensusFeature::iterator cf_it;
    ice.extractChannels(exp, cm_out);

    TEST_EQUAL(cm_out.size(), 5)
    ABORT_IF(cm_out.size() != 5)

    TEST_EQUAL(cm_out[0].size(), 4)
    cf_it = cm_out[0].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 643005.56)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 458708.97)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0) // is 182238.38 < 200.000
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 206543.3)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[0].end())

    cf_it = cm_out[1].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 851248.38)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 875994.77)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 322173.1)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 308647)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[1].end())

    cf_it = cm_out[2].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 898583.7)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 977466.23)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 406220.4)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 341145)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[2].end())

    cf_it = cm_out[3].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 593009)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 661448.27)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 249740.1)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 0.0) // is 188482 < 200.000
    ++cf_it;
    ABORT_IF(cf_it != cm_out[3].end())

    cf_it = cm_out[4].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 648863)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 632090)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 229391)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 236024)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[4].end())
  }
  {
    // load test data
    MSExperiment<Peak1D> exp;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp);

    // add some more information to the quant method
    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "");
    p.setValue("min_reporter_intensity", 200000.0);
    p.setValue("discard_low_intensity_quantifications", "true");
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ConsensusFeature::iterator cf_it;
    ice.extractChannels(exp, cm_out);

    TEST_EQUAL(cm_out.size(), 3)
    ABORT_IF(cm_out.size() != 3)

    cf_it = cm_out[0].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 851248.38)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 875994.77)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 322173.1)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 308647)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[0].end())

    cf_it = cm_out[1].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 898583.7)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 977466.23)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 406220.4)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 341145)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[1].end())

    cf_it = cm_out[2].begin();
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 648863)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 632090)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 229391)
    ++cf_it;
    TEST_REAL_SIMILAR(cf_it->getIntensity(), 236024)
    ++cf_it;
    ABORT_IF(cf_it != cm_out[2].end())
  }
  {
    // check precursor purity computation
    // - tested purities were validated manually
    // - dataset contains 2 ms1 and 5 ms2 spectra
    //   with the purity values listed below

    MSExperiment<Peak1D> exp_purity;
    MzMLFile mzmlfile;
    mzmlfile.load(OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML"), exp_purity);

    Param pItraq = q_method->getParameters();
    pItraq.setValue("channel_114_description", "ref");
    pItraq.setValue("channel_115_description", "something");
    pItraq.setValue("channel_116_description", "else");
    q_method->setParameters(pItraq);

    IsobaricChannelExtractor ice(q_method);

    // disable activation filtering
    Param p = ice.getParameters();
    p.setValue("select_activation", "");
    ice.setParameters(p);

    // extract channels
    ConsensusMap cm_out;
    ice.extractChannels(exp_purity, cm_out);

    TEST_EQUAL(cm_out.size(), 5)
    ABORT_IF(cm_out.size() != 5)

    // check results
    TEST_REAL_SIMILAR(cm_out[0].getMetaValue("precursor_purity"), 1.0)
    TEST_REAL_SIMILAR(cm_out[1].getMetaValue("precursor_purity"), 0.65472)
    TEST_REAL_SIMILAR(cm_out[2].getMetaValue("precursor_purity"), 0.775739)
    TEST_REAL_SIMILAR(cm_out[3].getMetaValue("precursor_purity"), 0.72009)
    TEST_REAL_SIMILAR(cm_out[4].getMetaValue("precursor_purity"), 1.0)

    // now filter by purity
    p.setValue("min_precursor_purity", 0.75);
    ice.setParameters(p);

    ConsensusMap cm_filtered;
    ice.extractChannels(exp_purity, cm_filtered);

    TEST_EQUAL(cm_filtered.size(), 3)
    ABORT_IF(cm_filtered.size() != 3)

    // check results
    TEST_REAL_SIMILAR(cm_filtered[0].getMetaValue("precursor_purity"), 1.0)
    TEST_REAL_SIMILAR(cm_filtered[1].getMetaValue("precursor_purity"), 0.775739)
    TEST_REAL_SIMILAR(cm_filtered[2].getMetaValue("precursor_purity"), 1.0)
  }
}

END_SECTION

delete q_method;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
