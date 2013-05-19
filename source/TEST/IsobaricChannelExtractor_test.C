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

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
///////////////////////////

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(IsobaricChannelExtractor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsobaricChannelExtractor* ptr = 0;
IsobaricChannelExtractor* null_ptr = 0;
IsobaricQuantitationMethod* q_method = new ItraqFourPlexQuantitationMethod();

START_SECTION((IsobaricChannelExtractor(const IsobaricQuantitationMethod *const quant_method)))
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

START_SECTION((IsobaricChannelExtractor& operator=(const IsobaricChannelExtractor &rhs)))
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

START_SECTION((void extractChannels(const MSExperiment< Peak1D > &ms_exp_data, ConsensusMap &consensus_map)))
{
  // load test data
  MzDataFile mz_data_file;
  MSExperiment<Peak1D > exp;
  mz_data_file.load(OPENMS_GET_TEST_DATA_PATH("ItraqChannelExtractor.mzData"),exp);

  // add some more information to the quant method
  Param pItraq = q_method->getParameters();
  pItraq.setValue("channel_114_description", "ref");
  pItraq.setValue("channel_115_description", "something");
  pItraq.setValue("channel_116_description", "else");
  q_method->setParameters(pItraq);

  IsobaricChannelExtractor ice(q_method);

  // disable activation filtering
  Param p = ice.getParameters();
  p.setValue("select_activation","");
  ice.setParameters(p);

  // extract channels
  ConsensusMap cm_out;
  ice.extractChannels(exp, cm_out);
  
  // compare results
  ConsensusXMLFile cm_file;
  String cm_file_out;
  NEW_TMP_FILE(cm_file_out);
  cm_file.store(cm_file_out,cm_out);
  WHITELIST("<?xml-stylesheet");
  TEST_FILE_SIMILAR(cm_file_out,OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor.consensusXML"));
}
END_SECTION

delete q_method;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



