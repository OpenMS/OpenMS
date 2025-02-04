// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>
///////////////////////////
//
#include <OpenMS/DATASTRUCTURES/Matrix.h>

using namespace OpenMS;
using namespace std;

START_TEST(ItraqEightPlexQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ItraqEightPlexQuantitationMethod* ptr = nullptr;
ItraqEightPlexQuantitationMethod* null_ptr = nullptr;
START_SECTION(ItraqEightPlexQuantitationMethod())
{
	ptr = new ItraqEightPlexQuantitationMethod();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ItraqEightPlexQuantitationMethod())
{
	delete ptr;
}
END_SECTION

START_SECTION((const String& getMethodName() const ))
{
  ItraqEightPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getMethodName(), "itraq8plex")
}
END_SECTION

START_SECTION((const IsobaricChannelList& getChannelInformation() const ))
{
  ItraqEightPlexQuantitationMethod quant_meth;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = quant_meth.getChannelInformation();
  
  TEST_EQUAL(channel_list.size(), 8)
  ABORT_IF(channel_list.size() != 8)
  
  // descriptions are empty by default
  TEST_STRING_EQUAL(channel_list[0].description, "")
  TEST_STRING_EQUAL(channel_list[1].description, "")
  TEST_STRING_EQUAL(channel_list[2].description, "")
  TEST_STRING_EQUAL(channel_list[3].description, "")
  TEST_STRING_EQUAL(channel_list[4].description, "")
  TEST_STRING_EQUAL(channel_list[5].description, "")
  TEST_STRING_EQUAL(channel_list[6].description, "")
  TEST_STRING_EQUAL(channel_list[7].description, "")              
    
  // check masses&co
  TEST_EQUAL(channel_list[0].name, 113)
  TEST_EQUAL(channel_list[0].id, 0)
  TEST_EQUAL(channel_list[0].center, 113.1078)
        
  TEST_EQUAL(channel_list[1].name, 114)
  TEST_EQUAL(channel_list[1].id, 1)
  TEST_EQUAL(channel_list[1].center, 114.1112)

  TEST_EQUAL(channel_list[2].name, 115)
  TEST_EQUAL(channel_list[2].id, 2)
  TEST_EQUAL(channel_list[2].center, 115.1082)

  TEST_EQUAL(channel_list[3].name, 116)
  TEST_EQUAL(channel_list[3].id, 3)
  TEST_EQUAL(channel_list[3].center, 116.1116)

  TEST_EQUAL(channel_list[4].name, 117)
  TEST_EQUAL(channel_list[4].id, 4)
  TEST_EQUAL(channel_list[4].center, 117.1149)

  TEST_EQUAL(channel_list[5].name, 118)
  TEST_EQUAL(channel_list[5].id, 5)
  TEST_EQUAL(channel_list[5].center, 118.1120)

  TEST_EQUAL(channel_list[6].name, 119)
  TEST_EQUAL(channel_list[6].id, 6)
  TEST_EQUAL(channel_list[6].center, 119.1153)

  TEST_EQUAL(channel_list[7].name, 121)
  TEST_EQUAL(channel_list[7].id, 7)
  TEST_EQUAL(channel_list[7].center, 121.1220)
}
END_SECTION

START_SECTION((Size getNumberOfChannels() const ))
{
  ItraqEightPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getNumberOfChannels(), 8)
}
END_SECTION

START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const ))
{
  ItraqEightPlexQuantitationMethod quant_meth;
  
  // we only check the default matrix here
  Matrix<double> m = quant_meth.getIsotopeCorrectionMatrix();
  TEST_EQUAL(m.rows(), 8)
  TEST_EQUAL(m.cols(), 8)
    
  ABORT_IF(m.rows() != 8)
  ABORT_IF(m.cols() != 8)
  
  double real_m[8][8] = {
    { 0.9289, 0.0094, 0, 0, 0, 0, 0, 0 },
    { 0.0689, 0.93, 0.0188, 0, 0, 0, 0, 0 },  
    { 0.0022, 0.059, 0.9312, 0.0282, 0.0006, 0, 0, 0 }, 
    { 0, 0.0016, 0.049, 0.9321, 0.0377, 0.0009, 0, 0 }, 
    { 0, 0, 0.001, 0.039, 0.9318, 0.0471, 0.0014, 0 }, 
    { 0, 0, 0, 0.0007, 0.0299, 0.9332, 0.0566, 0 }, 
    { 0, 0, 0, 0, 0, 0.0188, 0.9333, 0.0027 }, 
    { 0, 0, 0, 0, 0, 0, 0, 0.9211 }
  };
  
  for(size_t i = 0; i < m.rows(); ++i)
  {
    for(size_t j = 0; j < m.cols(); ++j)
    {
      TEST_REAL_SIMILAR(m(i,j), real_m[i][j])
    }
  }
}
END_SECTION

START_SECTION((Size getReferenceChannel() const ))
{
  ItraqEightPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getReferenceChannel(), 0)
    
  Param p;
  p.setValue("reference_channel",115);
  quant_meth.setParameters(p);
  
  TEST_EQUAL(quant_meth.getReferenceChannel(), 2)
    
  p.setValue("reference_channel",121);
  quant_meth.setParameters(p);

  TEST_EQUAL(quant_meth.getReferenceChannel(), 7)
}
END_SECTION

START_SECTION((ItraqEightPlexQuantitationMethod(const ItraqEightPlexQuantitationMethod &other)))
{
  ItraqEightPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_114_description", "new_description");
  p.setValue("reference_channel", 116);
  qm.setParameters(p);
  
  ItraqEightPlexQuantitationMethod qm2(qm);
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 3)
}
END_SECTION

START_SECTION((ItraqEightPlexQuantitationMethod& operator=(const ItraqEightPlexQuantitationMethod &rhs)))
{
  ItraqEightPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_114_description", "new_description");
  p.setValue("reference_channel", 116);
  qm.setParameters(p);
  
  ItraqEightPlexQuantitationMethod qm2 = qm;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
