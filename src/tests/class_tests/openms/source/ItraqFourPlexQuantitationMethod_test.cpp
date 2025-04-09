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
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Matrix.h>

using namespace OpenMS;
using namespace std;

START_TEST(ItraqFourPlexQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ItraqFourPlexQuantitationMethod* ptr = nullptr;
ItraqFourPlexQuantitationMethod* null_ptr = nullptr;
START_SECTION(ItraqFourPlexQuantitationMethod())
{
	ptr = new ItraqFourPlexQuantitationMethod();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ItraqFourPlexQuantitationMethod())
{
	delete ptr;
}
END_SECTION

START_SECTION((const String& getMethodName() const ))
{
  ItraqFourPlexQuantitationMethod quant_meth;
  TEST_STRING_EQUAL(quant_meth.getMethodName(), "itraq4plex")
}
END_SECTION
    
START_SECTION((ItraqFourPlexQuantitationMethod(const ItraqFourPlexQuantitationMethod &other)))
{
  ItraqFourPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_114_description", "new_description");
  p.setValue("reference_channel", 116);
  qm.setParameters(p);
  
  ItraqFourPlexQuantitationMethod qm2(qm);
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[0].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 2)    
}
END_SECTION

START_SECTION((ItraqFourPlexQuantitationMethod& operator=(const ItraqFourPlexQuantitationMethod &rhs)))
{
  ItraqFourPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_114_description", "new_description");
  p.setValue("reference_channel", 116);
  qm.setParameters(p);
  
  ItraqFourPlexQuantitationMethod qm2;
  qm2 = qm;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[0].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 2)
}
END_SECTION


START_SECTION((const IsobaricChannelList& getChannelInformation() const ))
{
  ItraqFourPlexQuantitationMethod quant_meth;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = quant_meth.getChannelInformation();
  
  TEST_EQUAL(channel_list.size(), 4)
  ABORT_IF(channel_list.size() != 4)
  
  // descriptions are empty by default
  TEST_STRING_EQUAL(channel_list[0].description, "")
  TEST_STRING_EQUAL(channel_list[1].description, "")
  TEST_STRING_EQUAL(channel_list[2].description, "")
  TEST_STRING_EQUAL(channel_list[3].description, "")
    
  // check masses&co
  TEST_EQUAL(channel_list[0].name, 114)
  TEST_EQUAL(channel_list[0].id, 0)
  TEST_EQUAL(channel_list[0].center, 114.1112)

  TEST_EQUAL(channel_list[1].name, 115)
  TEST_EQUAL(channel_list[1].id, 1)
  TEST_EQUAL(channel_list[1].center, 115.1082)

  TEST_EQUAL(channel_list[2].name, 116)
  TEST_EQUAL(channel_list[2].id, 2)
  TEST_EQUAL(channel_list[2].center, 116.1116)

  TEST_EQUAL(channel_list[3].name, 117)
  TEST_EQUAL(channel_list[3].id, 3)
  TEST_EQUAL(channel_list[3].center, 117.1149)
}
END_SECTION

START_SECTION((Size getNumberOfChannels() const ))
{
  ItraqFourPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getNumberOfChannels(), 4)
}
END_SECTION

START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const ))
{
  ItraqFourPlexQuantitationMethod quant_meth;
  
  // we only check the default matrix here
  Matrix<double> m = quant_meth.getIsotopeCorrectionMatrix();
  TEST_EQUAL(m.rows(), 4)
  TEST_EQUAL(m.cols(), 4)
    
  ABORT_IF(m.rows() != 4)
  ABORT_IF(m.cols() != 4)

  /**  
    0.929   0.02      0      0 
    0.059  0.923   0.03  0.001 
    0.002  0.056  0.924   0.04 
        0  0.001  0.045  0.923
  */
  double real_m[4][4] = {{0.929, 0.02, 0, 0},
    {0.059, 0.923, 0.03, 0.001},
    {0.002, 0.056, 0.924, 0.04},
    {0, 0.001, 0.045, 0.923}};
  
  for (size_t i = 0; i < m.rows(); ++i)
  {
    for (size_t j = 0; j < m.cols(); ++j)
    {
      TEST_REAL_SIMILAR(m(i,j), real_m[i][j])
    }
  }
}
END_SECTION

START_SECTION((virtual Size getReferenceChannel() const ))
{
  ItraqFourPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getReferenceChannel(), 0)
    
  Param p;
  p.setValue("reference_channel",115);
  quant_meth.setParameters(p);
  
  TEST_EQUAL(quant_meth.getReferenceChannel(), 1)
}
END_SECTION
  
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
