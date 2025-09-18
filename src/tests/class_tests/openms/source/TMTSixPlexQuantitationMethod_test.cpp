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
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Matrix.h>

using namespace OpenMS;
using namespace std;

START_TEST(TMTSixPlexQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TMTSixPlexQuantitationMethod* ptr = nullptr;
TMTSixPlexQuantitationMethod* null_ptr = nullptr;
START_SECTION(TMTSixPlexQuantitationMethod())
{
	ptr = new TMTSixPlexQuantitationMethod();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TMTSixPlexQuantitationMethod())
{
	delete ptr;
}
END_SECTION

START_SECTION((const String& getMethodName() const ))
{
  TMTSixPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getMethodName(), "tmt6plex")
}
END_SECTION

START_SECTION((const IsobaricChannelList& getChannelInformation() const ))
{
  TMTSixPlexQuantitationMethod quant_meth;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = quant_meth.getChannelInformation();
  
  TEST_EQUAL(channel_list.size(), 6)
  ABORT_IF(channel_list.size() != 6)
  
  // descriptions are empty by default
  TEST_STRING_EQUAL(channel_list[0].description, "")
  TEST_STRING_EQUAL(channel_list[1].description, "")
  TEST_STRING_EQUAL(channel_list[2].description, "")
  TEST_STRING_EQUAL(channel_list[3].description, "")
  TEST_STRING_EQUAL(channel_list[4].description, "")
  TEST_STRING_EQUAL(channel_list[5].description, "")
    
  // check masses&co
  TEST_EQUAL(channel_list[0].name, 126)
  TEST_EQUAL(channel_list[0].id, 0)
  TEST_EQUAL(channel_list[0].center, 126.127725)
        
  TEST_EQUAL(channel_list[1].name, 127)
  TEST_EQUAL(channel_list[1].id, 1)
  TEST_EQUAL(channel_list[1].center, 127.124760)

  TEST_EQUAL(channel_list[2].name, 128)
  TEST_EQUAL(channel_list[2].id, 2)
  TEST_EQUAL(channel_list[2].center, 128.134433)

  TEST_EQUAL(channel_list[3].name, 129)
  TEST_EQUAL(channel_list[3].id, 3)
  TEST_EQUAL(channel_list[3].center, 129.131468)

  TEST_EQUAL(channel_list[4].name, 130)
  TEST_EQUAL(channel_list[4].id, 4)
  TEST_EQUAL(channel_list[4].center, 130.141141)

  TEST_EQUAL(channel_list[5].name, 131)
  TEST_EQUAL(channel_list[5].id, 5)
  TEST_EQUAL(channel_list[5].center, 131.138176)

}
END_SECTION

START_SECTION((Size getNumberOfChannels() const ))
{
  TMTSixPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getNumberOfChannels(), 6)
}
END_SECTION

START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const ))
{
  TMTSixPlexQuantitationMethod quant_meth;
  
  // we only check the default matrix here which is an identity matrix 
  // for tmt6plex
  Matrix<double> m = quant_meth.getIsotopeCorrectionMatrix();
  TEST_EQUAL(m.rows(), 6)
  TEST_EQUAL(m.cols(), 6)
    
  ABORT_IF(m.rows() != 6)
  ABORT_IF(m.cols() != 6)
  
  for(size_t i = 0; i < m.rows(); ++i)
  {
    for(size_t j = 0; j < m.cols(); ++j)
    {
      if (i == j) { TEST_TRUE(m(i,j) > 0.5) } // diagonal entries should be largest
      else { TEST_TRUE(m(i,j) < 0.5) }
    }
  }
}
END_SECTION

START_SECTION((Size getReferenceChannel() const ))
{
  TMTSixPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getReferenceChannel(), 0)
    
  Param p;
  p.setValue("reference_channel",128);
  quant_meth.setParameters(p);
  
  TEST_EQUAL(quant_meth.getReferenceChannel(), 2)
}
END_SECTION

START_SECTION((TMTSixPlexQuantitationMethod(const TMTSixPlexQuantitationMethod &other)))
{
  TMTSixPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127_description", "new_description");
  p.setValue("reference_channel", 129);
  qm.setParameters(p);
  
  TMTSixPlexQuantitationMethod qm2(qm);
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 3)

}
END_SECTION

START_SECTION((TMTSixPlexQuantitationMethod& operator=(const TMTSixPlexQuantitationMethod &rhs)))
{
  TMTSixPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127_description", "new_description");
  p.setValue("reference_channel", 129);
  qm.setParameters(p);
  
  TMTSixPlexQuantitationMethod qm2 = qm;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
