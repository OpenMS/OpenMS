// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $, Julianus Pfeuffer $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/TMTThirtyFivePlexQuantitationMethod.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Matrix.h>

using namespace OpenMS;
using namespace std;

START_TEST(TMTThirtyFivePlexQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
TMTThirtyFivePlexQuantitationMethod* ptr = nullptr;
TMTThirtyFivePlexQuantitationMethod* null_ptr = nullptr;

START_SECTION(TMTThirtyFivePlexQuantitationMethod())
{
  ptr = new TMTThirtyFivePlexQuantitationMethod();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TMTThirtyFivePlexQuantitationMethod())
{
  delete ptr;
}
END_SECTION

START_SECTION((const String& getMethodName() const ))
{
  TMTThirtyFivePlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getMethodName(), "tmt35plex")
}
END_SECTION

START_SECTION((const IsobaricChannelList& getChannelInformation() const ))
{
  TMTThirtyFivePlexQuantitationMethod quant_meth;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = quant_meth.getChannelInformation();

  TEST_EQUAL(channel_list.size(), 35)
  ABORT_IF(channel_list.size() != 35)

  // descriptions are empty by default
  TEST_STRING_EQUAL(channel_list[0].description, "")
  TEST_STRING_EQUAL(channel_list[1].description, "")
  TEST_STRING_EQUAL(channel_list[2].description, "")
  TEST_STRING_EQUAL(channel_list[3].description, "")
  TEST_STRING_EQUAL(channel_list[4].description, "")
  TEST_STRING_EQUAL(channel_list[5].description, "")
  TEST_STRING_EQUAL(channel_list[6].description, "")
  TEST_STRING_EQUAL(channel_list[7].description, "")
  TEST_STRING_EQUAL(channel_list[8].description, "")
  TEST_STRING_EQUAL(channel_list[9].description, "")
  TEST_STRING_EQUAL(channel_list[10].description, "")
  TEST_STRING_EQUAL(channel_list[11].description, "")
  TEST_STRING_EQUAL(channel_list[12].description, "")
  TEST_STRING_EQUAL(channel_list[13].description, "")
  TEST_STRING_EQUAL(channel_list[14].description, "")
  TEST_STRING_EQUAL(channel_list[15].description, "")
  TEST_STRING_EQUAL(channel_list[16].description, "")
  TEST_STRING_EQUAL(channel_list[17].description, "")
  TEST_STRING_EQUAL(channel_list[18].description, "")
  TEST_STRING_EQUAL(channel_list[19].description, "")
  TEST_STRING_EQUAL(channel_list[20].description, "")
  TEST_STRING_EQUAL(channel_list[21].description, "")
  TEST_STRING_EQUAL(channel_list[22].description, "")
  TEST_STRING_EQUAL(channel_list[23].description, "")
  TEST_STRING_EQUAL(channel_list[24].description, "")
  TEST_STRING_EQUAL(channel_list[25].description, "")
  TEST_STRING_EQUAL(channel_list[26].description, "")
  TEST_STRING_EQUAL(channel_list[27].description, "")
  TEST_STRING_EQUAL(channel_list[28].description, "")
  TEST_STRING_EQUAL(channel_list[29].description, "")
  TEST_STRING_EQUAL(channel_list[30].description, "")
  TEST_STRING_EQUAL(channel_list[31].description, "")
  TEST_STRING_EQUAL(channel_list[32].description, "")
  TEST_STRING_EQUAL(channel_list[33].description, "")
  TEST_STRING_EQUAL(channel_list[34].description, "")

  // check masses&co
  TEST_EQUAL(channel_list[0].name, "126")
  TEST_EQUAL(channel_list[0].id, 0)
  TEST_EQUAL(channel_list[0].center, 126.127726)

  TEST_EQUAL(channel_list[1].name, "127N")
  TEST_EQUAL(channel_list[1].id, 1)
  TEST_EQUAL(channel_list[1].center, 127.124761)

  TEST_EQUAL(channel_list[2].name, "127C")
  TEST_EQUAL(channel_list[2].id, 2)
  TEST_EQUAL(channel_list[2].center, 127.131081)

  TEST_EQUAL(channel_list[3].name, "127D")
  TEST_EQUAL(channel_list[3].id, 3)
  TEST_EQUAL(channel_list[3].center, 127.134003)

  TEST_EQUAL(channel_list[4].name, "128N")
  TEST_EQUAL(channel_list[4].id, 4)
  TEST_EQUAL(channel_list[4].center, 128.128116)

  TEST_EQUAL(channel_list[5].name, "128C")
  TEST_EQUAL(channel_list[5].id, 5)
  TEST_EQUAL(channel_list[5].center, 128.134436)
  
  TEST_EQUAL(channel_list[6].name, "128ND")
  TEST_EQUAL(channel_list[6].id, 6)
  TEST_EQUAL(channel_list[6].center, 128.131038)

  TEST_EQUAL(channel_list[7].name, "128CD")
  TEST_EQUAL(channel_list[7].id, 7)
  TEST_EQUAL(channel_list[7].center, 128.137358)

  TEST_EQUAL(channel_list[8].name, "129N")
  TEST_EQUAL(channel_list[8].id, 8)
  TEST_EQUAL(channel_list[8].center, 129.131471)

  TEST_EQUAL(channel_list[9].name, "129C")
  TEST_EQUAL(channel_list[9].id, 9)
  TEST_EQUAL(channel_list[9].center, 129.137790)

  TEST_EQUAL(channel_list[10].name, "129ND")
  TEST_EQUAL(channel_list[10].id, 10)
  TEST_EQUAL(channel_list[10].center, 129.134393)

  TEST_EQUAL(channel_list[11].name, "129CD")
  TEST_EQUAL(channel_list[11].id, 11)
  TEST_EQUAL(channel_list[11].center, 129.140713)

  TEST_EQUAL(channel_list[12].name, "130N")
  TEST_EQUAL(channel_list[12].id, 12)
  TEST_EQUAL(channel_list[12].center, 130.134825)

  TEST_EQUAL(channel_list[13].name, "130C")
  TEST_EQUAL(channel_list[13].id, 13)
  TEST_EQUAL(channel_list[13].center, 130.141145)
  
  TEST_EQUAL(channel_list[14].name, "130ND")
  TEST_EQUAL(channel_list[14].id, 14)
  TEST_EQUAL(channel_list[14].center, 130.137748)

  TEST_EQUAL(channel_list[15].name, "130CD")
  TEST_EQUAL(channel_list[15].id, 15)
  TEST_EQUAL(channel_list[15].center, 130.144068)

  TEST_EQUAL(channel_list[16].name, "131N")
  TEST_EQUAL(channel_list[16].id, 16)
  TEST_EQUAL(channel_list[16].center, 131.138180)

  TEST_EQUAL(channel_list[17].name, "131C")
  TEST_EQUAL(channel_list[17].id, 17)
  TEST_EQUAL(channel_list[17].center, 131.144500)

  TEST_EQUAL(channel_list[18].name, "131ND")
  TEST_EQUAL(channel_list[18].id, 18)
  TEST_EQUAL(channel_list[18].center, 131.141103)

  TEST_EQUAL(channel_list[19].name, "131CD")
  TEST_EQUAL(channel_list[19].id, 19)
  TEST_EQUAL(channel_list[19].center, 131.147423)

  TEST_EQUAL(channel_list[20].name, "132N")
  TEST_EQUAL(channel_list[20].id, 20)
  TEST_EQUAL(channel_list[20].center, 132.141535)

  TEST_EQUAL(channel_list[21].name, "132C")
  TEST_EQUAL(channel_list[21].id, 21)
  TEST_EQUAL(channel_list[21].center, 132.147855)
  
  TEST_EQUAL(channel_list[22].name, "132ND")
  TEST_EQUAL(channel_list[22].id, 22)
  TEST_EQUAL(channel_list[22].center, 132.144458)

  TEST_EQUAL(channel_list[23].name, "132CD")
  TEST_EQUAL(channel_list[23].id, 23)
  TEST_EQUAL(channel_list[23].center, 132.150778)

  TEST_EQUAL(channel_list[24].name, "133N")
  TEST_EQUAL(channel_list[24].id, 24)
  TEST_EQUAL(channel_list[24].center, 133.144890)

  TEST_EQUAL(channel_list[25].name, "133C")
  TEST_EQUAL(channel_list[25].id, 25)
  TEST_EQUAL(channel_list[25].center, 133.151210)

  TEST_EQUAL(channel_list[26].name, "133ND")
  TEST_EQUAL(channel_list[26].id, 26)
  TEST_EQUAL(channel_list[26].center, 133.147813)

  TEST_EQUAL(channel_list[27].name, "133CD")
  TEST_EQUAL(channel_list[27].id, 27)
  TEST_EQUAL(channel_list[27].center, 133.154133)

  TEST_EQUAL(channel_list[28].name, "134N")
  TEST_EQUAL(channel_list[28].id, 28)
  TEST_EQUAL(channel_list[28].center, 134.148245)

  TEST_EQUAL(channel_list[29].name, "134C")
  TEST_EQUAL(channel_list[29].id, 29)
  TEST_EQUAL(channel_list[29].center, 134.154566)

  TEST_EQUAL(channel_list[30].name, "134ND")
  TEST_EQUAL(channel_list[30].id, 30)
  TEST_EQUAL(channel_list[30].center, 134.151171)
  
  TEST_EQUAL(channel_list[31].name, "134CD")
  TEST_EQUAL(channel_list[31].id, 31)
  TEST_EQUAL(channel_list[31].center, 134.157491)

  TEST_EQUAL(channel_list[32].name, "135N")
  TEST_EQUAL(channel_list[32].id, 32)
  TEST_EQUAL(channel_list[32].center, 135.151601)

  TEST_EQUAL(channel_list[33].name, "135ND")
  TEST_EQUAL(channel_list[33].id, 33)
  TEST_EQUAL(channel_list[33].center, 135.154526)

  TEST_EQUAL(channel_list[34].name, "135CD")
  TEST_EQUAL(channel_list[34].id, 34)
  TEST_EQUAL(channel_list[34].center, 135.160846)

  for(const auto& channel : channel_list)
  {
    TEST_EQUAL(channel.affected_channels.size(), 14)
  }
}
END_SECTION

START_SECTION((Size getNumberOfChannels() const ))
{
  TMTThirtyFivePlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getNumberOfChannels(), 35)
}
END_SECTION

// the test_matrix has to be defined

// START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const)){

//   double test_matrix[32][32] ={};

//   Matrix<double> test_Matrix;
//   test_Matrix.setMatrix<double,32,32>(test_matrix);

//   TMTThirtyFivePlexQuantitationMethod quant_meth;

//   // we only check the default matrix here which is an identity matrix
//   // for tmt32plex
//   Matrix<double> m = quant_meth.getIsotopeCorrectionMatrix();

//   TEST_EQUAL(m.rows(), 32)
//   TEST_EQUAL(m.cols(), 32)

//   ABORT_IF(m.rows() != 32)
//   ABORT_IF(m.cols() != 32)

//   for(size_t i = 0; i < m.rows(); ++i)
//   {
//     for(size_t j = 0; j < m.cols(); ++j)
//     {
//       TEST_REAL_SIMILAR(m(i,j), test_Matrix(i,j))
//     }
//   }
// }
// END_SECTION

START_SECTION((Size getReferenceChannel() const ))
{
  TMTThirtyFivePlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getReferenceChannel(), 0)

  Param p;
  p.setValue("reference_channel","128N");
  quant_meth.setParameters(p);

  TEST_EQUAL(quant_meth.getReferenceChannel(), 4)
}
END_SECTION

START_SECTION((TMTThirtyFivePlexQuantitationMethod(const TMTThirtyFivePlexQuantitationMethod &other)))
{
  TMTThirtyFivePlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127N_description", "new_description");
  p.setValue("reference_channel", "129C");
  qm.setParameters(p);

  TMTThirtyFivePlexQuantitationMethod qm2(qm);
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 9)
}
END_SECTION

START_SECTION((TMTThirtyFivePlexQuantitationMethod& operator=(const TMTThirtyFivePlexQuantitationMethod &rhs)))
{
  TMTThirtyFivePlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127N_description", "new_description");
  p.setValue("reference_channel", "129C");
  qm.setParameters(p);

  TMTThirtyFivePlexQuantitationMethod qm2;
  qm2 = qm;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 9)
}
END_SECTION

END_TEST