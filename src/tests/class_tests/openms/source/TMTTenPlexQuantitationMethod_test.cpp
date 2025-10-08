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
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Matrix.h>

using namespace OpenMS;
using namespace std;

START_TEST(TMTTenPlexQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TMTTenPlexQuantitationMethod* ptr = nullptr;
TMTTenPlexQuantitationMethod* null_ptr = nullptr;
START_SECTION(TMTTenPlexQuantitationMethod())
{
	ptr = new TMTTenPlexQuantitationMethod();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TMTTenPlexQuantitationMethod())
{
	delete ptr;
}
END_SECTION

START_SECTION((const String& getMethodName() const ))
{
  TMTTenPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getMethodName(), "tmt10plex")
}
END_SECTION

START_SECTION((const IsobaricChannelList& getChannelInformation() const ))
{
  TMTTenPlexQuantitationMethod quant_meth;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = quant_meth.getChannelInformation();

  TEST_EQUAL(channel_list.size(), 10)
  ABORT_IF(channel_list.size() != 10)

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

  // check masses&co
  TEST_EQUAL(channel_list[0].name, "126")
  TEST_EQUAL(channel_list[0].id, 0)
  TEST_EQUAL(channel_list[0].center, 126.127726)
  TEST_EQUAL(channel_list[0].affected_channels[0], -1)
  TEST_EQUAL(channel_list[0].affected_channels[1], -1)
  TEST_EQUAL(channel_list[0].affected_channels[2], 2)
  TEST_EQUAL(channel_list[0].affected_channels[3], 4)

  TEST_EQUAL(channel_list[1].name, "127N")
  TEST_EQUAL(channel_list[1].id, 1)
  TEST_EQUAL(channel_list[1].center, 127.124761)
  TEST_EQUAL(channel_list[1].affected_channels[0], -1)
  TEST_EQUAL(channel_list[1].affected_channels[1], -1)
  TEST_EQUAL(channel_list[1].affected_channels[2], 3)
  TEST_EQUAL(channel_list[1].affected_channels[3], 5)

  TEST_EQUAL(channel_list[2].name, "127C")
  TEST_EQUAL(channel_list[2].id, 2)
  TEST_EQUAL(channel_list[2].center, 127.131081)
  TEST_EQUAL(channel_list[2].affected_channels[0], -1)
  TEST_EQUAL(channel_list[2].affected_channels[1], 0)
  TEST_EQUAL(channel_list[2].affected_channels[2], 4)
  TEST_EQUAL(channel_list[2].affected_channels[3], 6)

  TEST_EQUAL(channel_list[3].name, "128N")
  TEST_EQUAL(channel_list[3].id, 3)
  TEST_EQUAL(channel_list[3].center, 128.128116)
  TEST_EQUAL(channel_list[3].affected_channels[0], -1)
  TEST_EQUAL(channel_list[3].affected_channels[1], 1)
  TEST_EQUAL(channel_list[3].affected_channels[2], 5)
  TEST_EQUAL(channel_list[3].affected_channels[3], 7)

  TEST_EQUAL(channel_list[4].name, "128C")
  TEST_EQUAL(channel_list[4].id, 4)
  TEST_EQUAL(channel_list[4].center, 128.134436)
  TEST_EQUAL(channel_list[4].affected_channels[0], 0)
  TEST_EQUAL(channel_list[4].affected_channels[1], 2)
  TEST_EQUAL(channel_list[4].affected_channels[2], 6)
  TEST_EQUAL(channel_list[4].affected_channels[3], 8)

  TEST_EQUAL(channel_list[5].name, "129N")
  TEST_EQUAL(channel_list[5].id, 5)
  TEST_EQUAL(channel_list[5].center, 129.131471)
  TEST_EQUAL(channel_list[5].affected_channels[0], 1)
  TEST_EQUAL(channel_list[5].affected_channels[1], 3)
  TEST_EQUAL(channel_list[5].affected_channels[2], 7)
  TEST_EQUAL(channel_list[5].affected_channels[3], 9)

  TEST_EQUAL(channel_list[6].name, "129C")
  TEST_EQUAL(channel_list[6].id, 6)
  TEST_EQUAL(channel_list[6].center, 129.137790)
  TEST_EQUAL(channel_list[6].affected_channels[0], 2)
  TEST_EQUAL(channel_list[6].affected_channels[1], 4)
  TEST_EQUAL(channel_list[6].affected_channels[2], 8)
  TEST_EQUAL(channel_list[6].affected_channels[3], -1)

  TEST_EQUAL(channel_list[7].name, "130N")
  TEST_EQUAL(channel_list[7].id, 7)
  TEST_EQUAL(channel_list[7].center, 130.134825)
  TEST_EQUAL(channel_list[7].affected_channels[0], 3)
  TEST_EQUAL(channel_list[7].affected_channels[1], 5)
  TEST_EQUAL(channel_list[7].affected_channels[2], 9)
  TEST_EQUAL(channel_list[7].affected_channels[3], -1)

  TEST_EQUAL(channel_list[8].name, "130C")
  TEST_EQUAL(channel_list[8].id, 8)
  TEST_EQUAL(channel_list[8].center, 130.141145)
  TEST_EQUAL(channel_list[8].affected_channels[0], 4)
  TEST_EQUAL(channel_list[8].affected_channels[1], 6)
  TEST_EQUAL(channel_list[8].affected_channels[2], -1)
  TEST_EQUAL(channel_list[8].affected_channels[3], -1)

  TEST_EQUAL(channel_list[9].name, "131")
  TEST_EQUAL(channel_list[9].id, 9)
  TEST_EQUAL(channel_list[9].center, 131.138180)
  TEST_EQUAL(channel_list[9].affected_channels[0], 5)
  TEST_EQUAL(channel_list[9].affected_channels[1], 7)
  TEST_EQUAL(channel_list[9].affected_channels[2], -1)
  TEST_EQUAL(channel_list[9].affected_channels[3], -1)
}
END_SECTION

START_SECTION((Size getNumberOfChannels() const ))
{
  TMTTenPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getNumberOfChannels(), 10)
}
END_SECTION

START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const ))
{

  double test_matrix[10][10] = {{0.9491,0,0.0037,0,0.0008,0,0,0,0,0},
                                {0,0.9448,0,0.0065,0,0.0001,0,0,0,0},
                                {0.0509,0,0.9412,0,0.0049,0,0,0,0,0},
                                {0,0.0527,0,0.9508,0,0.0071,0,0.0002,0,0},
                                {0,0,0.0536,0,0.9637,0,0.0132,0,0.0003,0},
                                {0,0,0,0.0417,0,0.9621,0,0.0128,0,0.0008},
                                {0,0,0.0015,0,0.0306,0,0.9606,0,0.0208,0},
                                {0,0,0,0.001,0,0.0307,0,0.9342,0,0.0199},
                                {0,0,0,0,0,0,0.0262,0,0.9566,0},
                                {0,0,0,0,0,0,0,0.0275,0,0.9628}};

  Matrix<double> test_Matrix;
	test_Matrix.setMatrix<double,10,10>(test_matrix);

  TMTTenPlexQuantitationMethod quant_meth;

  // we only check the default matrix here which is an identity matrix
  // for tmt10plex
  Matrix<double> m = quant_meth.getIsotopeCorrectionMatrix();
  TEST_EQUAL(m.rows(), 10)
  TEST_EQUAL(m.cols(), 10)

  ABORT_IF(m.rows() != 10)
  ABORT_IF(m.cols() != 10)

  for (size_t i = 0; i < m.rows(); ++i)
  {
    for (size_t j = 0; j < m.cols(); ++j)
    {
      if (i == j) TEST_REAL_SIMILAR(m(i,j), test_Matrix(i,j))
      else TEST_REAL_SIMILAR(m(i,j), test_Matrix(i,j))
    }
  }
}
END_SECTION

START_SECTION((Size getReferenceChannel() const ))
{
  TMTTenPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getReferenceChannel(), 0)

  Param p;
  p.setValue("reference_channel","128N");
  quant_meth.setParameters(p);

  TEST_EQUAL(quant_meth.getReferenceChannel(), 3)
}
END_SECTION

START_SECTION((TMTTenPlexQuantitationMethod(const TMTTenPlexQuantitationMethod &other)))
{
  TMTTenPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127N_description", "new_description");
  p.setValue("reference_channel", "129C");
  qm.setParameters(p);

  TMTTenPlexQuantitationMethod qm2(qm);
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 6)

}
END_SECTION

START_SECTION((TMTTenPlexQuantitationMethod& operator=(const TMTTenPlexQuantitationMethod &rhs)))
{
  TMTTenPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127N_description", "new_description");
  p.setValue("reference_channel", "130C");
  qm.setParameters(p);

  TMTTenPlexQuantitationMethod qm2 = qm;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 8)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
