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
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixteenPlexQuantitationMethod.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Matrix.h>

using namespace OpenMS;
using namespace std;

START_TEST(TMTSixteenPlexQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TMTSixteenPlexQuantitationMethod* ptr = nullptr;
TMTSixteenPlexQuantitationMethod* null_ptr = nullptr;
START_SECTION(TMTSixteenPlexQuantitationMethod())
{
  ptr = new TMTSixteenPlexQuantitationMethod();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TMTSixteenPlexQuantitationMethod())
{
	delete ptr;
}
END_SECTION

START_SECTION((const String& getMethodName() const ))
{
  TMTSixteenPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getMethodName(), "tmt16plex")
}
END_SECTION

START_SECTION((const IsobaricChannelList& getChannelInformation() const ))
{
  TMTSixteenPlexQuantitationMethod quant_meth;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = quant_meth.getChannelInformation();

  TEST_EQUAL(channel_list.size(), 16)
  ABORT_IF(channel_list.size() != 16)

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

  TEST_EQUAL(channel_list[3].name, "128N")
  TEST_EQUAL(channel_list[3].id, 3)
  TEST_EQUAL(channel_list[3].center, 128.128116)

  TEST_EQUAL(channel_list[4].name, "128C")
  TEST_EQUAL(channel_list[4].id, 4)
  TEST_EQUAL(channel_list[4].center, 128.134436)

  TEST_EQUAL(channel_list[5].name, "129N")
  TEST_EQUAL(channel_list[5].id, 5)
  TEST_EQUAL(channel_list[5].center, 129.131471)

  TEST_EQUAL(channel_list[6].name, "129C")
  TEST_EQUAL(channel_list[6].id, 6)
  TEST_EQUAL(channel_list[6].center, 129.137790)

  TEST_EQUAL(channel_list[7].name, "130N")
  TEST_EQUAL(channel_list[7].id, 7)
  TEST_EQUAL(channel_list[7].center, 130.134825)

  TEST_EQUAL(channel_list[8].name, "130C")
  TEST_EQUAL(channel_list[8].id, 8)
  TEST_EQUAL(channel_list[8].center, 130.141145)

  TEST_EQUAL(channel_list[9].name, "131N")
  TEST_EQUAL(channel_list[9].id, 9)
  TEST_EQUAL(channel_list[9].center, 131.138180)

  TEST_EQUAL(channel_list[10].name, "131C")
  TEST_EQUAL(channel_list[10].id, 10)
  TEST_EQUAL(channel_list[10].center, 131.144500)

  TEST_EQUAL(channel_list[11].name, "132N")
  TEST_EQUAL(channel_list[11].id, 11)
  TEST_EQUAL(channel_list[11].center, 132.141535)

  TEST_EQUAL(channel_list[12].name, "132C")
  TEST_EQUAL(channel_list[12].id, 12)
  TEST_EQUAL(channel_list[12].center, 132.147855)

  TEST_EQUAL(channel_list[13].name, "133N")
  TEST_EQUAL(channel_list[13].id, 13)
  TEST_EQUAL(channel_list[13].center, 133.144890)

  TEST_EQUAL(channel_list[14].name, "133C")
  TEST_EQUAL(channel_list[14].id, 14)
  TEST_EQUAL(channel_list[14].center, 133.151210)

  TEST_EQUAL(channel_list[15].name, "134N")
  TEST_EQUAL(channel_list[15].id, 15)
  TEST_EQUAL(channel_list[15].center, 134.148245)

  for (const auto& channel : channel_list)
  {
    TEST_EQUAL(channel.affected_channels.size(), 8)
  }
}
END_SECTION

START_SECTION((Size getNumberOfChannels() const ))
{
  TMTSixteenPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getNumberOfChannels(), 16)
}
END_SECTION

START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const ))
{

  double test_matrix[16][16] = {{0.9026, 0.0078, 0.0093, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
                                 {0.0031, 0.8948, 0, 0.0082, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
                                 {0.0909, 0, 0.8981, 0.0065, 0.0147, 0, 0.0013, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
                                 {0.0002, 0.0941, 0.0035, 0.9014, 0, 0.0146, 0, 0.0013, 0, 0, 0, 0, 0, 0, 0, 0, },
                                 {0.0032, 0, 0.0863, 0, 0.9113, 0.0128, 0.0259, 0, 0.0004, 0, 0, 0, 0, 0, 0, 0, },
                                 {0, 0.0033, 0.0001, 0.0813, 0.0034, 0.9025, 0, 0.0241, 0, 0.0003, 0, 0, 0, 0, 0, 0, },
                                 {0, 0, 0.0027, 0, 0.0691, 0, 0.907, 0.0027, 0.031, 0, 0.0008, 0, 0, 0, 0, 0, },
                                 {0, 0, 0, 0.0026, 0, 0.0686, 0.0032, 0.9151, 0, 0.0278, 0, 0.0015, 0, 0, 0, 0, },
                                 {0, 0, 0, 0, 0.0015, 0, 0.0607, 0, 0.9154, 0.0063, 0.039, 0.0001, 0.0011, 0, 0, 0, },
                                 {0, 0, 0, 0, 0, 0.0015, 0.001, 0.0558, 0.0042, 0.9187, 0, 0.0358, 0, 0.0007, 0, 0, },
                                 {0, 0, 0, 0, 0, 0, 0.0009, 0, 0.0482, 0, 0.9194, 0.0072, 0.0455, 0.0001, 0.0022, 0, },
                                 {0, 0, 0, 0, 0, 0, 0, 0.001, 0.0002, 0.0457, 0.0047, 0.9374, 0, 0.0314, 0, 0.003, },
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0.0006, 0, 0.0357, 0, 0.9305, 0.0073, 0.0496, 0.0003, },
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0012, 0, 0.018, 0.0043, 0.9262, 0, 0.0549, },
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0004, 0, 0.0186, 0, 0.9345, 0.0062, },
                                 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.034, 0.0034, 0.9242, }
  };

  Matrix<double> test_Matrix;
	test_Matrix.setMatrix<double, 16, 16>(test_matrix);

  TMTSixteenPlexQuantitationMethod quant_meth;

  // we only check the default matrix here which is an identity matrix
  // for tmt16plex
  Matrix<double> m = quant_meth.getIsotopeCorrectionMatrix();

  TEST_EQUAL(m.rows(), 16)
  TEST_EQUAL(m.cols(), 16)

  ABORT_IF(m.rows() != 16)
  ABORT_IF(m.cols() != 16)

  for (size_t i = 0; i < m.rows(); ++i)
  {
    for (size_t j = 0; j < m.cols(); ++j)
    {
      TEST_REAL_SIMILAR(m(i,j), test_Matrix(i,j))
    }
  }
}
END_SECTION

START_SECTION((Size getReferenceChannel() const ))
{
  TMTSixteenPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getReferenceChannel(), 0)

  Param p;
  p.setValue("reference_channel","128N");
  quant_meth.setParameters(p);

  TEST_EQUAL(quant_meth.getReferenceChannel(), 3)
}
END_SECTION

START_SECTION((TMTSixteenPlexQuantitationMethod(const TMTSixteenPlexQuantitationMethod &other)))
{
  TMTSixteenPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127N_description", "new_description");
  p.setValue("reference_channel", "129C");
  qm.setParameters(p);

  TMTSixteenPlexQuantitationMethod qm2(qm);
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 6)

}
END_SECTION

START_SECTION((TMTSixteenPlexQuantitationMethod& operator=(const TMTSixteenPlexQuantitationMethod &rhs)))
{
  TMTSixteenPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127N_description", "new_description");
  p.setValue("reference_channel", "130C");
  qm.setParameters(p);

  TMTSixteenPlexQuantitationMethod qm2 = qm;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 8)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
