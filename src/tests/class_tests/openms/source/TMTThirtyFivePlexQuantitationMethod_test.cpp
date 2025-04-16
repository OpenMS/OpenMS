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

  String channel_names[35] = {
    "126", "127N", "127C", "127D", "128N", 
    "128C","128ND","128CD","129N", "129C", 
    "129ND","129CD","130N","130C", "130ND",
    "130CD","131N", "131C","131ND","131CD", 
    "132N", "132C", "132ND","132CD","133N",
    "133C", "133ND","133CD","134N", "134C",
    "134ND","134CD","135N", "135ND","135CD"
  };

  double channel_centers[35] = {
    126.127726, 127.124761, 127.131081, 127.134003, 128.128116,
    128.134436, 128.131038, 128.137358, 129.131471, 129.137790,
    129.134393, 129.140713, 130.134825, 130.141145, 130.137748,
    130.144068, 131.138180, 131.144500, 131.141103, 131.147423,
    132.141535, 132.147855, 132.144458, 132.150778, 133.144890,
    133.151210, 133.147813, 133.154133, 134.148245, 134.154566,
    134.151171, 134.157491, 135.151601,135.154526 ,135.160846
   };

   // check that all channels have the correct name and id
   // and that the default descriptions are empty by default

  // descriptions are empty by default
  for(int i = 0; i< 35; i++){
    TEST_STRING_EQUAL(channel_list[i].description, "")
    TEST_EQUAL(channel_list[i].name, channel_names[i])
    TEST_EQUAL(channel_list[i].center, channel_centers[i])
    TEST_EQUAL(channel_list[i].id, i)
  }

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


START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const)){
  double test_matrix[35][35] ={
    {0.9026,0.0078,0.0093,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,},
    {0,0.8948,0,0,0.0082,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0.8981,0,0.0065,0.0147,0,0,0,0.0013,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0.0941,0.0035,0.8958,0,0,0,0,0.0146,0,0,0,0.0013,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0.0863,0,0.9014,0,0,0,0.0128,0.0259,0,0,0,0.0004,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0.0033,0.0001,0,0.0813,0.9113,0,0,0,0,0,0,0.0241,0,0,0,0.0003,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0.0027,0,0,0.0691,0.8915,0,0,0,0,0,0.0027,0.031,0,0,0,0.0008,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0.0026,0,0,0.895,0.0686,0.0032,0,0,0,0,0,0,0.0278,0,0,0,0.0015,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0.0015,0,0,0.9025,0.0607,0,0,0,0,0,0,0.0063,0.039,0,0,0.0001,0.0011,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0.0015,0.907,0,0,0.0558,0.0042,0,0,0,0,0,0,0.0358,0,0,0,0.0007,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0.0009,0.8839,0,0,0.0482,0,0,0,0,0,0,0.0072,0.0455,0,0,0.0001,0.0022,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0.9033,0.001,0.0002,0,0,0.0457,0.0047,0,0,0,0,0,0,0.0314,0,0,0,0.003,0.0003,0,0,0.0019,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0.9151,0.0006,0,0,0,0.0357,0,0,0,0,0,0,0.0073,0.0496,0,0,0.0003,0,0,0,0.0002,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0.9154,0,0,0.0012,0,0,0,0.018,0.0043,0,0,0,0,0,0,0.0549,0.0548,0,0,0.0542,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.8933,0,0,0.0004,0,0,0,0.0186,0,0,0,0,0,0,0.0062,0,0,0,0.0036,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9,0,0,0,0,0,0,0,0,0.034,0.0034,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9187,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9194,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.8955,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.8852,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9374,0,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9305,0,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.8909,0,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9061,0,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9262,0,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9345,0,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9234,0,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9166,0,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9242,0,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9421,0,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.922,0,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.917,0,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9401,0,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9178,0,}, 
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0.9068,}, 
  };

  Matrix<double> test_Matrix;
  test_Matrix.setMatrix<double,35,35>(test_matrix);

  TMTThirtyFivePlexQuantitationMethod quant_meth;

  // for tmt35plex
  Matrix<double> m = quant_meth.getIsotopeCorrectionMatrix();

  TEST_EQUAL(m.rows(), 35)
  TEST_EQUAL(m.cols(), 35)

  ABORT_IF(m.rows() != 35)
  ABORT_IF(m.cols() != 35)

  for(size_t i = 0; i < m.rows(); ++i)
  {
    for(size_t j = 0; j < m.cols(); ++j)
    {
      TEST_REAL_SIMILAR(m(i,j), test_Matrix(i,j))
    }
  }
}
END_SECTION

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