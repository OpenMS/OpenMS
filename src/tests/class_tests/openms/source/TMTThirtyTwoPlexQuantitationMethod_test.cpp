// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow, Julianus Pfeuffer $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/TMTThirtyTwoPlexQuantitationMethod.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Matrix.h>

using namespace OpenMS;
using namespace std;

START_TEST(TMTThirtyTwoPlexQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
TMTThirtyTwoPlexQuantitationMethod* ptr = nullptr;
TMTThirtyTwoPlexQuantitationMethod* null_ptr = nullptr;

START_SECTION(TMTThirtyTwoPlexQuantitationMethod())
{
  ptr = new TMTThirtyTwoPlexQuantitationMethod();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TMTThirtyTwoPlexQuantitationMethod())
{
  delete ptr;
}
END_SECTION

START_SECTION((const String& getMethodName() const ))
{
  TMTThirtyTwoPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getMethodName(), "tmt32plex")
}
END_SECTION

START_SECTION((const IsobaricChannelList& getChannelInformation() const ))
{
  TMTThirtyTwoPlexQuantitationMethod quant_meth;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = quant_meth.getChannelInformation();

  TEST_EQUAL(channel_list.size(), 32)
  ABORT_IF(channel_list.size() != 32)

  String channel_names[32] = {
    "126", "127N", "127C", "127D", "128N", "128C", "128ND", "128CD",
    "129N", "129C", "129ND", "129CD", "130N", "130C", "130ND", "130CD",
    "131N", "131C", "131ND", "131CD", "132N", "132C", "132ND", "132CD",
    "133N", "133C", "133ND", "133CD", "134N", "134ND", "134CD", "135ND"
  };

  double channel_centers[32] = {
    126.127726, 127.124761, 127.131081, 127.134003, 128.128116, 128.134436,
    128.131038, 128.137358, 129.131471, 129.137790, 129.134393, 129.140713,
    130.134825, 130.141145, 130.137748, 130.144068, 131.138180, 131.144500,
    131.141103, 131.147423, 132.141535, 132.147855, 132.144458, 132.150778,
    133.144890, 133.151210, 133.147813, 133.154133, 134.148245, 134.151171,
    134.157491, 135.154526
  };

  // descriptions are empty by default
  for(int i = 0; i < 32; ++i)
  {
    TEST_STRING_EQUAL(channel_list[i].description, "")
    TEST_EQUAL(channel_list[i].name, channel_names[i])
    TEST_EQUAL(channel_list[i].center, channel_centers[i])
    TEST_EQUAL(channel_list[i].id, i)
  }

  for(const auto& channel : channel_list)
  {
    // Each TMT 32-plex channel is affected by 14 other channels due to isotope overlap
    TEST_EQUAL(channel.affected_channels.size(), 14)
  }
}
END_SECTION

START_SECTION((Size getNumberOfChannels() const ))
{
  TMTThirtyTwoPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getNumberOfChannels(), 32)
}
END_SECTION


START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const)){
  TMTThirtyTwoPlexQuantitationMethod quant_meth;

  // The default correction matrix is derived from the Thermo TMTpro deuterated reagent
  // set Certificate of Analysis (product A40000817). The 16 non-deuterated channels match
  // the TMTpro datasheet (identical to TMT 16-/18-plex), so this is not an identity matrix.
  Matrix<double> m = quant_meth.getIsotopeCorrectionMatrix();

  TEST_EQUAL(m.rows(), 32)
  TEST_EQUAL(m.cols(), 32)

  ABORT_IF(m.rows() != 32)
  ABORT_IF(m.cols() != 32)

  // Channel "126" loses 9.74% to isotope impurities (0.31 + 9.09 + 0.02 + 0.32), so its
  // self-contribution (diagonal) is 1 - 0.0974 = 0.9026 and the +13C share (9.09%) leaks
  // into its neighbour 127C (channel index 2).
  TEST_REAL_SIMILAR(m(0,0), 0.9026)
  TEST_REAL_SIMILAR(m(2,0), 0.0909)

  // each channel column conserves intensity (self-contribution + leakages = 100%);
  // for "126" all impurity shares map to valid in-range neighbours
  double col0 = 0.0;
  for (Size i = 0; i < m.rows(); ++i) col0 += m(i, 0);
  TEST_REAL_SIMILAR(col0, 1.0)

  // a deuterated channel (127D, index 3) uses its 14-column row directly: it loses
  // 0.82 + 0.30 + 8.71 + 0.33 + 0.26 = 10.42% to impurities (self = 0.8958) and its
  // +13C share (8.71%) leaks into 128CD (channel index 7).
  TEST_REAL_SIMILAR(m(3,3), 0.8958)
  TEST_REAL_SIMILAR(m(7,3), 0.0871)

  // the default is no longer an identity matrix
  Size off_diagonal_nonzero = 0;
  for (Size i = 0; i < m.rows(); ++i)
    for (Size j = 0; j < m.cols(); ++j)
      if (i != j && m(i, j) > 0.0) ++off_diagonal_nonzero;
  TEST_EQUAL(off_diagonal_nonzero > 0, true)

  // a user-supplied all-NA correction matrix is honored and falls back to the identity
  // (16 non-deuterated rows with 8 columns + 16 deuterated rows with 14 columns)
  Param p = quant_meth.getParameters();
  p.setValue("correction_matrix", std::vector<std::string>(16, "NA/NA/NA/NA/NA/NA/NA/NA"));
  p.setValue("correction_matrix_deuterated", std::vector<std::string>(16, "NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA/NA"));
  quant_meth.setParameters(p);
  Matrix<double> id = quant_meth.getIsotopeCorrectionMatrix();
  for (Size i = 0; i < id.rows(); ++i)
    for (Size j = 0; j < id.cols(); ++j)
      TEST_REAL_SIMILAR(id(i,j), (i == j) ? 1.0 : 0.0)
}
END_SECTION

START_SECTION((Size getReferenceChannel() const ))
{
  TMTThirtyTwoPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getReferenceChannel(), 0)

  Param p;
  p.setValue("reference_channel","128N");
  quant_meth.setParameters(p);

  TEST_EQUAL(quant_meth.getReferenceChannel(), 4)
}
END_SECTION

START_SECTION((TMTThirtyTwoPlexQuantitationMethod(const TMTThirtyTwoPlexQuantitationMethod &other)))
{
  TMTThirtyTwoPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127N_description", "new_description");
  p.setValue("reference_channel", "129C");
  qm.setParameters(p);

  TMTThirtyTwoPlexQuantitationMethod qm2(qm);
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 9)
}
END_SECTION

START_SECTION((TMTThirtyTwoPlexQuantitationMethod& operator=(const TMTThirtyTwoPlexQuantitationMethod &rhs)))
{
  TMTThirtyTwoPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127N_description", "new_description");
  p.setValue("reference_channel", "129C");
  qm.setParameters(p);

  TMTThirtyTwoPlexQuantitationMethod qm2;
  qm2 = qm;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 9)
}
END_SECTION

END_TEST
