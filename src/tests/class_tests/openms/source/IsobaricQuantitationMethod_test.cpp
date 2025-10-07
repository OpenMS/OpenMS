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
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

class TestQuantitationMethod :
  public IsobaricQuantitationMethod
{
public:
  IsobaricChannelList channel_list;
  String name;
  StringList correction_list;

  TestQuantitationMethod()
  {
    setName("TestQuantitationMethod");
    channel_list.push_back(IsobaricChannelInformation("114", 0, "", 114.1112, {-1, -1, 1, 2}));
    channel_list.push_back(IsobaricChannelInformation("115", 1, "", 115.1082, {-1, 0, 2, 3}));
    channel_list.push_back(IsobaricChannelInformation("116", 2, "", 116.1116, {0, 1, 3, -1}));
    channel_list.push_back(IsobaricChannelInformation("117", 3, "", 117.1149, {1, 2, -1, -1}));
    name = "TestQuantitationMethod";
  }

  ~TestQuantitationMethod() override = default;

  const String& getMethodName() const override
  {
    return name;
  }

  const IsobaricChannelList& getChannelInformation() const override
  {
    return channel_list;
  }

  Size getNumberOfChannels() const override
  {
    return 4;
  }

  Matrix<double> getIsotopeCorrectionMatrix() const override
  {
    return stringListToIsotopeCorrectionMatrix_(correction_list);
  }

  Size getReferenceChannel() const override
  {
    return 0;
  }
};


START_TEST(IsobaricQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IsobaricQuantitationMethod* ptr = nullptr;
IsobaricQuantitationMethod* null_ptr = nullptr;
START_SECTION(IsobaricQuantitationMethod())
{
	ptr = new TestQuantitationMethod();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IsobaricQuantitationMethod())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual const String& getName() const =0))
{
  IsobaricQuantitationMethod* quant_method = new TestQuantitationMethod();
  TEST_STRING_EQUAL(quant_method->getName(), "TestQuantitationMethod")
  delete quant_method;
}
END_SECTION

START_SECTION((virtual const IsobaricChannelList& getChannelInformation() const =0))
{
  IsobaricQuantitationMethod* quant_method = new TestQuantitationMethod();
  IsobaricQuantitationMethod::IsobaricChannelList cl = quant_method->getChannelInformation();
  TEST_EQUAL(cl.size(), 4)
  ABORT_IF(cl.size() != 4)

  TEST_STRING_EQUAL(cl[0].description, "")
  TEST_EQUAL(cl[0].name, 114)
  TEST_EQUAL(cl[0].id, 0)
  TEST_EQUAL(cl[0].center, 114.1112)
  delete quant_method;
}
END_SECTION

START_SECTION((virtual Size getNumberOfChannels() const =0))
{
  IsobaricQuantitationMethod* quant_method = new TestQuantitationMethod();
  TEST_EQUAL(quant_method->getNumberOfChannels(), 4)
  delete quant_method;
}
END_SECTION

START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const =0 with Exception))
{
  auto* quant_method = new TestQuantitationMethod();
  // missing entry
  quant_method->correction_list = ListUtils::create<String>("0.0/1.0/5.9/0.2,0.0/2.0/    0.1,0.0/3.0/4.5/0.1,0.1/4.0/3.5/0.1");
  TEST_EXCEPTION(Exception::InvalidValue, quant_method->getIsotopeCorrectionMatrix())
  delete quant_method;
}
END_SECTION

START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const =0))
{
  auto* quant_method = new TestQuantitationMethod();
  quant_method->correction_list = ListUtils::create<String>("0.0/1.0/5.9/0.2,0.0/2.0/5.6/0.1,0.0/3.0/4.5/0.1,0.1/4.0/3.5/0.1");
  Matrix<double> m = quant_method->getIsotopeCorrectionMatrix();

  ABORT_IF(m.rows() != 4)
  ABORT_IF(m.cols() != 4)

  double real_m[4][4] = {{0.929, 0.02, 0, 0},
    {0.059, 0.923, 0.03, 0.001},
    {0.002, 0.056, 0.924, 0.04},
    {0, 0.001, 0.045, 0.923}};

  for (size_t i = 0; i < m.rows(); ++i)
  {
    for (size_t j = 0; j < m.cols(); ++j)
    {
      TEST_REAL_SIMILAR(real_m[i][j], m(i,j))
    }
  }

  quant_method->correction_list = ListUtils::create<String>("0.0/1.0/10.9/0.2,0.0/2.0/5.6/0.6,0.0/10.0/4.5/0.1,0.1/4.0/3.5/0.1");
  m = quant_method->getIsotopeCorrectionMatrix();

  ABORT_IF(m.rows() != 4)
  ABORT_IF(m.cols() != 4)

  double real_m2[4][4] = {{0.879, 0.02, 0, 0},
    {0.109, 0.918, 0.1, 0.001},
    {0.002, 0.056, 0.854, 0.04},
    {0, 0.006, 0.045, 0.923}};

  for(size_t i = 0; i < m.rows(); ++i)
  {
    for(size_t j = 0; j < m.cols(); ++j)
    {
      TEST_REAL_SIMILAR(real_m2[i][j], m(i,j))
    }
  }
  delete quant_method;
}
END_SECTION

START_SECTION((virtual Size getReferenceChannel() const =0))
{
  IsobaricQuantitationMethod* quant_method = new TestQuantitationMethod();
  TEST_EQUAL(quant_method->getReferenceChannel(), 0)
  delete quant_method;
}
END_SECTION

START_SECTION(([IsobaricQuantitationMethod::IsobaricChannelInformation] IsobaricChannelInformation(const Int name, const Int id, const String &description, const Peak2D::CoordinateType &center)))
{
  IsobaricQuantitationMethod::IsobaricChannelInformation cI(114, 0, "", 114.1112, {-1, -1, -1, -1});
  TEST_STRING_EQUAL(cI.description, "")
  TEST_EQUAL(cI.name, 114)
  TEST_EQUAL(cI.id, 0)
  TEST_EQUAL(cI.center, 114.1112)

  TEST_EQUAL(cI.affected_channels[0], -1)
  TEST_EQUAL(cI.affected_channels[1], -1)
  TEST_EQUAL(cI.affected_channels[2], -1)
  TEST_EQUAL(cI.affected_channels[3], -1)

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
