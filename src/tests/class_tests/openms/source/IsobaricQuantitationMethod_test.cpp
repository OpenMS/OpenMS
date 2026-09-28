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
  StringList correction_list;

  // a generic test stub; not a real method, so it reports MethodType::UNKNOWN to the base c'tor
  TestQuantitationMethod() :
    IsobaricQuantitationMethod(MethodType::UNKNOWN)
  {
    setName("TestQuantitationMethod");
    channel_list.push_back(IsobaricChannelInformation("114", 0, "", 114.1112, {-1, -1, 1, 2}));
    channel_list.push_back(IsobaricChannelInformation("115", 1, "", 115.1082, {-1, 0, 2, 3}));
    channel_list.push_back(IsobaricChannelInformation("116", 2, "", 116.1116, {0, 1, 3, -1}));
    channel_list.push_back(IsobaricChannelInformation("117", 3, "", 117.1149, {1, 2, -1, -1}));
  }

  ~TestQuantitationMethod() override = default;

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
START_SECTION(IsobaricQuantitationMethod(MethodType method_type))
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

START_SECTION((virtual const std::string& getName() const =0))
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
  quant_method->correction_list = ListUtils::create<std::string>("0.0/1.0/5.9/0.2,0.0/2.0/    0.1,0.0/3.0/4.5/0.1,0.1/4.0/3.5/0.1");
  TEST_EXCEPTION(Exception::InvalidValue, quant_method->getIsotopeCorrectionMatrix())
  delete quant_method;
}
END_SECTION

START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const =0))
{
  auto* quant_method = new TestQuantitationMethod();
  quant_method->correction_list = ListUtils::create<std::string>("0.0/1.0/5.9/0.2,0.0/2.0/5.6/0.1,0.0/3.0/4.5/0.1,0.1/4.0/3.5/0.1");
  Matrix<double> m = quant_method->getIsotopeCorrectionMatrix();

  ABORT_IF(m.rows() != 4)
  ABORT_IF(m.cols() != 4)

  double real_m[4][4] = {{0.929, 0.02, 0, 0},
    {0.059, 0.923, 0.03, 0.001},
    {0.002, 0.056, 0.924, 0.04},
    {0, 0.001, 0.045, 0.923}};

  for (Size i = 0; i < m.rows(); ++i)
  {
    for (Size j = 0; j < m.cols(); ++j)
    {
      TEST_REAL_SIMILAR(real_m[i][j], m(i,j))
    }
  }

  quant_method->correction_list = ListUtils::create<std::string>("0.0/1.0/10.9/0.2,0.0/2.0/5.6/0.6,0.0/10.0/4.5/0.1,0.1/4.0/3.5/0.1");
  m = quant_method->getIsotopeCorrectionMatrix();

  ABORT_IF(m.rows() != 4)
  ABORT_IF(m.cols() != 4)

  double real_m2[4][4] = {{0.879, 0.02, 0, 0},
    {0.109, 0.918, 0.1, 0.001},
    {0.002, 0.056, 0.854, 0.04},
    {0, 0.006, 0.045, 0.923}};

  for(Size i = 0; i < m.rows(); ++i)
  {
    for(Size j = 0; j < m.cols(); ++j)
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

START_SECTION(([IsobaricQuantitationMethod::IsobaricChannelInformation] IsobaricChannelInformation(const Int name, const Int id, const std::string &description, const Peak2D::CoordinateType &center)))
{
  IsobaricQuantitationMethod::IsobaricChannelInformation cI(StringUtils::toStr(114), 0, "", 114.1112, {-1, -1, -1, -1});
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

START_SECTION((static std::unique_ptr<IsobaricQuantitationMethod> create(MethodType mt)))
{
  using MT = IsobaricQuantitationMethod::MethodType;

  // UNKNOWN is the disabled/none sentinel and yields a null pointer (no exception).
  TEST_EQUAL(IsobaricQuantitationMethod::create(MT::UNKNOWN) == nullptr, true)

  // Every concrete method must be instantiable and report back the MethodType it was created from.
  for (int i = static_cast<int>(MT::UNKNOWN) + 1; i < static_cast<int>(MT::SIZE_OF_METHODTYPE); ++i)
  {
    const MT mt = static_cast<MT>(i);
    auto method = IsobaricQuantitationMethod::create(mt);
    TEST_EQUAL(method != nullptr, true)
    ABORT_IF(method == nullptr)
    TEST_EQUAL(static_cast<int>(method->getMethodType()), i)
    TEST_STRING_EQUAL(method->getMethodName(), std::string(IsobaricQuantitationMethod::methodTypeName(mt)))
  }

  // The new TMT 32/35-plex methods must instantiate with the correct channel counts (#9460).
  TEST_EQUAL(IsobaricQuantitationMethod::create(MT::TMT_32PLEX)->getNumberOfChannels(), 32)
  TEST_EQUAL(IsobaricQuantitationMethod::create(MT::TMT_35PLEX)->getNumberOfChannels(), 35)

  // The terminator and any out-of-range value are programming errors and must throw.
  TEST_EXCEPTION(Exception::IllegalArgument, IsobaricQuantitationMethod::create(MT::SIZE_OF_METHODTYPE))
  TEST_EXCEPTION(Exception::IllegalArgument, IsobaricQuantitationMethod::create(static_cast<MT>(static_cast<int>(MT::SIZE_OF_METHODTYPE) + 1)))
}
END_SECTION

START_SECTION((const std::string& getMethodName() const))
{
  using MT = IsobaricQuantitationMethod::MethodType;
  auto m = IsobaricQuantitationMethod::create(MT::TMT_6PLEX);
  TEST_STRING_EQUAL(m->getMethodName(), "tmt6plex")
}
END_SECTION

START_SECTION((MethodType getMethodType() const))
{
  using MT = IsobaricQuantitationMethod::MethodType;
  auto m = IsobaricQuantitationMethod::create(MT::ITRAQ_4PLEX);
  TEST_EQUAL(m->getMethodType() == MT::ITRAQ_4PLEX, true)
}
END_SECTION

START_SECTION((static std::string_view methodTypeName(MethodType mt)))
{
  using MT = IsobaricQuantitationMethod::MethodType;
  TEST_STRING_EQUAL(std::string(IsobaricQuantitationMethod::methodTypeName(MT::UNKNOWN)), "unknown")
  TEST_STRING_EQUAL(std::string(IsobaricQuantitationMethod::methodTypeName(MT::TMT_6PLEX)), "tmt6plex")
  TEST_STRING_EQUAL(std::string(IsobaricQuantitationMethod::methodTypeName(MT::ITRAQ_8PLEX)), "itraq8plex")
  TEST_STRING_EQUAL(std::string(IsobaricQuantitationMethod::methodTypeName(MT::TMT_32PLEX)), "tmt32plex")
  TEST_STRING_EQUAL(std::string(IsobaricQuantitationMethod::methodTypeName(MT::TMT_35PLEX)), "tmt35plex")
  // out-of-range values (the SIZE_OF_METHODTYPE terminator and beyond) are programming errors and must throw
  TEST_EXCEPTION(Exception::IllegalArgument, IsobaricQuantitationMethod::methodTypeName(MT::SIZE_OF_METHODTYPE))
  TEST_EXCEPTION(Exception::IllegalArgument, IsobaricQuantitationMethod::methodTypeName(static_cast<MT>(static_cast<int>(MT::SIZE_OF_METHODTYPE) + 1)))
}
END_SECTION

START_SECTION((static std::string_view methodDisplayName(MethodType mt)))
{
  using MT = IsobaricQuantitationMethod::MethodType;
  TEST_STRING_EQUAL(std::string(IsobaricQuantitationMethod::methodDisplayName(MT::UNKNOWN)), "none")
  TEST_STRING_EQUAL(std::string(IsobaricQuantitationMethod::methodDisplayName(MT::TMT_6PLEX)), "TMT 6-plex")
  TEST_STRING_EQUAL(std::string(IsobaricQuantitationMethod::methodDisplayName(MT::TMT_32PLEX)), "TMT 32-plex")
  TEST_STRING_EQUAL(std::string(IsobaricQuantitationMethod::methodDisplayName(MT::TMT_35PLEX)), "TMT 35-plex")
  // out-of-range values (the SIZE_OF_METHODTYPE terminator and beyond) are programming errors and must throw
  TEST_EXCEPTION(Exception::IllegalArgument, IsobaricQuantitationMethod::methodDisplayName(MT::SIZE_OF_METHODTYPE))
  TEST_EXCEPTION(Exception::IllegalArgument, IsobaricQuantitationMethod::methodDisplayName(static_cast<MT>(static_cast<int>(MT::SIZE_OF_METHODTYPE) + 1)))
}
END_SECTION

START_SECTION((static MethodType methodTypeFromName(std::string_view name)))
{
  using MT = IsobaricQuantitationMethod::MethodType;
  TEST_EQUAL(IsobaricQuantitationMethod::methodTypeFromName("tmt6plex") == MT::TMT_6PLEX, true)
  TEST_EQUAL(IsobaricQuantitationMethod::methodTypeFromName("itraq4plex") == MT::ITRAQ_4PLEX, true)
  TEST_EQUAL(IsobaricQuantitationMethod::methodTypeFromName("tmt32plex") == MT::TMT_32PLEX, true)
  TEST_EQUAL(IsobaricQuantitationMethod::methodTypeFromName("tmt35plex") == MT::TMT_35PLEX, true)
  // an unrecognized name returns the UNKNOWN sentinel (string lookup; does not throw)
  TEST_EQUAL(IsobaricQuantitationMethod::methodTypeFromName("not_a_method") == MT::UNKNOWN, true)

  // round-trip: name <-> type is bijective for every concrete method, so a newly added
  // method whose name mapping was forgotten (e.g. TMT 32/35-plex) would fail here.
  for (int i = static_cast<int>(MT::UNKNOWN) + 1; i < static_cast<int>(MT::SIZE_OF_METHODTYPE); ++i)
  {
    const MT mt = static_cast<MT>(i);
    TEST_EQUAL(IsobaricQuantitationMethod::methodTypeFromName(IsobaricQuantitationMethod::methodTypeName(mt)) == mt, true)
  }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
