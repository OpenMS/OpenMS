// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/BaseFeature.h>
///////////////////////////

#include <OpenMS/METADATA/PeptideIdentification.h>

START_TEST(BaseFeature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

typedef BaseFeature::QualityType QualityType;
typedef BaseFeature::WidthType WidthType;
BaseFeature* feat_ptr = nullptr;
BaseFeature* feat_nullPointer = nullptr;

START_SECTION((BaseFeature()))
{
  feat_ptr = new BaseFeature;
  TEST_NOT_EQUAL(feat_ptr, feat_nullPointer);
}
END_SECTION

START_SECTION((~BaseFeature()))
{
  delete feat_ptr;
}
END_SECTION

START_SECTION((QualityType getQuality() const))
  BaseFeature p;
  TEST_REAL_SIMILAR(p.getQuality(), 0.0)
  // continued in "setQuality" test
END_SECTION

START_SECTION((void setQuality(QualityType q)))
  BaseFeature p;
  p.setQuality((QualityType) 123.456);
  TEST_REAL_SIMILAR(p.getQuality(), 123.456)
  p.setQuality((QualityType)-0.12345);
  TEST_REAL_SIMILAR(p.getQuality(), -0.12345)
  p.setQuality((QualityType)0.0);
  TEST_REAL_SIMILAR(p.getQuality(), 0.0)
END_SECTION

START_SECTION((WidthType getWidth() const))
  BaseFeature p;
  TEST_REAL_SIMILAR(p.getWidth(), 0.0)
  // continued in "setWidth" test
END_SECTION

START_SECTION((void setWidth(WidthType fwhm)))
  BaseFeature p;
  p.setWidth((WidthType) 123.456);
  TEST_REAL_SIMILAR(p.getWidth(), (WidthType) 123.456)
  p.setWidth((WidthType) -0.12345);
  TEST_REAL_SIMILAR(p.getWidth(), (WidthType) -0.12345)
  p.setWidth((WidthType) 0.0);
  TEST_REAL_SIMILAR(p.getWidth(), (WidthType) 0.0)
END_SECTION

START_SECTION([EXTRA](IntensityType getIntensity() const))
  const BaseFeature p;
  TEST_REAL_SIMILAR(p.getIntensity(), 0.0)
END_SECTION

START_SECTION([EXTRA](const PositionType& getPosition() const))
  const BaseFeature	p;
  TEST_REAL_SIMILAR(p.getPosition()[0], 0.0)
  TEST_REAL_SIMILAR(p.getPosition()[1], 0.0)
END_SECTION

START_SECTION([EXTRA](IntensityType& getIntensity()))
  BaseFeature p;
  TEST_REAL_SIMILAR(p.getIntensity(), 0.0f)
  p.setIntensity(123.456f);
  TEST_REAL_SIMILAR(p.getIntensity(), 123.456f)
  p.setIntensity(-0.12345f);
  TEST_REAL_SIMILAR(p.getIntensity(), -0.12345f)
  p.setIntensity(0.0f);
  TEST_REAL_SIMILAR(p.getIntensity(), 0.0f)
END_SECTION

START_SECTION([EXTRA](PositionType& getPosition()))
  BaseFeature::PositionType pos;
  BaseFeature p;
  pos = p.getPosition();
  TEST_REAL_SIMILAR(pos[0], 0.0)
  TEST_REAL_SIMILAR(pos[1], 0.0)
  pos[0] = 1.0;
  pos[1] = 2.0;
  p.setPosition(pos);
  BaseFeature::PositionType pos2(p.getPosition());
  TEST_REAL_SIMILAR(pos2[0], 1.0)
  TEST_REAL_SIMILAR(pos2[1], 2.0)
END_SECTION

START_SECTION((const ChargeType& getCharge() const))
{
  BaseFeature const tmp;
  TEST_EQUAL(tmp.getCharge(),0);
  // continued in "setCharge" test
}
END_SECTION

START_SECTION((void setCharge(const ChargeType &ch)))
{
  BaseFeature tmp;
  TEST_EQUAL(tmp.getCharge(),0);
  tmp.setCharge(17);
  TEST_EQUAL(tmp.getCharge(),17);
}
END_SECTION

START_SECTION((BaseFeature(const BaseFeature &feature)))
{
  BaseFeature::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  BaseFeature p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  p.setMetaValue("cluster_id", 4711);
  p.setQuality((QualityType)0.9);

  BaseFeature copy_of_p(p);
  BaseFeature::PositionType pos2 = copy_of_p.getPosition();
  BaseFeature::IntensityType i2 = copy_of_p.getIntensity();
  BaseFeature::QualityType q2 = copy_of_p.getQuality();

  TEST_REAL_SIMILAR(i2, 123.456)
  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
  TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"), DataValue(4711));
  TEST_REAL_SIMILAR(q2, 0.9)
}
END_SECTION

START_SECTION((BaseFeature(BaseFeature &&feature)))
{
  // Ensure that BaseFeature has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(BaseFeature(std::declval<BaseFeature&&>())), true)

  BaseFeature::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  BaseFeature p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  p.setMetaValue("cluster_id",4711);
  p.setQuality((QualityType)0.9);

  BaseFeature orig = p;
  BaseFeature copy_of_p(std::move(p));

  BaseFeature::PositionType pos2 = copy_of_p.getPosition();
  BaseFeature::IntensityType i2 = copy_of_p.getIntensity();
  BaseFeature::QualityType q2 = copy_of_p.getQuality();

  TEST_REAL_SIMILAR(i2, 123.456)
  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
  TEST_EQUAL(copy_of_p.getMetaValue("cluster_id"), DataValue(4711));
  TEST_REAL_SIMILAR(q2, 0.9)
}
END_SECTION

START_SECTION((BaseFeature(const Peak2D& point)))
{
  Peak2D point;
  point.setRT(1.23);
  point.setMZ(4.56);
  point.setIntensity(OpenMS::Peak2D::IntensityType(7.89));

  BaseFeature copy(point);
  TEST_REAL_SIMILAR(copy.getRT(), 1.23);
  TEST_REAL_SIMILAR(copy.getMZ(), 4.56);
  TEST_REAL_SIMILAR(copy.getIntensity(), 7.89);
  TEST_EQUAL(copy.getQuality(), 0.0);
  TEST_EQUAL(copy.getCharge(), 0);
  TEST_EQUAL(copy.getWidth(), 0.0);
  TEST_EQUAL(copy.getIDMatches().empty(), true);
}
END_SECTION

START_SECTION((BaseFeature(const RichPeak2D& point)))
{
  RichPeak2D point;
  point.setRT(1.23);
  point.setMZ(4.56);
  point.setIntensity(OpenMS::Peak2D::IntensityType(7.89));
  point.setMetaValue("meta", "test");

  BaseFeature copy(point);
  TEST_REAL_SIMILAR(copy.getRT(), 1.23);
  TEST_REAL_SIMILAR(copy.getMZ(), 4.56);
  TEST_REAL_SIMILAR(copy.getIntensity(), 7.89);
  TEST_EQUAL(copy.getMetaValue("meta"), "test");
  TEST_EQUAL(copy.getQuality(), 0.0);
  TEST_EQUAL(copy.getCharge(), 0);
  TEST_EQUAL(copy.getWidth(), 0.0);
  TEST_EQUAL(copy.getIDMatches().empty(), true);
}
END_SECTION

START_SECTION((BaseFeature& operator=(const BaseFeature& rhs)))
  BaseFeature::PositionType pos;
  pos[0] = 21.21;
  pos[1] = 22.22;
  BaseFeature p;
  p.setIntensity(123.456f);
  p.setPosition(pos);
  p.setQuality((QualityType)0.9);

  BaseFeature copy_of_p;
  copy_of_p = p;

  BaseFeature::PositionType pos2 = copy_of_p.getPosition();
  BaseFeature::IntensityType i2 = copy_of_p.getIntensity();
  BaseFeature::QualityType q2 =  copy_of_p.getQuality();

  TEST_REAL_SIMILAR(i2, 123.456)
  TEST_REAL_SIMILAR(pos2[0], 21.21)
  TEST_REAL_SIMILAR(pos2[1], 22.22)
  TEST_REAL_SIMILAR(q2, 0.9)
END_SECTION

START_SECTION((bool operator==(const BaseFeature &rhs) const))
{
  BaseFeature p1;
  BaseFeature p2(p1);
  TEST_TRUE(p1 == p2)

  p1.setIntensity(5.0f);
  p1.setQuality((QualityType)0.9);
  TEST_EQUAL(p1 == p2, false)
  p2.setIntensity(5.0f);
  p2.setQuality((QualityType)0.9);
  TEST_TRUE(p1 == p2)

  p1.getPosition()[0] = 5;
  TEST_EQUAL(p1 == p2, false)
  p2.getPosition()[0] = 5;
  TEST_TRUE(p1 == p2)

  p1.addIDMatch({"run", IdentificationData::MatchId{1}});
  TEST_EQUAL(p1 == p2, false);
  p2.addIDMatch({"run", IdentificationData::MatchId{1}});
  TEST_TRUE(p1 == p2);
}
END_SECTION

START_SECTION((bool operator!=(const BaseFeature& rhs) const))
  BaseFeature p1;
  BaseFeature p2(p1);
  TEST_EQUAL(p1 != p2, false)

  p1.setIntensity(5.0f);
  TEST_FALSE(p1 == p2)
  p2.setIntensity(5.0f);
  TEST_EQUAL(p1 != p2, false)

  p1.getPosition()[0] = 5;
  TEST_FALSE(p1 == p2)
  p2.getPosition()[0] = 5;
  TEST_EQUAL(p1 != p2, false)

  p1.addIDMatch({"run", IdentificationData::MatchId{1}});
  TEST_FALSE(p1 == p2);
  p2.addIDMatch({"run", IdentificationData::MatchId{1}});
  TEST_EQUAL(p1 != p2, false);
END_SECTION

START_SECTION(([EXTRA]meta info with copy constructor))
  BaseFeature p;
  p.setMetaValue(2,std::string("bla"));
  BaseFeature p2(p);
  TEST_EQUAL(p.getMetaValue(2), "bla")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
  p.setMetaValue(2,std::string("bluff"));
  TEST_EQUAL(p.getMetaValue(2), "bluff")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

START_SECTION(([EXTRA]meta info with assignment))
  BaseFeature p;
  p.setMetaValue(2,std::string("bla"));
  BaseFeature p2 = p;
  TEST_EQUAL(p.getMetaValue(2), "bla")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
  p.setMetaValue(2,std::string("bluff"));
  TEST_EQUAL(p.getMetaValue(2), "bluff")
  TEST_EQUAL(p2.getMetaValue(2), "bla")
END_SECTION

START_SECTION(([BaseFeature::QualityLess] bool operator()(const BaseFeature &left, const BaseFeature &right) const ))
  BaseFeature f1, f2;
  f1.setQuality((QualityType)0.94);
  f2.setQuality((QualityType)0.78);
  BaseFeature::QualityLess oql;

  TEST_EQUAL(oql(f1, f2), 0);
  TEST_EQUAL(oql(f2, f1), 1);
END_SECTION

START_SECTION(([BaseFeature::QualityLess] bool operator()(const BaseFeature &left, const QualityType &right) const ))
  BaseFeature f1, f2;
  f1.setQuality((QualityType)0.94);
  f2.setQuality((QualityType)0.78);
  BaseFeature::QualityType rhs = f1.getQuality();
  BaseFeature::QualityLess oql;

  TEST_EQUAL(oql(f1, rhs), 0);
  TEST_EQUAL(oql(f2, rhs), 1);
END_SECTION

START_SECTION(([BaseFeature::QualityLess] bool operator()(const QualityType& left, const BaseFeature& right) const))
  BaseFeature f1, f2;
  f1.setQuality((QualityType)0.94);
  f2.setQuality((QualityType)0.78);
  BaseFeature::QualityType lhs = f2.getQuality();
  BaseFeature::QualityLess oql;

  TEST_EQUAL(oql(lhs,f2), 0);
  TEST_EQUAL(oql(lhs,f1), 1);
END_SECTION

START_SECTION(([BaseFeature::QualityLess] bool operator()(const QualityType& left, const QualityType& right) const ))
  BaseFeature f1, f2;
  f1.setQuality((QualityType)0.94);
  f2.setQuality((QualityType)0.78);
  BaseFeature::QualityType lhs = f1.getQuality();
  BaseFeature::QualityType rhs = f2.getQuality();
  BaseFeature::QualityLess oql;

  TEST_EQUAL(oql(lhs,rhs), 0);
  TEST_EQUAL(oql(rhs,lhs), 1);
END_SECTION


START_SECTION((static AnnotationState annotationState(const PeptideIdentificationList& peptides)))
  PeptideIdentificationList ids;

  TEST_EQUAL(BaseFeature::annotationState(ids), BaseFeature::AnnotationState::FEATURE_ID_NONE);
  ids.resize(1);
  TEST_EQUAL(BaseFeature::annotationState(ids), BaseFeature::AnnotationState::FEATURE_ID_NONE);

  PeptideHit hit;
  hit.setSequence(AASequence::fromString("ABCDE"));
  ids[0].setHits(std::vector<PeptideHit>(1, hit));
  TEST_EQUAL(BaseFeature::annotationState(ids), BaseFeature::AnnotationState::FEATURE_ID_SINGLE);

  ids.resize(2);
  ids[1].setHits(std::vector<PeptideHit>(1, hit)); // same as first hit
  TEST_EQUAL(BaseFeature::annotationState(ids), BaseFeature::AnnotationState::FEATURE_ID_MULTIPLE_SAME);

  hit.setSequence(AASequence::fromString("KRGH"));
  ids[1].setHits(std::vector<PeptideHit>(1, hit)); // different to first hit
  TEST_EQUAL(BaseFeature::annotationState(ids), BaseFeature::AnnotationState::FEATURE_ID_MULTIPLE_DIVERGENT);
END_SECTION

START_SECTION((AnnotationState getAnnotationState(const IdentificationData& data) const))
{
  IdentificationData data;
  auto& run = data.addRun("search");
  IdentificationData::ScoreDefinition score;
  score.name = "score";
  run.setPrimaryScore(run.addScore(score));
  const auto source = run.addSource({});
  const auto add = [&](const std::string& sequence) {
    IdentificationData::MatchData match;
    match.representation = sequence;
    return run.addMatch(run.addIdentification(source, {}), match, {1.0});
  };
  const auto first = add("ABCDE"), same = add("ABCDE"), other = add("KRGH");

  BaseFeature feature;
  TEST_EQUAL(feature.getAnnotationState(data), BaseFeature::AnnotationState::FEATURE_ID_NONE);
  feature.addIDMatch({run.getUuid(), first});
  TEST_EQUAL(feature.getAnnotationState(data), BaseFeature::AnnotationState::FEATURE_ID_SINGLE);
  feature.addIDMatch({run.getUuid(), same});
  TEST_EQUAL(feature.getAnnotationState(data), BaseFeature::AnnotationState::FEATURE_ID_MULTIPLE_SAME);
  feature.addIDMatch({run.getUuid(), other});
  TEST_EQUAL(feature.getAnnotationState(data), BaseFeature::AnnotationState::FEATURE_ID_MULTIPLE_DIVERGENT);
  // a link to a match that the data does not have
  feature.addIDMatch({run.getUuid(), IdentificationData::MatchId{1000}});
  TEST_EXCEPTION(Exception::MissingInformation, feature.getAnnotationState(data));
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
