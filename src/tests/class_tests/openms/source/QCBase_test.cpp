// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/QC/QCBase.h>
#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////

using namespace OpenMS;

START_TEST(SpectraMap, "$Id$")

  QCBase::SpectraMap* ptr = nullptr;
  QCBase::SpectraMap* nulpt = nullptr;
  START_SECTION(QCBase::SpectraMap())
    {
      ptr = new QCBase::SpectraMap();
      TEST_NOT_EQUAL(ptr, nulpt)
    }
  END_SECTION

  START_SECTION(~QCBase::SpectraMap())
    {
      delete ptr;
    }
  END_SECTION
  
  MSExperiment exp;
  MSSpectrum spec1;
  spec1.setNativeID("XTandem::0");
  MSSpectrum spec2;
  spec2.setNativeID("XTandem::1");
  MSSpectrum spec3;
  spec3.setNativeID("XTandem::2");
  exp.setSpectra({spec1,spec2,spec3});
  
  START_SECTION(QCBase::SpectraMap::calculateMap(const MSExperiment& exp))
    QCBase::SpectraMap spec_map;
    spec_map.calculateMap(exp);
    ABORT_IF(spec_map.size() != 3);
    TEST_EQUAL(spec_map.at("XTandem::0"), 0);
    TEST_EQUAL(spec_map.at("XTandem::1"), 1);
    TEST_EQUAL(spec_map.at("XTandem::2"), 2);
    TEST_EXCEPTION(Exception::ElementNotFound, spec_map.at("XTandem::15"));
  END_SECTION

  START_SECTION(QCBase::SpectraMap::SpectraMap(const MSExperiment& exp))
    QCBase::SpectraMap spec_map(exp);
    TEST_EQUAL(spec_map.size(), 3);
  END_SECTION
  
  START_SECTION(QCBase::SpectraMap::empty())
    QCBase::SpectraMap spec_map;
    TEST_EQUAL(spec_map.empty(),true);
  END_SECTION

  START_SECTION(QCBase::SpectraMap::clear())
    QCBase::SpectraMap spec_map;
    spec_map.calculateMap(exp);
    TEST_EQUAL(spec_map.empty(),false);
    spec_map.clear();
    TEST_EQUAL(spec_map.empty(),true);
  END_SECTION
  
  START_SECTION(QCBase::SpectraMap::at(const String& identifier))
    NOT_TESTABLE;
  END_SECTION
  
  START_SECTION(QCBase::SpectraMap::size())
    NOT_TESTABLE;
  END_SECTION
  
END_TEST

