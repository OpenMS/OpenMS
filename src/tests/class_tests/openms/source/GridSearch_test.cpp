// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ML/GRIDSEARCH/GridSearch.h>

using namespace OpenMS;
using namespace std;

START_TEST(GridSearch, "$Id$")

    START_SECTION(GridSearch lambda)
    {
        auto evaluator = [](double i, const std::string& j, double k, double l)
        {
          return i+j.length()+k+l;
        };

        GridSearch<double, std::string, double, double> gs({1,3,5,2},{"foo","barz"},{2},{3});
        std::array<size_t, 4> bestParamIdx{{0u,0u,0u,0u}};
        TEST_EQUAL(gs.getNrCombos(), 8)
        gs.evaluate(evaluator, -1.0, bestParamIdx);
        TEST_EQUAL(get<0>(bestParamIdx),2)
        TEST_EQUAL(get<1>(bestParamIdx),1)
        TEST_EQUAL(get<2>(bestParamIdx),0)
        TEST_EQUAL(get<3>(bestParamIdx),0)

    }
    END_SECTION

    START_SECTION(GridSearch Functor)
        {
          struct Evaluator
          {
            double operator()(double i, const std::string& j, double k, double l)
            {
              return i+j.length()+k+l;
            };
          };

          GridSearch<double, std::string, double, double> gs({1,3,5,2},{"foo","barz"},{2},{3});
          std::array<size_t, 4> bestParamIdx{{0u,0u,0u,0u}};
          TEST_EQUAL(gs.getNrCombos(), 8)
          gs.evaluate(Evaluator(), -1.0, bestParamIdx);
          TEST_EQUAL(get<0>(bestParamIdx),2)
          TEST_EQUAL(get<1>(bestParamIdx),1)
          TEST_EQUAL(get<2>(bestParamIdx),0)
          TEST_EQUAL(get<3>(bestParamIdx),0)

        }
    END_SECTION


END_TEST
