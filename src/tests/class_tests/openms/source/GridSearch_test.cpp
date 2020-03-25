// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/MATH/MISC/GridSearch.h>

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
