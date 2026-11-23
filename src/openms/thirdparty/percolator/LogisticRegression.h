/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/

#ifndef LOGISTICREGRESSION_H_
#define LOGISTICREGRESSION_H_

#include <vector>

#include "BaseSpline.h"
#include "Array.h"

namespace OpenMS { namespace Internal { namespace Percolator {


class LogisticRegression : public BaseSpline {
  public:
    LogisticRegression(){};
    virtual ~LogisticRegression(){};
    void predict(const std::vector<double>& x, std::vector<double>& predict) {
      return BaseSpline::predict(x, predict);
    }
    // The knots are the scores xx after the transform of BaseSpline::setData;
    // equal knots are merged, and a warning is printed if some are nearly
    // equal (see countNearTies).
    void setData(const std::vector<double>& xx, const std::vector<double>& yy,
                 const std::vector<double>& mm);
  protected:
    // Merges adjacent equal knots (ascending, after the transform of
    // setData), with their y and m summed, unless fewer than 4 knots would
    // remain.
    static void mergeEqualKnots(std::vector<double>& xx,
                                std::vector<double>& yy,
                                std::vector<double>& mm);
    // The number of pairs of neighbouring knots (ascending) that are less
    // than 1e-9 times the median of the gaps around them apart: up to 5 gaps
    // ending at the lower knot and up to 5 starting at the upper one (of an
    // even number of gaps the larger of the two in the middle).
    static std::size_t countNearTies(const std::vector<double>& xx);
    virtual void calcPZW();
    virtual void initg();
    virtual void limitg();
    std::vector<double> y, m;
    static const double gRange;
    Array<double> p;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*LOGISTICREGRESSION_H_*/
