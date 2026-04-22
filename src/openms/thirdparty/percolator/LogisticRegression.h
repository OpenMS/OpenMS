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

class LogisticRegression : public BaseSpline {
  public:
    LogisticRegression(){};
    virtual ~LogisticRegression(){};
    void predict(const std::vector<double>& x, std::vector<double>& predict) {
      return BaseSpline::predict(x, predict);
    }
    void setData(const std::vector<double>& xx, const std::vector<double>& yy,
                 const std::vector<double>& mm) {
      BaseSpline::setData(xx);
      y = yy;
      m = mm;
    }
  protected:
    virtual void calcPZW();
    virtual void initg();
    virtual void limitg();
    virtual void limitgamma();
    std::vector<double> y, m;
    static const double gRange;
    Array<double> p;
};

#endif /*LOGISTICREGRESSION_H_*/
