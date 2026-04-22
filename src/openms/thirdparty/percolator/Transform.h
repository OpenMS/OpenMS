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
#ifndef TRANSFORM_H_
#define TRANSFORM_H_

#include <cmath>

class Transform {
  public:
    Transform(double deltL = 0.0, double deltH = 0.0, bool dLt = false,
              bool dL = false) :
      deltaLow(deltL), deltaHigh(deltH), doLogit(dLt), doLog(dL) {
      ;
    }
    ~Transform() {
      ;
    }
    double operator()(double xx) {
      if ((!doLogit) && (!doLog)) {
        return xx;
      }
      if ((deltaLow > .0) || (deltaHigh > .0)) {
        if (doLogit) {
          xx *= (1. - deltaHigh - deltaLow);
        }
        xx += deltaLow;
      }
      if (doLogit) {
        return log(xx / (1. - xx));
      }
      return log(xx);
    }
  private:
    double deltaLow, deltaHigh;
    bool doLogit, doLog;
};

#endif /*TRANSFORM_H_*/
