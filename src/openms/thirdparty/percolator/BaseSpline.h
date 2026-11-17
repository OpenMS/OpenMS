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

#ifndef BASESPLINE_H_
#define BASESPLINE_H_

#include <assert.h>
#include "Transform.h"
#include "PackedVector.h"

namespace OpenMS { namespace Internal { namespace Percolator {


class BaseSpline {
  public:
    BaseSpline(){};
    virtual ~BaseSpline(){};
    double splineEval(double xx);
    static double convergeEpsilon;
    static double stepEpsilon;
    static double weightSlope;
    static double scaleAlpha;
    void roughnessPenaltyIRLS();
    // Penalized IRLS fit for a given alpha, from the current gnew. Returns
    // false if the step size was still above stepEpsilon after 20 steps, or
    // if a step failed (see reinschStep); g is then the last iterate.
    bool iterativeReweightedLeastSquares(double alpha);
    void predict(const vector<double>& x, vector<double>& predict);
    void setData(const vector<double>& x);
    double predict(double xx) {
      return splineEval(xx);
    }
  protected:
    virtual void calcPZW() {}
    virtual void initg() {
      int n = static_cast<int>(x.size());
      g = PackedVector(n);
      gnew = PackedVector(n);
      w = PackedVector(n);
      z = PackedVector(n,0.5);
      gamma = PackedVector(n-2);
    }
    virtual void limitg() {}
    // One step of the Reinsch algorithm: solves
    // (R + Q' alpha W^-1 Q) gamma = Q' z and sets gnew = z - alpha W^-1 Q gamma.
    // Returns false, and changes nothing, if two knots are equal or a result
    // is not finite.
    bool reinschStep(double alpha);
    double evaluateSlope(double alpha);
    double alphaLinearSearchBA(double min_p, double max_p,
                               double p1, double p2,
                               double cv1, double cv2);
    Transform transf;

    PackedVector gnew, w, z;
    PackedVector g, gamma;
    vector<double> x;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*BASESPLINE_H_*/
