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
#ifndef MASSHANDLER_H_
#define MASSHANDLER_H_

#include<string>
using namespace std;

class MassHandler {
  public:
    MassHandler();
    virtual ~MassHandler();
    static void setMonoisotopicMass(bool mi) {
      monoisotopic = mi;
    }
    static double massDiff(double observedMass, double calculatedMass, unsigned int charge);
    static bool monoisotopic;
};

#endif /*MASSHANDLER_H_*/
