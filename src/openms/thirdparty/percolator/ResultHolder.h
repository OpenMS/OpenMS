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
#ifndef RESULTHOLDER_H_
#define RESULTHOLDER_H_
#include <cmath>

using namespace std;

class ResultHolder {
   public:
    ResultHolder();
    ResultHolder(const double score, const double q, const double po,
                 const string& i, const string& pe = "", const string& p = "", const string& fn = "");
    virtual ~ResultHolder();
    double score, q, posterior;
    string id, pepSeq, prot, fileName;
    bool outputRT;
    double retentionTime = nan("");
};

bool operator>(const ResultHolder& one, const ResultHolder& other);
bool operator<(const ResultHolder& one, const ResultHolder& other);
ostream& operator<<(ostream& out, const ResultHolder& obj);

#endif /*RESULTHOLDER_H_*/
