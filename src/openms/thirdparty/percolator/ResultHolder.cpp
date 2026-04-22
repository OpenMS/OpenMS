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
#include <iostream>
#include <string>
using namespace std;
#include "ResultHolder.h"

ResultHolder::ResultHolder() : score(0.0), q(1.0), posterior(1.0), pepSeq(""), prot(""), outputRT(false) {
}

ResultHolder::ResultHolder(const double sc, const double qq,
                           const double po, const string& i,
                           const string& pe, const string& p, const string& fn) : score(sc), q(qq), posterior(po), id(i), pepSeq(pe), prot(p), fileName(fn),  outputRT(false) {
}

ResultHolder::~ResultHolder() {
}


bool operator>(const ResultHolder& one, const ResultHolder& other) {
    return (one.score > other.score);
}

bool operator<(const ResultHolder& one, const ResultHolder& other) {
    return (one.score < other.score);
}

ostream& operator<<(ostream& out, const ResultHolder& obj) {
    out << obj.id << "\t";
    if (obj.fileName.size() > 0) {
        out << obj.fileName << "\t";
    }
    out << obj.score << "\t" << obj.q << "\t";
    out << obj.posterior << "\t" << obj.pepSeq << "\t" << obj.prot;
    if (obj.outputRT) {
        out << "\t" << obj.retentionTime;
    }
    return out;
}
