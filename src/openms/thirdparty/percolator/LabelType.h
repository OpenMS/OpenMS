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
#ifndef LABEL_TYPE_H_
#define LABEL_TYPE_H_

#include <iostream>

enum class LabelType {
  DECOY = -1,
  UNDEFINED = 0,
  TARGET = 1,
  PSEUDO_TARGET = 2  // used for RESET algorithm
};

inline std::ostream& operator<<(std::ostream& os, LabelType label) {
  os << static_cast<int>(label);
  return os;
}

#endif /*LABEL_TYPE_H_*/