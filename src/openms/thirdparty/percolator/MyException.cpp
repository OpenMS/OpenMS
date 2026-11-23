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

#include "MyException.h"

#include <iostream>

namespace OpenMS { namespace Internal { namespace Percolator {


MyException::MyException(const string& ss)  : msg(ss) 
{
  
}

MyException::MyException(const ostream& ss)
{
   ostringstream x;
   x << ss.rdbuf();
   msg = std::string(x.str());
}


MyException::~MyException() throw()
{

}

const char* MyException::what() const throw()
{
  return msg.c_str(); 
}

void closeOutputFile(std::ofstream& out, const std::string& fileName) {
  out.close();
  if (out.fail()) {
    throw MyException("ERROR: Could not write the file " + fileName +
                      ". Check that the disk is not full and that the file "
                      "is writable.\n");
  }
}

void flushStandardOutput() {
  std::cout.flush();
  if (std::cout.fail()) {
    throw MyException("ERROR: Could not write the results to standard "
                      "output.\n");
  }
}

}}}  // namespace OpenMS::Internal::Percolator
