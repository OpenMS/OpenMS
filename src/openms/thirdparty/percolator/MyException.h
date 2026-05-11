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

#ifndef MYEXCEPTION_HH
#define MYEXCEPTION_HH

#include <string>
#include <exception>
#include <sstream>

namespace OpenMS { namespace Internal { namespace Percolator {

using namespace std;

class MyException : public std::exception {
 public:
  MyException(const std::string &ss);
  MyException(const std::ostream &ss);
  ~MyException() throw();
  const char* what() const throw();
 
 protected:
  std::string msg;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif
